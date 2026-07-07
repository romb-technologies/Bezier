#include "Bezier/bezier.h"
#include "Bezier/coefficients.h"
#include "Bezier/utils.h"

#include <algorithm>
#include <limits>
#include <numeric>

#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/Polynomials>

using namespace Bezier;
namespace bu = Bezier::Utils;
namespace bc = Bezier::Coefficients;

std::vector<unsigned> Bezier::Utils::visvalingamWyatt(const PointVector& polyline)
{
  // Vector of indices sorted by contribution to the polyline shape (first and last contribute the most by default)
  // Initialized with all indices, taking care to put the first and the last at the start
  std::vector<unsigned> by_contribution(polyline.size());
  std::iota(by_contribution.begin(), by_contribution.end(), -1);
  by_contribution.front() += polyline.size();

  // Helper structure to keep track of each Points neigbours and contribution to the polyline shape
  struct Vertex
  {
    size_t prev, next;
    double contribution;
  };
  std::vector<Vertex> vertices(polyline.size());

  // Visvalingam-Whyatt measures contribution as an area between 3 consecutive Points
  auto area = [&polyline](unsigned idx1, unsigned idx2, unsigned idx3) {
    return std::fabs(bu::cross(polyline[idx2] - polyline[idx1], polyline[idx3] - polyline[idx1])) / 2;
  };
  auto cmp = [&vertices](unsigned idx1, unsigned idx2) {
    return vertices[idx1].contribution < vertices[idx2].contribution;
  };

  // Initialize vertices (NONE marks an endpoint's missing neighbour; never read back)
  constexpr size_t NONE = 0;
  vertices.front() = {NONE, 1, 0.0};
  for (unsigned k = 1; k + 1 < polyline.size(); k++)
    vertices[k] = {k - 1, k + 1, area(k - 1, k, k + 1)};
  vertices.back() = {polyline.size() - 2, NONE, 0.0};

  // Smallest contribution will be at the end of the vector
  for (auto it = by_contribution.rbegin(); it != by_contribution.rend() - 2; ++it)
  {
    // Select and move a Point with smallest current contribution
    std::iter_swap(it, std::min_element(it, by_contribution.rend() - 2, cmp));

    // Update previous and next Vertices (neighbours and contributions)
    auto prev = vertices[*it].prev;
    auto next = vertices[*it].next;
    vertices[prev].next = next;
    vertices[next].prev = prev;
    vertices[prev].contribution = area(vertices[prev].prev, prev, next);
    vertices[next].contribution = area(prev, next, vertices[next].next);
  }

  return by_contribution;
}

std::vector<double> Bezier::Utils::solvePolynomial(const Eigen::VectorXd& polynomial)
{
  // Trim trailing zero coefficients from the polynomial
  unsigned idx = polynomial.size();
  while (idx && std::fabs(polynomial(idx - 1)) < bu::epsilon)
    --idx;
  if (idx < 2) // Polynomial is a constant
    return {};

  std::vector<double> roots;
  Eigen::PolynomialSolver<double, Eigen::Dynamic>(polynomial.head(idx)).realRoots(roots);
  roots.erase(std::remove_if(roots.begin(), roots.end(), [](double t) { return t < 0 || t > 1; }), roots.end());
  return roots;
}

// VarPro fitter: softmax-reparameterised footpoints, ridge-regularised control-point solve, analytic Jacobian.
Curve Bezier::Utils::fitBezier(const PointVector& points, unsigned order)
{
  const unsigned M = points.size();
  const unsigned N = std::min<unsigned>(order + 1, M);
  if (N < 3)
    return Curve(PointVector{points.front(), points.back()}); // a 2-point (order-1) fit is just the chord

  Eigen::MatrixX2d P(M, 2);
  for (unsigned k{}; k < M; k++)
    P.row(k) = points[k];

  // Ridge VarPro: the control-point solve is regularised, min ||Phi C - Y||^2 + lambda ||C-Cref||^2,
  // by augmenting the system with [Phi; sqrt(lambda) I] and [Y; sqrt(lambda) Cref]. This bounds the
  // surplus high-order control points (the over-order Runge blow-up) and makes Phi full column rank.
  // Residuals therefore have 2*(M + NF) rows (M geometric + NF regularisation, per coordinate).
  struct CostFunctor : Eigen::DenseFunctor<double>
  {
    const Eigen::MatrixX2d& P;
    const unsigned N, M, NF;
    mutable Eigen::VectorXd cu_, t_; // cached u (M-2), derived t (M)
    mutable Eigen::MatrixXd Phi_;
    mutable Eigen::MatrixX2d Y_;
    mutable Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_;
    double sqlam_;          // sqrt(ridge lambda)
    Eigen::MatrixX2d Cref_; // straight-chord reference for the interior control points (NF x 2)

    CostFunctor(const Eigen::MatrixX2d& P_, unsigned N_)
        : DenseFunctor<double>(P_.rows() - 2, 2 * (P_.rows() + N_ - 2)), P(P_), N(N_), M(P_.rows()), NF(N_ - 2)
    {
      t_.resize(M);
      t_(0) = 0.0;
      Phi_.resize(M + NF, NF);
      Y_.resize(M + NF, 2);
      Cref_.resize(NF, 2);
      for (unsigned k{}; k < NF; k++)
        Cref_.row(k) = P.row(0) + double(k + 1) / (N - 1) * (P.row(M - 1) - P.row(0));
      sqlam_ = std::sqrt(1e-5); // lambda=1e-5: keeps Phi full-rank and damps Runge blow-up without biasing the fit
      qr_.setThreshold(1e-7);
    }

    void prepare(const Eigen::VectorXd& u) const
    {
      if (cu_.size() == u.size() && cu_ == u)
        return;
      // M-1 interval log-gaps = [u (free), 0 (anchored last)]; softmax for overflow-safety.
      const double lmax = std::max(u.maxCoeff(), 0.0);
      Eigen::VectorXd g(M - 1);
      g << (u.array() - lmax).exp().matrix(), std::exp(-lmax);     // softmax-normalised
      std::partial_sum(g.data(), g.data() + M - 1, t_.data() + 1); // cumulative gaps
      t_ /= t_(M - 1);                                             // normalise -> t_(M-1) = 1 exactly
      Eigen::MatrixXd A = bu::powMatrix(t_, N) * bc::bernstein(N);
      // Augment with the ridge rows: Phi_ = [Phi; sqrt(lambda) I], Y_ = [Y; sqrt(lambda) Cref].
      Phi_ << A.middleCols(1, NF), sqlam_ * Eigen::MatrixXd::Identity(NF, NF);
      Y_ << P - A.col(0) * P.row(0) - A.col(N - 1) * P.row(M - 1), sqlam_ * Cref_;
      qr_.compute(Phi_);
      cu_ = u;
    }

    int operator()(const Eigen::VectorXd& u, Eigen::VectorXd& fvec) const
    {
      prepare(u);
      // (I - P_U) Y_aug residuals, x/y interleaved: top M geometric, bottom NF ridge
      Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>>(fvec.data(), M + NF, 2) =
          Phi_ * qr_.solve(Y_) - Y_;
      return 0;
    }

    int df(const Eigen::VectorXd& u, Eigen::MatrixXd& jac) const
    {
      prepare(u);
      Eigen::MatrixX2d cp(N, 2);
      cp << P.row(0), qr_.solve(Y_), P.row(M - 1); // ridge-regularised interior CPs
      Curve curve(cp);

      // Geometric Jacobian w.r.t. the interior t, augmented (2*(M+NF) x (M-2)). Only the top M rows
      // depend on t directly (via C'(t_l)); the (I-P_U) of the augmented system carries the ridge.
      Eigen::MatrixXd Jt(values(), M - 2);
      Eigen::MatrixXd Ur = qr_.householderQ() * Eigen::MatrixXd::Identity(M + NF, qr_.rank());
      for (unsigned l{1}; l + 1 < M; l++)
      {
        Eigen::VectorXd proj = -Ur * Ur.row(l).transpose(); // (I - P_U) e_l in the augmented space
        proj(l) += 1.0;
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>>(Jt.col(l - 1).data(), M + NF, 2) =
            proj * curve.derivativeAt(t_(l)).transpose();
      }

      // Chain rule over the M-2 free gaps (j = 1..M-2; the anchored last gap is not a variable):
      // W[l-1,j-1] = (t_j - t_{j-1})([j<=l] - t_l), then J_u = Jt * W  (2*(M+NF) x (M-2)).
      Eigen::MatrixXd W(M - 2, M - 2);
      for (unsigned l{1}; l + 1 < M; l++)
        for (unsigned j{1}; j + 1 < M; j++)
          W(l - 1, j - 1) = (t_(j) - t_(j - 1)) * ((j <= l ? 1.0 : 0.0) - t_(l));

      jac = Jt * W;
      return 0;
    }
  };

  // Centripetal initialisation: u_j = 0.5 * log(|dP_{j+1}| / |dP_last|), relative to the anchored last interval.
  const double d_last = std::max((P.row(M - 1) - P.row(M - 2)).norm(), 1e-12);
  Eigen::VectorXd u(M - 2);
  for (Eigen::Index j{}; j + 2 < M; j++)
    u(j) = 0.5 * std::log(std::max((P.row(j + 1) - P.row(j)).norm(), 1e-12) / d_last);

  CostFunctor functor(P, N);
  Eigen::LevenbergMarquardt<CostFunctor> lm(functor);
  lm.minimize(u);
  functor.prepare(u);

  Eigen::MatrixX2d cp(N, 2);
  cp << P.row(0), functor.qr_.solve(functor.Y_), P.row(M - 1);
  return Curve(cp);
}
