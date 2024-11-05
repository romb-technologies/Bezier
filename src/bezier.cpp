#include "Bezier/bezier.h"
#include "Bezier/coefficients.h"
#include "Bezier/utils.h"

#include <numeric>

#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/NumericalDiff>

using namespace Bezier;
namespace bu = Bezier::Utils;
namespace bc = Bezier::Coefficients;

Curve::Curve(Eigen::MatrixX2d points) : N_(points.rows()), control_points_(std::move(points)) {}

Curve::Curve(const std::vector<Point>& points) : N_(points.size()), control_points_(N_, 2)
{
  for (unsigned k{}; k < N_; k++)
    control_points_.row(k) = points[k];
}

Curve::Curve(const Curve& curve) : Curve(curve.control_points_) {}

Curve& Curve::operator=(const Curve& curve)
{
  control_points_ = curve.control_points_;
  clearCache();
  return *this;
}

unsigned Curve::order() const { return N_ - 1; }

std::vector<Point> Curve::controlPoints() const
{
  std::vector<Point> points(N_);
  for (unsigned k{}; k < N_; k++)
    points[k] = control_points_.row(k);
  return points;
}

Point Curve::controlPoint(unsigned idx) const { return control_points_.row(idx); }

std::pair<Point, Point> Curve::endPoints() const { return {control_points_.row(0), control_points_.row(N_ - 1)}; }

std::vector<Point> Curve::polyline() const
{
  return polyline(boundingBox().diagonal().norm() / 1000);
}

std::vector<Point> Curve::polyline(double flatness) const
{
  if (cached_polyline_ && std::fabs(cached_polyline_flatness_ - flatness) < bu::epsilon)
    return cached_polyline_.value();

  cached_polyline_flatness_ = flatness;
  cached_polyline_t_.emplace({0.0});
  cached_polyline_.emplace({Point(control_points_.row(0))});

  std::vector<std::tuple<Eigen::MatrixX2d, double, double>> subcurves;
  subcurves.emplace_back(control_points_, 0.0, 1.0);

  while (!subcurves.empty())
  {
    auto [cp, t1, t2] = std::move(subcurves.back());
    subcurves.pop_back();

    if (bu::maxDeviation(cp) <= flatness)
    {
      cached_polyline_t_->emplace_back(t2);
      cached_polyline_->emplace_back(cp.row(N_ - 1));
    }
    else
    {
      subcurves.emplace_back(bc::rightSplit(N_) * cp, (t1 + t2) / 2, t2);
      subcurves.emplace_back(bc::leftSplit(N_) * cp, t1, (t1 + t2) / 2);
    }
  }

  return cached_polyline_.value();
}

double Curve::length() const { return length(1.0); }

double Curve::length(double t) const
{
  if (t < 0.0 || t > 1.0)
    throw std::logic_error{"Length can only be calculated for t within [0.0, 1.0] range."};

  if (cached_chebyshev_polynomial_)
    return bu::evaluateChebyshev(t, cached_chebyshev_polynomial_.value());
  auto& chebyshev = cached_chebyshev_polynomial_.emplace();

  constexpr unsigned START_LOG_N = 10;
  unsigned log_n = START_LOG_N;
  Eigen::VectorXd derivative_cache(bu::exp2(log_n) + 1);
  derivative_cache.head(2) << derivativeAt(1.0).norm(), derivativeAt(0.0).norm();

  auto updateDerivativeCache = [this, &derivative_cache](unsigned n) {
    auto derFunc = [n, this](int k) { return derivativeAt((1 + std::cos((2 * k + 1) * M_PI / n)) / 2).norm(); };
    derivative_cache.conservativeResize(n + 1);
    derivative_cache.tail(n / 2) = Eigen::VectorXd::NullaryExpr(n / 2, derFunc);
  };

  for (unsigned k{1}; k < log_n; k++)
    updateDerivativeCache(bu::exp2(k));

  do
  {
    unsigned n{bu::exp2(log_n++)};
    Eigen::VectorXd coeff(2 * n);
    coeff(0) = derivative_cache(0);
    coeff(n) = derivative_cache(1);
    updateDerivativeCache(n);

    for (unsigned k{1}; k < log_n; k++)
    {
      Eigen::ArrayXi lin_spaced = Eigen::ArrayXi::LinSpaced(bu::exp2(k - 1), 0, bu::exp2(k - 1) - 1);
      Eigen::ArrayXi index_c = bu::exp2(log_n - (k + 1)) + lin_spaced * bu::exp2(log_n - k);
      Eigen::ArrayXi index_dc = bu::exp2(k - 1) + 1 + lin_spaced;
      for (unsigned i{}; i < lin_spaced.size(); i++)
        coeff(index_c(i)) = coeff(2 * n - index_c(i)) = derivative_cache(index_dc(i)) / n;
    }

    Eigen::VectorXcd fft_out;
    Eigen::FFT<double>().fwd(fft_out, coeff);
    chebyshev.resize(n);
    chebyshev << 0, (fft_out.real().head(n - 1) - fft_out.real().segment(2, n - 1))
                        .cwiseQuotient(Eigen::VectorXd::LinSpaced(n - 1, 4, 4 * (n - 1)));
  } while (std::fabs(chebyshev.tail(1).value()) > bu::epsilon);

  // Trim trailing zero coefficients
  unsigned idx = chebyshev.size();
  while (std::fabs(chebyshev(idx - 1)) < bu::epsilon)
    --idx;
  chebyshev.conservativeResize(idx);

  chebyshev(0) -= bu::evaluateChebyshev(0, chebyshev);
  return bu::evaluateChebyshev(t, chebyshev);
}

double Curve::length(double t1, double t2) const { return length(t2) - length(t1); }

double Curve::step(double t, double dS) const
{
  if (std::fabs(dS) < bu::epsilon) // no-op
    return t;

  double s_t = length(t);

  std::pair<double, double> lbracket, rbracket, guess{t, 0.0};
  if (dS < 0)
  {
    lbracket = {0.0, -s_t};
    if (dS < lbracket.second + bu::epsilon) // out-of-scope
      return 0.0;
    rbracket = guess;
  }
  else // s > 0
  {
    rbracket = {1.0, length() - s_t};
    if (dS > rbracket.second - bu::epsilon) // out-of-scope
      return 1.0;
    lbracket = guess;
  }

  while (std::fabs(guess.second - dS) > bu::epsilon)
  {
    // Halley's method
    double f = guess.second - dS;
    double f_d = derivativeAt(guess.first).norm();
    double f_d2 = derivativeAt(2, guess.first).norm();
    guess.first -= (2 * f * f_d) / (2 * f_d * f_d - f * f_d2);

    // root bracketing, if not in bounds, use bisection method
    if (guess.first <= lbracket.first || guess.first >= rbracket.first)
      guess.first = (lbracket.first + rbracket.first) / 2;

    if (rbracket.first - lbracket.first < bu::epsilon)
      break;

    guess.second = length(guess.first) - s_t;
    (guess.second < dS ? lbracket : rbracket) = guess;
  }

  return guess.first;
}

void Curve::reverse()
{
  control_points_ = control_points_.colwise().reverse().eval();
  clearCache();
}

void Curve::setControlPoint(unsigned idx, const Point& point)
{
  control_points_.row(idx) = point;
  clearCache();
}

void Curve::raiseOrder()
{
  control_points_ = bc::raiseOrder(N_) * control_points_;
  clearCache();
}

void Curve::lowerOrder()
{
  if (N_ == 2)
    throw std::logic_error{"Cannot further reduce the order of curve."};
  control_points_ = bc::lowerOrder(N_) * control_points_;
  clearCache();
}

Point Curve::valueAt(double t) const
{
  return N_ == 0 ? Point(0, 0) : (bu::powVector(t, N_) * bc::bernstein(N_) * control_points_).transpose();
}

Eigen::MatrixX2d Curve::valueAt(const std::vector<double>& t_vector) const
{
  auto t_map = Eigen::Map<const Eigen::VectorXd>(t_vector.data(), t_vector.size());
  return bu::powMatrix(t_map, N_) * bc::bernstein(N_) * control_points_;
}

double Curve::curvatureAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);

  return bu::cross(d1, d2) / bu::pow(d1.norm(), 3);
}

double Curve::curvatureDerivativeAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);
  Vector d3 = derivativeAt(3, t);

  return (d1.squaredNorm() * bu::cross(d1, d3) - 3 * d1.dot(d2) * bu::cross(d1, d2)) / bu::pow(d1.norm(), 5);
}

Vector Curve::tangentAt(double t, bool normalize) const
{
  Vector d = derivativeAt(t);
  return normalize && d.norm() > 0 ? d.normalized() : d;
}

Vector Curve::normalAt(double t, bool normalize) const
{
  Vector tangent = tangentAt(t, normalize);
  return {-tangent.y(), tangent.x()};
}

const Curve& Curve::derivative() const
{
  if (!cached_derivative_)
    cached_derivative_ = std::make_unique<const Curve>(
        (N_ - 1) * (control_points_.bottomRows(N_ - 1) - control_points_.topRows(N_ - 1)));
  return *cached_derivative_;
}

const Curve& Curve::derivative(unsigned n) const
{
  auto* nth_derivative = this;
  for (unsigned k{}; k < n; k++)
    nth_derivative = &nth_derivative->derivative();
  return *nth_derivative;
}

Vector Curve::derivativeAt(double t) const { return derivative().valueAt(t); }

Vector Curve::derivativeAt(unsigned n, double t) const { return derivative(n).valueAt(t); }

std::vector<double> Curve::roots() const
{
  if (cached_roots_)
    return cached_roots_.value();

  Eigen::MatrixX2d bezier_polynomial = bc::bernstein(N_) * control_points_;
  return cached_roots_.emplace(bu::concatenate(bu::solvePolynomial(bezier_polynomial.col(0)), //
                                               bu::solvePolynomial(bezier_polynomial.col(1))));
}

std::vector<double> Curve::extrema() const { return derivative().roots(); }

BoundingBox Curve::boundingBox() const
{
  if (cached_bounding_box_)
    return cached_bounding_box_.value();

  auto extremes = valueAt(extrema());
  extremes.conservativeResize(extremes.rows() + 2, Eigen::NoChange);
  extremes.bottomRows<2>() << control_points_.row(0), control_points_.row(N_ - 1);
  return cached_bounding_box_.emplace(extremes.colwise().minCoeff(), extremes.colwise().maxCoeff());
}

std::vector<Curve> Curve::splitCurve(const std::vector<double>& t) const
{
  auto sorted_t = t;
  std::sort(sorted_t.begin(), sorted_t.end());
  std::vector<Curve> subcurves;
  subcurves.reserve(sorted_t.size() + 1);
  auto leftover_cp = control_points_;
  for (unsigned k{}; k < sorted_t.size(); k++)
  {
    subcurves.emplace_back(bc::leftSplit(N_, sorted_t[k]) * leftover_cp);
    leftover_cp = bc::rightSplit(N_, sorted_t[k]) * leftover_cp;
    std::for_each(sorted_t.begin() + k + 1, sorted_t.end(), [t = sorted_t[k]](double& x) { x = (x - t) / (1 - t); });
  }
  subcurves.emplace_back(std::move(leftover_cp));
  return subcurves;
}

std::vector<Curve> Curve::splitCurve(double t) const
{
  return {Curve(bc::leftSplit(N_, t) * control_points_), Curve(bc::rightSplit(N_, t) * control_points_)};
}

std::vector<Point> Curve::intersections(const Curve& curve) const
{
  std::vector<std::pair<Eigen::MatrixX2d, Eigen::MatrixX2d>> cp_pairs;
  if (!control_points_.isApprox(curve.control_points_))
    cp_pairs.emplace_back(control_points_, curve.control_points_);
  else
  {
    // If self-similar, split curve into subcurves at extremas and compare each pair of distinct subcurves
    auto subcurves = splitCurve(extrema());
    for(unsigned k{}; k < subcurves.size(); k++)
      for(unsigned i{k + 1}; i < subcurves.size(); i++)
        cp_pairs.emplace_back(subcurves[k].control_points_, subcurves[i].control_points_);
  }

  auto insertPairs = [&cp_pairs](const auto& scp1, const auto& scp2) {
    for (const auto& cp1 : scp1)
      for (const auto& cp2 : scp2)
        cp_pairs.emplace_back(cp1, cp2);
  };

  auto splitCP = [](const Eigen::MatrixX2d& cp) -> std::array<Eigen::MatrixX2d, 2> {
    return {bc::leftSplit(cp.rows()) * cp, bc::rightSplit(cp.rows()) * cp};
  };

  std::vector<Point> intersections;
  auto insertIntersection = [&intersections](const Eigen::MatrixX2d& cp1, const Eigen::MatrixX2d& cp2) {
    // Intersection of two line segments (Victor Lecomte - Handbook of geometry for competitive programmers)
    auto a1 = cp1.row(0), a2 = cp1.bottomRows<1>();
    auto b1 = cp2.row(0), b2 = cp2.bottomRows<1>();
    double oa = bu::cross(b2 - b1, a1 - b1);
    double ob = bu::cross(b2 - b1, a2 - b1);
    double oc = bu::cross(a2 - a1, b1 - a1);
    double od = bu::cross(a2 - a1, b2 - a1);

    // If intersection exists, insert it into solution vector
    if (oa * ob < 0 && oc * od < 0)
      intersections.emplace_back((a1 * ob - a2 * oa) / (ob - oa));
  };

  while (!cp_pairs.empty())
  {
    auto [cp1, cp2] = std::move(cp_pairs.back());
    cp_pairs.pop_back();

    BoundingBox bbox1(cp1.colwise().minCoeff(), cp1.colwise().maxCoeff());
    BoundingBox bbox2(cp2.colwise().minCoeff(), cp2.colwise().maxCoeff());
    if (!bbox1.intersects(bbox2))
      continue; // no intersection, cheap check

    // Split each curve until both are flat enough to be represented as line segment
    if (bu::maxDeviation(cp1) < bu::epsilon && bu::maxDeviation(cp2) < bu::epsilon)
      insertIntersection(cp1, cp2);
    else if (bbox1.diagonal().norm() < bu::epsilon)
      insertPairs(std::array{cp1}, splitCP(cp2));
    else if (bbox2.diagonal().norm() < bu::epsilon)
      insertPairs(splitCP(cp1), std::array{cp2});
    else
      insertPairs(splitCP(cp1), splitCP(cp2));
  }

  return intersections;
}

double Curve::projectPoint(const Point& point) const
{
  if (!cached_projection_polynomial_const_ || !cached_projection_polynomial_derivative_)
  {
    Eigen::MatrixX2d curve_polynomial = bc::bernstein(N_) * control_points_;
    Eigen::MatrixX2d derivate_polynomial = bc::bernstein(N_ - 1) * derivative().control_points_;

    Eigen::VectorXd polynomial_part = Eigen::VectorXd::Zero(2 * N_ - 2);
    for (unsigned k{}; k < N_; k++)
      polynomial_part.middleRows(k, N_ - 1) += derivate_polynomial * curve_polynomial.row(k).transpose();

    cached_projection_polynomial_const_.emplace(std::move(polynomial_part));
    cached_projection_polynomial_derivative_.emplace(std::move(derivate_polynomial));
  }

  Eigen::VectorXd polynomial = cached_projection_polynomial_const_.value();
  polynomial.topRows(N_ - 1) -= cached_projection_polynomial_derivative_.value() * point.transpose();

  double min_t{0.0}, min_dist{bu::dist(point, valueAt(0.0))};

  for (double t : bu::concatenate(bu::solvePolynomial(polynomial), {1.0}))
    if (double dist = bu::dist(point, valueAt(t)); dist < min_dist)
      std::tie(min_t, min_dist) = std::make_pair(t, dist);
  return min_t;
}

double Curve::distance(const Point& point) const { return bu::dist(point, valueAt(projectPoint(point))); }

void Curve::applyContinuity(const Curve& curve, const std::vector<double>& beta_coeffs)
{
  unsigned c_order = beta_coeffs.size();
  Eigen::Map<const Eigen::VectorXd> beta(beta_coeffs.data(), beta_coeffs.size());

  // pascal triangle matrix (binomial coefficients) - rowwise
  Eigen::MatrixXd pascal_matrix = Eigen::MatrixXd::Zero(c_order + 1, c_order + 1);
  pascal_matrix.diagonal(1).setLinSpaced(1, c_order);
  pascal_matrix = pascal_matrix.exp();

  // inverse of pascal matrix, i.e., pascal matrix with alternating signs - colwise
  Eigen::MatrixXd pascal_alternating_matrix = pascal_matrix.transpose().inverse();

  // https://en.wikipedia.org/wiki/Bell_polynomials -> equivalent to equations of geometric continuity
  Eigen::MatrixXd bell_matrix = Eigen::MatrixXd::Zero(c_order + 1, c_order + 1);
  bell_matrix(0, c_order) = 1;
  for (unsigned k{}; k < c_order; k++)
    bell_matrix.block(1, c_order - k - 1, k + 1, 1) = bell_matrix.block(0, c_order - k, k + 1, k + 1) *
                                                      pascal_matrix.col(k).head(k + 1).cwiseProduct(beta.head(k + 1));

  // derivatives of given curve
  Eigen::Matrix2Xd derivatives(2, c_order + 1);
  for (unsigned k{}; k < c_order + 1; k++)
    derivatives.col(k) = curve.derivative(k).control_points_.bottomRows<1>().transpose();

  // based on the beta coefficients and geometric continuity equations, calculate new derivatives
  Eigen::MatrixXd new_derivatives = (derivatives * bell_matrix).rowwise().reverse().transpose();

  // diagonal: (N-1)! / (N-k-1)!
  std::function<double(int)> permFunc = [x = 1. / N_, N = N_](int k) mutable { return x *= N - k; };
  Eigen::MatrixXd permutation_matrix = Eigen::VectorXd::NullaryExpr(c_order + 1, permFunc).asDiagonal();

  // calculate new control points
  control_points_.topRows(c_order + 1) = (permutation_matrix * pascal_alternating_matrix).inverse() * new_derivatives;
  clearCache();
}

Curve Curve::offsetCurve(const Curve& curve, double offset, unsigned order)
{
  std::vector<Point> offset_polyline = curve.polyline();
  for (unsigned k{}; k < offset_polyline.size(); k++)
    offset_polyline[k] += offset * curve.normalAt((*curve.cached_polyline_t_)[k]);
  return fromPolyline(offset_polyline, order ? order : curve.order() + 1);
}

Curve Curve::joinCurves(const Curve& curve1, const Curve& curve2, unsigned order)
{
  if (order == 1)
    return Curve(std::vector<Point>{curve1.control_points_.row(0), curve2.control_points_.row(curve2.N_ - 1)});
  return fromPolyline(bu::concatenate(curve1.polyline(), curve2.polyline()),
                      order ? order : curve1.order() + curve2.order());
}

Curve Curve::fromPolyline(const std::vector<Point>& polyline, unsigned order)
{
  const unsigned N = std::min(order ? order + 1 : polyline.size(), polyline.size());

  if (polyline.size() < 2)
    throw std::logic_error{"Polyline must have at least two points."};
  if (N == 2)
    return Curve(std::vector{polyline.front(), polyline.back()});

  // Sort the polyline points by their contribution to the Visvalingam-Whyatt
  // simplification algorithm, and keep the N most contributing points in original order.
  auto vw = bu::visvalingamWyatt(polyline);
  std::sort(vw.begin(), vw.begin() + N);

  // Divide polyline into subparts where the simplified polyline points are located.
  std::vector<std::vector<Point>> subpolylines;
  subpolylines.reserve(N - 1);
  subpolylines.emplace_back(std::vector{polyline.front()});
  for (unsigned k{1}; k + 1 < polyline.size(); k++)
  {
    subpolylines.back().emplace_back(polyline[k]);
    if (std::binary_search(vw.begin(), vw.begin() + N, k))
      subpolylines.emplace_back(std::vector{polyline[k]});
  }
  subpolylines.back().emplace_back(polyline.back());

  // Initialize vector t where each element represents a normalized cumulative
  // distance between consecutive simplified points along the simplified polyline.
  Eigen::VectorXd t(N);
  Eigen::MatrixX2d P(N, 2);
  for (unsigned k{}; k < N; k++)
  {
    P.row(k) = polyline[vw[k]];
    t(k) = k == 0 ? 0 : t(k - 1) + bu::dist(P.row(k), P.row(k - 1));
  }
  t /= t(N - 1);

  // Compute the control points for a Bezier curve such that it passes through
  // the simplified polyline points at parameter t.
  auto getCurve = [&P, M = bc::bernstein(N)](const Eigen::VectorXd& t) {
    Eigen::MatrixXd T = bu::powMatrix(t, t.size());
    return Curve(M.inverse() * (T.transpose() * T).inverse() * T.transpose() * P);
  };

  // Cost functor calculates RMSD and length difference for each subcurve/subpolyline
  // divided at parameter t, where C(t_i) = P_i.
  struct CostFunctor : public Eigen::DenseFunctor<double>
  {
    using GetCurveFun = std::function<Curve(const Eigen::VectorXd&)>;
    GetCurveFun getCurve;
    std::vector<std::vector<Point>> subpolylines;

    CostFunctor(int N, GetCurveFun getCurve, std::vector<std::vector<Point>> subpolylines)
        : DenseFunctor<double>(N - 2, 2 * N - 2), getCurve(std::move(getCurve)), subpolylines(std::move(subpolylines))
    {
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
    {
      auto curve = getCurve((Eigen::VectorXd(inputs() + 2) << 0, x, 1).finished());
      auto subcurves = curve.splitCurve(std::vector<double>(x.data(), x.data() + inputs()));
      for (int k = 0; k <= inputs(); k++)
      {
        auto polyline = subcurves[k].polyline();
        auto fun = [&](double acc, const Point& p) { return acc + bu::pow(bu::dist(subpolylines[k], p), 2); };
        fvec(k) = std::sqrt(std::accumulate(polyline.begin(), polyline.end(), 0.0, fun) / polyline.size());
        fvec(values() / 2 + k) = std::fabs(bu::polylineLength(polyline) - bu::polylineLength(subpolylines[k]));
      }
      return 0;
    }
  };

  // Use Levenberg-Marquardt optimization to find the control points that minimize
  // the RMSD and length difference of the Bezier curve and the simplified polyline.
  CostFunctor costFun(N, getCurve, subpolylines);
  Eigen::NumericalDiff<CostFunctor> num_diff(costFun);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<CostFunctor>> lm(num_diff);
  Eigen::VectorXd x = t.segment(1, N - 2);
  lm.minimize(x);

  return getCurve((Eigen::VectorXd(N) << 0, x, 1).finished());
}

void Curve::clearCache()
{
  N_ = control_points_.rows();
  cached_derivative_.reset();
  cached_roots_.reset();
  cached_bounding_box_.reset();
  cached_polyline_.reset();
  cached_polyline_t_.reset();
  cached_projection_polynomial_const_.reset();
  cached_projection_polynomial_derivative_.reset();
  cached_chebyshev_polynomial_.reset();
}
