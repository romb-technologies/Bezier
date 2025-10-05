#include "Bezier/bezier.h"
#include "Bezier/coefficients.h"
#include "Bezier/declarations.h"
#include "Bezier/utils.h"

#include <numeric>

#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Bezier;
namespace bu = Bezier::Utils;
namespace bc = Bezier::Coefficients;

///// Curve::Curve

Curve::Curve(Eigen::MatrixX2d points) : N_(points.rows()), control_points_(std::move(points)) {}

Curve::Curve(const PointVector& points)
    : control_points_(Eigen::Index(points.size()), Eigen::Index(2)), N_(points.size())
{
  for (unsigned k = 0; k < N_; k++)
    control_points_.row(k) = points[k];
}

Curve::Curve(const Curve& curve) : Curve(curve.control_points_) {}

Curve& Curve::operator=(const Curve& curve)
{
  N_ = control_points_.rows();
  control_points_ = curve.control_points_;
  cache.clear();
  return *this;
}

unsigned Curve::order() const { return N_ - 1; }

PointVector Curve::controlPoints() const
{
  PointVector points(N_);
  for (unsigned k = 0; k < N_; k++)
    points[k] = control_points_.row(k);
  return points;
}

Point Curve::controlPoint(unsigned idx) const { return control_points_.row(idx); }

std::pair<Point, Point> Curve::endPoints() const { return {control_points_.row(0), control_points_.row(N_ - 1)}; }

PointVector Curve::polyline() const { return polyline(boundingBox().diagonal().norm() / 1000); }

PointVector Curve::polyline(double flatness) const
{
  if (cache.polyline && std::fabs(cache.polyline_flatness - flatness) < bu::epsilon)
    return cache.polyline.value();

  cache.polyline_flatness = flatness;
  cache.polyline_t.emplace({0.0});
  cache.polyline.emplace({Point(control_points_.row(0))});

  std::vector<std::tuple<Eigen::MatrixX2d, double, double>> subcurves;
  subcurves.emplace_back(control_points_, 0.0, 1.0);

  while (!subcurves.empty())
  {
    auto [cp, t1, t2] = std::move(subcurves.back());
    subcurves.pop_back();

    if (bu::maxDeviation(cp) <= flatness)
    {
      cache.polyline_t->emplace_back(t2);
      cache.polyline->emplace_back(cp.row(N_ - 1));
    }
    else
    {
      subcurves.emplace_back(bc::rightSplit(N_) * cp, (t1 + t2) / 2, t2);
      subcurves.emplace_back(bc::leftSplit(N_) * cp, t1, (t1 + t2) / 2);
    }
  }

  return cache.polyline.value();
}

ParamVector Curve::polylineParams() const { return polylineParams(boundingBox().diagonal().norm() / 1000); }

ParamVector Curve::polylineParams(double flatness) const
{
  polyline(flatness);
  return cache.polyline_t.value();
}

double Curve::length() const { return length(1.0); }

double Curve::length(double t) const
{
  if (t < 0.0 || t > 1.0)
    throw std::logic_error{"Length can only be calculated for t within [0.0, 1.0] range."};

  if (!cache.chebyshev_polynomial)
  {
    constexpr unsigned START_LOG_N = 10;
    unsigned log_n = START_LOG_N - 1;
    unsigned n = bu::exp2(START_LOG_N - 1);

    Eigen::VectorXd derivative_cache(2 * n + 1);
    auto updateDerivativeCache = [this, &derivative_cache](double n) {
      derivative_cache.conservativeResize(n + 1);
      derivative_cache.tail(n / 2) =
          ((1 + Eigen::cos(Eigen::ArrayXd::LinSpaced(n / 2, 1, n - 1) * M_PI / n)) / 2).unaryExpr([this](double t) {
            return derivativeAt(t).norm();
          });
    };

    derivative_cache.head(2) << derivativeAt(1.0).norm(), derivativeAt(0.0).norm();
    for (unsigned k = 2; k <= n; k *= 2)
      updateDerivativeCache(k);

    Eigen::VectorXd chebyshev;
    Eigen::FFT<double> fft;
    Eigen::VectorXcd fft_out;
    do
    {
      n *= 2;
      log_n++;
      updateDerivativeCache(n);

      unsigned N = 2 * n;
      Eigen::VectorXd coeff(N);
      coeff(0) = derivative_cache(0);
      coeff(n) = derivative_cache(1);

      for (unsigned k = 1; k <= log_n; k++)
      {
        auto lin_spaced = Eigen::ArrayXi::LinSpaced(bu::exp2(k - 1), 0, bu::exp2(k - 1) - 1);
        auto index_c = bu::exp2(log_n + 1 - (k + 1)) + lin_spaced * bu::exp2(log_n + 1 - k);
        auto index_dc = bu::exp2(k - 1) + 1 + lin_spaced;
        // TODO: make use of slicing & indexing in Eigen3.4
        // coeff(index_c) = coeff(N - index_c) = derivative_cache(index_dc) / n;
        for (unsigned i = 0; i < lin_spaced.size(); i++)
          coeff(index_c(i)) = coeff(N - index_c(i)) = derivative_cache(index_dc(i)) / n;
      }

      fft.fwd(fft_out, coeff);
      chebyshev = (fft_out.real().head(n - 1) - fft_out.real().segment(2, n - 1)).array() /
                  Eigen::ArrayXd::LinSpaced(n - 1, 4, 4 * (n - 1));
    } while (std::fabs(chebyshev.tail<1>()[0]) > bu::epsilon * 1e-2);

    unsigned cut = 0;
    while (std::fabs(chebyshev(cut)) > bu::epsilon * 1e-2)
      cut++;
    cache.chebyshev_polynomial.emplace(cut + 1);
    (*cache.chebyshev_polynomial) << 0, chebyshev.head(cut);
    (*cache.chebyshev_polynomial)(0) = -bu::evaluateChebyshev(0, *cache.chebyshev_polynomial);
  }
  return bu::evaluateChebyshev(t, *cache.chebyshev_polynomial);
}

double Curve::length(double t1, double t2) const { return length(t2) - length(t1); }

double Curve::step(double t, double ds) const
{
  if (std::fabs(ds) < bu::epsilon) // no-op
    return t;

  double s_t = length(t);

  struct RootState
  {
    double t; // curve parameter
    double s; // arc length offset
  } lbracket, rbracket, guess{t, 0.0};

  if (ds < 0)
  {
    lbracket = {0.0, -s_t};
    if (ds < lbracket.s + bu::epsilon) // out-of-scope
      return 0.0;
    rbracket = guess;
  }
  else // ds > 0
  {
    rbracket = {1.0, length() - s_t};
    if (ds > rbracket.s - bu::epsilon) // out-of-scope
      return 1.0;
    lbracket = guess;
  }

  while (std::fabs(guess.s - ds) > bu::epsilon)
  {
    // Halley's method
    double f = guess.s - ds;
    Vector dC = derivativeAt(guess.t), dC2 = derivativeAt(2, guess.t);
    double df = dC.norm(), df2 = dC.dot(dC2) / df;
    guess.t -= (2 * f * df) / (2 * df * df - f * df2);

    // root bracketing, if not in bounds, use bisection method
    if (guess.t <= lbracket.t || guess.t >= rbracket.t)
      guess.t = (lbracket.t + rbracket.t) / 2;

    if (rbracket.t - lbracket.t < bu::epsilon)
      break;

    guess.s = length(guess.t) - s_t;
    (guess.s < ds ? lbracket : rbracket) = guess;
  }

  return guess.t;
}

void Curve::reverse()
{
  control_points_ = control_points_.colwise().reverse().eval();
  cache.clear();
}

void Curve::setControlPoint(unsigned idx, const Point& point)
{
  control_points_.row(idx) = point;
  cache.clear();
}

void Curve::raiseOrder()
{
  control_points_ = bc::raiseOrder(N_++) * control_points_;
  cache.clear();
}

void Curve::lowerOrder()
{
  if (N_ == 2)
    throw std::logic_error{"Cannot further reduce the order of curve."};
  control_points_ = bc::lowerOrder(N_--) * control_points_;
  cache.clear();
}

Point Curve::valueAt(double t) const
{
  if (N_ == 0)
    return {0, 0};
  return (bu::powVector(t, N_) * bc::bernstein(N_) * control_points_).transpose();
}

Eigen::MatrixX2d Curve::valueAt(const ParamVector& t_vector) const
{
  Eigen::VectorXd t_map = Eigen::Map<const Eigen::VectorXd>(t_vector.data(), t_vector.size());
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
  Vector p(derivativeAt(t));
  if (normalize && p.norm() > 0)
    p.normalize();
  return p;
}

Vector Curve::normalAt(double t, bool normalize) const
{
  Vector tangent = tangentAt(t, normalize);
  return {-tangent.y(), tangent.x()};
}

const Curve& Curve::derivative() const
{
  if (!cache.derivative)
  {
    if (N_ == 1)
      cache.derivative = std::make_unique<const Curve>(PointVector{Point(0, 0)});
    else
      cache.derivative = std::make_unique<const Curve>(
          (N_ - 1) * (control_points_.bottomRows(N_ - 1) - control_points_.topRows(N_ - 1)));
  }
  return *cache.derivative;
}

const Curve& Curve::derivative(unsigned n) const
{
  auto nth_derivative = this;
  for (unsigned k = 0; k < n; k++)
    nth_derivative = &nth_derivative->derivative();
  return *nth_derivative;
}

Vector Curve::derivativeAt(double t) const { return derivative().valueAt(t); }

Vector Curve::derivativeAt(unsigned n, double t) const { return derivative(n).valueAt(t); }

ParamVector Curve::roots() const
{
  if (!cache.roots)
  {
    cache.roots.emplace();
    if (N_ > 1)
    {
      Eigen::MatrixXd bezier_polynomial = bc::bernstein(N_) * control_points_;
      cache.roots =
          bu::concatenate(bu::solvePolynomial(bezier_polynomial.col(0)), bu::solvePolynomial(bezier_polynomial.col(1)));
    }
  }
  return *cache.roots;
}

ParamVector Curve::extrema() const { return derivative().roots(); }

BoundingBox Curve::boundingBox() const
{
  if (!cache.bounding_box)
  {
    auto extremes = valueAt(extrema());
    extremes.conservativeResize(extremes.rows() + 2, Eigen::NoChange);
    extremes.row(extremes.rows() - 1) = control_points_.row(0);
    extremes.row(extremes.rows() - 2) = control_points_.row(N_ - 1);

    cache.bounding_box.emplace(Point(extremes.col(0).minCoeff(), extremes.col(1).minCoeff()),
                               Point(extremes.col(0).maxCoeff(), extremes.col(1).maxCoeff()));
  }
  return *cache.bounding_box;
}

std::vector<Curve> Curve::splitCurve(const ParamVector& t) const
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

std::pair<Curve, Curve> Curve::splitCurve(double t) const
{
  return {Curve(bc::leftSplit(N_, t) * control_points_), Curve(bc::rightSplit(N_, t) * control_points_)};
}

PointVector Curve::intersections(const Curve& curve) const
{
  std::vector<std::pair<Eigen::MatrixX2d, Eigen::MatrixX2d>> cp_pairs;
  if (!control_points_.isApprox(curve.control_points_))
    cp_pairs.emplace_back(control_points_, curve.control_points_);
  else
  {
    // If self-similar, split curve into subcurves at extremas and compare each pair of distinct subcurves
    auto subcurves = splitCurve(extrema());
    for (unsigned k{}; k < subcurves.size(); k++)
      for (unsigned i{k + 1}; i < subcurves.size(); i++)
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

  PointVector intersections;
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
  if (!cache.projection_polynomial_const || !cache.projection_polynomial_der)
  {
    Eigen::MatrixXd curve_polynomial = (bc::bernstein(N_) * control_points_);
    Eigen::MatrixX2d derivate_polynomial = (bc::bernstein(N_ - 1) * derivative().control_points_);

    Eigen::VectorXd polynomial_part = Eigen::VectorXd::Zero(curve_polynomial.rows() + derivate_polynomial.rows() - 1);
    for (unsigned k = 0; k < curve_polynomial.rows(); k++)
      polynomial_part.middleRows(k, derivate_polynomial.rows()) +=
          derivate_polynomial * curve_polynomial.row(k).transpose();

    cache.projection_polynomial_const.emplace(std::move(polynomial_part));
    cache.projection_polynomial_der.emplace(std::move(derivate_polynomial));
  }

  Eigen::VectorXd polynomial = cache.projection_polynomial_const.value();
  polynomial.topRows(N_ - 1) -= cache.projection_polynomial_der.value() * point;

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

  // pascal triangle matrix (binomial coefficients) - rowwise
  Eigen::MatrixXd pascal_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  pascal_matrix.row(0).setOnes();
  for (unsigned k = 1; k <= c_order; k++)
    for (unsigned i = 1; i <= k; i++)
      pascal_matrix(i, k) = pascal_matrix(i - 1, k - 1) + pascal_matrix(i, k - 1);

  // inverse of pascal matrix, i.e., pascal matrix with alternating signs - colwise
  Eigen::MatrixXd pascal_alternating_matrix = pascal_matrix.transpose().inverse();

  // https://en.wikipedia.org/wiki/Bell_polynomials -> equivalent to equations of geometric continuity
  Eigen::MatrixXd bell_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  bell_matrix(0, c_order) = 1;
  for (unsigned k = 0; k < c_order; k++)
    bell_matrix.block(1, c_order - k - 1, k + 1, 1) =
        bell_matrix.block(0, c_order - k, k + 1, k + 1) *
        pascal_matrix.block(0, k, k + 1, 1)
            .cwiseProduct(Eigen::Map<const Eigen::MatrixXd>(beta_coeffs.data(), k + 1, 1));

  // diagonal: (N-1)! / (N-k-1)!
  Eigen::MatrixXd factorial_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  factorial_matrix(0, 0) = 1;
  for (unsigned k = 1; k <= c_order; k++)
    factorial_matrix(k, k) = factorial_matrix(k - 1, k - 1) * (N_ - k);

  // derivatives of given curve
  Eigen::Matrix2Xd derivatives(Eigen::Index(2), Eigen::Index(c_order + 1));
  for (unsigned k = 0; k < c_order + 1; k++)
    derivatives.col(k) = curve.derivative(k).control_points_.bottomRows(1).transpose();

  // based on the beta coefficients and geometric continuity equations, calculate new derivatives
  Eigen::MatrixXd new_derivatives = (derivatives * bell_matrix).rowwise().reverse().transpose();

  // calculate new control points
  control_points_.topRows(c_order + 1) = (factorial_matrix * pascal_alternating_matrix).inverse() * new_derivatives;
  cache.clear();
}

void Curve::Cache::clear()
{
  derivative.reset();
  roots.reset();
  bounding_box.reset();
  polyline.reset();
  polyline_t.reset();
  projection_polynomial_const.reset();
  projection_polynomial_der.reset();
  chebyshev_polynomial.reset();
  polyline_flatness = 0.0;
}

