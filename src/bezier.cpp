#include "Bezier/bezier.h"
#include "Bezier/coefficients.h"
#include "Bezier/declarations.h"
#include "Bezier/utils.h"

#include <unsupported/Eigen/FFT>

#include <numeric>

using namespace Bezier;
namespace bu = Bezier::Utils;
namespace bc = Bezier::Coefficients;

Curve::Curve(Eigen::MatrixX2d points) : N_(points.rows()), control_points_(std::move(points)) {}

Curve::Curve(const PointVector& points) : N_(points.size()), control_points_(N_, 2)
{
  for (unsigned k{}; k < N_; k++)
    control_points_.row(k) = points[k];
}

Curve::Curve(const Curve& curve) : N_(curve.N_), control_points_(curve.control_points_)
{
  std::lock_guard lock{curve.cache_mutex_};
  cache_ = curve.cache_;
}

Curve::Curve(Curve&& curve) noexcept
    : N_(curve.N_), control_points_(std::move(curve.control_points_)), cache_(std::move(curve.cache_))
{
}

Curve& Curve::operator=(const Curve& curve)
{
  Curve tmp{curve}; // self-assignment safe
  swap(tmp);
  return *this;
}

Curve& Curve::operator=(Curve&& curve) noexcept
{
  swap(curve);
  return *this;
}

void Curve::swap(Curve& other) noexcept
{
  std::swap(N_, other.N_);
  control_points_.swap(other.control_points_);
  std::swap(cache_, other.cache_);
}

Curve::Cache::Cache(const Cache& other)
    : derivative(other.derivative ? std::make_unique<const Curve>(*other.derivative) : nullptr), roots(other.roots),
      bounding_box(other.bounding_box), polyline(other.polyline), polyline_t(other.polyline_t),
      polyline_flatness(other.polyline_flatness), projection_polynomial_const(other.projection_polynomial_const),
      projection_polynomial_der(other.projection_polynomial_der), chebyshev_polynomial(other.chebyshev_polynomial)
{
}

Curve::Cache& Curve::Cache::operator=(const Cache& other)
{
  Cache tmp{other};
  *this = std::move(tmp);
  return *this;
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

unsigned Curve::order() const { return N_ - 1; }

PointVector Curve::controlPoints() const
{
  PointVector points(N_);
  for (unsigned k{}; k < N_; k++)
    points[k] = control_points_.row(k);
  return points;
}

Point Curve::controlPoint(unsigned idx) const { return control_points_.row(idx); }

std::pair<Point, Point> Curve::endPoints() const { return {control_points_.row(0), control_points_.row(N_ - 1)}; }

PointVector Curve::polyline() const { return polyline(boundingBox().diagonal().norm() / 1000); }

PointVector Curve::polyline(double flatness) const
{
  std::lock_guard lock{cache_mutex_};
  if (cache_.polyline && std::fabs(cache_.polyline_flatness - flatness) < bu::epsilon)
    return *cache_.polyline;

  cache_.polyline_flatness = flatness;
  cache_.polyline_t.emplace({0.0});
  cache_.polyline.emplace({Point(control_points_.row(0))});

  std::vector<std::tuple<Eigen::MatrixX2d, double, double>> subcurves;
  subcurves.emplace_back(control_points_, 0.0, 1.0);

  while (!subcurves.empty())
  {
    auto [cp, t1, t2] = std::move(subcurves.back());
    subcurves.pop_back();

    if (bu::maxDeviation(cp) <= flatness)
    {
      cache_.polyline_t->emplace_back(t2);
      cache_.polyline->emplace_back(cp.row(N_ - 1));
    }
    else
    {
      subcurves.emplace_back(bc::rightSplit(N_) * cp, (t1 + t2) / 2, t2);
      subcurves.emplace_back(bc::leftSplit(N_) * cp, t1, (t1 + t2) / 2);
    }
  }

  return *cache_.polyline;
}

ParamVector Curve::polylineParams() const { return polylineParams(boundingBox().diagonal().norm() / 1000); }

ParamVector Curve::polylineParams(double flatness) const
{
  std::lock_guard lock{cache_mutex_};
  polyline(flatness);
  return *cache_.polyline_t;
}

double Curve::length() const { return length(1.0); }

double Curve::length(double t) const
{
  if (t < 0.0 || t > 1.0)
    throw std::logic_error{"Length can only be calculated for t within [0.0, 1.0] range."};

  std::lock_guard lock{cache_mutex_};
  if (cache_.chebyshev_polynomial)
    return bu::evaluateChebyshev(t, *cache_.chebyshev_polynomial);
  auto& chebyshev = cache_.chebyshev_polynomial.emplace();

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

  // Trim trailing zero coefficients (keep at least 2)
  unsigned idx = chebyshev.size();
  while (idx > 2 && std::fabs(chebyshev(idx - 1)) < bu::epsilon)
    --idx;
  chebyshev.conservativeResize(idx);

  chebyshev(0) -= bu::evaluateChebyshev(0, chebyshev);
  return bu::evaluateChebyshev(t, chebyshev);
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
    if (!(guess.t > lbracket.t && guess.t < rbracket.t))
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
  cache_.clear();
}

void Curve::setControlPoint(unsigned idx, const Point& point)
{
  control_points_.row(idx) = point;
  cache_.clear();
}

void Curve::raiseOrder()
{
  control_points_ = bc::raiseOrder(N_++) * control_points_;
  cache_.clear();
}

void Curve::lowerOrder()
{
  if (N_ == 2)
    throw std::logic_error{"Cannot further reduce the order of curve."};
  control_points_ = bc::lowerOrder(N_--) * control_points_;
  cache_.clear();
}

Point Curve::valueAt(double t) const
{
  return N_ ? (bu::powVector(t, N_) * bc::bernstein(N_) * control_points_).transpose() : Point(0, 0);
}

Eigen::MatrixX2d Curve::valueAt(const ParamVector& t_vector) const
{
  auto t_map = Eigen::Map<const Eigen::VectorXd>(t_vector.data(), t_vector.size());
  return bu::powMatrix(t_map, N_) * bc::bernstein(N_) * control_points_;
}

double Curve::curvatureAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);
  if (double d1n = d1.norm(); d1n >= bu::epsilon)
    return bu::cross(d1, d2) / bu::pow(d1n, 3); // regular point (a straight curve gives 0)

  // Cusp: curvature is unbounded -> signed infinity (radius of curvature -> 0)
  // Stationary: curvature is zero (radius of curvature -> inf)
  for (unsigned a{2}; a <= order(); a++)
    if (Vector da = derivativeAt(a, t); da.squaredNorm() > bu::epsilon)
    {
      double turn = bu::cross(da, derivativeAt(a + 1, t));
      return std::fabs(turn) > bu::epsilon ? std::copysign(std::numeric_limits<double>::infinity(), turn) : 0.0;
    }
  return 0.0;
}

double Curve::curvatureDerivativeAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);
  Vector d3 = derivativeAt(3, t);
  double d1n = d1.norm();
  return d1n < bu::epsilon
             ? 0.0
             : (d1.squaredNorm() * bu::cross(d1, d3) - 3 * d1.dot(d2) * bu::cross(d1, d2)) / bu::pow(d1n, 5);
}

Vector Curve::tangentAt(double t) const
{
  // tangent direction is the first non-vanishing derivative (zero only if the curve is a point)
  for (unsigned k{1}; k <= order(); k++)
    if (Vector d = derivativeAt(k, t); d.squaredNorm() > bu::epsilon)
      return d.normalized();
  return Vector(0, 0);
}

Vector Curve::normalAt(double t) const
{
  Vector tangent = tangentAt(t);
  return {-tangent.y(), tangent.x()};
}

const Curve& Curve::derivative() const
{
  std::lock_guard lock{cache_mutex_};
  if (!cache_.derivative)
    cache_.derivative = N_ == 1 ? std::make_unique<const Curve>(PointVector{Point(0, 0)})
                                : std::make_unique<const Curve>((N_ - 1) * (control_points_.bottomRows(N_ - 1) -
                                                                            control_points_.topRows(N_ - 1)));
  return *cache_.derivative;
}

const Curve& Curve::derivative(unsigned n) const
{
  auto nth_derivative = this;
  for (unsigned k{}; k < n; k++)
    nth_derivative = &nth_derivative->derivative();
  return *nth_derivative;
}

Vector Curve::derivativeAt(double t) const { return derivative().valueAt(t); }

Vector Curve::derivativeAt(unsigned n, double t) const { return derivative(n).valueAt(t); }

ParamVector Curve::roots() const
{
  std::lock_guard lock{cache_mutex_};
  if (cache_.roots)
    return *cache_.roots;

  Eigen::MatrixX2d bezier_polynomial = bc::bernstein(N_) * control_points_;
  return cache_.roots.emplace(bu::concatenate(bu::solvePolynomial(bezier_polynomial.col(0)), //
                                              bu::solvePolynomial(bezier_polynomial.col(1))));
}

ParamVector Curve::extrema() const { return derivative().roots(); }

BoundingBox Curve::boundingBox() const
{
  std::lock_guard lock{cache_mutex_};
  if (cache_.bounding_box)
    return *cache_.bounding_box;

  auto extremes = valueAt(extrema());
  extremes.conservativeResize(extremes.rows() + 2, Eigen::NoChange);
  extremes.bottomRows<2>() << control_points_.row(0), control_points_.row(N_ - 1);
  return cache_.bounding_box.emplace(extremes.colwise().minCoeff(), extremes.colwise().maxCoeff());
}

std::vector<Curve> Curve::splitCurve(const ParamVector& t_vector) const
{
  auto t_sorted = t_vector;
  std::sort(t_sorted.begin(), t_sorted.end());
  std::vector<Curve> subcurves;
  subcurves.reserve(t_sorted.size() + 1);
  auto leftover_cp = control_points_;
  for (unsigned k{}; k < t_sorted.size(); k++)
  {
    subcurves.emplace_back(bc::leftSplit(N_, t_sorted[k]) * leftover_cp);
    leftover_cp = bc::rightSplit(N_, t_sorted[k]) * leftover_cp;
    std::for_each(t_sorted.begin() + k + 1, t_sorted.end(), [t = t_sorted[k]](double& x) { x = (x - t) / (1 - t); });
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
  if (N_ != curve.N_ || !control_points_.isApprox(curve.control_points_))
    cp_pairs.emplace_back(control_points_, curve.control_points_);
  else
  {
    // If self-similar, split curve into subcurves at extremas and compare each pair of non-adjacent subcurves
    auto t = extrema();
    std::sort(t.begin(), t.end());
    t.erase(std::unique(t.begin(), t.end(), [](double a, double b) { return std::abs(a - b) < bu::epsilon; }),
            t.end()); // collapse repeated extrema (cusps)
    auto subcurves = splitCurve(t);
    for (unsigned k{}; k < subcurves.size(); k++)
      for (unsigned i{k + 2}; i < subcurves.size(); i++)
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

    // If intersection is detected and unique, insert it
    if (oa * ob <= 0 && oc * od <= 0 && ob != oa)
    {
      Point p = (a1 * ob - a2 * oa) / (ob - oa);
      if (std::none_of(intersections.begin(), intersections.end(),
                       [&p](const Point& q) { return (p - q).norm() < bu::epsilon; }))
        intersections.emplace_back(p);
    }
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
  std::lock_guard lock{cache_mutex_};
  if (!cache_.projection_polynomial_const || !cache_.projection_polynomial_der)
  {
    Eigen::MatrixX2d curve_polynomial = bc::bernstein(N_) * control_points_;
    Eigen::MatrixX2d derivate_polynomial = bc::bernstein(N_ - 1) * derivative().control_points_;

    Eigen::VectorXd polynomial_part = Eigen::VectorXd::Zero(2 * N_ - 2);
    for (unsigned k{}; k < N_; k++)
      polynomial_part.middleRows(k, N_ - 1) += derivate_polynomial * curve_polynomial.row(k).transpose();

    cache_.projection_polynomial_const.emplace(std::move(polynomial_part));
    cache_.projection_polynomial_der.emplace(std::move(derivate_polynomial));
  }

  Eigen::VectorXd polynomial = *cache_.projection_polynomial_const;
  polynomial.topRows(N_ - 1) -= *cache_.projection_polynomial_der * point;

  double min_t{}, min_dist{std::numeric_limits<double>::infinity()};
  for (double t : bu::concatenate(bu::solvePolynomial(polynomial), {0.0, 1.0}))
    if (double dist = bu::dist(point, valueAt(t)); dist < min_dist)
      std::tie(min_t, min_dist) = std::make_pair(t, dist);
  return min_t;
}

double Curve::distance(const Point& point) const { return bu::dist(point, valueAt(projectPoint(point))); }

void Curve::applyContinuity(const Curve& curve, const std::vector<double>& beta_coeffs)
{
  unsigned c_order = beta_coeffs.size();

  // Raise this curve's order if needed so it has enough control points to
  // satisfy the requested continuity order (shape is preserved).
  while (N_ <= c_order)
    raiseOrder();
  Eigen::Map<const Eigen::VectorXd> beta(beta_coeffs.data(), beta_coeffs.size());

  // pascal triangle matrix (binomial coefficients) - rowwise
  Eigen::MatrixXd pascal_matrix = Eigen::MatrixXd::Zero(c_order + 1, c_order + 1);
  for (unsigned k = 0; k <= c_order; k++)
  {
    double c = 1.0; // C(k,0)
    for (unsigned i = 0; i <= k; i++)
    {
      pascal_matrix(i, k) = c;
      c = c * (k - i) / (i + 1); // C(k,i+1)
    }
  }

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
  Eigen::VectorXd perm = Eigen::VectorXd::NullaryExpr(c_order + 1, permFunc);

  // calculate new control points
  control_points_.topRows(c_order + 1) = pascal_matrix.transpose() * perm.cwiseInverse().asDiagonal() * new_derivatives;
  cache_.clear();
}

Curve Curve::offsetCurve(const Curve& curve, double offset, unsigned order)
{
  PointVector offset_polyline = curve.polyline();
  ParamVector polyline_t = curve.polylineParams();
  for (unsigned k{}; k < offset_polyline.size(); k++)
    offset_polyline[k] += offset * curve.normalAt(polyline_t[k]);
  return fromPolyline(offset_polyline, order);
}

Curve Curve::joinCurves(const Curve& curve1, const Curve& curve2, unsigned order)
{
  if (order == 1)
    return Curve(PointVector{curve1.control_points_.row(0), curve2.control_points_.row(curve2.N_ - 1)});
  return fromPolyline(bu::concatenate(curve1.polyline(), curve2.polyline()), order);
}

Curve Curve::fromPolyline(const PointVector& polyline, unsigned order)
{
  if (polyline.size() < 2)
    throw std::logic_error{"Polyline must have at least two points."};
  if (order == 1 || polyline.size() == 2)
    return Curve(PointVector{polyline.front(), polyline.back()});

  const unsigned M = polyline.size();

  // Lazily compute Visvalingam-Whyatt order on first call; return K most-significant points in original order.
  std::vector<unsigned> vw;
  auto reducedPolyline = [&](unsigned K) -> PointVector {
    if (K >= M)
      return polyline;
    if (vw.empty())
      vw = bu::visvalingamWhyatt(polyline);
    std::vector<unsigned> idx(vw.begin(), vw.begin() + K);
    std::sort(idx.begin(), idx.end());
    PointVector out;
    out.reserve(K);
    for (unsigned i : idx)
      out.push_back(polyline[i]);
    return out;
  };

  // Fixed order: fit its 3N most-significant points.
  if (order != 0)
    return bu::fitBezier(reducedPolyline(3 * (order + 1)), order);

  // Automatic order selection: for each candidate order fit its 3N most-significant points (Visvalingam-Whyatt)
  // and score a floored BIC (M*ln(RSS/M) + 2*(n-1)*ln(M)) over the full polyline; pick the minimum (cap 12).
  // RSS floor: M*(0.05% of bbox diagonal)^2 so numerically-perfect data does not chase zero.
  constexpr unsigned MAX_AUTO_ORDER = 12;
  const double rss_floor = [&] {
    Point lo = polyline.front(), hi = polyline.front();
    for (const Point& p : polyline)
    {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }
    return M * bu::pow(5e-4 * (hi - lo).norm(), 2);
  }();

  std::optional<Curve> best;
  double best_bic = std::numeric_limits<double>::max();
  for (unsigned n{1}; n <= MAX_AUTO_ORDER && n < M; n++)
  {
    Curve curve = bu::fitBezier(reducedPolyline(3 * (n + 1)), n);
    double rss{};
    for (const Point& p : polyline)
      rss += bu::pow(curve.distance(p), 2);
    if (rss <= rss_floor)
      return curve; // numerically perfect; higher orders would only add complexity
    double bic = M * std::log(rss / M) + 2.0 * (n - 1) * std::log(double(M));
    if (bic < best_bic)
    {
      best_bic = bic;
      best = std::move(curve);
    }
  }
  return *best;
}
