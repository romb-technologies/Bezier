#include "Bezier/bezier.h"

#include <numeric>

#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/Polynomials>

using namespace Bezier;

///// Additional declarations

#ifndef __cpp_lib_make_unique
namespace std
{
template <typename T, typename... Args> inline std::unique_ptr<T> make_unique(Args&&... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
} // namespace std
#endif

struct _PolynomialRoots : public std::vector<double>
{
  explicit _PolynomialRoots(unsigned reserve) { std::vector<double>::reserve(reserve); }
  void clear(){};          // no-op so that PolynomialSolver::RealRoots() doesn't clear it
  void push_back(double t) // only allow valid roots
  {
    if (t >= 0 && t <= 1)
      std::vector<double>::push_back(t);
  }
};

inline unsigned _exp2(unsigned exp) { return 1 << exp; }

inline double _pow(double base, unsigned exp)
{
  double result = exp & 1 ? base : 1;
  while (exp >>= 1)
  {
    base *= base;
    if (exp & 1)
      result *= base;
  }
  return result;
}

inline Eigen::RowVectorXd _powSeries(double base, unsigned exp)
{
  Eigen::RowVectorXd power_series(exp);
  power_series(0) = 1;
  for (unsigned k = 1; k < exp; k++)
    power_series(k) = power_series(k - 1) * base;
  return power_series;
}

inline Eigen::VectorXd _trimZeroes(const Eigen::VectorXd& vec)
{
  auto idx = vec.size();
  while (idx && std::abs(vec(idx - 1)) < _epsilon)
    --idx;
  return vec.head(idx);
}

///// Curve::Curve

Curve::Curve(Eigen::MatrixX2d points) : control_points_(std::move(points)), N_(control_points_.rows()) {}

Curve::Curve(const PointVector& points)
    : control_points_(Eigen::Index(points.size()), Eigen::Index(2)), N_(points.size())
{
  for (unsigned k = 0; k < N_; k++)
    control_points_.row(k) = points[k];
}

Curve::Curve(const Curve& curve) : Curve(curve.control_points_) {}

Curve& Curve::operator=(const Curve& curve)
{
  control_points_ = curve.control_points_;
  resetCache();
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

PointVector Curve::polyline(double flatness) const
{
  if (!cached_polyline_ || cached_polyline_flatness_ != flatness)
  {
    cached_polyline_flatness_ = flatness;
    cached_polyline_ = std::make_unique<PointVector>();
    cached_polyline_->emplace_back(control_points_.row(0));

    std::vector<Eigen::MatrixX2d> subcurves;
    subcurves.emplace_back(control_points_);

    // we calculate in squared distances
    flatness *= flatness;

    double coeff{1};
    if (N_ < 10)
    {
      // for N_ == 10, coeff is 0.9922, so we ignore it for higher orders
      coeff -= std::exp2(2. - N_);
      coeff *= coeff;
    }

    while (!subcurves.empty())
    {
      Eigen::MatrixX2d cp(std::move(subcurves.back()));
      subcurves.pop_back();
      const Point& p1 = cp.row(0);
      const Point& p2 = cp.row(N_ - 1);
      Vector u = p2 - p1;

      auto deviation = [&p1, &p2, &u](double x, double y) {
        Point q(x, y);
        Vector v = q - p1;
        double t = u.dot(v) / u.squaredNorm();
        if (t < 0)
          return v.squaredNorm();
        if (t > 1)
          return (q - p2).squaredNorm();
        return (p1 + t * u - q).squaredNorm();
      };

      if (coeff * cp.rowwise().redux(deviation).maxCoeff() <= flatness)
        cached_polyline_->emplace_back(cp.row(N_ - 1));
      else
      {
        subcurves.emplace_back(splittingCoeffsRight(N_) * cp);
        subcurves.emplace_back(splittingCoeffsLeft(N_) * cp);
      }
    }
  }

  return *cached_polyline_;
}

double Curve::length() const { return length(1.0); }

double Curve::length(double t) const
{
  if (t < 0.0 || t > 1.0)
    throw std::logic_error{"Length can only be calculated for t within [0.0, 1.0] range."};

  auto evaluate_chebyshev = [](double t, const Eigen::VectorXd& coeff) {
    t = 2 * t - 1;
    double tn{t}, tn_1{1}, res{coeff(0) + coeff(1) * t};
    for (unsigned k = 2; k < coeff.size(); k++)
    {
      std::swap(tn_1, tn);
      tn = 2 * t * tn_1 - tn;
      res += coeff(k) * tn;
    }
    return res;
  };

  if (!cached_chebyshev_coeffs_)
  {
    constexpr unsigned START_LOG_N = 10;
    unsigned log_n = START_LOG_N - 1;
    unsigned n = _exp2(START_LOG_N - 1);

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
        auto lin_spaced = Eigen::ArrayXi::LinSpaced(_exp2(k - 1), 0, _exp2(k - 1) - 1);
        auto index_c = _exp2(log_n + 1 - (k + 1)) + lin_spaced * _exp2(log_n + 1 - k);
        auto index_dc = _exp2(k - 1) + 1 + lin_spaced;
        // TODO: make use of slicing & indexing in Eigen3.4
        // coeff(index_c) = coeff(N - index_c) = derivative_cache(index_dc) / n;
        for (unsigned i = 0; i < lin_spaced.size(); i++)
          coeff(index_c(i)) = coeff(N - index_c(i)) = derivative_cache(index_dc(i)) / n;
      }

      fft.fwd(fft_out, coeff);
      chebyshev = (fft_out.real().head(n - 1) - fft_out.real().segment(2, n - 1)).array() /
                  Eigen::ArrayXd::LinSpaced(n - 1, 4, 4 * (n - 1));
    } while (std::fabs(chebyshev.tail<1>()[0]) > _epsilon * 1e-2);

    unsigned cut = 0;
    while (std::fabs(chebyshev(cut)) > _epsilon * 1e-2)
      cut++;
    cached_chebyshev_coeffs_ = std::make_unique<Eigen::VectorXd>(cut + 1);
    (*cached_chebyshev_coeffs_) << 0, chebyshev.head(cut);
    (*cached_chebyshev_coeffs_)(0) = -evaluate_chebyshev(0, *cached_chebyshev_coeffs_);
  }
  return evaluate_chebyshev(t, *cached_chebyshev_coeffs_);
}

double Curve::length(double t1, double t2) const { return length(t2) - length(t1); }

double Curve::iterateByLength(double t, double s) const
{
  if (std::fabs(s) < _epsilon) // no-op
    return t;

  double s_t = length(t);

  std::pair<double, double> lbracket, rbracket, guess{t, 0.0};
  if (s < 0)
  {
    lbracket = {0.0, -s_t};
    if (s < lbracket.second + _epsilon) // out-of-scope
      return 0.0;
    rbracket = guess;
  }
  else // s > 0
  {
    rbracket = {1.0, length() - s_t};
    if (s > rbracket.second - _epsilon) // out-of-scope
      return 1.0;
    lbracket = guess;
  }

  while (std::fabs(guess.second - s) > _epsilon)
  {
    // Halley's method
    double f = guess.second - s;
    double f_d = derivativeAt(guess.first).norm();
    double f_d2 = derivativeAt(2, guess.first).norm();
    guess.first -= (2 * f * f_d) / (2 * f_d * f_d - f * f_d2);

    // root bracketing, if not in bounds, use bisection method
    if (guess.first <= lbracket.first || guess.first >= rbracket.first)
      guess.first = (lbracket.first + rbracket.first) / 2;

    if (rbracket.first - lbracket.first < _epsilon)
      break;

    guess.second = length(guess.first) - s_t;
    (guess.second < s ? lbracket : rbracket) = guess;
  }

  return guess.first;
}

void Curve::reverse()
{
  control_points_ = control_points_.colwise().reverse().eval();
  resetCache();
}

void Curve::setControlPoint(unsigned idx, const Point& point)
{
  control_points_.row(idx) = point;
  resetCache();
}

void Curve::manipulateCurvature(double t, const Point& point)
{
  if (N_ < 3 || N_ > 4)
    throw std::logic_error{"Only quadratic and cubic curves can be manipulated"};

  double r = std::fabs((_pow(t, N_ - 1) + _pow(1 - t, N_ - 1) - 1) / (_pow(t, N_ - 1) + _pow(1 - t, N_ - 1)));
  double u = _pow(1 - t, N_ - 1) / (_pow(t, N_ - 1) + _pow(1 - t, N_ - 1));
  Point C = u * control_points_.row(0) + (1 - u) * control_points_.row(N_ - 1);
  const Point& B = point;
  Point A = B - (C - B) / r;

  switch (N_)
  {
  case 3:
    control_points_.row(1) = A;
    break;
  case 4:
    Point e1 = control_points_.row(0) * _pow(1 - t, 2) + control_points_.row(1) * 2 * t * (1 - t) +
               control_points_.row(2) * _pow(t, 2);
    Point e2 = control_points_.row(1) * _pow(1 - t, 2) + control_points_.row(2) * 2 * t * (1 - t) +
               control_points_.row(3) * _pow(t, 2);
    e1 = B + e1 - valueAt(t);
    e2 = B + e2 - valueAt(t);
    Point v1 = A - (A - e1) / (1 - t);
    Point v2 = A + (e2 - A) / t;
    control_points_.row(1).noalias() = control_points_.row(0) + (v1.transpose() - control_points_.row(0)) / t;
    control_points_.row(2).noalias() = control_points_.row(3) - (control_points_.row(3) - v2.transpose()) / (1 - t);
  }
  resetCache();
}

void Curve::elevateOrder()
{
  control_points_ = elevateOrderCoeffs(N_) * control_points_;
  resetCache();
}

void Curve::lowerOrder()
{
  if (N_ == 2)
    throw std::logic_error{"Cannot further reduce the order of curve."};
  control_points_ = lowerOrderCoeffs(N_) * control_points_;
  resetCache();
}

Point Curve::valueAt(double t) const
{
  if (N_ == 0)
    return {0, 0};
  return (_powSeries(t, N_) * bernsteinCoeffs(N_) * control_points_).transpose();
}

Eigen::MatrixX2d Curve::valueAt(const std::vector<double>& t_vector) const
{
  Eigen::MatrixXd power_basis(t_vector.size(), N_);
  for (unsigned k = 0; k < t_vector.size(); k++)
    power_basis.row(k) = _powSeries(t_vector[k], N_);
  return power_basis * bernsteinCoeffs(N_) * control_points_;
}

double Curve::curvatureAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);

  return (d1.x() * d2.y() - d1.y() * d2.x()) / _pow(d1.norm(), 3);
}

double Curve::curvatureDerivativeAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);
  Vector d3 = derivativeAt(3, t);

  return (d1.x() * d3.y() - d1.y() * d3.x()) / _pow(d1.norm(), 3) -
         3 * d1.dot(d2) * (d1.x() * d2.y() - d1.y() * d2.x()) / _pow(d1.norm(), 5);
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
  if (!cached_derivative_)
  {
    cached_derivative_ = N_ == 1 ? std::make_unique<const Curve>(PointVector{Point(0, 0)})
                                 : std::make_unique<const Curve>((N_ - 1) * (control_points_.bottomRows(N_ - 1) -
                                                                             control_points_.topRows(N_ - 1)));
  }
  return *cached_derivative_;
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

std::vector<double> Curve::roots() const
{
  if (!cached_roots_)
  {
    cached_roots_ = std::make_unique<std::vector<double>>();
    if (N_ > 1)
    {
      Eigen::MatrixXd bezier_polynomial = bernsteinCoeffs(N_) * control_points_;
      Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
      auto trimmed_x = _trimZeroes(bezier_polynomial.col(0));
      auto trimmed_y = _trimZeroes(bezier_polynomial.col(1));
      _PolynomialRoots roots(trimmed_x.size() + trimmed_y.size());
      if (trimmed_x.size() > 1)
      {
        poly_solver.compute(trimmed_x);
        poly_solver.realRoots(roots);
      }
      if (trimmed_y.size() > 1)
      {
        poly_solver.compute(trimmed_y);
        poly_solver.realRoots(roots);
      }
      std::swap(roots, *cached_roots_);
    }
  }
  return *cached_roots_;
}

std::vector<double> Curve::extrema() const { return derivative().roots(); }

BoundingBox Curve::boundingBox() const
{
  if (!cached_bounding_box_)
  {
    auto extremes = valueAt(extrema());
    extremes.conservativeResize(extremes.rows() + 2, Eigen::NoChange);
    extremes.row(extremes.rows() - 1) = control_points_.row(0);
    extremes.row(extremes.rows() - 2) = control_points_.row(N_ - 1);

    cached_bounding_box_ = std::make_unique<BoundingBox>(Point(extremes.col(0).minCoeff(), extremes.col(1).minCoeff()),
                                                         Point(extremes.col(0).maxCoeff(), extremes.col(1).maxCoeff()));
  }
  return *cached_bounding_box_;
}

std::pair<Curve, Curve> Curve::splitCurve(double t) const
{
  return {Curve(splittingCoeffsLeft(N_, t) * control_points_), Curve(splittingCoeffsRight(N_, t) * control_points_)};
}

PointVector Curve::intersections(const Curve& curve) const
{
  PointVector intersections;
  auto addIntersection = [&intersections](Point new_point) {
    // check if not already found, and add new point
    if (std::none_of(intersections.begin(), intersections.end(),
                     [&new_point](const Point& point) { return (point - new_point).norm() < _epsilon; }))
      intersections.emplace_back(std::move(new_point));
  };

  std::vector<std::pair<Eigen::MatrixX2d, Eigen::MatrixX2d>> subcurve_pairs;

  if (this != &curve)
    subcurve_pairs.emplace_back(control_points_, curve.control_points_);
  else
  {
    // for self intersections divide curve into subcurves at extrema
    auto t = extrema();
    std::sort(t.begin(), t.end());
    std::vector<Eigen::MatrixX2d> subcurves;
    subcurves.emplace_back(control_points_);
    for (unsigned k = 0; k < t.size(); k++)
    {
      Eigen::MatrixX2d new_cp = std::move(subcurves.back());
      subcurves.pop_back();
      subcurves.emplace_back(splittingCoeffsLeft(N_, t[k] - _epsilon / 2) * new_cp);
      subcurves.emplace_back(splittingCoeffsRight(N_, t[k] + _epsilon / 2) * new_cp);

#if __cpp_init_captures
      std::for_each(t.begin() + k + 1, t.end(), [t = t[k]](double& x) { x = (x - t) / (1 - t); });
#else
      std::for_each(t.begin() + k + 1, t.end(), [&t, k](double& x) { x = (x - t[k]) / (1 - t[k]); });
#endif
    }

    // create all pairs of subcurves
    for (unsigned k = 0; k < subcurves.size(); k++)
      for (unsigned i = k + 1; i < subcurves.size(); i++)
        subcurve_pairs.emplace_back(subcurves[k], subcurves[i]);
  }

  while (!subcurve_pairs.empty())
  {
#if __cpp_structured_bindings
    auto [cp_a, cp_b] = std::move(subcurve_pairs.back());
#else
    Eigen::MatrixX2d cp_a, cp_b;
    std::tie(cp_a, cp_b) = std::move(subcurve_pairs.back());
#endif
    subcurve_pairs.pop_back();

    BoundingBox bbox1(Point(cp_a.col(0).minCoeff(), cp_a.col(1).minCoeff()),
                      Point(cp_a.col(0).maxCoeff(), cp_a.col(1).maxCoeff()));
    BoundingBox bbox2(Point(cp_b.col(0).minCoeff(), cp_b.col(1).minCoeff()),
                      Point(cp_b.col(0).maxCoeff(), cp_b.col(1).maxCoeff()));

    if (!bbox1.intersects(bbox2))
      ; // no intersection
    else if (bbox1.diagonal().norm() < _epsilon)
      addIntersection(bbox1.center());
    else if (bbox2.diagonal().norm() < _epsilon)
      addIntersection(bbox2.center());
    else
    {
      // intersection exists, but segments are still too large
      // - divide both segments in half
      // - insert all combinations for next iteration
      // - last pair is one where both subcurves have smallest t ranges
      Eigen::MatrixX2d subcurve_a_1(splittingCoeffsRight(N_) * cp_a);
      Eigen::MatrixX2d subcurve_a_2(splittingCoeffsLeft(N_) * cp_a);
      Eigen::MatrixX2d subcurve_b_1(splittingCoeffsRight(cp_b.rows()) * cp_b);
      Eigen::MatrixX2d subcurve_b_2(splittingCoeffsLeft(cp_b.rows()) * cp_b);
      subcurve_pairs.emplace_back(subcurve_a_1, subcurve_b_1);
      subcurve_pairs.emplace_back(subcurve_a_2, std::move(subcurve_b_1));
      subcurve_pairs.emplace_back(std::move(subcurve_a_1), subcurve_b_2);
      subcurve_pairs.emplace_back(std::move(subcurve_a_2), std::move(subcurve_b_2));
    }
  }

  return intersections;
}

double Curve::projectPoint(const Point& point) const
{
  if (!cached_projection_polynomial_part_)
  {
    Eigen::MatrixXd curve_polynomial = (bernsteinCoeffs(N_) * control_points_);
    Eigen::MatrixXd derivate_polynomial = (bernsteinCoeffs(N_ - 1) * derivative().control_points_);

    Eigen::VectorXd polynomial_part = Eigen::VectorXd::Zero(curve_polynomial.rows() + derivate_polynomial.rows() - 1);
    for (unsigned k = 0; k < curve_polynomial.rows(); k++)
      polynomial_part.middleRows(k, derivate_polynomial.rows()) +=
          derivate_polynomial * curve_polynomial.row(k).transpose();

    cached_projection_polynomial_part_ = std::make_unique<Eigen::VectorXd>(std::move(polynomial_part));
    cached_projection_polynomial_derivative_ = std::move(derivate_polynomial);
  }

  Eigen::VectorXd polynomial = *cached_projection_polynomial_part_;
  polynomial.topRows(cached_projection_polynomial_derivative_.rows()) -=
      cached_projection_polynomial_derivative_ * point;

  auto trimmed = _trimZeroes(polynomial);
  _PolynomialRoots candidates(trimmed.size());
  if (trimmed.size() > 1)
  {
    Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver(trimmed);
    poly_solver.realRoots(candidates);
  }

  double min_t = (point - valueAt(0.0)).norm() < (point - valueAt(1.0)).norm() ? 0.0 : 1.0;
  double min_dist = (point - valueAt(min_t)).norm();

  return std::accumulate(candidates.begin(), candidates.end(), std::make_pair(min_t, min_dist),
                         [this, &point](std::pair<double, double> min, double t) {
                           double dist = (point - valueAt(t)).norm();
                           return dist < min.second ? std::make_pair(t, dist) : min;
                         })
      .first;
}

double Curve::distance(const Point& point) const { return (point - valueAt(projectPoint(point))).norm(); }

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
  resetCache();
}

void Curve::resetCache()
{
  N_ = control_points_.rows();
  cached_derivative_.reset();
  cached_roots_.reset();
  cached_bounding_box_.reset();
  cached_polyline_.reset();
  cached_projection_polynomial_part_.reset();
  cached_chebyshev_coeffs_.reset();
}

Curve::CoeffsMap Curve::bernstein_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_left_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_right_ = CoeffsMap();
Curve::CoeffsMap Curve::elevate_order_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::lower_order_coeffs_ = CoeffsMap();

Curve::Coeffs Curve::bernsteinCoeffs(unsigned n)
{
  if (!bernstein_coeffs_.count(n))
  {
    bernstein_coeffs_.insert({n, Coeffs::Zero(n, n)});
    bernstein_coeffs_[n].diagonal(-1).setLinSpaced(-1, -static_cast<int>(n - 1));
    bernstein_coeffs_[n] = bernstein_coeffs_[n].exp();
    for (unsigned k = 0, binomial = 1; k < n; binomial = binomial * (n - k - 1) / (k + 1), k++)
      bernstein_coeffs_[n].row(k) *= binomial;
  }
  return bernstein_coeffs_[n];
}

Curve::Coeffs Curve::splittingCoeffsLeft(unsigned n, double t)
{
  if (t == 0.5)
  {
    if (!splitting_coeffs_left_.count(n))
    {
      splitting_coeffs_left_.insert({n, Coeffs::Zero(n, n)});
      splitting_coeffs_left_[n].diagonal() = _powSeries(0.5, n);
      splitting_coeffs_left_[n] = bernsteinCoeffs(n).inverse() * splitting_coeffs_left_[n] * bernsteinCoeffs(n);
    }
    return splitting_coeffs_left_[n];
  }

  Curve::Coeffs coeffs(Coeffs::Zero(n, n));
  coeffs.diagonal() = _powSeries(t, n);
  return bernsteinCoeffs(n).inverse() * coeffs * bernsteinCoeffs(n);
}

Curve::Coeffs Curve::splittingCoeffsRight(unsigned n, double t)
{
  if (t == 0.5)
  {
    if (!splitting_coeffs_right_.count(n))
    {
      splitting_coeffs_right_.insert({n, Coeffs::Zero(n, n)});
      Curve::Coeffs temp_splitting_coeffs_left = splittingCoeffsLeft(n);
      for (unsigned k = 0; k < n; k++)
        splitting_coeffs_right_[n].block(0, n - 1 - k, n - k, 1) =
            temp_splitting_coeffs_left.diagonal(-static_cast<int>(k)).reverse();
    }
    return splitting_coeffs_right_[n];
  }

  Curve::Coeffs coeffs(Coeffs::Zero(n, n));
  Curve::Coeffs temp_splitting_coeffs_left = splittingCoeffsLeft(n, t);
  for (unsigned k = 0; k < n; k++)
    coeffs.block(0, n - 1 - k, n - k, 1) = temp_splitting_coeffs_left.diagonal(-static_cast<int>(k)).reverse();
  return coeffs;
}

Curve::Coeffs Curve::elevateOrderCoeffs(unsigned n)
{
  if (!elevate_order_coeffs_.count(n))
  {
    elevate_order_coeffs_.insert({n, Coeffs::Zero(n + 1, n)});
    elevate_order_coeffs_[n].diagonal().setLinSpaced(1, 1 - (n - 1.) / n);
    elevate_order_coeffs_[n].diagonal(-1).setLinSpaced(1. / n, 1);
  }
  return elevate_order_coeffs_[n];
}

Curve::Coeffs Curve::lowerOrderCoeffs(unsigned n)
{
  if (!lower_order_coeffs_.count(n))
  {
    lower_order_coeffs_.insert({n, Coeffs::Zero(n - 1, n)});
    lower_order_coeffs_[n].noalias() = (elevateOrderCoeffs(n - 1).transpose() * elevateOrderCoeffs(n - 1)).inverse() *
                                       elevateOrderCoeffs(n - 1).transpose();
  }
  return lower_order_coeffs_[n];
}
