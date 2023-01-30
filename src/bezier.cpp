#include "Bezier/bezier.h"
#include "Bezier/legendre_gauss.h"

#include <numeric>

#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/Polynomials>

using namespace Bezier;

Eigen::VectorXd trimZeroes(const Eigen::VectorXd& vec)
{
  auto idx = vec.size();
  while (idx && vec(idx - 1) == 0.0)
    --idx;
  return vec.head(idx);
}

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

double Curve::length() const { return length(0.0, 1.0); }

double Curve::length(double t) const { return length(0.0, t); }

double Curve::length(double t1, double t2) const
{
#if __cpp_lib_clamp
  auto N = std::clamp(N_ * static_cast<unsigned>(std::ceil(std::fabs(t2 - t1) / 0.2)), 0u, 63u);
#else
  auto N = N_ * static_cast<unsigned>(std::ceil(std::fabs(t2 - t1) / 0.2));
  if (N > 63)
    N = 63;
#endif

  return std::accumulate(LegendreGauss::coefficients[N].begin(), LegendreGauss::coefficients[N].end(), 0.0,
                         [&](double sum, const std::pair<double, double>& coeff) {
                           return sum + std::get<1>(coeff) *
                                            derivativeAt(std::get<0>(coeff) * (t2 - t1) / 2 + (t1 + t2) / 2).norm();
                         }) *
         (t2 - t1) / 2;
}

double Curve::iterateByLength(double t, double s) const
{
  struct DistanceFunctor
  {
    const Curve& c_;
    double t_, s_;

    DistanceFunctor(const Curve& c, double t, double s) : c_(c), t_(t), s_(s) {}
    int values() const { return 1; }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
    {
      fvec(0) = c_.length(t_, x(0)) - s_;
      return 0;
    }

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fvec) const
    {
      fvec(0) = c_.derivativeAt(x(0)).norm();
      return 0;
    }

  } df(*this, t, s);

  // do not use "Eigen::HybridNonLinearSolver", because it gives wrong result near curve inflections
  Eigen::LevenbergMarquardt<DistanceFunctor> lm(df);
  Eigen::VectorXd x(1);
  x << t;
  lm.minimize(x);

  return x(0);
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

  double r =
      std::fabs((std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1) - 1) / (std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1)));
  double u = std::pow(1 - t, N_ - 1) / (std::pow(t, N_ - 1) + std::pow(1 - t, N_ - 1));
  Point C = u * control_points_.row(0) + (1 - u) * control_points_.row(N_ - 1);
  const Point& B = point;
  Point A = B - (C - B) / r;

  switch (N_)
  {
  case 3:
    control_points_.row(1) = A;
    break;
  case 4:
    Point e1 = control_points_.row(0) * std::pow(1 - t, 2) + control_points_.row(1) * 2 * t * (1 - t) +
               control_points_.row(2) * std::pow(t, 2);
    Point e2 = control_points_.row(1) * std::pow(1 - t, 2) + control_points_.row(2) * 2 * t * (1 - t) +
               control_points_.row(3) * std::pow(t, 2);
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
  return (Eigen::pow(t, Eigen::ArrayXd::LinSpaced(N_, 0, N_ - 1)).matrix().transpose() * bernsteinCoeffs(N_) *
          control_points_)
      .transpose();
}

Eigen::MatrixX2d Curve::valueAt(const std::vector<double>& t_vector) const
{
  auto t_matrix = Eigen::Map<const Eigen::VectorXd>(t_vector.data(), t_vector.size()).replicate(1, N_);
  auto p_matrix = Eigen::ArrayXd::LinSpaced(N_, 0, N_ - 1).transpose().replicate(t_vector.size(), 1);
  Eigen::MatrixXd power_basis = Eigen::pow(t_matrix.array(), p_matrix.array());
  return power_basis * bernsteinCoeffs(N_) * control_points_;
}

double Curve::curvatureAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);

  return (d1.x() * d2.y() - d1.y() * d2.x()) / std::pow(d1.norm(), 3);
}

double Curve::curvatureDerivativeAt(double t) const
{
  Vector d1 = derivativeAt(t);
  Vector d2 = derivativeAt(2, t);
  Vector d3 = derivativeAt(3, t);

  return (d1.x() * d3.y() - d1.y() * d3.x()) / std::pow(d1.norm(), 3) -
         3 * d1.dot(d2) * (d1.x() * d2.y() - d1.y() * d2.x()) / std::pow(d1.norm(), 5);
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
      std::vector<double> roots_X, roots_Y;
      Eigen::MatrixXd bezier_polynomial = bernsteinCoeffs(N_) * control_points_;
      Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver;
      auto trimmed = trimZeroes(bezier_polynomial.col(0));
      if (trimmed.size() > 1)
      {
        poly_solver.compute(trimmed);
        poly_solver.realRoots(roots_X);
      }
      trimmed = trimZeroes(bezier_polynomial.col(1));
      if (trimmed.size() > 1)
      {
        poly_solver.compute(trimmed);
        poly_solver.realRoots(roots_Y);
      }
      cached_roots_->reserve(roots_X.size() + roots_Y.size());
      std::copy_if(std::make_move_iterator(roots_X.begin()), std::make_move_iterator(roots_X.end()),
                   std::back_inserter(*cached_roots_), [](double t) { return t >= 0 && t <= 1; });
      std::copy_if(std::make_move_iterator(roots_Y.begin()), std::make_move_iterator(roots_Y.end()),
                   std::back_inserter(*cached_roots_), [](double t) { return t >= 0 && t <= 1; });
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

std::pair<Curve, Curve> Curve::splitCurve(double z) const
{
  return {Curve(splittingCoeffsLeft(N_, z) * control_points_), Curve(splittingCoeffsRight(N_, z) * control_points_)};
}

PointVector Curve::intersections(const Curve& curve) const
{
  const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

  PointVector intersections;
  auto addIntersection = [&intersections, epsilon](Point new_point) {
    // check if not already found, and add new point
    if (std::none_of(intersections.begin(), intersections.end(),
                     [&new_point, epsilon](const Point& point) { return (point - new_point).norm() < epsilon; }))
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
      subcurves.emplace_back(splittingCoeffsLeft(N_, t[k] - epsilon / 2) * new_cp);
      subcurves.emplace_back(splittingCoeffsRight(N_, t[k] + epsilon / 2) * new_cp);

#if __cpp_init_captures
      std::for_each(t.begin() + k + 1, t.end(), [t = t[k]](double& x) { x = (x - t) / (1 - t); });
#else
      std::for_each(t.begin() + k + 1, t.end(), [t, k](double& x) { x = (x - t[k]) / (1 - t[k]); });
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
    else if (bbox1.diagonal().norm() < epsilon)
      addIntersection(bbox1.center());
    else if (bbox2.diagonal().norm() < epsilon)
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

  std::vector<double> candidates;
  auto trimmed = trimZeroes(polynomial);
  if (trimmed.size() > 1)
  {
    Eigen::PolynomialSolver<double, Eigen::Dynamic> poly_solver(trimmed);
    poly_solver.realRoots(candidates);
  }

  double min_t = (point - valueAt(0.0)).norm() < (point - valueAt(1.0)).norm() ? 0.0 : 1.0;
  double min_dist = (point - valueAt(min_t)).norm();

  return std::accumulate(candidates.begin(), candidates.end(), std::make_pair(min_t, min_dist),
                         [&](std::pair<double, double> min, double t) {
                           if (t < 0 || t > 1)
                             return min;
                           double dist = (point - valueAt(t)).norm();
                           return dist < min.second ? std::make_pair(t, dist) : min;
                         })
      .first;
}

double Curve::distance(const Point& point) const { return (point - valueAt(projectPoint(point))).norm(); }

void Curve::applyContinuity(const Curve& source_curve, const std::vector<double>& beta_coeffs)
{
  unsigned c_order = beta_coeffs.size();

  Eigen::MatrixXd pascal_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  Eigen::MatrixXd pascal_alterating_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  pascal_alterating_matrix.diagonal(-1).setLinSpaced(-1, -static_cast<int>(c_order));
  pascal_alterating_matrix = pascal_alterating_matrix.exp();
  pascal_matrix = pascal_alterating_matrix.cwiseAbs().transpose();

  Eigen::MatrixXd bell_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  bell_matrix(0, c_order) = 1;

  for (unsigned k = 0; k < c_order; k++)
    bell_matrix.block(1, c_order - k - 1, k + 1, 1) =
        bell_matrix.block(0, c_order - k, k + 1, k + 1) *
        pascal_matrix.block(0, k, k + 1, 1)
            .cwiseProduct(Eigen::Map<const Eigen::MatrixXd>(beta_coeffs.data(), k + 1, 1));

  Eigen::MatrixXd factorial_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));

  // factorial(k) = std::tgamma(k+1)
  factorial_matrix.diagonal() = Eigen::ArrayXd::LinSpaced(c_order + 1, 0, c_order).unaryExpr([this](unsigned k) {
    return std::tgamma(N_) / std::tgamma(N_ - k);
  });

  Eigen::Matrix2Xd derivatives(Eigen::Index(2), Eigen::Index(c_order + 1));
  for (unsigned k = 0; k < c_order + 1; k++)
    derivatives.col(k) = source_curve.derivative(k).control_points_.bottomRows(1).transpose();

  Eigen::MatrixXd derivatives_wanted = (derivatives * bell_matrix).rowwise().reverse().transpose();

  control_points_.topRows(c_order + 1) = (factorial_matrix * pascal_alterating_matrix).inverse() * derivatives_wanted;
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
}

Curve::CoeffsMap Curve::bernstein_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_left_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_right_ = CoeffsMap();
Curve::CoeffsMap Curve::elevate_order_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::lower_order_coeffs_ = CoeffsMap();

Curve::Coeffs Curve::bernsteinCoeffs(unsigned n)
{
  if (bernstein_coeffs_.find(n) == bernstein_coeffs_.end())
  {
    auto binomial = [](unsigned n, unsigned k) {
      return std::exp(std::lgamma(n + 1) - std::lgamma(n - k + 1) - std::lgamma(k + 1));
    };
    bernstein_coeffs_.insert({n, Coeffs::Zero(n, n)});
    bernstein_coeffs_[n].diagonal(-1).setLinSpaced(-1, -static_cast<int>(n - 1));
    bernstein_coeffs_[n] = bernstein_coeffs_[n].exp();
    for (unsigned k = 0; k < n; k++)
      bernstein_coeffs_[n].row(k) *= binomial(n - 1, k);
  }
  return bernstein_coeffs_[n];
}

Curve::Coeffs Curve::splittingCoeffsLeft(unsigned n, double z)
{
  if (z == 0.5)
  {
    if (splitting_coeffs_left_.find(n) == splitting_coeffs_left_.end())
    {
      splitting_coeffs_left_.insert({n, Coeffs::Zero(n, n)});
      splitting_coeffs_left_[n].diagonal() = Eigen::pow(0.5, Eigen::ArrayXd::LinSpaced(n, 0, n - 1));
      splitting_coeffs_left_[n] = bernsteinCoeffs(n).inverse() * splitting_coeffs_left_[n] * bernsteinCoeffs(n);
    }
    return splitting_coeffs_left_[n];
  }

  Curve::Coeffs coeffs(Coeffs::Zero(n, n));
  coeffs.diagonal() = Eigen::pow(z, Eigen::ArrayXd::LinSpaced(n, 0, n - 1));
  return bernsteinCoeffs(n).inverse() * coeffs * bernsteinCoeffs(n);
}

Curve::Coeffs Curve::splittingCoeffsRight(unsigned n, double z)
{
  if (z == 0.5)
  {
    if (splitting_coeffs_right_.find(n) == splitting_coeffs_right_.end())
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
  Curve::Coeffs temp_splitting_coeffs_left = splittingCoeffsLeft(n, z);
  for (unsigned k = 0; k < n; k++)
    coeffs.block(0, n - 1 - k, n - k, 1) = temp_splitting_coeffs_left.diagonal(-static_cast<int>(k)).reverse();
  return coeffs;
}

Curve::Coeffs Curve::elevateOrderCoeffs(unsigned n)
{
  if (elevate_order_coeffs_.find(n) == elevate_order_coeffs_.end())
  {
    elevate_order_coeffs_.insert({n, Coeffs::Zero(n + 1, n)});
    elevate_order_coeffs_[n].diagonal().setLinSpaced(1, 1 - (n - 1.) / n);
    elevate_order_coeffs_[n].diagonal(-1).setLinSpaced(1. / n, 1);
  }
  return elevate_order_coeffs_[n];
}

Curve::Coeffs Curve::lowerOrderCoeffs(unsigned n)
{
  if (lower_order_coeffs_.find(n) == lower_order_coeffs_.end())
  {
    lower_order_coeffs_.insert({n, Coeffs::Zero(n - 1, n)});
    lower_order_coeffs_[n].noalias() = (elevateOrderCoeffs(n - 1).transpose() * elevateOrderCoeffs(n - 1)).inverse() *
                                       elevateOrderCoeffs(n - 1).transpose();
  }
  return lower_order_coeffs_[n];
}
