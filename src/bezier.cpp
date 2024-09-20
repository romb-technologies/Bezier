#include "Bezier/bezier.h"
#include "Bezier/utils.h"

#include <numeric>

#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/Polynomials>

using namespace Bezier;
using namespace Bezier::Utils;

Curve::Curve(Eigen::MatrixX2d points) : control_points_(std::move(points)), N_(control_points_.rows()) {}

Curve::Curve(const PointVector& points)
    : control_points_(Eigen::Index(points.size()), Eigen::Index(2)), N_(points.size())
{
  for (unsigned k{}; k < N_; k++)
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
  for (unsigned k{}; k < N_; k++)
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
    cached_polyline_t_ = std::make_unique<std::vector<double>>();
    cached_polyline_t_->emplace_back(0.0);

    std::vector<std::tuple<Eigen::MatrixX2d, double, double>> subcurves;
    subcurves.emplace_back(control_points_, 0.0, 1.0);

    // we calculate in squared distances
    flatness *= flatness;

    // for N_ == 10, coeff is 0.9922, so we ignore it for higher orders
    const double coeff{N_ >= 10 ? 1 : _pow(1 - std::exp2(2. - N_), 2)};

    while (!subcurves.empty())
    {
      auto [cp, t1, t2] = std::move(subcurves.back());
      subcurves.pop_back();
      const Point& p1 = cp.row(0);
      const Point& p2 = cp.row(N_ - 1);
      const Vector u = p2 - p1;

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
      {
        cached_polyline_->emplace_back(cp.row(N_ - 1));
        cached_polyline_t_->emplace_back(t2);
      }
      else
      {
        subcurves.emplace_back(splittingCoeffsRight(N_) * cp, (t1 + t2) / 2, t2);
        subcurves.emplace_back(splittingCoeffsLeft(N_) * cp, t1, (t1 + t2) / 2);
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
    for (unsigned k{2}; k < coeff.size(); k++)
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
    for (unsigned k{2}; k <= n; k *= 2)
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

      for (unsigned k{1}; k <= log_n; k++)
      {
        auto lin_spaced = Eigen::ArrayXi::LinSpaced(_exp2(k - 1), 0, _exp2(k - 1) - 1);
        auto index_c = _exp2(log_n + 1 - (k + 1)) + lin_spaced * _exp2(log_n + 1 - k);
        auto index_dc = _exp2(k - 1) + 1 + lin_spaced;
        // TODO: make use of slicing & indexing in Eigen3.4
        // coeff(index_c) = coeff(N - index_c) = derivative_cache(index_dc) / n;
        for (unsigned i{}; i < lin_spaced.size(); i++)
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
  for (unsigned k{}; k < t_vector.size(); k++)
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
    cached_derivative_ = N_ == 1 ? std::make_unique<const Curve>(PointVector{Point(0, 0)})
                                 : std::make_unique<const Curve>((N_ - 1) * (control_points_.bottomRows(N_ - 1) -
                                                                             control_points_.topRows(N_ - 1)));
  return *cached_derivative_;
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
    for (unsigned k{}; k < t.size(); k++)
    {
      Eigen::MatrixX2d new_cp = std::move(subcurves.back());
      subcurves.pop_back();
      subcurves.emplace_back(splittingCoeffsLeft(N_, t[k] - _epsilon / 2) * new_cp);
      subcurves.emplace_back(splittingCoeffsRight(N_, t[k] + _epsilon / 2) * new_cp);

      std::for_each(t.begin() + k + 1, t.end(), [t = t[k]](double& x) { x = (x - t) / (1 - t); });
    }

    // create all pairs of subcurves
    for (unsigned k{}; k < subcurves.size(); k++)
      for (unsigned i{k + 1}; i < subcurves.size(); i++)
        subcurve_pairs.emplace_back(subcurves[k], subcurves[i]);
  }

  while (!subcurve_pairs.empty())
  {
    auto [cp_a, cp_b] = std::move(subcurve_pairs.back());
    subcurve_pairs.pop_back();

    BoundingBox bbox1(Point(cp_a.col(0).minCoeff(), cp_a.col(1).minCoeff()),
                      Point(cp_a.col(0).maxCoeff(), cp_a.col(1).maxCoeff()));
    BoundingBox bbox2(Point(cp_b.col(0).minCoeff(), cp_b.col(1).minCoeff()),
                      Point(cp_b.col(0).maxCoeff(), cp_b.col(1).maxCoeff()));

    if (!bbox1.intersects(bbox2))
      continue; // no intersection
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
    for (unsigned k{}; k < curve_polynomial.rows(); k++)
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
    Eigen::PolynomialSolver<double, Eigen::Dynamic>(trimmed).realRoots(candidates);

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
  for (unsigned k{1}; k <= c_order; k++)
    for (unsigned i{1}; i <= k; i++)
      pascal_matrix(i, k) = pascal_matrix(i - 1, k - 1) + pascal_matrix(i, k - 1);

  // inverse of pascal matrix, i.e., pascal matrix with alternating signs - colwise
  Eigen::MatrixXd pascal_alternating_matrix = pascal_matrix.transpose().inverse();

  // https://en.wikipedia.org/wiki/Bell_polynomials -> equivalent to equations of geometric continuity
  Eigen::MatrixXd bell_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  bell_matrix(0, c_order) = 1;
  for (unsigned k{}; k < c_order; k++)
    bell_matrix.block(1, c_order - k - 1, k + 1, 1) =
        bell_matrix.block(0, c_order - k, k + 1, k + 1) *
        pascal_matrix.block(0, k, k + 1, 1)
            .cwiseProduct(Eigen::Map<const Eigen::MatrixXd>(beta_coeffs.data(), k + 1, 1));

  // diagonal: (N-1)! / (N-k-1)!
  Eigen::MatrixXd factorial_matrix(Eigen::MatrixXd::Zero(c_order + 1, c_order + 1));
  factorial_matrix(0, 0) = 1;
  for (unsigned k{1}; k <= c_order; k++)
    factorial_matrix(k, k) = factorial_matrix(k - 1, k - 1) * (N_ - k);

  // derivatives of given curve
  Eigen::Matrix2Xd derivatives(Eigen::Index(2), Eigen::Index(c_order + 1));
  for (unsigned k{}; k < c_order + 1; k++)
    derivatives.col(k) = curve.derivative(k).control_points_.bottomRows(1).transpose();

  // based on the beta coefficients and geometric continuity equations, calculate new derivatives
  Eigen::MatrixXd new_derivatives = (derivatives * bell_matrix).rowwise().reverse().transpose();

  // calculate new control points
  control_points_.topRows(c_order + 1) = (factorial_matrix * pascal_alternating_matrix).inverse() * new_derivatives;
  resetCache();
}

Curve Curve::offsetCurve(const Curve& curve, double offset, unsigned order)
{
  PointVector offset_polyline;
  if (!curve.cached_polyline_)
    curve.polyline();
  offset_polyline.reserve(curve.cached_polyline_->size());
  for (unsigned k{}; k < curve.cached_polyline_->size(); k++)
    offset_polyline.emplace_back((*curve.cached_polyline_)[k] +
                                 offset * curve.normalAt((*curve.cached_polyline_t_)[k]));
  return fromPolyline(offset_polyline, order ? order : curve.order() + 1);
}

Curve Curve::joinCurves(const Curve& curve1, const Curve& curve2, unsigned int order)
{
  if (order == 1)
    return Curve(PointVector{curve1.control_points_.row(0), curve2.control_points_.row(curve2.N_ - 1)});

  auto polyline = curve1.polyline();
  auto polyline2 = curve2.polyline();
  polyline.reserve(polyline.size() + polyline2.size());
  polyline.insert(polyline.end(), std::make_move_iterator(polyline2.begin()), std::make_move_iterator(polyline2.end()));

  return fromPolyline(polyline, order ? order : curve1.order() + curve2.order());
}

Curve Curve::fromPolyline(const PointVector& polyline, unsigned int order)
{
  const unsigned N = std::min(order ? order + 1 : polyline.size(), polyline.size());

  if (polyline.size() < 2)
    throw std::logic_error{"Polyline must have at least two points."};
  if (N == 2)
    return Curve(PointVector{polyline.front(), polyline.back()});

  // Select N points that most influence the shape of the polyline,
  // either based on a specified order or using the full polyline.
  auto simplified = _polylineSimplify(polyline, N);

  // Initialize vector t where each element represents a normalized cumulative
  // distance between consecutive simplified points along the simplified polyline.
  Eigen::VectorXd t(N);
  Eigen::MatrixXd P(N, 2), M = bernsteinCoeffs(N);
  for (unsigned k{}; k < N; k++)
  {
    P.row(k) = simplified[k];
    t(k) = k == 0 ? 0 : t(k - 1) + (P.row(k) - P.row(k - 1)).norm();
  }
  t /= t(N - 1);

  // Compute the control points for a Bezier curve such that C(t_i) = P_i for all i,
  // by solving the matrix form of the Bezier curve equation.
  auto getCurve = [&M, &P](const Eigen::VectorXd& t) {
    Eigen::MatrixXd T(t.size(), t.size());
    for (unsigned k{}; k < t.size(); k++)
      T.row(k) = _powSeries(t(k), t.size());
    return Curve(M.inverse() * (T.transpose() * T).inverse() * T.transpose() * P);
  };

  // Calculate the root mean square distance (RMSD) between the polyline
  // representation of the curve and the original polyline points.
  auto rmsd = [&polyline](const Curve& c) {
    double rmsd{};
    auto polyline_c = c.polyline();
    for (const auto& p : polyline_c)
      rmsd += _pow(_polylineDist(polyline, p), 2);
    return std::sqrt(rmsd / polyline_c.size());
  };

  // Calculate the absolute difference in length between the polyline representation
  // of the curve and the original polyline points.
  const double polyline_length = _polylineLength(polyline);
  auto length_diff = [&polyline_length](const Curve& c) {
    return std::fabs(_polylineLength(c.polyline()) - polyline_length);
  };

  struct CostFunctor : public Eigen::DenseFunctor<double>
  {
    using GetCurveFun = std::function<Curve(Eigen::VectorXd)>;
    using CostFun = std::function<double(Curve)>;
    GetCurveFun getCurve;
    CostFun f1, f2;

    CostFunctor(int N, GetCurveFun getCurve, CostFun f1, CostFun f2)
        : DenseFunctor<double>(N - 2, 2), getCurve(getCurve), f1(f1), f2(f2)
    {
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
    {
      auto c = getCurve((Eigen::VectorXd(inputs() + 2) << 0, x, 1).finished());
      fvec(0) = f1(c);
      fvec(1) = f2(c);
      return 0;
    }
  };

  CostFunctor costFun(N, getCurve, rmsd, length_diff);
  Eigen::NumericalDiff<CostFunctor> num_diff(costFun);

  Eigen::MatrixXd J(2, N - 2);
  Eigen::VectorXd fvec(2), x = t.segment(1, N - 2);
  costFun(x, fvec);
  num_diff.df(x, J);

  struct ParetoData
  {
    Eigen::VectorXd f, x;
    Eigen::MatrixXd J;
    double alpha;
  };

  constexpr double alpha_init = 0.1;
  constexpr unsigned pareto_max_size = 10;
  std::vector<ParetoData> pareto_front{{fvec, x, J, alpha_init}};
  pareto_front.reserve(2 * pareto_max_size);

  for (bool finished{false}; !finished;)
  {
    finished = true;

    for (auto& sample : pareto_front)
    {
      if (sample.alpha < std::sqrt(_epsilon))
        continue;
      finished = false;

      do
      {
        x = sample.x - sample.alpha * (sample.J.row(0) + sample.J.row(1)).transpose().normalized();
        sample.alpha /= 2;
      } while (x.minCoeff() < 0.0 || x.maxCoeff() > 1.0 || (x.array().head(N - 3) >= x.array().tail(N - 3)).any());

      costFun(x, fvec);
      if (std::any_of(pareto_front.begin(), pareto_front.end(),
                      [&fvec](const auto& p) { return (p.f.array() <= fvec.array()).all(); }))
        continue;

      num_diff.df(x, J);
      pareto_front.push_back({fvec, x, J, 2 * sample.alpha});
    }

    std::sort(pareto_front.begin(), pareto_front.end(), [](const auto& p1, const auto& p2) {
      return p1.f(0) == p2.f(0) ? p1.f(1) < p2.f(1) : p1.f(0) < p2.f(0);
    });

    double best_f1 = std::numeric_limits<double>::max();
    auto erase_it = std::remove_if(pareto_front.begin(), pareto_front.end(), [&best_f1](const auto& p) {
      return p.f(1) <= best_f1 ? (best_f1 = p.f(1), false) : true;
    });

    if (std::distance(erase_it, pareto_front.begin() + pareto_max_size) < 0)
      erase_it = pareto_front.begin() + pareto_max_size;
    pareto_front.erase(erase_it, pareto_front.end());
  }

  return getCurve((Eigen::VectorXd(N) << 0, pareto_front.front().x, 1).finished());
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
    for (unsigned k{}, binomial = 1; k < n; binomial = binomial * (n - k - 1) / (k + 1), k++)
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
      for (unsigned k{}; k < n; k++)
        splitting_coeffs_right_[n].block(0, n - 1 - k, n - k, 1) =
            temp_splitting_coeffs_left.diagonal(-static_cast<int>(k)).reverse();
    }
    return splitting_coeffs_right_[n];
  }

  Curve::Coeffs coeffs(Coeffs::Zero(n, n));
  Curve::Coeffs temp_splitting_coeffs_left = splittingCoeffsLeft(n, t);
  for (unsigned k{}; k < n; k++)
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
