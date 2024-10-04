#include "Bezier/bezier.h"
#include "Bezier/utils.h"

#include <unsupported/Eigen/FFT>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/NumericalDiff>

using namespace Bezier;
namespace bu = Bezier::Utils;

Curve::Curve(Eigen::MatrixX2d points) : control_points_(std::move(points)), N_(control_points_.rows()) {}

Curve::Curve(const PointVector& points) : control_points_(points.size(), 2), N_(points.size())
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
      subcurves.emplace_back(splittingCoeffsRight(N_) * cp, (t1 + t2) / 2, t2);
      subcurves.emplace_back(splittingCoeffsLeft(N_) * cp, t1, (t1 + t2) / 2);
    }
  }

  return cached_polyline_.value();
}

double Curve::length() const { return length(1.0); }

double Curve::length(double t) const
{
  if (t < 0.0 || t > 1.0)
    throw std::logic_error{"Length can only be calculated for t within [0.0, 1.0] range."};

  if (cached_chebyshev_coeffs_)
    return bu::evaluateChebyshev(t, cached_chebyshev_coeffs_.value());
  auto& chebyshev = cached_chebyshev_coeffs_.emplace();

  constexpr unsigned START_LOG_N = 10;
  unsigned log_n = START_LOG_N - 1;
  unsigned n = bu::exp2(START_LOG_N - 1);

  Eigen::VectorXd derivative_cache(2 * n + 1);
  auto updateDerivativeCache = [this, &derivative_cache](double n) {
    derivative_cache.conservativeResize(n + 1);
    derivative_cache.tail(n / 2) = Eigen::VectorXd::NullaryExpr(
        n / 2, [n, this](int k) { return derivativeAt((1 + std::cos((2 * k + 1) * M_PI / n)) / 2).norm(); });
  };

  derivative_cache.head(2) << derivativeAt(1.0).norm(), derivativeAt(0.0).norm();
  for (unsigned k{2}; k <= n; k *= 2)
    updateDerivativeCache(k);

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
      auto lin_spaced = Eigen::ArrayXi::LinSpaced(bu::exp2(k - 1), 0, bu::exp2(k - 1) - 1);
      auto index_c = bu::exp2(log_n + 1 - (k + 1)) + lin_spaced * bu::exp2(log_n + 1 - k);
      auto index_dc = bu::exp2(k - 1) + 1 + lin_spaced;
      for (unsigned i{}; i < lin_spaced.size(); i++)
        coeff(index_c(i)) = coeff(N - index_c(i)) = derivative_cache(index_dc(i)) / n;
    }

    Eigen::VectorXcd fft_out;
    Eigen::FFT<double>().fwd(fft_out, coeff);
    chebyshev.resize(n);
    chebyshev << 0, (fft_out.real().head(n - 1) - fft_out.real().segment(2, n - 1)).array() /
                        Eigen::ArrayXd::LinSpaced(n - 1, 4, 4 * (n - 1));
  } while (std::fabs(chebyshev.tail<1>()[0]) > bu::epsilon);

  while (std::fabs(chebyshev.tail<1>()[0]) < bu::epsilon)
    chebyshev.conservativeResize(chebyshev.size() - 1);
  chebyshev(0) -= bu::evaluateChebyshev(0, chebyshev);

  return bu::evaluateChebyshev(t, chebyshev);
}

double Curve::length(double t1, double t2) const { return length(t2) - length(t1); }

double Curve::iterateByLength(double t, double s) const
{
  if (std::fabs(s) < bu::epsilon) // no-op
    return t;

  double s_t = length(t);

  std::pair<double, double> lbracket, rbracket, guess{t, 0.0};
  if (s < 0)
  {
    lbracket = {0.0, -s_t};
    if (s < lbracket.second + bu::epsilon) // out-of-scope
      return 0.0;
    rbracket = guess;
  }
  else // s > 0
  {
    rbracket = {1.0, length() - s_t};
    if (s > rbracket.second - bu::epsilon) // out-of-scope
      return 1.0;
    lbracket = guess;
  }

  while (std::fabs(guess.second - s) > bu::epsilon)
  {
    // Halley's method
    double f = guess.second - s;
    double f_d = derivativeAt(guess.first).norm();
    double f_d2 = derivativeAt(2, guess.first).norm();
    guess.first -= (2 * f * f_d) / (2 * f_d * f_d - f * f_d2);

    // root bracketing, if not in bounds, use bisection method
    if (guess.first <= lbracket.first || guess.first >= rbracket.first)
      guess.first = (lbracket.first + rbracket.first) / 2;

    if (rbracket.first - lbracket.first < bu::epsilon)
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
  return N_ == 0 ? Point(0, 0) : (bu::powVector(t, N_) * bernsteinCoeffs(N_) * control_points_).transpose();
}

Eigen::MatrixX2d Curve::valueAt(const std::vector<double>& t_vector) const
{
  auto t_map = Eigen::Map<const Eigen::VectorXd>(t_vector.data(), t_vector.size());
  return bu::powMatrix(t_map, N_) * bernsteinCoeffs(N_) * control_points_;
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
  if (cached_roots_)
    return cached_roots_.value();

  Eigen::MatrixXd bezier_polynomial = bernsteinCoeffs(N_) * control_points_;
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
  extremes.bottomRows(2) << control_points_.row(0), control_points_.row(N_ - 1);
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
    subcurves.emplace_back(splittingCoeffsLeft(N_, sorted_t[k]) * leftover_cp);
    leftover_cp = splittingCoeffsRight(N_, sorted_t[k]) * leftover_cp;
    std::for_each(sorted_t.begin() + k + 1, sorted_t.end(), [t = sorted_t[k]](double& x) { x = (x - t) / (1 - t); });
  }
  subcurves.emplace_back(std::move(leftover_cp));
  return subcurves;
}

std::vector<Curve> Curve::splitCurve(double t) const
{
  return {Curve(splittingCoeffsLeft(N_, t) * control_points_), Curve(splittingCoeffsRight(N_, t) * control_points_)};
}

PointVector Curve::intersections(const Curve& curve) const
{
  using CP = Eigen::MatrixX2d;

  std::vector<std::pair<CP, CP>> cp_pairs;
  if (this != &curve)
    cp_pairs.emplace_back(control_points_, curve.control_points_);
  else
  {
    // For self intersections divide curve into subcurves at extrema and create all pairs of subcurves
    // - leave epsilon sized gap between subcurves to avoid false positives
    auto t = extrema();
    t.resize(t.size() * 2);
    for (unsigned k{}; k < t.size() / 2; k++)
      t[k + t.size() / 2] = (t[k] -= bu::epsilon / 2) + bu::epsilon;
    std::sort(t.begin(), t.end());
    auto subcurves = splitCurve(t);
    for (auto c1 = subcurves.begin(); c1 < subcurves.end(); c1 += 2)
      for (auto c2 = c1 + 2; c2 < subcurves.end(); c2 += 2)
        cp_pairs.emplace_back(c1->control_points_, c2->control_points_);
  }

  auto insertPairs = [&cp_pairs](auto&& scp1, auto&& scp2) {
    for (auto cp1 : scp1)
      for (auto cp2 : scp2)
        cp_pairs.emplace_back(cp1, cp2);
  };

  auto splitCP = [](const CP& cp) -> std::array<CP, 2> {
    return {splittingCoeffsRight(cp.rows()) * cp, splittingCoeffsLeft(cp.rows()) * cp};
  };

  PointVector intersections;
  auto insertIntersection = [&intersections](const CP& cp1, const CP& cp2) {
    /// Intersection of two segments (Victor Lecomte - Handbook of geometry for competitive programmers)
    auto a1 = cp1.row(0), a2 = cp1.row(cp1.rows() - 1);
    auto b1 = cp2.row(0), b2 = cp2.row(cp2.rows() - 1);
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

    // Split each curve until is it flat enough to be represented as segment
    bool cp1_done{bu::maxDeviation(cp1) < bu::epsilon};
    bool cp2_done{bu::maxDeviation(cp2) < bu::epsilon};
    if (cp1_done && cp2_done)
      insertIntersection(cp1, cp2);
    else if (cp1_done)
      insertPairs(std::array{cp1}, splitCP(cp2));
    else if (cp2_done)
      insertPairs(splitCP(cp1), std::array{cp2});
    else
      insertPairs(splitCP(cp1), splitCP(cp2));
  }

  return intersections;
}

double Curve::projectPoint(const Point& point) const
{
  if (!cached_projection_polynomial_part_ || !cached_projection_polynomial_derivative_)
  {
    Eigen::MatrixXd curve_polynomial = (bernsteinCoeffs(N_) * control_points_);
    Eigen::MatrixXd derivate_polynomial = (bernsteinCoeffs(N_ - 1) * derivative().control_points_);

    Eigen::VectorXd polynomial_part = Eigen::VectorXd::Zero(curve_polynomial.rows() + derivate_polynomial.rows() - 1);
    for (unsigned k{}; k < curve_polynomial.rows(); k++)
      polynomial_part.middleRows(k, derivate_polynomial.rows()) +=
          derivate_polynomial * curve_polynomial.row(k).transpose();

    cached_projection_polynomial_part_.emplace(std::move(polynomial_part));
    cached_projection_polynomial_derivative_.emplace(std::move(derivate_polynomial));
  }

  Eigen::VectorXd polynomial = cached_projection_polynomial_part_.value();
  polynomial.topRows(cached_projection_polynomial_derivative_->rows()) -=
      cached_projection_polynomial_derivative_.value() * point;

  double min_t{0.0}, min_dist{bu::dist(point, valueAt(0.0))};

  for (auto t : bu::concatenate(bu::solvePolynomial(polynomial), {1.0}))
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
  Eigen::Matrix2Xd derivatives(2, c_order + 1);
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

Curve Curve::joinCurves(const Curve& curve1, const Curve& curve2, unsigned order)
{
  if (order == 1)
    return Curve(PointVector{curve1.control_points_.row(0), curve2.control_points_.row(curve2.N_ - 1)});
  return fromPolyline(bu::concatenate(curve1.polyline(), curve2.polyline()),
                      order ? order : curve1.order() + curve2.order());
}

Curve Curve::fromPolyline(const PointVector& polyline, unsigned order)
{
  const unsigned N = std::min(order ? order + 1 : polyline.size(), polyline.size());

  if (polyline.size() < 2)
    throw std::logic_error{"Polyline must have at least two points."};
  if (N == 2)
    return Curve(PointVector{polyline.front(), polyline.back()});

  // Sort the polyline points by their contribution to the Visvalingam-Whyatt
  // simplification algorithm, and keep the N most contributing points in original order.
  auto vw = bu::visvalingamWyatt(polyline);
  std::sort(vw.begin(), vw.begin() + N);

  // Divide polyline into subparts where the simplified polyline points are located.
  std::vector<PointVector> subpolylines;
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
  auto getCurve = [&P, M = bernsteinCoeffs(N)](const Eigen::VectorXd& t) {
    Eigen::MatrixXd T = bu::powMatrix(t, t.size());
    return Curve(M.inverse() * (T.transpose() * T).inverse() * T.transpose() * P);
  };

  // Cost functor calculates RMSD and length difference for each subcurve/subpolyline
  // divided at parameter t, where C(t_i) = P_i.
  struct CostFunctor : public Eigen::DenseFunctor<double>
  {
    using GetCurveFun = std::function<Curve(const Eigen::VectorXd&)>;
    GetCurveFun getCurve;
    std::vector<PointVector> subpolylines;

    CostFunctor(int N, GetCurveFun getCurve, const std::vector<PointVector>& subpolylines)
        : DenseFunctor<double>(N - 2, 2 * N - 2), getCurve(getCurve), subpolylines(subpolylines)
    {
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
    {
      auto rmsd = [](const Curve& curve, const PointVector& polyline) {
        double rmsd{};
        auto polyline_c = curve.polyline();
        for (const auto& p : polyline_c)
          rmsd += bu::pow(bu::dist(polyline, p), 2);
        return std::sqrt(rmsd / polyline_c.size());
      };

      auto length_diff = [](const Curve& curve, const PointVector& polyline) {
        return std::fabs(bu::polylineLength(curve.polyline()) - bu::polylineLength(polyline));
      };

      auto curve = getCurve((Eigen::VectorXd(inputs() + 2) << 0, x, 1).finished());
      auto subcurves = curve.splitCurve(std::vector<double>(x.data(), x.data() + inputs()));
      for (unsigned k = 0; k < subcurves.size(); k++)
      {
        fvec(k) = rmsd(subcurves[k], subpolylines[k]);
        fvec(values() / 2 + k) = length_diff(subcurves[k], subpolylines[k]);
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

void Curve::resetCache()
{
  N_ = control_points_.rows();
  cached_derivative_.reset();
  cached_roots_.reset();
  cached_bounding_box_.reset();
  cached_polyline_.reset();
  cached_polyline_t_.reset();
  cached_projection_polynomial_part_.reset();
  cached_projection_polynomial_derivative_.reset();
  cached_chebyshev_coeffs_.reset();
}

Curve::CoeffsMap Curve::bernstein_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_left_ = CoeffsMap();
Curve::CoeffsMap Curve::splitting_coeffs_right_ = CoeffsMap();
Curve::CoeffsMap Curve::elevate_order_coeffs_ = CoeffsMap();
Curve::CoeffsMap Curve::lower_order_coeffs_ = CoeffsMap();

Curve::Coeffs Curve::bernsteinCoeffs(unsigned n)
{
  auto coeffs = [n]() -> Coeffs {
    Coeffs coeffs = Coeffs::Zero(n, n);
    coeffs.diagonal(-1).setLinSpaced(-1, -static_cast<int>(n - 1));
    coeffs = coeffs.exp();
    coeffs.array().colwise() *= coeffs.row(n - 1).transpose().array().abs();
    return coeffs;
  };
  return bernstein_coeffs_.try_emplace(n, bu::lazyFunctor(coeffs)).first->second;
}

Curve::Coeffs Curve::splittingCoeffsLeft(unsigned n, double t)
{
  auto scl = [n](double t) -> Coeffs {
    return bernsteinCoeffs(n).inverse() * bu::powVector(t, n).asDiagonal() * bernsteinCoeffs(n);
  };
  return t == 0.5 ? splitting_coeffs_left_.try_emplace(n, bu::lazyFunctor(scl, 0.5)).first->second : scl(t);
}

Curve::Coeffs Curve::splittingCoeffsRight(unsigned n, double t)
{
  auto scr = [n](double t) -> Coeffs {
    Coeffs coeffs = splittingCoeffsLeft(n, t);
    for (unsigned k{}; k < n; k++)
      coeffs.col(n - 1 - k).head(n - k) = coeffs.diagonal(-static_cast<int>(k)).reverse();
    coeffs.triangularView<Eigen::StrictlyLower>().setZero();
    return coeffs;
  };
  return t == 0.5 ? splitting_coeffs_right_.try_emplace(n, bu::lazyFunctor(scr, 0.5)).first->second : scr(t);
}

Curve::Coeffs Curve::elevateOrderCoeffs(unsigned n)
{
  auto eoc = [n]() -> Coeffs {
    Coeffs coeffs = Coeffs::Zero(n + 1, n);
    coeffs.diagonal(-1) = coeffs.diagonal().setLinSpaced(1, 1. / n).reverse();
    return coeffs;
  };
  return elevate_order_coeffs_.try_emplace(n, bu::lazyFunctor(eoc)).first->second;
}

Curve::Coeffs Curve::lowerOrderCoeffs(unsigned n)
{
  auto loc = [n]() -> Coeffs { return elevateOrderCoeffs(n - 1).completeOrthogonalDecomposition().pseudoInverse(); };
  return lower_order_coeffs_.try_emplace(n, bu::lazyFunctor(loc)).first->second;
}
