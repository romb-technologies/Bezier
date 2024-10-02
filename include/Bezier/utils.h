#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <vector>

#include "declarations.h"

namespace Bezier
{
namespace Utils
{

/// Precision for numerical methods
const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

/// 2D analogue cross operation
inline double cross(const Vector& u, const Vector& v) { return u.x() * v.y() - u.y() * v.x(); }

/// Calculate unsigned power of number 2
inline unsigned exp2(unsigned exp) { return 1 << exp; }

/// Calculate power for integer exponents
template <typename T> inline T pow(T base, unsigned exp)
{
  T result = exp & 1 ? base : 1;
  while (exp >>= 1)
  {
    base *= base;
    if (exp & 1)
      result *= base;
  }
  return result;
}

/// Compute binomial coefficient
inline Eigen::RowVectorXd powSeries(double base, unsigned exp)
{
  Eigen::RowVectorXd power_series(exp);
  power_series(0) = 1;
  for (unsigned k = 1; k < exp; k++)
    power_series(k) = power_series(k - 1) * base;
  return power_series;
}

/// Concatenate two vectors
template <typename T> inline std::vector<T> concatenate(std::vector<T> v1, std::vector<T> v2)
{
  v1.reserve(v1.size() + v2.size());
  v1.insert(v1.end(), std::make_move_iterator(v2.begin()), std::make_move_iterator(v2.end()));
  return v1;
}

/// Evaluate Chebyshev polynomial at point t
inline double evaluateChebyshev(double t, const Eigen::VectorXd& coeffs)
{
  t = 2 * t - 1;
  double tn{t}, tn_1{1}, res{coeffs(0) + coeffs(1) * t};
  for (unsigned k{2}; k < coeffs.size(); k++)
  {
    std::swap(tn_1, tn);
    tn = 2 * t * tn_1 - tn;
    res += coeffs(k) * tn;
  }
  return res;
}

/// Distance between two points
inline double dist(const Point& p1, const Point& p2) { return (p1 - p2).norm(); }

/// Distance between a point and a segment
inline double dist(const Point& p1, const Point& p2, const Point& point)
{
  Vector u = p2 - p1;
  Vector v = point - p1;
  double t = u.dot(v) / u.squaredNorm();
  if (t < 0)
    return v.norm();
  if (t > 1)
    return (point - p2).norm();
  return (p1 + t * u - point).norm();
}

/// Distance between a point and a polyline
inline double dist(const PointVector& polyline, const Point& point)
{
  double min_dist = std::numeric_limits<double>::max();
  for (size_t k{1}; k < polyline.size(); k++)
    min_dist = std::min(min_dist, dist(polyline[k - 1], polyline[k], point));
  return min_dist;
}

/// Upper limit for curve deviation from the line segment connecting first and last point
inline double maxDeviation(const Eigen::MatrixX2d& cp)
{
  // for N_ == 10, coeff is 0.9922, so we ignore it for higher orders
  const double coeff{cp.rows() >= 10 ? 1 : 1 - std::exp2(2. - cp.rows())};
  auto deviationFunc = [&cp](int k) { return dist(cp.row(0), cp.row(cp.rows() - 1), cp.row(k + 1)); };
  return cp.rows() <= 2 ? 0 : coeff * Eigen::VectorXd::NullaryExpr(cp.rows() - 2, deviationFunc).maxCoeff();
}

/// Sort indices of polyline points by their contribution to the polyline shape
std::vector<unsigned> visvalingamWyatt(const PointVector& polyline);

/// Simplify polyline to N points
PointVector polylineSimplify(const PointVector& polyline, unsigned N);

/// Length of a polyline
inline double polylineLength(const PointVector& polyline)
{
  double length{};
  for (size_t k{1}; k < polyline.size(); k++)
    length += (polyline[k] - polyline[k - 1]).norm();
  return length;
}

/// Find solutions to polynomial equation (limited to [0, 1])
std::vector<double> solvePolynomial(const Eigen::VectorXd& polynomial);

} // namespace Utils
} // namespace Bezier

#endif // UTILS_H
