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
template <typename T> inline std::vector<T> concatenate(std::vector<T>&& v1, std::vector<T>&& v2)
{
  v1.reserve(v1.size() + v2.size());
  v1.insert(v1.end(), std::make_move_iterator(v2.begin()), std::make_move_iterator(v2.end()));
  return std::move(v1);
}

/// Length of a polyline
inline double polylineLength(const PointVector& polyline)
{
  double length{};
  for (size_t k{1}; k < polyline.size(); k++)
    length += (polyline[k] - polyline[k - 1]).norm();
  return length;
}

/// Distance between a point and a polyline
inline double polylineDist(const PointVector& polyline, const Point& point)
{
  auto distSeg = [&point](const Point& p1, const Point& p2) {
    Vector u = p2 - p1;
    Vector v = point - p1;
    double t = u.dot(v) / u.squaredNorm();
    if (t < 0)
      return v.squaredNorm();
    if (t > 1)
      return (point - p2).squaredNorm();
    return (p1 + t * u - point).squaredNorm();
  };

  double dist = std::numeric_limits<double>::max();
  for (size_t k{1}; k < polyline.size(); k++)
    dist = std::min(dist, distSeg(polyline[k - 1], polyline[k]));
  return dist;
}

/// Sort indices of polyline points by their contribution to the polyline shape
std::vector<unsigned> visvalingamWyatt(const PointVector& polyline);

/// Simplify polyline to N points
inline PointVector polylineSimplify(const PointVector& polyline, unsigned N)
{
  if (polyline.size() < 2)
    throw std::logic_error{"Polyline must have at least two points."};
  if (polyline.size() < N)
    return PointVector(polyline);
  if (N == 2)
    return PointVector{polyline.front(), polyline.back()};

  auto by_contribution = visvalingamWyatt(polyline);
  std::sort(by_contribution.begin(), by_contribution.begin() + N);
  PointVector simplified(N);
  for (size_t k{0}; k < N; k++)
    simplified[k] = polyline[by_contribution[k]];
  return simplified;
}

/// Find solutions to polynomial equation (limited to [0, 1])
std::vector<double> solvePolynomial(const Eigen::VectorXd& polynomial);

} // namespace Utils
} // namespace Bezier

#endif // UTILS_H
