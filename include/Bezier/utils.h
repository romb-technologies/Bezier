#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <vector>

#include "declarations.h"

namespace Bezier
{
namespace Utils
{
/*!
 * \brief Precision for numerical methods
 */
const double _epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

inline unsigned _exp2(unsigned exp) { return 1 << exp; }

template <typename T> inline T _pow(T base, unsigned exp)
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

inline Eigen::RowVectorXd _powSeries(double base, unsigned exp)
{
  Eigen::RowVectorXd power_series(exp);
  power_series(0) = 1;
  for (unsigned k = 1; k < exp; k++)
    power_series(k) = power_series(k - 1) * base;
  return power_series;
}

template <typename T> inline std::vector<T> _concatenate(std::vector<T>&& v1, std::vector<T>&& v2)
{
  v1.reserve(v1.size() + v2.size());
  v1.insert(v1.end(), std::make_move_iterator(v2.begin()), std::make_move_iterator(v2.end()));
  return std::move(v1);
}

inline double _polylineLength(const PointVector& polyline)
{
  double length{};
  for (size_t k{1}; k < polyline.size(); k++)
    length += (polyline[k] - polyline[k - 1]).norm();
  return length;
}

inline double _polylineDist(const PointVector& polyline, const Point& point)
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

PointVector _polylineSimplify(const PointVector& polyline, unsigned N);

std::vector<double> _solvePolynomial(const Eigen::VectorXd& polynomial);

} // namespace Utils
} // namespace Bezier

#endif // UTILS_H
