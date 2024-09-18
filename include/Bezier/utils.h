#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Dense>
#include <vector>

#include "declarations.h"

#ifndef __cpp_lib_make_unique
#include <memory>
namespace std
{
template <typename T, typename... Args> inline std::unique_ptr<T> make_unique(Args&&... args)
{
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
} // namespace std
#endif

namespace Bezier
{
namespace Utils
{
/*!
 * \brief Precision for numerical methods
 */
const double _epsilon = std::sqrt(std::numeric_limits<double>::epsilon());

struct _PolynomialRoots : public std::vector<double>
{
  explicit _PolynomialRoots(unsigned reserve) { std::vector<double>::reserve(reserve); }
  void clear() {}          // no-op so that PolynomialSolver::RealRoots() doesn't clear it
  void push_back(double t) // only allow valid roots
  {
    if (t >= 0 && t <= 1)
      std::vector<double>::push_back(t);
  }
};

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

inline Eigen::VectorXd _trimZeroes(const Eigen::VectorXd& vec)
{
  auto idx = vec.size();
  while (idx && std::abs(vec(idx - 1)) < _epsilon)
    --idx;
  return vec.head(idx);
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

} // namespace Utils
} // namespace Bezier

#endif // UTILS_H
