/*
 * Copyright 2024 Mirko Kokot
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef UTILS_H
#define UTILS_H

#include "Bezier/declarations.h"

#include <vector>

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
    if (base *= base; exp & 1)
      result *= base;
  return result;
}

/// A vector of powers of a given base
inline Eigen::RowVectorXd powVector(double base, unsigned exp)
{
  std::function<double(int)> powFunc = [x = 1., base](int k) mutable { return k ? x *= base : x; };
  return Eigen::RowVectorXd::NullaryExpr(exp, powFunc);
}

/// A matrix of powers of given bases
inline Eigen::MatrixXd powMatrix(const Eigen::VectorXd& base, unsigned exp)
{
  Eigen::MatrixXd power_matrix(base.size(), exp);
  power_matrix.col(0).setOnes();
  for (unsigned k = 1; k < exp; k++)
    power_matrix.col(k) = power_matrix.col(k - 1).cwiseProduct(base);
  return power_matrix;
}

/// Concatenate two vectors
template <typename T> inline std::vector<T> concatenate(std::vector<T>&& v1, std::vector<T>&& v2)
{
  v1.reserve(v1.size() + v2.size());
  v1.insert(v1.end(), std::make_move_iterator(v2.begin()), std::make_move_iterator(v2.end()));
  return std::move(v1);
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
  double t = std::clamp(u.dot(v) / u.squaredNorm(), 0., 1.);
  return (p1 + t * u - point).norm();
}

/// Distance between a point and a polyline
inline double dist(const std::vector<Point>& polyline, const Point& point)
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
std::vector<unsigned> visvalingamWyatt(const std::vector<Point>& polyline);

/// Simplify polyline to N points
std::vector<Point> polylineSimplify(const std::vector<Point>& polyline, unsigned N);

/// Length of a polyline
inline double polylineLength(const std::vector<Point>& polyline)
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
