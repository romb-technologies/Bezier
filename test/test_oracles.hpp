#ifndef BEZIER_TEST_ORACLES_HPP
#define BEZIER_TEST_ORACLES_HPP

// Independent reference implementations used as oracles in tests.
// None of these may call into the library under test — they exist so the
// test suite survives the planned v0.4.0 algorithm rework.

#include <cmath>

#include <Eigen/Dense>

#include "Bezier/declarations.h"

namespace Bezier
{
namespace Oracles
{

/// Tolerance guide:
/// - algebraically exact operations (control-point round trips, derivative
///   control points): use EXPECT_EQ / EXPECT_DOUBLE_EQ, no constant needed
/// - kGeom: geometric identities (de Casteljau agreement, split continuity,
///   reverse mapping); absolute, adequate for coordinates of magnitude ~1e2
/// - kFit:  approximation methods (fromPolyline/offsetCurve/joinCurves),
///   expressed as a fraction of the bounding-box diagonal
constexpr double kGeom = 1e-9;
constexpr double kFit = 0.05;

/// Pointwise de Casteljau evaluation of a Bezier curve given its control points.
inline Point deCasteljau(PointVector cp, double t)
{
  for (size_t k = cp.size(); k > 1; k--)
    for (size_t i = 0; i + 1 < k; i++)
      cp[i] = (1 - t) * cp[i] + t * cp[i + 1];
  return cp[0];
}

/// Arc length of the curve segment [t1, t2] via dense chord sampling.
/// Relative accuracy ~1/n^2; with n = 1e5 good to ~1e-9 relative.
inline double chordLength(const PointVector& cp, double t1 = 0.0, double t2 = 1.0, size_t n = 100000)
{
  double length{};
  Point prev = deCasteljau(cp, t1);
  for (size_t i = 1; i <= n; i++)
  {
    Point next = deCasteljau(cp, t1 + (t2 - t1) * static_cast<double>(i) / static_cast<double>(n));
    length += (next - prev).norm();
    prev = next;
  }
  return length;
}

/// Central-difference derivative of a point-valued function of t.
template <typename F> inline Vector centralDiff(F&& f, double t, double h = 1e-6)
{
  return (f(t + h) - f(t - h)) / (2 * h);
}

/// Control points of the derivative (hodograph): n * (P_{i+1} - P_i).
inline PointVector derivativeControlPoints(const PointVector& cp)
{
  PointVector d;
  if (cp.size() < 2)
    return d;
  const double n = static_cast<double>(cp.size() - 1);
  d.reserve(cp.size() - 1);
  for (size_t i = 0; i + 1 < cp.size(); i++)
    d.emplace_back(n * (cp[i + 1] - cp[i]));
  return d;
}

/// Signed curvature kappa(t) = (x'y'' - y'x'') / (x'^2 + y'^2)^{3/2}, from exact
/// Bezier derivatives — no finite-difference noise.
inline double curvature(const PointVector& cp, double t)
{
  PointVector d1 = derivativeControlPoints(cp);
  PointVector d2 = derivativeControlPoints(d1);
  Vector r1 = deCasteljau(d1, t), r2 = deCasteljau(d2, t);
  return (r1.x() * r2.y() - r1.y() * r2.x()) / std::pow(r1.squaredNorm(), 1.5);
}

/// d kappa / dt, closed form: (a'·s - 1.5·a·s') / s^{5/2} with a = x'y''-y'x'',
/// s = ‖r'‖², a' = x'y'''-y'x''', s' = 2(x'x''+y'y'').
inline double curvatureDerivative(const PointVector& cp, double t)
{
  PointVector d1 = derivativeControlPoints(cp);
  PointVector d2 = derivativeControlPoints(d1);
  PointVector d3 = derivativeControlPoints(d2);
  Vector r1 = deCasteljau(d1, t), r2 = deCasteljau(d2, t), r3 = deCasteljau(d3, t);
  double s = r1.squaredNorm();
  double a = r1.x() * r2.y() - r1.y() * r2.x();
  double ap = r1.x() * r3.y() - r1.y() * r3.x();
  double sp = 2 * (r1.x() * r2.x() + r1.y() * r2.y());
  return (ap * s - 1.5 * a * sp) / std::pow(s, 2.5);
}

} // namespace Oracles
} // namespace Bezier

#endif // BEZIER_TEST_ORACLES_HPP
