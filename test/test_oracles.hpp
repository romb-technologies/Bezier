#ifndef BEZIER_TEST_ORACLES_HPP
#define BEZIER_TEST_ORACLES_HPP

// Independent reference implementations used as oracles in tests.
// None of these may call into the library under test — they exist so the
// test suite survives the planned v0.4.0 algorithm rework.

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
/// - kChar: characterization anchors pinned to current-algorithm output
///   (matches Utils::epsilon = sqrt(machine epsilon))
/// - kFit:  approximation methods (fromPolyline/offsetCurve/joinCurves),
///   expressed as a fraction of the bounding-box diagonal
constexpr double kGeom = 1e-9;
constexpr double kChar = 1.5e-8;
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

/// Parameter of the closest curve point to p, found by dense sampling.
inline double denseArgminDistance(const PointVector& cp, const Point& p, size_t n = 10000)
{
  double best_t{}, best_dist{(deCasteljau(cp, 0.0) - p).norm()};
  for (size_t i = 1; i <= n; i++)
  {
    double t = static_cast<double>(i) / static_cast<double>(n);
    double dist = (deCasteljau(cp, t) - p).norm();
    if (dist < best_dist)
    {
      best_dist = dist;
      best_t = t;
    }
  }
  return best_t;
}

} // namespace Oracles
} // namespace Bezier

#endif // BEZIER_TEST_ORACLES_HPP
