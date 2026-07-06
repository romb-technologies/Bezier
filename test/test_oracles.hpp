#ifndef BEZIER_TEST_ORACLES_HPP
#define BEZIER_TEST_ORACLES_HPP

// Independent reference implementations used as oracles in tests.
// None of these may call into the library under test — they exist so the
// test suite survives the planned v0.4.0 algorithm rework.

#include <cmath>
#include <limits>
#include <vector>

#include <Eigen/Dense>

#include "Bezier/declarations.h"

namespace Bezier
{
namespace Oracles
{

/// Tolerance policy: kGeom (machine precision) is the default. Every library
/// method — geometric OR iterative — converges to Utils::epsilon (= sqrt of the
/// machine epsilon) internally, so it is the right bound for almost everything:
///   - length() grows its Chebyshev expansion until the tail coeff < epsilon
///   - step() runs Halley's method until |s - ds| < epsilon
///   - projectPoint()/distance() are bounded by root-finder / projection precision
/// Depart from kGeom only with a reason:
///   - exact: control-point copies, orders, sizes -> EXPECT_EQ / EXPECT_DOUBLE_EQ
///   - kAlgebraic: FINER than kGeom. Pure arithmetic / matrix identities with no
///                 geometric cancellation — Utils::pow, Chebyshev evaluation,
///                 Bernstein / split / order-change matrices. A few ULP above eps.
///   - kGeom:      closed-form geometric identities AND converged iterative
///                 methods (length/step/projection/distance). Stacked-eps sites
///                 (e.g. step calls length) use a small named multiple of kGeom,
///                 not a new tier.
/// Approximation methods (fromPolyline/offsetCurve/joinCurves) are NOT machine
/// exact for interior fit; their bounds live locally in approx_test.cpp, tied to
/// input resolution or the measured impl floor — never a bare percentage.
constexpr double kAlgebraic = 1e-9;
constexpr double kGeom = std::sqrt(std::numeric_limits<double>::epsilon()); // ~1.5e-8

/// Sample-count tiers for parameter sweeps (see sampleParams). Named by purpose
/// so tests stop scattering ad-hoc step counts.
constexpr unsigned kCoarseSamples = 20;    // general property sweep
constexpr unsigned kDenseSamples = 100;    // fine resolution to catch a worst-case sample
constexpr unsigned kFlatnessSamples = 1000; // matched to span/1000 geometric bounds

/// Inclusive parameter samples 0 = t_0 < ... < t_steps = 1 (steps+1 points).
/// Index-based, so unlike `for (double t{}; t <= 1.0; t += step)` it always
/// samples both endpoints exactly, regardless of float round-off.
inline std::vector<double> sampleParams(unsigned steps)
{
  std::vector<double> ts(steps + 1);
  for (unsigned i = 0; i <= steps; i++)
    ts[i] = static_cast<double>(i) / steps;
  return ts;
}

/// Pointwise de Casteljau evaluation of a Bezier curve given its control points.
/// Precondition: cp is non-empty -- every call site in this suite passes real
/// curve control points, so this is never exercised with cp = {}.
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

/// Control points of the derivative (hodograph): n * (P_{i+1} - P_i).
inline PointVector derivativeControlPoints(const PointVector& cp)
{
  PointVector d;
  if (cp.size() < 2)
    return d;
  double n = static_cast<double>(cp.size() - 1);
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
