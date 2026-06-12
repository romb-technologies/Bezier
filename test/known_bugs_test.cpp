#include "test_data.hpp"
#include "test_oracles.hpp"

#include <cmath>

#include <gtest/gtest.h>

#include "Bezier/bezier.h"
#include "Bezier/utils.h"

namespace Bezier
{

// Known bugs in the current devel-v040 state, written as DISABLED_ tests.
// Each cites the corresponding item in REVIEW-v4-devel.md and asserts the
// intended post-fix behavior. Enabling these tests (and keeping them green)
// is the acceptance gate for the v0.4.0 algorithm rework.

// REVIEW-v4-devel §2.1: Curve::intersections compares control-point matrices
// with isApprox without checking dimensions first; intersecting curves of
// different order triggers an Eigen size-mismatch assert/crash.
// Post-fix: returns the intersection points, no crash.
TEST(KnownBugsTest, DISABLED_IntersectionsDifferentOrder)
{
  Curve cubic{curvePointsAsMatrix()};
  Curve quadratic{PointVector{{84, 100}, {150, 200}, {180, 120}}};
  PointVector intersections;
  ASSERT_NO_THROW(intersections = cubic.intersections(quadratic));
  EXPECT_GE(intersections.size(), 1u);
}

// REVIEW-v4-devel §2.2: the segment-intersection predicate uses strict
// inequalities (oa*ob < 0 && oc*od < 0), so intersections exactly at a
// shared endpoint or tangential touch are missed.
// Post-fix: two curves sharing an endpoint report that point.
TEST(KnownBugsTest, DISABLED_IntersectionAtSharedEndpoint)
{
  Curve first{curvePointsAsMatrix()};
  // Second curve starts exactly where the first one ends
  Point shared = first.endPoints().second;
  Curve second{PointVector{shared, {shared.x() + 50, shared.y() + 80}, {shared.x() + 120, shared.y() + 10},
                           {shared.x() + 150, shared.y() + 90}}};

  PointVector intersections = first.intersections(second);
  ASSERT_GE(intersections.size(), 1u);
  bool found_shared = false;
  for (const Point& p : intersections)
    found_shared |= (p - shared).norm() < Utils::epsilon;
  EXPECT_TRUE(found_shared) << "Shared endpoint not reported as intersection";
}

// REVIEW-v4-devel §2.3: length() trims trailing Chebyshev coefficients
// without an index guard; for a degenerate (zero-length) curve the trim
// underflows / leaves too few coefficients for evaluateChebyshev.
// Post-fix: a degenerate curve has length 0.
TEST(KnownBugsTest, DISABLED_LengthDegenerateCurve)
{
  Curve degenerate{PointVector{{10, 10}, {10, 10}, {10, 10}, {10, 10}}};
  double length{};
  ASSERT_NO_THROW(length = degenerate.length());
  EXPECT_NEAR(length, 0.0, Utils::epsilon);
}

// REVIEW-v4-devel §2.6: at a cusp the derivative vanishes, Halley iteration
// in step() divides by zero and the NaN never satisfies the convergence
// test — infinite loop (this test relies on the ctest TIMEOUT to fail
// instead of hanging). Post-fix: step returns a finite parameter.
TEST(KnownBugsTest, DISABLED_StepAtCusp)
{
  // Doubled control point creates a cusp at t = 0.5
  Curve cusped{PointVector{{0, 0}, {100, 100}, {100, 100}, {200, 0}}};
  double t = cusped.step(0.5, 10.0);
  EXPECT_TRUE(std::isfinite(t));
  EXPECT_GE(t, 0.0);
  EXPECT_LE(t, 1.0);
}

// At zero-derivative points the normalized normal is NaN, which silently
// poisons offsetCurve output. Post-fix: normalAt returns a finite vector
// (or the API defines an explicit contract) at cusps.
TEST(KnownBugsTest, DISABLED_NormalAtCuspIsFinite)
{
  Curve cusped{PointVector{{0, 0}, {100, 100}, {100, 100}, {200, 0}}};
  Vector normal = cusped.normalAt(0.5);
  EXPECT_TRUE(std::isfinite(normal.x()));
  EXPECT_TRUE(std::isfinite(normal.y()));
}

// REVIEW-v4-devel §2.6: applyContinuity does not validate the requested
// continuity order against the curve order; beta_coeffs.size() >= N_
// produces out-of-range derivative access instead of a clear error.
// Post-fix: throws a logic error.
TEST(KnownBugsTest, DISABLED_ApplyContinuityOrderValidation)
{
  Curve curve{curvePointsAsMatrix()};
  Curve locked{rootPointsAsMatrix()};
  // Cubic has N_ = 4; requesting continuity order 4 needs derivatives past
  // what an order-3 curve can provide
  EXPECT_THROW(curve.applyContinuity(locked, {1.0, 1.0, 1.0, 1.0}), std::logic_error);
}

// Found while writing approx_test.cpp: fromPolyline with order = 0 uses
// polyline.size() - 1 as the order, so a realistic polyline (~29 points)
// dispatches an order-28 Levenberg-Marquardt fit that effectively never
// finishes (caught by the ctest TIMEOUT). Post-fix: automatic order
// selection completes in reasonable time on real polylines.
TEST(KnownBugsTest, DISABLED_FromPolylineAutoOrderHangsOnRealisticPolyline)
{
  Curve curve{curvePointsAsMatrix()};
  PointVector polyline = curve.polyline();
  Curve fitted = Curve::fromPolyline(polyline);
  EXPECT_GE(fitted.order(), 1u);
}

} // namespace Bezier
