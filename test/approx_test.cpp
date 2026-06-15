#include "test_data.hpp"
#include "test_oracles.hpp"

#include <gtest/gtest.h>

#include "Bezier/bezier.h"
#include "Bezier/utils.h"

namespace Bezier
{

// Loose property tests for the approximation methods. They assert the
// contract (fit quality relative to the bounding-box diagonal), not the
// current Levenberg-Marquardt implementation, so they survive the planned
// rewrite of these methods.

class ApproxTest : public ::testing::Test
{
public:
  ApproxTest() : curve_{curvePointsAsMatrix()}, fit_tolerance_{Oracles::kFit * curve_.boundingBox().diagonal().norm()}
  {
  }

protected:
  Curve curve_;
  double fit_tolerance_;
};

TEST_F(ApproxTest, FromPolylineFitsSourcePolyline)
{
  PointVector polyline = curve_.polyline();
  Curve fitted = Curve::fromPolyline(polyline, curve_.order());
  EXPECT_EQ(fitted.order(), curve_.order());

  // Sampled fitted points stay close to the source polyline
  for (double t{}; t <= 1.0; t += 0.02)
    EXPECT_LE(Utils::dist(polyline, fitted.valueAt(t)), fit_tolerance_) << "fit too far at t=" << t;

  // First and last polyline points are preserved
  EXPECT_LE((fitted.valueAt(0.0) - polyline.front()).norm(), fit_tolerance_);
  EXPECT_LE((fitted.valueAt(1.0) - polyline.back()).norm(), fit_tolerance_);
}

TEST_F(ApproxTest, FromPolylineEdgeCases)
{
  // order = 0 -> automatic order selection returns a usable curve.
  PointVector short_polyline = Utils::polylineSimplified(curve_.polyline(), 6);
  Curve automatic = Curve::fromPolyline(short_polyline);
  EXPECT_GE(automatic.order(), 1u);
  for (double t{}; t <= 1.0; t += 0.1)
    EXPECT_LE(Utils::dist(short_polyline, automatic.valueAt(t)), fit_tolerance_);

  // Fewer than two points is an error
  EXPECT_THROW(Curve::fromPolyline(PointVector{}), std::logic_error);
  EXPECT_THROW(Curve::fromPolyline(PointVector{{1, 1}}), std::logic_error);

  // Two points produce the order-1 curve through both
  Curve segment = Curve::fromPolyline(PointVector{{1, 2}, {5, 6}});
  EXPECT_EQ(segment.order(), 1u);
  EXPECT_EQ(segment.endPoints().first, Point(1, 2));
  EXPECT_EQ(segment.endPoints().second, Point(5, 6));
}

TEST_F(ApproxTest, FromPolylineAutoOrderOnRealisticPolyline)
{
  // Guards the auto-order cap: full-resolution input must not hang (FromPolyline
  // EdgeCases only exercises a 6-point polyline, too small to trigger it).
  Curve fitted = Curve::fromPolyline(curve_.polyline());
  EXPECT_GE(fitted.order(), 1u);
}

TEST_F(ApproxTest, OffsetCurveStraightLine)
{
  // Analytic case: offsetting a straight line shifts it by exactly |d|
  Curve line{PointVector{{0, 0}, {100, 0}}};
  for (double offset : {5.0, -5.0})
  {
    Curve offset_curve = Curve::offsetCurve(line, offset);
    for (double t{}; t <= 1.0; t += 0.05)
      EXPECT_NEAR(line.distance(offset_curve.valueAt(t)), std::fabs(offset), 1e-6)
          << "offset=" << offset << " t=" << t;
  }
}

TEST_F(ApproxTest, OffsetCurveGentleArc)
{
  // Gentle arc with small offset: sampled distance to the source curve stays
  // within a few percent of |d|
  Curve arc{PointVector{{0, 0}, {50, 10}, {100, 0}}};
  const double offset{2.0};
  Curve offset_curve = Curve::offsetCurve(arc, offset);
  for (double t{}; t <= 1.0; t += 0.05)
    EXPECT_NEAR(arc.distance(offset_curve.valueAt(t)), offset, 0.05 * offset) << "t=" << t;
}

TEST_F(ApproxTest, OffsetCurveRespectsRequestedOrder)
{
  EXPECT_EQ(Curve::offsetCurve(curve_, 5.0, 4).order(), 4u);
  // Default order is source order + 1
  EXPECT_EQ(Curve::offsetCurve(curve_, 5.0).order(), curve_.order() + 1);
}

TEST_F(ApproxTest, JoinCurvesOrderOneIsExactLine)
{
  Curve other{rootPointsAsMatrix()};
  Curve joined = Curve::joinCurves(curve_, other, 1);
  EXPECT_EQ(joined.order(), 1u);
  EXPECT_EQ(joined.endPoints().first, curve_.endPoints().first);
  EXPECT_EQ(joined.endPoints().second, other.endPoints().second);
}

TEST_F(ApproxTest, JoinCurvesSmoothPair)
{
  // Two C0-chained cubics (translated copies, same trick as the PolyCurve fixture)
  Eigen::MatrixX2d cp = curvePointsAsMatrix();
  Eigen::RowVector2d delta = cp.row(cp.rows() - 1) - cp.row(0);
  Curve first{cp};
  Curve second{Eigen::MatrixX2d(cp.rowwise() + delta)};

  Curve joined = Curve::joinCurves(first, second);

  // Endpoints of the result match the outer endpoints of the pair
  EXPECT_LE((joined.endPoints().first - first.endPoints().first).norm(), fit_tolerance_);
  EXPECT_LE((joined.endPoints().second - second.endPoints().second).norm(), fit_tolerance_);

  // Sampled result stays near the original pair
  for (double t{}; t <= 1.0; t += 0.02)
  {
    Point p = joined.valueAt(t);
    EXPECT_LE(std::min(first.distance(p), second.distance(p)), fit_tolerance_) << "joined curve too far at t=" << t;
  }
}

} // namespace Bezier
