#include "test_data.hpp"
#include "test_oracles.hpp"

#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>

#include "Bezier/bezier.h"
#include "Bezier/utils.h"

namespace Bezier
{

// Property tests for the approximation methods. Parts that are exact by
// construction — interpolated endpoints (fromPolyline pins t=0 and t=1), orders,
// and the analytic straight-line offset — are asserted at machine precision. The
// genuinely approximate interior fit of a curved shape is bounded by the input
// polyline's own resolution or the measured implementation floor, never a bare
// percentage of the bounding box. The measured floors are pinned to the current
// Levenberg-Marquardt implementation and get revisited when it is reworked.

class ApproxTest : public ::testing::Test
{
public:
  ApproxTest() : curve_{curvePointsAsMatrix()}, polyline_flatness_{curve_.boundingBox().diagonal().norm() / 1000}
  {
  }

protected:
  Curve curve_;
  double polyline_flatness_; // resolution of curve_.polyline(): default flatness is 0.1% of the bbox diagonal
};

TEST_F(ApproxTest, FromPolylineFitsSourcePolyline)
{
  PointVector polyline = curve_.polyline();
  Curve fitted = Curve::fromPolyline(polyline, curve_.order());
  EXPECT_EQ(fitted.order(), curve_.order());

  // A same-order refit tracks the source polyline to within its own resolution
  // (the fit cannot be asked to beat its input); measured worst ~0.9x flatness.
  constexpr double kFitFactor = 1.5;
  for (double t : Oracles::sampleParams(Oracles::kDenseSamples))
    EXPECT_LE(Utils::dist(polyline, fitted.valueAt(t)), kFitFactor * polyline_flatness_) << "fit too far at t=" << t;

  // Endpoints are interpolated (t=0/1 pinned) -> match to solver precision (~1e-11)
  EXPECT_NEAR((fitted.valueAt(0.0) - polyline.front()).norm(), 0.0, Oracles::kGeom);
  EXPECT_NEAR((fitted.valueAt(1.0) - polyline.back()).norm(), 0.0, Oracles::kGeom);
}

TEST_F(ApproxTest, FromPolylineEdgeCases)
{
  // order = 0 -> automatic order selection returns a usable curve.
  // 6 most shape-significant points of the curve's polyline (short enough not to hit the auto-order cap).
  PointVector dense = curve_.polyline();
  std::vector<unsigned> vw = Utils::visvalingamWhyatt(dense);
  std::vector<unsigned> idx(vw.begin(), vw.begin() + 6);
  std::sort(idx.begin(), idx.end());
  PointVector short_polyline;
  for (unsigned i : idx)
    short_polyline.push_back(dense[i]);
  Curve automatic = Curve::fromPolyline(short_polyline);
  EXPECT_GE(automatic.order(), 1u);

  // The auto fit follows the deliberately coarse 6-point simplified polyline; the
  // bound is the measured floor (~2.0) for this low-resolution input.
  constexpr double kAutoFitTol = 3.0;
  for (double t : Oracles::sampleParams(Oracles::kCoarseSamples))
    EXPECT_LE(Utils::dist(short_polyline, automatic.valueAt(t)), kAutoFitTol);

  // Fewer than two points is an error
  EXPECT_THROW(Curve::fromPolyline(PointVector{}), std::logic_error);
  EXPECT_THROW(Curve::fromPolyline(PointVector{{1, 1}}), std::logic_error);

  // Two points produce the order-1 curve through both (exact copy)
  Curve segment = Curve::fromPolyline(PointVector{{1, 2}, {5, 6}});
  EXPECT_EQ(segment.order(), 1u);
  EXPECT_EQ(segment.endPoints().first, Point(1, 2));
  EXPECT_EQ(segment.endPoints().second, Point(5, 6));
}

TEST_F(ApproxTest, FromPolylineAutoOrderOnRealisticPolyline)
{
  // Auto-order fromPolyline on a full polyline completes and returns a usable curve.
  Curve fitted = Curve::fromPolyline(curve_.polyline());
  EXPECT_GE(fitted.order(), 1u);
}

TEST_F(ApproxTest, OffsetCurveStraightLine)
{
  // Analytic case: offsetting a line yields an exact parallel segment (fromPolyline
  // of two points), so the distance to the source is exactly |d| to projection precision.
  Curve line{PointVector{{0, 0}, {100, 0}}};
  for (double offset : {5.0, -5.0})
  {
    Curve offset_curve = Curve::offsetCurve(line, offset);
    for (double t : Oracles::sampleParams(Oracles::kCoarseSamples))
      EXPECT_NEAR(line.distance(offset_curve.valueAt(t)), std::fabs(offset), Oracles::kGeom)
          << "offset=" << offset << " t=" << t;
  }
}

TEST_F(ApproxTest, OffsetCurveGentleArc)
{
  // Offsetting a curved arc by a constant normal distance is only an approximate
  // locus (curvature bends the true offset off the constant-distance curve), plus
  // fit error. Bound is the measured deviation (~0.056, ~3% of |d|) of the impl.
  constexpr double kOffsetArcTol = 0.07;
  Curve arc{PointVector{{0, 0}, {50, 10}, {100, 0}}};
  double offset{2.0};
  Curve offset_curve = Curve::offsetCurve(arc, offset);
  for (double t : Oracles::sampleParams(Oracles::kCoarseSamples))
    EXPECT_NEAR(arc.distance(offset_curve.valueAt(t)), offset, kOffsetArcTol) << "t=" << t;
}

TEST_F(ApproxTest, OffsetCurveRespectsRequestedOrder)
{
  // An explicitly requested order is honoured exactly.
  EXPECT_EQ(Curve::offsetCurve(curve_, 5.0, 4).order(), 4u);
  // Default order = 0 delegates to fromPolyline's automatic order selection.
  EXPECT_GE(Curve::offsetCurve(curve_, 5.0).order(), 1u);
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

  // Endpoints are interpolated by fromPolyline (t=0/1 pinned); the default order-6
  // solve is ill-conditioned, so they match only to the measured Bernstein/
  // Vandermonde floor (~1.3e-6) -- still 100x tighter than the old 1e-4.
  constexpr double kApproxEndpoint = 1e-5;
  EXPECT_NEAR((joined.endPoints().first - first.endPoints().first).norm(), 0.0, kApproxEndpoint);
  EXPECT_NEAR((joined.endPoints().second - second.endPoints().second).norm(), 0.0, kApproxEndpoint);

  // The joined order-6 curve approximates the two-cubic chain to a measured floor (~2.2).
  constexpr double kJoinFitTol = 3.0;
  for (double t : Oracles::sampleParams(Oracles::kDenseSamples))
  {
    Point p = joined.valueAt(t);
    EXPECT_LE(std::min(first.distance(p), second.distance(p)), kJoinFitTol) << "joined curve too far at t=" << t;
  }
}

TEST_F(ApproxTest, FromPolylineOverOrderStaysBounded)
{
  // The ridge regularisation must bound the over-order (Runge) blow-up: fitting a cubic's polyline at
  // orders well above 3 stays near the data instead of exploding (control points -> 1e3..1e9, curve
  // oscillating between the samples). Without the ridge these deviations are 1e2..1e9% of the diagonal.
  PointVector polyline = curve_.polyline();
  // The ridge keeps the over-order fit within the input polyline's own resolution
  // (measured worst ~0.9x flatness, same as the same-order refit) instead of exploding.
  constexpr double kFitFactor = 1.5;
  for (unsigned order : {8u, 9u, 10u})
  {
    Curve fitted = Curve::fromPolyline(polyline, order);
    for (double t : Oracles::sampleParams(Oracles::kDenseSamples))
    {
      Point p = fitted.valueAt(t);
      ASSERT_TRUE(std::isfinite(p.x()) && std::isfinite(p.y())) << "non-finite at order=" << order << " t=" << t;
      EXPECT_LE(Utils::dist(polyline, p), kFitFactor * polyline_flatness_)
          << "over-order fit strayed at order=" << order << " t=" << t;
    }
  }
}

TEST_F(ApproxTest, FromPolylineAutoOrderRecoversParsimoniousOrder)
{
  // order = 0 must recover a parsimonious order on cleanly representable input, not the cap or 1.
  // Auto-order recovers each cubic well; the bound is a measured floor (worst ~0.12)
  // covering both fixtures -- their bbox diagonals differ, so a single
  // polyline_flatness_ does not apply.
  constexpr double kAutoOrderTol = 0.2;
  for (const auto& cp : {curvePointsAsMatrix(), rootPointsAsMatrix()})
  {
    Curve source{cp}; // both fixtures are cubics
    PointVector polyline = source.polyline();
    Curve fitted = Curve::fromPolyline(polyline); // order = 0 -> automatic selection
    EXPECT_GE(fitted.order(), 2u);
    EXPECT_LE(fitted.order(), source.order() + 1) << "auto-order did not stay parsimonious";
    for (double t : Oracles::sampleParams(Oracles::kDenseSamples))
      EXPECT_LE(Utils::dist(polyline, fitted.valueAt(t)), kAutoOrderTol);
  }
}

TEST_F(ApproxTest, FromPolylineDenseInputFixedOrder)
{
  // A fixed-order fit on dense input is trimmed to ~3N points internally, so it stays accurate and
  // keeps the requested order without the optimiser's variable count (and cost) scaling with M.
  // The suite-wide 60s ctest TIMEOUT is the "does not hang / no O(M^3) blow-up" guard.
  PointVector dense = curve_.polyline(curve_.boundingBox().diagonal().norm() / 100000);
  ASSERT_GT(dense.size(), 100u) << "expected a dense polyline (M >> 3N)";
  Curve fitted = Curve::fromPolyline(dense, 3);
  EXPECT_EQ(fitted.order(), 3u);
  // Input is ~exact (100x finer than default), so an order-3 fit of a cubic is bounded by
  // the solver floor, not input resolution (measured worst ~6e-3).
  constexpr double kDenseFitTol = 1e-2;
  for (double t : Oracles::sampleParams(Oracles::kDenseSamples))
    EXPECT_LE(Utils::dist(dense, fitted.valueAt(t)), kDenseFitTol);
}

TEST_F(ApproxTest, FromPolylineDegenerateInputsStayFinite)
{
  // Order exceeding the available points clamps to a valid lower order without throwing.
  Curve clamped = Curve::fromPolyline(curvePointsAsVector(), 8); // 4 points, order 8 requested
  EXPECT_LE(clamped.order(), 3u);

  // Duplicate-consecutive and fully-collinear polylines must not produce NaNs (the 1e-12 init floors
  // on the centripetal log-gaps guard against zero-length segments).
  PointVector duplicates{{0, 0}, {0, 0}, {50, 50}, {50, 50}, {100, 100}, {100, 100}};
  PointVector collinear{{0, 0}, {25, 0}, {50, 0}, {75, 0}, {100, 0}};
  for (const PointVector& poly : {duplicates, collinear})
  {
    Curve fitted = Curve::fromPolyline(poly, 3);
    for (const Point& cp : fitted.controlPoints())
      EXPECT_TRUE(std::isfinite(cp.x()) && std::isfinite(cp.y()));
  }
}

} // namespace Bezier
