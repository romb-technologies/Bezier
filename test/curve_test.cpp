#include "test_data.hpp"
#include "test_oracles.hpp"

#include <cmath>

#include <gtest/gtest.h>

#include "Bezier/bezier.h"
#include "Bezier/declarations.h"
#include "Bezier/utils.h"

namespace Bezier
{

class BezierTest : public ::testing::Test
{
public:
  BezierTest() : curve_{curvePointsAsMatrix()}, curve_roots_{rootPointsAsMatrix()} {}

protected:
  Curve curve_;
  Curve curve_roots_;
};

TEST_F(BezierTest, CurveOrderTest) { EXPECT_EQ(curve_.order(), 3) << "Curve order differs from expected."; }

TEST_F(BezierTest, CurveControlPointsTest)
{
  PointVector points_expected = curvePointsAsVector();
  PointVector points_curve = curve_.controlPoints();
  EXPECT_EQ(points_expected, points_curve) << "Curve control points differ from expected";
}

TEST_F(BezierTest, CurveControlPointTest)
{
  Point point_expected = curvePointsAsVector()[2];
  Point point_curve = curve_.controlPoint(2);
  EXPECT_EQ(point_expected, point_curve) << "Curve control point at idx differs from expected";
}

TEST_F(BezierTest, CurveEndPointsTest)
{
  std::pair<Point, Point> end_points_expected{curvePointsAsVector().front(), curvePointsAsVector().back()};
  std::pair<Point, Point> point_curve = curve_.endPoints();
  EXPECT_EQ(end_points_expected, point_curve) << "Curve end points differ from expected";
}

TEST_F(BezierTest, CurvePolylineFlatnessTest)
{
  // Finer flatness yields more vertices (correctness is CurvePolylineContractTest)
  PointVector polyline = curve_.polyline();
  EXPECT_LT(curve_.polyline(1.0).size(), polyline.size());
  EXPECT_GT(curve_.polyline(0.01).size(), polyline.size());
}

TEST_F(BezierTest, CurveLengthMonotonicTest)
{
  // length(t) starts at 0 and strictly increases; full length equals length(1)
  // (absolute value is checked against the chord-length oracle elsewhere)
  EXPECT_NEAR(curve_.length(0.0), 0.0, Oracles::kGeom);
  double prev = 0.0;
  for (double t : Oracles::sampleParams(Oracles::kCoarseSamples))
  {
    double len = curve_.length(t);
    if (t > 0.0)
      EXPECT_GT(len, prev) << "length not strictly increasing at t=" << t;
    prev = len;
  }
  EXPECT_DOUBLE_EQ(curve_.length(), curve_.length(1.0));
}

TEST_F(BezierTest, CurveSetControlPointTest)
{
  Point new_control_point{40, 40};
  EXPECT_NO_THROW(curve_.setControlPoint(2, new_control_point)) << "Curve set control point failed";
  EXPECT_EQ(curve_.controlPoint(2), new_control_point) << "Failed to set control point";
}

TEST_F(BezierTest, CurveValueAtMultipleParamsTest)
{
  // Batch overload must agree with the scalar one (de Casteljau oracle)
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  Eigen::MatrixX2d points = curve_.valueAt(t_vals);
  ASSERT_EQ(points.rows(), static_cast<Eigen::Index>(t_vals.size()));
  PointVector cp = curve_.controlPoints();
  for (size_t i = 0; i < t_vals.size(); i++)
  {
    Point expected = Oracles::deCasteljau(cp, t_vals[i]);
    EXPECT_NEAR((points.row(i).transpose() - expected).norm(), 0.0, Oracles::kGeom) << "t=" << t_vals[i];
  }
}

TEST_F(BezierTest, CurveCurvatureMatchesOracle)
{
  PointVector cp = curve_.controlPoints();
  for (double t : {0.0, 0.25, 0.5, 0.75, 1.0})
  {
    EXPECT_NEAR(curve_.curvatureAt(t), Oracles::curvature(cp, t), Oracles::kGeom) << "kappa at t=" << t;
    EXPECT_NEAR(curve_.curvatureDerivativeAt(t), Oracles::curvatureDerivative(cp, t), Oracles::kGeom)
        << "kappa' at t=" << t;
  }
}

TEST_F(BezierTest, CurveDerivativeTest)
{
  PointVector expected = TestData::toPointVector(TestData::kExpectedFirstDerivative);
  PointVector derivative_control_points = curve_.derivative().controlPoints();
  ASSERT_EQ(expected.size(), derivative_control_points.size());
  for (size_t i = 0; i < derivative_control_points.size(); i++)
  {
    EXPECT_EQ(expected[i], derivative_control_points[i]);
  }
}

TEST_F(BezierTest, CurveNthDerivativeTest)
{
  int N = 2;
  PointVector expected = TestData::toPointVector(TestData::kExpectedSecondDerivative);
  PointVector derivative_control_points = curve_.derivative(N).controlPoints();
  ASSERT_EQ(expected.size(), derivative_control_points.size());
  for (size_t i = 0; i < derivative_control_points.size(); i++)
  {
    EXPECT_EQ(expected[i], derivative_control_points[i]);
  }
}

TEST_F(BezierTest, CurveRootsTest)
{
  std::vector<double> curve_roots = curve_roots_.roots();

  for (double root : curve_roots)
  {
    EXPECT_GE(root, 0.0);
    EXPECT_LE(root, 1.0);

    // a root lies on one of the axes
    Point p = curve_roots_.valueAt(root);
    EXPECT_TRUE(std::abs(p.x()) < Oracles::kGeom || std::abs(p.y()) < Oracles::kGeom)
        << "root t=" << root << " not on an axis: (" << p.x() << ", " << p.y() << ")";
  }

  EXPECT_EQ(curve_roots.size(), 3) << "this fixture has 3 roots";
}

TEST_F(BezierTest, CurveExtremaTest)
{
  std::vector<double> curve_extrema = curve_roots_.extrema();

  for (double t : curve_extrema)
  {
    EXPECT_GE(t, 0.0);
    EXPECT_LE(t, 1.0);

    // one derivative component vanishes at an axis-aligned extremum
    Vector deriv = curve_roots_.derivativeAt(t);
    EXPECT_TRUE(std::abs(deriv.x()) < Oracles::kGeom || std::abs(deriv.y()) < Oracles::kGeom)
        << "extremum t=" << t << " derivative not axis-aligned: (" << deriv.x() << ", " << deriv.y() << ")";
  }

  EXPECT_EQ(curve_extrema.size(), 2) << "this fixture has 2 extrema";
}

TEST_F(BezierTest, CurveBoundingBoxTest)
{
  BoundingBox bbox = curve_.boundingBox();

  auto endpoints = curve_.endPoints();
  EXPECT_GE(endpoints.first.x(), bbox.min().x());
  EXPECT_LE(endpoints.first.x(), bbox.max().x());
  EXPECT_GE(endpoints.first.y(), bbox.min().y());
  EXPECT_LE(endpoints.first.y(), bbox.max().y());
  EXPECT_GE(endpoints.second.x(), bbox.min().x());
  EXPECT_LE(endpoints.second.x(), bbox.max().x());
  EXPECT_GE(endpoints.second.y(), bbox.min().y());
  EXPECT_LE(endpoints.second.y(), bbox.max().y());

  for (double t : curve_.extrema())
  {
    Point p = curve_.valueAt(t);
    EXPECT_GE(p.x(), bbox.min().x() - Oracles::kGeom) << "extremum outside bbox";
    EXPECT_LE(p.x(), bbox.max().x() + Oracles::kGeom) << "extremum outside bbox";
    EXPECT_GE(p.y(), bbox.min().y() - Oracles::kGeom) << "extremum outside bbox";
    EXPECT_LE(p.y(), bbox.max().y() + Oracles::kGeom) << "extremum outside bbox";
  }

  EXPECT_LT(bbox.min().x(), bbox.max().x());
  EXPECT_LT(bbox.min().y(), bbox.max().y());

  // every sampled point lies inside the box
  for (double t : Oracles::sampleParams(Oracles::kCoarseSamples))
  {
    Point p = curve_.valueAt(t);
    EXPECT_GE(p.x(), bbox.min().x() - Oracles::kGeom) << "Point at t=" << t << " outside bbox";
    EXPECT_LE(p.x(), bbox.max().x() + Oracles::kGeom) << "Point at t=" << t << " outside bbox";
    EXPECT_GE(p.y(), bbox.min().y() - Oracles::kGeom) << "Point at t=" << t << " outside bbox";
    EXPECT_LE(p.y(), bbox.max().y() + Oracles::kGeom) << "Point at t=" << t << " outside bbox";
  }
}

TEST_F(BezierTest, CurveIntersectionsTest)
{
  // Crossing curve of a different order (quartic) than the cubic fixture
  Curve curve_with_intersections{intersectionPointsAsMatrix()};
  PointVector intersections = curve_.intersections(curve_with_intersections);

  EXPECT_GE(intersections.size(), 1) << "Too few intersections found";
  EXPECT_LE(intersections.size(), 9) << "Too many intersections found";

  // Reported intersection points come from bounding-box subdivision, so they
  // carry more error than a closed-form result; intersections are reworked in
  // the follow-up PR, where this local bound gets revisited.
  constexpr double kIntersectTol = 1e-4;
  for (const Point& p : intersections)
  {
    EXPECT_NEAR(curve_.distance(p), 0.0, kIntersectTol);
    EXPECT_NEAR(curve_with_intersections.distance(p), 0.0, kIntersectTol);
  }
}

TEST_F(BezierTest, CurveSplitTest)
{
  std::pair<Curve, Curve> split_curves = curve_.splitCurve(0.5);

  // halves meet at the split point and keep the original order
  Point left_end = split_curves.first.endPoints().second;
  Point right_start = split_curves.second.endPoints().first;
  EXPECT_NEAR((left_end - right_start).norm(), 0.0, Oracles::kGeom);

  Point original_mid = curve_.valueAt(0.5);
  EXPECT_NEAR((left_end - original_mid).norm(), 0.0, Oracles::kGeom);

  EXPECT_EQ(split_curves.first.order(), curve_.order());
  EXPECT_EQ(split_curves.second.order(), curve_.order());

  auto orig_endpoints = curve_.endPoints();
  auto left_start = split_curves.first.endPoints().first;
  auto right_end = split_curves.second.endPoints().second;
  EXPECT_NEAR((left_start - orig_endpoints.first).norm(), 0.0, Oracles::kGeom);
  EXPECT_NEAR((right_end - orig_endpoints.second).norm(), 0.0, Oracles::kGeom);

  // each half reparametrizes the matching slice of the original
  auto p1 = split_curves.first.valueAt(0.0);
  auto p2 = curve_.valueAt(0.0);
  EXPECT_NEAR((p1 - p2).norm(), 0.0, Oracles::kGeom);

  p1 = split_curves.first.valueAt(0.5);
  p2 = curve_.valueAt(0.25);
  EXPECT_NEAR((p1 - p2).norm(), 0.0, Oracles::kGeom);

  p1 = split_curves.second.valueAt(0.5);
  p2 = curve_.valueAt(0.75);
  EXPECT_NEAR((p1 - p2).norm(), 0.0, Oracles::kGeom);

  p1 = split_curves.second.valueAt(1.0);
  p2 = curve_.valueAt(1.0);
  EXPECT_NEAR((p1 - p2).norm(), 0.0, Oracles::kGeom);
}

TEST_F(BezierTest, CurveProjectPointTest)
{
  Point point{100, 150};
  double t_projected = curve_.projectPoint(point);
  EXPECT_GE(t_projected, 0.0);
  EXPECT_LE(t_projected, 1.0);

  // No sampled point is closer than the projection
  double min_dist = (point - curve_.valueAt(t_projected)).norm();
  for (double t : Oracles::sampleParams(Oracles::kDenseSamples))
    EXPECT_GE((point - curve_.valueAt(t)).norm(), min_dist - Oracles::kGeom) << "closer point at t=" << t;
}

TEST_F(BezierTest, CurveDistanceTest)
{
  Point point{100, 150};
  double distance = curve_.distance(point);
  EXPECT_GE(distance, 0.0);

  // distance == ‖point − valueAt(projectPoint)‖
  double expected_distance = (point - curve_.valueAt(curve_.projectPoint(point))).norm();
  EXPECT_NEAR(distance, expected_distance, Oracles::kGeom);
}

// ---------------------------------------------------------------------------
// Oracle-based and property tests (refactor-proof; see test_oracles.hpp)
// ---------------------------------------------------------------------------

// Deterministic control points for an arbitrary-order test curve
static PointVector makeControlPoints(unsigned order)
{
  PointVector cp;
  cp.reserve(order + 1);
  for (unsigned i = 0; i <= order; i++)
    cp.emplace_back(100 * std::cos(1.7 * i + order), 100 * std::sin(2.3 * i) + 10 * i);
  return cp;
}

TEST_F(BezierTest, CurvePolylineContractTest)
{
  PointVector polyline = curve_.polyline();
  ParamVector params = curve_.polylineParams();

  ASSERT_EQ(polyline.size(), params.size()) << "polyline and polylineParams must correspond 1:1";
  EXPECT_DOUBLE_EQ(params.front(), 0.0);
  EXPECT_DOUBLE_EQ(params.back(), 1.0);

  for (size_t i = 0; i < params.size(); i++)
  {
    Point on_curve = curve_.valueAt(params[i]);
    EXPECT_NEAR((on_curve - polyline[i]).norm(), 0.0, Oracles::kGeom) << "Polyline vertex " << i << " not on curve";
  }

  // Dense samples of the curve stay within flatness of the polyline
  double flatness = curve_.boundingBox().diagonal().norm() / 1000;
  for (double t : Oracles::sampleParams(Oracles::kFlatnessSamples))
    EXPECT_LE(Utils::dist(polyline, curve_.valueAt(t)), flatness + Oracles::kGeom) << "Curve too far at t=" << t;
}

TEST(CurveOracleTests, ValueAtMatchesDeCasteljau)
{
  for (unsigned order = 1; order <= 6; order++)
  {
    PointVector cp = makeControlPoints(order);
    Curve curve{cp};
    for (double t : Oracles::sampleParams(Oracles::kCoarseSamples))
    {
      Point expected = Oracles::deCasteljau(cp, t);
      Point actual = curve.valueAt(t);
      EXPECT_NEAR((actual - expected).norm(), 0.0, Oracles::kGeom) << "order=" << order << " t=" << t;
    }
  }
}

TEST_F(BezierTest, CurveLengthMatchesChordLengthOracle)
{
  // Fixture cubic
  double oracle = Oracles::chordLength(curve_.controlPoints());
  EXPECT_NEAR(curve_.length(), oracle, Oracles::kGeom * oracle);

  // Partial lengths
  EXPECT_NEAR(curve_.length(0.3), Oracles::chordLength(curve_.controlPoints(), 0.0, 0.3), Oracles::kGeom * oracle);
  EXPECT_NEAR(curve_.length(0.2, 0.8), Oracles::chordLength(curve_.controlPoints(), 0.2, 0.8), Oracles::kGeom * oracle);

  // Straight line: length equals the chord exactly
  Curve line{PointVector{{0, 0}, {3, 4}}};
  EXPECT_NEAR(line.length(), 5.0, Oracles::kGeom);

  // Symmetric quadratic — regression guard for the master-branch Chebyshev
  // truncation bug (even coefficients vanish); v4 must get this right
  Curve symmetric{PointVector{{0, 0}, {1, 2}, {2, 0}}};
  double symmetric_oracle = Oracles::chordLength(symmetric.controlPoints());
  EXPECT_NEAR(symmetric.length(), symmetric_oracle, Oracles::kGeom * symmetric_oracle);
}

TEST_F(BezierTest, CurveStepLengthRoundTripTest)
{
  for (double t : {0.1, 0.42, 0.7})
    for (double ds : {5.0, 35.0, -5.0, -15.0})
    {
      double new_t = curve_.step(t, ds);
      EXPECT_NEAR(curve_.length(t, new_t), ds, Oracles::kGeom) << "Round trip failed for t=" << t << " ds=" << ds;
    }

  // Out-of-range distances clamp to the curve ends
  EXPECT_DOUBLE_EQ(curve_.step(0.5, 1e6), 1.0);
  EXPECT_DOUBLE_EQ(curve_.step(0.5, -1e6), 0.0);
}

TEST_F(BezierTest, CurveReverseTest)
{
  Curve reversed{curve_};
  reversed.reverse();
  for (double t : Oracles::sampleParams(Oracles::kCoarseSamples))
  {
    Point expected = curve_.valueAt(1.0 - t);
    Point actual = reversed.valueAt(t);
    EXPECT_NEAR((actual - expected).norm(), 0.0, Oracles::kGeom) << "reverse mismatch at t=" << t;
  }
}

TEST_F(BezierTest, CurveSplitMultiTest)
{
  ParamVector t_split{0.25, 0.6};
  std::vector<Curve> pieces = curve_.splitCurve(t_split);
  ASSERT_EQ(pieces.size(), 3u);

  // Pieces map linearly onto the original parameter ranges
  std::array<std::pair<double, double>, 3> ranges{{{0.0, 0.25}, {0.25, 0.6}, {0.6, 1.0}}};
  for (size_t k = 0; k < pieces.size(); k++)
    for (double s : Oracles::sampleParams(Oracles::kCoarseSamples))
    {
      Point expected = curve_.valueAt(ranges[k].first + s * (ranges[k].second - ranges[k].first));
      Point actual = pieces[k].valueAt(s);
      EXPECT_NEAR((actual - expected).norm(), 0.0, Oracles::kGeom) << "piece " << k << " s=" << s;
    }

  // C0 continuity between adjacent pieces
  for (size_t k = 0; k + 1 < pieces.size(); k++)
  {
    Point left = pieces[k].endPoints().second;
    Point right = pieces[k + 1].endPoints().first;
    EXPECT_NEAR((left - right).norm(), 0.0, Oracles::kGeom);
  }

  // Unsorted input produces the same pieces
  std::vector<Curve> pieces_unsorted = curve_.splitCurve(ParamVector{0.6, 0.25});
  ASSERT_EQ(pieces_unsorted.size(), 3u);
  for (size_t k = 0; k < pieces.size(); k++)
    EXPECT_EQ(pieces[k].controlPoints(), pieces_unsorted[k].controlPoints());

  // Empty parameter vector returns the whole curve as a single piece
  std::vector<Curve> whole = curve_.splitCurve(ParamVector{});
  ASSERT_EQ(whole.size(), 1u);
  EXPECT_EQ(whole.front().controlPoints(), curve_.controlPoints());
}

TEST_F(BezierTest, CurveRaiseOrderPreservesShape)
{
  Curve raised{curve_};
  raised.raiseOrder();
  EXPECT_EQ(raised.order(), curve_.order() + 1);
  for (double t : Oracles::sampleParams(Oracles::kCoarseSamples))
  {
    Point expected = curve_.valueAt(t);
    Point actual = raised.valueAt(t);
    EXPECT_NEAR((actual - expected).norm(), 0.0, Oracles::kGeom) << "raiseOrder changed shape at t=" << t;
  }

  // Lowering the raised curve recovers the original control points
  raised.lowerOrder();
  ASSERT_EQ(raised.order(), curve_.order());
  PointVector original = curve_.controlPoints();
  PointVector round_trip = raised.controlPoints();
  for (size_t i = 0; i < original.size(); i++)
  {
    EXPECT_NEAR((round_trip[i] - original[i]).norm(), 0.0, Oracles::kGeom);
  }
}

TEST(CurveOrderTests, LowerOrderThrowsOnFirstOrderCurve)
{
  Curve line{PointVector{{0, 0}, {1, 1}}};
  EXPECT_THROW(line.lowerOrder(), std::logic_error);
}

TEST_F(BezierTest, CurveDerivativeAtMatchesHodograph)
{
  // Compare against the exact derivative (hodograph evaluated by de Casteljau),
  // not a finite difference: closed-form vs closed-form holds to machine precision.
  PointVector cp = curve_.controlPoints();
  PointVector hodograph = Oracles::derivativeControlPoints(cp);
  for (double t : {0.1, 0.3, 0.5, 0.7, 0.9})
  {
    Vector expected = Oracles::deCasteljau(hodograph, t);
    Vector actual = curve_.derivativeAt(t);
    EXPECT_NEAR((actual - expected).norm(), 0.0, Oracles::kGeom) << "derivative mismatch at t=" << t;
  }
}

TEST_F(BezierTest, CurveTangentNormalPropertiesTest)
{
  for (double t : {0.0, 0.2, 0.5, 0.8, 1.0})
  {
    Vector tangent = curve_.tangentAt(t);
    Vector normal = curve_.normalAt(t);

    EXPECT_NEAR(tangent.norm(), 1.0, Oracles::kGeom) << "tangent not unit at t=" << t;
    EXPECT_NEAR(normal.norm(), 1.0, Oracles::kGeom) << "normal not unit at t=" << t;
    EXPECT_NEAR(tangent.dot(normal), 0.0, Oracles::kGeom) << "tangent/normal not orthogonal at t=" << t;

    // normal is the tangent rotated +90 degrees
    EXPECT_NEAR((normal - Vector(-tangent.y(), tangent.x())).norm(), 0.0, Oracles::kGeom);

    // tangent points along the derivative (same direction, not just parallel)
    Vector derivative = curve_.derivativeAt(t);
    EXPECT_NEAR(tangent.x() * derivative.y() - tangent.y() * derivative.x(), 0.0, Oracles::kGeom * derivative.norm());
    EXPECT_GT(tangent.dot(derivative), 0.0) << "tangent points against the curve at t=" << t;
  }
}

TEST_F(BezierTest, CurveProjectPointRecoveryTest)
{
  // Projecting a point that lies on the curve recovers its parameter
  for (double t : {0.1, 0.3, 0.5, 0.7, 0.9})
  {
    double projected = curve_.projectPoint(curve_.valueAt(t));
    EXPECT_NEAR(projected, t, Oracles::kGeom) << "did not recover parameter t=" << t;
  }

  // Points beyond the endpoints project to the endpoints
  Point before_start = curve_.valueAt(0.0) - 10 * curve_.tangentAt(0.0);
  Point past_end = curve_.valueAt(1.0) + 10 * curve_.tangentAt(1.0);
  EXPECT_DOUBLE_EQ(curve_.projectPoint(before_start), 0.0);
  EXPECT_DOUBLE_EQ(curve_.projectPoint(past_end), 1.0);
}

TEST_F(BezierTest, CurveApplyContinuityC1Test)
{
  Curve continued{curve_roots_};
  continued.applyContinuity(curve_, {1.0});

  // C1: position and first derivative match the locked curve's end
  Point pos_expected = curve_.valueAt(1.0);
  Point pos_actual = continued.valueAt(0.0);
  EXPECT_NEAR((pos_actual - pos_expected).norm(), 0.0, Oracles::kGeom);

  Vector der_expected = curve_.derivativeAt(1.0);
  Vector der_actual = continued.derivativeAt(0.0);
  EXPECT_NEAR((der_actual - der_expected).norm(), 0.0, Oracles::kGeom);
}

TEST_F(BezierTest, CurveApplyContinuityG1Test)
{
  constexpr double beta = 2.5;
  Curve continued{curve_roots_};
  continued.applyContinuity(curve_, {beta});

  // G1: position matches; first derivative is scaled by beta
  Point pos_expected = curve_.valueAt(1.0);
  Point pos_actual = continued.valueAt(0.0);
  EXPECT_NEAR((pos_actual - pos_expected).norm(), 0.0, Oracles::kGeom);

  Vector der_expected = beta * curve_.derivativeAt(1.0);
  Vector der_actual = continued.derivativeAt(0.0);
  EXPECT_NEAR((der_actual - der_expected).norm(), 0.0, Oracles::kGeom);
}

TEST_F(BezierTest, CurveCopySemanticsTest)
{
  double original_length = curve_.length();
  PointVector original_points = curve_.controlPoints();

  Curve copy{curve_};
  EXPECT_EQ(copy.controlPoints(), original_points);

  // Mutating the copy must not affect the original or its cached data
  copy.setControlPoint(1, Point{0, 0});
  EXPECT_EQ(curve_.controlPoints(), original_points);
  EXPECT_DOUBLE_EQ(curve_.length(), original_length);
  EXPECT_NE(copy.length(), original_length);

  Curve assigned{curve_roots_};
  assigned = curve_;
  EXPECT_EQ(assigned.controlPoints(), original_points);
  assigned.setControlPoint(0, Point{-1, -1});
  EXPECT_EQ(curve_.controlPoints(), original_points);
}

TEST_F(BezierTest, CurveLengthReversedArgsContractTest)
{
  // Pinned contract: Curve::length(t1, t2) == length(t2) - length(t1),
  // so swapped arguments yield a negative value. A future contract change
  // (e.g. throwing or taking the absolute value) must edit this test.
  EXPECT_DOUBLE_EQ(curve_.length(0.8, 0.2), -curve_.length(0.2, 0.8));
  EXPECT_LT(curve_.length(0.8, 0.2), 0.0);
}

// ---------------------------------------------------------------------------
// Degenerate, cusp and contract edge cases
// ---------------------------------------------------------------------------

TEST_F(BezierTest, CurveIntersectionsSharedEndpointTest)
{
  // Only transversal crossings are reported; a contact exactly at a shared
  // endpoint is not returned (see Curve::intersections).
  Point shared = curve_.endPoints().second;
  Curve second{PointVector{shared,
                           {shared.x() + 50, shared.y() + 80},
                           {shared.x() + 120, shared.y() + 10},
                           {shared.x() + 150, shared.y() + 90}}};
  for (const Point& p : curve_.intersections(second))
    EXPECT_GT((p - shared).norm(), Oracles::kGeom) << "Endpoint contact should not be reported";
}

TEST_F(BezierTest, CurveApplyContinuityRaisesOrderTest)
{
  // A continuity order >= the curve order raises the curve's own order to make
  // room; derivatives beyond the locked curve's degree are zero vectors.
  Curve continued{curve_roots_};
  ASSERT_NO_THROW(continued.applyContinuity(curve_, {1.0, 1.0, 1.0, 1.0}));
  EXPECT_GE(continued.order(), 4u);
  EXPECT_TRUE(continued.valueAt(0.0).isApprox(curve_.valueAt(1.0), Oracles::kGeom));
}

TEST(CurveLengthTests, DegenerateCurveHasZeroLength)
{
  Curve degenerate{PointVector{{10, 10}, {10, 10}, {10, 10}, {10, 10}}};
  double length{};
  ASSERT_NO_THROW(length = degenerate.length());
  EXPECT_NEAR(length, 0.0, Oracles::kGeom);
}

TEST(CurveCuspTests, StepReturnsFiniteParameter)
{
  // Genuine cusp at t = 0.5 (velocity is (300(1-2t)^2, 300(1-2t)))
  Curve cusped{PointVector{{0, 0}, {100, 100}, {0, 100}, {100, 0}}};
  double t = cusped.step(0.5, 10.0);
  EXPECT_TRUE(std::isfinite(t));
  EXPECT_GE(t, 0.0);
  EXPECT_LE(t, 1.0);
}

TEST(CurveCuspTests, NormalIsFiniteUnitVector)
{
  // At a cusp the direction comes from the first non-vanishing derivative.
  Curve cusped{PointVector{{0, 0}, {100, 100}, {0, 100}, {100, 0}}};
  Vector normal = cusped.normalAt(0.5);
  ASSERT_TRUE(std::isfinite(normal.x()));
  ASSERT_TRUE(std::isfinite(normal.y()));
  EXPECT_NEAR(normal.norm(), 1.0, Oracles::kGeom);
}

} // namespace Bezier
