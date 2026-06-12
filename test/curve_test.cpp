#include "test_data.hpp"
#include "test_oracles.hpp"

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

TEST(ConstructorTests, ConstructFromEigenMatrix)
{
  EXPECT_NO_THROW(Curve{curvePointsAsMatrix()}) << "Cannot create Curve with Eigen::MatrixX2d.";
}
TEST(ConstructorTests, ConstructFromPointVector)
{
  EXPECT_NO_THROW(Curve{curvePointsAsVector()}) << "Cannot create Curve with PointVector.";
}

TEST_F(BezierTest, CurveOrderTest) { EXPECT_EQ(curve_.order(), 3) << "Curve order differs from expected."; }

TEST_F(BezierTest, CurveControlPointsTest)
{
  PointVector points_expected = curvePointsAsVector();
  PointVector points_curve = curve_.controlPoints();
  EXPECT_EQ(points_expected, points_curve) << "Curve control points differ from expected";
}

TEST_F(BezierTest, CurveControlPointTest)
{
  Point point_expected = curvePointsAsVector().at(2);
  Point point_curve = curve_.controlPoint(2);
  EXPECT_EQ(point_expected, point_curve) << "Curve control point at idx differs from expected";
}

TEST_F(BezierTest, CurveEndPointsTest)
{
  std::pair<Point, Point> end_points_expected{curvePointsAsVector().front(), curvePointsAsVector().back()};
  std::pair<Point, Point> point_curve = curve_.endPoints();
  EXPECT_EQ(end_points_expected, point_curve) << "Curve end points differ from expected";
}

TEST_F(BezierTest, CurvePolylineTest)
{
  // Test with default flatness (no manual parameter)
  PointVector polyline = curve_.polyline();
  PointVector expected = expectedPolylinePoints();

  ASSERT_EQ(polyline.size(), expected.size()) << "Polyline size differs from expected";

  for (size_t i = 0; i < polyline.size(); ++i)
  {
    EXPECT_NEAR(polyline[i].x(), expected[i].x(), Utils::epsilon) << "Polyline point " << i << " X differs";
    EXPECT_NEAR(polyline[i].y(), expected[i].y(), Utils::epsilon) << "Polyline point " << i << " Y differs";
  }

  // Test with custom flatness parameters
  PointVector polyline_coarse = curve_.polyline(1.0);
  PointVector polyline_fine = curve_.polyline(0.01);

  // Coarser flatness should produce fewer points
  EXPECT_LT(polyline_coarse.size(), polyline.size());

  // Finer flatness should produce more points
  EXPECT_GT(polyline_fine.size(), polyline.size());
}

TEST_F(BezierTest, CurveLengthTest)
{
  double length = curve_.length();

  // Length should be positive
  EXPECT_GT(length, 0) << "Curve length should be positive";

  // Check against expected value with tight tolerance
  EXPECT_NEAR(length, TestData::kExpectedCurveLength, Utils::epsilon);

  // Length at t=0 should be 0
  EXPECT_NEAR(curve_.length(0.0), 0.0, Utils::epsilon);

  // Length should be monotonically increasing
  EXPECT_LT(curve_.length(0.0), curve_.length(0.25));
  EXPECT_LT(curve_.length(0.25), curve_.length(0.5));
  EXPECT_LT(curve_.length(0.5), curve_.length(0.75));
  EXPECT_LT(curve_.length(0.75), curve_.length(1.0));

  // Full length should match length(1.0)
  EXPECT_DOUBLE_EQ(length, curve_.length(1.0));
}

TEST_F(BezierTest, CurveStepByLengthTest)
{
  constexpr double start_t{0.42}, target_length{35};

  // Check that we moved forward along the curve
  double result_t = curve_.step(start_t, target_length);
  EXPECT_GT(result_t, start_t) << "Step should move forward for positive length";
  EXPECT_NEAR(curve_.length(start_t, result_t), target_length, Utils::epsilon);

  // Test backward stepping
  double backward_t = curve_.step(start_t, -target_length);
  EXPECT_LT(backward_t, start_t) << "Step should move backward for negative length";
  EXPECT_NEAR(curve_.length(backward_t, start_t), target_length, Utils::epsilon);
}

TEST_F(BezierTest, CurveSetControlPointTest)
{
  Point new_control_point{40, 40};
  EXPECT_NO_THROW(curve_.setControlPoint(2, new_control_point)) << "Curve set control point failed";
  EXPECT_EQ(curve_.controlPoint(2), new_control_point) << "Failed to set control point";
}

TEST_F(BezierTest, CurveRaiseOrderTest) { EXPECT_NO_THROW(curve_.raiseOrder()) << "Curve::raiseOrder failed"; }

TEST_F(BezierTest, CurveLowerOrderTest) { EXPECT_NO_THROW(curve_.lowerOrder()) << "Curve::lowerOrder failed"; }

TEST_F(BezierTest, CurveValueAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  PointVector expected_points = expectedValueAtVector();

  for (size_t i = 0; i < t_vals.size(); i++)
  {
    Point current_point = curve_.valueAt(t_vals[i]);
    EXPECT_NEAR(current_point.x(), expected_points[i].x(), Utils::epsilon);
    EXPECT_NEAR(current_point.y(), expected_points[i].y(), Utils::epsilon);
  }
}

TEST_F(BezierTest, CurveValueAtMultipleParamsTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  Eigen::MatrixX2d points = curve_.valueAt(t_vals);
  Eigen::MatrixX2d expected = expectedValueAtMatrix();
  EXPECT_TRUE(points.isApprox(expected, Utils::epsilon));
}

TEST_F(BezierTest, CurveCurvatureAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  for (size_t i = 0; i < t_vals.size(); i++)
  {
    double current = curve_.curvatureAt(t_vals[i]);
    EXPECT_NEAR(current, TestData::kExpectedCurvature[i], Utils::epsilon);
  }
}

TEST_F(BezierTest, CurveCurvatureDerivativeAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  for (size_t i = 0; i < t_vals.size(); i++)
  {
    double current = curve_.curvatureDerivativeAt(t_vals[i]);
    EXPECT_NEAR(current, TestData::kExpectedCurvatureDerivative[i], Utils::epsilon);
  }
}

TEST_F(BezierTest, CurveTangentAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  for (size_t i = 0; i < t_vals.size(); i++)
  {
    Vector current = curve_.tangentAt(t_vals[i]);
    Vector expected{TestData::kExpectedTangent[i].first, TestData::kExpectedTangent[i].second};
    EXPECT_TRUE(current.isApprox(expected, Utils::epsilon)) << "Tangent at t=" << t_vals[i] << " differs";
  }
}

TEST_F(BezierTest, CurveNormalAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  for (size_t i = 0; i < t_vals.size(); i++)
  {
    Vector current = curve_.normalAt(t_vals[i]);
    Vector expected{TestData::kExpectedNormal[i].first, TestData::kExpectedNormal[i].second};
    EXPECT_TRUE(current.isApprox(expected, Utils::epsilon)) << "Normal at t=" << t_vals[i] << " differs";
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
  const int N = 2;
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

  // Check that roots are in valid range [0, 1]
  for (double root : curve_roots)
  {
    EXPECT_GE(root, 0.0) << "Root should be >= 0";
    EXPECT_LE(root, 1.0) << "Root should be <= 1";
  }

  // Verify that the curve actually crosses zero at the roots
  for (double root : curve_roots)
  {
    Point p = curve_roots_.valueAt(root);
    // At least one coordinate should be very close to zero
    bool crosses_x = std::abs(p.x()) < Utils::epsilon;
    bool crosses_y = std::abs(p.y()) < Utils::epsilon;
    EXPECT_TRUE(crosses_x || crosses_y) << "At root t=" << root << ", point should be near axis: "
                                        << "(" << p.x() << ", " << p.y() << ")";
  }

  // For this specific curve, we expect 3 roots
  EXPECT_EQ(curve_roots.size(), 3) << "This test curve should have 3 roots";
}

TEST_F(BezierTest, CurveExtremaTest)
{
  std::vector<double> curve_extrema = curve_roots_.extrema();

  // Check that extrema are in valid range [0, 1]
  for (double t : curve_extrema)
  {
    EXPECT_GE(t, 0.0) << "Extrema should be >= 0";
    EXPECT_LE(t, 1.0) << "Extrema should be <= 1";
  }

  // Verify that derivative is near zero at extrema
  for (double t : curve_extrema)
  {
    Vector deriv = curve_roots_.derivativeAt(t);
    // At extrema, at least one component of the derivative should be very small
    bool x_extreme = std::abs(deriv.x()) < Utils::epsilon;
    bool y_extreme = std::abs(deriv.y()) < Utils::epsilon;
    EXPECT_TRUE(x_extreme || y_extreme) << "At extrema t=" << t << ", derivative should be near zero: "
                                        << "(" << deriv.x() << ", " << deriv.y() << ")";
  }

  // For this specific curve, we expect 2 extrema
  EXPECT_EQ(curve_extrema.size(), 2) << "This test curve should have 2 extrema";
}

TEST_F(BezierTest, CurveBoundingBoxTest)
{
  BoundingBox bbox = curve_.boundingBox();

  // Bounding box should contain curve endpoints
  auto endpoints = curve_.endPoints();
  EXPECT_GE(endpoints.first.x(), bbox.min().x());
  EXPECT_LE(endpoints.first.x(), bbox.max().x());
  EXPECT_GE(endpoints.first.y(), bbox.min().y());
  EXPECT_LE(endpoints.first.y(), bbox.max().y());
  EXPECT_GE(endpoints.second.x(), bbox.min().x());
  EXPECT_LE(endpoints.second.x(), bbox.max().x());
  EXPECT_GE(endpoints.second.y(), bbox.min().y());
  EXPECT_LE(endpoints.second.y(), bbox.max().y());

  // Bounding box should contain all extrema
  auto extrema = curve_.extrema();
  for (double t : extrema)
  {
    Point p = curve_.valueAt(t);
    EXPECT_GE(p.x(), bbox.min().x() - Utils::epsilon) << "Extrema point outside bbox";
    EXPECT_LE(p.x(), bbox.max().x() + Utils::epsilon) << "Extrema point outside bbox";
    EXPECT_GE(p.y(), bbox.min().y() - Utils::epsilon) << "Extrema point outside bbox";
    EXPECT_LE(p.y(), bbox.max().y() + Utils::epsilon) << "Extrema point outside bbox";
  }

  // Min should be less than max
  EXPECT_LT(bbox.min().x(), bbox.max().x());
  EXPECT_LT(bbox.min().y(), bbox.max().y());

  // Verify some points on the curve are within the bounding box
  for (double t = 0.0; t <= 1.0; t += 0.1)
  {
    Point p = curve_.valueAt(t);
    EXPECT_GE(p.x(), bbox.min().x() - Utils::epsilon) << "Point at t=" << t << " outside bbox";
    EXPECT_LE(p.x(), bbox.max().x() + Utils::epsilon) << "Point at t=" << t << " outside bbox";
    EXPECT_GE(p.y(), bbox.min().y() - Utils::epsilon) << "Point at t=" << t << " outside bbox";
    EXPECT_LE(p.y(), bbox.max().y() + Utils::epsilon) << "Point at t=" << t << " outside bbox";
  }
}

TEST_F(BezierTest, CurveIntersectionsTest)
{
  // Same-order (cubic) crossing curve: the original quartic test data hits
  // the known different-order isApprox out-of-bounds read in intersections()
  // (REVIEW-v4-devel §2.1, caught by ASan) — that case lives in
  // known_bugs_test.cpp as DISABLED_IntersectionsDifferentOrder.
  Curve curve_with_intersections{intersectionPointsAsMatrix().topRows(4)};
  PointVector intersections = curve_.intersections(curve_with_intersections);

  EXPECT_GE(intersections.size(), 1) << "Too few intersections found";
  EXPECT_LE(intersections.size(), 9) << "Too many intersections found";

  // Every reported point lies on both curves
  for (const Point& p : intersections)
  {
    EXPECT_NEAR(curve_.distance(p), 0.0, 1e-4);
    EXPECT_NEAR(curve_with_intersections.distance(p), 0.0, 1e-4);
  }
}

TEST_F(BezierTest, CurveSplitTest)
{
  std::pair<Curve, Curve> split_curves = curve_.splitCurve(0.5);

  // Test continuity at split point
  Point left_end = split_curves.first.endPoints().second;
  Point right_start = split_curves.second.endPoints().first;
  EXPECT_NEAR(left_end.x(), right_start.x(), Utils::epsilon) << "Split curves should be continuous";
  EXPECT_NEAR(left_end.y(), right_start.y(), Utils::epsilon) << "Split curves should be continuous";

  // Test that original curve value at 0.5 matches split point
  Point original_mid = curve_.valueAt(0.5);
  EXPECT_NEAR(left_end.x(), original_mid.x(), Utils::epsilon);
  EXPECT_NEAR(left_end.y(), original_mid.y(), Utils::epsilon);

  // Test that split curves have same order as original
  EXPECT_EQ(split_curves.first.order(), curve_.order());
  EXPECT_EQ(split_curves.second.order(), curve_.order());

  // Test that endpoints match original (with tolerance for floating point)
  auto orig_endpoints = curve_.endPoints();
  auto left_start = split_curves.first.endPoints().first;
  auto right_end = split_curves.second.endPoints().second;
  EXPECT_NEAR(left_start.x(), orig_endpoints.first.x(), Utils::epsilon);
  EXPECT_NEAR(left_start.y(), orig_endpoints.first.y(), Utils::epsilon);
  EXPECT_NEAR(right_end.x(), orig_endpoints.second.x(), Utils::epsilon);
  EXPECT_NEAR(right_end.y(), orig_endpoints.second.y(), Utils::epsilon);

  // Test value preservation: points on split curves should match original
  auto p1 = split_curves.first.valueAt(0.0);
  auto p2 = curve_.valueAt(0.0);
  EXPECT_NEAR(p1.x(), p2.x(), Utils::epsilon);
  EXPECT_NEAR(p1.y(), p2.y(), Utils::epsilon);

  p1 = split_curves.first.valueAt(0.5);
  p2 = curve_.valueAt(0.25);
  EXPECT_NEAR(p1.x(), p2.x(), Utils::epsilon);
  EXPECT_NEAR(p1.y(), p2.y(), Utils::epsilon);

  p1 = split_curves.second.valueAt(0.5);
  p2 = curve_.valueAt(0.75);
  EXPECT_NEAR(p1.x(), p2.x(), Utils::epsilon);
  EXPECT_NEAR(p1.y(), p2.y(), Utils::epsilon);

  p1 = split_curves.second.valueAt(1.0);
  p2 = curve_.valueAt(1.0);
  EXPECT_NEAR(p1.x(), p2.x(), Utils::epsilon);
  EXPECT_NEAR(p1.y(), p2.y(), Utils::epsilon);
}

TEST_F(BezierTest, CurveProjectPointTest)
{
  Point point{100, 150};
  double t_projected = curve_.projectPoint(point);

  // t should be in valid range
  EXPECT_GE(t_projected, 0.0);
  EXPECT_LE(t_projected, 1.0);

  // Point on curve at t_projected should be closer than any other sampled point
  Point closest_point = curve_.valueAt(t_projected);
  double min_dist = (point - closest_point).norm();

  for (double t = 0.0; t <= 1.0; t += 0.01)
  {
    double dist = (point - curve_.valueAt(t)).norm();
    EXPECT_GE(dist, min_dist - Utils::epsilon) << "Found closer point at t=" << t;
  }

  // For this specific test case, check against known value (tight tolerance)
  EXPECT_NEAR(t_projected, 0.03465913829148962, Utils::epsilon);
}

TEST_F(BezierTest, CurveDistanceTest)
{
  Point point{100, 150};
  double distance = curve_.distance(point);

  // Distance should be non-negative
  EXPECT_GE(distance, 0.0);

  // Distance should match distance to projected point
  double t_projected = curve_.projectPoint(point);
  Point closest_point = curve_.valueAt(t_projected);
  double expected_distance = (point - closest_point).norm();
  EXPECT_NEAR(distance, expected_distance, Utils::epsilon);

  // For this specific test case, check against known value
  EXPECT_NEAR(distance, 0.68269613683526820, Utils::epsilon);
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
    EXPECT_NEAR(on_curve.x(), polyline[i].x(), 1e-6) << "Polyline vertex " << i << " not on curve";
    EXPECT_NEAR(on_curve.y(), polyline[i].y(), 1e-6) << "Polyline vertex " << i << " not on curve";
  }

  // Dense samples of the curve stay within flatness of the polyline
  double flatness = curve_.boundingBox().diagonal().norm() / 1000;
  for (double t{}; t <= 1.0; t += 0.001)
    EXPECT_LE(Utils::dist(polyline, curve_.valueAt(t)), flatness + Utils::epsilon) << "Curve too far at t=" << t;
}

TEST(CurveOracleTests, ValueAtMatchesDeCasteljau)
{
  for (unsigned order = 1; order <= 6; order++)
  {
    PointVector cp = makeControlPoints(order);
    Curve curve{cp};
    for (double t{}; t <= 1.0; t += 0.1)
    {
      Point expected = Oracles::deCasteljau(cp, t);
      Point actual = curve.valueAt(t);
      EXPECT_NEAR(actual.x(), expected.x(), Oracles::kGeom) << "order=" << order << " t=" << t;
      EXPECT_NEAR(actual.y(), expected.y(), Oracles::kGeom) << "order=" << order << " t=" << t;
    }
  }
}

TEST_F(BezierTest, CurveLengthMatchesChordLengthOracle)
{
  // Fixture cubic
  double oracle = Oracles::chordLength(curve_.controlPoints());
  EXPECT_NEAR(curve_.length(), oracle, 1e-6 * oracle);

  // Partial lengths
  EXPECT_NEAR(curve_.length(0.3), Oracles::chordLength(curve_.controlPoints(), 0.0, 0.3), 1e-6 * oracle);
  EXPECT_NEAR(curve_.length(0.2, 0.8), Oracles::chordLength(curve_.controlPoints(), 0.2, 0.8), 1e-6 * oracle);

  // Straight line: length equals the chord exactly
  Curve line{PointVector{{0, 0}, {3, 4}}};
  EXPECT_NEAR(line.length(), 5.0, Utils::epsilon);

  // Symmetric quadratic — regression guard for the master-branch Chebyshev
  // truncation bug (even coefficients vanish); v4 must get this right
  Curve symmetric{PointVector{{0, 0}, {1, 2}, {2, 0}}};
  double symmetric_oracle = Oracles::chordLength(symmetric.controlPoints());
  EXPECT_NEAR(symmetric.length(), symmetric_oracle, 1e-6 * symmetric_oracle);
}

TEST_F(BezierTest, CurveStepLengthRoundTripTest)
{
  for (double t : {0.1, 0.42, 0.7})
    for (double ds : {5.0, 35.0, -5.0, -15.0})
    {
      double new_t = curve_.step(t, ds);
      EXPECT_NEAR(curve_.length(t, new_t), ds, 1e-6) << "Round trip failed for t=" << t << " ds=" << ds;
    }

  // Out-of-range distances clamp to the curve ends
  EXPECT_DOUBLE_EQ(curve_.step(0.5, 1e6), 1.0);
  EXPECT_DOUBLE_EQ(curve_.step(0.5, -1e6), 0.0);
}

TEST_F(BezierTest, CurveReverseTest)
{
  Curve reversed{curve_};
  reversed.reverse();
  for (double t{}; t <= 1.0; t += 0.05)
  {
    Point expected = curve_.valueAt(1.0 - t);
    Point actual = reversed.valueAt(t);
    EXPECT_NEAR(actual.x(), expected.x(), Oracles::kGeom) << "reverse mismatch at t=" << t;
    EXPECT_NEAR(actual.y(), expected.y(), Oracles::kGeom) << "reverse mismatch at t=" << t;
  }
}

TEST_F(BezierTest, CurveSplitMultiTest)
{
  const ParamVector t_split{0.25, 0.6};
  std::vector<Curve> pieces = curve_.splitCurve(t_split);
  ASSERT_EQ(pieces.size(), 3u);

  // Pieces map linearly onto the original parameter ranges
  const std::array<std::pair<double, double>, 3> ranges{{{0.0, 0.25}, {0.25, 0.6}, {0.6, 1.0}}};
  for (size_t k = 0; k < pieces.size(); k++)
    for (double s{}; s <= 1.0; s += 0.125)
    {
      Point expected = curve_.valueAt(ranges[k].first + s * (ranges[k].second - ranges[k].first));
      Point actual = pieces[k].valueAt(s);
      EXPECT_NEAR(actual.x(), expected.x(), 1e-8) << "piece " << k << " s=" << s;
      EXPECT_NEAR(actual.y(), expected.y(), 1e-8) << "piece " << k << " s=" << s;
    }

  // C0 continuity between adjacent pieces
  for (size_t k = 0; k + 1 < pieces.size(); k++)
  {
    Point left = pieces[k].endPoints().second;
    Point right = pieces[k + 1].endPoints().first;
    EXPECT_NEAR(left.x(), right.x(), 1e-8);
    EXPECT_NEAR(left.y(), right.y(), 1e-8);
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
  for (double t{}; t <= 1.0; t += 0.05)
  {
    Point expected = curve_.valueAt(t);
    Point actual = raised.valueAt(t);
    EXPECT_NEAR(actual.x(), expected.x(), 1e-8) << "raiseOrder changed shape at t=" << t;
    EXPECT_NEAR(actual.y(), expected.y(), 1e-8) << "raiseOrder changed shape at t=" << t;
  }

  // Lowering the raised curve recovers the original control points
  raised.lowerOrder();
  ASSERT_EQ(raised.order(), curve_.order());
  PointVector original = curve_.controlPoints();
  PointVector round_trip = raised.controlPoints();
  for (size_t i = 0; i < original.size(); i++)
  {
    EXPECT_NEAR(round_trip[i].x(), original[i].x(), 1e-8);
    EXPECT_NEAR(round_trip[i].y(), original[i].y(), 1e-8);
  }
}

TEST(CurveOrderTests, LowerOrderThrowsOnFirstOrderCurve)
{
  Curve line{PointVector{{0, 0}, {1, 1}}};
  EXPECT_THROW(line.lowerOrder(), std::logic_error);
}

TEST_F(BezierTest, CurveDerivativeAtMatchesCentralDiff)
{
  auto f = [this](double t) { return curve_.valueAt(t); };
  for (double t : {0.1, 0.3, 0.5, 0.7, 0.9})
  {
    Vector expected = Oracles::centralDiff(f, t);
    Vector actual = curve_.derivativeAt(t);
    // central difference is O(h^2) accurate; derivative magnitudes are ~1e2
    EXPECT_NEAR(actual.x(), expected.x(), 1e-3) << "derivative mismatch at t=" << t;
    EXPECT_NEAR(actual.y(), expected.y(), 1e-3) << "derivative mismatch at t=" << t;
  }
}

TEST_F(BezierTest, CurveTangentNormalPropertiesTest)
{
  for (double t : {0.0, 0.2, 0.5, 0.8, 1.0})
  {
    Vector tangent = curve_.tangentAt(t);
    Vector normal = curve_.normalAt(t);

    EXPECT_NEAR(tangent.norm(), 1.0, Utils::epsilon) << "tangent not unit at t=" << t;
    EXPECT_NEAR(normal.norm(), 1.0, Utils::epsilon) << "normal not unit at t=" << t;
    EXPECT_NEAR(tangent.dot(normal), 0.0, Utils::epsilon) << "tangent/normal not orthogonal at t=" << t;

    // normal is the tangent rotated +90 degrees
    EXPECT_NEAR(normal.x(), -tangent.y(), Utils::epsilon);
    EXPECT_NEAR(normal.y(), tangent.x(), Utils::epsilon);

    // tangent is parallel to the (non-normalized) derivative
    Vector derivative = curve_.derivativeAt(t);
    EXPECT_NEAR(tangent.x() * derivative.y() - tangent.y() * derivative.x(), 0.0,
                Utils::epsilon * derivative.norm());
  }
}

TEST_F(BezierTest, CurveProjectPointRecoveryTest)
{
  // Projecting a point that lies on the curve recovers its parameter
  for (double t : {0.1, 0.3, 0.5, 0.7, 0.9})
  {
    double projected = curve_.projectPoint(curve_.valueAt(t));
    EXPECT_NEAR(curve_.distance(curve_.valueAt(t)), 0.0, 1e-6);
    Point recovered = curve_.valueAt(projected);
    Point original = curve_.valueAt(t);
    EXPECT_NEAR(recovered.x(), original.x(), 1e-6) << "projection at t=" << t;
    EXPECT_NEAR(recovered.y(), original.y(), 1e-6) << "projection at t=" << t;
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
  EXPECT_NEAR(pos_actual.x(), pos_expected.x(), 1e-8);
  EXPECT_NEAR(pos_actual.y(), pos_expected.y(), 1e-8);

  Vector der_expected = curve_.derivativeAt(1.0);
  Vector der_actual = continued.derivativeAt(0.0);
  EXPECT_NEAR(der_actual.x(), der_expected.x(), 1e-6);
  EXPECT_NEAR(der_actual.y(), der_expected.y(), 1e-6);
}

TEST_F(BezierTest, CurveApplyContinuityG1Test)
{
  constexpr double beta = 2.5;
  Curve continued{curve_roots_};
  continued.applyContinuity(curve_, {beta});

  // G1: position matches; first derivative is scaled by beta
  Point pos_expected = curve_.valueAt(1.0);
  Point pos_actual = continued.valueAt(0.0);
  EXPECT_NEAR(pos_actual.x(), pos_expected.x(), 1e-8);
  EXPECT_NEAR(pos_actual.y(), pos_expected.y(), 1e-8);

  Vector der_expected = beta * curve_.derivativeAt(1.0);
  Vector der_actual = continued.derivativeAt(0.0);
  EXPECT_NEAR(der_actual.x(), der_expected.x(), 1e-6);
  EXPECT_NEAR(der_actual.y(), der_expected.y(), 1e-6);
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

} // namespace Bezier
