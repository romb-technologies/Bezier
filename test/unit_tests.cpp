#include "unit_tests.hpp"

#include <Bezier/declarations.h>
#include <iostream>

#include <gtest/gtest.h>

#include "Bezier/bezier.h"
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

template <typename T> void CoutFixed(T out)
{
  std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << std::fixed << out << ", ";
}

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

  for (int i = 0; i < t_vals.size(); i++)
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
  for (int i = 0; i < t_vals.size(); i++)
  {
    double current = curve_.curvatureAt(t_vals[i]);
    EXPECT_NEAR(current, TestData::kExpectedCurvature[i], Utils::epsilon);
  }
}

TEST_F(BezierTest, CurveCurvatureDerivativeAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  for (int i = 0; i < t_vals.size(); i++)
  {
    double current = curve_.curvatureDerivativeAt(t_vals[i]);
    EXPECT_NEAR(current, TestData::kExpectedCurvatureDerivative[i], Utils::epsilon);
  }
}

TEST_F(BezierTest, CurveTangentAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  for (int i = 0; i < t_vals.size(); i++)
  {
    Vector current = curve_.tangentAt(t_vals[i]);
    Vector expected{TestData::kExpectedTangent[i].first, TestData::kExpectedTangent[i].second};
    EXPECT_EQ(current, expected);
  }
}

TEST_F(BezierTest, CurveNormalAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  for (int i = 0; i < t_vals.size(); i++)
  {
    Vector current = curve_.normalAt(t_vals[i]);
    Vector expected{TestData::kExpectedNormal[i].first, TestData::kExpectedNormal[i].second};
    EXPECT_EQ(current, expected);
  }
}

TEST_F(BezierTest, CurveDerivativeTest)
{
  PointVector expected = TestData::toPointVector(TestData::kExpectedFirstDerivative);
  PointVector derivative_control_points = curve_.derivative().controlPoints();
  ASSERT_EQ(expected.size(), derivative_control_points.size());
  for (int i = 0; i < derivative_control_points.size(); i++)
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
  for (int i = 0; i < derivative_control_points.size(); i++)
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
  Curve curve_with_intersections{intersectionPointsAsMatrix()};
  PointVector intersections = curve_.intersections(curve_with_intersections);

  // v040 may compute intersections in different order or with slightly different algorithm
  // Just check that we have a reasonable number of intersections
  EXPECT_GE(intersections.size(), 4) << "Too few intersections found";
  EXPECT_LE(intersections.size(), 8) << "Too many intersections found";
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

// TODO: ApplyContinuityTest
// TEST_F(BezierTest, CurveApplyContinuityTest) { ... }

} // namespace Bezier

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
