#include "unit_tests.hpp"

#include <iostream>

#include <gtest/gtest.h>

#include "Bezier/bezier.h"

namespace Bezier
{

class BezierTest : public ::testing::Test
{
public:
  BezierTest() : curve_{getCurvePointsAsEigenMatrix()}, curve_roots_{getRootsCurvePointsAsEigenMatrix()} {}

protected:
  Curve curve_;
  // Curve that passes through axes
  Curve curve_roots_;
};

template <typename T> void CoutFixed(T out)
{
  std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << std::fixed << out << ", ";
}

TEST(ConstructorTests, ConstructFromEigenMatrix)
{
  EXPECT_NO_THROW(Curve{getCurvePointsAsEigenMatrix()}) << "Cannot create Curve with Eigen::MatrixX2d.";
}
TEST(ConstructorTests, ConstructFromPointVector)
{
  EXPECT_NO_THROW(Curve{getCurvePointsAsPointVector()}) << "Cannot create Curve with PointVector.";
}

TEST_F(BezierTest, CurveOrderTest) { EXPECT_EQ(curve_.order(), 3) << "Curve order differs from expected."; }

TEST_F(BezierTest, CurveControlPointsTest)
{
  PointVector points_expected = getCurvePointsAsPointVector();
  PointVector points_curve = curve_.controlPoints();
  EXPECT_EQ(points_expected, points_curve) << "Curve control points differ from expected";
}

TEST_F(BezierTest, CurveControlPointTest)
{
  Point point_expected = getCurvePointsAsPointVector().at(2);
  Point point_curve = curve_.controlPoint(2);
  EXPECT_EQ(point_expected, point_curve) << "Curve control point at idx differs from expected";
}

TEST_F(BezierTest, CurveEndPointsTest)
{
  std::pair<Point, Point> end_points_expected{getCurvePointsAsPointVector().at(0),
                                              getCurvePointsAsPointVector().back()};
  std::pair<Point, Point> point_curve = curve_.endPoints();
  EXPECT_EQ(end_points_expected, point_curve) << "Curve end points differ from expected";
}

TEST_F(BezierTest, CurvePolylineTest)
{
  PointVector polyline = curve_.polyline();
  PointVector expected_polyline = getExpectedPolylinePoints();
  ASSERT_EQ(polyline.size(), expected_polyline.size()) << "Polyline size differs from expected";
  for (int i = 0; i < polyline.size(); i++)
  {
    EXPECT_EQ(polyline[i], expected_polyline[i]) << "Polyline differs from expected";
  }
}

TEST_F(BezierTest, CurveLengthTest)
{
  EXPECT_NEAR(curve_.length(), 189.11484365930446, 1e-8) << "Curve length differs from expected";
}

TEST_F(BezierTest, CurveIterateByLengthTest)
{
  EXPECT_NEAR(curve_.iterateByLength(0.2, 35), 0.5513946586820645, 1e-8)
      << "Parameter t returned from the function differs from expected";
}

TEST_F(BezierTest, CurveSetControlPointTest)
{
  Point new_control_point{40, 40};
  EXPECT_NO_THROW(curve_.setControlPoint(2, new_control_point)) << "Curve set control point failed";
  EXPECT_EQ(curve_.controlPoint(2), new_control_point) << "Failed to set control point";
}

TEST_F(BezierTest, CurveElevateOrderTest) { EXPECT_NO_THROW(curve_.elevateOrder()) << "Curve::elevateOrder failed"; }

TEST_F(BezierTest, CurveLowerOrderTest) { EXPECT_NO_THROW(curve_.lowerOrder()) << "Curve::lowerOrder failed"; }

TEST_F(BezierTest, CurveValueAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  PointVector expected_points{{83.99999999999998579, 161.99999999999997158},
                              {148.78124999999994316, 115.90625000000000000},
                              {143.24999999999994316, 133.75000000000000000},
                              {132.09374999999991473, 152.71875000000000000},
                              {180.00000000000002842, 109.99999999999997158}};
  for (int i = 0; i < t_vals.size(); i++)
  {
    Point current_point = curve_.valueAt(t_vals[i]);
    EXPECT_EQ(current_point, expected_points[i]) << "Curve Point at t = " << t_vals[i] << " differs from expected";
  }
}

TEST_F(BezierTest, CurveValueAtMultipleParamsTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  Eigen::MatrixX2d points = curve_.valueAt(t_vals);
  EXPECT_EQ(points, getExpectedValueAt()) << "Matrix differs from expected";
}

TEST_F(BezierTest, CurveCurvatureAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  std::vector<double> expected{0.00052864283054623, 0.13400265200820813, 0.00339166466643027, -0.83667680576075631,
                               -0.00024618707208090};
  for (int i = 0; i < t_vals.size(); i++)
  {
    double current = curve_.curvatureAt(t_vals[i]);
    EXPECT_NEAR(current, expected[i], 1e-12) << "Curvature values differ from expected";
  }
}

TEST_F(BezierTest, CurveCurvatureDerivativeAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  std::vector<double> expected{0.00620125036871112, 6.26701004231499148, -0.06892991355501528, 86.75145336780214222,
                               0.00362325743277114};
  for (int i = 0; i < t_vals.size(); i++)
  {
    double current = curve_.curvatureDerivativeAt(t_vals[i]);
    EXPECT_NEAR(current, expected[i], 1e-12) << "Curvature derivative differs from expected";
  }
}

TEST_F(BezierTest, CurveTangentAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  std::vector<Vector> expected{{0.77523498551728931, -0.63167295116223576},
                               {0.98169156767741872, -0.19047746835532003},
                               {-0.55219905682528458, 0.83371230148131203},
                               {0.98328200498446017, -0.18208926018230745},
                               {0.72335554414357217, -0.69047574668250067}};
  for (int i = 0; i < t_vals.size(); i++)
  {
    Vector current = curve_.tangentAt(t_vals[i]);
    EXPECT_EQ(current, expected[i]) << "Curve Tangent at t = " << t_vals[i] << " differs from expected";
  }
}

TEST_F(BezierTest, CurveNormalAtTest)
{
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  std::vector<Vector> expected{{0.63167295116223576, 0.77523498551728931},
                               {0.19047746835532003, 0.98169156767741872},
                               {-0.83371230148131203, -0.55219905682528458},
                               {0.18208926018230745, 0.98328200498446017},
                               {0.69047574668250067, 0.72335554414357217}};
  for (int i = 0; i < t_vals.size(); i++)
  {
    Vector current = curve_.normalAt(t_vals[i]);
    EXPECT_EQ(current, expected[i]) << "Curve Normal at t = " << t_vals[i] << " differs from expected";
  }
}

TEST_F(BezierTest, CurveDerivativeTest)
{
  PointVector expected_control_points{{486., -396.}, {-594., 618.}, {396., -378.}};
  PointVector derivative_control_points = curve_.derivative().controlPoints();
  ASSERT_EQ(expected_control_points.size(), derivative_control_points.size())
      << "Curve derivative control points size differs from expected";
  for (int i = 0; i < derivative_control_points.size(); i++)
  {
    EXPECT_EQ(expected_control_points[i], derivative_control_points[i])
        << "Derivative control point differs from expected";
  }
}

TEST_F(BezierTest, CurveNthDerivativeTest)
{
  const int N = 2;
  PointVector expected_control_points{{-2160., 2028.}, {1980., -1992.}};
  PointVector derivative_control_points = curve_.derivative(N).controlPoints();
  ASSERT_EQ(expected_control_points.size(), derivative_control_points.size())
      << "Second curve derivative control points size differs from expected";
  for (int i = 0; i < derivative_control_points.size(); i++)
  {
    EXPECT_EQ(expected_control_points[i], derivative_control_points[i])
        << "Second derivative control point differs from expected";
  }
}

TEST_F(BezierTest, CurveRootsTest)
{
  std::vector<double> expected_roots{0.15960764467040239, 0.20703566777685056, 0.81790457456613741};
  std::vector<double> curve_roots = curve_roots_.roots();
  ASSERT_EQ(expected_roots.size(), curve_roots.size()) << "Roots vector size differs from expected";
  for (int i = 0; i < curve_roots.size(); i++)
  {
    EXPECT_EQ(expected_roots[i], curve_roots[i]) << "Curve roots differ from expected";
  }
}

TEST_F(BezierTest, CurveExtremaTest)
{
  std::vector<double> expected_extrema{0.69733039467846381, 0.51985862621226775};
  std::vector<double> curve_extrema = curve_roots_.extrema();
  ASSERT_EQ(expected_extrema.size(), curve_extrema.size()) << "Extrema vector size differs from expected";
  for (int i = 0; i < curve_extrema.size(); i++)
  {
    EXPECT_EQ(expected_extrema[i], curve_extrema[i]) << "Curve extrema differs from expected";
  }
}

TEST_F(BezierTest, CurveBoundingBoxTest)
{
  BoundingBox curve_bounding_box = curve_.boundingBox();
  Point expected_min = {84, 110};
  Point expected_max = {180, 162};
  EXPECT_EQ(curve_bounding_box.min(), expected_min) << "Bounding box minimum differs from expected";
  EXPECT_EQ(curve_bounding_box.max(), expected_max) << "Bounding box maximum differs from expected";
}

TEST_F(BezierTest, CurveIntersectionsTest)
{
  Curve curve_with_intersections{getIntersectionCurvePointsAsEigenMatrix()};
  PointVector expected_intersections{
      {143.28410318501778420, 118.12287028981754133}, {150.87261159971575353, 120.47985696319850035},
      {179.99996184822515488, 110.00038146265423222}, {166.28506725972277991, 123.01852075448888968},
      {166.28401126949671607, 123.01855001520399924}, {128.65721095296441945, 127.43178611658311183}};
  PointVector intersections = curve_.intersections(curve_with_intersections);
  ASSERT_EQ(expected_intersections.size(), intersections.size()) << "Intersections vector size differs from expected";
  for (int i = 0; i < intersections.size(); i++)
  {
    EXPECT_EQ(expected_intersections[i], intersections[i]) << "Curve intersections differs from expected";
  }
}

TEST_F(BezierTest, CurveSplitTest)
{
  std::pair<Curve, Curve> split_curves = curve_.splitCurve();
  std::pair<PointVector, PointVector> split_curves_control_points{split_curves.first.controlPoints(),
                                                                  split_curves.second.controlPoints()};
  std::pair<PointVector, PointVector> expected_control_points = {{{84.00000000000001421, 161.99999999999997158},
                                                                  {165.00000000000000000, 95.99999999999997158},
                                                                  {156.00000000000000000, 114.49999999999997158},
                                                                  {143.25000000000000000, 133.74999999999997158}},
                                                                 {{143.25000000000000000, 133.74999999999997158},
                                                                  {130.49999999999997158, 153.00000000000000000},
                                                                  {114.00000000000000000, 172.99999999999997158},
                                                                  {179.99999999999997158, 109.99999999999997158}}};

  // Check first curve
  auto& curve = split_curves_control_points.first;
  auto& expected_curve = expected_control_points.first;
  ASSERT_EQ(curve.size(), expected_curve.size()) << "First curve control points size differs from expected";
  for (int i = 0; i < curve.size(); i++)
  {
    EXPECT_EQ(curve[i], expected_curve[i]) << "First curve control points differ from expected";
  }
  // Check second curve
  curve = split_curves_control_points.second;
  expected_curve = expected_control_points.second;
  ASSERT_EQ(curve.size(), expected_curve.size()) << "Second curve control points size differs from expected";
  for (int i = 0; i < curve.size(); i++)
  {
    EXPECT_EQ(curve[i], expected_curve[i]) << "Second curve control points differ from expected";
  }
}

TEST_F(BezierTest, CurveProjectPointTest)
{
  Point point{100, 150};
  double t_closest_to_point = curve_.projectPoint(point);
  EXPECT_NEAR(t_closest_to_point, 0.03465913829148962, 1e-12) << "Calculated `t` differs from expected";
}

TEST_F(BezierTest, CurveDistanceTest)
{
  Point point{100, 150};
  double distance_to_point = curve_.distance(point);
  EXPECT_EQ(distance_to_point, 0.68269613683526820) << "Calculated distance differs from expected";
}

TEST_F(BezierTest, CurveManipulateCurvatureTest)
{
  Point point{150, 150};
  // Manipulate curvature, so it passes through `point` at t = 0.5
  ASSERT_NO_THROW(curve_.manipulateCurvature(0.5, point));
  const PointVector& curve_control_points{curve_.controlPoints()};
  const PointVector expected_control_points{{84.00000000000000000, 162.00000000000000000},
                                            {255.00000000000022737, 51.66666666666668561},
                                            {57.00000000000022737, 257.66666666666668561},
                                            {180.00000000000000000, 110.00000000000000000}};
  ASSERT_EQ(expected_control_points.size(), curve_control_points.size())
      << "Curve size after manipulation differs from expected";
  for (int i = 0; i < curve_control_points.size(); i++)
  {
    EXPECT_EQ(expected_control_points[i], curve_control_points[i])
        << "Curve control points after manipulation differ from expected";
  }
}

// TODO:
/*
TEST_F(BezierTest, CurveApplyContinuityTest) {
}
*/

} // namespace Bezier

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
