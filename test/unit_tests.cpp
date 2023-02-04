#include "unit_tests.hpp"

#include <iostream>

#include "Bezier/bezier.h"
#include "Bezier/declarations.h"
#include "Eigen/Dense"
#include "gtest/gtest.h"

namespace Bezier {

class BezierTest : public ::testing::Test {
 public:
  BezierTest() : curve_1_{getCurvePointsAsEigenMatrix()} {}

 protected:
  Curve curve_1_;
};

template <typename T>
void CoutFixed(T out) {
  std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
            << std::fixed << out << ", ";
}

TEST(ConstructorTests, ConstructFromEigenMatrix) {
  EXPECT_NO_THROW(Curve{getCurvePointsAsEigenMatrix()})
      << "Cannot create Curve with Eigen::MatrixX2d.";
}
TEST(ConstructorTests, ConstructFromPointVector) {
  EXPECT_NO_THROW(Curve{getCurvePointsAsPointVector()})
      << "Cannot create Curve with PointVector.";
}

TEST_F(BezierTest, CurveOrderTest) {
  EXPECT_EQ(curve_1_.order(), 3) << "Curve order differs from expected.";
}

TEST_F(BezierTest, CurveControlPointsTest) {
  PointVector points_expected = getCurvePointsAsPointVector();
  PointVector points_curve = curve_1_.controlPoints();
  EXPECT_EQ(points_expected, points_curve)
      << "Curve control points differ from expected";
}

TEST_F(BezierTest, CurveControlPointTest) {
  Point point_expected = getCurvePointsAsPointVector().at(2);
  Point point_curve = curve_1_.controlPoint(2);
  EXPECT_EQ(point_expected, point_curve)
      << "Curve control point at idx differs from expected";
}

TEST_F(BezierTest, CurveEndPointsTest) {
  std::pair<Point, Point> end_points_expected{
      getCurvePointsAsPointVector().at(0),
      getCurvePointsAsPointVector().back()};
  std::pair<Point, Point> point_curve = curve_1_.endPoints();
  EXPECT_EQ(end_points_expected, point_curve)
      << "Curve end points differ from expected";
}

TEST_F(BezierTest, CurvePolylineTest) {
  PointVector polyline = curve_1_.polyline();
  PointVector expected_polyline = getExpectedPolylinePoints();
  ASSERT_EQ(polyline.size(), expected_polyline.size())
      << "polyline size differs from expected";
  for (int i = 0; i < polyline.size(); i++) {
    EXPECT_EQ(polyline[i], expected_polyline[i]);
  }
}

TEST_F(BezierTest, CurveLengthTest) {
  EXPECT_NEAR(curve_1_.length(), 189.11484365930446, 1e-6)
      << "Curve length differs from expected";
}

TEST_F(BezierTest, CurveIterateByLengthTest) {
  EXPECT_NEAR(curve_1_.iterateByLength(0.2, 35), 0.5513946586820645, 1e-6)
      << "Parameter t returned from the function differs from expected";
}

TEST_F(BezierTest, CurveSetControlPointTest) {
  Point new_control_point{40, 40};
  EXPECT_NO_THROW(curve_1_.setControlPoint(2, new_control_point))
      << "Curve set control point failed";
  EXPECT_EQ(curve_1_.controlPoint(2), new_control_point)
      << "Failed to set control point";
}

TEST_F(BezierTest, CurveElevateOrderTest) {
  EXPECT_NO_THROW(curve_1_.elevateOrder()) << "Curve::elevateOrder failed";
}

TEST_F(BezierTest, CurveLowerOrderTest) {
  EXPECT_NO_THROW(curve_1_.lowerOrder()) << "Curve::lowerOrder failed";
}

TEST_F(BezierTest, CurveValueAtTest) {
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  PointVector expected_points{getExpectedValueAtPoints()};
  for (int i = 0; i < t_vals.size(); i++) {
    Point current_point = curve_1_.valueAt(t_vals[i]);
    EXPECT_EQ(current_point, expected_points[i]);
  }
}

TEST_F(BezierTest, CurveValueAtMultipleParamsTest) {
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  PointVector expected_points{getExpectedValueAtPoints()};
  Eigen::MatrixX2d points = curve_1_.valueAt(t_vals);
  EXPECT_EQ(points, getExpectedValueAtMatrix());
}

TEST_F(BezierTest, CurveCurvatureAtTest) {
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  std::vector<double> expected{0.00052864283054623, 0.13400265200820813,
                               0.00339166466643027, -0.83667680576075631,
                               -0.00024618707208090};
  for (int i = 0; i < t_vals.size(); i++) {
    double current = curve_1_.curvatureAt(t_vals[i]);
    EXPECT_NEAR(current, expected[i], 1e-12);
  }
}

TEST_F(BezierTest, CurveCurvatureDerivativeAtTest) {
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  std::vector<double> expected{0.00620125036871112, 6.26701004231499148,
                               -0.06892991355501528, 86.75145336780214222,
                               0.00362325743277114};
  for (int i = 0; i < t_vals.size(); i++) {
    double current = curve_1_.curvatureDerivativeAt(t_vals[i]);
    EXPECT_NEAR(current, expected[i], 1e-12);
  }
}

TEST_F(BezierTest, CurveTangentAtTest) {
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  std::vector<Vector> expected{{0.77523498551728931, -0.63167295116223576},
                               {0.98169156767741872, -0.19047746835532003},
                               {-0.55219905682528458, 0.83371230148131203},
                               {0.98328200498446017, -0.18208926018230745},
                               {0.72335554414357217, -0.69047574668250067}};
  for (int i = 0; i < t_vals.size(); i++) {
    Vector current = curve_1_.tangentAt(t_vals[i]);
    EXPECT_EQ(current, expected[i]);
  }
}

TEST_F(BezierTest, CurveNormalAtTest) {
  std::vector<double> t_vals{0., 0.25, 0.5, 0.75, 1};
  std::vector<Vector> expected{{0.63167295116223576, 0.77523498551728931},
                               {0.19047746835532003, 0.98169156767741872},
                               {-0.83371230148131203, -0.55219905682528458},
                               {0.18208926018230745, 0.98328200498446017},
                               {0.69047574668250067, 0.72335554414357217}};
  for (int i = 0; i < t_vals.size(); i++) {
    Vector current = curve_1_.normalAt(t_vals[i]);
    EXPECT_EQ(current, expected[i]);
  }
}

}  // namespace Bezier

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// TODO: manipulateCurvature test
/*
std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
          << std::fixed << "{" << current_point.x() << ", "
          << current_point.y() << "}, ";
*/