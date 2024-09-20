#ifndef BEZIER_UNIT_TESTS_HPP
#define BEZIER_UNIT_TESTS_HPP

#include "Bezier/bezier.h"
#include "Bezier/declarations.h"
#include "Eigen/Dense"

namespace Bezier {
static Eigen::MatrixX2d getCurvePointsAsEigenMatrix() {
  Eigen::MatrixX2d curve_pts;
  curve_pts.resize(4, 2);
  curve_pts << 84, 162, 246, 30, 48, 236, 180, 110;
  return curve_pts;
}

static Eigen::MatrixX2d getRootsCurvePointsAsEigenMatrix() {
  Eigen::MatrixX2d curve_pts;
  curve_pts.resize(4, 2);
  curve_pts << -50, -50, 75, 48, 64, 65, 50, -50;
  return curve_pts;
}

static Eigen::MatrixX2d getIntersectionCurvePointsAsEigenMatrix() {
  Eigen::MatrixX2d curve_pts;
  curve_pts.resize(5, 2);
  curve_pts << 180, 110, 175, 160, 60, 48, 164, 165, 124, 134;
  return curve_pts;
}

static PointVector getCurvePointsAsPointVector() {
  PointVector curve_pts{{84, 162}, {246, 30}, {48, 236}, {180, 110}};
  return curve_pts;
}

static PointVector getExpectedPolylinePoints() {
  // clang-format off
  return {{84, 162}, {110.32470703124997, 141.04736328124994}, {129.22265624999997, 127.03515624999991}, {141.70458984374997, 118.98193359374994}, {148.78124999999994, 115.90624999999991}, {150.60845947265619, 115.92828369140619}, {151.46337890624991, 116.82666015624993}, {150.76171874999994, 120.76171874999994}, {143.24999999999994, 133.74999999999991}, {134.33203124999994, 147.01953124999994}, {131.87255859374994, 151.30615234374994}, {132.09374999999991, 152.71874999999991}, {144.62109374999994, 142.99609374999994}, {179.99999999999991, 109.99999999999991}};
  // clang-format on
}

static Eigen::MatrixX2d getExpectedValueAt() {
  Eigen::MatrixX2d mat;
  mat.resize(5, 2);
  mat << 83.99999999999998579, 161.99999999999997158, 148.78124999999997158,
      115.90625000000000000, 143.24999999999994316, 133.75000000000000000,
      132.09374999999994316, 152.71875000000005684, 180.00000000000002842,
      109.99999999999997158;
  return mat;
}

}  // namespace Bezier

#endif  // BEZIER_UNIT_TESTS_HPP