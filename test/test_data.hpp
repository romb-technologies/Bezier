#ifndef BEZIER_TEST_DATA_HPP
#define BEZIER_TEST_DATA_HPP

#include <Eigen/Dense>

#include "Bezier/declarations.h"

namespace Bezier
{
namespace TestData
{

// Test curve control points
constexpr std::array<std::pair<double, double>, 4> kCurvePts{
    {{84.000, 162.000}, {246.000, 30.000}, {48.000, 236.000}, {180.000, 110.000}}};

constexpr std::array<std::pair<double, double>, 4> kRootsPts{
    {{-50.000, -50.000}, {75.000, 48.000}, {64.000, 65.000}, {50.000, -50.000}}};

constexpr std::array<std::pair<double, double>, 5> kIntersectPts{
    {{180.000, 110.000}, {175.000, 160.000}, {60.000, 48.000}, {164.000, 165.000}, {124.000, 134.000}}};

// First derivative control points
constexpr std::array<std::pair<double, double>, 3> kExpectedFirstDerivative{
    {{486.000, -396.000}, {-594.000, 618.000}, {396.000, -378.000}}};

// Second derivative control points
constexpr std::array<std::pair<double, double>, 2> kExpectedSecondDerivative{
    {{-2160.000, 2028.000}, {1980.000, -1992.000}}};

// Helper to convert array to PointVector
template <size_t N> inline PointVector toPointVector(const std::array<std::pair<double, double>, N>& pts)
{
  PointVector v;
  v.reserve(N);
  for (const auto& p : pts)
    v.emplace_back(p.first, p.second);
  return v;
}

// Helper to convert array to Eigen matrix
template <size_t N> inline Eigen::MatrixX2d toEigenMatrix(const std::array<std::pair<double, double>, N>& pts)
{
  Eigen::MatrixX2d m(N, 2);
  for (size_t i = 0; i < N; ++i)
  {
    m(static_cast<Eigen::Index>(i), 0) = pts[i].first;
    m(static_cast<Eigen::Index>(i), 1) = pts[i].second;
  }
  return m;
}

} // namespace TestData

// Convenience accessors
inline PointVector curvePointsAsVector() { return TestData::toPointVector(TestData::kCurvePts); }
inline Eigen::MatrixX2d curvePointsAsMatrix() { return TestData::toEigenMatrix(TestData::kCurvePts); }
inline Eigen::MatrixX2d rootPointsAsMatrix() { return TestData::toEigenMatrix(TestData::kRootsPts); }
inline Eigen::MatrixX2d intersectionPointsAsMatrix() { return TestData::toEigenMatrix(TestData::kIntersectPts); }

} // namespace Bezier

#endif // BEZIER_TEST_DATA_HPP
