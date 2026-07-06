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
    {{84.0, 162.0}, {246.0, 30.0}, {48.0, 236.0}, {180.0, 110.0}}};

constexpr std::array<std::pair<double, double>, 4> kRootsPts{
    {{-50.0, -50.0}, {75.0, 48.0}, {64.0, 65.0}, {50.0, -50.0}}};

constexpr std::array<std::pair<double, double>, 5> kIntersectPts{
    {{180.0, 110.0}, {175.0, 160.0}, {60.0, 48.0}, {164.0, 165.0}, {124.0, 134.0}}};

// First derivative control points
constexpr std::array<std::pair<double, double>, 3> kExpectedFirstDerivative{
    {{486.0, -396.0}, {-594.0, 618.0}, {396.0, -378.0}}};

// Second derivative control points
constexpr std::array<std::pair<double, double>, 2> kExpectedSecondDerivative{
    {{-2160.0, 2028.0}, {1980.0, -1992.0}}};

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
