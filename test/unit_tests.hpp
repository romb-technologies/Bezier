#ifndef BEZIER_UNIT_TESTS_HPP
#define BEZIER_UNIT_TESTS_HPP

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

// Expected values for validation (full precision from actual computation)
constexpr std::array<std::pair<double, double>, 5> kExpectedValueAt{{{83.999999999999986, 161.99999999999997},
                                                                     {148.78124999999997, 115.90625000000001},
                                                                     {143.24999999999997, 133.75000000000003},
                                                                     {132.09374999999991, 152.71875000000006},
                                                                     {180.00000000000003, 110.00000000000007}}};

// Expected length of the test curve (current algorithm precision)
constexpr double kExpectedCurveLength = 189.086275311931;

// Curvature values at specific t values
constexpr std::array<double, 5> kExpectedCurvature{0.00052864283054623, 0.13400265200820813, 0.00339166466643027,
                                                   -0.83667680576075631, -0.00024618707208090};

// Curvature derivative values at specific t values
constexpr std::array<double, 5> kExpectedCurvatureDerivative{
    0.00620125036871112, 6.26701004231499148, -0.06892991355501528, 86.75145336780214222, 0.00362325743277114};

// Tangent vectors at specific t values
constexpr std::array<std::pair<double, double>, 5> kExpectedTangent{{{0.77523498551728931, -0.63167295116223576},
                                                                     {0.98169156767741872, -0.19047746835532003},
                                                                     {-0.55219905682528458, 0.83371230148131203},
                                                                     {0.98328200498446017, -0.18208926018230745},
                                                                     {0.72335554414357217, -0.69047574668250067}}};

// Normal vectors at specific t values
constexpr std::array<std::pair<double, double>, 5> kExpectedNormal{{{0.63167295116223576, 0.77523498551728931},
                                                                    {0.19047746835532003, 0.98169156767741872},
                                                                    {-0.83371230148131203, -0.55219905682528458},
                                                                    {0.18208926018230745, 0.98328200498446017},
                                                                    {0.69047574668250067, 0.72335554414357217}}};

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
inline PointVector expectedValueAtVector() { return TestData::toPointVector(TestData::kExpectedValueAt); }
inline Eigen::MatrixX2d expectedValueAtMatrix() { return TestData::toEigenMatrix(TestData::kExpectedValueAt); }
inline PointVector expectedPolylinePoints()
{
  // clang-format off
  return { {84, 162}, {98.153869628906222, 150.59478759765619}, {110.32470703124994, 141.04736328124991}, {120.63885498046869, 133.23504638671869}, {129.22265624999991, 127.0351562499999}, {136.20245361328119, 122.32501220703119}, {141.70458984374991, 118.98193359374991}, {145.85540771484366, 116.88323974609368}, {148.78124999999989, 115.90624999999989}, {149.82428741455072, 115.80005645751945}, {150.60845947265616, 115.92828369140616}, {151.14955902099604, 116.27559661865226}, {151.46337890624989, 116.82666015624989}, {151.47235107421869, 118.47869873046868}, {150.76171874999989, 120.7617187499999}, {149.45782470703114, 123.55303955078119}, {147.68701171874986, 126.72998046874991}, {143.24999999999986, 133.74999999999991}, {134.33203124999994, 147.01953124999997}, {131.87255859374994, 151.30615234374994}, {131.58489990234369, 152.43304443359364}, {131.73186492919911, 152.68871307373033}, {132.09374999999986, 152.71874999999983}, {133.52545166015619, 152.04058837890619}, {136.00634765624991, 150.27587890624991}, {144.62109374999991, 142.99609374999991}, {158.94873046874994, 129.89794921874991}, {179.99999999999989, 109.99999999999989} };
  // clang-format on
}

} // namespace Bezier

#endif // BEZIER_UNIT_TESTS_HPP
