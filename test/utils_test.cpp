#include "test_data.hpp"
#include "test_oracles.hpp"

#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>

#include "Bezier/utils.h"

namespace Bezier
{

TEST(UtilsTest, PowMatchesStdPow)
{
  for (double base : {0.0, 0.5, 1.0, -1.7, 3.0})
    for (unsigned exp = 0; exp <= 10; exp++)
      EXPECT_NEAR(Utils::pow(base, exp), std::pow(base, exp), 1e-12 * std::fabs(std::pow(base, exp)))
          << "base=" << base << " exp=" << exp;

  // Integer instantiation is exact
  EXPECT_EQ(Utils::pow(3, 7u), 2187);
}

TEST(UtilsTest, Exp2)
{
  for (unsigned exp = 0; exp <= 16; exp++)
    EXPECT_EQ(Utils::exp2(exp), 1u << exp);
}

TEST(UtilsTest, PowVectorMatchesPlainLoop)
{
  // Regression net for the planned replacement of the stateful NullaryExpr
  for (double base : {0.0, 0.37, 1.0, -2.5})
    for (unsigned exp = 1; exp <= 12; exp++)
    {
      Eigen::RowVectorXd actual = Utils::powVector(base, exp);
      ASSERT_EQ(actual.size(), exp);
      double power = 1.0;
      for (unsigned k = 0; k < exp; k++, power *= base)
        EXPECT_DOUBLE_EQ(actual(k), power) << "base=" << base << " exp=" << exp << " k=" << k;
    }
}

TEST(UtilsTest, PowMatrixMatchesPlainLoop)
{
  Eigen::VectorXd bases(4);
  bases << 0.0, 0.25, 1.0, -1.5;
  const unsigned exp = 8;
  Eigen::MatrixXd actual = Utils::powMatrix(bases, exp);
  ASSERT_EQ(actual.rows(), bases.size());
  ASSERT_EQ(actual.cols(), exp);
  for (unsigned i = 0; i < bases.size(); i++)
  {
    double power = 1.0;
    for (unsigned k = 0; k < exp; k++, power *= bases(i))
      EXPECT_DOUBLE_EQ(actual(i, k), power) << "row=" << i << " col=" << k;
  }
}

TEST(UtilsTest, EvaluateChebyshevMatchesNaiveSum)
{
  // T_k evaluated via the trigonometric definition; argument mapped from [0,1]
  Eigen::VectorXd coeffs(6);
  coeffs << 1.5, -0.3, 0.7, 0.05, -0.2, 0.11;
  for (double t : {0.0, 0.1, 0.35, 0.5, 0.82, 1.0})
  {
    double x = 2 * t - 1;
    double expected{};
    for (unsigned k = 0; k < coeffs.size(); k++)
      expected += coeffs(k) * std::cos(k * std::acos(std::clamp(x, -1.0, 1.0)));
    EXPECT_NEAR(Utils::evaluateChebyshev(t, coeffs), expected, 1e-12) << "t=" << t;
  }
}

TEST(UtilsTest, DistOverloads)
{
  // point - point
  EXPECT_DOUBLE_EQ(Utils::dist(Point{0, 0}, Point{3, 4}), 5.0);

  // point - segment: perpendicular foot inside the segment
  EXPECT_DOUBLE_EQ(Utils::dist(Point{0, 0}, Point{10, 0}, Point{5, 2}), 2.0);
  // clamped to segment start / end
  EXPECT_DOUBLE_EQ(Utils::dist(Point{0, 0}, Point{10, 0}, Point{-3, 4}), 5.0);
  EXPECT_DOUBLE_EQ(Utils::dist(Point{0, 0}, Point{10, 0}, Point{13, 4}), 5.0);

  // point - polyline: nearest of the segments wins
  PointVector polyline{{0, 0}, {10, 0}, {10, 10}};
  EXPECT_DOUBLE_EQ(Utils::dist(polyline, Point{5, 1}), 1.0);
  EXPECT_DOUBLE_EQ(Utils::dist(polyline, Point{12, 5}), 2.0);
}

TEST(UtilsTest, PolylineLength)
{
  PointVector polyline{{0, 0}, {3, 4}, {3, 14}};
  EXPECT_DOUBLE_EQ(Utils::polylineLength(polyline), 15.0);
  EXPECT_DOUBLE_EQ(Utils::polylineLength(PointVector{{1, 1}}), 0.0);
  EXPECT_DOUBLE_EQ(Utils::polylineLength(PointVector{}), 0.0);
}

TEST(UtilsTest, Concatenate)
{
  std::vector<int> a{1, 2, 3}, b{4, 5};
  std::vector<int> result = Utils::concatenate(std::move(a), std::move(b));
  EXPECT_EQ(result, (std::vector<int>{1, 2, 3, 4, 5}));
}

TEST(UtilsTest, MaxDeviationIsUpperBound)
{
  // maxDeviation must bound the true deviation of the curve from its chord
  for (unsigned order : {3u, 4u, 5u})
  {
    Eigen::MatrixX2d cp(order + 1, 2);
    for (unsigned i = 0; i <= order; i++)
      cp.row(i) << 100 * std::cos(1.3 * i + order), 100 * std::sin(2.1 * i);

    PointVector cp_vec;
    for (unsigned i = 0; i <= order; i++)
      cp_vec.emplace_back(cp(i, 0), cp(i, 1));

    double bound = Utils::maxDeviation(cp);
    double true_deviation{};
    for (double t{}; t <= 1.0; t += 0.001)
      true_deviation =
          std::max(true_deviation, Utils::dist(cp_vec.front(), cp_vec.back(), Oracles::deCasteljau(cp_vec, t)));
    EXPECT_GE(bound, true_deviation - 1e-12) << "order=" << order;
  }

  // Two points: the "curve" is its own chord
  Eigen::MatrixX2d segment(2, 2);
  segment << 0, 0, 5, 5;
  EXPECT_DOUBLE_EQ(Utils::maxDeviation(segment), 0.0);
}

TEST(UtilsTest, VisvalingamWyattOrdering)
{
  // Zig-zag polyline with one nearly collinear vertex (index 3)
  PointVector polyline{{0, 0}, {1, 5}, {2, -3}, {3, 0.001}, {4, 6}, {5, -2}, {6, 0}};
  std::vector<unsigned> order = Utils::visvalingamWyatt(polyline);

  // Output is a permutation of all indices
  ASSERT_EQ(order.size(), polyline.size());
  std::vector<unsigned> sorted = order;
  std::sort(sorted.begin(), sorted.end());
  for (unsigned k = 0; k < sorted.size(); k++)
    EXPECT_EQ(sorted[k], k) << "not a permutation";

  // Endpoints contribute the most: first two entries are {n-1, 0}
  EXPECT_EQ(order[0], polyline.size() - 1);
  EXPECT_EQ(order[1], 0u);

  // The nearly collinear vertex is eliminated first (appears last)
  EXPECT_EQ(order.back(), 3u);
}

TEST(UtilsTest, PolylineSimplified)
{
  PointVector polyline{{0, 0}, {1, 5}, {2, -3}, {3, 0.001}, {4, 6}, {5, -2}, {6, 0}};

  for (unsigned N : {2u, 3u, 5u})
  {
    PointVector simplified = Utils::polylineSimplified(polyline, N);
    ASSERT_EQ(simplified.size(), N) << "N=" << N;

    // Endpoints preserved
    EXPECT_EQ(simplified.front(), polyline.front());
    EXPECT_EQ(simplified.back(), polyline.back());

    // Result is an ordered subset of the input
    size_t search_from = 0;
    for (const Point& p : simplified)
    {
      auto it = std::find_if(polyline.begin() + search_from, polyline.end(),
                             [&p](const Point& q) { return q == p; });
      ASSERT_NE(it, polyline.end()) << "simplified point not found in order";
      search_from = (it - polyline.begin()) + 1;
    }
  }

  // N >= size returns the input unchanged
  EXPECT_EQ(Utils::polylineSimplified(polyline, 7), polyline);
  EXPECT_EQ(Utils::polylineSimplified(polyline, 100), polyline);
}

TEST(UtilsTest, SolvePolynomial)
{
  // Coefficients in increasing powers of x

  // (x - 0.3)(x - 0.7) = 0.21 - x + x^2 -> both roots inside [0, 1]
  Eigen::VectorXd quadratic(3);
  quadratic << 0.21, -1.0, 1.0;
  std::vector<double> roots = Utils::solvePolynomial(quadratic);
  ASSERT_EQ(roots.size(), 2u);
  std::sort(roots.begin(), roots.end());
  EXPECT_NEAR(roots[0], 0.3, 1e-12);
  EXPECT_NEAR(roots[1], 0.7, 1e-12);

  // (x - 0.5)(x - 2) = 1 - 2.5x + x^2 -> the root outside [0, 1] is excluded
  Eigen::VectorXd with_outside(3);
  with_outside << 1.0, -2.5, 1.0;
  roots = Utils::solvePolynomial(with_outside);
  ASSERT_EQ(roots.size(), 1u);
  EXPECT_NEAR(roots[0], 0.5, 1e-12);

  // Constant polynomials have no roots
  Eigen::VectorXd constant(1);
  constant << 5.0;
  EXPECT_TRUE(Utils::solvePolynomial(constant).empty());

  // Trailing zero coefficients are trimmed before solving
  Eigen::VectorXd padded(5);
  padded << 0.21, -1.0, 1.0, 0.0, 0.0;
  roots = Utils::solvePolynomial(padded);
  ASSERT_EQ(roots.size(), 2u);
  std::sort(roots.begin(), roots.end());
  EXPECT_NEAR(roots[0], 0.3, 1e-12);
  EXPECT_NEAR(roots[1], 0.7, 1e-12);

  // A constant padded with zeros is still a constant
  Eigen::VectorXd padded_constant(3);
  padded_constant << 5.0, 0.0, 0.0;
  EXPECT_TRUE(Utils::solvePolynomial(padded_constant).empty());
}

} // namespace Bezier
