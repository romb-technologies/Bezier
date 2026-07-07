#include "test_data.hpp"
#include "test_oracles.hpp"

#include <future>
#include <vector>

#include <gtest/gtest.h>

#include "Bezier/coefficients.h"
#include "Bezier/utils.h"

namespace Bezier
{

namespace
{

double binomial(unsigned n, unsigned k)
{
  double result = 1.0;
  for (unsigned i = 1; i <= k; i++)
    result *= static_cast<double>(n - k + i) / i;
  return result;
}

PointVector toPoints(const Eigen::MatrixX2d& cp)
{
  PointVector points;
  points.reserve(cp.rows());
  for (Eigen::Index i = 0; i < cp.rows(); i++)
    points.emplace_back(cp(i, 0), cp(i, 1));
  return points;
}

} // namespace

// These are the contracts the planned closed-form coefficient rewrite must
// preserve — asserted against the math, not the current matrix construction.

TEST(CoefficientsTest, BernsteinMatchesBinomialBasis)
{
  // powVector(t, n) * bernstein(n) yields the Bernstein basis functions
  // B_{k,m}(t) = C(m,k) t^k (1-t)^(m-k), with m = n - 1
  for (unsigned n = 2; n <= 8; n++)
  {
    unsigned m = n - 1;
    for (double t : {0.0, 0.2, 0.5, 0.77, 1.0})
    {
      Eigen::RowVectorXd basis = Utils::powVector(t, n) * Coefficients::bernstein(n);
      ASSERT_EQ(basis.size(), n);
      for (unsigned k = 0; k < n; k++)
      {
        double expected = binomial(m, k) * std::pow(t, k) * std::pow(1 - t, m - k);
        EXPECT_NEAR(basis(k), expected, Oracles::kAlgebraic) << "n=" << n << " k=" << k << " t=" << t;
      }
    }
  }
}

TEST(CoefficientsTest, SplitMatricesMatchDeCasteljau)
{
  Eigen::MatrixX2d cp = curvePointsAsMatrix();
  PointVector cp_points = toPoints(cp);
  unsigned n = cp.rows();

  // t = 0.5 exercises the cached path, t = 0.3 the uncached one
  for (double t_split : {0.5, 0.3})
  {
    PointVector left = toPoints(Coefficients::leftSplit(n, t_split) * cp);
    PointVector right = toPoints(Coefficients::rightSplit(n, t_split) * cp);

    for (double s : Oracles::sampleParams(Oracles::kCoarseSamples))
    {
      Point left_expected = Oracles::deCasteljau(cp_points, s * t_split);
      Point left_actual = Oracles::deCasteljau(left, s);
      EXPECT_NEAR(left_actual.x(), left_expected.x(), Oracles::kGeom) << "left t_split=" << t_split << " s=" << s;
      EXPECT_NEAR(left_actual.y(), left_expected.y(), Oracles::kGeom) << "left t_split=" << t_split << " s=" << s;

      Point right_expected = Oracles::deCasteljau(cp_points, t_split + s * (1 - t_split));
      Point right_actual = Oracles::deCasteljau(right, s);
      EXPECT_NEAR(right_actual.x(), right_expected.x(), Oracles::kGeom) << "right t_split=" << t_split << " s=" << s;
      EXPECT_NEAR(right_actual.y(), right_expected.y(), Oracles::kGeom) << "right t_split=" << t_split << " s=" << s;
    }

    // Both pieces share the split point
    EXPECT_NEAR((left.back() - right.front()).norm(), 0.0, Oracles::kGeom);
  }
}

TEST(CoefficientsTest, LowerOrderIsLeftInverseOfRaiseOrder)
{
  for (unsigned n = 2; n <= 8; n++)
  {
    Eigen::MatrixXd product = Coefficients::lowerOrder(n + 1) * Coefficients::raiseOrder(n);
    EXPECT_TRUE(product.isApprox(Eigen::MatrixXd::Identity(n, n), Oracles::kAlgebraic)) << "n=" << n;
  }
}

TEST(CoefficientsTest, CacheIdempotence)
{
  for (unsigned n = 2; n <= 6; n++)
  {
    EXPECT_TRUE(Coefficients::bernstein(n).isApprox(Coefficients::bernstein(n)));
    EXPECT_TRUE(Coefficients::leftSplit(n).isApprox(Coefficients::leftSplit(n)));
    EXPECT_TRUE(Coefficients::rightSplit(n).isApprox(Coefficients::rightSplit(n)));
    EXPECT_TRUE(Coefficients::raiseOrder(n).isApprox(Coefficients::raiseOrder(n)));
    EXPECT_TRUE(Coefficients::lowerOrder(n).isApprox(Coefficients::lowerOrder(n)));
  }
}

TEST(CoefficientsTest, ConcurrentFactoryAccess)
{
  // Hammer the cached factories from parallel tasks; meaningful failures
  // surface under the sanitizer CI job (data races, torn cache writes).
  // kTasks comfortably oversubscribes typical CI core counts to maximize
  // contention; kRepsPerTask repeats enough per task to catch a race without
  // slowing the sanitizer job down; kOrderRange (orders 2..8) cycles every
  // small order so multiple tasks land on the same cache entry concurrently.
  constexpr unsigned kTasks = 16;
  constexpr unsigned kRepsPerTask = 50;
  constexpr unsigned kOrderRange = 7;
  std::vector<std::future<bool>> tasks;
  for (unsigned task = 0; task < kTasks; task++)
    tasks.push_back(std::async(std::launch::async, [task]() {
      for (unsigned rep = 0; rep < kRepsPerTask; rep++)
      {
        unsigned n = 2 + (task + rep) % kOrderRange;
        if (Coefficients::bernstein(n).rows() != n || Coefficients::leftSplit(n).rows() != n ||
            Coefficients::rightSplit(n).rows() != n || Coefficients::raiseOrder(n).rows() != n + 1 ||
            Coefficients::lowerOrder(n).rows() != n - 1)
          return false;
      }
      return true;
    }));
  for (auto& task : tasks)
    EXPECT_TRUE(task.get());
}

} // namespace Bezier
