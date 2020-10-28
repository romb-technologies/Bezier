#ifndef STURM_H
#define STURM_H

#include <Eigen/Dense>
#include <vector>

namespace Sturm
{

/*!
 * \brief Local shape of function at root
 */
enum RootType
{
  All, /*! Find all roots */
  Convex, /*! Find all convex roots */
  Concave /*! Find all concave roots */
};

/*!
 * \brief chain
 * \param polynomial Coefficients of polynomial equation (highest degree first)
 * \param epsilon Precision
 * \return
 */
inline Eigen::MatrixXd chain(const Eigen::VectorXd& polynomial, double epsilon = 0.001)
{
  Eigen::MatrixXd sturm_chain = Eigen::MatrixXd::Zero(polynomial.size(), polynomial.size() + 2);

  sturm_chain.row(0).head(polynomial.size()) = polynomial;
  for (uint j = 1; j < polynomial.size(); j++)
    sturm_chain(1, j) = (polynomial.size() - j) * sturm_chain(0, j - 1);

  for (uint i = 2; i < polynomial.size(); i++)
  {
    const Eigen::VectorXd& d2 = sturm_chain.row(i - 2).tail(polynomial.size() + 2 - i + 2);
    const Eigen::VectorXd& d1 = sturm_chain.row(i - 1).tail(polynomial.size() + 2 - i + 1);

    if (std::fabs(d1.norm() - std::fabs(d1(d1.size() - 3))) < epsilon)
      return sturm_chain.block(0, 0, i, polynomial.size());

    if (std::fabs(d1(0)) > epsilon)
    {
      double T, M;
      T = d2(0) / d1(0);
      M = (d2(1) - T * d1(1)) / d1(0);
      for (uint j = 0; j < polynomial.size() - i; j++)
        sturm_chain(i, i + j) = -(d2(j + 2) - M * d1(j + 1) - T * d1(j + 2));
    }
    else
    {
      auto Trim = [epsilon](const Eigen::VectorXd p) {
        uint k = 0;
        while (std::fabs(p(k)) < epsilon)
          k++;
        return p.tail(p.size() - k);
      };
      auto Inflate = [polynomial](const Eigen::VectorXd p) {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(polynomial.size() + 2);
        v.segment(polynomial.size() - p.size(), p.size()) = p;
        return v;
      };

      Eigen::VectorXd a = Trim(d2.head(d2.size() - 2));
      Eigen::VectorXd b = Trim(d1.head(d1.size() - 2));

      Eigen::VectorXd r = a;
      while (r.size() && r.size() >= b.size())
      {
        double L = r(0) / b(0);
        uint X = static_cast<uint>(r.size() - b.size());

        for (uint k = 0; k < b.size(); k++)
          r(r.size() - (X + b.size()) + k) -= L * b(k);
        r = Trim(r);
      }
      sturm_chain.row(i) = Inflate(r);
    }
  }
  return sturm_chain.leftCols(polynomial.size());
};

/*!
 * \brief Counts the number of roots in a given interval
 * \param sturm_chain Series of Sturm functions
 * \param t1 start of interval
 * \param t2 end of interval
 * \return Number of roots
 */
inline int interval(const Eigen::MatrixXd& sturm_chain, double t1, double t2)
{
  Eigen::VectorXd power_basis_1 =
      Eigen::pow(t1, Eigen::ArrayXd::LinSpaced(sturm_chain.cols(), sturm_chain.cols() - 1, 0));
  Eigen::VectorXd power_basis_2 =
      Eigen::pow(t2, Eigen::ArrayXd::LinSpaced(sturm_chain.cols(), sturm_chain.cols() - 1, 0));

  Eigen::VectorXd signcount1 = sturm_chain * power_basis_1;
  Eigen::VectorXd signcount2 = sturm_chain * power_basis_2;

  int count1{0}, count2{0};
  for (uint k = 1; k < signcount1.size(); k++)
  {
    if (signcount1(k - 1) * signcount1(k) <= 0 && signcount1(k) != 0)
      count1++;
    if (signcount2(k - 1) * signcount2(k) <= 0 && signcount2(k) != 0)
      count2++;
  }
  return count1 - count2;
};

/*!
 * \brief Find all roots in the [0, 1] interval
 * \param polynomial Coefficients of polynomial equation (highest degree first)
 * \param root_type Type of roots to find
 * \param epsilon Precision of resulting roots
 * \return A vector of roots
 */
inline std::vector<double> roots(const Eigen::VectorXd& polynomial, RootType root_type = All, double epsilon = 0.001)
{
  auto sturm_chain = chain(polynomial);
  std::vector<std::tuple<double, double, bool>> intervals;
  std::vector<double> roots;

  auto iterate = [&](std::tuple<double, double, bool> item) {
    double a = std::get<0>(item);
    double b = std::get<1>(item);
    bool flag = std::get<2>(item);
    double a_b = (a + b) / 2;
    uint root_num = static_cast<uint>(interval(sturm_chain, std::get<0>(item), std::get<1>(item)));
    switch (root_num)
    {
    case 0:
      return;
    case 1:
      if (a_b - a < epsilon)
      {
        roots.emplace_back((a + a_b) / 2);
        return;
      }
      else if (root_type != All && !flag)
      {
        Eigen::VectorXd power_basis_a =
            Eigen::pow(a, Eigen::ArrayXd::LinSpaced(sturm_chain.cols(), sturm_chain.cols() - 1, 0)).eval();
        Eigen::VectorXd power_basis_b =
            Eigen::pow(b, Eigen::ArrayXd::LinSpaced(sturm_chain.cols(), sturm_chain.cols() - 1, 0)).eval();
        double g_a = sturm_chain.row(0) * power_basis_a;
        double g_b = sturm_chain.row(0) * power_basis_b;
        switch (root_type)
        {
        case Convex:
          if (g_a <= 0 && g_b > 0)
            flag = true;
          else
            return;

          break;
        case Concave:
          if (g_a > 0 && g_b <= 0)
            flag = true;
          else
            return;
          break;
        case All:;
        }
      }
      [[clang::fallthrough]];
    default:
      intervals.emplace_back(a, a_b, flag);
      intervals.emplace_back(a_b, b, flag);
    }
  };

  iterate({0.0, 1.0, false});

  while (!intervals.empty())
  {
    double a, b, a_b;
    bool flag;
    std::tie(a, b, flag) = intervals.back();
    intervals.pop_back();
    a_b = (a + b) / 2;

    iterate({a, a_b, flag});
    iterate({a_b, b, flag});
  }

  return roots;
}
} // namespace Sturm
#endif // STURM_H
