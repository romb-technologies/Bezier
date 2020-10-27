#ifndef STURM_H
#define STURM_H

#include <Eigen/Dense>
#include <vector>

namespace Sturm
{

enum RootType
{
  All,
  Min,
  Max
};

inline Eigen::MatrixXd chain(const Eigen::VectorXd& poly, double epsilon = 0.001)
{
  Eigen::MatrixXd sturm = Eigen::MatrixXd::Zero(poly.size(), poly.size() + 2);

  sturm.row(0).head(poly.size()) = poly;
  for (uint j = 1; j < poly.size(); j++)
    sturm(1, j) = (poly.size() - j) * sturm(0, j - 1);

  for (uint i = 2; i < poly.size(); i++)
  {
    const Eigen::VectorXd& d2 = sturm.row(i - 2).tail(poly.size() + 2 - i + 2);
    const Eigen::VectorXd& d1 = sturm.row(i - 1).tail(poly.size() + 2 - i + 1);

    if (std::fabs(d1.norm() - std::fabs(d1(d1.size() - 3))) < epsilon)
      return sturm.block(0, 0, i, poly.size());

    if (std::fabs(d1(0)) > epsilon)
    {
      double T, M;
      T = d2(0) / d1(0);
      M = (d2(1) - T * d1(1)) / d1(0);
      for (uint j = 0; j < poly.size() - i; j++)
        sturm(i, i + j) = -(d2(j + 2) - M * d1(j + 1) - T * d1(j + 2));
    }
    else
    {
      auto Trim = [epsilon](const Eigen::VectorXd p) {
        uint k = 0;
        while (std::fabs(p(k)) < epsilon)
          k++;
        return p.tail(p.size() - k);
      };
      auto Inflate = [poly](const Eigen::VectorXd p) {
        Eigen::VectorXd v = Eigen::VectorXd::Zero(poly.size() + 2);
        v.segment(poly.size() - p.size(), p.size()) = p;
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
      sturm.row(i) = Inflate(r);
    }
  }
  return sturm.leftCols(poly.size());
};

inline int interval(const Eigen::MatrixXd& sturm_chain, double t1, double t2)
{
  Eigen::VectorXd power_basis_1 =
      Eigen::pow(t1, Eigen::ArrayXd::LinSpaced(sturm_chain.cols(), sturm_chain.cols() - 1, 0));
  Eigen::VectorXd power_basis_2 =
      Eigen::pow(t2, Eigen::ArrayXd::LinSpaced(sturm_chain.cols(), sturm_chain.cols() - 1, 0));

  Eigen::VectorXd signcount1 = sturm_chain * power_basis_1;
  Eigen::VectorXd signcount2 = sturm_chain * power_basis_2;

  int count1{0}, count2(0);
  for (uint k = 1; k < signcount1.size(); k++)
  {
    if (signcount1(k - 1) * signcount1(k) <= 0 && signcount1(k) != 0)
      count1++;
    if (signcount2(k - 1) * signcount2(k) <= 0 && signcount2(k) != 0)
      count2++;
  }
  return count1 - count2;
};

inline std::vector<double> roots(const Eigen::VectorXd& poly, RootType root_type = All, double epsilon = 0.001)
{
  auto sturm_chain = chain(poly);
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
            Eigen::pow(a, Eigen::ArrayXd::LinSpaced(sturm_chain.cols(), sturm_chain.cols() - 1, 0)).eval();
        double g_a = sturm_chain.row(0) * power_basis_a;
        double g_b = sturm_chain.row(0) * power_basis_b;
        switch (root_type)
        {
        case Min:
          if (g_a <= 0 || g_b > 0)
            flag = true;
          break;
        case Max:
          if (g_a > 0 || g_b <= 0)
            flag = true;
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
