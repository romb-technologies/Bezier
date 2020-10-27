#ifndef STURM_H
#define STURM_H

#include <Eigen/Dense>
#include <vector>

inline Eigen::MatrixXd SturmChain(const Eigen::VectorXd& poly, double epsilon = 0.001)
{
  Eigen::MatrixXd sturm = Eigen::MatrixXd::Zero(poly.size(), poly.size() + 2);

  sturm.row(0).head(poly.size()) = poly;
  for (uint j = 1; j < poly.size(); j++)
    sturm(1, j) = (poly.size() - j) * sturm(0, j - 1);

  for (uint i = 2; i < poly.size(); i++)
  {
    const Eigen::VectorXd& d2 = sturm.row(i - 2).tail(poly.size() + 2 - i + 2);
    const Eigen::VectorXd& d1 = sturm.row(i - 1).tail(poly.size() + 2 - i + 1);

    if(std::fabs(d1.norm() - std::fabs(d1(d1.size()-3))) < epsilon)
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

inline int SturmInterval(const Eigen::MatrixXd& sturm_chain, double t1, double t2)
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

inline std::vector<double> SturmRoots(const Eigen::VectorXd& poly, double epsilon = 0.001)
{
  auto sturm_chain = SturmChain(poly);
  std::vector<std::tuple<double, double, uint>> intervals;
  std::vector<std::tuple<double, double, double>> root_candidates;
  std::vector<double> roots;
  uint temp_res = static_cast<uint>(SturmInterval(sturm_chain, 0.0, 1.0));

  switch (temp_res)
  {
  case 0:
    return std::vector<double>();
  case 1:
    root_candidates.emplace_back(0.0, 1.0, 0.5);
    break;
  default:
    intervals.reserve(temp_res);
    root_candidates.reserve(temp_res);
    roots.reserve(temp_res);
    intervals.emplace_back(0.0, 1.0, temp_res);
  }

  while (!intervals.empty())
  {
    double a = std::get<0>(intervals.back());
    double b = std::get<1>(intervals.back());
    double a_b = (a + b) / 2;
    uint count = std::get<2>(intervals.back());
    intervals.pop_back();

    uint left = static_cast<uint>(SturmInterval(sturm_chain, a, a_b));
    uint right = count - left;

    switch (left)
    {
    case 0:
      break;
    case 1:
      root_candidates.emplace_back(a, a_b, (a + a_b) / 2);
      break;
    default:
      intervals.emplace_back(a, a_b, left);
    }

    switch (right)
    {
    case 0:
      break;
    case 1:
      root_candidates.emplace_back(a_b, b, (a_b + b) / 2);
      break;
    default:
      intervals.emplace_back(a_b, b, right);
    }
  }

  // fine search
  Eigen::Matrix<double, 3, Eigen::Dynamic> derivates;
  derivates.resize(3, poly.size());
  derivates.setZero();
  derivates.topRows(2) = sturm_chain.topRows(2);
  for (uint j = 2; j < poly.size(); j++)
    derivates(2, j) = (poly.size() - j) * derivates(1, j - 1);
  for (auto candidate : root_candidates)
  {
    double a, b, t;
    std::tie(a, b, t) = candidate;
    Eigen::VectorXd power_basis = Eigen::pow(t, Eigen::ArrayXd::LinSpaced(poly.size(), poly.size() - 1, 0));
    Eigen::Vector3d f_fd_fdd = derivates * power_basis;
    while (std::fabs(f_fd_fdd(0)) > epsilon)
    {
      double t_old = t;
      t -= 2 * f_fd_fdd(0) * f_fd_fdd(1) / (2 * f_fd_fdd(1) * f_fd_fdd(1) - f_fd_fdd(0) * f_fd_fdd(2));
      if (t < a && t_old < a)
        t = b;
      if (t > b && t_old > b)
        t = a;
      power_basis = Eigen::pow(t, Eigen::ArrayXd::LinSpaced(poly.size(), poly.size() - 1, 0));
      f_fd_fdd = derivates * power_basis;
    }
    roots.emplace_back(t);
  }

  return roots;
}

#endif // STURM_H
