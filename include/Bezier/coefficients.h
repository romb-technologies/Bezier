#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include "Bezier/utils.h"

#include <unordered_map>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace Bezier
{
namespace Coefficients
{

template <class Func, class... Args> struct lazyFunctor
{
  Func& fun_;
  std::tuple<Args...> args_;
  lazyFunctor(Func& fun, Args... args) : fun_(fun), args_(std::make_tuple(args...)) {}
  template <class Out> operator Out() { return invoke(std::index_sequence_for<Args...>{}); }
  template <std::size_t... I> auto invoke(std::index_sequence<I...>) { return fun_(std::get<I>(args_)...); }
};

inline Eigen::MatrixXd bernstein(unsigned n)
{
  static std::unordered_map<unsigned, Eigen::MatrixXd> cache;
  auto fun = [n]() -> Eigen::MatrixXd {
    Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(n, n);
    coeffs.diagonal(-1).setLinSpaced(-1, -static_cast<int>(n - 1));
    coeffs = coeffs.exp();
    coeffs.array().colwise() *= coeffs.row(n - 1).transpose().array().abs();
    return coeffs;
  };
  return cache.try_emplace(n, lazyFunctor(fun)).first->second;
}

inline Eigen::MatrixXd leftSplit(unsigned n, double t = 0.5)
{
  static std::unordered_map<unsigned, Eigen::MatrixXd> cache;
  auto fun = [n](double t) -> Eigen::MatrixXd {
    return bernstein(n).inverse() * Bezier::Utils::powVector(t, n).asDiagonal() * bernstein(n);
  };
  return t == 0.5 ? cache.try_emplace(n, lazyFunctor(fun, 0.5)).first->second : fun(t);
}

inline Eigen::MatrixXd rightSplit(unsigned n, double t = 0.5)
{
  static std::unordered_map<unsigned, Eigen::MatrixXd> cache;
  auto fun = [n](double t) -> Eigen::MatrixXd {
    Eigen::MatrixXd coeffs = leftSplit(n, t);
    for (unsigned k{}; k < n; k++)
      coeffs.col(n - 1 - k).head(n - k) = coeffs.diagonal(-static_cast<int>(k)).reverse();
    coeffs.triangularView<Eigen::StrictlyLower>().setZero();
    return coeffs;
  };
  return t == 0.5 ? cache.try_emplace(n, lazyFunctor(fun, 0.5)).first->second : fun(t);
}

inline Eigen::MatrixXd raiseOrder(unsigned n)
{
  static std::unordered_map<unsigned, Eigen::MatrixXd> cache;
  auto fun = [n]() -> Eigen::MatrixXd {
    Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(n + 1, n);
    coeffs.diagonal(-1) = coeffs.diagonal().setLinSpaced(1, 1. / n).reverse();
    return coeffs;
  };
  return cache.try_emplace(n, lazyFunctor(fun)).first->second;
}

inline Eigen::MatrixXd lowerOrder(unsigned n)
{
  static std::unordered_map<unsigned, Eigen::MatrixXd> cache;
  auto fun = [n]() -> Eigen::MatrixXd { return raiseOrder(n - 1).completeOrthogonalDecomposition().pseudoInverse(); };
  return cache.try_emplace(n, lazyFunctor(fun)).first->second;
}

} // namespace Coefficients
} // namespace Bezier

#endif // COEFFICIENTS_H
