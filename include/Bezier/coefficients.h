/*
 * Copyright 2024 Mirko Kokot
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H

#include "Bezier/utils.h"

#include <mutex>
#include <unordered_map>

namespace Bezier
{
namespace Coefficients
{

template <class Func> struct lazyFunctor
{
  Func& fun_;
  template <class Out> operator Out() { return fun_(); }
};

template <class Compute> inline const Eigen::MatrixXd& memoize(unsigned n, Compute compute)
{
  static std::unordered_map<unsigned, Eigen::MatrixXd> cache;
  static std::mutex cache_mutex;

  std::lock_guard<std::mutex> lock(cache_mutex);
  return cache.try_emplace(n, lazyFunctor<Compute>{compute}).first->second;
}

inline const Eigen::MatrixXd& bernstein(unsigned n)
{
  return memoize(n, [n] {
    Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(n, n);
    double rn = 1.0; // C(n-1, k) row factor
    for (unsigned k = 0; k < n; ++k)
    {
      double c = 1.0; // C(k,0)
      for (unsigned i = 0; i <= k; ++i)
      {
        coeffs(k, i) = ((k - i) & 1u ? -1.0 : 1.0) * c * rn;
        c = c * (k - i) / (i + 1); // C(k,i+1)
      }
      rn = rn * (n - 1 - k) / (k + 1); // C(n-1, k+1)
    }
    return coeffs;
  });
}

inline Eigen::MatrixXd leftSplit(unsigned n, double t)
{
  Eigen::RowVectorXd tp = Utils::powVector(t, n);     // [t^0 .. t^{n-1}]
  Eigen::RowVectorXd ct = Utils::powVector(1 - t, n); // [(1-t)^0 .. (1-t)^{n-1}]
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(n, n);
  for (unsigned k = 0; k < n; ++k)
  {
    double c = 1.0; // C(k,0)
    for (unsigned i = 0; i <= k; ++i)
    {
      L(k, i) = c * tp(i) * ct(k - i);
      c = c * (k - i) / (i + 1); // C(k,i+1)
    }
  }
  return L;
}

inline const Eigen::MatrixXd& leftSplit(unsigned n)
{
  return memoize(n, [n] { return leftSplit(n, 0.5); });
}

inline Eigen::MatrixXd rightSplit(unsigned n, double t) { return leftSplit(n, 1.0 - t).reverse(); }

inline const Eigen::MatrixXd& rightSplit(unsigned n)
{
  return memoize(n, [n] { return leftSplit(n).reverse().eval(); });
}

inline const Eigen::MatrixXd& raiseOrder(unsigned n)
{
  return memoize(n, [n] {
    Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(n + 1, n);
    coeffs.diagonal(-1) = coeffs.diagonal().setLinSpaced(1, 1. / n).reverse();
    return coeffs;
  });
}

inline const Eigen::MatrixXd& lowerOrder(unsigned n)
{
  return memoize(n, [n] { return raiseOrder(n - 1).completeOrthogonalDecomposition().pseudoInverse().eval(); });
}

} // namespace Coefficients
} // namespace Bezier

#endif // COEFFICIENTS_H
