#include "Bezier/utils.h"

#include <unsupported/Eigen/Polynomials>

using namespace Bezier;
namespace bu = Bezier::Utils;

std::vector<double> Bezier::Utils::solvePolynomial(const Eigen::VectorXd& polynomial)
{
  // Trim trailing zero coefficients from the polynomial
  unsigned idx = polynomial.size();
  while (idx && std::fabs(polynomial(idx - 1)) < bu::epsilon)
    --idx;
  if (idx < 2) // Polynomial is a constant
    return {};

  struct PolynomialRoots : public std::vector<double>
  {
    PolynomialRoots(unsigned size) : std::vector<double>() { reserve(size); }
    void push_back(double t) // only allow valid roots
    {
      if (t >= 0 && t <= 1)
        std::vector<double>::push_back(t);
    }
  } roots(idx);
  Eigen::PolynomialSolver<double, Eigen::Dynamic>(polynomial.head(idx)).realRoots(roots);
  return roots;
}
