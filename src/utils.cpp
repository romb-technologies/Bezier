#include "Bezier/utils.h"

#include <limits>
#include <numeric>

#include <unsupported/Eigen/Polynomials>

using namespace Bezier;
namespace bu = Bezier::Utils;

std::vector<unsigned> Bezier::Utils::visvalingamWyatt(const PointVector& polyline)
{
  // Vector of indices sorted by contribution to the polyline shape (first and last contribute the most by default)
  // Initialized with all indices, taking care to put the first and the last at the start
  std::vector<unsigned> by_contribution(polyline.size());
  std::iota(by_contribution.begin(), by_contribution.end(), -1);
  by_contribution.front() += polyline.size();

  // Helper structure to keep track of each Points neigbours and contribution to the polyline shape
  struct Vertex
  {
    size_t prev, next;
    double contribution;
  };
  std::vector<Vertex> vertices(polyline.size());

  // Visvalingam-Whyatt measures contribution as an area between 3 consecutive Points
  auto area = [&polyline](unsigned idx1, unsigned idx2, unsigned idx3) {
    return std::fabs(bu::cross(polyline[idx2] - polyline[idx1], polyline[idx3] - polyline[idx1])) / 2;
  };
  auto cmp = [&vertices](unsigned idx1, unsigned idx2) {
    return vertices[idx1].contribution < vertices[idx2].contribution;
  };

  // Initialize vertices
  constexpr size_t NaN = std::numeric_limits<size_t>::quiet_NaN();
  vertices.front() = {NaN, 1, 0.0};
  for (unsigned k = 1; k + 1 < polyline.size(); k++)
    vertices[k] = {k - 1, k + 1, area(k - 1, k, k + 1)};
  vertices.back() = {polyline.size() - 2, NaN, 0.0};

  // Smallest contribution will be at the end of the vector
  for (auto it = by_contribution.rbegin(); it != by_contribution.rend() - 2; ++it)
  {
    // Select and move a Point with smallest current contribution
    std::iter_swap(it, std::min_element(it, by_contribution.rend() - 2, cmp));

    // Update previous and next Vertices (neighbours and contributions)
    auto prev = vertices[*it].prev;
    auto next = vertices[*it].next;
    vertices[prev].next = next;
    vertices[next].prev = prev;
    vertices[prev].contribution = area(vertices[prev].prev, prev, next);
    vertices[next].contribution = area(prev, next, vertices[next].next);
  }

  return by_contribution;
}

PointVector Bezier::Utils::polylineSimplified(const PointVector& polyline, unsigned int N)
{
  if (polyline.size() < N)
    return polyline;
  if (N == 2)
    return std::vector{polyline.front(), polyline.back()};

  std::vector<Point> simplified(N);
  auto vw = visvalingamWyatt(polyline);
  std::sort(vw.begin(), vw.begin() + N);
  std::transform(vw.begin(), vw.begin() + N, simplified.begin(), [&polyline](unsigned k) { return polyline[k]; });
  return simplified;
}

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
