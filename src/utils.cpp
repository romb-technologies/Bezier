#include "Bezier/utils.h"

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

  // Visvalingam-Whyatt measures contribution as an area between 3 consecutive Points
  auto area = [&polyline](unsigned id1, unsigned id2, unsigned id3) {
    const auto& A = polyline[id1];
    const auto& B = polyline[id2];
    const auto& C = polyline[id3];
    return std::fabs((B.x() - A.x()) * (C.y() - A.y()) - (C.x() - A.x()) * (B.y() - A.y()));
  };

  // Initialize vertices
  std::vector<Vertex> vertices(polyline.size());
  vertices.front() = {0, 1, 0.0};
  for (unsigned k = 1; k + 1 < polyline.size(); k++)
    vertices[k] = {k - 1, k + 1, area(k - 1, k, k + 1)};
  vertices.back() = {polyline.size() - 2, polyline.size(), 0.0};

  // Smallest contribution will be at the end of the vector
  for (auto it = by_contribution.rbegin(); it != by_contribution.rend() - 2; ++it)
  {
    // Select and move a Point with smallest current contribution
    std::iter_swap(it, std::min_element(it, by_contribution.rend() - 2, [&vertices](unsigned idx1, unsigned idx2) {
                     return vertices[idx1].contribution < vertices[idx2].contribution;
                   }));

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

PointVector Bezier::Utils::polylineSimplify(const PointVector &polyline, unsigned int N)
{
  if (polyline.size() < 2)
    throw std::logic_error{"Polyline must have at least two points."};
  if (polyline.size() < N)
    return PointVector(polyline);
  if (N == 2)
    return PointVector{polyline.front(), polyline.back()};

  auto by_contribution = visvalingamWyatt(polyline);
  std::sort(by_contribution.begin(), by_contribution.begin() + N);
  PointVector simplified(N);
  for (size_t k{0}; k < N; k++)
    simplified[k] = polyline[by_contribution[k]];
  return simplified;
}

std::vector<double> Bezier::Utils::solvePolynomial(const Eigen::VectorXd& polynomial)
{
  // Trim leading zeros from the polynomial
  auto idx = polynomial.size();
  while (idx && std::abs(polynomial(idx - 1)) < bu::epsilon)
    --idx;
  // Polynomial is constant
  if (idx < 2)
    return {};

  struct PolynomialRoots : public std::vector<double>
  {
    void push_back(double t) // only allow valid roots
    {
      if (t >= 0 && t <= 1)
        std::vector<double>::push_back(t);
    }
  } roots;
  roots.reserve(idx);
  Eigen::PolynomialSolver<double, Eigen::Dynamic>(polynomial.head(idx)).realRoots(roots);
  return roots;
}
