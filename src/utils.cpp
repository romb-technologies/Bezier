#include "Bezier/utils.h"

#include <numeric>

#include <unsupported/Eigen/Polynomials>

using namespace Bezier;
namespace bu = Bezier::Utils;

std::vector<unsigned> Bezier::Utils::visvalingamWyatt(const PointVector& polyline)
{
  // Vector of indices sorted by contribution to the polyline shape
  std::vector<unsigned> by_contribution;
  by_contribution.reserve(polyline.size());

  // Helper structure to keep track of each Points neigbours and contribution to the polyline shape
  struct Vertex
  {
    size_t prev, next;
    double contribution;
  };
  std::vector<Vertex> vertices;
  vertices.reserve(polyline.size());

  // Set is used to keep track of sorted contributions
  auto cmp = [&vertices](unsigned idx1, unsigned idx2) {
    return vertices[idx1].contribution > vertices[idx2].contribution;
  };
  std::vector<unsigned> min_heap(polyline.size() - 2);
  std::iota(min_heap.begin(), min_heap.end(), 1);

  // Visvalingam-Whyatt measures contribution as an area between 3 consecutive Points
  auto area = [&polyline](unsigned id1, unsigned id2, unsigned id3) {
    const auto& A = polyline[id1];
    const auto& B = polyline[id2];
    const auto& C = polyline[id3];
    return std::fabs((B.x() - A.x()) * (C.y() - A.y()) - (C.x() - A.x()) * (B.y() - A.y()));
  };

  // Initialize structures used for the algorithm
  vertices.push_back({0, 1, 0.0});
  for (unsigned k = 1; k + 1 < polyline.size(); k++)
    vertices.push_back({k - 1, k + 1, area(k - 1, k, k + 1)});
  vertices.push_back({polyline.size() - 2, polyline.size(), 0.0});

  while (!min_heap.empty())
  {
    // Select and erase a Point with smallest current contribution to polyline shape
    std::make_heap(min_heap.begin(), min_heap.end(), cmp);
    std::pop_heap(min_heap.begin(), min_heap.end(), cmp);
    by_contribution.push_back(min_heap.back());
    min_heap.pop_back();

    // Update previous and next Vertex:
    // - update neighbours neighbours
    // - update neighbours contribution
    auto prev = vertices[by_contribution.back()].prev;
    auto next = vertices[by_contribution.back()].next;
    vertices[prev].next = vertices[by_contribution.back()].next;
    vertices[next].prev = vertices[by_contribution.back()].prev;
    vertices[prev].contribution = area(vertices[prev].prev, prev, vertices[prev].next);
    vertices[next].contribution = area(vertices[next].prev, next, vertices[next].next);
  }

  by_contribution.push_back(0);
  by_contribution.push_back(polyline.size() - 1);
  std::reverse(by_contribution.begin(), by_contribution.end());

  return by_contribution;
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
