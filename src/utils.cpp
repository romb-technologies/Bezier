#include "Bezier/utils.h"

#include <numeric>
#include <iostream>

using namespace Bezier;

PointVector Bezier::_polylineSimplify(const PointVector &polyline, unsigned int N)
{
  if (polyline.size() < 2)
    throw std::logic_error{"Polyline must have at least two points."};
  if (polyline.size() < N)
    return PointVector(polyline);
  if (N == 2)
    return PointVector{polyline.front(), polyline.back()};

  // Simplification is done using the Visvalingam-Whyatt algorithm
  // Helper structure to keep track of each Points neigbours and contribution to the polyline shape
  struct Vertex {
    size_t prev, next;
    double contribution;
  };
  std::vector<Vertex> vertices;
  vertices.reserve(polyline.size());

  // Set is used to keep track of sorted contributions
  auto cmp=[&vertices](size_t idx1, size_t idx2) {
    return vertices[idx1].contribution > vertices[idx2].contribution;
  };
  std::vector<size_t> by_contribution(polyline.size()-2);
  std::iota(by_contribution.begin(), by_contribution.end(), 1);

  // Visvalingam-Whyatt measures contribution as an area between 3 consecutive Points
  auto area = [&polyline](size_t id1, size_t id2, size_t id3)
  {
    const auto& A = polyline[id1];
    const auto& B = polyline[id2];
    const auto& C = polyline[id3];
    return std::fabs((B.x() - A.x()) * (C.y() - A.y()) - (C.x() - A.x()) * (B.y() - A.y()));
  };

  // Initialize structures used for the algorithm
  vertices.push_back({0, 1, 0.0});
  for(size_t k=1; k+1<polyline.size();k++)
    vertices.push_back({k-1, k+1, area(k-1, k, k+1)});
  vertices.push_back({polyline.size()-2, polyline.size(), 0.0});

  // Simplify polyline until N points are left (-2 for start/end points)
  while(by_contribution.size() > N - 2)
  {
    // Select and erase a Point with smallest current contribution to polyline shape
    std::make_heap(by_contribution.begin(), by_contribution.end(), cmp);
    std::pop_heap(by_contribution.begin(), by_contribution.end(), cmp);
    auto curr = by_contribution.back();
    by_contribution.pop_back();
    std::cout << vertices[curr].contribution << " ";

    // Update previous and next Vertex:
    // - update neighbours neighbours
    // - update neighbours contribution
    auto prev = vertices[curr].prev;
    auto next = vertices[curr].next;
    vertices[prev].next = vertices[curr].next;
    vertices[next].prev = vertices[curr].prev;
    vertices[prev].contribution = area(vertices[prev].prev, prev, vertices[prev].next);
    vertices[next].contribution = area(vertices[next].prev, next, vertices[next].next);
  }
std::cout << std::endl;
  // Reconstruct simplified polyline
  PointVector simplified;
  simplified.reserve(N);
  for(size_t k = 0; k < polyline.size(); k = vertices[k].next)
    simplified.push_back(polyline[k]);
  return simplified;
}
