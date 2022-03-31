#include "Bezier/polycurve.h"

#include <execution>
#include <functional>
#include <mutex>
#include <numeric>
#include <utility>

enum Policy
{
  seq,
  par
};

template <class F> auto maybe_parallel(Policy p, F f)
{
  switch (p)
  {
  case seq:
    return f(std::execution::seq);
  case par:
    return f(std::execution::par_unseq);
  }
}

using namespace Bezier;

PolyCurve::PolyCurve(std::deque<Curve> curves) : curves_(std::move(curves)) {}

void PolyCurve::insertAt(uint idx, Curve curve) { curves_.emplace(curves_.begin() + idx, curve); }

void PolyCurve::insertFront(Curve curve) { curves_.emplace_front(curve); }

void PolyCurve::insertBack(Curve curve) { curves_.emplace_back(curve); }

void PolyCurve::removeAt(uint idx) { curves_.erase(curves_.begin() + idx); }

void PolyCurve::removeFirst() { curves_.pop_front(); }

void PolyCurve::removeBack() { curves_.pop_back(); }

uint PolyCurve::size() const { return static_cast<uint>(curves_.size()); }

uint PolyCurve::curveIdx(double t) const
{
  uint idx = static_cast<uint>(t);
  return idx - (idx == size());
}

Curve& PolyCurve::curve(uint idx) { return curves_[idx]; }

const Curve& PolyCurve::curve(uint idx) const { return curves_[idx]; }

std::deque<Curve>& PolyCurve::curves() { return curves_; }

const std::deque<Curve>& PolyCurve::curves() const { return curves_; }

PointVector PolyCurve::polyline(double flatness) const
{
  PointVector polyline;
  for (uint idx = 0; idx < size(); idx++)
  {
    auto new_poly = curves_[idx].polyline(flatness);
    polyline.reserve(polyline.size() + new_poly.size() - (idx ? 1 : 0));
    polyline.insert(polyline.end(), new_poly.begin() + (idx ? 1 : 0), new_poly.end());
  }
  return polyline;
}

double PolyCurve::length() const { return length(0, size()); }

double PolyCurve::length(double t) const { return length(0, t); }

double PolyCurve::length(double t1, double t2) const
{
  uint idx1 = curveIdx(t1);
  uint idx2 = curveIdx(t2);

  if (idx1 == idx2)
    return curves_[idx1].length(t1 - idx1, t2 - idx2);
  if (idx1 + 1 == idx2)
    return curves_[idx1].length(t1 - idx1, 1.0) + curves_[idx2].length(0.0, t2 - idx2);

  return maybe_parallel(size() > 300 ? par : seq, [&](auto& pol) {
    return std::reduce(
        pol, begin(curves_) + idx1 + 1, begin(curves_) + idx2,
        curves_[idx1].length(t1 - idx1, 1.0) + curves_[idx2].length(0.0, t2 - idx2), [](const auto& a, const auto& b) {
          if constexpr (std::is_same<decltype(a), const Curve&>() && std::is_same<decltype(b), const Curve&>())
            return a.length() + b.length();
          if constexpr (std::is_same<decltype(a), const double&>() && std::is_same<decltype(b), const double&>())
            return a + b;
          if constexpr (std::is_same<decltype(a), const double&>() && std::is_same<decltype(b), const Curve&>())
            return a + b.length();
          if constexpr (std::is_same<decltype(a), const Curve&>() && std::is_same<decltype(b), const double&>())
            return a.length() + b;
        });
  });
}

double PolyCurve::iterateByLength(double t, double s, double epsilon) const
{
  double s_t = length(t);

  if (s_t + s < 0)
    return 0;
  if (s_t + s > length())
    return size();

  uint idx = curveIdx(t);
  t -= idx;
  s_t -= length(idx);

  bool first_iteration = true;
  while (curves_[idx].length() - s_t < s)
  {
    s -= curve(idx++).length() - s_t;
    if (first_iteration)
    {
      first_iteration = false;
      s_t = 0;
      t = 0;
    }
  }

  return idx + curves_[idx].iterateByLength(t, s, epsilon);
}

std::pair<Point, Point> PolyCurve::endPoints() const
{
  return std::make_pair(curves_[0].endPoints().first, curves_[size() - 1].endPoints().second);
}

PointVector PolyCurve::controlPoints() const
{
  PointVector cp;
  for (const auto& curve : curves_)
  {
    auto cp_c = curve.controlPoints();
    cp.reserve(cp.size() + cp_c.size());
    cp.insert(cp.end(), cp_c.begin(), cp_c.end());
  }
  return cp;
}

void PolyCurve::moveControlPoint(uint idx, const Point& point)
{
  for (auto& curve : curves_)
    if (idx <= curve.order())
    {
      curve.moveControlPoint(idx, point);
      break;
    }
    else
      --idx -= curve.order();
}

Point PolyCurve::valueAt(double t) const
{
  uint idx = curveIdx(t);
  return curves_[idx].valueAt(t - idx);
}

PointVector PolyCurve::valueAt(const std::vector<double>& t_vector) const
{
  PointVector points(t_vector.size());
  std::transform(t_vector.begin(), t_vector.end(), points.begin(), [this](double t) { return valueAt(t); });
  return points;
}

double PolyCurve::curvatureAt(double t) const
{
  uint idx = curveIdx(t);
  return curves_[idx].curvatureAt(t - idx);
}

double PolyCurve::curvatureDerivativeAt(double t) const
{
  uint idx = curveIdx(t);
  return curves_[idx].curvatureDerivativeAt(t - idx);
}

Vector PolyCurve::tangentAt(double t, bool normalize) const
{
  uint idx = curveIdx(t);
  return curves_[idx].tangentAt(t - idx, normalize);
}

Vector PolyCurve::normalAt(double t, bool normalize) const
{
  uint idx = curveIdx(t);
  return curves_[idx].normalAt(t - idx, normalize);
}

Point PolyCurve::derivativeAt(double t) const
{
  uint idx = curveIdx(t);
  return curves_[idx].derivativeAt(t - idx);
}

Point PolyCurve::derivativeAt(uint n, double t) const
{
  uint idx = curveIdx(t);
  return curves_[idx].derivativeAt(n, t - idx);
}

BoundingBox PolyCurve::boundingBox() const
{
  BoundingBox bbox;
  maybe_parallel(size() > 10 ? par : seq, [&](auto& pol) {
    std::for_each(pol, curves_.begin(), curves_.end(),
                  [&bbox](const Curve& curve) { bbox.extend(curve.boundingBox()); });
  });
  return bbox;
}

template <> PointVector PolyCurve::intersections<Curve>(const Curve& curve, double epsilon) const
{
  PointVector points;
  std::mutex m;
  std::for_each(std::execution::par_unseq, curves_.begin(), curves_.end(), [&](const Curve& curve_) {
    auto new_points = curve_.intersections(curve, epsilon);
    std::scoped_lock l(m);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
  });
  return points;
}

template <> PointVector PolyCurve::intersections<PolyCurve>(const PolyCurve& poly_curve, double epsilon) const
{
  PointVector points;
  std::mutex m;
  std::for_each(std::execution::par_unseq, curves_.begin(), curves_.end(), [&](const Curve& curve) {
    auto new_points = poly_curve.intersections(curve, epsilon);
    std::scoped_lock l(m);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
  });
  return points;
}

double PolyCurve::projectPoint(const Point& point) const
{
  auto getPair = [&point](const Curve& curve) -> std::pair<double, double> {
    double t = curve.projectPoint(point);
    return {t, (point - curve.valueAt(t)).norm()};
  };

  return maybe_parallel(size() > 20 ? par : seq, [&](auto& pol) {
    return std::reduce(pol, curves_.begin() + 1, curves_.end(), getPair(curves_[0]),
                       [&getPair](const auto& a, const auto& b) {
                         if constexpr (std::is_same<decltype(a), const Curve&>() &&
                                       std::is_same<decltype(b), const Curve&>())
                         {
                           auto a_ = getPair(a);
                           auto b_ = getPair(b);
                           return a_.second < b_.second ? a_ : b_;
                         }
                         if constexpr (std::is_same<decltype(a), const std::pair<double, double>&>() &&
                                       std::is_same<decltype(b), const std::pair<double, double>&>())
                         {
                           return a.second < b.second ? a : b;
                         }
                         if constexpr (std::is_same<decltype(a), const std::pair<double, double>&>() &&
                                       std::is_same<decltype(b), const Curve&>())
                         {
                           auto b_ = getPair(b);
                           return a.second < b_.second ? a : b_;
                         }
                         if constexpr (std::is_same<decltype(a), const Curve&>() &&
                                       std::is_same<decltype(b), const std::pair<double, double>&>())
                         {
                           auto a_ = getPair(a);
                           return a_.second < b.second ? a_ : b;
                         }
                       })
        .first;
  });
}

std::vector<double> PolyCurve::projectPoint(const PointVector& point_vector) const
{
  std::vector<double> t_vector(point_vector.size());
  std::transform(point_vector.begin(), point_vector.end(), t_vector.begin(),
                 [this](const Point& point) { return projectPoint(point); });
  return t_vector;
}

double PolyCurve::distance(const Point& point) const { return (point - valueAt(projectPoint(point))).norm(); }

std::vector<double> PolyCurve::distance(const PointVector& point_vector) const
{
  std::vector<double> dist_vector(point_vector.size());
  std::transform(point_vector.begin(), point_vector.end(), dist_vector.begin(),
                 [this](const Point& point) { return distance(point); });
  return dist_vector;
}
