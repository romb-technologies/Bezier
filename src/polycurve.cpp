#include "BezierCpp/polycurve.h"
#include "BezierCpp/bezier.h"

#include <numeric>

inline double binomial(uint n, uint k) { return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1)); }

namespace Bezier
{
PolyCurve::PolyCurve(const std::deque<CurvePtr>& part_list) : curves_(part_list) {}

PolyCurve::PolyCurve(CurvePtr& curve) { curves_.push_back(curve); }

PolyCurve::PolyCurve(std::vector<CurvePtr>& curve_list)
{
  for (auto&& curve : curve_list)
    insertBack(curve);
}

void PolyCurve::insertAt(uint idx, CurvePtr& curve)
{
  Point s_1, s_2, e_1, e_2;
  std::tie(s_1, e_1) = curve->getEndPoints();
  if (idx > 0) // check with curve before
  {
    std::tie(s_2, e_2) = curves_.at(idx - 1)->getEndPoints();

    double s_e = (s_1 - e_2).norm();
    double e_e = (e_1 - e_2).norm();

    if (e_e < s_e) // we need to reverse the curve
    {
      curve->reverse();
      std::tie(s_1, e_1) = curve->getEndPoints();
    }

    curve->manipulateControlPoint(0, (s_1 + e_2) / 2);
    curves_.at(idx - 1)->manipulateControlPoint(curves_.at(idx - 1)->getOrder(), (s_1 + e_2) / 2);
  }
  if (idx + 1 < getSize()) // check with curve after
  {
    std::tie(s_2, e_2) = curves_.at(idx)->getEndPoints();

    double e_s = (e_1 - s_2).norm();
    double s_s = (s_1 - s_2).norm();

    if (s_s < e_s) // we ned to reverse the curve
    {
      curve->reverse();
      std::tie(s_1, e_1) = curve->getEndPoints();
    }

    curve->manipulateControlPoint(curve->getOrder(), (e_1 + s_2) / 2);
    curves_.at(idx + 1)->manipulateControlPoint(0, (e_1 + s_2) / 2);
  }

  curves_.insert(curves_.begin() + idx, curve);
}

void PolyCurve::insertFront(CurvePtr& curve) { insertAt(0, curve); }

void PolyCurve::insertBack(CurvePtr& curve) { insertAt(getSize(), curve); }

void PolyCurve::removeAt(uint idx)
{
  if (idx == 0)
    removeFirst();
  else if (idx == getSize() - 1)
    removeBack();
  else
  {
    Point s_1, s_2, e_1, e_2;
    std::tie(s_1, e_1) = curves_.at(idx - 1)->getEndPoints();
    std::tie(s_2, e_2) = curves_.at(idx + 1)->getEndPoints();
    curves_.at(idx - 1)->manipulateControlPoint(curves_.at(idx - 1)->getOrder(), (e_1 + s_2) / 2);
    curves_.at(idx + 1)->manipulateControlPoint(0, (e_1 + s_2) / 2);
    curves_.erase(curves_.begin() + idx);
  }
}

void PolyCurve::removeFirst() { curves_.pop_front(); }

void PolyCurve::removeBack() { curves_.pop_back(); }

PolyCurve PolyCurve::getSubPolyCurve(uint idx_l, uint idx_r) const
{
  return PolyCurve(std::deque<CurvePtr>(curves_.begin() + idx_l, curves_.begin() + idx_r));
}

uint PolyCurve::getSize() const { return static_cast<uint>(curves_.size()); }

std::pair<uint, double> PolyCurve::getCurveIdx(double t) const
{
  uint idx = static_cast<uint>(t);
  if (idx == getSize())
    --idx; // for the last point of last curve
  return {idx, t - idx};
}

CurvePtr PolyCurve::getCurvePtr(uint idx) const { return curves_.at(idx); }

std::vector<CurvePtr> PolyCurve::getCurveList() const { return std::vector<CurvePtr>(curves_.begin(), curves_.end()); }

PointVector PolyCurve::getPolyline(double smoothness, double precision) const
{
  PointVector polyline;
  for (uint k = 0; k < getSize(); k++)
  {
    auto new_poly = getCurvePtr(k)->getPolyline(smoothness, precision);
    polyline.reserve(polyline.size() + new_poly.size() - (k > 0 ? 1 : 0));
    polyline.insert(polyline.end(), new_poly.begin() + (k > 0 ? 1 : 0), new_poly.end());
  }
  return polyline;
}

double PolyCurve::getLength() const
{
  double l = 0;
  for (auto& curve : curves_)
    l += curve->getLength();
  return l;
}

double PolyCurve::getLength(double t) const
{
  uint idx;
  std::tie(idx, t) = getCurveIdx(t);

  return std::accumulate(begin(curves_), begin(curves_) + idx, curves_.at(idx)->getLength(t),
                         [](double sum, Bezier::ConstCurvePtr curve) { return sum + curve->getLength(); });
}

double PolyCurve::getLength(double t1, double t2) const { return getLength(t2) - getLength(t1); }

double PolyCurve::iterateByLength(double t, double s, double epsilon, std::size_t max_iter) const
{
  double s_t = getLength(t);
  if (s_t + s < 0 || s_t + s > getLength())
    throw std::out_of_range{"Resulting parameter t not in [0, n] range."};

  uint idx;
  std::tie(idx, t) = getCurveIdx(t);
  s_t -= getLength(idx);

  while (getCurvePtr(idx)->getLength() - s_t < s)
  {
    s -= getCurvePtr(idx++)->getLength() - s_t;
    s_t = 0;
    t = 0;
  }

  return idx + getCurvePtr(idx)->iterateByLength(t, s, epsilon, max_iter);
}

std::pair<Point, Point> PolyCurve::getEndPoints() const
{
  return std::make_pair(curves_.front()->getEndPoints().first, curves_.back()->getEndPoints().second);
}

PointVector PolyCurve::getControlPoints() const
{
  PointVector cp;
  for (auto&& curve : curves_)
  {
    auto cp_c = curve->getControlPoints();
    cp.reserve(cp.size() + cp_c.size());
    cp.insert(cp.end(), cp_c.begin(), cp_c.end());
  }
  return cp;
}

void PolyCurve::manipulateControlPoint(uint idx, const Point& point)
{
  for (auto&& curve : curves_)
    if (idx <= curve->getOrder())
    {
      curve->manipulateControlPoint(idx, point);
      break;
    }
    else
      --idx -= curve->getOrder();
}

Point PolyCurve::valueAt(double t) const
{
  uint idx;
  std::tie(idx, t) = getCurveIdx(t);
  return getCurvePtr(idx)->valueAt(t);
}

double PolyCurve::curvatureAt(double t) const
{
  uint idx;
  std::tie(idx, t) = getCurveIdx(t);
  return getCurvePtr(idx)->curvatureAt(t);
}

Vec2 PolyCurve::tangentAt(double t, bool normalize) const
{
  uint idx;
  std::tie(idx, t) = getCurveIdx(t);
  return getCurvePtr(idx)->tangentAt(t, normalize);
}

Vec2 PolyCurve::normalAt(double t, bool normalize) const
{
  uint idx;
  std::tie(idx, t) = getCurveIdx(t);
  return getCurvePtr(idx)->normalAt(t, normalize);
}

BBox PolyCurve::getBBox(bool use_roots) const
{
  BBox bbox;
  for (uint k = 0; k < getSize(); k++)
    bbox.extend(curves_.at(k)->getBBox(use_roots));
  return bbox;
}

template <>
PointVector PolyCurve::getPointsOfIntersection<Curve>(const Curve& curve, bool stop_at_first, double epsilon) const
{
  PointVector points;
  for (uint k = 0; k < getSize(); k++)
  {
    auto new_points = curves_.at(k)->getPointsOfIntersection(curve, stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (!points.empty() && stop_at_first)
      break;
  }
  return points;
}

template <>
PointVector PolyCurve::getPointsOfIntersection<PolyCurve>(const PolyCurve& poly_curve, bool stop_at_first,
                                                          double epsilon) const
{
  PointVector points;
  for (uint k = 0; k < getSize(); k++)
  {
    auto new_points = poly_curve.getPointsOfIntersection(*curves_.at(k), stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (!points.empty() && stop_at_first)
      break;
  }
  return points;
}

double PolyCurve::projectPoint(const Point& point, double step, double epsilon) const
{
  double min_t = curves_.front()->projectPoint(point, step, epsilon);
  double min_dist = (point - curves_.front()->valueAt(min_t)).norm();

  for (uint k = 1; k < getSize(); k++)
  {
    double t = curves_.at(k)->projectPoint(point, step, epsilon);
    double dist = (point - curves_.at(k)->valueAt(t)).norm();
    if (dist < min_dist)
    {
      min_dist = dist;
      min_t = k + t;
    }
  }
  return min_t;
}
} // namespace Bezier
