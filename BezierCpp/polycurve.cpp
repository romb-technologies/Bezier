#include "polycurve.h"

namespace Bezier
{
PolyCurve::PolyCurve(const std::deque<PolyCurve::Part>& part_list) : parts_(part_list) {}

PolyCurve::PolyCurve(CurvePtr& curve) { parts_.push_back(Part(curve, false)); }

PolyCurve::PolyCurve(std::vector<CurvePtr>& curve_list) : PolyCurve(curve_list.front())
{
  // first curve is inserted with constructor delegation
  bool first = true;
  for (auto&& curve : curve_list)
    if (first)
      first = false;
    else
      addCurve(curve);
}

void PolyCurve::addCurve(CurvePtr& curve)
{
  // TODO: LinkType isn't implemented yet

  Point s_1, s_2, e_1, e_2;
  std::tie(s_1, e_1) = curve->getEndPoints();
  std::tie(s_2, e_2) = this->getEndPoints();

  double s_s = (s_1 - s_2).norm();
  double s_e = (s_1 - e_2).norm();
  double e_s = (e_1 - s_2).norm();
  double e_e = (e_1 - e_2).norm();

  if (e_s <= s_s && e_s <= s_e && e_s <= e_e)
  {
    // insert in front without reversing
    parts_.push_front(Part(curve, false));
  }
  else if (s_s <= s_e && s_s <= e_e)
  {
    // insert in front with reversing
    parts_.push_front(Part(curve, true));
  }
  else if (s_e <= e_e)
  {
    // insert at back without reversing
    parts_.push_back(Part(curve, false));
  }
  else
  {
    // insert at back with reversing
    parts_.push_back(Part(curve, true));
  }
}

void PolyCurve::removeFirst() { parts_.pop_front(); }

void PolyCurve::removeLast() { parts_.pop_back(); }

PolyCurve PolyCurve::getSubPolyCurve(uint idx_l, uint idx_r)
{
  return PolyCurve(std::deque<Part>(parts_.begin() + idx_l, parts_.begin() + idx_r));
}

uint PolyCurve::getSize() const { return static_cast<uint>(parts_.size()); }

uint PolyCurve::getCurveIdx(const CurvePtr& curve) const
{
  uint idx = 0;
  for (; idx < getSize(); idx++)
    if (curve == parts_.at(idx).curve)
      break;
  return idx;
}

CurvePtr PolyCurve::getCurvePtr(uint idx) const { return parts_.at(idx).curve; }

PointVector PolyCurve::getPolyline(double smoothness, double precision) const
{
  PointVector polyline;
  for (uint k = 0; k < getSize(); k++)
  {
    auto new_poly = getCurvePtr(k)->getPolyline(smoothness, precision);
    polyline.reserve(polyline.size() + new_poly.size());
    polyline.insert(polyline.end(), new_poly.begin(), new_poly.end());
  }
  return polyline;
}

std::pair<Point, Point> PolyCurve::getEndPoints() const
{
  std::pair<Point, Point> end_points;
  if (parts_.front().reversed)
    end_points.first = parts_.front().curve->getEndPoints().second;
  else
    end_points.first = parts_.front().curve->getEndPoints().first;

  if (parts_.back().reversed)
    end_points.second = parts_.front().curve->getEndPoints().first;
  else
    end_points.second = parts_.front().curve->getEndPoints().second;
  return end_points;
}

Point PolyCurve::valueAt(double t) const
{
  uint idx = static_cast<uint>(t);
  if (idx == getSize()) --idx;
  double t_new = parts_.at(idx).reversed ? 1 - (t - idx) : t - idx;
  return getCurvePtr(idx)->valueAt(t_new);
}

double PolyCurve::curvatureAt(double t) const
{
  uint idx = static_cast<uint>(t);
  if (idx == getSize()) --idx;
  double t_new = parts_.at(idx).reversed ? 1 - (t - idx) : t - idx;
  return getCurvePtr(idx)->curvatureAt(t_new);
}

Vec2 PolyCurve::tangentAt(double t, bool normalize) const
{
  uint idx = static_cast<uint>(t);
  if (idx == getSize()) --idx;
  double t_new = parts_.at(idx).reversed ? 1 - (t - idx) : t - idx;
  return (parts_.at(idx).reversed ? -1 : 1) * getCurvePtr(idx)->tangentAt(t_new, normalize);
}

Vec2 PolyCurve::normalAt(double t, bool normalize) const { return tangentAt(t, normalize); }

BBox PolyCurve::getBBox(bool use_roots) const
{
  BBox bbox;
  for (uint k = 0; k < getSize(); k++)
    bbox.extend(getCurvePtr(k)->getBBox(use_roots));
  return bbox;
}

PointVector PolyCurve::getPointsOfIntersection(const Curve& curve, bool stop_at_first, double epsilon) const
{
  PointVector points;
  for (uint k = 0; k < getSize(); k++)
  {
    auto new_points = getCurvePtr(k)->getPointsOfIntersection(curve, stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (!points.empty() && stop_at_first)
      break;
  }
  return points;
}

PointVector PolyCurve::getPointsOfIntersection(const PolyCurve& poly_curve, bool stop_at_first, double epsilon) const
{
  PointVector points;
  for (uint k = 0; k < poly_curve.getSize(); k++)
  {
    auto new_points = getPointsOfIntersection(*poly_curve.getCurvePtr(k), stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (!points.empty() && stop_at_first)
      break;
  }
  return points;
}

double PolyCurve::projectPoint(const Point& point, double step) const
{
  double min_t = parts_.front().curve->projectPoint(point, step);
  double min_dist = (point - parts_.front().curve->valueAt(min_t)).norm();
  min_t = parts_.front().reversed ? 1 - min_t : min_t;

  for (uint k = 1; k < getSize(); k++)
  {
    double t = getCurvePtr(k)->projectPoint(point, step);
    double dist = (point - getCurvePtr(k)->valueAt(t)).norm();
    if (dist < min_dist)
    {
      min_dist = dist;
      min_t = k + (parts_.at(k).reversed ? 1 - t : t);
    }
  }
  return min_t;
}
}
