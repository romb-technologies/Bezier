#include "polycurve.h"

namespace Bezier
{
PolyCurve::PolyCurve(CurvePtr& curve1, CurvePtr& curve2, LinkType type)
{
  // we need to connect two curves, but there are four possible combinations
  // find and pick combination for two nearest endpoints
  // TODO: LinkType isn't implemented yet
  Point s_1, s_2, e_1, e_2;
  std::tie(s_1, e_1) = curve1->getEndPoints();
  std::tie(s_1, e_1) = curve2->getEndPoints();

  double s_s = (s_1 - s_2).norm();
  double s_e = (s_1 - e_2).norm();
  double e_s = (e_1 - s_2).norm();
  double e_e = (e_1 - e_2).norm();

  Link::LinkPoint point;

  if (s_s <= s_e && s_s <= e_s && s_s <= e_e)
  {
    curve1->manipulateControlPoint(0, (s_1 + s_2) / 2);
    curve2->manipulateControlPoint(0, (s_1 + s_2) / 2);
    point = Link::LinkPoint::s_s;
  }
  else if (s_e <= e_s && s_e <= e_e)
  {
    curve1->manipulateControlPoint(0, (s_1 + e_2) / 2);
    curve2->manipulateControlPoint(curve2->getOrder(), (s_1 + e_2) / 2);
    point = Link::LinkPoint::s_e;
  }
  else if (e_s <= e_e)
  {
    curve1->manipulateControlPoint(curve1->getOrder(), (e_1 + s_2) / 2);
    curve2->manipulateControlPoint(0, (e_1 + s_2) / 2);
    point = Link::LinkPoint::e_s;
  }
  else
  {
    curve1->manipulateControlPoint(curve1->getOrder(), (e_1 + e_2) / 2);
    curve2->manipulateControlPoint(curve2->getOrder(), (e_1 + e_2) / 2);
    point = Link::LinkPoint::e_e;
  }

  chain_.push_back(Link(std::make_pair(curve1, curve2), type, point));
}

PolyCurve::PolyCurve(std::vector<CurvePtr>& curve_list, const std::vector<PolyCurve::LinkType>& type_list)
    : PolyCurve(curve_list[0], curve_list[1], type_list.empty() ? C0 : type_list[0])
{
  // first pair is created by constructor delegation
  // for adding new curve, we only have to check free end point of last curve in chain
  for (uint k = 1; k < curve_list.size() - 1; k++)
  {
    Point s_1, s_2, e_1, e_2;
    std::tie(s_1, e_1) = curve_list[k]->getEndPoints();
    std::tie(s_1, e_1) = curve_list[k + 1]->getEndPoints();
    Link::LinkPoint point;
    if (chain_.at(k - 1).link_point == Link::LinkPoint::s_s || chain_.at(k - 1).link_point == Link::LinkPoint::e_s)
    {
      double e_s = (e_1 - s_2).norm();
      double e_e = (e_1 - e_2).norm();
      if (e_s <= e_e)
      {
        curve_list[k]->manipulateControlPoint(curve_list[k]->getOrder(), (e_1 + s_2) / 2);
        curve_list[k + 1]->manipulateControlPoint(0, (e_1 + s_2) / 2);
        point = Link::LinkPoint::e_s;
      }
      else
      {
        curve_list[k]->manipulateControlPoint(curve_list[k]->getOrder(), (e_1 + e_2) / 2);
        curve_list[k + 1]->manipulateControlPoint(curve_list[k + 1]->getOrder(), (e_1 + e_2) / 2);
        point = Link::LinkPoint::e_e;
      }
    }
    else
    {
      double s_s = (s_1 - s_2).norm();
      double s_e = (s_1 - e_2).norm();
      if (s_s <= s_e)
      {
        curve_list[k]->manipulateControlPoint(0, (s_1 + s_2) / 2);
        curve_list[k + 1]->manipulateControlPoint(0, (s_1 + s_2) / 2);
        point = Link::LinkPoint::s_s;
      }
      else
      {
        curve_list[k]->manipulateControlPoint(0, (s_1 + e_2) / 2);
        curve_list[k + 1]->manipulateControlPoint(curve_list[k + 1]->getOrder(), (s_1 + e_2) / 2);
        point = Link::LinkPoint::s_e;
      }
    }

    chain_.push_back(Link(std::make_pair(curve_list[k], curve_list[k + 1]), type_list.empty() ? C0 : type_list[k], point));
  }
}

void PolyCurve::addCurve(CurvePtr& curve)
{
  // TODO: check which end of polycurve is nearer
}

void PolyCurve::removeFirst()
{
  chain_.pop_front();
}

void PolyCurve::removeLast()
{
  chain_.pop_back();
}

PolyCurve PolyCurve::getSubPolyCurve(uint idx_l, uint idx_r)
{
  std::vector<CurvePtr> poly;
  std::vector<LinkType> type;
  for(; idx_l < idx_r; idx_l++)
  {
    poly.push_back(getCurvePtr(idx_l));
    type.push_back(chain_.at(idx_l).link_type);
  }
  poly.push_back(getCurvePtr(idx_r));
  return PolyCurve(poly, type);
}

uint PolyCurve::getSize() const { return static_cast<uint>(chain_.size() + 1); }

int PolyCurve::getCurveIdx(const CurvePtr& curve) const
{
  int idx = -1;
  for (uint k = 0; k < getSize(); k++)
    if (curve == getCurvePtr(k))
    {
      idx = static_cast<int>(k);
      break;
    }
  return idx;
}

CurvePtr PolyCurve::getCurvePtr(uint idx) const
{
  if (idx < chain_.size())
    return chain_.at(idx).curve_pair.first;
  if (idx == chain_.size() && !chain_.empty())
    return chain_.back().curve_pair.first;
  return CurvePtr();
}

PointVector PolyCurve::getPolyline(double smoothness, double precision) const
{
  PointVector polyline;
  for(uint k = 0; k < getSize(); k++)
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
  auto ep_s = getCurvePtr(0)->getEndPoints();
  auto ep_e = getCurvePtr(getSize()-1)->getEndPoints();
  if (chain_.front().link_point == Link::LinkPoint::s_e || chain_.front().link_point == Link::LinkPoint::s_s)
    end_points.first = ep_s.first;
  else
    end_points.first = ep_s.second;

  if (chain_.back().link_point == Link::LinkPoint::s_s || chain_.back().link_point == Link::LinkPoint::e_s)
    end_points.second = ep_e.second;
  else
    end_points.second = ep_e.first;
  return std::pair<Point, Point>();
}

Point PolyCurve::valueAt(double t) const
{
  uint idx = static_cast<uint>(t);
  return getCurvePtr(idx)->valueAt(t - idx);
}

double PolyCurve::curvatureAt(double t) const
{
  uint idx = static_cast<uint>(t);
  return getCurvePtr(idx)->curvatureAt(t - idx);
}

Vec2 PolyCurve::tangentAt(double t) const
{
  uint idx = static_cast<uint>(t);
  return getCurvePtr(idx)->tangentAt(t - idx);
}

Vec2 PolyCurve::normalAt(double t) const
{
  uint idx = static_cast<uint>(t);
  return getCurvePtr(idx)->normalAt(t - idx);
}

BBox PolyCurve::getBBox(bool use_roots) const
{
  BBox bbox;
  for (uint k = 0; k < getSize(); k++)
    bbox.extend(getCurvePtr(k)->getBBox(use_roots));
  return bbox;
}

std::vector<Point> PolyCurve::getPointsOfIntersection(const Curve& curve, bool stop_at_first, double epsilon) const
{
  std::vector<Point> points;
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

std::vector<Point> PolyCurve::getPointsOfIntersection(const PolyCurve& poly_curve, bool stop_at_first,
                                                      double epsilon) const
{
  std::vector<Point> points;
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
  double min_t = getCurvePtr(0)->projectPoint(point, step);
  double min_dist = (point - getCurvePtr(0)->valueAt(min_t)).norm();

  for (uint k = 1; k < getSize(); k++)
  {
    double t = getCurvePtr(k)->projectPoint(point, step);
    double dist = (point - getCurvePtr(k)->valueAt(t)).norm();
    if (dist < min_dist)
    {
      min_dist = dist;
      min_t = k + t;
    }
  }
  return min_t;
}
}
