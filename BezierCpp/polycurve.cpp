#include "polycurve.h"

namespace Bezier {
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

  chain_.push_back(Link(curve1, curve2, type, point));
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
    if (chain_[k - 1].link_point == Link::LinkPoint::s_s || chain_[k - 1].link_point == Link::LinkPoint::e_s)
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

    chain_.push_back(Link(curve_list[k], curve_list[k + 1], type_list.empty() ? C0 : type_list[k], point));
  }
}

void PolyCurve::addCurve(CurvePtr& curve)
{
  // TODO: check which end of polycurve is nearer
}

void PolyCurve::removeFirst()
{
  if (!chain_.empty())
    chain_.pop_front();
}

void PolyCurve::removeLast()
{
  if (!chain_.empty())
    chain_.pop_back();
}

uint PolyCurve::getSize() const { return static_cast<uint>(chain_.size() + 1); }

int PolyCurve::getCurveIdx(const CurvePtr& curve) const
{
  int idx = -1;
  for (uint k = 0; k < chain_.size(); k++)
    if (curve == chain_[k].c[0])
    {
      idx = static_cast<int>(k);
      break;
    }
  if (idx < 0)
    if (curve == chain_.back().c[1])
      idx = static_cast<int>(chain_.size());
  return idx;
}

CurvePtr PolyCurve::getCurvePtr(uint idx) const
{
  if (idx < chain_.size())
    return chain_[idx].c[0];
  if (idx == chain_.size() && !chain_.empty())
    return chain_.back().c[1];
  return CurvePtr();
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
  if (!chain_.empty())
  {
    bbox = chain_[0].c[0]->getBBox(use_roots);
    for (uint k = 0; k < chain_.size(); k++)
      bbox.extend(chain_[k].c[1]->getBBox(use_roots));
  }
  return bbox;
}

std::vector<Point> PolyCurve::getPointsOfIntersection(const Curve& curve, bool stop_at_first, double epsilon) const
{
  std::vector<Point> points;
  if (!chain_.empty())
  {
    auto new_points = chain_[0].c[0]->getPointsOfIntersection(curve, stop_at_first, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), new_points.begin(), new_points.end());
    if (points.empty() || !stop_at_first)
      for (uint k = 0; k < chain_.size(); k++)
      {
        new_points = chain_[k].c[1]->getPointsOfIntersection(curve, stop_at_first, epsilon);
        points.reserve(points.size() + new_points.size());
        points.insert(points.end(), new_points.begin(), new_points.end());
        if (!points.empty() && stop_at_first)
          break;
      }
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
  double min_t = 0;
  if (!chain_.empty())
  {
    double min_t = chain_[0].c[0]->projectPoint(point, step);
    double min_dist = (point - chain_[0].c[0]->valueAt(min_t)).norm();
    for (uint k = 0; k < chain_.size(); k++)
    {
      double t = chain_[k].c[1]->projectPoint(point, step);
      double dist = (point - chain_[k].c[1]->valueAt(t)).norm();
      if (dist < min_dist)
      {
        min_dist = dist;
        min_t = 1 + k + t;
      }
    }
  }
  return min_t;
}
}
