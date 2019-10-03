#include "polycurve.h"
#include "bezier.h"

inline double binomial(uint n, uint k) { return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1)); }

namespace Bezier
{
PolyCurve::PolyCurve(const std::deque<CurvePtr>& part_list) : curves_(part_list) {}

void PolyCurve::prepareForContinuity(uint idx)
{
  // raise this curve enough for that only one end point can be affected
  auto curve = getCurvePtr(idx);

  while (curve->getOrder() < continuity_[idx].order + 1)
    curve->elevateOrder();

  if (idx > 0)
  {
    // backward
    // no need to check order (it is always checked in forward direction)
    getCurvePtr(idx - 1)->reverse();
    getCurvePtr(idx)->reverse();
    applyContinuity(getCurvePtr(idx), getCurvePtr(idx - 1), continuity_.at(idx - 1));
    getCurvePtr(idx - 1)->reverse();
    getCurvePtr(idx)->reverse();
  }
  if (idx < getSize() - 1)
  {
    // forward
    // raise next curve enough for that only one end point can be affected
    auto curve_f = getCurvePtr(idx + 1);
    while (curve_f->getOrder() < continuity_[idx].order + 1)
      curve_f->elevateOrder();
    applyContinuity(getCurvePtr(idx), getCurvePtr(idx + 1), continuity_.at(idx));
  }
}

void PolyCurve::applyContinuity(CurvePtr from, CurvePtr to, Continuity con)
{
  PointVector N;
  PointVector P = from->getControlPoints();
  PointVector Q = to->getControlPoints();

  uint n = from->getOrder();
  uint m = to->getOrder();

  auto F = [n, m](uint k)
  {
    // n! == tgamma(n +1)
    return std::tgamma(m - k + 1) / std::tgamma(m + 1) *
           std::tgamma(n + 1) / std::tgamma(n - k + 1);
  };

  auto R = [from, to, con](uint k)
  {
    if (k == 0 || con.type == 'C') return 1.0;
    return to->getDerivative()->valueAt(0).norm() / from->getDerivative()->valueAt(1).norm();
//    ConstCurvePtr der_q = from;
//    ConstCurvePtr der_p = to;
//    for(uint d = 0; d < k; d++)
//    {
//      der_q = der_q->getDerivative();
//      der_p = der_p->getDerivative();
//    }
//    return der_p->valueAt(0).norm() / der_q->valueAt(1).norm();
  };

  for (uint x = 0; x <= con.order; x++)
  {
    Point new_point;
    new_point << 0, 0;
    for (uint i = 0; i <= x; i++)
    {
      double inner_sum = 0.0;
      for (uint k = i; k <= x; k++)
        inner_sum += binomial(x-i, k-i)*F(k) * R(k);
      new_point += std::pow(-1, i) * binomial(x, i) * P.at(n - i) * inner_sum;
    }
    N.push_back(new_point);
  }

  for (uint x = 0; x <= con.order; x++)
    to->manipulateControlPoint(x, N.at(x));
}

PolyCurve::PolyCurve(CurvePtr& curve) { curves_.push_back(curve); }

PolyCurve::PolyCurve(std::vector<CurvePtr>& curve_list)
{
  for (auto&& curve : curve_list)
    insertBack(curve);
}

void PolyCurve::insertAt(uint idx, CurvePtr& curve)
{
  // TODO: Continuity isn't implemented yet
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
  continuity_.resize(getSize() - 1);
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

void PolyCurve::setContinuity(uint idx, Continuity c) {
  continuity_[idx] = c;
  prepareForContinuity(idx);
}

PolyCurve PolyCurve::getSubPolyCurve(uint idx_l, uint idx_r)
{
  return PolyCurve(std::deque<CurvePtr>(curves_.begin() + idx_l, curves_.begin() + idx_r));
}

uint PolyCurve::getSize() const { return static_cast<uint>(curves_.size()); }

uint PolyCurve::getCurveIdx(const CurvePtr& curve) const
{
  uint idx = 0;
  for (; idx < getSize(); idx++)
    if (curve == curves_.at(idx))
      break;
  return idx;
}

CurvePtr PolyCurve::getCurvePtr(uint idx) const { return curves_.at(idx); }

std::vector<CurvePtr> PolyCurve::getCurveList() const
{
  return std::vector<CurvePtr>(curves_.begin(), curves_.end());
}

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
  for(auto& curve : curves_)
    l += curve->getLength();
  return l;
}

double PolyCurve::getLength(double t) const
{
  double l = 0;
  uint idx = static_cast<uint>(t);
  if (idx == getSize()) // for the last point of last curve
    --idx;
  for (uint k = 0; k < idx; k++)
    l += curves_.at(k)->getLength();
  l += curves_.at(idx)->getLength(t - idx);
  return l;
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
      prepareForContinuity(getCurveIdx(curve));
      break;
    }
    else
      --idx -= curve->getOrder();
}

Point PolyCurve::valueAt(double t) const
{
  uint idx = static_cast<uint>(t);
  if (idx == getSize()) // for the last point of last curve
    --idx;
  return getCurvePtr(idx)->valueAt(t - idx);
}

double PolyCurve::curvatureAt(double t) const
{
  uint idx = static_cast<uint>(t);
  if (idx == getSize()) // for the last point of last curve
    --idx;
  return getCurvePtr(idx)->curvatureAt(t - idx);
}

Vec2 PolyCurve::tangentAt(double t, bool normalize) const
{
  uint idx = static_cast<uint>(t);
  if (idx == getSize()) // for the last point of last curve
    --idx;
  return getCurvePtr(idx)->tangentAt(t - idx, normalize);
}

Vec2 PolyCurve::normalAt(double t, bool normalize) const
{
  uint idx = static_cast<uint>(t);
  if (idx == getSize()) // for the last point of last curve
    --idx;
  return getCurvePtr(idx)->normalAt(t - idx, normalize);
}

BBox PolyCurve::getBBox(bool use_roots) const
{
  BBox bbox;
  for (uint k = 0; k < getSize(); k++)
    bbox.extend(curves_.at(k)->getBBox(use_roots));
  return bbox;
}

PointVector PolyCurve::getPointsOfIntersection(const Curve& curve, bool stop_at_first, double epsilon) const
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

PointVector PolyCurve::getPointsOfIntersection(const PolyCurve& poly_curve, bool stop_at_first, double epsilon) const
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
}
