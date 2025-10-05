#include "Bezier/polycurve.h"
#include "Bezier/utils.h"

#include <numeric>

using namespace Bezier;
namespace bu = Bezier::Utils;

PolyCurve::PolyCurve(std::deque<Curve> curves) : curves_(std::move(curves)) {}

void PolyCurve::insertAt(unsigned idx, Curve curve) { curves_.emplace(curves_.begin() + idx, curve); }

void PolyCurve::insertFront(Curve curve) { curves_.emplace_front(curve); }

void PolyCurve::insertBack(Curve curve) { curves_.emplace_back(curve); }

void PolyCurve::removeAt(unsigned idx) { curves_.erase(curves_.begin() + idx); }

void PolyCurve::removeFirst() { curves_.pop_front(); }

void PolyCurve::removeBack() { curves_.pop_back(); }

unsigned PolyCurve::size() const { return curves_.size(); }

unsigned PolyCurve::curveIdx(double t) const
{
  unsigned idx = t;
  return idx - (idx == size());
}

Curve& PolyCurve::curve(unsigned idx) { return curves_[idx]; }

const Curve& PolyCurve::curve(unsigned idx) const { return curves_[idx]; }

std::deque<Curve>& PolyCurve::curves() { return curves_; }

const std::deque<Curve>& PolyCurve::curves() const { return curves_; }

PointVector PolyCurve::polyline() const
{
  double flatness = std::numeric_limits<double>::infinity();
  for (const auto& curve : curves_)
    flatness = std::min(flatness, curve.boundingBox().diagonal().norm() / 1000);
  return polyline(flatness);
}

PointVector PolyCurve::polyline(double flatness) const
{
  if (curves_.empty())
    return {};
  std::vector<Point> polyline = curves_[0].polyline(flatness);
  for (unsigned k{1}; k < curves_.size(); k++)
  {
    polyline.pop_back(); // remove duplicate point
    polyline = bu::concatenate(std::move(polyline), curves_[k].polyline(flatness));
  }
  return polyline;
}

ParamVector PolyCurve::polylineParams() const
{
  double flatness = std::numeric_limits<double>::infinity();
  for (const auto& curve : curves_)
    flatness = std::min(flatness, curve.boundingBox().diagonal().norm() / 1000);
  return polylineParams(flatness);
}

ParamVector PolyCurve::polylineParams(double flatness) const
{
  ParamVector params;
  for (unsigned k{}; k < curves_.size(); k++)
  {
    auto subcurve_params = curves_[k].polylineParams(flatness);
    std::transform(subcurve_params.begin() + (k != 0), subcurve_params.end(), std::back_inserter(params),
                   [k](double t) { return k + t; });
  }
  return params;
}

double PolyCurve::length() const { return length(0.0, size()); }

double PolyCurve::length(double t) const { return length(0.0, t); }

double PolyCurve::length(double t1, double t2) const
{
  int sign{1};
  if (t2 < t1)
  {
    sign = -1;
    std::swap(t1, t2);
  }

  unsigned idx1 = curveIdx(t1);
  unsigned idx2 = curveIdx(t2);

  if (idx1 == idx2)
    return sign * curves_[idx1].length(t1 - idx1, t2 - idx2);
  if (idx1 + 1 == idx2)
    return sign * curves_[idx1].length(t1 - idx1, 1.0) + curves_[idx2].length(t2 - idx2);

  return sign * std::accumulate(curves_.begin() + idx1 + 1, curves_.begin() + idx2,
                                curves_[idx1].length(t1 - idx1, 1.0) + curves_[idx2].length(t2 - idx2),
                                [](double sum, const Curve& curve) { return sum + curve.length(); });
}

double PolyCurve::step(double t, double dS) const
{
  if (std::fabs(dS) < bu::epsilon) // no-op
    return t;

  double s_t = length(t);

  if (dS < -s_t + bu::epsilon) // out-of-scope
    return 0.0;

  if (dS > length() - s_t - bu::epsilon) // out-of-scope
    return size();

  unsigned idx = curveIdx(t);
  t -= idx;

  s_t = dS < 0 ? s_t - length(idx) : length(idx + 1) - s_t;

  while (-s_t > dS + bu::epsilon)
  {
    dS += s_t;
    s_t = curves_[--idx].length();
    t = 1.0;
  }
  while (s_t < dS - bu::epsilon)
  {
    dS -= s_t;
    s_t = curves_[++idx].length();
    t = 0.0;
  }

  return idx + curves_[idx].step(t, dS);
}

std::pair<Point, Point> PolyCurve::endPoints() const
{
  return std::make_pair(curves_[0].endPoints().first, curves_[size() - 1].endPoints().second);
}

PointVector PolyCurve::controlPoints() const
{
  std::vector<Point> cp;
  for (const auto& curve : curves_)
    cp = bu::concatenate(std::move(cp), curve.controlPoints());
  return cp;
}

void PolyCurve::setControlPoint(unsigned idx, const Point& point)
{
  for (auto it = curves_.begin(); it != curves_.end(); idx -= it++->order() + 1)
    if (idx <= it->order())
    {
      it->setControlPoint(idx, point);
      break;
    }
}

Point PolyCurve::valueAt(double t) const
{
  unsigned idx = curveIdx(t);
  return curves_[idx].valueAt(t - idx);
}

PointVector PolyCurve::valueAt(const ParamVector& t_vector) const
{
  std::vector<Point> points(t_vector.size());
  std::transform(t_vector.begin(), t_vector.end(), points.begin(), [this](double t) { return valueAt(t); });
  return points;
}

double PolyCurve::curvatureAt(double t) const
{
  unsigned idx = curveIdx(t);
  return curves_[idx].curvatureAt(t - idx);
}

double PolyCurve::curvatureDerivativeAt(double t) const
{
  unsigned idx = curveIdx(t);
  return curves_[idx].curvatureDerivativeAt(t - idx);
}

Vector PolyCurve::tangentAt(double t, bool normalize) const
{
  unsigned idx = curveIdx(t);
  return curves_[idx].tangentAt(t - idx, normalize);
}

Vector PolyCurve::normalAt(double t, bool normalize) const
{
  unsigned idx = curveIdx(t);
  return curves_[idx].normalAt(t - idx, normalize);
}

Point PolyCurve::derivativeAt(double t) const
{
  unsigned idx = curveIdx(t);
  return curves_[idx].derivativeAt(t - idx);
}

Point PolyCurve::derivativeAt(unsigned n, double t) const
{
  unsigned idx = curveIdx(t);
  return curves_[idx].derivativeAt(n, t - idx);
}

BoundingBox PolyCurve::boundingBox() const
{
  BoundingBox bbox;
  for (const auto& curve : curves_)
    bbox.extend(curve.boundingBox());
  return bbox;
}

template <> PointVector PolyCurve::intersections<Curve>(const Curve& curve) const
{
  std::vector<Point> points;
  for (const auto& curve2 : curves_)
    points = bu::concatenate(std::move(points), curve2.intersections(curve));
  return points;
}

template <> PointVector PolyCurve::intersections<PolyCurve>(const PolyCurve& poly_curve) const
{
  std::vector<Point> points;
  for (const auto& curve : curves_)
    points = bu::concatenate(std::move(points), poly_curve.intersections(curve));
  return points;
}

double PolyCurve::projectPoint(const Point& point) const
{
  double min_t{}, min_dist{std::numeric_limits<double>::infinity()};
  for (unsigned k = 0; k < curves_.size(); k++)
  {
    double t = curves_[k].projectPoint(point);
    double dist = (point - curves_[k].valueAt(t)).norm();
    if (dist < min_dist)
      std::tie(min_t, min_dist) = std::make_tuple(k + t, dist);
  }
  return min_t;
}

ParamVector PolyCurve::projectPoint(const PointVector& point_vector) const
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
