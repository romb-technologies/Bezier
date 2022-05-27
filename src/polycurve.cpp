#include "Bezier/polycurve.h"

#include <numeric>

using namespace Bezier;

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

PointVector PolyCurve::polyline(double flatness) const
{
  PointVector polyline;
  for (unsigned idx = 0; idx < size(); idx++)
  {
    auto new_poly = curves_[idx].polyline(flatness);
    polyline.reserve(polyline.size() + new_poly.size() - (idx ? 1 : 0));
    polyline.insert(polyline.end(), std::make_move_iterator(new_poly.begin() + (idx ? 1 : 0)),
                    std::make_move_iterator(new_poly.end()));
  }
  return polyline;
}

double PolyCurve::length() const { return length(0, size()); }

double PolyCurve::length(double t) const { return length(0, t); }

double PolyCurve::length(double t1, double t2, double epsilon) const
{
  unsigned idx1 = curveIdx(t1);
  unsigned idx2 = curveIdx(t2);

  if (idx1 == idx2)
    return curves_[idx1].length(t1 - idx1, t2 - idx2, epsilon);
  if (idx1 + 1 == idx2)
    return curves_[idx1].length(t1 - idx1, 1.0, epsilon) + curves_[idx2].length(0.0, t2 - idx2, epsilon);

  return std::accumulate(begin(curves_) + idx1 + 1, begin(curves_) + idx2,
                         curves_[idx1].length(t1 - idx1, 1.0, epsilon) + curves_[idx2].length(0.0, t2 - idx2, epsilon),
                         [&epsilon](double sum, const Curve& curve) { return sum + curve.length(0, 1, epsilon); });
}

double PolyCurve::iterateByLength(double t, double s, double epsilon) const
{
  double s_t = length(t);

  if (s_t + s < 0)
    return 0;
  if (s_t + s > length())
    return size();

  unsigned idx = curveIdx(t);
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
    cp.insert(cp.end(), std::make_move_iterator(cp_c.begin()), std::make_move_iterator(cp_c.end()));
  }
  return cp;
}

void PolyCurve::setControlPoint(unsigned idx, const Point& point)
{
  for (auto& curve : curves_)
    if (idx <= curve.order())
    {
      curve.setControlPoint(idx, point);
      break;
    }
    else
      --idx -= curve.order();
}

Point PolyCurve::valueAt(double t) const
{
  unsigned idx = curveIdx(t);
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

// namespace is a workaround for a bug on old gcc versions:
// https://stackoverflow.com/questions/25311512/specialization-of-template-in-different-namespace
namespace Bezier
{

template <> PointVector PolyCurve::intersections<Curve>(const Curve& curve, double epsilon) const
{
  PointVector points;
  for (const auto& curve2 : curves_)
  {
    auto new_points = curve2.intersections(curve, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), std::make_move_iterator(new_points.begin()), std::make_move_iterator(new_points.end()));
  }
  return points;
}

template <> PointVector PolyCurve::intersections<PolyCurve>(const PolyCurve& poly_curve, double epsilon) const
{
  PointVector points;
  for (const auto& curve : curves_)
  {
    auto new_points = poly_curve.intersections(curve, epsilon);
    points.reserve(points.size() + new_points.size());
    points.insert(points.end(), std::make_move_iterator(new_points.begin()), std::make_move_iterator(new_points.end()));
  }
  return points;
}

} // namespace Bezier

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
