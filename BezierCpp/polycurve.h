#ifndef POLYCURVE_H
#define POLYCURVE_H

#include <bezier.h>
#include <deque>

namespace Bezier
{
class PolyCurve
{
public:
  enum LinkType // functionality not yet implemented
  {
    C0,
    C1,
    G1 // etc
  };

private:
  struct Link
  {
    enum LinkPoint
    {
      s_s,
      s_e,
      e_s,
      e_e
    } link_point;
    CurvePtr c[2];
    LinkType link_type;
    Link(CurvePtr& curve1, CurvePtr& curve2, LinkType type, LinkPoint point)
    {
      c[0] = curve1;
      c[1] = curve2;
      link_type = type;
      link_point = point;
    }
  };

  std::deque<Link> chain_;

public:
  PolyCurve(CurvePtr& curve1, CurvePtr& curve2, LinkType type = C0);
  PolyCurve(std::vector<CurvePtr>& curve_list, const std::vector<LinkType>& type_list = std::vector<LinkType>());

  void addCurve(CurvePtr& curve);
  void removeFirst();
  void removeLast();
  uint getSize() const;
  int getCurveIdx(const CurvePtr& curve) const;
  CurvePtr getCurvePtr(uint idx) const;

  Point valueAt(double t) const;
  double curvatureAt(double t) const;
  Vec2 tangentAt(double t) const;
  Vec2 normalAt(double t) const;

  BBox getBBox(bool use_roots = true) const;

  std::vector<Point> getPointsOfIntersection(const Curve& curve, bool stop_at_first = false,
                                             double epsilon = 0.001) const;
  std::vector<Point> getPointsOfIntersection(const PolyCurve& poly_curve, bool stop_at_first = false,
                                             double epsilon = 0.001) const;
  double projectPoint(const Point& point, double step = 0.01) const;
};
}
#endif // POLYCURVE_H
