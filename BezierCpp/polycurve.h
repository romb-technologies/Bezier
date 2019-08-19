#ifndef POLYCURVE_H
#define POLYCURVE_H

#include "declarations.h"

namespace Bezier
{

struct Continuity
{
  char type;
  uint order;

  Continuity(char type = 'C', uint order = 0) : type(type), order(order){}
};

class PolyCurve
{
private:

  std::deque<CurvePtr> curves_;
  std::vector<Continuity> continuity_;
  // TODO: Continuity isn't implemented yet

  PolyCurve(const std::deque<CurvePtr>& curve_list);

  void applyContinuity(uint idx);
  void calculateContinuity(CurvePtr from, CurvePtr to, Continuity con);

public:
  PolyCurve(CurvePtr& curve1);
  PolyCurve(std::vector<CurvePtr>& curve_list);

  void insertAt(uint idx, CurvePtr& curve);
  void insertFront(CurvePtr& curve);
  void insertBack(CurvePtr& curve);

  void removeAt(uint idx);
  void removeFirst();
  void removeBack();

  void setContinuity(uint idx, Continuity c = Continuity());

  PolyCurve getSubPolyCurve(uint idx_l, uint idx_r);
  uint getSize() const;
  uint getCurveIdx(const CurvePtr& curve) const;
  CurvePtr getCurvePtr(uint idx) const;

  PointVector getPolyline(double smoothness = 1.0001, double precision = 1.0) const;
  std::pair<Point, Point> getEndPoints() const;
  /*!
   * \brief Get the control points
   * \return A vector of control points
   */
  PointVector getControlPoints() const;

  /*!
   * \brief Set the new coordinates to a control point
   * \param index Index of chosen control point
   * \param point New control point
   */
  void manipulateControlPoint(uint idx, const Point& point);

  Point valueAt(double t) const;
  double curvatureAt(double t) const;
  Vec2 tangentAt(double t, bool normalize = true) const;
  Vec2 normalAt(double t, bool normalize = true) const;

  BBox getBBox(bool use_roots = true) const;

  std::vector<Point> getPointsOfIntersection(const Curve& curve, bool stop_at_first = false,
                                             double epsilon = 0.001) const;
  std::vector<Point> getPointsOfIntersection(const PolyCurve& poly_curve, bool stop_at_first = false,
                                             double epsilon = 0.001) const;
  double projectPoint(const Point& point, double step = 0.01) const;
};
}
#endif // POLYCURVE_H