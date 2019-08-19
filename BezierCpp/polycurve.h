/*
 * Copyright 2019 Mirko Kokot
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef POLYCURVE_H
#define POLYCURVE_H

#include "declarations.h"

namespace Bezier
{

/*!
 * \brief A structure for declaring continuity type
 * \var type Type of order (C - parametric, G - geometric)
 * \var order Order of continuity
 */
struct Continuity
{
  char type;
  uint order;

  Continuity(char type = 'C', uint order = 0) : type(type), order(order){}
};

/*!
 * \brief A Bezier polycurve class
 *
 * A class for linking multiple Bezier curves with at least
 * C0 continuity. It allows subcurve and continuity manipulation.
 * Both parametric and geometric continuity are supported.
 *
 * \warning Range of parameter 't' depends on number of subcurves.
 * To access n-th subcurve, t has to be in range [n-1, n>
 */
class PolyCurve
{
private:
  /// Structure for holding underlying Bezier curves
  std::deque<CurvePtr> curves_;
  /// Structure containing information about continuity of all curve neighbours
  std::vector<Continuity> continuity_;

  /*!
   * \brief Constructor for easier creation of sub-polycurve
   * \param curve_list A list of continuus sub-curves
   */
  PolyCurve(const std::deque<CurvePtr>& curve_list);

  /*!
   * \brief Preparation checks for applying continuity
   * \param idx Index of subcurve
   */
  void prepareForContinuity(uint idx);

  /*!
   * \brief Calculate and apply continuity from one curve to its neighbour
   * \param from Origin curve (curve on which calculations are based on)
   * \param to Recipient curve (curve to which calculations are applied)
   * \param con Continuity to apply
   */
  void applyContinuity(CurvePtr from, CurvePtr to, Continuity con);

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /*!
   * \brief Create the Bezier polycurve with only one subcurve
   * \param curve A single curve
   */
  PolyCurve(CurvePtr& curve);

  /*!
   * \brief Create the Bezier polycurve from vector of curves
   * \param curve_list A list of curves
   */
  PolyCurve(std::vector<CurvePtr>& curve_list);

  /*!
   * \brief Insert new curve into polycurve
   * \param idx Index where to insert new curve
   * \param curve A curve to insert
   */
  void insertAt(uint idx, CurvePtr& curve);

  /*!
   * \brief Insert new curve at the beginning of polycurve
   * \param curve A curve to insert
   */
  void insertFront(CurvePtr& curve);

  /*!
   * \brief Insert new curve at the end of polycurve
   * \param curve A curve to insert
   */
  void insertBack(CurvePtr& curve);

  /*!
   * \brief Remove a subcurve from polycurve
   * \param idx Index of subcurve to remove
   */
  void removeAt(uint idx);

  /*!
   * \brief Remove a subcurve from the beginning of polycurve
   */
  void removeFirst();

  /*!
   * \brief Remove a subcurve from the end of polycurve
   */
  void removeBack();

  /*!
   * \brief Set continuity between two subcurves
   * \param idx Index of subcurve joint
   * \param c Type of continuity
   */
  void setContinuity(uint idx, Continuity c = Continuity());

  /*!
   * \brief Get sub-polycurve
   * \param idx_l Index of first subcurve (start)
   * \param idx_r Index of last subcurve (end)
   */
  PolyCurve getSubPolyCurve(uint idx_l, uint idx_r);

  /*!
   * \brief Get number of subcurves
   * \return Number of subcurves
   */
  uint getSize() const;

  /*!
   * \brief Get index of subcurve
   * \param curve A subcurve pointer
   * \return An index of subcurve
   */
  uint getCurveIdx(const CurvePtr& curve) const;

  /*!
   * \brief Get pointer of subcurve
   * \param idx Subcurve index
   * \return A shared pointer
   */
  CurvePtr getCurvePtr(uint idx) const;

  /*!
   * \brief Get list of all subcurves
   * \return A vector of pointers
   */
  std::vector<CurvePtr> getCurveList() const;

  /*!
   * \brief Get a polyline representation of polycurve as a vector of points on curve
   * \param smoothness Smoothness factor > 1 (more resulting points when closer to 1)
   * \param precision Minimal distance between two subsequent points
   * \return A vector of polyline vertices
   */
  PointVector getPolyline(double smoothness = 1.0001, double precision = 1.0) const;

  /*!
   * \brief Get first and last control points
   * \return A pair of end points
   */
  std::pair<Point, Point> getEndPoints() const;

  /*!
   * \brief Get the control points of all subcurves
   * \return A vector of control points
   */
  PointVector getControlPoints() const;

  /*!
   * \brief Set the new coordinates to a control point
   * \param index Index of chosen control point
   * \param point New control point
   */
  void manipulateControlPoint(uint idx, const Point& point);

  /*!
   * \brief Get the point on polycurve for a given t
   * \param t Curve parameter
   * \return Point on a polycurve for a given t
   */
  Point valueAt(double t) const;

  /*!
   * \brief Get curvature of polycurve for a given t
   * \param t Curve parameter
   * \return Curvature of a polycurve for a given t
   */
  double curvatureAt(double t) const;

  /*!
   * \brief Get the tangent of polycurve for a given t
   * \param t Curve parameter
   * \param normalize If the resulting tangent should be normalized
   * \return Tangent of a polycurve for a given t
   */
  Vec2 tangentAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the normal of polycurve for a given t
   * \param t Curve parameter
   * \param normalize If the resulting normal should be normalized
   * \return Normal of a polycurve for given t
   */
  Vec2 normalAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the bounding box of polycurve
   * \param use_roots If algorithm should use roots
   * \return Bounding box (if use_roots is false, returns the bounding box of control points)
   */
  BBox getBBox(bool use_roots = true) const;

  /*!
   * \brief Get the points of intersection with another curve
   * \param curve Curve to intersect with
   * \param stop_at_first If first point of intersection is enough
   * \param epsilon Precision of resulting intersection
   * \return A vector af points of intersection between curves
   *
   * \warning subcurve self-intersection not yet implemented
   */
  std::vector<Point> getPointsOfIntersection(const Curve& curve, bool stop_at_first = false,
                                             double epsilon = 0.001) const;

  /*!
   * \brief Get the points of intersection with another polycurve
   * \param poly_curve Polycurve to intersect with
   * \param stop_at_first If first point of intersection is enough
   * \param epsilon Precision of resulting intersection
   * \return A vector af points of intersection between curves
   *
   * \warning subcurve self-intersection not yet implemented
   */
  std::vector<Point> getPointsOfIntersection(const PolyCurve& poly_curve, bool stop_at_first = false,
                                             double epsilon = 0.001) const;

  /*!
   * \brief Get the parameter t where polycurve is closest to given point
   * \param point Point to project on polycurve
   * \param step Size of step in coarse search
   * \param epsilon Precision of resulting projection
   * \return Parameter t
   */
  double projectPoint(const Point& point, double step = 0.01) const;
};
}
#endif // POLYCURVE_H
