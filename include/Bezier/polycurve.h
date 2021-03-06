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

#include <deque>

#include "declarations.h"

namespace Bezier
{

/*!
 * \brief A Bezier polycurve class
 *
 * A class for linking multiple Bezier curves with at least
 * C0 continuity. It allows subcurve manipulation.
 *
 * \warning Range of parameter 't' depends on number of subcurves.
 * To access n-th subcurve, t has to be in range [n-1, n>
 */
class PolyCurve
{
public:
  /*!
   * \brief Create the empty Bezier polycurve
   */
  PolyCurve() = default;

  /*!
   * \brief Create the Bezier polycurve with only one subcurve
   * \param curve A single curve
   */
  PolyCurve(std::shared_ptr<Curve>& curve);

  /*!
   * \brief Create the Bezier polycurve from vector of curves
   * \param curve_list A list of curves
   */
  PolyCurve(std::vector<std::shared_ptr<Curve>>& curve_list);

  /*!
   * \brief Create a copy of Bezier polycurve
   * \param polycurve A Bezier polycurve to copy
   */
  PolyCurve(const PolyCurve& poly_curve);

  /*!
   * \brief Insert new curve into polycurve
   * \param idx Index where to insert new curve
   * \param curve A curve to insert
   */
  void insertAt(uint idx, std::shared_ptr<Curve>& curve);

  /*!
   * \brief Insert new curve at the beginning of polycurve
   * \param curve A curve to insert
   */
  void insertFront(std::shared_ptr<Curve>& curve);

  /*!
   * \brief Insert new curve at the end of polycurve
   * \param curve A curve to insert
   */
  void insertBack(std::shared_ptr<Curve>& curve);

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
   * \brief Get sub-polycurve
   * \param idx_l Index of first subcurve (start)
   * \param idx_r Index of last subcurve (end)
   */
  PolyCurve subPolyCurve(uint idx_l, uint idx_r) const;

  /*!
   * \brief Get number of subcurves
   * \return Number of subcurves
   */
  uint size() const;

  /*!
   * \brief Resolve polycurve parameter to subcurve index
   * \param t A polycurve parameter
   * \return An index of of subcurve where parameter t is
   */
  uint curveIdx(double t) const;

  /*!
   * \brief Get pointer of a subcurve
   * \param idx Subcurve index
   * \return A shared pointer
   */
  std::shared_ptr<Curve> curvePtr(uint idx) const;

  /*!
   * \brief Get list of all subcurves
   * \return A vector of pointers
   */
  std::vector<std::shared_ptr<Curve>> curveList() const;

  /*!
   * \brief Get a polyline representation of polycurve as a vector of points on curve
   * \param smoothness Smoothness factor > 1 (more resulting points when closer to 1)
   * \param precision Minimal distance between two subsequent points
   * \return A vector of polyline vertices
   */
  PointVector polyline(double smoothness = 1.0001, double precision = 1.0) const;

  /*!
   * \brief Compute exaxt arc length with Legendre-Gauss quadrature
   * \return Arc length
   * \warning Precision depends on value of LEGENDRE_GAUSS_N at compile time
   */
  double length() const;

  /*!
   * \brief Compute exact arc length with Legendre-Gauss quadrature
   * \param t A Polyurve parameter to which length is computed
   * \return Arc length from start to parameter t
   * \warning Precision depends on value of LEGENDRE_GAUSS_N at compile time
   */
  double length(double t) const;

  /*!
   * \brief Compute exact arc length with Legendre-Gauss quadrature
   * \param t1 A Polyurve parameter from which length is computed
   * \param t2 A Polyurve parameter to which length is computed
   * \return Arc length between paramaters t1 and t2
   * \warning Precision depends on value of LEGENDRE_GAUSS_N at compile time
   */
  double length(double t1, double t2) const;

  /*!
   * \brief Compute parameter t which is S distance from given t
   * \param t Curve parameter
   * \param s Distance to iterate
   * \param epsilon Precision of resulting t
   * \param max_iter Maximum number of iterations for Newton-Rhapson
   * \return New parameter t
   */
  double iterateByLength(double t, double s, double epsilon = 0.001, std::size_t max_iter = 15) const;

  /*!
   * \brief Get first and last control points
   * \return A pair of end points
   */
  std::pair<Point, Point> endPoints() const;

  /*!
   * \brief Get the control points of all subcurves
   * \return A vector of control points
   */
  PointVector controlPoints() const;

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
   * \param t A Polyurve parameter
   * \return Curvature of a polycurve for a given t
   */
  double curvatureAt(double t) const;

  /*!
   * \brief Get curvature derivative of curve for a given t
   * \param t Curve parameter
   * \return Curvature derivative of a curve for a given t
   */
  double curvatureDerivativeAt(double t) const;

  /*!
   * \brief Get the tangent of polycurve for a given t
   * \param t A Polyurve parameter
   * \param normalize If the resulting tangent should be normalized
   * \return Tangent of a polycurve for a given t
   */
  Vector tangentAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the normal of polycurve for a given t
   * \param t A Polyurve parameter
   * \param normalize If the resulting normal should be normalized
   * \return Normal of a polycurve for given t
   */
  Vector normalAt(double t, bool normalize = true) const;

  /*!
   * \brief Get value of a derivative for a given t
   * \param t Curve parameter
   * \return Derivative curve
   */
  Point derivativeAt(double t) const;

  /*!
   * \brief Get value of an nth derivative for a given t
   * \param n Desired number of derivative
   * \param t Curve parameter
   * \return Derivative curve
   */
  Point derivativeAt(uint n, double t) const;

  /*!
   * \brief Get the bounding box of polycurve
   * \param use_roots If algorithm should use roots
   * \return Bounding box (if use_roots is false, returns the bounding box of control points)
   */
  BoundingBox boundingBox(bool use_roots = true) const;

  /*!
   * \brief Get the points of intersection with another curve or polycurve
   * \param curve Curve to intersect with
   * \param stop_at_first If first point of intersection is enough
   * \param epsilon Precision of resulting intersection
   * \return A vector af points of intersection between curves
   */
  template <typename Curve_PolyCurve>
  std::vector<Point> pointsOfIntersection(const Curve_PolyCurve& curve, bool stop_at_first = false,
                                             double epsilon = 0.001) const;

  /*!
   * \brief Get the parameter t where polycurve is closest to given point
   * \param point Point to project on polycurve
   * \param step Size of step in coarse search
   * \param epsilon Precision of resulting projection
   * \return Parameter t
   */
  double projectPoint(const Point& point, double step = 0.01, double epsilon = 0.001) const;

private:
  /// Structure for holding underlying Bezier curves
  std::deque<std::shared_ptr<Curve>> curves_;

  /*!
   * \brief Constructor for easier creation of sub-polycurve
   * \param curve_list A list of continuus sub-curves
   */
  PolyCurve(std::deque<std::shared_ptr<Curve>>  curve_list);
};

} // namespace Bezier
#endif // POLYCURVE_H
