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

#include "Bezier/declarations.h"
#include "Bezier/bezier.h"

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
  PolyCurve() = default;
  ~PolyCurve() = default;

  /*!
   * \brief Create the Bezier polycurve from deque of curves
   * \param curve_list A list of curves
   */
  PolyCurve(std::deque<Curve> curves);

  PolyCurve(const PolyCurve&) = default;
  PolyCurve(PolyCurve&&) = default;
  PolyCurve& operator=(const PolyCurve&) = default;
  PolyCurve& operator=(PolyCurve&&) = default;

  /*!
   * \brief Insert new curve into polycurve
   * \param idx Index where to insert new curve
   * \param curve A curve to insert
   */
  void insertAt(uint idx, Curve curve);

  /*!
   * \brief Insert new curve at the beginning of the polycurve
   * \param curve A curve to insert
   */
  void insertFront(Curve curve);

  /*!
   * \brief Insert new curve at the end of the polycurve
   * \param curve A curve to insert
   */
  void insertBack(Curve curve);

  /*!
   * \brief Remove a subcurve from the polycurve
   * \param idx Index of subcurve to remove
   */
  void removeAt(uint idx);

  /*!
   * \brief Remove a subcurve from the beginning of the polycurve
   */
  void removeFirst();

  /*!
   * \brief Remove a subcurve from the end of the polycurve
   */
  void removeBack();

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

  ///@{
  /*!
   * \brief Get pointer of a subcurve
   * \param idx Subcurve index
   * \return A reference to curve
   */
  Curve& curve(uint idx);
  const Curve& curve(uint idx) const;
  ///@}

  ///@{
  /*!
   * \brief Get list of all subcurves
   * \return A vector of curve references
   */
  std::deque<Curve>& curves();
  const std::deque<Curve>& curves() const;
  ///@}

  /*!
   * \brief Get a polyline representation of the polycurve as a vector of points on curve
   * \param flatness Error tolerance of approximation
   * \return A vector of polyline vertices
   */
  PointVector polyline(double flatness = 0.5) const;

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
   * \return New parameter t
   */
  double iterateByLength(double t, double s, double epsilon = 0.001) const;

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
  void setControlPoint(uint idx, const Point& point);

  /*!
   * \brief Get the point on polycurve for a given t
   * \param t Curve parameter
   * \return Point on a polycurve for a given t
   */
  Point valueAt(double t) const;

  /*!
   * \brief Get the point vector on polycurve for given parameters
   * \param t_vector Curve parameters
   * \return Vector of points on a polycurve for given parameters
   */
  PointVector valueAt(const std::vector<double>& t_vector) const;

  /*!
   * \brief Get curvature of the polycurve for a given t
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
   * \brief Get the tangent of the polycurve for a given t
   * \param t A Polyurve parameter
   * \param normalize If the resulting tangent should be normalized
   * \return Tangent of a polycurve for a given t
   */
  Vector tangentAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the normal of the polycurve for a given t
   * \param t A Polyurve parameter
   * \param normalize If the resulting normal should be normalized
   * \return Normal of a polycurve for given t
   */
  Vector normalAt(double t, bool normalize = true) const;

  /*!
   * \brief Get value of a derivative for a given t
   * \param t Curve parameter
   * \return Curve derivative at t
   */
  Vector derivativeAt(double t) const;

  /*!
   * \brief Get value of an nth derivative for a given t
   * \param n Desired number of derivative
   * \param t Curve parameter
   * \return nth curve derivative at t
   */
  Vector derivativeAt(uint n, double t) const;

  /*!
   * \brief Get the bounding box of the polycurve
   * \return Bounding box
   */
  BoundingBox boundingBox() const;

  /*!
   * \brief Get the points of intersection with another curve or polycurve
   * \param curve Curve to intersect with
   * \param epsilon Precision of resulting intersection
   * \return A vector af points of intersection between curves
   */
  template <typename Curve_PolyCurve>
  PointVector intersections(const Curve_PolyCurve& curve, double epsilon = 0.001) const;

  /*!
   * \brief Get the parameter t where polycurve is closest to given point
   * \param point Point to project on polycurve
   * \return double t
   */
  double projectPoint(const Point& point) const;

  /*!
   * \brief Get the parameter t vector where polycurve is closest to given points
   * \param point_vector Points to project on polycurve
   * \return Vector of parameters t
   */
  std::vector<double> projectPoint(const PointVector& point_vector) const;

  /*!
   * \brief Get distance of the point to the polycurve
   * \param point Point to project on the polycurve
   * \return Distance to the curve
   */
  double distance(const Point& point) const;

  /*!
   * \brief Get the distance vector of points to the polycurve
   * \param point_vector Points to project on the polycurve
   * \return Vector of distances
   */
  std::vector<double> distance(const PointVector& point_vector) const;

protected:
  /// Structure for holding underlying Bezier curves
  std::deque<Curve> curves_;
};

} // namespace Bezier
#endif // POLYCURVE_H
