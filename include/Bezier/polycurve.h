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

#include "Bezier/bezier.h"
#include "Bezier/declarations.h"

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
  void insertAt(unsigned idx, Curve curve);

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
  void removeAt(unsigned idx);

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
  unsigned size() const;

  /*!
   * \brief Resolve polycurve parameter to subcurve index
   * \param t A polycurve parameter
   * \return An index of of subcurve where parameter t is
   */
  unsigned curveIdx(double t) const;

  ///@{
  /*!
   * \brief Get pointer of a subcurve
   * \param idx Subcurve index
   * \return A reference to curve
   */
  Curve& curve(unsigned idx);
  const Curve& curve(unsigned idx) const;
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
   * \brief Get a polyline representation of the polycurve as a vector of points on polycurve
   * \return A vector of polyline vertices
   * \note Default flatness parameter is calculated as 0.1% for smallest bounding box diagonal of all subcurves
   */
  std::vector<Point> polyline() const;

  /*!
   * \brief Get a polyline representation of the polycurve as a vector of points on polycurve
   * \param flatness Error tolerance of approximation
   * \return A vector of polyline vertices
   */
  std::vector<Point> polyline(double flatness) const;

  /*!
   * \brief Compute exact arc length using Chebyshev polynomials
   * \return Arc length
   */
  double length() const;

  /*!
   * \brief Compute exact arc length using Chebyshev polynomials
   * \param t A Polycurve parameter to which length is computed
   * \return Arc length from start to parameter t
   */
  double length(double t) const;

  /*!
   * \brief Compute exact arc length using Chebyshev polynomials
   * \param t1 A Polycurve parameter from which length is computed
   * \param t2 A Polycurve parameter to which length is computed
   * \return Arc length between paramaters t1 and t2
   */
  double length(double t1, double t2) const;

  /*!
   * \brief Compute parameter t which is dS distance from given t
   * \param t Curve parameter
   * \param dS Distance to move
   * \return New parameter t
   */
  double step(double t, double dS) const;

  /*!
   * \brief Get first and last control points
   * \return A pair of end points
   */
  std::pair<Point, Point> endPoints() const;

  /*!
   * \brief Get the control points of all subcurves
   * \return A vector of control points
   */
  std::vector<Point> controlPoints() const;

  /*!
   * \brief Set the new coordinates to a control point
   * \param index Index of chosen control point
   * \param point New control point
   */
  void setControlPoint(unsigned idx, const Point& point);

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
  std::vector<Point> valueAt(const std::vector<double>& t_vector) const;

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
  Vector derivativeAt(unsigned n, double t) const;

  /*!
   * \brief Get the bounding box of the polycurve
   * \return Bounding box
   */
  BoundingBox boundingBox() const;

  /*!
   * \brief Get the points of intersection with another curve or polycurve
   * \param curve Curve to intersect with
   * \return A vector af points of intersection between curves
   */
  template <typename Curve_PolyCurve> std::vector<Point> intersections(const Curve_PolyCurve& curve) const;

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
  std::vector<double> projectPoint(const std::vector<Point>& point_vector) const;

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
  std::vector<double> distance(const std::vector<Point>& point_vector) const;

protected:
  /// Structure for holding underlying Bezier curves
  std::deque<Curve> curves_;
};

} // namespace Bezier
#endif // POLYCURVE_H
