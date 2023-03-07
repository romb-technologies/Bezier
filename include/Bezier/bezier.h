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

#ifndef BEZIER_H
#define BEZIER_H

#include <map>
#include <memory>

#include "declarations.h"

namespace Bezier
{
/*!
 * \brief A Bezier curve class
 *
 * A class for storing and using any-order Bezier curve.
 * It uses private and static caching for storing often accessed data.
 * Private caching is used for data concerning individual curve, while
 * static caching is used for common data (coefficient matrices)
 */
class Curve
{
public:
  ~Curve() = default;

  /*!
   * \brief Create the Bezier curve
   * \param points Nx2 matrix where each row is one of N control points that define the curve
   */
  Curve(Eigen::MatrixX2d points);

  /*!
   * \brief Create the Bezier curve
   * \param points A vector of control points that define the curve
   */
  Curve(const PointVector& points);

  Curve(const Curve& curve);
  Curve(Curve&&) = default;
  Curve& operator=(const Curve&);
  Curve& operator=(Curve&&) = default;

  /*!
   * \brief Get order of the curve (Nth order curve is described with N+1 points);
   * \return Order of curve
   */
  unsigned order() const;

  /*!
   * \brief Get a vector of control points
   * \return A vector of control points
   */
  PointVector controlPoints() const;

  /*!
   * \brief Get the control point at index idx
   * \param idx Index of chosen control point
   * \return A vector of control points
   */
  Point controlPoint(unsigned idx) const;

  /*!
   * \brief Get first and last control points
   * \return A pair of end points
   */
  std::pair<Point, Point> endPoints() const;

  /*!
   * \brief Get a polyline representation of the curve as a vector of points on curve
   * \param flatness Error tolerance of approximation
   * \return A vector of polyline vertices
   */
  PointVector polyline(double flatness = 0.5) const;

  /*!
   * \brief Compute exact arc length using Chebyshev polynomials
   * \return Arc length
   */
  double length() const;

  /*!
   * \brief Compute exact arc length using Chebyshev polynomials
   * \param t Curve parameter to which length is computed
   * \return Arc length from start to parameter t
   */
  double length(double t) const;

  /*!
   * \brief Compute exact arc length using Chebyshev polynomials
   * \param t1 Curve parameter from which length is computed
   * \param t2 Curve parameter to which length is computed
   * \return Arc length between paramaters t1 and t2
   */
  double length(double t1, double t2) const;

  /*!
   * \brief Compute parameter t which is S distance from given t
   * \param t Curve parameter
   * \param s Distance to iterate
   * \return New parameter t
   */
  double iterateByLength(double t, double s) const;

  /*!
   * \brief Reverse order of control points
   */
  void reverse();

  /*!
   * \brief Set the new coordinates to a control point
   * \param idx Index of chosen control point
   * \param point New control point
   */
  void setControlPoint(unsigned idx, const Point& point);

  /*!
   * \brief Manipulate the curve so that it passes through wanted point for given 't'
   * \param t Curve parameter
   * \param point Point where curve should pass for a given t
   *
   * \warning CAN THROW: Only works for quadratic and cubic curves
   * \warning Resets cached data
   */
  void manipulateCurvature(double t, const Point& point);

  /*!
   * \brief Raise the curve order by 1
   *
   * Curve will always retain its shape
   * \warning Resets cached data
   */
  void elevateOrder();

  /*!
   * \brief Lower the curve order by 1
   *
   * If current shape cannot be described by lower order, it will be best aproximation
   * \warning CAN THROW: Cannot be called for curves of 1st order
   * \warning Resets cached data
   */
  void lowerOrder();

  /*!
   * \brief Get the point on curve for a given t
   * \param t Curve parameter
   * \return Point on a curve for a given t
   */
  Point valueAt(double t) const;

  /*!
   * \brief Get the point vector on curve for given parameters
   * \param t_vector Curve parameters
   * \return Matrix of points on a curve for given parameters
   */
  Eigen::MatrixX2d valueAt(const std::vector<double>& t_vector) const;

  /*!
   * \brief Get curvature of the curve for a given t
   * \param t Curve parameter
   * \return Curvature of a curve for a given t
   */
  double curvatureAt(double t) const;

  /*!
   * \brief Get curvature derivative of the curve for a given t
   * \param t Curve parameter
   * \return Curvature derivative of a curve for a given t
   */
  double curvatureDerivativeAt(double t) const;

  /*!
   * \brief Get the tangent of the curve for a given t
   * \param t Curve parameter
   * \param normalize If the resulting tangent should be normalized
   * \return Tangent of a curve for a given t
   */
  Vector tangentAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the normal of the curve for a given t
   * \param t Curve parameter
   * \param normalize If the resulting normal should be normalized
   * \return Normal of a curve for given t
   */
  Vector normalAt(double t, bool normalize = true) const;

  /*!
   * \brief Get the derivative of a curve
   * \return Derivative curve
   */
  const Curve& derivative() const;

  /*!
   * \brief Get the nth derivative of a curve
   * \param n Desired number of derivative
   * \return Derivative curve
   * \warning double n cannot be zero
   */
  const Curve& derivative(unsigned n) const;

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
   * \brief Get roots of the curve on both axes
   * \return A vector of parameters where curve passes through axes
   */
  std::vector<double> roots() const;

  /*!
   * \brief Get all extrema of the curve
   * \return A vector of parameters where extrema are
   */
  std::vector<double> extrema() const;

  /*!
   * \brief Get the bounding box of curve
   * \return Bounding box (if use_roots is false, returns the bounding box of control points)
   */
  BoundingBox boundingBox() const;

  /*!
   * \brief Split the curve into two subcurves
   * \param t Curve parameter at which to split the curve
   * \return Pair of two subcurves
   */
  std::pair<Curve, Curve> splitCurve(double t = 0.5) const;

  /*!
   * \brief Get the points of intersection with another curve
   * \param curve Curve to intersect with
   * \return A vector af points of intersection between curves
   */
  PointVector intersections(const Curve& curve) const;

  /*!
   * \brief Get the parameter t where curve is closest to given point
   * \param point Point to project on curve
   * \return double t
   */
  double projectPoint(const Point& point) const;

  /*!
   * \brief Get distance of the point to the curve
   * \param point Point to project on curve
   * \return Distance to the curve
   */
  double distance(const Point& point) const;

  /*!
   * \brief applyContinuity Apply geometric continuity based on the another curve.
   * \param locked_curve Curve on which calculation are based.
   * \param beta_coeffs Beta-constraints used to calculate continuity. Size defines continuity order.
   */
  void applyContinuity(const Curve& source_curve, const std::vector<double>& beta_coeffs);

protected:
  /*!
   * \brief N x 2 matrix where each row corresponds to control Point
   * \warning Any changes made to control_points_ require a call to resetCache() funtion!
   */
  Eigen::MatrixX2d control_points_;

  /// Reset all privately cached data
  inline void resetCache();

private:
  /// Number of control points (order + 1)
  unsigned N_{};

  /*!
   * \brief Coefficients for matrix operations
   */
  using Coeffs = Eigen::MatrixXd;

  /*!
   * \brief Map of different coefficient matrices, depending on the order of the curve
   */
  using CoeffsMap = std::map<unsigned, Coeffs>;

  // private caching
  mutable std::unique_ptr<const Curve> cached_derivative_;    /*! If generated, stores derivative for later use */
  mutable std::unique_ptr<std::vector<double>> cached_roots_; /*! If generated, stores roots for later use */
  mutable std::unique_ptr<BoundingBox> cached_bounding_box_;  /*! If generated, stores bounding box for later use */
  mutable std::unique_ptr<PointVector> cached_polyline_;      /*! If generated, stores polyline for later use */
  mutable double cached_polyline_flatness_{};                 /*! Flatness of cached polyline */
  mutable std::unique_ptr<Eigen::VectorXd>
      cached_projection_polynomial_part_; /*! Constant part of point projection polynomial */
  mutable Eigen::MatrixXd
      cached_projection_polynomial_derivative_; /*! Polynomial representation of the curve derivative */
  mutable std::unique_ptr<Eigen::VectorXd> cached_chebyshev_coeffs_; /*!  If generated, stores chebyshev coefficients
                                                                        for calculating the length of the curve */

  // static caching
  static CoeffsMap bernstein_coeffs_;       /*! Map of Bernstein coefficients */
  static CoeffsMap splitting_coeffs_left_;  /*! Map of coefficients to get subcurve for t = [0, 0.5] */
  static CoeffsMap splitting_coeffs_right_; /*! Map of coefficients to get subcurve for t = [0.5, 1] */
  static CoeffsMap elevate_order_coeffs_;   /*! Map of coefficients for elevating the order of the curve */
  static CoeffsMap lower_order_coeffs_;     /*! Map of coefficients for lowering the order of the curve */

  /// Static getter function for Bernstein coefficients
  static Coeffs bernsteinCoeffs(unsigned n);
  /// Static getter function for coefficients to get a subcurve [0, z];
  static Coeffs splittingCoeffsLeft(unsigned n, double t = 0.5);
  /// Static getter function for coefficients to get a subcurve [z, 1];
  static Coeffs splittingCoeffsRight(unsigned n, double t = 0.5);
  /// Static getter function for coefficients to elevate order of curve
  static Coeffs elevateOrderCoeffs(unsigned n);
  /// Static getter function for coefficients to lower order of curve
  static Coeffs lowerOrderCoeffs(unsigned n);
};

} // namespace Bezier

#endif // BEZIER_H
