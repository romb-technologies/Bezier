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

#ifndef POLYLINE_H
#define POLYLINE_H

#include "declarations.h"

namespace Bezier
{
/*!
 * \brief A polyline class
 *
 * A class for operations on polyline, mainly used for
 * polyline representation of curves. Underlying structure
 * is a vector of points.
 */
class PolyLine : public PointVector
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /*!
   * \brief Create the polyline
   */
  PolyLine();
  /*!
   * \brief Create the polyline
   * \param vector A vector of points that define the polyline
   */
  PolyLine(PointVector vector) : PointVector(vector) {}

  /*!
   * \brief Get the length of polycurve
   */
  double getLength() const;

  /*!
   * \brief Get the percentage of length where point is on a polyline
   * \param point A point for which we calculate percentage
   */
  double getPercentage(const Point& point) const;

  /*!
   * \brief Get the point on polyline corresponding to length percentage
   * \param percent Percentage of length for which we get point on polyline
   */
  Point getPoint(double percent) const;
};
}

#endif // POLYLINE_H
