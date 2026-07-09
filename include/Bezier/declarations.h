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

#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <Eigen/Dense>

/*!
 * Nominal namespace containing forward declarations and typedefs
 */
namespace Bezier
{
/*!
 * \brief A Bezier curve class
 */
class Curve;

/*!
 * \brief A Bezier polycurve class
 */
class PolyCurve;

/*!
 * \brief Point in xy plane
 */
using Point = Eigen::Vector2d;

/*!
 * \brief A vector of Points
 */
using PointVector = std::vector<Point>;

/*!
 * \brief A vector of curve parameters
 */
using ParamVector = std::vector<double>;

/*!
 * \brief A Vector in xy plane
 */
using Vector = Eigen::Vector2d;

/*!
 * \brief Bounding box class
 */
using BoundingBox = Eigen::AlignedBox2d;
} // namespace Bezier
#endif // DECLARATIONS_H
