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
