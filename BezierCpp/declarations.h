#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <memory>

#ifdef BEZIER_POLYCURVE
#include <deque>
#endif

namespace Bezier
{
class Curve;

/*!
 * \brief Shared pointer of Curve
 */
typedef std::shared_ptr<Curve> CurvePtr;
/*!
 * \brief Shared pointer of const Curve
 */
typedef std::shared_ptr<const Curve> ConstCurvePtr;
/*!
 * \brief Point in xy plane
 */
typedef Eigen::Vector2d Point;
/*!
 * \brief A Vector in xy plane
 */
typedef Eigen::Vector2d Vec2;
/*!
 * \brief A vector of Points
 */
typedef std::vector<Point> PointVector;
/*!
 * \brief Bounding box
 */
typedef Eigen::AlignedBox2d BBox;

#ifdef BEZIER_POLYLINE
class PolyLine;
#else
typedef PointVector PolyLine;
#endif

#ifdef BEZIER_POLYCURVE
class PolyCurve;
typedef std::shared_ptr<PolyCurve> PolyCurvePtr;
typedef std::shared_ptr<const PolyCurve> ConstPolyCurvePtr;
#endif
}

#endif // DECLARATIONS_H
