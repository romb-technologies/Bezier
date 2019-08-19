#ifndef DECLARATIONS_H
#define DECLARATIONS_H

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <memory>
#include <deque>

/*!
 * Nominal namespace containing class pre-definitions and typedefs
 */
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
class Curve;

/*!
 * \brief A polyline class
 *
 * A class for operations on polyline, mainly used for
 * polyline representation of curves. Underlying structure
 * is a vector of points.
 */
class PolyLine;

/*!
 * \brief A Bezier polycurve class
 *
 * A class for linking multiple Bezier curves with at least
 * C0 continuity. It allows subcurve and continuity manipulation.
 * Both parametric and geometric continuity are supported.
 */
class PolyCurve;

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
/*!
 * \brief Shared pointer of Polycurve
 */
typedef std::shared_ptr<PolyCurve> PolyCurvePtr;
/*!
 * \brief Shared pointer of const Polycurve
 */
typedef std::shared_ptr<const PolyCurve> ConstPolyCurvePtr;
}
#endif // DECLARATIONS_H
