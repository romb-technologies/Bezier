#ifndef POLYLINE_H
#define POLYLINE_H

#include "declarations.h"

namespace Bezier
{
class PolyLine : public PointVector
{
public:
  PolyLine();

  double getLength() const;

  double getPercentage(const Point& point) const;
};
}

#endif // POLYLINE_H
