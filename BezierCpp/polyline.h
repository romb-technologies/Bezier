#ifndef POLYLINE_H
#define POLYLINE_H

#include "declarations.h"

namespace Bezier {
class PolyLine : public PointVector
{
public:
  PolyLine();

  double getLength();

  double getPercentage(Point point);
};
}

#endif // POLYLINE_H
