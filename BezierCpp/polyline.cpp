#include "polyline.h"

namespace Bezier
{
PolyLine::PolyLine() {}

double PolyLine::getLength() const
{
  double length = 0.0;
  for (uint k = 1; k < size(); k++)
    length += (at(k - 1) - at(k)).norm();
  return length;
}

double PolyLine::getPercentage(const Point& point) const
{
  double lenght_to_point = 0.0;
  for (uint k = 1; k < size(); k++)
  {
    auto p1 = at(k) - at(k-1);
    auto p2 = point - at(k-1);
    double x_r = p2.x() / p1.x();
    double y_r = p2.y() / p1.y();
    if ( x_r >= 0 && x_r <= 1 && y_r >= 0 && y_r <= 1)
    {
      lenght_to_point += p2.norm();
      break;
    }
    else
    {
      lenght_to_point += p1.norm();
    }
  }

  return lenght_to_point / getLength();
}
}
