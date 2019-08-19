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
    auto p1 = at(k) - at(k - 1);
    auto p2 = point - at(k - 1);
    double k1 = p1.y() / p1.x();
    double k2 = p2.y() / p2.x();
    // check if point is between current pair of control points:
    // point is closer to first cp than cp
    // AND they are near collinear (slope method)
    if (p2.norm() <= p1.norm() && std::fabs(k1 - k2) < 0.1)
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

Point PolyLine::getPoint(double percent) const
{
  double length = percent * getLength();
  Point point;

  for (uint k = 1; k < size(); k++)
  {
    double p_dist = (at(k) - at(k - 1)).norm();
    if (p_dist < length)
    {
      length -= p_dist;
    }
    else
    {
      point = at(k - 1) + (at(k) - at(k - 1)) * length / p_dist;
      break;
    }
  }

  return point;
}
}
