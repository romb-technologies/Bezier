#include "polyline.h"

namespace Bezier {
PolyLine::PolyLine()
{

}

double PolyLine::getLength()
{
  double length = 0.0;
  for(uint k = 1; k < size(); k++)
    length += (at(k-1) - at(k)).norm();
  return length;
}

double PolyLine::getPercentage(Point point)
{
  double lenght_to_point = 0.0;
  double lenght = 0.0;
  bool found = false;
  for(uint k = 1; k < size(); k++)
  {
    double d1 = (at(k-1) - point).norm();
    double d2 = (at(k) - point).norm();
    double l = (at(k-1) - at(k)).norm();
    lenght += l;
    if(!found)
    {
      if(d1 < l && d2 <l)
      {
        lenght_to_point = d1;
        found = true;
      }
      else {
        lenght_to_point += l;
      }
    }
  }

  return lenght_to_point / lenght;
}
}

