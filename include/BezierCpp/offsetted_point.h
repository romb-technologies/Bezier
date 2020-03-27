#ifndef OFFSETTED_POINT_H
#define OFFSETTED_POINT_H

#include "declarations.h"

namespace Bezier {

template<typename Curve_PolyCurve>
Point offsettedPoint(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
Point offsettedPointDerivation_1(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
Point offsettedPointDerivation_2(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
Point offsettedPointDerivation_3(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
double offsettedPointKappa(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
double offsettedPointKappaDerived(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
double offsettedPointPathLenght(const Curve_PolyCurve& curve, Bezier::Point offset);

template<typename Curve_PolyCurve>
double offsettedPointPathLenght(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
double offsettedPointPathLenght(const Curve_PolyCurve& curve, double t1, double t2, Bezier::Point offset);

}

#endif // OFFSETTED_POINT_H
