#ifndef OFFSETTED_POINT_H
#define OFFSETTED_POINT_H

#include "declarations.h"

namespace Bezier {

template<typename Curve_PolyCurve>
Point getOffsettedPoint(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
Point getOffsettedPointDerivation_1(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
Point getOffsettedPointDerivation_2(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
Point getOffsettedPointDerivation_3(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
double getOffsettedPointKappa(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
double getOffsettedPointKappaDerived(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
double getOffsettedPointPathLenght(const Curve_PolyCurve& curve, Bezier::Point offset);

template<typename Curve_PolyCurve>
double getOffsettedPointPathLenght(const Curve_PolyCurve& curve, double t, Bezier::Point offset);

template<typename Curve_PolyCurve>
double getOffsettedPointPathLenght(const Curve_PolyCurve& curve, double t1, double t2, Bezier::Point offset);

}

#endif // OFFSETTED_POINT_H
