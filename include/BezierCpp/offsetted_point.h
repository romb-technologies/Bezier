#ifndef OFFSETTED_POINT_H
#define OFFSETTED_POINT_H

#include "declarations.h"

namespace Bezier {

template<class T_curve>
Point getOffsettedPoint(const T_curve& curve, double t, Bezier::Point offset);

template<class T_curve>
Point getOffsettedPointDerivation_1(const T_curve& curve, double t, Bezier::Point offset);

template<class T_curve>
Point getOffsettedPointDerivation_2(const T_curve& curve, double t, Bezier::Point offset);

template<class T_curve>
Point getOffsettedPointDerivation_3(const T_curve& curve, double t, Bezier::Point offset);

template<class T_curve>
double getOffsettedPointKappa(const T_curve& curve, double t, Bezier::Point offset);

template<class T_curve>
double getOffsettedPointKappaDerived(const T_curve& curve, double t, Bezier::Point offset);

}

#endif // OFFSETTED_POINT_H
