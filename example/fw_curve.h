#ifndef FW_CURVE_H
#define FW_CURVE_H

#include "bezier.h"

Bezier::Point getFrontWheelPosition(const Bezier::Curve& curve, double t);

Bezier::Point getFrontWheelCurveDerivation_1(const Bezier::Curve& curve, double t);

Bezier::Point getFrontWheelCurveDerivation_2(const Bezier::Curve& curve, double t);

Bezier::Point getFrontWheelCurveDerivation_3(const Bezier::Curve& curve, double t);

double getFrontWheelCurveKappa(const Bezier::Curve& curve, double t);

double getFrontWheelCurveKappaDerived(const Bezier::Curve& curve, double t);

#endif // FW_CURVE_H
