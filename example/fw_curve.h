#ifndef FW_CURVE_H
#define FW_CURVE_H

#include "bezier.h"

constexpr double L = 15;
constexpr double d = 0;

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

inline Bezier::Point getFrontWheelPosition(const Bezier::Curve& curve, double t)
{
  Bezier::Point p_fw;

  auto P = curve.valueAt(t);
  double x = P.x();
  double y = P.y();

  auto dP = curve.getDerivative()->valueAt(t);
  double dx = dP.x();
  double dy = dP.y();

  auto ddP = curve.getDerivative()->getDerivative()->valueAt(t);
  double ddx = ddP.x();
  double ddy = ddP.y();

  auto helper_1 = dx * ddx * 2.0 + dy * ddy * 2.0;
  auto helper_2a = ddx / dP.norm() - helper_1 * dx / pow(dP.norm(), 3) / 2.0;
  auto helper_2 = pow(fabs(helper_2a), 2.0);
  auto helper_3a = ddy / dP.norm() - helper_1 * dy / pow(dP.norm(), 3) / 2.0;
  auto helper_3 = pow(fabs(helper_3a), 2.0);

  p_fw.x() = x + L * 1.0 / dP.norm() * dx + d * helper_2a / sqrt(helper_2 + helper_3);

  p_fw.y() = y + L * 1.0 / dP.norm() * dy + d * helper_3a / sqrt(helper_2 + helper_3);

  return p_fw;
}

#endif // FW_CURVE_H
