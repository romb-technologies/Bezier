#ifndef FW_CURVE_H
#define FW_CURVE_H

#include "bezier.h"

constexpr double L = 15;
constexpr double d = 0;

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
template <typename T> int dirac(T val) { return val == T(0) ? INFINITY : T(0); }

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

  auto helper_1 = dx * ddx * 2 + dy * ddy * 2;
  auto helper_2 = ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2;
  auto helper_3 = ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2;

  p_fw.x() = x + L * dx / dP.norm() + d * helper_2 / sqrt(helper_2 * helper_2 + helper_3 * helper_3);
  p_fw.y() = y + L * dy / dP.norm() + d * helper_3 / sqrt(helper_2 * helper_2 + helper_3 * helper_3);

  return p_fw;
}

inline Bezier::Point getFrontWheelTangent(const Bezier::Curve& curve, double t)
{
  Bezier::Point t_fw;

  auto P = curve.valueAt(t);
  double x = P.x();
  double y = P.y();

  auto dP = curve.getDerivative()->valueAt(t);
  double dx = dP.x();
  double dy = dP.y();

  auto ddP = curve.getDerivative()->getDerivative()->valueAt(t);
  double ddx = ddP.x();
  double ddy = ddP.y();

  auto dddP = curve.getDerivative()->getDerivative()->getDerivative()->valueAt(t);
  double dddx = dddP.x();
  double dddy = dddP.y();

  auto helper_1 = dx * ddx * 2 + dy * ddy * 2;
  auto helper_2 = dx * dddx * 2 + dy * dddy * 2;

  auto helper_3 = ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2;
  auto helper_4 = ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2;

  auto helper_5 =
      dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx;
  auto helper_6 =
      dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy;
  auto helper_7 = dx / pow(dP.norm(), 3) *
                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + fabs(dx) * dirac(dx) * pow(ddx, 2.0) * 4.0 +
                   fabs(dy) * dirac(dy) * pow(ddy, 2.0) * 4.0 + helper_2) /
                  2;


  t_fw.x() = dx + d / sqrt(pow(helper_3, 2.0) + pow(helper_4, 2.0)) * (helper_5 - helper_7) + L * ddx / dP.norm() -
             d * (helper_3 * (helper_5 - helper_7) * 2.0 + helper_4 * (helper_6 - helper_7) * 2.0) *
                 pow(pow(helper_3, 2.0) + pow(helper_4, 2.0), 3.0 / 2.0) * helper_3 / 2 -
             L * helper_1 / pow(dP.norm(), 3) * dx / 2;

  t_fw.y() = dy + d / sqrt(pow(helper_3, 2.0) + pow(helper_4, 2.0)) * (helper_6 - helper_7) + L / dP.norm() * ddy -
             d *
                 (helper_3 * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) / helper_3)) *
                      (helper_5 - helper_7) * 2.0 +
                  helper_4 * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) / helper_4)) *
                      (helper_6 - helper_7) * 2.0) *
                 pow(pow(helper_3, 2.0) + pow(helper_4, 2.0), 3.0 / 2.0) *
                 (ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) / 2 -
             L * helper_1 / pow(dP.norm(), 3) * dy / 2;

  return t_fw;
}

#endif // FW_CURVE_H
