#ifndef FW_CURVE_H
#define FW_CURVE_H

#include "bezier.h"

constexpr double L = 15;
constexpr double d = 0;

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

inline Bezier::Point getFrontWheelCurveDerivation_1(const Bezier::Curve& curve, double t)
{
  Bezier::Point d1_fw;

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
  auto helper_7 = dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + helper_2) / 2;

  d1_fw.x() = dx + d / sqrt(pow(helper_3, 2.0) + pow(helper_4, 2.0)) * (helper_5 - helper_7) + L * ddx / dP.norm() -
              d * (helper_3 * (helper_5 - helper_7) * 2.0 + helper_4 * (helper_6 - helper_7) * 2.0) *
                  pow(pow(helper_3, 2.0) + pow(helper_4, 2.0), 3.0 / 2.0) * helper_3 / 2 -
              L * helper_1 / pow(dP.norm(), 3) * dx / 2;

  d1_fw.y() = dy + d / sqrt(pow(helper_3, 2.0) + pow(helper_4, 2.0)) * (helper_6 - helper_7) + L / dP.norm() * ddy -
              d *
                  (helper_3 * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) / helper_3)) *
                       (helper_5 - helper_7) * 2.0 +
                   helper_4 * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) / helper_4)) *
                       (helper_6 - helper_7) * 2.0) *
                  pow(pow(helper_3, 2.0) + pow(helper_4, 2.0), 3.0 / 2.0) *
                  (ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) / 2 -
              L * helper_1 / pow(dP.norm(), 3) * dy / 2;

  return d1_fw;
}

inline Bezier::Point getFrontWheelCurveDerivation_2(const Bezier::Curve& curve, double t)
{
  Bezier::Point d2_fw;

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

  auto ddddP = curve.getDerivative()->getDerivative()->getDerivative()->getDerivative()->valueAt(t);
  double ddddx = ddddP.x();
  double ddddy = ddddP.y();

  auto helper_1 = dx * ddx * 2 + dy * ddy * 2;

  auto helper_3 = ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2;
  auto helper_4 = ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2;

  d2_fw
      .x() = -d / sqrt(pow(helper_3, 2) + pow(helper_4, 2)) *
                 (-ddddx / dP.norm() -
                  pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                  helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                  dx / pow(dP.norm(), 3) *
                      (pow((dx / fabs(dx)), 2.0) * ddx * dddx * 6.0 + dx * ddddx * 2.0 +
                       pow((dy / fabs(dy)), 2.0) * ddy * dddy * 6.0 + dy * ddddy * 2.0) /
                      2 +
                  ddx / pow(dP.norm(), 3) *
                      (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                       pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                       dx * dddx * 2.0 + dy * dddy * 2.0) *
                      (3.0 / 2.0) +
                  pow(dx * ddx * 2.0 + dy * ddy * 2.0, 3.0) / pow(pow(dx, 2) + pow(dy, 2), 7.0 / 2.0) * dx *
                      (1.5E1 / 8.0) -
                  helper_1 / pow(dP.norm(), 5) * dx *
                      (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                       pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                       dx * dddx * 2.0 + dy * dddy * 2.0) *
                      (9.0 / 4.0)) +
             L * dddx / dP.norm() + ddx -
             L / pow(dP.norm(), 3) * dx *
                 (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +
                  dx * dddx * 2.0 + dy * dddy * 2.0) *
                 0.5 -
             d / pow(pow(helper_3, 2) + pow(helper_4, 2), 3.0 / 2.0) *
                 (ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) *
                 (pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) /
                       fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2)),
                      2.0) *
                      pow(dddx / dP.norm() +
                              pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddx -
                              dx / pow(dP.norm(), 3) *
                                  (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                                   pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                                   dx * dddx * 2.0 + dy * dddy * 2.0) *
                                  0.5,
                          2.0) *
                      2.0 +
                  pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) /
                       fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2)),
                      2.0) *
                      pow(dddy / dP.norm() +
                              pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddy -
                              dy /
                                  pow(dP.norm(), 3) *
                                  (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                                   pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                                   dx *
                                       dddx * 2.0 +
                                   dy *
                                       dddy * 2.0) *
                                  0.5,
                          2.0) *
                      2.0 -
                  fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) *
                      (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2))) *
                      (-ddddx / dP.norm() -
                       pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                       helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                       dx /
                           pow(dP.norm(), 3) *
                           (pow((dx / fabs(dx)), 2.0) * ddx * dddx * 6.0 + dx * ddddx * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                           0.5 +
                       ddx / pow(dP.norm(), 3) *
                           (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                            dx *
                                dddx * 2.0 +
                            dy *
                                dddy * 2.0) *
                           (3.0 / 2.0) +
                       pow(dx * ddx * 2.0 + dy * ddy * 2.0, 3.0) / pow(pow(dx, 2) + pow(dy, 2), 7.0 / 2.0) * dx *
                           (1.5E1 / 8.0) -
                       helper_1 / pow(dP.norm(), 5) * dx *
                           (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                            dx *
                                dddx * 2.0 +
                            dy *
                                dddy * 2.0) *
                           (9.0 / 4.0)) *
                      2.0 -
                  fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) *
                      (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) /
                        fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2))) *
                      (-ddddy / dP.norm() -
                       pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) +
                       helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
                       dy /
                           pow(dP.norm(), 3) *
                           (pow((dx / fabs(dx)), 2.0) * ddx * dddx * 6.0 + dx * ddddx * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                           0.5 +
                       ddy / pow(dP.norm(), 3) *
                           (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (3.0 / 2.0) +
                       pow(dx * ddx * 2.0 + dy * ddy * 2.0, 3.0) / pow(pow(dx, 2) + pow(dy, 2), 7.0 / 2.0) * dy *
                           (1.5E1 / 8.0) -
                       helper_1 / pow(dP.norm(), 5) * dy *
                           (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (9.0 / 4.0)) *
                      2.0) /
                 2 -
             d * (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2))) * (dddx / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) / 2) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2))) * (dddy / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) / 2) * 2.0) / pow(pow(helper_3, 2) + pow(helper_4, 2), 3.0 / 2.0) *
                 (dddx / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                  helper_1 / pow(dP.norm(), 3) * ddx -
                  dx / pow(dP.norm(), 3) *
                      (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                       pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                       dx * dddx * 2.0 + dy * dddy * 2.0) *
                      0.5) +
             d * pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2))) * (dddx / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) / 2) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2))) * (dddy / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) / 2) * 2.0, 2.0) *
                 pow(pow(helper_3, 2) + pow(helper_4, 2), 5.0 / 2.0) *
                 (ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) * (3.0 / 4.0) +
             L * pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
             L * helper_1 / pow(dP.norm(), 3) * ddx;

  d2_fw
      .y() = -d / sqrt(pow(helper_3, 2) + pow(helper_4, 2)) *
                 (-ddddy / dP.norm() -
                  pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) +
                  helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
                  dy / pow(dP.norm(), 3) *
                      (pow((dx / fabs(dx)), 2.0) * ddx * dddx * 6.0 + dx * ddddx * 2.0 +
                       pow((dy / fabs(dy)), 2.0) * ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                      0.5 +
                  ddy / pow(dP.norm(), 3) *
                      (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                       pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                      (3.0 / 2.0) +
                  pow(dx * ddx * 2.0 + dy * ddy * 2.0, 3.0) / pow(pow(dx, 2) + pow(dy, 2), 7.0 / 2.0) * dy *
                      (1.5E1 / 8.0) -
                  helper_1 / pow(dP.norm(), 5) * dy *
                      (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                       pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                      (9.0 / 4.0)) +
             L * dddy / dP.norm() + ddy -
             L / pow(dP.norm(), 3) * dy *
                 (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +
                  dx * dddx * 2.0 + dy * dddy * 2.0) *
                 0.5 -
             d / pow(pow(helper_3, 2) + pow(helper_4, 2), 3.0 / 2.0) *
                 (ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) *
                 (pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) /
                       fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2)),
                      2.0) *
                      pow(dddx / dP.norm() +
                              pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddx -
                              dx / pow(dP.norm(), 3) *
                                  (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                                   pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                                   dx * dddx * 2.0 + dy * dddy * 2.0) *
                                  0.5,
                          2.0) *
                      2.0 +
                  pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) /
                       fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2)),
                      2.0) *
                      pow(dddy / dP.norm() +
                              pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddy -
                              dy /
                                  pow(dP.norm(), 3) *
                                  (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                                   pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                                   dx *
                                       dddx * 2.0 +
                                   dy *
                                       dddy * 2.0) *
                                  0.5,
                          2.0) *
                      2.0 -
                  fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) *
                      (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2))) *
                      (-ddddx / dP.norm() -
                       pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                       helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                       dx /
                           pow(dP.norm(), 3) *
                           (pow((dx / fabs(dx)), 2.0) * ddx * dddx * 6.0 + dx * ddddx * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                           0.5 +
                       ddx / pow(dP.norm(), 3) *
                           (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                            dx *
                                dddx * 2.0 +
                            dy *
                                dddy * 2.0) *
                           (3.0 / 2.0) +
                       pow(dx * ddx * 2.0 + dy * ddy * 2.0, 3.0) / pow(pow(dx, 2) + pow(dy, 2), 7.0 / 2.0) * dx *
                           (1.5E1 / 8.0) -
                       helper_1 / pow(dP.norm(), 5) * dx *
                           (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                            dx *
                                dddx * 2.0 +
                            dy *
                                dddy * 2.0) *
                           (9.0 / 4.0)) *
                      2.0 -
                  fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) *
                      (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) /
                        fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2))) *
                      (-ddddy / dP.norm() -
                       pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) +
                       helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
                       dy /
                           pow(dP.norm(), 3) *
                           (pow((dx / fabs(dx)), 2.0) * ddx * dddx * 6.0 + dx * ddddx * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                           0.5 +
                       ddy / pow(dP.norm(), 3) *
                           (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (3.0 / 2.0) +
                       pow(dx * ddx * 2.0 + dy * ddy * 2.0, 3.0) / pow(pow(dx, 2) + pow(dy, 2), 7.0 / 2.0) * dy *
                           (1.5E1 / 8.0) -
                       helper_1 / pow(dP.norm(), 5) * dy *
                           (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                            pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (9.0 / 4.0)) *
                      2.0) *
                 0.5 -
             d * (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2))) * (dddx / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) / 2) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2))) * (dddy / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) / 2) * 2.0) / pow(pow(helper_3, 2) + pow(helper_4, 2), 3.0 / 2.0) *
                 (dddy / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                  helper_1 / pow(dP.norm(), 3) * ddy -
                  dy / pow(dP.norm(), 3) *
                      (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 +
                       pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 +

                       dx * dddx * 2.0 + dy * dddy * 2.0) *
                      0.5) +
             d * pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx / 2))) * (dddx / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) / 2) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2))) * (dddy / dP.norm() + pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow((dx / fabs(dx)), 2.0) * pow(ddx, 2.0) * 2.0 + pow((dy / fabs(dy)), 2.0) * pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) / 2) * 2.0, 2.0) *
                 pow(pow(helper_3, 2) + pow(helper_4, 2), 5.0 / 2.0) *
                 (ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy / 2) * (3.0 / 4.0) +
             L * pow(dx * ddx * 2.0 + dy * ddy * 2.0, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
             L * helper_1 / pow(dP.norm(), 3) * ddy;

  return d2_fw;
}

#endif // FW_CURVE_H
