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

inline Bezier::Point getFrontWheelCurveDerivation_3(const Bezier::Curve& curve, double t)
{
  Bezier::Point d3_fw;

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

  auto dddddP = curve.getDerivative()->getDerivative()->getDerivative()->getDerivative()->getDerivative()->valueAt(t);
  double dddddx = dddddP.x();
  double dddddy = dddddP.y();

  auto helper_1 = dx * ddx * 2 + dy * ddy * 2;

  d3_fw.x() = d /
                  sqrt(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                       pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0)) *
                  (dddddx / dP.norm() - pow(helper_1, 3.0) / pow(dP.norm(), 7) * ddx * (1.5E1 / 2.0) +
                   dx / pow(dP.norm(), 5) *
                       pow(pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0, 2.0) *
                       (9.0 / 4.0) -
                   dx / pow(dP.norm(), 3) *
                       (pow(dddx, 2.0) * 6.0 + pow(dddy, 2.0) * 6.0 + ddx * ddddx * 8.0 + dx * dddddx * 2.0 +
                        ddy * ddddy * 8.0 + dy * dddddy * 2.0) *
                       0.5 +
                   pow(helper_1, 2.0) / pow(dP.norm(), 5) * dddx * (9.0 / 2.0) -
                   helper_1 / pow(dP.norm(), 3) * ddddx * 2.0 +
                   pow(helper_1, 4.0) / pow(dP.norm(), 9) * dx * (1.05E2 / 1.6E1) -
                   ddx / pow(dP.norm(), 3) *
                       (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 2.0 -
                   1.0 / pow(dP.norm(), 3) * dddx *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 3.0 -
                   pow(helper_1, 2.0) / pow(dP.norm(), 7) * dx *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (4.5E1 / 4.0) +
                   helper_1 / pow(dP.norm(), 5) * dx *
                       (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 3.0 +
                   helper_1 / pow(dP.norm(), 5) * ddx *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 9.0) +
              L / dP.norm() * ddddx +
              d *
                  pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                              (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                                fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                              (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                               helper_1 / pow(dP.norm(), 3) * ddx -
                               dx / pow(dP.norm(), 3) *
                                   (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                                   0.5) *
                              2.0 +
                          fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) *
                              (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                                fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) *
                              (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                               helper_1 / pow(dP.norm(), 3) * ddy -
                               dy / pow(dP.norm(), 3) *
                                   (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                                   0.5) *
                              2.0,
                      2.0) *
                  1.0 /
                  pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                          pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0),
                      5.0 / 2.0) *
                  (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                   helper_1 / pow(dP.norm(), 3) * ddx -
                   dx / pow(dP.norm(), 3) *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) *
                  (9.0 / 4.0) -
              d / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 3.0 / 2.0) *
                  (ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                  (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                       (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                         fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                       (dddddx / dP.norm() - pow(helper_1, 3.0) / pow(dP.norm(), 7) * ddx * (1.5E1 / 2.0) +
                        dx / pow(dP.norm(), 5) *
                            pow(pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0, 2.0) *
                            (9.0 / 4.0) -
                        dx / pow(dP.norm(), 3) *
                            (pow(dddx, 2.0) * 6.0 + pow(dddy, 2.0) * 6.0 + ddx * ddddx * 8.0 + dx * dddddx * 2.0 +
                             ddy * ddddy * 8.0 + dy * dddddy * 2.0) *
                            0.5 +
                        pow(helper_1, 2.0) / pow(dP.norm(), 5) * dddx * (9.0 / 2.0) -
                        helper_1 / pow(dP.norm(), 3) * ddddx * 2.0 +
                        pow(helper_1, 4.0) / pow(dP.norm(), 9) * dx * (1.05E2 / 1.6E1) -
                        ddx / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 2.0 -
                        1.0 / pow(dP.norm(), 3) * dddx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 3.0 -
                        pow(helper_1, 2.0) / pow(dP.norm(), 7) * dx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (4.5E1 / 4.0) +
                        helper_1 / pow(dP.norm(), 5) * dx *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 3.0 +
                        helper_1 / pow(dP.norm(), 5) * ddx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 9.0) *
                       2.0 -
                   pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5)),
                       2.0) *
                       (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                        helper_1 / pow(dP.norm(), 3) * ddx -
                        dx / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) *
                       (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                        helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                        dx / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                        ddx / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (3.0 / 2.0) +
                        pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) -
                        helper_1 / pow(dP.norm(), 5) * dx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (9.0 / 4.0)) *
                       6.0 -
                   pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                        fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)),
                       2.0) *
                       (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                        helper_1 / pow(dP.norm(), 3) * ddy -
                        dy / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) *
                       (-ddddy / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) +
                        helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
                        dy / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                        ddy / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (3.0 / 2.0) +
                        pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) -
                        helper_1 / pow(dP.norm(), 5) * dy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (9.0 / 4.0)) *
                       6.0 +
                   fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) *
                       (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                         fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) *
                       (dddddy / dP.norm() - pow(helper_1, 3.0) / pow(dP.norm(), 7) * ddy * (1.5E1 / 2.0) +
                        dy / pow(dP.norm(), 5) *
                            pow(pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0, 2.0) *
                            (9.0 / 4.0) -
                        dy / pow(dP.norm(), 3) *
                            (pow(dddx, 2.0) * 6.0 + pow(dddy, 2.0) * 6.0 + ddx * ddddx * 8.0 + dx * dddddx * 2.0 +
                             ddy * ddddy * 8.0 + dy * dddddy * 2.0) *
                            0.5 +
                        pow(helper_1, 2.0) / pow(dP.norm(), 5) * dddy * (9.0 / 2.0) -
                        helper_1 / pow(dP.norm(), 3) * ddddy * 2.0 +
                        pow(helper_1, 4.0) / pow(dP.norm(), 9) * dy * (1.05E2 / 1.6E1) -
                        ddy / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 2.0 -
                        1.0 / pow(dP.norm(), 3) * dddy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 3.0 -
                        pow(helper_1, 2.0) / pow(dP.norm(), 7) * dy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (4.5E1 / 4.0) +
                        helper_1 / pow(dP.norm(), 5) * dy *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 3.0 +
                        helper_1 / pow(dP.norm(), 5) * ddy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 9.0) *
                       2.0) *
                  0.5 +
              L * pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) -
              L * helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
              d * (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0) /
                  pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                          pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0),
                      3.0 / 2.0) *
                  (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                   helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                   dx / pow(dP.norm(), 3) *
                       (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                   ddx / pow(dP.norm(), 3) *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) +
                   pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) -
                   helper_1 / pow(dP.norm(), 5) * dx *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0)) *
                  (3.0 / 2.0) -
              L / pow(dP.norm(), 3) * dx * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                  0.5 -
              L / pow(dP.norm(), 3) * ddx * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) - d / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 3.0 / 2.0) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * (pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5)), 2.0) * pow(dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5, 2.0) * 2.0 + pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)), 2.0) * pow(dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5, 2.0) * 2.0 - fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) + helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) + dx / pow(dP.norm(), 3) * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 + ddx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) + pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) - helper_1 / pow(dP.norm(), 5) * dx * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0)) * 2.0 - fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (-ddddy / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) + helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) + dy / pow(dP.norm(), 3) * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 + ddy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) + pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) - helper_1 / pow(dP.norm(), 5) * dy * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0)) * 2.0) * (3.0 / 2.0) + dddx - d * pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0, 3.0) / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 7.0 / 2.0) * (ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (1.5E1 / 8.0) -
              L * pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) +
              L * helper_1 / pow(dP.norm(), 5) * dx *
                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0) +
              d * (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0) /
                  pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                          pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0),
                      5.0 / 2.0) *
                  (ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                  (pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5)),
                       2.0) *
                       pow(dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                               helper_1 / pow(dP.norm(), 3) * ddx -
                               dx / pow(dP.norm(), 3) *
                                   (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                                   0.5,
                           2.0) *
                       2.0 +
                   pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                        fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)),
                       2.0) *
                       pow(dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                               helper_1 / pow(dP.norm(), 3) * ddy -
                               dy / pow(dP.norm(), 3) *
                                   (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                                   0.5,
                           2.0) *
                       2.0 -
                   fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                       (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                         fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                       (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                        helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                        dx / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                        ddx / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (3.0 / 2.0) +
                        pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) -
                        helper_1 / pow(dP.norm(), 5) * dx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (9.0 / 4.0)) *
                       2.0 -
                   fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) *
                       (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                         fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) *
                       (-ddddy / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) +
                        helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
                        dy / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                        ddy / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (3.0 / 2.0) +
                        pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) -
                        helper_1 / pow(dP.norm(), 5) * dy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (9.0 / 4.0)) *
                       2.0) *
                  (9.0 / 4.0);

  d3_fw.y() = d /
                  sqrt(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                       pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0)) *
                  (dddddy / dP.norm() - pow(helper_1, 3.0) / pow(dP.norm(), 7) * ddy * (1.5E1 / 2.0) +
                   dy / pow(dP.norm(), 5) *
                       pow(pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0, 2.0) *
                       (9.0 / 4.0) -
                   dy / pow(dP.norm(), 3) *
                       (pow(dddx, 2.0) * 6.0 + pow(dddy, 2.0) * 6.0 + ddx * ddddx * 8.0 + dx * dddddx * 2.0 +
                        ddy * ddddy * 8.0 + dy * dddddy * 2.0) *
                       0.5 +
                   pow(helper_1, 2.0) / pow(dP.norm(), 5) * dddy * (9.0 / 2.0) -
                   helper_1 / pow(dP.norm(), 3) * ddddy * 2.0 +
                   pow(helper_1, 4.0) / pow(dP.norm(), 9) * dy * (1.05E2 / 1.6E1) -
                   ddy / pow(dP.norm(), 3) *
                       (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 2.0 -
                   1.0 / pow(dP.norm(), 3) * dddy *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 3.0 -
                   pow(helper_1, 2.0) / pow(dP.norm(), 7) * dy *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (4.5E1 / 4.0) +
                   helper_1 / pow(dP.norm(), 5) * dy *
                       (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 3.0 +
                   helper_1 / pow(dP.norm(), 5) * ddy *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 9.0) +
              L / dP.norm() * ddddy +
              d *
                  pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                              (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                                fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                              (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                               helper_1 / pow(dP.norm(), 3) * ddx -
                               dx / pow(dP.norm(), 3) *
                                   (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                                   0.5) *
                              2.0 +
                          fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) *
                              (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                                fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) *
                              (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                               helper_1 / pow(dP.norm(), 3) * ddy -
                               dy / pow(dP.norm(), 3) *
                                   (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                                   0.5) *
                              2.0,
                      2.0) *
                  1.0 /
                  pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                          pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0),
                      5.0 / 2.0) *
                  (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                   helper_1 / pow(dP.norm(), 3) * ddy -
                   dy / pow(dP.norm(), 3) *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) *
                  (9.0 / 4.0) -
              d / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 3.0 / 2.0) *
                  (ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) *
                  (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                       (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                         fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                       (dddddx / dP.norm() - pow(helper_1, 3.0) / pow(dP.norm(), 7) * ddx * (1.5E1 / 2.0) +
                        dx / pow(dP.norm(), 5) *
                            pow(pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0, 2.0) *
                            (9.0 / 4.0) -
                        dx / pow(dP.norm(), 3) *
                            (pow(dddx, 2.0) * 6.0 + pow(dddy, 2.0) * 6.0 + ddx * ddddx * 8.0 + dx * dddddx * 2.0 +
                             ddy * ddddy * 8.0 + dy * dddddy * 2.0) *
                            0.5 +
                        pow(helper_1, 2.0) / pow(dP.norm(), 5) * dddx * (9.0 / 2.0) -
                        helper_1 / pow(dP.norm(), 3) * ddddx * 2.0 +
                        pow(helper_1, 4.0) / pow(dP.norm(), 9) * dx * (1.05E2 / 1.6E1) -
                        ddx / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 2.0 -
                        1.0 / pow(dP.norm(), 3) * dddx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 3.0 -
                        pow(helper_1, 2.0) / pow(dP.norm(), 7) * dx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (4.5E1 / 4.0) +
                        helper_1 / pow(dP.norm(), 5) * dx *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 3.0 +
                        helper_1 / pow(dP.norm(), 5) * ddx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 9.0) *
                       2.0 -
                   pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5)),
                       2.0) *
                       (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                        helper_1 / pow(dP.norm(), 3) * ddx -
                        dx / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) *
                       (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                        helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                        dx / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                        ddx / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (3.0 / 2.0) +
                        pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) -
                        helper_1 / pow(dP.norm(), 5) * dx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (9.0 / 4.0)) *
                       6.0 -
                   pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                        fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)),
                       2.0) *
                       (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                        helper_1 / pow(dP.norm(), 3) * ddy -
                        dy / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) *
                       (-ddddy / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) +
                        helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
                        dy / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                        ddy / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (3.0 / 2.0) +
                        pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) -
                        helper_1 / pow(dP.norm(), 5) * dy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (9.0 / 4.0)) *
                       6.0 +
                   fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) *
                       (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                         fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) *
                       (dddddy / dP.norm() - pow(helper_1, 3.0) / pow(dP.norm(), 7) * ddy * (1.5E1 / 2.0) +
                        dy / pow(dP.norm(), 5) *
                            pow(pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0, 2.0) *
                            (9.0 / 4.0) -
                        dy / pow(dP.norm(), 3) *
                            (pow(dddx, 2.0) * 6.0 + pow(dddy, 2.0) * 6.0 + ddx * ddddx * 8.0 + dx * dddddx * 2.0 +
                             ddy * ddddy * 8.0 + dy * dddddy * 2.0) *
                            0.5 +
                        pow(helper_1, 2.0) / pow(dP.norm(), 5) * dddy * (9.0 / 2.0) -
                        helper_1 / pow(dP.norm(), 3) * ddddy * 2.0 +
                        pow(helper_1, 4.0) / pow(dP.norm(), 9) * dy * (1.05E2 / 1.6E1) -
                        ddy / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 2.0 -
                        1.0 / pow(dP.norm(), 3) * dddy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 3.0 -
                        pow(helper_1, 2.0) / pow(dP.norm(), 7) * dy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (4.5E1 / 4.0) +
                        helper_1 / pow(dP.norm(), 5) * dy *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 3.0 +
                        helper_1 / pow(dP.norm(), 5) * ddy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 9.0) *
                       2.0) *
                  0.5 +
              L * pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) -
              L * helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
              d * (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0) /
                  pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                          pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0),
                      3.0 / 2.0) *
                  (-ddddy / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) +
                   helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
                   dy / pow(dP.norm(), 3) *
                       (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                   ddy / pow(dP.norm(), 3) *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) +
                   pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) -
                   helper_1 / pow(dP.norm(), 5) * dy *
                       (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0)) *
                  (3.0 / 2.0) -
              L / pow(dP.norm(), 3) * dy * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                  0.5 -
              L / pow(dP.norm(), 3) * ddy * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) - d / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 3.0 / 2.0) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * (pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5)), 2.0) * pow(dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5, 2.0) * 2.0 + pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)), 2.0) * pow(dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5, 2.0) * 2.0 - fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) + helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) + dx / pow(dP.norm(), 3) * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 + ddx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) + pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) - helper_1 / pow(dP.norm(), 5) * dx * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0)) * 2.0 - fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (-ddddy / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) + helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) + dy / pow(dP.norm(), 3) * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 + ddy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) + pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) - helper_1 / pow(dP.norm(), 5) * dy * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0)) * 2.0) * (3.0 / 2.0) + dddy - d * pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0, 3.0) / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 7.0 / 2.0) * (ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (1.5E1 / 8.0) -
              L * pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) +
              L * helper_1 / pow(dP.norm(), 5) * dy *
                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0) +
              d * (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0) /
                  pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                          pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0),
                      5.0 / 2.0) *
                  (ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) *
                  (pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5)),
                       2.0) *
                       pow(dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                               helper_1 / pow(dP.norm(), 3) * ddx -
                               dx / pow(dP.norm(), 3) *
                                   (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                                   0.5,
                           2.0) *
                       2.0 +
                   pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                        fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)),
                       2.0) *
                       pow(dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                               helper_1 / pow(dP.norm(), 3) * ddy -
                               dy / pow(dP.norm(), 3) *
                                   (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                                   0.5,
                           2.0) *
                       2.0 -
                   fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                       (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                         fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                       (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                        helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                        dx / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                        ddx / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (3.0 / 2.0) +
                        pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) -
                        helper_1 / pow(dP.norm(), 5) * dx *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (9.0 / 4.0)) *
                       2.0 -
                   fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) *
                       (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                         fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) *
                       (-ddddy / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) +
                        helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) +
                        dy / pow(dP.norm(), 3) *
                            (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 +
                        ddy / pow(dP.norm(), 3) *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (3.0 / 2.0) +
                        pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) -
                        helper_1 / pow(dP.norm(), 5) * dy *
                            (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                            (9.0 / 4.0)) *
                       2.0) *
                  (9.0 / 4.0);

  return d3_fw;
}

inline double getFrontWheelCurveKappa(const Bezier::Curve& curve, double t)
{
  auto d1 = getFrontWheelCurveDerivation_1(curve, t);
  auto d2 = getFrontWheelCurveDerivation_2(curve, t);

  return std::fabs((d1.x() * d2.y() - d1.y() * d2.x()) / pow(d1.norm(), 3));
}

inline double getFrontWheelCurveKappaDerived(const Bezier::Curve& curve, double t)
{
  auto d1 = getFrontWheelCurveDerivation_1(curve, t);
  auto d2 = getFrontWheelCurveDerivation_2(curve, t);
  auto d3 = getFrontWheelCurveDerivation_3(curve, t);

  return d1.x() * d3.y() + d2.x() * d2.y() - d1.y() * d3.x() / pow(d1.norm(), 3) - d2.x() * d2.y() / pow(d1.norm(), 3) +
         (d1.x() * d2.x() * 2.0 + d1.y() * d2.y() * 2.0) / pow(d1.norm(), 5) * d1.y() * d2.x() * (3.0 / 2.0);
}

#endif // FW_CURVE_H
