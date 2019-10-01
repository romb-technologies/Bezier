#ifndef FW_CURVE2_H
#define FW_CURVE2_H

#include "fw_curve.h"

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

  d3_fw
      .x() = d /
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
                              dy /
                                  pow(dP.norm(), 3) *
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
                  dx /
                      pow(dP.norm(), 3) *
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
                            ddy * ddddy * 8.0 + dy * dddddy * 2.0 +) *
                           0.5 +
                       pow(helper_1, 2.0) / pow(dP.norm(), 5) * dddx * (9.0 / 2.0) -
                       helper_1 / pow(dP.norm(), 3) * ddddx * 2.0 +
                       pow(helper_1, 4.0) / pow(dP.norm(), 9) * dx * (1.05E2 / 1.6E1) -
                       ddx / pow(dP.norm(), 3) *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 2.0 -
                       1.0 / pow(dP.norm(), 3) *
                           dddx * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           3.0 -
                       pow(helper_1, 2.0) / pow(dP.norm(), 7) * dx *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (4.5E1 / 4.0) +
                       helper_1 / pow(dP.norm(), 5) * dx *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 3.0 +
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
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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
                            ddy * ddddy * 8.0 + dy * dddddy * 2.0 +) *
                           0.5 +
                       pow(helper_1, 2.0) / pow(dP.norm(), 5) * dddy * (9.0 / 2.0) -
                       helper_1 / pow(dP.norm(), 3) * ddddy * 2.0 +
                       pow(helper_1, 4.0) / pow(dP.norm(), 9) * dy * (1.05E2 / 1.6E1) -
                       ddy / pow(dP.norm(), 3) *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 2.0 -
                       1.0 / pow(dP.norm(), 3) *
                           dddy * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           3.0 -
                       pow(helper_1, 2.0) / pow(dP.norm(), 7) * dy *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (4.5E1 / 4.0) +
                       helper_1 / pow(dP.norm(), 5) * dy *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 3.0 +
                       helper_1 / pow(dP.norm(), 5) * ddy *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 9.0) *
                      2.0) *
                 0.5 +
             L * pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) -
             L * helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) + d * (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0) / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 3.0 / 2.0) * (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) + helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) + dx / pow(dP.norm(), 3) * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 + ddx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) + pow(helper_1, 3.0) / pow(dP.norm(), 7) * dx * (1.5E1 / 8.0) - helper_1 / pow(dP.norm(), 5) * dx * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0)) * (3.0 / 2.0) -
             L / pow(dP.norm(), 3) * dx * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                 0.5 -
             L / pow(dP.norm(), 3) * ddx *
                 (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) -
             d / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 3.0 / 2.0) *
                 (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                  helper_1 / pow(dP.norm(), 3) * ddx -
                  dx /
                      pow(dP.norm(), 3) *
                      (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) *
                 (pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                       fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5)),
                      2.0) *
                      pow(dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddx -
                              dx /
                                  pow(dP.norm(), 3) *
                                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5,
                          2.0) *
                      2.0 +
                  pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                       fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)),
                      2.0) *
                      pow(dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddy -
                              dy /
                                  pow(dP.norm(), 3) *
                                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5,
                          2.0) *
                      2.0 -
                  fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                      (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                      (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                       helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                       dx /
                           pow(dP.norm(), 3) *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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
                       dy /
                           pow(dP.norm(), 3) *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
                       ddy / pow(dP.norm(), 3) *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (3.0 / 2.0) +
                       pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) -
                       helper_1 / pow(dP.norm(), 5) * dy *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (9.0 / 4.0)) *
                      2.0) *
                 (3.0 / 2.0) +
             dddx - d * pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0, 3.0) / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 7.0 / 2.0) * (ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (1.5E1 / 8.0) -
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
                                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5,
                          2.0) *
                      2.0 +
                  pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                       fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)),
                      2.0) *
                      pow(dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddy -
                              dy / pow(dP.norm(), 3) *
                                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5,
                          2.0) *
                      2.0 -
                  fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                      (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                      (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                       helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                       dx / pow(dP.norm(), 3) *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
                       ddy / pow(dP.norm(), 3) *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (3.0 / 2.0) +
                       pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) -
                       helper_1 / pow(dP.norm(), 5) * dy *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (9.0 / 4.0)) *
                      2.0) *
                 (9.0 / 4.0);

  d3_fw
      .y() = d /
                 sqrt(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) +
                      pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0)) *
                 (dddddy / dP.norm() - pow(helper_1, 3.0) / pow(dP.norm(), 7) * ddy * (1.5E1 / 2.0) +
                  dy / pow(dP.norm(), 5) *
                      pow(pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0, 2.0) *
                      (9.0 / 4.0) -
                  dy / pow(dP.norm(), 3) *
                      (pow(dddx, 2.0) * 6.0 + pow(dddy, 2.0) * 6.0 + ddx * ddddx * 8.0 + dx * dddddx * 2.0 +
                       ddy * ddddy * 8.0 + dy * dddddy * 2.0 +) *
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
                              dy /
                                  pow(dP.norm(), 3) *
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
                  dy /
                      pow(dP.norm(), 3) *
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
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 2.0 -
                       1.0 / pow(dP.norm(), 3) *
                           dddx * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           3.0 -
                       pow(helper_1, 2.0) / pow(dP.norm(), 7) * dx *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (4.5E1 / 4.0) +
                       helper_1 / pow(dP.norm(), 5) * dx *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 3.0 +
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
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 2.0 -
                       1.0 / pow(dP.norm(), 3) *
                           dddy * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           3.0 -
                       pow(helper_1, 2.0) / pow(dP.norm(), 7) * dy *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (4.5E1 / 4.0) +
                       helper_1 / pow(dP.norm(), 5) * dy *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 3.0 +
                       helper_1 / pow(dP.norm(), 5) * ddy *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 9.0) *
                      2.0) *
                 0.5 +
             L * pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) -
             L * helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) + d * (fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0) / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 3.0 / 2.0) * (-ddddy / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddy * (9.0 / 4.0) + helper_1 / pow(dP.norm(), 3) * dddy * (3.0 / 2.0) + dy / pow(dP.norm(), 3) * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) * 0.5 + ddy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) + pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) - helper_1 / pow(dP.norm(), 5) * dy * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (9.0 / 4.0)) * (3.0 / 2.0) -
             L / pow(dP.norm(), 3) * dy * (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0) *
                 0.5 -
             L / pow(dP.norm(), 3) * ddy *
                 (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * (3.0 / 2.0) -
             d / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 3.0 / 2.0) *
                 (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                  helper_1 / pow(dP.norm(), 3) * ddy -
                  dy /
                      pow(dP.norm(), 3) *
                      (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) *
                 (pow(((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                       fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5)),
                      2.0) *
                      pow(dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddx -
                              dx /
                                  pow(dP.norm(), 3) *
                                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5,
                          2.0) *
                      2.0 +
                  pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                       fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)),
                      2.0) *
                      pow(dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddy -
                              dy /
                                  pow(dP.norm(), 3) *
                                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5,
                          2.0) *
                      2.0 -
                  fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                      (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                      (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                       helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                       dx /
                           pow(dP.norm(), 3) *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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
                       dy /
                           pow(dP.norm(), 3) *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
                       ddy / pow(dP.norm(), 3) *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (3.0 / 2.0) +
                       pow(helper_1, 3.0) / pow(dP.norm(), 7) * dy * (1.5E1 / 8.0) -
                       helper_1 / pow(dP.norm(), 5) * dy *
                           (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) *
                           (9.0 / 4.0)) *
                      2.0) *
                 (3.0 / 2.0) +
             dddy - d * pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) * (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) / fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) * (dddx / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dx * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddx - dx / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0 + fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) / fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5))) * (dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) - helper_1 / pow(dP.norm(), 3) * ddy - dy / pow(dP.norm(), 3) * (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5) * 2.0, 3.0) / pow(pow(fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5), 2.0) + pow(fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5), 2.0), 7.0 / 2.0) * (ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) * (1.5E1 / 8.0) -
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
                                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5,
                          2.0) *
                      2.0 +
                  pow(((ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5) /
                       fabs(ddy / dP.norm() - helper_1 / pow(dP.norm(), 3) * dy * 0.5)),
                      2.0) *
                      pow(dddy / dP.norm() + pow(helper_1, 2.0) / pow(dP.norm(), 5) * dy * (3.0 / 4.0) -
                              helper_1 / pow(dP.norm(), 3) * ddy -
                              dy / pow(dP.norm(), 3) *
                                  (pow(ddx, 2.0) * 2.0 + pow(ddy, 2.0) * 2.0 + dx * dddx * 2.0 + dy * dddy * 2.0) * 0.5,
                          2.0) *
                      2.0 -
                  fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) *
                      (((ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5) /
                        fabs(ddx / dP.norm() - helper_1 / pow(dP.norm(), 3) * dx * 0.5))) *
                      (-ddddx / dP.norm() - pow(helper_1, 2.0) / pow(dP.norm(), 5) * ddx * (9.0 / 4.0) +
                       helper_1 / pow(dP.norm(), 3) * dddx * (3.0 / 2.0) +
                       dx / pow(dP.norm(), 3) *
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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
                           (ddx * dddx * 6.0 + dx * ddddx * 2.0 + ddy * dddy * 6.0 + dy * ddddy * 2.0 +) * 0.5 +
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

#endif // FW_CURVE2_H
