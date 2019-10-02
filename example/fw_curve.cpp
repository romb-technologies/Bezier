#include "fw_curve.h"

#define LL 20
#define dd 7

inline double LDx(Bezier::Point p) {
  auto asd = LL * p.x();
  auto asd2 = dd * p.y();
  return  asd - asd2;
}

inline double LDy(Bezier::Point p) {
  return LL * p.y() + dd * p.x(); }


Bezier::Point getFrontWheelPosition(const Bezier::Curve &curve, double t)
{
  Bezier::Point p_fw;

  auto P = curve.valueAt(t);
  auto d1 = curve.getDerivative()->valueAt(t);

  auto asd = LDx(d1);
  auto asd2 = LDy(d1);

  p_fw.x() = P.x() + LDx(d1) / d1.norm();
  p_fw.y() = P.y() + LDy(d1) / d1.norm();

  return p_fw;
}

Bezier::Point getFrontWheelCurveDerivation_1(const Bezier::Curve &curve, double t)
{
  Bezier::Point d1_fw;

  auto d1 = curve.getDerivative()->valueAt(t);
  auto d2 = curve.getDerivative()->getDerivative()->valueAt(t);

  auto n1 = d1.x() * d2.x() + d1.y() * d2.y();

  d1_fw.x() = d1.x() + LDx(d2) / d1.norm() - n1 * LDx(d1) / pow(d1.norm(), 3);
  d1_fw.y() = d1.y() + LDy(d2) / d1.norm() - n1 * LDy(d1) / pow(d1.norm(), 3);

  return d1_fw;
}

Bezier::Point getFrontWheelCurveDerivation_2(const Bezier::Curve &curve, double t)
{
  Bezier::Point d2_fw;

  auto d1 = curve.getDerivative()->valueAt(t);
  auto d2 = curve.getDerivative()->getDerivative()->valueAt(t);
  auto d3 = curve.getDerivative()->getDerivative()->getDerivative()->valueAt(t);

  auto n1 = d2.x() * d2.x() + d2.y() * d2.y() + d1.x() * d3.x() + d1.y() * d3.y();
  auto n2 = d1.x() * d2.x() + d1.y() * d2.y();

  d2_fw.x() = d2.x() + LDx(d3) / d1.norm() - (2 * n2 * LDx(d2) + n1 * LDx(d1)) / std::pow(d1.norm(), 3) +
      3 * n2 * n2 * LDx(d1) / std::pow(d1.norm(), 5);

  d2_fw.y() = d2.y() + LDy(d3) / d1.norm() - (2 * n2 * LDy(d2) + n1 * LDy(d1)) / std::pow(d1.norm(), 3) +
      3 * n2 * n2 * LDy(d1) / std::pow(d1.norm(), 5);

  return d2_fw;
}

Bezier::Point getFrontWheelCurveDerivation_3(const Bezier::Curve &curve, double t)
{
  Bezier::Point d3_fw;

  auto d1 = curve.getDerivative()->valueAt(t);
  auto d2 = curve.getDerivative()->getDerivative()->valueAt(t);
  auto d3 = curve.getDerivative()->getDerivative()->getDerivative()->valueAt(t);
  auto d4 = curve.getDerivative()->getDerivative()->getDerivative()->getDerivative()->valueAt(t);

  auto n1 = (d2.x() * d2.x() + d2.y() * d2.y()) + (d3.x() * d1.x() + d3.y() * d1.y());
  auto n2 = 3 * (d3.x() * d2.x() + d3.y() * d2.y()) + (d4.x() * d1.x() + d4.y() * d1.y());
  auto n3 = d2.x() * d1.x() + d2.y() * d1.y();

  d3_fw.x() = d3.x() +
      LDx(d4) / d1.norm() -
      (3 * n3 * LDx(d3) + 3 * n1 * LDx(d2) + n2 * LDx(d1)) / std::pow(d1.norm(), 3) +
      (9 * n3 * n3 * LDx(d2) + 9 * n3 * n1 * LDx(d1)) / std::pow(d1.norm(), 5) -
      (15 * n3 * n3 * n3 * LDx(d1)) / std::pow(d1.norm(), 7);

  d3_fw.y() = d3.y() +
      LDy(d4) / d1.norm() -
      (3 * n3 * LDy(d3) + 3 * n1 * LDy(d2) + n2 * LDy(d1)) / std::pow(d1.norm(), 3) +
      (9 * n3 * n3 * LDy(d2) + 9 * n3 * n1 * LDy(d1)) / std::pow(d1.norm(), 5) -
      (15 * n3 * n3 * n3 * LDy(d1)) / std::pow(d1.norm(), 7);

  return d3_fw;
}

double getFrontWheelCurveKappa(const Bezier::Curve &curve, double t)
{
  auto d1 = getFrontWheelCurveDerivation_1(curve, t);
  auto d2 = getFrontWheelCurveDerivation_2(curve, t);

  return std::fabs((d1.x() * d2.y() - d1.y() * d2.x()) / pow(d1.norm(), 3));
}

double getFrontWheelCurveKappaDerived(const Bezier::Curve &curve, double t)
{
  auto d1 = getFrontWheelCurveDerivation_1(curve, t);
  auto d2 = getFrontWheelCurveDerivation_2(curve, t);
  auto d3 = getFrontWheelCurveDerivation_3(curve, t);

  return d1.x() * d3.y() + d2.x() * d2.y() - d1.y() * d3.x() / pow(d1.norm(), 3) - d2.x() * d2.y() / pow(d1.norm(), 3) +
      (d1.x() * d2.x() + d1.y() * d2.y()) / pow(d1.norm(), 5) * d1.y() * d2.x() * 3.0;
}
