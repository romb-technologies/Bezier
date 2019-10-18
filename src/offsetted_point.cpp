#include "BezierCpp/offsetted_point.h"
#include "BezierCpp/bezier.h"
#include "BezierCpp/polycurve.h"

inline double LDx(Bezier::Point p, Bezier::Point offset) { return offset.x() * p.x() - offset.y() * p.y(); }

inline double LDy(Bezier::Point p, Bezier::Point offset) { return offset.x() * p.y() + offset.y() * p.x(); }

namespace Bezier
{

template <> Bezier::Point getOffsettedPoint<Bezier::Curve>(const Bezier::Curve& curve, double t, Bezier::Point offset)
{
  Bezier::Point p_off;

  auto P = curve.valueAt(t);
  auto d1 = curve.getDerivativeAt(t);

  p_off.x() = P.x() + LDx(d1, offset) / d1.norm();
  p_off.y() = P.y() + LDy(d1, offset) / d1.norm();

  return p_off;
}

template <>
Bezier::Point getOffsettedPoint<Bezier::PolyCurve>(const Bezier::PolyCurve& poly_curve, double t, Bezier::Point offset)
{
  uint idx = static_cast<uint>(t);
  if (idx == poly_curve.getSize()) // for the last point of last curve
    --idx;

  return getOffsettedPoint<Bezier::Curve>(*poly_curve.getCurvePtr(idx), t - idx, offset);
}

template <>
Bezier::Point getOffsettedPointDerivation_1<Bezier::Curve>(const Bezier::Curve& curve, double t, Bezier::Point offset)
{
  Bezier::Point d1_off;

  auto d1 = curve.getDerivativeAt(t);
  auto d2 = curve.getDerivativeAt(2, t);

  auto n1 = d1.x() * d2.x() + d1.y() * d2.y();

  d1_off.x() = d1.x() + LDx(d2, offset) / d1.norm() - n1 * LDx(d1, offset) / pow(d1.norm(), 3);
  d1_off.y() = d1.y() + LDy(d2, offset) / d1.norm() - n1 * LDy(d1, offset) / pow(d1.norm(), 3);

  return d1_off;
}

template <>
Bezier::Point getOffsettedPointDerivation_1<Bezier::PolyCurve>(const Bezier::PolyCurve& poly_curve, double t,
                                                               Bezier::Point offset)
{
  uint idx = static_cast<uint>(t);
  if (idx == poly_curve.getSize()) // for the last point of last curve
    --idx;

  return getOffsettedPointDerivation_1<Bezier::Curve>(*poly_curve.getCurvePtr(idx), t - idx, offset);
}

template <>
Bezier::Point getOffsettedPointDerivation_2<Bezier::Curve>(const Bezier::Curve& curve, double t, Bezier::Point offset)
{
  Bezier::Point d2_off;

  auto d1 = curve.getDerivativeAt(t);
  auto d2 = curve.getDerivativeAt(2, t);
  auto d3 = curve.getDerivativeAt(3, t);

  auto n1 = d2.x() * d2.x() + d2.y() * d2.y() + d1.x() * d3.x() + d1.y() * d3.y();
  auto n2 = d1.x() * d2.x() + d1.y() * d2.y();

  d2_off.x() = d2.x() + LDx(d3, offset) / d1.norm() -
               (2 * n2 * LDx(d2, offset) + n1 * LDx(d1, offset)) / std::pow(d1.norm(), 3) +
               3 * n2 * n2 * LDx(d1, offset) / std::pow(d1.norm(), 5);

  d2_off.y() = d2.y() + LDy(d3, offset) / d1.norm() -
               (2 * n2 * LDy(d2, offset) + n1 * LDy(d1, offset)) / std::pow(d1.norm(), 3) +
               3 * n2 * n2 * LDy(d1, offset) / std::pow(d1.norm(), 5);

  return d2_off;
}

template <>
Bezier::Point getOffsettedPointDerivation_2<Bezier::PolyCurve>(const Bezier::PolyCurve& poly_curve, double t,
                                                               Bezier::Point offset)
{
  uint idx = static_cast<uint>(t);
  if (idx == poly_curve.getSize()) // for the last point of last curve
    --idx;

  return getOffsettedPointDerivation_2<Bezier::Curve>(*poly_curve.getCurvePtr(idx), t - idx, offset);
}

template <>
Bezier::Point getOffsettedPointDerivation_3<Bezier::Curve>(const Bezier::Curve& curve, double t, Bezier::Point offset)
{
  Bezier::Point d3_off;

  auto d1 = curve.getDerivativeAt(t);
  auto d2 = curve.getDerivativeAt(2, t);
  auto d3 = curve.getDerivativeAt(3, t);
  auto d4 = curve.getDerivativeAt(4, t);

  auto n1 = (d2.x() * d2.x() + d2.y() * d2.y()) + (d3.x() * d1.x() + d3.y() * d1.y());
  auto n2 = 3 * (d3.x() * d2.x() + d3.y() * d2.y()) + (d4.x() * d1.x() + d4.y() * d1.y());
  auto n3 = d2.x() * d1.x() + d2.y() * d1.y();

  d3_off.x() = d3.x() + LDx(d4, offset) / d1.norm() -
               (3 * n3 * LDx(d3, offset) + 3 * n1 * LDx(d2, offset) + n2 * LDx(d1, offset)) / std::pow(d1.norm(), 3) +
               (9 * n3 * n3 * LDx(d2, offset) + 9 * n3 * n1 * LDx(d1, offset)) / std::pow(d1.norm(), 5) -
               (15 * n3 * n3 * n3 * LDx(d1, offset)) / std::pow(d1.norm(), 7);

  d3_off.y() = d3.y() + LDy(d4, offset) / d1.norm() -
               (3 * n3 * LDy(d3, offset) + 3 * n1 * LDy(d2, offset) + n2 * LDy(d1, offset)) / std::pow(d1.norm(), 3) +
               (9 * n3 * n3 * LDy(d2, offset) + 9 * n3 * n1 * LDy(d1, offset)) / std::pow(d1.norm(), 5) -
               (15 * n3 * n3 * n3 * LDy(d1, offset)) / std::pow(d1.norm(), 7);

  return d3_off;
}

template <>
Bezier::Point getOffsettedPointDerivation_3<Bezier::PolyCurve>(const Bezier::PolyCurve& poly_curve, double t,
                                                               Bezier::Point offset)
{
  uint idx = static_cast<uint>(t);
  if (idx == poly_curve.getSize()) // for the last point of last curve
    --idx;

  return getOffsettedPointDerivation_3<Bezier::Curve>(*poly_curve.getCurvePtr(idx), t - idx, offset);
}

template <> double getOffsettedPointKappa<Bezier::Curve>(const Bezier::Curve& curve, double t, Bezier::Point offset)
{
  auto d1 = getOffsettedPointDerivation_1(curve, t, offset);
  auto d2 = getOffsettedPointDerivation_2(curve, t, offset);

  return std::fabs((d1.x() * d2.y() - d1.y() * d2.x()) / pow(d1.norm(), 3));
}

template <>
double getOffsettedPointKappa<Bezier::PolyCurve>(const Bezier::PolyCurve& poly_curve, double t, Bezier::Point offset)
{
  uint idx = static_cast<uint>(t);
  if (idx == poly_curve.getSize()) // for the last point of last curve
    --idx;

  return getOffsettedPointKappa<Bezier::Curve>(*poly_curve.getCurvePtr(idx), t - idx, offset);
}

template <>
double getOffsettedPointKappaDerived<Bezier::Curve>(const Bezier::Curve& curve, double t, Bezier::Point offset)
{
  auto d1 = getOffsettedPointDerivation_1(curve, t, offset);
  auto d2 = getOffsettedPointDerivation_2(curve, t, offset);
  auto d3 = getOffsettedPointDerivation_3(curve, t, offset);

  return d1.x() * d3.y() + d2.x() * d2.y() - d1.y() * d3.x() / pow(d1.norm(), 3) - d2.x() * d2.y() / pow(d1.norm(), 3) +
         (d1.x() * d2.x() + d1.y() * d2.y()) / pow(d1.norm(), 5) * d1.y() * d2.x() * 3.0;
}

template <>
double getOffsettedPointKappaDerived<Bezier::PolyCurve>(const Bezier::PolyCurve& poly_curve, double t,
                                                        Bezier::Point offset)
{
  uint idx = static_cast<uint>(t);
  if (idx == poly_curve.getSize()) // for the last point of last curve
    --idx;

  return getOffsettedPointKappaDerived<Bezier::Curve>(*poly_curve.getCurvePtr(idx), t - idx, offset);
}

} // namespace Bezier
