#ifndef QPOLYCURVE_H
#define QPOLYCURVE_H

#include <QGraphicsItem>

#include "BezierCpp/declarations.h"
#include "BezierCpp/polycurve.h"

class qPolyCurve : public QGraphicsItem, public Bezier::PolyCurve
{
private:
  bool draw_control_points = false;
  bool draw_curvature_radious = false;
public:
  qPolyCurve(std::vector<Bezier::CurvePtr>& curve_list) : QGraphicsItem(), Bezier::PolyCurve(curve_list) {}
  qPolyCurve(Bezier::CurvePtr& curve) : QGraphicsItem(), Bezier::PolyCurve(curve) {}
  int type() const Q_DECL_OVERRIDE;
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) Q_DECL_OVERRIDE;
  QRectF boundingRect() const Q_DECL_OVERRIDE;
  void prepareGeometryChange() { QGraphicsItem::prepareGeometryChange(); }

  bool getDraw_control_points() const;
  void setDraw_control_points(bool value);
  bool getDraw_curvature_radious() const;
  void setDraw_curvature_radious(bool value);
};

#endif // QPOLYCURVE_H
