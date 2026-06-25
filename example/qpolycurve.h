#ifndef QPOLYCURVE_H
#define QPOLYCURVE_H

#include <QGraphicsItem>

#include "Bezier/declarations.h"
#include "Bezier/polycurve.h"
#include "Bezier/bezier.h"

class qPolyCurve : public QGraphicsItem, public Bezier::PolyCurve
{
private:
  bool draw_curvature_radious = false;
public:
  qPolyCurve(const std::deque<Bezier::Curve>& curve_list) : QGraphicsItem(), Bezier::PolyCurve(curve_list) {}
  qPolyCurve(const Bezier::Curve& curve) : QGraphicsItem(), Bezier::PolyCurve(std::deque<Bezier::Curve>{curve}) {}
  int type() const Q_DECL_OVERRIDE;
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) Q_DECL_OVERRIDE;
  QRectF boundingRect() const Q_DECL_OVERRIDE;
  QVariant itemChange(GraphicsItemChange change, const QVariant& value) Q_DECL_OVERRIDE;
  void prepareGeometryChange() { QGraphicsItem::prepareGeometryChange(); }

  bool getDraw_curvature_radious() const;
  void setDraw_curvature_radious(bool value);
};

#endif // QPOLYCURVE_H
