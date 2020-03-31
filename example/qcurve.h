#ifndef QCURVE_H
#define QCURVE_H

#include <QGraphicsItem>

#include "Bezier/declarations.h"
#include "Bezier/bezier.h"

class qCurve : public QGraphicsItem, public Bezier::Curve
{
private:
  bool draw_control_points = false;
  bool draw_curvature_radious = false;
  bool locked = false;

public:
  qCurve(const Eigen::MatrixX2d& points) : QGraphicsItem(), Bezier::Curve(points) {}
  qCurve(const Bezier::Curve& curve) : QGraphicsItem(), Bezier::Curve(curve) {}
  qCurve(Bezier::Curve&& curve) : QGraphicsItem(), Bezier::Curve(curve) {}

  int type() const Q_DECL_OVERRIDE;
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) Q_DECL_OVERRIDE;
  QRectF boundingRect() const Q_DECL_OVERRIDE;

  void prepareGeometryChange() { QGraphicsItem::prepareGeometryChange(); }
  void setDraw_control_points(bool value);
  void setDraw_curvature_radious(bool value);
  bool getDraw_control_points() const;
  bool getDraw_curvature_radious() const;
  std::shared_ptr<Bezier::Curve> getSharedPtr();
  bool getLocked() const;
  void setLocked(bool value);
};

#endif // QCURVE_H
