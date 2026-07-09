#ifndef QCURVE_H
#define QCURVE_H

#include <QGraphicsItem>

#include "Bezier/declarations.h"
#include "Bezier/bezier.h"

class qCurve : public QGraphicsItem, public Bezier::Curve
{
private:
  bool draw_curvature_radius = false;
  Bezier::PointVector offset_polyline_;
  double offset_value_ = 0.0;

public:
  enum
  {
    Type = QGraphicsItem::UserType + 1
  };

  qCurve(const Eigen::MatrixX2d& points) : QGraphicsItem(), Bezier::Curve(points) {}
  qCurve(const Bezier::Curve& curve) : QGraphicsItem(), Bezier::Curve(curve) {}
  qCurve(Bezier::Curve&& curve) : QGraphicsItem(), Bezier::Curve(std::move(curve)) {}

  int type() const Q_DECL_OVERRIDE;
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) Q_DECL_OVERRIDE;
  QRectF boundingRect() const Q_DECL_OVERRIDE;
  QVariant itemChange(GraphicsItemChange change, const QVariant& value) Q_DECL_OVERRIDE;

  void prepareGeometryChange()
  {
    offset_polyline_.clear();
    QGraphicsItem::prepareGeometryChange();
  }
  const Bezier::PointVector& offsetPolyline(double offset);
  void setDraw_curvature_radius(bool value);
  bool getDraw_curvature_radius() const;
};

#endif // QCURVE_H
