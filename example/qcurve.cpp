#include "qcurve.h"

#include <QPainter>
#include <QPen>

void qCurve::setDraw_curvature_radius(bool value) { draw_curvature_radius = value; }

bool qCurve::getDraw_curvature_radius() const { return draw_curvature_radius; }

bool qCurve::getLocked() const
{
    return locked;
}

void qCurve::setLocked(bool value)
{
    locked = value;
}

int qCurve::type() const { return QGraphicsItem::UserType + 1; }

void qCurve::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
  Q_UNUSED(option)
  Q_UNUSED(widget)

  setFlag(GraphicsItemFlag::ItemIsSelectable, true);

  painter->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform, true);

  QPen pen;
  pen.setStyle(isSelected() ? Qt::DashDotLine : Qt::SolidLine);
  pen.setColor(getLocked() ? Qt::red : Qt::black);
  painter->setPen(pen);
  QPainterPath curve;
  auto poly = polyline();
  curve.moveTo(poly[0].x(), poly[0].y());
  for (uint k = 1; k < poly.size(); k++)
    curve.lineTo(poly[k].x(), poly[k].y());
  painter->drawPath(curve);

  if (isSelected()) // control points show while the curve is selected
  {
    const int d = 6;
    painter->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
    Bezier::PointVector points = controlPoints();
    for (uint k = 1; k < points.size(); k++)
    {
      painter->setPen(Qt::blue);
      painter->drawEllipse(QRectF(points[k - 1].x() - d / 2, points[k - 1].y() - d / 2, d, d));
      painter->setPen(QPen(QBrush(Qt::gray), 1, Qt::DotLine));
      painter->drawLine(QLineF(points[k - 1].x(), points[k - 1].y(), points[k].x(), points[k].y()));
    }
    painter->setPen(Qt::blue);
    painter->drawEllipse(QRectF(points.back().x() - d / 2, points.back().y() - d / 2, d, d));
  }

  if (draw_curvature_radius)
  {
    for (double t = 1.0 / 100; t <= 1.0; t += 1.0 / 200)
    {
      painter->setPen(QColor(abs(255 * (0.5 - t)), (int)(255 * t), (int)(255 * (1 - t))));
      auto p = valueAt(t);
      auto tangent = tangentAt(t);
      Bezier::Point normal(-tangent.y(), tangent.x());
      double kappa = curvatureAt(t);
      auto n1 = p + normal * kappa * 100;
      auto n2 = p - normal * kappa * 100;
      painter->drawLine(QLineF(n1.x(), n1.y(), n2.x(), n2.y()));
    }
  }
}

QRectF qCurve::boundingRect() const
{
  // Always cover the control points (even when not drawn) so toggling selection never leaves
  // ghost handles: the repaint region stays constant and large enough to erase them.
  auto bbox = boundingBox();
  QRectF rect(QPointF(bbox.min().x(), bbox.min().y()), QPointF(bbox.max().x(), bbox.max().y()));
  for (const auto& cp : controlPoints())
    rect |= QRectF(cp.x() - 3, cp.y() - 3, 6, 6);
  return rect;
}

QVariant qCurve::itemChange(GraphicsItemChange change, const QVariant& value)
{
  if (change == ItemSelectedHasChanged)
    update(); // selection toggles whether control points are drawn -> repaint
  return QGraphicsItem::itemChange(change, value);
}
