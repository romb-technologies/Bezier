#include "customscene.h"

void qCurve::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
  Q_UNUSED(option)
  Q_UNUSED(widget)

  painter->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform, true);

  QPainterPath curve;
  auto poly = getPolyline();
  curve.moveTo(poly[0].x(), poly[0].y());
  for (uint k = 1; k < poly.size(); k++)
    curve.lineTo(poly[k].x(), poly[k].y());
  painter->drawPath(curve);

  const int d = 6;
  painter->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
  Bezier::PointVector points = getControlPoints();
  for (uint k = 1; k < points.size(); k++)
  {
    painter->setPen(Qt::blue);
    painter->drawEllipse(points[k - 1].x() - d / 2, points[k - 1].y() - d / 2, d, d);
    painter->setPen(QPen(QBrush(Qt::gray), 1, Qt::DotLine));
    painter->drawLine(points[k - 1].x(), points[k - 1].y(), points[k].x(), points[k].y());
  }
  painter->setPen(Qt::blue);
  painter->drawEllipse(points.back().x() - d / 2, points.back().y() - d / 2, d, d);
}

QRectF qCurve::boundingRect() const
{
  auto bbox = getBBox(false);
  return QRectF(QPointF(bbox.min().x(), bbox.min().y()), QPointF(bbox.max().x(), bbox.max().y()));
}

void CustomScene::drawForeground(QPainter* painter, const QRectF& rect)
{
  Q_UNUSED(rect)

  if (draw_box_inter)
  {
    painter->setPen(Qt::blue);
    for (auto&& curve : curves)
    {
      auto bbox = curve->getBBox(true);
      painter->drawRect(bbox.min().x(), bbox.min().y(), bbox.max().x() - bbox.min().x(),
                        bbox.max().y() - bbox.min().y());
    }

    painter->setPen(Qt::red);
    painter->setBrush(QBrush(Qt::red, Qt::SolidPattern));
    for (int k = 0; k < curves.size() - 1; k++)
      for (int i = k + 1; i < curves.size(); i++)
      {
        auto inter = curves[k]->getPointsOfIntersection(*curves[i]);
        for (auto& dot : inter)
          painter->drawEllipse(QPointF(dot.x(), dot.y()), 3, 3);
      }
  }
}

void CustomScene::mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  const int sensitivity = 5;
  Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
  if (mouseEvent->button() == Qt::RightButton)
  {
    dot = addEllipse(QRectF(QPointF(p.x(), p.y()), QSizeF(6, 6)), QPen(Qt::yellow), QBrush(Qt::red, Qt::SolidPattern));
    for (auto&& curve : curves)
    {
      auto t1 = curve->projectPointOnCurve(p);
      auto p1 = curve->valueAt(t1);
      auto tan1 = curve->tangentAt(t1);
      line.push_back(addLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())), QPen(Qt::red)));
      tan.push_back(addLine(QLineF(QPointF(p1.x(), p1.y()) - 150 * QPointF(tan1.x(), tan1.y()),
                                   QPointF(p1.x(), p1.y()) + 150 * QPointF(tan1.x(), tan1.y())),
                            QPen(Qt::blue)));
    }
    show_projection = true;
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
    for (auto&& curve : curves)
    {
      auto pv = curve->getControlPoints();
      for (uint k = 0; k < pv.size(); k++)
        if ((pv[k] - p).norm() < sensitivity)
        {
          update_cp = true;
          cp_to_update = std::make_pair(curve, k);
        }
      if (update_cp)
        break;
      double t = curve->projectPointOnCurve(p);
      auto pt = curve->valueAt(t);
      if ((pt - p).norm() < 10)
      {
        update_curvature = true;
        t_to_update = std::make_pair(curve, t);
        break;
      }
    }
  }
  if (mouseEvent->button() == Qt::MiddleButton)
  {
    for (auto&& curve : curves)
    {
      auto t = curve->projectPointOnCurve(p);
      if ((curve->valueAt(t) - p).norm() < sensitivity)
      {
        this->removeItem(curve);
        qCurve* c1, *c2;
        auto split = curve->splitCurve(t);
        c1 = new qCurve(split.first);
        c2 = new qCurve(split.second);
        this->addItem(c1);
        this->addItem(c2);
        curves.push_back(c1);
        curves.push_back(c2);
        curves.removeOne(curve);
        update();
        break;
      }
    }
  }
  QGraphicsScene::mousePressEvent(mouseEvent);
}

void CustomScene::mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
  if (show_projection)
  {
    dot->setRect(QRectF(QPointF(p.x() - 3, p.y() - 3), QSizeF(6, 6)));
    for (int k = 0; k < curves.size(); k++)
    {
      auto t1 = curves[k]->projectPointOnCurve(Bezier::Point(p.x(), p.y()));
      auto p1 = curves[k]->valueAt(t1);
      auto tan1 = curves[k]->tangentAt(t1);
      line[k]->setLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())));
      tan[k]->setLine(QLineF(QPointF(p1.x(), p1.y()) - 500 * QPointF(tan1.x(), tan1.y()),
                             QPointF(p1.x(), p1.y()) + 500 * QPointF(tan1.x(), tan1.y())));
    }
  }
  if (update_cp)
  {
    cp_to_update.first->manipulateControlPoint(cp_to_update.second, p);
    update();
  }
  if (update_curvature)
  {
    try
    {
      t_to_update.first->manipulateCurvature(t_to_update.second, p);
      update();
    }
    catch (char const* err)
    {
      update_curvature = false;
      QMessageBox::warning(nullptr, "Warning", QString().sprintf("%s", err));
    }
  }
  QGraphicsScene::mouseMoveEvent(mouseEvent);
}

void CustomScene::mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  if (mouseEvent->button() == Qt::RightButton)
  {
    if (show_projection)
    {
      removeItem(dot);
      for (int k = 0; k < curves.size(); k++)
      {
        removeItem(line[k]);
        removeItem(tan[k]);
      }
      line.clear();
      tan.clear();
      show_projection = false;
    }
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
    update_cp = false;
    update_curvature = false;
  }
  QGraphicsScene::mouseReleaseEvent(mouseEvent);
}

void CustomScene::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
  if (mouseEvent->button() == Qt::LeftButton)
    for (auto&& curve : curves)
    {
      double t = curve->projectPointOnCurve(p);
      auto pt = curve->valueAt(t);
      if ((pt - p).norm() < 10)
      {
        curve->elevateOrder();
      }
    }
  if (mouseEvent->button() == Qt::RightButton)
    for (auto&& curve : curves)
    {
      double t = curve->projectPointOnCurve(p);
      auto pt = curve->valueAt(t);
      if ((pt - p).norm() < 10)
      {
        try
        {
          curve->lowerOrder();
        }
        catch (char const* err)
        {
          QMessageBox::warning(nullptr, "Warning", QString().sprintf("%s", err));
        }
      }
    }
  if (mouseEvent->button() == Qt::MiddleButton)
  {
    if (draw_box_inter)
      draw_box_inter = false;
    else
      draw_box_inter = true;
  }
  update();
  QGraphicsScene::mouseDoubleClickEvent(mouseEvent);
}
