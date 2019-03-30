#include "customscene.h"

void qCurve::setDraw_control_points(bool value) { draw_control_points = value; }

void qCurve::setDraw_curvature_radious(bool value) { draw_curvature_radious = value; }

bool qCurve::getDraw_control_points() const { return draw_control_points; }

void qCurve::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
{
  Q_UNUSED(option)
  Q_UNUSED(widget)

  setFlag(GraphicsItemFlag::ItemIsSelectable, true);

  painter->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform, true);

  painter->setPen(QPen(isSelected() ? Qt::DashDotLine : Qt::SolidLine));
  QPainterPath curve;
  auto poly = getPolyline();
  curve.moveTo(poly[0].x(), poly[0].y());
  for (uint k = 1; k < poly.size(); k++)
    curve.lineTo(poly[k].x(), poly[k].y());
  painter->drawPath(curve);

  if (draw_control_points)
  {
    const int d = 6;
    painter->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
    Bezier::PointVector points = getControlPoints();
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

  if (draw_curvature_radious)
  {
    painter->setPen(Qt::green);
    for (double t = 0; t <= 1.0; t += 1.0 / 500)
    {
      auto p = valueAt(t);
      auto n = p + normalAt(t) / curvatureAt(t);
      painter->drawLine(QLineF(p.x(), p.y(), n.x(), n.y()));
    }
  }
}

QRectF qCurve::boundingRect() const
{
  auto bbox = getBBox(false);
  return QRectF(QPointF(bbox.min().x(), bbox.min().y()), QPointF(bbox.max().x(), bbox.max().y()));
}

void CustomScene::drawForeground(QPainter* painter, const QRectF& rect)
{
  Q_UNUSED(rect)

  if (draw_box_)
  {
    painter->setPen(Qt::blue);
    for (auto&& curve : curves)
    {
      auto bbox = curve->getBBox(true);
      painter->drawRect(bbox.min().x(), bbox.min().y(), bbox.max().x() - bbox.min().x(),
                        bbox.max().y() - bbox.min().y());
    }
  }

  if (draw_inter_)
  {
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
      auto t1 = curve->projectPoint(p);
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
      if (mouseEvent->modifiers().testFlag(Qt::ControlModifier))
      {
        double t = curve->projectPoint(p);
        auto pt = curve->valueAt(t);
        if ((pt - p).norm() < 10)
          curve->setSelected(true);
      }
      else
      {
        for (auto&& item : selectedItems())
          item->setSelected(false);
        auto pv = curve->getControlPoints();
        for (uint k = 0; k < pv.size(); k++)
          if ((pv[k] - p).norm() < sensitivity && curve->getDraw_control_points())
          {
            update_cp = true;
            cp_to_update = std::make_pair(curve, k);
          }
        if (update_cp)
          break;
        double t = curve->projectPoint(p);
        auto pt = curve->valueAt(t);
        auto ep = curve->getEndPoints();
        if ((pt - p).norm() < 10 && (pt - ep.first).norm() > 20 && (pt - ep.second).norm() > 20)
        {
          update_curvature = true;
          t_to_update = std::make_pair(curve, t);
          break;
        }
      }
    }
  }
  if (mouseEvent->button() == Qt::MiddleButton)
  {
    for (auto&& curve : curves)
    {
      auto t = curve->projectPoint(p);
      if ((curve->valueAt(t) - p).norm() < sensitivity)
      {
        this->removeItem(curve);
        auto split = curve->splitCurve(t);
        delete curve;
        curves.removeOne(curve);
        qCurve* c1, *c2;
        c1 = new qCurve(split.first);
        c2 = new qCurve(split.second);
        this->addItem(c1);
        this->addItem(c2);
        curves.push_back(c1);
        curves.push_back(c2);
        update();
        break;
      }
    }
  }
}

void CustomScene::mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
  if (show_projection)
  {
    dot->setRect(QRectF(QPointF(p.x() - 3, p.y() - 3), QSizeF(6, 6)));
    for (int k = 0; k < curves.size(); k++)
    {
      auto t1 = curves[k]->projectPoint(Bezier::Point(p.x(), p.y()));
      auto p1 = curves[k]->valueAt(t1);
      auto tan1 = curves[k]->tangentAt(t1);
      line[k]->setLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())));
      tan[k]->setLine(QLineF(QPointF(p1.x(), p1.y()) - 500 * QPointF(tan1.x(), tan1.y()),
                             QPointF(p1.x(), p1.y()) + 500 * QPointF(tan1.x(), tan1.y())));
    }
  }
  if (update_cp)
  {
    cp_to_update.first->prepareGeometryChange();
    cp_to_update.first->manipulateControlPoint(cp_to_update.second, p);
    update();
  }
  if (update_curvature)
  {
    try
    {
      t_to_update.first->prepareGeometryChange();
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
}

void CustomScene::keyPressEvent(QKeyEvent* event)
{
  qDebug() << event->key();
  if (event->key() == 72) // key H
  {
    QMessageBox::information(nullptr, "Help", " Usage \n \
- starts with two Bezier curves (with 4 and 5 control points respectively) \n \
- Zoom in/out: *__Ctrl + mouse wheel__* \n \
- Manipulate control point or point on curve: *__Left mouse buttom__* \n \
- Project mouse pointer on all curves and show tangent: *__Right mouse buttom__* \n \
- Split curve at mouse point: *__Middle mouse buttom__* \n \
- Raise order of the curve: *__Double left click__* \n \
- Lower order of the curve *__Double right click__* \n \
- Toggle bounding boxes and curve intersections: *__Double middle click__*");
  }

  if (event->key() == 66) // key B
  {
    draw_box_ = !draw_box_;
    update();
  }
  if (event->key() == 73) // key I
  {
    draw_inter_ = !draw_inter_;
    update();
  }
  if (event->key() == 67) // key C
  {
    for (auto&& curve : selectedItems())
      static_cast<qCurve*>(curve)->setDraw_curvature_radious(true);
    if (selectedItems().empty())
      for (auto&& curve : curves)
        curve->setDraw_curvature_radious(false);
    update();
  }
  if (event->key() == 80) // key P
  {
    for (auto&& curve : selectedItems())
      static_cast<qCurve*>(curve)->setDraw_control_points(true);
    if (selectedItems().empty())
      for (auto&& curve : curves)
        curve->setDraw_control_points(false);
    update();
  }
  if (event->key() == 16777235) // key UP
  {
    for (auto&& curve : selectedItems())
      static_cast<qCurve*>(curve)->elevateOrder();
    update();
  }
  if (event->key() == 16777237) // key DOWN
  {
    for (auto&& curve : selectedItems())
      try
      {
          static_cast<qCurve*>(curve)->lowerOrder();
      }
      catch (char const* err)
      {
        QMessageBox::warning(nullptr, "Warning", QString().sprintf("%s", err));
      }
    update();
  }

  QGraphicsScene::keyPressEvent(event);
}
