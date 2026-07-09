#include "customscene.h"

#include <QGraphicsPathItem>
#include <QGraphicsSceneMouseEvent>
#include <QKeyEvent>
#include <QMessageBox>

CustomScene::CustomScene(QObject* parent) : QGraphicsScene(parent)
{
  // Track the order in which curves were selected (selectedItems() is unordered)
  connect(this, &QGraphicsScene::selectionChanged, this, [this] {
    auto selected = selectedItems();
    for (auto it = selection_order_.begin(); it != selection_order_.end();)
      it = selected.contains(*it) ? it + 1 : selection_order_.erase(it);
    for (auto* item : selected)
      if (!selection_order_.contains(item))
        selection_order_.append(item);
  });
}

void CustomScene::drawForeground(QPainter* painter, const QRectF& rect)
{
  Q_UNUSED(rect)

  if (draw_box_)
  {
    painter->setPen(Qt::blue);
    for (auto&& item : items())
      if (auto* curve = qgraphicsitem_cast<qCurve*>(item))
      {
        auto bbox = curve->boundingBox();
        painter->drawRect(bbox.min().x(), bbox.min().y(), bbox.max().x() - bbox.min().x(),
                          bbox.max().y() - bbox.min().y());
      }
  }

  if (draw_inter_)
  {
    painter->setPen(Qt::red);
    painter->setBrush(QBrush(Qt::red, Qt::SolidPattern));
    auto all = items();
    for (int k = 0; k < all.size(); k++)
      for (int i = k; i < all.size(); i++)
      {
        auto* curve1 = qgraphicsitem_cast<qCurve*>(all[k]);
        auto* curve2 = qgraphicsitem_cast<qCurve*>(all[i]);
        if (!curve1 || !curve2)
          continue;
        for (auto& dot : curve2->intersections(*curve1))
          painter->drawEllipse(QPointF(dot.x(), dot.y()), 3, 3);
      }
  }

  if (offset_ != 0.0)
  {
    painter->setPen(QPen(Qt::darkGreen, 0, Qt::DashLine));
    painter->setBrush(Qt::NoBrush);
    for (auto&& item : items())
      if (auto* curve = qgraphicsitem_cast<qCurve*>(item))
      {
        const auto& poly = curve->offsetPolyline(offset_);
        QPainterPath path;
        path.moveTo(poly[0].x(), poly[0].y());
        for (uint k = 1; k < poly.size(); k++)
          path.lineTo(poly[k].x(), poly[k].y());
        painter->drawPath(path);
      }
  }

  if (draw_extrema_)
  {
    painter->setPen(Qt::magenta);
    painter->setBrush(QBrush(Qt::magenta, Qt::SolidPattern));
    for (auto&& item : items())
      if (auto* curve = qgraphicsitem_cast<qCurve*>(item))
        for (double t : curve->extrema())
        {
          auto e = curve->valueAt(t);
          painter->drawEllipse(QPointF(e.x(), e.y()), 4, 4);
        }
  }

  if (input_mode_ == InputMode::PlaceControlPoints)
  {
    painter->setPen(Qt::blue);
    painter->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
    for (const auto& cp : draw_pts_)
      painter->drawEllipse(QPointF(cp.x(), cp.y()), 3, 3);
  }
}

void CustomScene::mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  const int sensitivity = 5;
  Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());

  if (input_mode_ == InputMode::DrawFreehand && mouseEvent->button() == Qt::LeftButton)
  {
    draw_pts_ = {p};
    updatePreview();
    return;
  }
  if (input_mode_ == InputMode::PlaceControlPoints)
  {
    if (mouseEvent->button() == Qt::LeftButton)
    {
      draw_pts_.push_back(p);
      updatePreview();
    }
    else if (mouseEvent->button() == Qt::RightButton)
      finalizeControlPoints();
    return;
  }
  if (mouseEvent->button() == Qt::RightButton)
  {
    dot = addEllipse(QRectF(QPointF(p.x(), p.y()), QSizeF(6, 6)), QPen(Qt::yellow), QBrush(Qt::red, Qt::SolidPattern));
    for (auto&& item : items())
    {
      auto* curve = qgraphicsitem_cast<qCurve*>(item);
      if (!curve)
        continue;
      auto t1 = curve->projectPoint(p);
      auto p1 = curve->valueAt(t1);
      auto tan1 = curve->tangentAt(t1);
      line.insert(item, addLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())), QPen(Qt::red)));
      tan.insert(item, addLine(QLineF(QPointF(p1.x(), p1.y()) - 150 * QPointF(tan1.x(), tan1.y()),
                                      QPointF(p1.x(), p1.y()) + 150 * QPointF(tan1.x(), tan1.y())),
                               QPen(Qt::blue)));
      auto t2 = curve->step(t1, 50);
      auto a = curve->valueAt(t2);
      byLength.insert(item, addEllipse(QRectF(QPointF(a.x() - 3, a.y() - 3), QSizeF(6, 6)), QPen(Qt::yellow),
                                       QBrush(Qt::red, Qt::SolidPattern)));
    }
    show_projection = true;
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
    // 1) Grab a control point of an already-selected curve (control points show on selection)
    for (auto&& item : items())
    {
      auto* curve = qgraphicsitem_cast<qCurve*>(item);
      if (!curve || !curve->isSelected())
        continue;
      Bezier::PointVector pv = curve->controlPoints();
      for (uint k = 0; k < pv.size(); k++)
        if ((pv[k] - p).norm() < sensitivity)
        {
          update_cp = true;
          cp_to_update = std::make_pair(item, k);
        }
      if (update_cp)
        break;
    }
    if (update_cp)
      return;

    // 2) Otherwise select the curve under the cursor; Ctrl/Shift adds to the selection
    bool additive = mouseEvent->modifiers().testFlag(Qt::ControlModifier) ||
                    mouseEvent->modifiers().testFlag(Qt::ShiftModifier);
    QGraphicsItem* hit = nullptr;
    for (auto&& item : items())
      if (auto* curve = qgraphicsitem_cast<qCurve*>(item))
        if ((curve->valueAt(curve->projectPoint(p)) - p).norm() < 10)
        {
          hit = item;
          break;
        }
    if (!additive)
      for (auto&& item : selectedItems())
        item->setSelected(false);
    if (hit)
      hit->setSelected(true);
  }
  if (mouseEvent->button() == Qt::MiddleButton)
  {
    for (auto&& item : items())
      if (auto* curve = qgraphicsitem_cast<qCurve*>(item))
      {
        auto t = curve->projectPoint(p);
        if ((curve->valueAt(t) - p).norm() < sensitivity)
        {
          removeItem(item);
          auto split = curve->splitCurve(t);
          delete curve;
          addItem(new qCurve(split.first));
          addItem(new qCurve(split.second));
          update();
          break;
        }
      }
  }
}

void CustomScene::mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());

  if (input_mode_ != InputMode::Normal)
  {
    if (input_mode_ == InputMode::DrawFreehand && (mouseEvent->buttons() & Qt::LeftButton))
    {
      draw_pts_.push_back(p);
      updatePreview();
    }
    QGraphicsScene::mouseMoveEvent(mouseEvent);
    return;
  }

  if (show_projection)
  {
    dot->setRect(QRectF(QPointF(p.x() - 3, p.y() - 3), QSizeF(6, 6)));
    for (auto&& item : items())
    {
      auto* curve = qgraphicsitem_cast<qCurve*>(item);
      if (!curve)
        continue;
      auto t1 = curve->projectPoint(p);
      auto p1 = curve->valueAt(t1);
      auto tan1 = curve->tangentAt(t1);
      line[item]->setLine(QLineF(QPointF(p.x(), p.y()), QPointF(p1.x(), p1.y())));
      tan[item]->setLine(QLineF(QPointF(p1.x(), p1.y()) - 500 * QPointF(tan1.x(), tan1.y()),
                                QPointF(p1.x(), p1.y()) + 500 * QPointF(tan1.x(), tan1.y())));
      auto t2 = curve->step(t1, 50);
      auto a = curve->valueAt(t2);
      byLength[item]->setRect(QRectF(QPointF(a.x() - 3, a.y() - 3), QSizeF(6, 6)));
    }
  }
  if (update_cp)
  {
    if (auto* curve = qgraphicsitem_cast<qCurve*>(cp_to_update.first))
    {
      curve->prepareGeometryChange();
      curve->setControlPoint(cp_to_update.second, p);
    }
    update();
  }
  QGraphicsScene::mouseMoveEvent(mouseEvent);
}

void CustomScene::mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent)
{
  if (input_mode_ == InputMode::DrawFreehand && mouseEvent->button() == Qt::LeftButton)
  {
    finalizeFreehand();
    return;
  }
  if (input_mode_ != InputMode::Normal)
    return;

  if (mouseEvent->button() == Qt::RightButton)
  {
    if (show_projection)
    {
      removeItem(dot);
      for (auto&& item : items())
        if (qgraphicsitem_cast<qCurve*>(item))
        {
          removeItem(line[item]);
          removeItem(tan[item]);
          removeItem(byLength[item]);
        }
      line.clear();
      tan.clear();
      byLength.clear();
      show_projection = false;
    }
  }
  if (mouseEvent->button() == Qt::LeftButton)
  {
    update_cp = false;
  }
}

void CustomScene::keyPressEvent(QKeyEvent* keyEvent)
{
  if (keyEvent->key() == Qt::Key_H)
    showHelp();

  if (keyEvent->key() == Qt::Key_B)
  {
    draw_box_ = !draw_box_;
    update();
  }
  if (keyEvent->key() == Qt::Key_I)
  {
    draw_inter_ = !draw_inter_;
    update();
  }
  if (keyEvent->key() == Qt::Key_E)
  {
    draw_extrema_ = !draw_extrema_;
    update();
  }
  if (keyEvent->key() == Qt::Key_C)
  {
    for (auto&& item : selectedItems())
      if (auto* curve = qgraphicsitem_cast<qCurve*>(item))
        curve->setDraw_curvature_radius(!curve->getDraw_curvature_radius());
    update();
  }
  if (keyEvent->key() == Qt::Key_Up)
  {
    for (auto&& item : selectedItems())
      if (auto* curve = qgraphicsitem_cast<qCurve*>(item))
      {
        curve->prepareGeometryChange();
        curve->raiseOrder();
      }
    update();
  }
  if (keyEvent->key() == Qt::Key_Down)
  {
    for (auto&& item : selectedItems())
      if (auto* curve = qgraphicsitem_cast<qCurve*>(item))
        try
        {
          curve->prepareGeometryChange();
          curve->lowerOrder();
        }
        catch (const std::exception& err)
        {
          QMessageBox::warning(nullptr, "Warning", err.what());
        }
    update();
  }
  if (keyEvent->key() == Qt::Key_O)
    o_held_ = true;
  if (o_held_ && (keyEvent->key() == Qt::Key_Plus || keyEvent->key() == Qt::Key_Equal))
  {
    offset_ += 5;
    update();
  }
  if (o_held_ && keyEvent->key() == Qt::Key_Minus)
  {
    offset_ -= 5;
    update();
  }
  if (keyEvent->key() == Qt::Key_J) // fit a single new curve through two selected curves
  {
    if (selection_order_.size() != 2)
      return;
    auto* curve1 = qgraphicsitem_cast<qCurve*>(selection_order_.first());
    auto* curve2 = qgraphicsitem_cast<qCurve*>(selection_order_.last());
    if (curve1 && curve2)
    {
      auto* joined = new qCurve(Bezier::Curve::joinCurves(*curve1, *curve2));
      removeItem(curve1);
      removeItem(curve2);
      addItem(joined);
      update();
    }
  }
  if (keyEvent->key() == Qt::Key_Delete)
  {
    for (auto&& item : selectedItems())
    {
      QRectF vacated = item->sceneBoundingRect();
      removeItem(item);
      delete item;
      update(vacated); // repaint the vacated region (bare update() misses the last item)
    }
  }
}

void CustomScene::keyReleaseEvent(QKeyEvent* keyEvent)
{
  if (keyEvent->key() == Qt::Key_O && !keyEvent->isAutoRepeat())
    o_held_ = false;
}

void CustomScene::setInputMode(InputMode mode)
{
  draw_pts_.clear();
  if (preview_)
  {
    removeItem(preview_);
    delete preview_;
    preview_ = nullptr;
  }
  input_mode_ = mode;
  update();
}

void CustomScene::updatePreview()
{
  if (draw_pts_.empty())
    return;
  QPainterPath path;
  path.moveTo(draw_pts_[0].x(), draw_pts_[0].y());
  for (uint k = 1; k < draw_pts_.size(); k++)
    path.lineTo(draw_pts_[k].x(), draw_pts_[k].y());

  if (preview_)
    preview_->setPath(path);
  else
    preview_ = addPath(path, QPen(Qt::darkGray, 0, Qt::DashLine));
  update();
}

void CustomScene::finalizeFreehand()
{
  if (preview_)
  {
    removeItem(preview_);
    delete preview_;
    preview_ = nullptr;
  }
  if (draw_pts_.size() >= 2)
  {
    auto* curve = new qCurve(Bezier::Curve::fromPolyline(draw_pts_));
    addItem(curve);
    for (auto&& item : selectedItems())
      item->setSelected(false);
    curve->setSelected(true);
    emit modeFinished();
  }
  draw_pts_.clear();
  update();
}

void CustomScene::finalizeControlPoints()
{
  if (preview_)
  {
    removeItem(preview_);
    delete preview_;
    preview_ = nullptr;
  }
  if (draw_pts_.size() >= 2)
  {
    auto* curve = new qCurve(Bezier::Curve(draw_pts_));
    addItem(curve);
    for (auto&& item : selectedItems())
      item->setSelected(false);
    curve->setSelected(true);
    emit modeFinished();
  }
  draw_pts_.clear();
  update();
}

void CustomScene::showHelp()
{
  QMessageBox::information(nullptr, "Help", "\
Mouse controls:\n\
Right click - project mouse pointer on all curves\n\
Middle click - split curve under cursor in two\n\
Ctrl + Scroll - zoom in/out\n\
Left click - select curve under cursor (Ctrl/Shift to add to selection)\n\
Left click empty space - deselect\n\
Left click + drag control point - manipulate control point (shown while selected)\n\
\n\
Drawing modes (toggle from the side panel):\n\
Free-hand draw - drag to sketch a stroke, released as a fitted curve\n\
Control points - left click to place points, right click to finish\n\
\n\
Side panel:\n\
Apply continuity - reshape the 2nd selected curve to continue from the 1st,\n\
with given beta coefficients (their count sets the continuity order)\n\
\n\
Keyboard shortcuts:\n\
H - display help\n\
B - toggle bounding box display\n\
I - toggle intersections display\n\
E - toggle extrema display\n\
C - toggle curvature display (of selected curves)\n\
Key Up - raise the order (of selected curves)\n\
Key Down - lower the order (of selected curves)\n\
J - fit a single new curve through two selected curves\n\
O + +/- - increase/decrease offset distance (shown for all curves when != 0)\n\
Delete - delete curve");
}
