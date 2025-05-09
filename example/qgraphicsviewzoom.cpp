#include "qgraphicsviewzoom.h"

#include <QApplication>
#include <QMouseEvent>
#include <QScrollBar>
#include <qmath.h>

QGraphicsViewZoom::QGraphicsViewZoom(QGraphicsView* view) : QObject(view), m_view(view)
{
  m_view->viewport()->installEventFilter(this);
  m_view->setMouseTracking(true);
  m_modifiers = Qt::ControlModifier;
  m_zoom_factor_base = 1.0005;
}

void QGraphicsViewZoom::gentle_zoom(double factor)
{
  if (m_view->transform().m11() * factor < 1.0 / 30)
    return; // zoom x1/30 is minimum
  if (m_view->transform().m11() * factor > 100)
    return; // zoom x100 is maximum
  m_view->scale(factor, factor);
  m_view->centerOn(target_scene_pos);
  QPointF delta_viewport_pos =
      target_viewport_pos - QPointF(m_view->viewport()->width() / 2.0, m_view->viewport()->height() / 2.0);
  QPointF viewport_center = m_view->mapFromScene(target_scene_pos) - delta_viewport_pos;
  m_view->centerOn(m_view->mapToScene(viewport_center.toPoint()));

  emit zoomed();
}

void QGraphicsViewZoom::set_modifiers(Qt::KeyboardModifiers modifiers) { m_modifiers = modifiers; }

void QGraphicsViewZoom::set_zoom_factor_base(double value) { m_zoom_factor_base = value; }

bool QGraphicsViewZoom::eventFilter(QObject* object, QEvent* event)
{
  if (event->type() == QEvent::MouseMove)
  {
    QMouseEvent* mouse_event = static_cast<QMouseEvent*>(event);
    QPointF delta = target_viewport_pos - mouse_event->pos();
    if (qAbs(delta.x()) > 5 || qAbs(delta.y()) > 5)
    {
      target_viewport_pos = mouse_event->pos();
      target_scene_pos = m_view->mapToScene(mouse_event->pos());
    }
  }
  else if (event->type() == QEvent::Wheel)
  {
    QWheelEvent* wheel_event = static_cast<QWheelEvent*>(event);
    if (QApplication::keyboardModifiers() == m_modifiers)
    {
      if (wheel_event->orientation() == Qt::Vertical)
      {
        double angle = wheel_event->angleDelta().y();
        double factor = qPow(m_zoom_factor_base, angle);
        gentle_zoom(factor);
        return true;
      }
    }
  }
  Q_UNUSED(object)
  return false;
}
