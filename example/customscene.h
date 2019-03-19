#ifndef CUSTOMSCENE_H
#define CUSTOMSCENE_H

#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QMessageBox>

#include "bezier.h"
#include "qgraphicsviewzoom.h"

namespace Ui
{
class MainWindow;
}

inline QTextStream& qStdOut()
{
  static QTextStream ts(stdout);
  return ts;
}

class qCurve : public QGraphicsItem, public Bezier::Curve
{
public:
  qCurve(const Eigen::MatrixX2d& points) : QGraphicsItem(), Bezier::Curve(points) {}
  qCurve(const Bezier::Curve& curve) : QGraphicsItem(), Bezier::Curve(curve) {}
  qCurve(Bezier::Curve&& curve) : QGraphicsItem(), Bezier::Curve(curve) {}
  void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) Q_DECL_OVERRIDE;
  QRectF boundingRect() const Q_DECL_OVERRIDE;
};

class CustomScene : public QGraphicsScene
{
private:
  QGraphicsEllipseItem* dot;
  QVector<QGraphicsLineItem*> line;
  QVector<QGraphicsLineItem*> tan;
  bool draw_box_inter = false;
  bool show_projection = false;
  bool update_curvature = false;
  std::pair<qCurve*, double> t_to_update;
  bool update_cp = false;
  std::pair<qCurve*, uint> cp_to_update;

protected:
  void drawForeground(QPainter *painter, const QRectF &rect) Q_DECL_OVERRIDE;

  void mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;

  void mouseMoveEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;

  void mouseReleaseEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;

  void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* mouseEvent) Q_DECL_OVERRIDE;

public:
  QVector<qCurve*> curves;
};

#endif // CUSTOMSCENE_H
