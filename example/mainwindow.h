#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QMessageBox>

#include "bezier.h"
#include "qgraphicsviewzoom.h"

namespace Ui {
class MainWindow;
}

inline QTextStream& qStdOut()
{
    static QTextStream ts( stdout );
    return ts;
}

class qCurve : public QGraphicsItem, public Bezier::Curve
{
public:
  qCurve(const Eigen::MatrixX2d &points) : QGraphicsItem(), Bezier::Curve(points){}
  qCurve(const Bezier::Curve &curve) : QGraphicsItem(), Bezier::Curve(curve){}
  qCurve(Bezier::Curve &&curve) : QGraphicsItem(), Bezier::Curve(curve){}
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) Q_DECL_OVERRIDE
  {
    painter->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform, true);

    QPainterPath curve;
    auto poly = getPolyline();
    curve.moveTo(poly[0].x(), poly[0].y());
    for(int k = 1; k < poly.size(); k++)
        curve.lineTo(poly[k].x(), poly[k].y());
    painter->drawPath(curve);

    const int d = 6;
    painter->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
    Bezier::PointVector points = getControlPoints();
    for(int k = 1; k < points.size(); k++)
    {
      painter->setPen(Qt::blue);
      painter->drawEllipse(points[k-1].x()-d/2, points[k-1].y()-d/2, d, d);
      painter->setPen(QPen(QBrush(Qt::gray), 1, Qt::DotLine));
      painter->drawLine(points[k-1].x(), points[k-1].y(), points[k].x(), points[k].y());
    }
    painter->setPen(Qt::blue);
    painter->drawEllipse(points.back().x()-d/2, points.back().y()-d/2, d, d);
  }

  QRectF boundingRect() const Q_DECL_OVERRIDE {
    auto bbox = getBBox(false);
    return QRectF(QPointF(bbox.min().x(), bbox.min().y()), QPointF(bbox.max().x(),bbox.max().y()));
  }

};

class CustomScene : public QGraphicsScene
{
private:
  // for projection
  QGraphicsEllipseItem *dot;
  QVector<QGraphicsLineItem *> line;
  QVector<QGraphicsLineItem *> tan;
  QVector<QGraphicsRectItem *> boxes;
  QVector<QGraphicsEllipseItem *> dots;
  bool draw_box_inter = false;
  bool show_projection = false;
  bool update_curvature = false;
  std::pair<qCurve *, double> t_to_update;
  bool update_cp = false;
  std::pair<qCurve *, int> cp_to_update;
public:
  QVector<qCurve *> curves;

  void updateBoxIter(){
    for(auto &box: boxes)
      removeItem(box);
    for(auto &dot: dots)
      removeItem(dot);
    if (draw_box_inter){
      for(auto &curve: curves){
        auto bbox = curve->getBBox(true);
        boxes.push_back(addRect(bbox.min().x(), bbox.min().y(),
                                bbox.max().x() - bbox.min().x(), bbox.max().y() - bbox.min().y(),
                                QPen(Qt::blue)));
      }
      for(int k = 0; k < curves.size() - 1; k++)
        for(int i = k + 1; i < curves.size(); i++){
          auto inter = curves[k]->getPointsOfIntersection(*curves[i]);
          for(auto &dot: inter)
            dots.push_back(addEllipse(dot.x()-3, dot.y()-3, 6, 6, QPen(Qt::red), QBrush(Qt::red,Qt::SolidPattern)));
        }
    }
  }

  void mousePressEvent(QGraphicsSceneMouseEvent *mouseEvent){
    const int sensitivity = 5;
    Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
    if(mouseEvent->button() == Qt::RightButton){
      dot = addEllipse(QRectF(QPointF(p.x(), p.y()),QSizeF(6, 6)),QPen(Qt::yellow), QBrush(Qt::red,Qt::SolidPattern));
      for(auto curve : curves){
        auto t1 = curve->projectPointOnCurve(p);
        auto p1 = curve->valueAt(t1);
        auto tan1 = curve->tangentAt(t1);
        line.push_back(addLine(QLineF(QPointF(p.x(),p.y()),QPointF(p1.x(),p1.y())),QPen(Qt::red)));
        tan.push_back(addLine(QLineF(QPointF(p1.x(),p1.y()) - 150*QPointF(tan1.x(),tan1.y()),
                                      QPointF(p1.x(),p1.y()) + 150*QPointF(tan1.x(),tan1.y())),QPen(Qt::blue)));
      }
      show_projection = true;
    }
    if(mouseEvent->button() == Qt::LeftButton){
      for(auto &curve: curves){
        auto pv = curve->getControlPoints();
        for(int k = 0; k < pv.size(); k++)
          if ((pv[k] - p).norm() < sensitivity){
            update_cp = true;
            cp_to_update = std::make_pair(curve, k);
          }
        if(update_cp) break;
        double t = curve->projectPointOnCurve(p);
        auto pt = curve->valueAt(t);
        if((pt - p).norm() < 10){
          update_curvature = true;
          t_to_update = std::make_pair(curve, t);
          break;
        }
      }
    }
    if(mouseEvent->button() == Qt::MiddleButton){
      for(auto &curve: curves){
        auto t = curve->projectPointOnCurve(p);
        if((curve->valueAt(t) - p).norm() < sensitivity){
          this->removeItem(curve);
          qCurve *c1, *c2;
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

  void mouseMoveEvent(QGraphicsSceneMouseEvent *mouseEvent){
    Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
    if(show_projection){
      dot->setRect(QRectF(QPointF(p.x()-3,p.y()-3),QSizeF(6, 6)));
      for(int k = 0; k < curves.size(); k++){
        auto t1 = curves[k]->projectPointOnCurve(Bezier::Point(p.x(), p.y()));
        auto p1 = curves[k]->valueAt(t1);
        auto tan1 = curves[k]->tangentAt(t1);
        line[k]->setLine(QLineF(QPointF(p.x(),p.y()),QPointF(p1.x(),p1.y())));
        tan[k]->setLine(QLineF(QPointF(p1.x(),p1.y()) - 500*QPointF(tan1.x(),tan1.y()),
                        QPointF(p1.x(),p1.y()) + 500*QPointF(tan1.x(),tan1.y())));
      }
    }
    if(update_cp){
      cp_to_update.first->manipulateControlPoint(cp_to_update.second, p);
      updateBoxIter();
      update();
    }
    if(update_curvature){
      try{
        t_to_update.first->manipulateCurvature(t_to_update.second, p);
        update();
      }
      catch(char const* err)
      {
        update_curvature = false;
        QMessageBox::warning(nullptr,"Warning", QString().sprintf("%s",err));
      }
      updateBoxIter();
    }
    QGraphicsScene::mouseMoveEvent(mouseEvent);
  }

  void mouseReleaseEvent(QGraphicsSceneMouseEvent *mouseEvent){
    if(mouseEvent->button() == Qt::RightButton){
      if(show_projection){
        removeItem(dot);
        for(int k = 0; k < curves.size(); k++){
          removeItem(line[k]);
          removeItem(tan[k]);
        }
        line.clear();
        tan.clear();
        show_projection = false;
      }
    }
    if(mouseEvent->button() == Qt::LeftButton){
      update_cp = false;
      update_curvature = false;
    }
    QGraphicsScene::mouseReleaseEvent(mouseEvent);
  }

  void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *mouseEvent){
    Bezier::Point p(mouseEvent->scenePos().x(), mouseEvent->scenePos().y());
    if(mouseEvent->button() == Qt::LeftButton)
      for(auto &curve: curves){
        double t = curve->projectPointOnCurve(p);
        auto pt = curve->valueAt(t);
        if((pt - p).norm() < 10){
          curve->elevateOrder();
        }
    }
    if(mouseEvent->button() == Qt::RightButton)
      for(auto &curve: curves){
        double t = curve->projectPointOnCurve(p);
        auto pt = curve->valueAt(t);
        if((pt - p).norm() < 10){
          try{
            curve->lowerOrder();
          }
          catch(char const *err)
          {
            QMessageBox::warning(nullptr,"Warning", QString().sprintf("%s",err));
          }
        }
    }
    if(mouseEvent->button() == Qt::MiddleButton){
      if(draw_box_inter) draw_box_inter = false;
      else draw_box_inter = true;
    }
    updateBoxIter();
    update();
    QGraphicsScene::mouseDoubleClickEvent(mouseEvent);
  }
};

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

private:
  Ui::MainWindow *ui;
  CustomScene *scene;
};




#endif // MAINWINDOW_H
