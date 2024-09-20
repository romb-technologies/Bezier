#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent), ui(new Ui::MainWindow), scene(new CustomScene)
{
  ui->setupUi(this);

  ui->graphicsView->setScene(scene);
  new QGraphicsViewZoom(ui->graphicsView);

  Eigen::MatrixX2d cp1, cp2;
  cp1.resize(4, 2);
  cp2.resize(5, 2);
  cp1 << 84, 162, 246, 30, 48, 236, 180, 110;

  cp2 << 180, 110, 175, 160, 60, 48, 164, 165, 124, 134;

  scene->addItem(new qCurve(cp1 * 5));
  scene->addItem(new qCurve(cp2 * 5));

  ui->graphicsView->centerOn(scene->itemsBoundingRect().center());
}

MainWindow::~MainWindow() { delete ui; }
