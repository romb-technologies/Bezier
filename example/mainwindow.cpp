#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow),
  scene(new CustomScene)
{
  ui->setupUi(this);

  ui->graphicsView->setScene(scene);
  new QGraphicsViewZoom(ui->graphicsView);

  Eigen::MatrixX2d cp1, cp2;

/*
  cp1.resize(4,2);
  cp2.resize(4,2);
  cp1 << 10, 100,
         90, 30,
         164, 254,
         220, 220;

  cp2 << 101, 237,
         66, 121,
         217, 43,c1 = new qCurve(split.first);
         237, 182;/*

  cp1.resize(4,2);
  cp2.resize(4,2);
  cp1 << 34.8, 5.1,
      31.7, 6.2,
      31.1, 9.0,
      34.8, 10.0;

  cp2 << 33, 7.5,
      32.5, 8.4,
      31.2, 7.8,
      31.6, 6.8;
*/
  cp1.resize(4,2);
  cp2.resize(5,2);
  cp1 << 84, 162,
      246, 3,
      48, 236,
      180, 110;

  cp2 << 105, 79,
      275, 210,
      16, 48,
      164, 165,
      128, 128;

  scene->curves.push_back(new qCurve(cp1*5));
  scene->curves.push_back(new qCurve(cp2*5));

  scene->addItem(scene->curves[0]);
  scene->addItem(scene->curves[1]);
}

MainWindow::~MainWindow()
{
  delete ui;
}
