#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QDebug>
#include <chrono>


#include <iostream>
#include <iomanip>
#include "../include/Bezier/utils.h"

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

  using std::chrono::high_resolution_clock;
  using std::chrono::duration_cast;
  using std::chrono::duration;
  using std::chrono::microseconds;
  using std::chrono::milliseconds;

  // auto c1 = Bezier::Curve(cp1);
  // auto c2 = Bezier::Curve(cp2);
  // int N = 1000;

  // Eigen::ArrayXd l = Eigen::ArrayXd::Random(N).abs() * c2.length();
  // c2.step(0.0, c2.length()/2);
  // auto t1 = high_resolution_clock::now();
  // for (int k{}; k < N; k++)
  //   volatile auto x = c2.step2(0, l(k));
  // auto t2 = high_resolution_clock::now();
  // for (int k{}; k < N; k++)
  //   volatile auto x = c2.step(0, l(k));
  // auto t3 = high_resolution_clock::now();

  // std::cout << duration_cast<microseconds>(t2 - t1).count() << std::endl;
  // std::cout << duration_cast<microseconds>(t3 - t2).count() << std::endl;

  ui->graphicsView->centerOn(scene->itemsBoundingRect().center());
}

MainWindow::~MainWindow() { delete ui; }
