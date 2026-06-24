#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFrame>
#include <QLabel>
#include <QVBoxLayout>

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent), ui(new Ui::MainWindow), scene(new CustomScene)
{
  ui->setupUi(this);

  ui->graphicsView->setScene(scene);
  new QGraphicsViewZoom(ui->graphicsView);

  buildDashboard();
  connect(scene, &QGraphicsScene::selectionChanged, this, &MainWindow::refreshDashboard);
  connect(scene, &QGraphicsScene::changed, this, &MainWindow::refreshDashboard);

  Eigen::MatrixX2d cp1, cp2;
  cp1.resize(4, 2);
  cp2.resize(5, 2);
  cp1 << 84, 162,
      246, 30,
      48, 236,
      180, 110;

  cp2 << 180, 110,
      175, 160,
      60, 48,
      164, 165,
      124, 134;

  scene->addItem(new qCurve(cp1 * 5));
  scene->addItem(new qCurve(cp2 * 5));

  ui->graphicsView->centerOn(scene->itemsBoundingRect().center());
  refreshDashboard();
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::buildDashboard()
{
  auto* panel = new QFrame;
  panel->setFrameShape(QFrame::StyledPanel);
  panel->setMinimumWidth(200);
  panel->setMaximumWidth(240);

  auto* vbox = new QVBoxLayout(panel);
  vbox->addWidget(lbl_curves_ = new QLabel);
  vbox->addWidget(lbl_polycurves_ = new QLabel);
  vbox->addWidget(lbl_selected_ = new QLabel);
  vbox->addWidget(lbl_order_ = new QLabel);
  vbox->addWidget(lbl_length_ = new QLabel);
  vbox->addWidget(lbl_offset_ = new QLabel);
  vbox->addStretch(1);
  // Buttons (draw mode, control-point mode, help) will be added here.

  ui->gridLayout->addWidget(panel, 0, 1);
  ui->gridLayout->setColumnStretch(0, 1);
}

void MainWindow::refreshDashboard()
{
  int n_curves = 0, n_polycurves = 0;
  for (auto* item : scene->items())
  {
    if (item->type() == QGraphicsItem::UserType + 1)
      n_curves++;
    else if (item->type() == QGraphicsItem::UserType + 2)
      n_polycurves++;
  }

  auto selected = scene->selectedItems();
  lbl_curves_->setText(QString("Curves: %1").arg(n_curves));
  lbl_polycurves_->setText(QString("Polycurves: %1").arg(n_polycurves));
  lbl_selected_->setText(QString("Selected: %1").arg(selected.size()));

  if (selected.size() == 1 && selected.first()->type() == QGraphicsItem::UserType + 1)
  {
    auto* curve = static_cast<qCurve*>(selected.first());
    lbl_order_->setText(QString("Order: %1").arg(curve->order()));
    lbl_length_->setText(QString("Length: %1").arg(curve->length(), 0, 'f', 1));
  }
  else
  {
    lbl_order_->setText("Order: –");
    lbl_length_->setText("Length: –");
  }

  lbl_offset_->setText(QString("Offset: %1").arg(scene->offset(), 0, 'f', 0));
}
