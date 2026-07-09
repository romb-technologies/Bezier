#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QFrame>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
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
  vbox->addWidget(lbl_selected_ = new QLabel);
  vbox->addWidget(lbl_order_ = new QLabel);
  vbox->addWidget(lbl_length_ = new QLabel);
  vbox->addWidget(lbl_offset_ = new QLabel);
  vbox->addStretch(1);

  vbox->addWidget(new QLabel("Continuity betas:"));
  beta_edit_ = new QLineEdit("1, 0, 0");
  beta_edit_->setToolTip("Comma-separated beta coefficients; their count sets the continuity order");
  btn_continuity_ = new QPushButton("Apply continuity");
  btn_continuity_->setToolTip("Reshapes the 2nd selected curve to continue from the 1st");
  vbox->addWidget(beta_edit_);
  vbox->addWidget(btn_continuity_);

  btn_draw_ = new QPushButton("Free-hand draw");
  btn_draw_->setCheckable(true);
  btn_cp_ = new QPushButton("Control points");
  btn_cp_->setCheckable(true);
  auto* btn_help = new QPushButton("Help");
  vbox->addWidget(btn_draw_);
  vbox->addWidget(btn_cp_);
  vbox->addWidget(btn_help);

  connect(btn_draw_, &QPushButton::toggled, this, [this](bool on) {
    if (on)
    {
      btn_cp_->setChecked(false);
      scene->setInputMode(CustomScene::InputMode::DrawFreehand);
    }
    else if (!btn_cp_->isChecked())
      scene->setInputMode(CustomScene::InputMode::Normal);
  });
  connect(btn_cp_, &QPushButton::toggled, this, [this](bool on) {
    if (on)
    {
      btn_draw_->setChecked(false);
      scene->setInputMode(CustomScene::InputMode::PlaceControlPoints);
    }
    else if (!btn_draw_->isChecked())
      scene->setInputMode(CustomScene::InputMode::Normal);
  });
  connect(btn_continuity_, &QPushButton::clicked, this, &MainWindow::applyContinuity);
  connect(btn_help, &QPushButton::clicked, this, [this] { scene->showHelp(); });
  connect(scene, &CustomScene::modeFinished, this, [this] {
    btn_draw_->setChecked(false);
    btn_cp_->setChecked(false);
  });

  ui->gridLayout->addWidget(panel, 0, 1);
  ui->gridLayout->setColumnStretch(0, 1);
}

void MainWindow::refreshDashboard()
{
  int n_curves = 0;
  for (auto* item : scene->items())
    if (qgraphicsitem_cast<qCurve*>(item))
      n_curves++;

  auto selected = scene->selectedItems();
  lbl_curves_->setText(QString("Curves: %1").arg(n_curves));
  lbl_selected_->setText(QString("Selected: %1").arg(selected.size()));
  btn_continuity_->setEnabled(selected.size() == 2);

  auto* single = selected.size() == 1 ? qgraphicsitem_cast<qCurve*>(selected.first()) : nullptr;
  if (single)
  {
    lbl_order_->setText(QString("Order: %1").arg(single->order()));
    lbl_length_->setText(QString("Length: %1").arg(single->length(), 0, 'f', 1));
  }
  else
  {
    lbl_order_->setText("Order: –");
    lbl_length_->setText("Length: –");
  }

  lbl_offset_->setText(QString("Offset: %1").arg(scene->offset(), 0, 'f', 0));
}

void MainWindow::applyContinuity()
{
  auto order = scene->selectionOrder();
  if (order.size() != 2)
    return;
  auto* source = qgraphicsitem_cast<qCurve*>(order.first());
  auto* target = qgraphicsitem_cast<qCurve*>(order.last());
  if (!source || !target)
    return;

  std::vector<double> betas;
  for (const auto& part : beta_edit_->text().split(',', Qt::SkipEmptyParts))
  {
    bool ok;
    betas.push_back(part.trimmed().toDouble(&ok));
    if (!ok)
    {
      QMessageBox::warning(this, "Warning", "Invalid beta values; enter comma-separated numbers.");
      return;
    }
  }
  if (betas.empty())
    return;

  target->prepareGeometryChange();
  target->applyContinuity(*source, betas);
  scene->update();
}
