#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "customscene.h"

class QLabel;
class QLineEdit;
class QPushButton;

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  explicit MainWindow(QWidget* parent = nullptr);
  ~MainWindow();

private:
  void buildDashboard();
  void refreshDashboard();
  void applyContinuity();

  Ui::MainWindow* ui;
  CustomScene* scene;

  QLabel* lbl_curves_;
  QLabel* lbl_selected_;
  QLabel* lbl_order_;
  QLabel* lbl_length_;
  QLabel* lbl_offset_;

  QLineEdit* beta_edit_;
  QPushButton* btn_continuity_;
  QPushButton* btn_draw_;
  QPushButton* btn_cp_;
};

#endif // MAINWINDOW_H
