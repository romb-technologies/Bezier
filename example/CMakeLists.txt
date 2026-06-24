find_package(Qt5 REQUIRED COMPONENTS Core Gui Widgets)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

add_executable(bezier_example
  main.cpp
  mainwindow.cpp
  qgraphicsviewzoom.cpp
  customscene.cpp
  qcurve.cpp
  qpolycurve.cpp
  mainwindow.ui
)

target_link_libraries(bezier_example PRIVATE
  Bezier::Bezier
  Qt5::Core
  Qt5::Gui
  Qt5::Widgets
)
