#-------------------------------------------------
#
# Project created by QtCreator 2019-03-15T10:50:10
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = bezier_example
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++17
QMAKE_CXXFLAGS += -O3 -DNDEBUG

INCLUDEPATH += /usr/include/eigen3 \
               ../include

SOURCES += \
        main.cpp \
        mainwindow.cpp \
        qgraphicsviewzoom.cpp \
        customscene.cpp \
        qcurve.cpp \
        qpolycurve.cpp \
        ../src/utils.cpp \
        ../src/bezier.cpp \
        ../src/polycurve.cpp \

HEADERS += \
        ../include/Bezier/coefficients.h \
        mainwindow.h \
        qgraphicsviewzoom.h \
        customscene.h \
        qcurve.h \
        qpolycurve.h \
        ../include/Bezier/utils.h \
        ../include/Bezier/bezier.h \
        ../include/Bezier/polycurve.h \
        ../include/Bezier/declarations.h

FORMS += \
        mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
