#-------------------------------------------------
#
# Project created by QtCreator 2014-03-21T13:16:29
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = pkexgtract_gui
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h

FORMS    += mainwindow.ui

#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../vis_studio/release/ -limageclasses
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../vis_studio/debug/ -limageclasses
#else:unix: LIBS += -L$$PWD/../../vis_studio/ -limageclasses

#INCLUDEPATH += $$PWD/../../src/imageclasses
#DEPENDPATH += $$PWD/../../src/imageclasses

#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../vis_studio/install/lib/ -lgdal_i
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../vis_studio/install/lib/ -lgdal_i
#else:unix: LIBS += -L$$PWD/../../../../vis_studio/install/lib/ -lgdal_i

#INCLUDEPATH += $$PWD/../../../../vis_studio/install/include
#DEPENDPATH += $$PWD/../../../../vis_studio/install/include

#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../vis_studio/release/ -lbase
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../vis_studio/debug/ -lbase
#else:unix: LIBS += -L$$PWD/../../vis_studio/ -lbase

#INCLUDEPATH += $$PWD/../../src
#DEPENDPATH += $$PWD/../../src/base
