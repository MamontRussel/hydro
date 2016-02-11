#-------------------------------------------------
#
# Project created by QtCreator 2016-02-08T18:44:25
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = hydro
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    qcustomplot.cpp \
    input.cpp \
    output.cpp \
    time_integration.cpp \
    art_heat.cpp \
    art_visc.cpp \
    single_step.cpp \
    av_vel.cpp \
    density.cpp \
    direct_find.cpp \
    ext_force.cpp \
    hsml.cpp \
    int_force.cpp \
    virt_part.cpp \
    viscosity.cpp \
    eos.cpp \
    kernel.cpp

HEADERS  += mainwindow.h \
    qcustomplot.h \
    art_heat.h \
    art_visc.h \
    av_vel.h \
    density.h \
    direct_find.h \
    EOS.h \
    ext_force.h \
    hsml.h \
    int_force.h \
    kernel.h \
    modul.h \
    output.h \
    single_step.h \
    time_integration.h \
    virt_part.h \
    viscosity.h \
    input.h

FORMS    += mainwindow.ui
