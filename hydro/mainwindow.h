#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "ui_mainwindow.h"
#include <iostream>
#include <QTimer>
#include "modul.h"
#include "input.h"
#include "output.h"
#include "time_integration.h"

using namespace std;

// This is a SPH code, the followings are the
// basic parameters needed in this code or calculated by this code
// mass-- mass of particles [in]
// ntotal-- total particle number ues [in]
// dt Time step used in the time integration [in]
// itype-- types of particles [in]
// x-- coordinates of particles [in/out]
// vx-- velocities of particles [in/out]
// rho-- dnesities of particles [in/out]
// p-- pressure of particles [in/out]
// u-- internal energy of particles [in/out]
// hsml-- smoothing lengths of particles [in/out]
// c-- sound velocity of particles [out]
// s-- entropy of particles [out]
// e-- total energy of particles [out]

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    void paintEvent(QPaintEvent *);
    QPainter painter;
public slots:
    void updateBar();
private slots:
    void on_pushButton_clicked();
    void on_radioButton_2_clicked();
    void on_radioButton_clicked();

private:
    Ui::MainWindow ui;
    void drawPlots();
    void init();
    float **x,**vx,*mass,*rho,*p,*u,*c,*s,*e,*hsml,dt;
    int ntotal, maxtimestep;
    bool bResult;
    int *itype;
};

#endif // MAINWINDOW_H
