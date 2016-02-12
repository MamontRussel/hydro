#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    ui.setupUi(this);
}

void MainWindow::drawPlots()
{
    ui.plot1->addGraph();
    QVector<double> index(ntotal),rhod(ntotal);
    for (int i=0; i<ntotal; ++i)
    {
        index[i] = i;
        rhod[i] = rho[i];
    }
    ui.plot1->graph(0)->setData(index, rhod);
    ui.plot1->graph(0)->rescaleAxes(true);
    ui.plot1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui.plot2->addGraph();
    QVector<double> vd(ntotal);
    for (int i=0; i<ntotal; ++i)
    {
        vd[i] = vx[1][i];
    }
    ui.plot2->graph(0)->setData(index, vd);
    ui.plot2->graph(0)->rescaleAxes(true);
    ui.plot2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui.plot3->addGraph();
    QVector<double> ud(ntotal);
    for (int i=0; i<ntotal; ++i)
    {
        ud[i] = u[i];
    }
    ui.plot3->graph(0)->setData(index, ud);
    ui.plot3->graph(0)->rescaleAxes(true);
    ui.plot3->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui.plot4->addGraph();
    QVector<double> pd(ntotal);
    for (int i=0; i<ntotal; ++i)
    {
        pd[i] = p[i];
    }
    ui.plot4->graph(0)->setData(index, pd);
    ui.plot4->graph(0)->rescaleAxes(true);
    ui.plot4->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

void MainWindow::init()
{
    itype = new int[maxn];
    x = new float*[dim+1];
    vx = new float*[dim+1];
    mass = new float[maxn];
    rho = new float[maxn];
    p = new float[maxn];
    u = new float[maxn];
    c = new float[maxn];
    s = new float[maxn];
    e = new float[maxn];
    hsml = new float[maxn];

    for (int i = 0; i <= dim; i++)
    {
        x[i] = new float[maxn];
        vx[i] = new float[maxn];
        for (int j = 0; j < maxn; j++)
        {
            x[i][j] = (float)NULL;
            vx[i][j] = (float)NULL;
            mass[j] = (float)NULL;
            rho[j] = (float)NULL;
            p[j] = (float)NULL;
            u[j] = (float)NULL;
            c[j] = (float)NULL;
            s[j] = (float)NULL;
            hsml[j] = (float)NULL;
         }
    }
}

void MainWindow::on_pushButton_clicked()
{
    init();
    ui.progressBar->setValue(0);
    if (shocktube) dt = 0.005;
    if (shearcavity) dt = 5.e-5;

    input(x, vx, mass, rho, p, u, itype, hsml, ntotal);
    maxtimestep=ui.lineEdit->text().toInt();
    time_integration(x, vx, mass, rho, p, u,itype, hsml, ntotal, maxtimestep, dt);
    output(x, vx, mass, rho, p, u, itype, hsml, ntotal);
    drawPlots();
    ui.progressBar->setValue(100);

    delete itype;
    delete mass;
    delete rho;
    delete p;
    delete u;
    delete c;
    delete s;
    delete e;
    delete hsml;
}
