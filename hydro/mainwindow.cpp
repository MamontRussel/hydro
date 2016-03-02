#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    ui.setupUi(this);
}

void MainWindow::drawPlots()
{
    double min=10000,max=-10000,min1,max1,min2,max2,max3,min3;
    min1=min2=min3=min;
    max1=max2=max3=max;
    ui.plot1->addGraph();
    QVector<double> index(ntotal),rhod(ntotal);
    for (int i=130; i<370; ++i)
    {
        index[i] = i;
        rhod[i] = rho[i];
        if(rho[i]>max)max=rho[i];
        if(rho[i]<min)min=rho[i];
    }
    ui.plot1->graph(0)->setData(index, rhod);
    ui.plot1->xAxis->setRange(230,350);
    ui.plot1->yAxis->setRange(min,max);
    ui.plot1->replot();
    //ui.plot1->graph(0)->rescaleAxes(true);
    //ui.plot1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui.plot2->addGraph();
    QVector<double> vd(ntotal);
    for (int i=130; i<370; ++i)
    {
        vd[i] = vx[1][i];
        if(vx[1][i]>max1)max1=vx[1][i];
        if(vx[1][i]<min1)min1=vx[1][i];
    }
    ui.plot2->graph(0)->setData(index, vd);
    ui.plot2->xAxis->setRange(230,350);
    ui.plot2->yAxis->setRange(min1,max1);
    ui.plot2->replot();
    //ui.plot2->graph(0)->rescaleAxes(true);
    //ui.plot2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui.plot3->addGraph();
    QVector<double> ud(ntotal);
    for (int i=130; i<370; ++i)
    {
        ud[i] = u[i];
        if(u[i]>max2)max2=u[i];
        if(u[i]<min2)min2=u[i];
    }
    ui.plot3->graph(0)->setData(index, ud);
    ui.plot3->xAxis->setRange(230,350);
    ui.plot3->yAxis->setRange(min2,max2);
    ui.plot3->replot();
    //ui.plot3->graph(0)->rescaleAxes(true);
    //ui.plot3->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui.plot4->addGraph();
    QVector<double> pd(ntotal);
    for (int i=130; i<370; ++i)
    {
        pd[i] = p[i];
        if(p[i]>max3)max3=p[i];
        if(p[i]<min3)min3=p[i];
    }
    ui.plot4->graph(0)->setData(index, pd);
    ui.plot4->xAxis->setRange(230,350);
    ui.plot4->yAxis->setRange(min3,max3);
    ui.plot4->replot();
//    ui.plot4->graph(0)->rescaleAxes(true);
//    ui.plot4->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
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
