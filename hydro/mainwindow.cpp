#include "mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    ui.setupUi(this);
    init();
}

void MainWindow::drawPlots()
{
    ui.plot1->addGraph();
    ui.plot1->graph(0)->setPen(QPen(Qt::blue)); // line color blue for first graph
    ui.plot1->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20))); // first graph will be filled with translucent blue
    ui.plot1->addGraph();
    ui.plot1->graph(1)->setPen(QPen(Qt::red)); // line color red for second graph
    // generate some points of data (y0 for first, y1 for second graph):
    QVector<double> x(250), y0(250), y1(250);
    for (int i=0; i<250; ++i)
    {
      x[i] = i;
      y0[i] = qExp(-i/150.0)*qCos(i/10.0); // exponentially decaying cosine
      y1[i] = qExp(-i/150.0);              // exponential envelope
    }
    // configure right and top axis to show ticks but no labels:
    // (see QCPAxisRect::setupFullAxesBox for a quicker method to do this)
    ui.plot1->xAxis2->setVisible(true);
    ui.plot1->xAxis2->setTickLabels(false);
    ui.plot1->yAxis2->setVisible(true);
    ui.plot1->yAxis2->setTickLabels(false);
    // make left and bottom axes always transfer their ranges to right and top axes:
    connect(ui.plot1->xAxis, SIGNAL(rangeChanged(QCPRange)), ui.plot1->xAxis2, SLOT(setRange(QCPRange)));
    connect(ui.plot1->yAxis, SIGNAL(rangeChanged(QCPRange)), ui.plot1->yAxis2, SLOT(setRange(QCPRange)));
    // pass data points to graphs:
    ui.plot1->graph(0)->setData(x, y0);
    ui.plot1->graph(1)->setData(x, y1);
    // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    ui.plot1->graph(0)->rescaleAxes();
    // same thing for graph 1, but only enlarge ranges (in case graph 1 is smaller than graph 0):
    ui.plot1->graph(1)->rescaleAxes(true);
    // Note: we could have also just called customPlot->rescaleAxes(); instead
    // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
    ui.plot1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

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

//    delete itype;
//    delete mass;
//    delete rho;
//    delete p;
//    delete u;
//    delete c;
//    delete s;
//    delete e;
//    delete hsml;
}

void MainWindow::on_pushButton_clicked()
{
    if (shocktube) dt = 0.005;
    if (shearcavity) dt = 5.e-5;

    input(x, vx, mass, rho, p, u, itype, hsml, ntotal);
    maxtimestep=ui.lineEdit->text().toInt();
    time_integration(x, vx, mass, rho, p, u,itype, hsml, ntotal, maxtimestep, dt);
    output(x, vx, mass, rho, p, u, itype, hsml, ntotal);
}
