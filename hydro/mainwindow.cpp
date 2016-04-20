#include "mainwindow.h"

// ============  Changeable PARAMETERS  =========================

// dim : Dimension of the problem (1, 2 or 3)
int dim=1;

// Switches for different senarios
// summation_density = .TRUE. : Use density summation model in the code,
// 			               .FALSE. : Use continuiity equation
// 	average_velocity = .TRUE. : Monaghan treatment on average velocity,
// 			               .FALSE.: No average treatment.
// 	config_input =.TRUE. : Load initial configuration data,
// 	            	.FALSE.: Generate initial configuration.
// 	virtual_part = .TRUE. : Use vritual particle,
// 			           .FALSE.: No use of vritual particle.
// 	vp_input = 	.TRUE. : Load virtual particle information,
// 		         	.FALSE.: Generate virtual particle information,
// 	visc = 		.true. : Consider viscosity,
// 		       	.false.: No viscosity.
//   ex_force =	.true. : Consider external force,
// 			        .false.: No external force.
// 	visc_artificial = .true. : Consider artificial viscosity,
// 			              .false.: No considering of artificial viscosity.
// 	heat_artificial = .true. : Consider artificial heating,
// 			              .false.: No considering of artificial heating,
// 	self_gravity = .true. : Considering self_gravity,
// 			           .false.: No considering of self_gravity
// 	nor_density = .true. : Density normalization by using CSPM,
// 			          .false.: No normalization.
bool summation_density=true, average_velocity=false, config_input=false;
bool virtual_part=false, vp_input=false, visc=false, ex_force=false, heat_artificial=false;
bool visc_artificial=true, self_gravity=false, nor_density=false;

// Simulation cases
// shocktube = .true. : carry out shock tube simulation
// shearcavity = .true. : carry out shear cavity simulation
bool shocktube = true, shearcavity = false;

// ==================================================

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent)
{
    ui.setupUi(this);
    ui.tabWidget_2->hide();
}

void MainWindow::drawPlots()
{
    if(dim==2)
    {
        drawIzolines();
        double rX,rY,minX=1000,minY=1000,maxX=-1000,maxY=-1000;
        QVector<double> x_,y_;
        for(int i=0;i<=41;i++)
        {
            int k=30 + (i-1)*40;
            if(i==0)
            {
                rX=0;
                rY=0;
            }
            else if(i==41)
            {
                rX=0.001;
                rY=0;
            }
            else
            {
                rX=x[1][k];
                rY=vx[2][k];
            }
            if(rX<minX)minX=rX;
            if(rX>maxX)maxX=rX;
            if(rY<minY)minY=rY;
            if(rY>maxY)maxY=rY;
            x_.push_back(rX*1000);
            y_.push_back(rY*1000);
        }

        ui.plot1_2->addGraph();
        ui.plot1_2->graph(0)->setData(x_, y_);
        ui.plot1_2->xAxis->setRange(minX*1000,maxX*1000);
        ui.plot1_2->yAxis->setRange(minY*1000,maxY*1000);
        ui.plot1_2->replot();
//        ui.plot1_2->graph(0)->rescaleAxes(true);
//        ui.plot1_2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);


        double minX2=1000,minY2=1000,maxX2=-1000,maxY2=-1000;
        QVector<double> x_2,y_2;
        for(int i=1;i<=41;i++)
        {
            if(i==41)
            {
                rX=0;
                rY=0;
            }
            else
            {
                rX=x[2][i+800];
                rY=vx[1][i+800];
            }
            if(rX<minX2)minX2=rX;
            if(rX>maxX2)maxX2=rX;
            if(rY<minY2)minY2=rY;
            if(rY>maxY2)maxY2=rY;
            x_2.push_back(rX*1000);
            y_2.push_back(rY*1000);
        }
        qSort(y_2);
        ui.plot2_2->addGraph();
        ui.plot2_2->graph(0)->setData(y_2,x_2);
        ui.plot2_2->xAxis->setRange(minX2*1000,maxX2*1000);
        ui.plot2_2->yAxis->setRange(minY2*1000,maxY2*1000);
        ui.plot2_2->replot();
        ui.plot2_2->graph(0)->rescaleAxes(true);
        ui.plot2_2->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    }
    else
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
}

void MainWindow::drawIzolines()
{


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

void MainWindow::on_radioButton_2_clicked()
{
    dim=1;
    summation_density=true;
    average_velocity=false;
    config_input=false;
    virtual_part=false;
    vp_input=false;
    visc=false;
    ex_force=false;
    visc_artificial=true;
    heat_artificial=false;
    self_gravity=false;
    nor_density=false;
    shocktube=true;
    shearcavity=false;
    ui.lineEdit->setText("22");
    ui.tabWidget->show();
    ui.tabWidget_2->hide();
}

void MainWindow::on_radioButton_clicked()
{
    dim=2;
    summation_density=true;
    average_velocity=true;
    config_input=false;
    virtual_part=true;
    vp_input=false;
    visc=true;
    ex_force=true;
    visc_artificial=false;
    heat_artificial=false;
    self_gravity=false;
    nor_density=true;
    shocktube=false;
    shearcavity=true;
    ui.lineEdit->setText("1000");
    ui.tabWidget->hide();
    ui.tabWidget_2->show();
}
