#include "art_heat.h"

void art_heat(int ntotal,float *hsml,float *mass,float **x,float **vx,int &niac,float *rho,
    float *u,float *c,int *pair_i,int *pair_j,float **dwdx,float *dedt)
{
    int i,j;
    float dx, rr, h, mrho, mhsml, hvcc, mui, muj, muij, rdwdx, g1, g2;
    float *dvx = new float[dim+1];
    float *vcc = new float[maxn];

    //Parameter for the artificial heat conduction:
    g1=0.1;
    g2=1.0;
    for(i=1;i<=ntotal;i++)
    {
        vcc[i] = 0.0;
        dedt[i] = 0.0;
    }

    for(int k=1;k<=niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        for(int d=1;d<=dim;d++)
            dvx[d] = vx[d][j] - vx[d][i];

        hvcc = dvx[1]*dwdx[1][k];
        for(int d=2;d<=dim;d++)
          hvcc = hvcc + dvx[d]*dwdx[d][k];

        vcc[i] = vcc[i] + mass[j]*hvcc/rho[j];
        vcc[j] = vcc[j] + mass[i]*hvcc/rho[i];
    }

    for(int k=1;k<=niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        mhsml= (hsml[i]+hsml[j])/2.;
        mrho = 0.5e0*(rho[i] + rho[j]);
        rr = 0.e0;
        rdwdx = 0.e0;
        for(int d=1;d<=dim;d++)
        {
            dx = x[d][i] - x[d][j];
            rr = rr + dx*dx;
            rdwdx = rdwdx + dx*dwdx[d][k];
        }
        mui=g1*hsml[i]*c[i] + g2*hsml[i]*hsml[i]*(fabs(vcc[i])-vcc[i]);//Посчитаем что это операция возведения в степень
        muj=g1*hsml[j]*c[j] + g2*hsml[j]*hsml[j]*(fabs(vcc[j])-vcc[j]);//только hsml. Т.к сама операция выполняется раньше других
        muij = 0.5*(mui+muj);
        h = muij/(mrho*(rr+0.01*mhsml*mhsml))*rdwdx;
        dedt[i] = dedt[i] + mass[j]*h*(u[i]-u[j]);
        dedt[j] = dedt[j] + mass[i]*h*(u[j]-u[i]);
    }

    for(i=1;i<=ntotal;i++)
        dedt[i] = 2.0e0*dedt[i];

    delete dvx;
    delete vcc;
}

