#include "av_vel.h"

void av_vel(int ntotal,float *mass,int &niac,int *pair_i,int *pair_j,float *w,float **vx,
  float *rho,float **av)
{
    int i,j;
    float *dvx = new float[dim+1];
    float epsilon;

    // epsilon a small constants chosen by experence,
    // may lead to instability,
    // for example, for the 1 dimensional shock tube problem, the E <= 0.3
    epsilon = 0.3;

    for(i=1;i<=ntotal;i++)
        for(int d=1;d<=dim;d++)
            av[d][i] = 0.;

    for(int k=1;k<=niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        for(int d=1;d<=dim;d++)
        {
            dvx[d] = vx[d][i] - vx[d][j];
            av[d][i] = av[d][i] - 2*mass[j]*dvx[d]/(rho[i]+rho[j])*w[k];
            av[d][j] = av[d][j] + 2*mass[i]*dvx[d]/(rho[i]+rho[j])*w[k];
        }
    }

    for(i=1;i<=ntotal;i++)
        for(int d=1;d<=dim;d++)
            av[d][i] = epsilon * av[d][i];

    delete[] dvx;
}
