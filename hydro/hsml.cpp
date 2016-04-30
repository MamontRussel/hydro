#include "hsml.h"

void h_upgrade(float dt, int ntotal, float *mass, float **vx, float *rho, int niac, int *pair_i,
    int *pair_j, float **dwdx, float *hsml)
{
    int i, j;
    float fac, hvcc;
    float *dvx = new float[dim+1];
    float *vcc = new float[maxn];
    float *dhsml = new float[maxn];

    if (sle==0 )//Keep smoothing length unchanged. 2D
        return;
    else if (sle==2)
    {
        //	dh/dt = (-1/dim)*(h/rho)*(drho/dt).
        for(i=1;i<=ntotal;i++)
            vcc[i] = 0.e0;

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

        for(i=1;i<=ntotal;i++)
        {
            dhsml[i] = (hsml[i]/dim)*vcc[i];
            hsml[i] = hsml[i] + dt*dhsml[i];
            if (hsml[i]<=0) hsml[i] = hsml[i] - dt*dhsml[i];
        }
    }
    else if(sle==1)
    {
        fac = 2.0;
        for(i=1;i<=ntotal;i++)
            hsml[i] = fac * pow(mass[i]/rho[i],1./dim);
    }

    delete[] dvx;
    delete[] vcc;
    delete[] dhsml;
}
