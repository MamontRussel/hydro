#include "density.h"

void sum_density(int ntotal, float *hsml, float *mass, int niac, int *pair_i,int *pair_j,
    float *w,float *rho)
{
    int i, j;
    float selfdens,r;
    float *hv = new float[dim+1];
    float *wi = new float[maxn];

    // wi(maxn) integration'of the kernel itself
    for(int d=1;d<=dim;d++)
        hv[d] = 0.e0;

    // Self density of each particle: Wii (Kernel for distance 0)
    // and take contribution of particle itself:
    r = 0;
    // Firstly calculate the integration of the kernel over the space
    for(i=1;i<=ntotal;i++)
    {
        kernel(r,hv,hsml[i],selfdens,hv);
        wi[i]=selfdens*mass[i]/rho[i];
    }

    for(int k=1;k<=niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        wi[i] = wi[i] + mass[j]/rho[j]*w[k];
        wi[j] = wi[j] + mass[i]/rho[i]*w[k];
    }

    // Secondly calculate the rho integration over the space
    for(i=1;i<=ntotal;i++)
    {
        kernel(r,hv,hsml[i],selfdens,hv);
        rho[i] = selfdens*mass[i];
    }

    // Calculate SPH sum for rho:
    for(int k=1;k<=niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        rho[i] = rho[i] + mass[j]*w[k];
        rho[j] = rho[j] + mass[i]*w[k];
    }

    // Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
    if (nor_density)
        for(i=1;i<=ntotal;i++)
            rho[i]=rho[i]/wi[i];

    delete hv;
    delete wi;
}

void con_density(int ntotal,float *mass,int niac,int *pair_i,int *pair_j,float **dwdx,
    float **vx,float *drhodt)
{
    int i,j;
    float vcc;
    float *dvx = new float[dim+1];
    for(i=1;i<=ntotal;i++)
        drhodt[i] = 0;

    for(int k=1;k<=niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        for(int d=1;d<=dim;d++)
          dvx[d] = vx[d][i] - vx[d][j];

        vcc = dvx[1]*dwdx[1][k];
        for(int d=2;d<=dim;d++)
          vcc = vcc + dvx[d]*dwdx[d][k];

        drhodt[i] = drhodt[i] + mass[j]*vcc;
        drhodt[j] = drhodt[j] + mass[i]*vcc;
    }

    delete dvx;
}
