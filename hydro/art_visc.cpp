#include "art_visc.h"

void art_visc(int ntotal,float *hsml,float *mass,float **x,float **vx,int niac,float *rho,
    float *c,int *pair_i,int *pair_j, float **dwdx,float **dvxdt,float *dedt)
{
    int i,j;
    float dx,piv,muv, vr, rr, h, mc, mrho, mhsml;
    float *dvx = new float[dim+1];
    // Parameter for the artificial viscosity:
    // Shear viscosity
    const int alpha = 1.e0;
    // Bulk viscosity
    const int beta = 1.e0;
    // Parameter to avoid singularities
    float etq = 0.1e0;

    for(i=1;i<=ntotal;i++)
    {
        for(int d=1;d<=dim;d++)
            dvxdt[d][i] = 0.e0;
        dedt[i] = 0.e0;
    }
    //Calculate SPH sum for artificial viscosity
    for(int k=1;k<=niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        mhsml= (hsml[i]+hsml[j])/2;
        vr = 0.e0;
        rr = 0.e0;
        for(int d=1;d<=dim;d++)
        {
            dvx[d] = vx[d][i] - vx[d][j];
            dx = x[d][i] - x[d][j];
            vr = vr + dvx[d]*dx;
            rr = rr + dx*dx;
        }
        //Artificial viscous force only if v_ij * r_ij < 0
        if (vr<0.e0)
        {
            //Calculate muv_ij = hsml v_lj * r_ij / ( r_ij"2 + hsml~2 etq~2 )
            muv = mhsml*vr/(rr + mhsml*mhsml*etq*etq);
            //Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij~2) / rho_ij
            mc = 0.5*(c[i] + c[j]);
            mrho = 0.5*(rho[i] + rho[j]);
            piv = (beta*muv - alpha*mc)*muv/mrho;
            //Calculate SPH sum for artificial viscous force
            for(int d=1;d<=dim;d++)
            {
                h = -piv*dwdx[d][k];
                dvxdt[d][i] = dvxdt[d][i] + mass[j]*h;
                dvxdt[d][j] = dvxdt[d][j] - mass[i]*h;
                dedt[i] = dedt[i] - mass[j]*dvx[d]*h;
                dedt[j] = dedt[j] - mass[i]*dvx[d]*h;
            }
        }
    }
    //Change of specific internal energy:
    for(i=1;i<=ntotal;i++)
        dedt[i] = 0.5*dedt[i];

    delete dvx;
}
