#include "ext_force.h"

void ext_force(int ntotal,float **x,int &niac,int *pair_i,int *pair_j,int *itype,float **dvxdt)
{
    int i, j;
    float *dx = new float[dim+1];
    float rr, f, rr0, dd, p1, p2;

    for(int i=1;i<=ntotal;i++)
        for(int d=1;d<=dim;d++)
            dvxdt[d][i] = 0.;

    // Consider self-gravity or not ?
    if (self_gravity)
        for(int i=1;i<=ntotal;i++)
            dvxdt[dim][i] = -9.8;

    // Boundary particle force and penalty anti-penetration force.
    rr0 = 1.25e-5;
    dd = 1.e-2;
    p1 = 12;
    p2 = 4;

    for(int k=1;k<=niac;k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        if(itype[i]>0&&itype[j]<0)
        {
            rr = 0.;
            for(int d=1;d<=dim;d++)
            {
                dx[d] = x[d][i] - x[d][j];
                rr = rr + dx[d]*dx[d];
            }
            rr = sqrt(rr);
            if(rr<rr0)
            {
                f = (pow(rr0/rr,p1)-(pow(rr0/rr,p2)))/rr*rr;
                for(int d=1;d<=dim;d++)
                    dvxdt[d][i] = dvxdt[d][i] + dd*dx[d]*f;
            }
        }
    }

    delete dx;
}
