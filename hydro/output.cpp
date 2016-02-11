#include "output.h"

void output(float **x,float **vx,float *mass,float *rho,
  float *p,float *u,int *itype,float *hsml,int ntotal)
{
    FILE *out1, *out2, *out3;

    out1 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/f_xv.dat", "w");
    out2 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/f_state.dat", "w");
    out3 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/f_other.dat", "w");

    fprintf(out1, "%d\n", ntotal);
    for( int i=1;i<=ntotal;i++)
    {
        fprintf(out1, "%d ", i);
        for (int d = 1; d <= dim;d++)
            fprintf(out1, "%f %f",x[d][i], vx[d][i]);
        fprintf(out1, "\n");
        fprintf(out2, "%d %f %f %f %f \n", i, mass[i], rho[i], p[i], u[i]);
        fprintf(out3, "%d %d %f \n", i,itype[i], hsml[i]);
    }
    fclose(out1);
    fclose(out2);
    fclose(out3);
}

