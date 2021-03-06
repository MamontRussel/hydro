#include "virt_part.h"

void virt_part(int itimestep, int ntotal, int &nvirt, float *hsml, float *mass, float **x,
    float **vx, float *rho, float *u, float *p, int *itype)
{
    FILE *in1, *in2, *in3;
    int i,im=0, mp;
    float x1, dx, v_inf;

    if (vp_input)//FALSE
    {
        in1 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/data/xv_vp.dat", "r");
        in2 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/data/state_vp.dat", "r");
        in3 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/data/other_vp.dat", "r");
        fscanf(in1, "%d", &nvirt);

        for (int j = 1; j <= nvirt; j++)
        {
            i = ntotal + j;
            fprintf(in1, "%d ", i);
            for (int d = 1; d <= dim; d++)
                fprintf(in1, "%f ",x[d][i]);
            for (int d = 1; d <= dim; d++)
                fprintf(in1, "%f ",vx[d][i]);
            fprintf(in1, "\n");
            fprintf(in2, "%d %f %f %f %f", im, mass[i], rho[i], p[i], u[i]);
            fprintf(in3, "%d %d %f", im, itype[i], hsml[i]);
        }
        fclose(in1);
        fclose(in2);
        fclose(in3);
    }
    else
    {
        nvirt = 0;
        mp = 40;
        x1 = 1.0e-3;
        dx = x1 / mp;
        v_inf = 1.e-3;

        //Monaghan type virtual particle on the Upper side
        for(int i=1;i<=2*mp+1;i++)
        {
            nvirt++;
            x[1][ntotal + nvirt] = (i-1)*dx/2;
            x[2][ntotal + nvirt] = x1;
            vx[1][ntotal + nvirt] = v_inf;
            vx[2][ntotal + nvirt] = 0.;
        }
        //Monaghan type virtual particle on the Lower side
        for(int i=1;i<=2*mp+1;i++)
        {
            nvirt++;
            x[1][ntotal + nvirt] = (i-1)*dx/2;
            x[2][ntotal + nvirt] = 0.;
            vx[1][ntotal + nvirt] = 0.;
            vx[2][ntotal + nvirt] = 0.;
        }
        //Monaghan type virtual particle on the Left side
        for(int i=1;i<=2*mp-1;i++)
        {
              nvirt = nvirt + 1;
              x[1][ntotal + nvirt] = 0.;
              x[2][ntotal + nvirt] = i*dx/2;
              vx[1][ntotal + nvirt] = 0.;
              vx[2][ntotal + nvirt] = 0.;
        }
        //Monaghan type virtual particle on the Right side
        for(int i=1;i<=2*mp-1;i++)
        {
              nvirt = nvirt + 1;
              x[1][ntotal + nvirt] = x1;
              x[2][ntotal + nvirt] = i*dx/2;
              vx[1][ntotal + nvirt] = 0.;
              vx[2][ntotal + nvirt] = 0.;
        }

        for(int i=1;i<=nvirt;i++)
        {
              rho[ntotal + i] = 1000.;
              mass[ntotal + i] = rho[ntotal + i] * dx * dx;
              p[ntotal + i] =0.;
              u[ntotal + i] = 357.1;
              itype[ntotal + i] = -2;
              hsml[ntotal + i] = dx;
        }
    }

    if ((itimestep%save_step)==0)
    {
        in1 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/data/xv_vp.dat", "w");
        in2 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/data/state_vp.dat", "w");
        in3 = fopen("/Users/Mamont/Documents/GitHub/hydro/hydro/data/other_vp.dat", "w");

        fprintf(in1, "%d\n", nvirt);
        for(int j=ntotal+1;j<=ntotal+nvirt;j++)
        {
            fprintf(in1, "%d ", j);
            for (int d = 1; d <= dim; d++)
                fprintf(in1, "%f ",x[d][j]);
            for (int d = 1; d <= dim; d++)
                fprintf(in1, "%f ",vx[d][j]);
            fprintf(in1, "\n");
            fprintf(in2, "%d %f %f %f %f \n", j, mass[j], rho[j], p[j], u[j]);
            fprintf(in3, "%d %d %f \n", j, itype[j], hsml[j]);
        }

        fclose(in1);
        fclose(in2);
        fclose(in3);
    }
    if ((itimestep%print_step)==0)
    {
        if (int_stat)
        {
            cout<<" Statistics: Virtual boundary particles:\n";
            cout<<" Number of virtual particles: "<<nvirt<<endl;
        }
    }
}
