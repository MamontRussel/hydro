// Function to determine the information of virtual particles
// Here only the Monaghan type virtual particles for the 2D shear cavity driven problem are generated.
// itimestep : Current time step [in]
// ntotal Number of particles
// nvirt Number of virtual particles
// hsml Smoothing Length
// mass Particle masses
// x Coordinates of all particles
// vx Velocities of all particles
// rho Density
// u internal energy
// itype: type of particles

#include <cstdio>

void virt_part(int &itimestep,int ntotal,int &nvirt,double *hsml,double *mass,double **x,
	double **vx,double *rho,double *u,double *p,int *itype)
{
	FILE *in1, *in2, *in3;
	int i,im, mp;
	double x1, dx, v_inf;

	if (vp_input)
	{
		in1 = fopen("xv_vp.dat", "r");
		in2 = fopen("state_vp.da", "r");
		in3 = fopen("other_vp.dat", "r");
		fscanf(in1, "%d", &nvirt);

		for (int j = 1; j <= nvirt; j++)
		{
			i = ntotal + j;
			for (int d = 1; d <= dim; d++)
				fprintf(in1, "%d %f %f ", im, x[d][i], vx[d][i]);
			fprintf(in2, "%d %f %f %f", im, mass[i], rho[i], p[i], u[i]);
			fprintf(in3, "%d %d %f ", im, itype[i], hsml[i]);
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
		for(int j=1;j<=2*mp+1;j++)
		{
			nvirt++;
			x[1][ntotal + nvirt] = (j-1)*dx/2;
			x[2][ntotal + nvirt] = x1;
			vx[1][ntotal + nvirt] = v_inf;
			vx[2][ntotal + nvirt] = 0.;
		}
		//Monaghan type virtual particle on the Lower side
		for(int j=1;j<=2*mp+1;j++)
		{
			nvirt++;
			x[1][ntotal + nvirt] = (j-1)*dx/2;
			x[2][ntotal + nvirt] = 0.;
			vx[1][ntotal + nvirt] = 0.;
			vx[2][ntotal + nvirt] = 0.;
		}
		//Monaghan type virtual particle on the Left side
		for(int j=1;j<=2*mp-1;j++)
		{
			  nvirt = nvirt + 1;
			  x[1][ntotal + nvirt] = 0.;
			  x[2][ntotal + nvirt] = j*dx/2;
			  vx[1][ntotal + nvirt] = 0.;
			  vx[2][ntotal + nvirt] = 0.;
		}
		//Monaghan type virtual particle on the Right side
		for(int j=1;j<=2*mp-1;j++)
		{
			  nvirt = nvirt + 1;
			  x[1][ntotal + nvirt] = x1;
			  x[2][ntotal + nvirt] = j*dx/2;
			  vx[1][ntotal + nvirt] = 0.;
			  vx[2][ntotal + nvirt] = 0.;
		}
		for(int j=1;j<=nvirt;j++)
		{
			  rho[ntotal + j] = 1000.;
			  mass[ntotal + j] = rho[ntotal + i] * dx * dx;
			  p[ntotal + j] =0.;
			  u[ntotal + j] = 357.1;
			  itype[ntotal + j] = -2;
			  hsml[ntotal + j] = dx;
		}
	}
	if ((itimestep%save_step)==0)
	{
		in1 = fopen("xv_vp.dat", "w");
		in2 = fopen("state_vp.da", "w");
		in3 = fopen("other_vp.dat", "w");

		fprintf(in1, "%d\n", nvirt);
		for(int j=ntotal+1;j<=ntotal+nvirt;j++)
		{
			for (int d = 1; d <= dim; d++)
					fprintf(in1, "%d %f %f ", i, x[d][i], vx[d][i]);
			fprintf(in2, "%d %f %f %f", i, mass[i], rho[i], p[i], u[i]);
			fprintf(in3, "%d %d %f ", i, itype[i], hsml[i]);
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
