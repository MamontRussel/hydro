// Function for loading or generating initial particle information
// x-- coordinates of particles [out]
// vx-- velocities of particles [out]
// mass-- mass of particles [out]
// rho-- dnesities of particles [out]
// p-- pressure of particles [out]
// u-- internal energy of particles [out]
// itype-- types of particles [out]
// hsml-- smoothing lengths of particles [out]
// ntotal-- total particle number [out]

#include <cstdio>

void shock_tube(double **x, double **vx, double *mass, double *rho,
	double *p, double *u, int *itype, double *hsml, int &ntotal);

void shear_cavity(double **x, double **vx, double *mass, double *rho,
	double *p, double *u, int *itype, double *hsml, int &ntotal);

void input(double **x, double **vx, double *mass, double *rho,
	double *p, double *u, int *itype, double *hsml, int &ntotal)
{

	FILE *in1, *in2, *in3;
	int im;
	// load initial particle information from external disk file
	if(config_input)  // FALSE!!!
	{
		in1 = fopen("f_xv.dat", "r");
		in2 = fopen("f_state.dat", "r");
		in3 = fopen("f_other.dat", "r");

		cout << "**********************************************\n";
		cout << "Loading initial particle configuration...\n";
		fscanf(in1, "%d", &ntotal);
		cout << "Total number of particles " << ntotal << endl;
		cout << "**********************************************\n";
		for(int i=1;i<=ntotal;i++)
		{
			for (int d = 1; d <= dim; d++)
				fscanf(in1, "%d %f %f ", &im, &x[d][i], &vx[d][i]);
			fscanf(in2, "%d %f %f %f %f",&im, &mass[i], rho[i],p[i],u[i]);
			fscanf(in3, "%d %d %f",&im, &itype[i], hsml[i]);
		}
	}
	else
	{   
		in1 = fopen("ini_xv.dat", "w");
		in2 = fopen("ini_state.dat", "w");
		in3 = fopen("ini_other.dat", "w");
		// shocktube = TRUE
		if (shocktube)shock_tube(x, vx, mass, rho, p, u, itype, hsml, ntotal);
		if (shearcavity)shear_cavity(x, vx, mass, rho, p, u, itype, hsml, ntotal);

		for(int i=1;i<=ntotal;i++)
		{
			for (int d = 1; d <= dim; d++)
				fprintf(in1, "%d %f %f \n", i, x[d][i], vx[d][i]);
			fprintf(in2, "%d %f %f %f %f \n", i, mass[i], rho[i], p[i], u[i]);
			fprintf(in3, "%d %d %f \n", i, itype[i], hsml[i]);
		}
		cout << "**********************************************\n";
		cout << "**Initial particle configuration generated ***\n";
		cout << "********  Total number of particles   " << ntotal << endl;
		cout << "**********************************************\n";
	}
	fclose(in1);
	fclose(in2);
	fclose(in3);
}

  // This Function is used to generate initial data for the 1d noh shock tube problem
  // x-- coordinates of particles
  // vx-- velocities of particles
  // mass-- mass of particles
  // rho-- dnesities of particles
  // p-- pressure of particles
  // u-- internal energy of particles
  // itype-- types of particles
  // =1 ideal gas
  // hsml-- smoothing lengths of particles
  // ntotal-- total particle number
void shock_tube(double **x, double **vx, double *mass, double *rho,
	double *p, double *u, int *itype, double *hsml, int &ntotal)
 {
  
	double space_x;
	ntotal=400;
	space_x=0.6/80.;

	for(int i=1;i<=ntotal;i++)
	{
		mass[i]=0.75/400.;
		hsml[i]=0.015;
		itype[i]=1;
		for(int d=1;d<=dim;d++)
		{
			x[d][i] = 0.;
			vx[d][i] = 0.;

		}
	}

	for(int i=1;i<=320;i++)
		x[1][i]=-0.6+space_x/4.*(i-1);

	for(int i=320+1;i<=ntotal;i++)
		x[1][i]=0.+space_x*(i-320);

	for(int i=1;i<=ntotal;i++)
	{
		if (x[1][i]<=1.e-8)
		{
			u[i]=2.5;
			rho[i]=1.;
			p[i]=1.;
		}
		if (x[1][i]>1.e-8)
		{
			u[i]=1.795;
			rho[i]=0.25;
			p[i]=0.1795;
		}
	}
}

  // This Function is used to generate
  // 2d shear driven cavity probem with
  // x-- coordinates of particles
  // vx-- velocities of particles
  // mass-- mass of particles
  // rho-- dnesities of particles
  // p-- pressure of particles
  // u-- internal energy of particles
  // itype-- types of particles
  // 	=2 water
  // h-- smoothing lengths of particles
  // ntotal-- total particle number
  // initial data for the Re = 1
void shear_cavity(double **x, double **vx, double *mass, double *rho,
	double *p, double *u, int *itype, double *hsml, int &ntotal)
{
	int  m, n, mp, np, k;
	double x1, y1, dx, dy;
	// Giving mass and smoothing length as well as other data,
	m = 41;
	n = 41;
	mp = m-1;
	np = n-1;
	ntotal = mp * np;
	x1 = 1.e-3;
	y1 = 1.e-3;
	dx = x1/mp;
	dy = y1/np;

	for(int i=1;i<=mp;i++)
		for(int j=1;j<=np;j++)
		{
			  k = j + (i-1)* np;
			  x[1][k] = (i-1)*dx + dx/2.;
			  x[2][k] = (j-1)*dy + dy/2.;
		}

	for(int i=1;i<=mp*np;i++)
	{
		vx[1][i] = 0.;
		vx[2][i] = 0.;
		rho[i] = 1000.;
		mass[i] = dx*dy*rho[i];
		p[i]= 0.;
		u[i]=357.1;
		itype[i] = 2;
		hsml[i] = dx;
	}
}
