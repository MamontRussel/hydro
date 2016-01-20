// Function for saving particle information to external disk file
// x-- coordinates of particles
// vx-- velocities of particles
// mass-- mass of particles
// rho-- densities of particles
// p-- pressure of particles
// u-- internal energy of particles
// c-- sound velocity of particles
// itype-- types of particles
// hsml-- smoothing lengths of particles
// ntotal-- total particle number

#include <cstdio>

void output(double **x,double **vx,double *mass,double *rho,
  double *p,double *u,double *c,int *itype,double *hsml,int ntotal)
{
	FILE *out1, *out2, *out3;

	out1 = fopen("f_xv.dat", "w");
	out2 = fopen("f_state.dat", "w");
	out3 = fopen("f_other.dat", "w");

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
