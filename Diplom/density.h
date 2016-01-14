// =================================================
// 	Subroutine to calculate the density with SPH summation algorithm.
// 	ntotal Number of particles [in]
// 	hsml Smoothing Length [in]
// 	mass Particle masses [in]
// 	niac Number of interaction pairs [in]
// 	pair_i List of first partner of interaction pair [in]
// 	pair_j List of second partner of interaction pair [in]
// 	w Kernel for all interaction pairs [in]
//  itype : type of particles [in]
//  x : Coordinates of all particles [in]
//  rho : Density [out]

#include "kernel.h"

void sum_density(int ntotal, double *hsml, double *mass, int niac, int *pair_i,int *pair_j, 
	double *w, int *itype, double *rho)
{
	int i, j;
	double selfdens,r;
	double *hv = new double[dim];
	double *wi = new double[maxn];

	// wi(maxn) integration'of the kernel itself
	for(int d=0;d<dim;d++)
		hv[d] = 0;

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

// ===============================
// Subroutine to calculate 'the density with SPH continuiity approach.
// ntotal Number of particles [in]
// mass Particle masses [in]
// niac Number of interaction pairs [in]
// pair_i List of first partner of interaction pair [in]
// pair_j List of second partner of interaction pair [in]
// dwdx derivation of Kernel for all interaction pairs [in]
// vx Velocities of all particles [in]
// itype : type of particles [in]
// x Coordinates of all particles [in]
// rho Density [in]
// drhodt Density change rate of each particle [out]

void con_density(int ntotal,double *mass,int niac,int *pair_i,int *pair_j,double **dwdx,
	double **vx, int *itype,double **x,double *rho,double *drhodt)
{
	int i,j;
	double vcc;
	double *dvx = new double[dim];
	for(i=1;i<=ntotal;i++)
		drhodt[i] = 0;

	for(int k=1;k<=niac;k++)
	{
		i = pair_i[k];
		j = pair_j[k];
		for(int d=0;d<dim;d++)
		  dvx[d] = vx[d][i] - vx[d][j];

		vcc = dvx[0]*dwdx[0][k];
		for(int d=0;d<dim;d++)
		  vcc = vcc + dvx[d]*dwdx[d][k];

		drhodt[i] = drhodt[i] + mass[j]*vcc;
		drhodt[j] = drhodt[j] + mass[i]*vcc;
	}

	delete dvx;
}
