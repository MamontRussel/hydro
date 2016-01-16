// Function to determine the right hand side of a differential
// equation in a single step for performing time integration
// In this routine and its subroutines the SPH algorithms are performed
// itimestep: Current timestep number [in]
// dt Timestep [in]
// ntotal Number of particles [in]
// hsml Smoothing Length [in]
// mass Particle masses [in]
// x  Particle position [in]
// vx Particle velocity [in]
// u Particle internal energy [in]
// s Particle entropy (not used here) [in]
// rho Density [in/out]
// P Pressure [out]
// t Temperature [in/out]
// tdsdt Production of viscous entropy t*ds/dt [out]
// dx dx = vx = dx/dt [out]
// dvx dvx = dvx/dt, force per unit mass [out]
// du du = du/dt [out]
// ds ds = ds/dt [out]
// drho drho = drho/dt [out]
// itype Type of particle [in]
// av Monaghan average velocity [out]

#include "virt_part.h"
#include "viscosity.h"
#include "art_visc.h"
#include "art_heat.h"
#include "av_vel.h"
#include "density.h"
#include "ext_force.h"
#include "int_force.h"
#include "hsml.h"
#include "direct_find.h"

void single_step(int &itimestep, double &dt, int &ntotal, double *hsml, double *mass,
  double **x,double **vx,double *u,double *s,double *rho,double *p,double *t,double *tdsdt,
  double **dx,double **dvx,double *du,double *ds,double *drho,int *itype,double **av)
{

	int nvirt = 0, niac = 0;
	int *pair_i = new int[max_interaction];
	int *pair_j = new int[max_interaction];
	int *ns = new int[maxn];
	double *w = new double[max_interaction];
	double **dwdx = new double*[dim+1];
	double **indvxdt = new double*[dim+1];
	double **exdvxdt = new double*[dim+1];
	double **ardvxdt = new double*[dim+1];
	double *avdudt = new double[maxn];
	double *ahdudt = new double[maxn];
	double *c = new double[maxn];
	double *eta = new double[maxn];

	for (int i = 0; i <= dim; i++)
	{
		dwdx[i] = new double[max_interaction];
		indvxdt[i] = new double[maxn];
		exdvxdt[i] = new double[maxn];
		ardvxdt[i] = new double[maxn];
		for (int j = 1; j < maxn; j++)
		{
			indvxdt[i][j] = NULL;
			exdvxdt[i][j] = NULL;
			ardvxdt[i][j] = NULL;;
		}
		for (int j = 0; j < max_interaction; j++)
			dwdx[i][j] = NULL;

	}

	for(int i=1;i<=ntotal;i++)
	{
		avdudt[i] = 0.;
		ahdudt[i] = 0.;
		for(int d=1;d<=dim;d++)
		{
    		indvxdt[d][i] = 0;
    		ardvxdt[d][i] = 0;
    		exdvxdt[d][i] = 0;
		}
	}

	//Positions of virtual (boundary) particles:
	nvirt = 0;
	if (virtual_part)virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,rho,u,p,itype);

	//Interaction parameters, calculating neighboring particles
	//and optimzing smoothing length
	if (nnps==1)direct_find(itimestep, ntotal + nvirt, hsml,x,niac, pair_i,pair_j,w,dwdx,ns);
	// в данной версии не поддерживается
	//else if (nnps==2)link_list(itimestep, ntotal+nvirt,hsml[1],x,niac,pair_i,pair_j,w,dwdx,ns);
	//else if (nnps==3)tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,pair_j,w,dwdx,ns);

	//Density approximation or change rate
	if (summation_density)sum_density(ntotal+nvirt,hsml,mass,niac,pair_i,pair_j,w,itype,rho);
	else con_density(ntotal+nvirt,mass,niac,pair_i,pair_j,dwdx,vx, itype,x,rho, drho);

	//Dynamic viscosity:
	if (visc)viscosity(ntotal+nvirt,itype,x,rho,eta);

	//Internal forces:
	int_force(itimestep,dt,ntotal+nvirt,hsml,mass,vx,niac,rho,eta,pair_i,pair_j,dwdx,
		u,itype,x,t,c,p,indvxdt,tdsdt,du);

	//Artificial viscosity:
	if (visc_artificial)art_visc(ntotal+nvirt,hsml,mass,x,vx,niac,rho,c,pair_i,pair_j,
		w,dwdx,ardvxdt, avdudt);

	//External forces:
	if (ex_force)ext_force(ntotal+nvirt,mass,x,niac,pair_i,pair_j,itype, hsml, exdvxdt);

	//Calculating the neighboring particles and undating HSML
	if (sle!=0)h_upgrade(dt,ntotal, mass, vx, rho, niac,pair_i, pair_j, dwdx, hsml);
	if (heat_artificial)art_heat(ntotal+nvirt,hsml,mass,x,vx,niac,rho,u, c,pair_i,pair_j,
		w,dwdx,ahdudt);

	//Calculating average velocity of each partile for avoiding penetration
	if (average_velocity)av_vel(ntotal,mass,niac,pair_i,pair_j, w, vx, rho, av);

	//Convert velocity, force, and energy to f and dfdt
	for(int i=1;i<=ntotal;i++)
	{
		for(int d=1;d<=dim;d++)
  			dvx[d][i] = indvxdt[d][i] + exdvxdt[d][i] + ardvxdt[d][i];

  		du[i] = du[i] + avdudt[i] + ahdudt[i];
	}
	if ((itimestep % print_step)==0)
	{
		cout<<"\n**** Information for particle **** "<<moni_particle<<endl;
		cout << "internal a=" << indvxdt[1][moni_particle] << " artifical a=" << ardvxdt[1][moni_particle]<<
			"\nexternal a=" << exdvxdt[1][moni_particle]<< " total a="<<dvx[1][moni_particle]<<endl;
	}

	delete pair_i;
	delete pair_j;
	delete ns;
	delete w;
	for (int i = 0; i <= dim; i++)
		delete dwdx[i];
	delete dwdx;
	for (int i = 0; i <= dim; i++)
		delete indvxdt[i];
	delete indvxdt;
	for (int i = 0; i <= dim; i++)
		delete exdvxdt[i];
	delete exdvxdt;
	for (int i = 0; i <= dim; i++)
		delete ardvxdt[i];
	delete ardvxdt;
	delete avdudt;
	delete ahdudt;
	delete c;
	delete eta;
}
