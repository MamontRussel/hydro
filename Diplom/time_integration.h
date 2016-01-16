// x-- coordinates of particles[ input/output]
// vx-- velocities of particles [input/output]
// mass-- mass of particles [input]
// rho-- dnesities of particles [input/output]
// p-- pressure of particles [input/output]
// u-- internal energy of particles
// c-- sound velocity of particles
// s-- entropy of particles, not used here
// e-- total energy of particles
// itype-- types of particles
// 	=1 ideal gas
// 	=2 water
// 	=3 tnt
// hsml-- smoothing lengths of particles
// ntotal-- total particle number
// maxtimestep-- maximum timesteps
// dt-- timestep

#include "single_step.h"

void time_integration(double **x, double **vx, double *mass,
  double *rho,double *p,double *u,double *c,double *s,
  double *e,int *itype,double *hsml,int &ntotal,int &maxtimestep,double &dt )
{

  int itimestep, current_ts=0, nstart=0;
  double **v_min = new double*[dim+1];
  double *u_min = new double[maxn];
  double *rho_min = new double[maxn];
  double **dx = new double*[dim+1];
  double **dvx = new double*[dim+1];
  double *du = new double[maxn];
  double *drho = new double[maxn];
  double **av = new double*[dim+1];
  double *ds = new double[maxn];
  double *t = new double[maxn];
  double *tdsdt = new double[maxn];
  double time = 0, temp_rho, temp_u;

  for (int i = 0; i <= dim; i++)
  {
	  v_min[i] = new double[maxn];
	  dx[i] = new double[maxn];
	  dvx[i] = new double[maxn];
	  av[i] = new double[maxn];
	  for (int j = 0; j < maxn; j++)
	  {
		  v_min[i][j] = NULL;
		  dx[i][j] = NULL;
		  dvx[i][j] = NULL;
		  av[i][j] = NULL;;
	  }
  }

  for(int i=1;i<=ntotal;i++)
    for(int d=1;d<=dim;d++)
    {
      av[d][i] = 0;
    }

  for (itimestep = nstart+1;itimestep<=nstart+maxtimestep;itimestep++)
  {
    current_ts=current_ts+1;
    if ((itimestep%print_step)==0)
    {
      cout<<"-------------------------------------------------\n";
      cout<<"current time step = "<<itimestep<<" current time= "<<(double)(time+dt)<<endl;
      cout<<"-------------------------------------------------\n";
    }

    //If not first time step, then update thermal energy, density and
    //velocity half a time step
    if (itimestep!=1)
    {
      for(int i=1;i<=ntotal;i++)
      {
        u_min[i] = u[i];
        temp_u=0;
        if (dim==1)temp_u= -nsym*p[i]*vx[1][i]/x[1][i]/rho[i];
        u[i]= u[i] + (dt/2.)* (du[i]+temp_u);
        if(u[i]<0) u[i] = 0;
        if (!summation_density)
        {
          rho_min[i] = rho[i];
          temp_rho=0;
          if (dim==1)temp_rho= -nsym*rho[i]*vx[1][i]/x[1][i];
          rho[i] = rho[i] +(dt/2.)*(drho[i]+ temp_rho);
        }
        for(int d=1;d<=dim;d++)
        {
          v_min[d][i] = vx[d][i];
          vx[d][i] = vx[d][i] + (dt/2.)*dvx[d][i];
        }
      }
    }

    // Definition of variables out of the function vector:
	single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u,
		s, rho, p, t, tdsdt, dx, dvx, du, ds, drho, itype, av);
    // ====================================================
    if (itimestep==1)
    {
      for(int i=1;i<=ntotal;i++)
      {
        temp_u=0;
        if (dim==1) temp_u=-nsym*p[i]*vx[1][i]/x[1][i]/rho[i];
        u[i] = u[i] + (dt/2.)*(du[i] + temp_u);
        if(u[i]<0) u[i] = 0;
        if(!summation_density )
        {
          temp_rho=0;
          if (dim==1) temp_rho=-nsym*rho[i]*vx[1][i]/x[1][i];
          rho[i] = rho[i] + (dt/2.)* (drho[i]+temp_rho);
        }
        for(int d=1;d<=dim;d++)
        {
          vx[d][i] = vx[d][i] + (dt/2.) * dvx[d][i] + av[d][i];
          x[d][i] = x[d][i] + dt * vx[d][i];
        }
      }
    }
    else
    {
      for(int i=1;i<=ntotal;i++)
      {
        temp_u=0;
        if (dim==1) temp_u=-nsym*p[i]*vx[1][i]/x[1][i]/rho[i];
        u[i] = u_min[i] + dt*(du[i]+temp_u);
        if(u[i]<0) u[i] = 0;

        if (!summation_density )
        {
          temp_rho=0;
          if (dim==1) temp_rho=-nsym*rho[i]*vx[1][i]/x[1][i];
          rho[i] = rho_min[i] + dt*(drho[i]+temp_rho);
        }

        for(int d=1;d<=dim;d++)
        {
          vx[d][i] = v_min[d][i] + dt * dvx[d][i] + av[d][i];
          x[d][i] = x[d][i] + dt * vx[d][i];
        }
      }
    }
    time = time + dt;

    if ((itimestep % save_step)==0)
    {
		output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal);
    }
    if ((itimestep & print_step)==0)
    {
		cout << "x  velocity  dvx\n";
		cout << x[1][moni_particle] << " " << vx[1][moni_particle] << " " << dvx[1][moni_particle]<<endl;
    }
  }

  //double **v_min = new double*[dim];
  for (int i = 0; i <= dim; i++)
	  delete v_min[i];
  delete v_min;
  delete u_min;
  delete rho_min;
  //double **dx = new double*[dim];
  //double **dvx = new double*[dim];
  delete du;
  delete drho;
  //double **av = new double*[dim];
  delete ds;
  delete t;
  delete tdsdt;

  nstart=current_ts;
}
