#include "time_integration.h"

void time_integration(float **x, float **vx, float *mass,
  float *rho,float *p,float *u,int *itype,float *hsml,int &ntotal,int &maxtimestep,float &dt )
{

  int itimestep, current_ts=0, nstart=0;
  float **v_min = new float*[dim+1];
  float *u_min = new float[maxn];
  float *rho_min = new float[maxn];
  float **dx = new float*[dim+1];
  float **dvx = new float*[dim+1];
  float *du = new float[maxn];
  float *drho = new float[maxn];
  float **av = new float*[dim+1];
  float *ds = new float[maxn];
  float *t = new float[maxn];
  float *tdsdt = new float[maxn];
  float time = 0, temp_rho, temp_u;

  for (int i = 0; i <= dim; i++)
  {
      v_min[i] = new float[maxn];
      dx[i] = new float[maxn];
      dvx[i] = new float[maxn];
      av[i] = new float[maxn];
      for (int j = 0; j < maxn; j++)
      {
          v_min[i][j] = (float)NULL;
          dx[i][j] = (float)NULL;
          dvx[i][j] = (float)NULL;
          av[i][j] = (float)NULL;;
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
      cout<<"current time step = "<<itimestep<<" current time= "<<(float)(time+dt)<<endl;
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
    single_step(itimestep, dt, ntotal, hsml, mass, x, vx, u,rho, p, tdsdt, dvx, du, drho, itype, av);
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
        output(x, vx, mass, rho, p, u,itype, hsml, ntotal);
    }
    if ((itimestep & print_step)==0)
    {
        cout << "x  velocity  dvx\n";
        cout << x[1][moni_particle] << " " << vx[1][moni_particle] << " " << dvx[1][moni_particle]<<endl;
    }
  }

  //float **v_min = new float*[dim];
  for (int i = 0; i <= dim; i++)
      delete v_min[i];
  delete v_min;
  delete u_min;
  delete rho_min;
  //float **dx = new float*[dim];
  //float **dvx = new float*[dim];
  delete du;
  delete drho;
  //float **av = new float*[dim];
  delete ds;
  delete t;
  delete tdsdt;

  nstart=current_ts;
}
