#include "single_step.h"

void single_step(int &itimestep, float &dt, int &ntotal, float *hsml, float *mass,
  float **x,float **vx,float *u,float *rho,float *p,float *tdsdt,
  float **dvx,float *du,float *drho,int *itype,float **av)
{

    int nvirt, niac = 0;
    int *pair_i = new int[max_interaction];
    int *pair_j = new int[max_interaction];
    int *ns = new int[maxn];
    float *w = new float[max_interaction];
    float **dwdx = new float*[dim+1];
    float **indvxdt = new float*[dim+1];
    float **exdvxdt = new float*[dim+1];
    float **ardvxdt = new float*[dim+1];
    float *avdudt = new float[maxn];
    float *ahdudt = new float[maxn];
    float *c = new float[maxn];
    float *eta = new float[maxn];

    for (int i = 0; i <= dim; i++)
    {
        dwdx[i] = new float[max_interaction];
        indvxdt[i] = new float[maxn];
        exdvxdt[i] = new float[maxn];
        ardvxdt[i] = new float[maxn];
        for (int j = 1; j < maxn; j++)
        {
            indvxdt[i][j] = (float)NULL;
            exdvxdt[i][j] = (float)NULL;
            ardvxdt[i][j] = (float)NULL;;
        }
        for (int j = 0; j < max_interaction; j++)
            dwdx[i][j] = (float)NULL;

    }

    for(int i=1;i<=ntotal;i++)
    {
        avdudt[i] = 0.;
        ahdudt[i] = 0.;
        for(int d=1;d<=dim;d++)
        {
            indvxdt[d][i] = 0.;
            ardvxdt[d][i] = 0.;
            exdvxdt[d][i] = 0.;
        }
    }

    //Positions of virtual (boundary) particles:
    nvirt = 0;
    if (virtual_part)virt_part(itimestep, ntotal,nvirt,hsml,mass,x,vx,rho,u,p,itype);

    //Interaction parameters, calculating neighboring particles
    //and optimzing smoothing length
    if (nnps==1)direct_find(itimestep, ntotal + nvirt, hsml,x,niac, pair_i,pair_j,w,dwdx,ns);

    //Density approximation or change rate
    if (summation_density)sum_density(ntotal+nvirt,hsml,mass,niac,pair_i,pair_j,w,rho);
    else con_density(ntotal+nvirt,mass,niac,pair_i,pair_j,dwdx,vx,drho);

    //Dynamic viscosity:
    if (visc)viscosity(ntotal+nvirt,itype,eta);

    //Internal forces:
    int_force(ntotal+nvirt,mass,vx,niac,rho,eta,pair_i,pair_j,dwdx,
        u,itype,c,p,indvxdt,tdsdt,du);

    //Artificial viscosity:
    if (visc_artificial)art_visc(ntotal+nvirt,hsml,mass,x,vx,niac,rho,c,pair_i,pair_j,dwdx,ardvxdt, avdudt);

    //External forces:
    if (ex_force)ext_force(ntotal+nvirt,x,niac,pair_i,pair_j,itype,exdvxdt);

    //Calculating the neighboring particles and undating HSML
    if (sle!=0)h_upgrade(dt,ntotal, mass, vx, rho, niac,pair_i, pair_j, dwdx, hsml);
    if (heat_artificial)art_heat(ntotal+nvirt,hsml,mass,x,vx,niac,rho,u, c,pair_i,pair_j,dwdx,ahdudt);

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
