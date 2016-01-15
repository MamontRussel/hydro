// Subroutine to evolve smoothing length
// dt time step[in]
// ntotal Number of particles[in]
// mass Particle masses [in]
// vx Velocities of all particles [in]
// rho Density [in]
// niac Number of interaction pairs  [in]
// pair_i List of first partner of interaction pair  [in]
// pair_j List of second partner of interaction pair [in]
// dwdx Derivative of kernel with respect to x, y and [in]
// hsml Smoothing Length [in/out]

void h_upgrade(double &dt,int ntotal,double *mass,double **vx,double *rho,int &niac,int *pair_i,
	int *pair_j,double **dwdx,double *hsml)
{
	int i, j;
	double fac, hvcc;
	double *dvx = new double[dim];
	double *vcc = new double[maxn];
	double *dhsml = new double[maxn];

	if (sle==0 )//Keep smoothing length unchanged.
		return;
	else if (sle==2)
	{
		//	dh/dt = (-1/dim)*(h/rho)*(drho/dt).
		for(i=1;i<=ntotal;i++)
			vcc[i] = 0;
		for(int k=1;k<=niac;k++)
		{
			i = pair_i[k];
			j = pair_j[k];
			for(int d=0;d<dim;d++)
				dvx[d] = vx[d][j] - vx[d][i];
			hvcc = dvx[0]*dwdx[0][k];
			for(int d=1;d<dim;d++)
				hvcc = hvcc + dvx[d]*dwdx[d][k];

			vcc[i] = vcc[i] + mass[j] *hvcc/rho[j];
			vcc[j] = vcc[j] + mass[i] *hvcc/rho[i];
		}
		for(i=1;i<=ntotal;i++)
		{
			dhsml[i] = (hsml[i]/dim)*vcc[i];
			hsml[i] = hsml[i] + dt*dhsml[i];
			if (hsml[i]<=0) hsml[i] = hsml[i] - dt*dhsml[i];
		}
	}
	else if(sle==1)
	{
		fac = 2.0;
		for(i=1;i<=ntotal;i++)
			hsml[i] = fac * pow(mass[i]/rho[i],1/dim);
	}

	delete dvx;
	delete vcc;
	delete dhsml;
}
