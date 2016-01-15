// 	Subroutine to calculate the artificial heat (Fulk, 1994, p, a-17)
// ntotal Number of particles [in]
// hsml Smoothing Length [in]
// mass Particle masses [in]
// x Coordinates of all particles [in]
// vx Velocities of all particles [in]
// rho Density [in]
// u specific internal energy [in]
// c Sound veolcity [in]
// niac Number of interaction pairs [in]
// pair_i List of first partner of interaction pair [in]
// pair_j List of second partner of interaction pair [in]
// w Kernel for all interaction pairs [in]
// dwdx Derivative of kernel with respect to x, y and z [in]
// dedt produced artificial heat, adding to energy Eq. [out]

void art_heat(int ntotal,double *hsml,double *mass,double **x,double **vx,int &niac,double *rho,
	double *u,double *c,int *pair_i,int *pair_j,double* w,double **dwdx,double *dedt)
{
	int i,j;
	double dx, rr, h, mrho, mhsml, hvcc, mui, muj, muij, rdwdx, g1, g2;
	double *dvx = new double[dim];
	double *vcc = new double[maxn];

	//Parameter for the artificial heat conduction:
	g1=0.1;
	g2=1.0;
	for(i=1;i<=ntotal;i++)
	{
  		vcc[i] = 0;
  		dedt[i] = 0;
	}
	for(int k=1;k<=niac;k++)
	{
  		i = pair_i[k];
  		j = pair_j[k];
		for(int d=0;d<dim;d++)
  			dvx[d] = vx[d][j] - vx[d][i];

  		hvcc = dvx[0]*dwdx[0][k];
  		for(int d=1;d<dim;d++)
		  hvcc = hvcc + dvx[d]*dwdx[d][k];

  		vcc[i] = vcc[i] + mass[j]*hvcc/rho[j];
  		vcc[j] = vcc[j] + mass[i]*hvcc/rho[i];
	}
	for(int k=1;k<=niac;k++)
	{
  		i = pair_i[k];
  		j = pair_j[k];
  		mhsml= (hsml[i]+hsml[j])/2;
  		mrho = 0.5e0*(rho[i] + rho[j]);
  		rr = 0;
  		rdwdx = 0;
  		for(int d=0;d<dim;d++)
		{
  			dx = x[d][i] - x[d][j];
  			rr = rr + dx*dx;
  			rdwdx = rdwdx + dx*dwdx[d][k];
  		}
  		// mui=g1*hsml[i]*c[i] + g2*hsml[i]**2*(abs(vcc[i])-vcc[i]);
  		// muj=g1*hsml[j]*c[j] + g2*hsml[j]**2*(abs(vcc[j])-vcc[j]);
		mui=g1*hsml[i]*c[i] + g2*hsml[i]*hsml[i]*(abs(vcc[i])-vcc[i]);//Посчитаем что это операция возведения в степень
  		muj=g1*hsml[j]*c[j] + g2*hsml[j]*hsml[j]*(abs(vcc[j])-vcc[j]);//только hsml. Т.к сама операция выполняется раньше других
  		muij = 0.5*(mui+muj);
  		//h = muij/(mrho*(rr+0.01*mhsml**2))*rdwdx;
		h = muij/(mrho*(rr+0.01*mhsml*mhsml))*rdwdx;
  		dedt[i] = dedt[i] + mass[j]*h*(u[i]-u[j]);
  		dedt[j] = dedt[j] + mass[i]*h*(u[j]-u[i]);
	}
	for(i=1;i<=ntotal;i++)
  		dedt[i] = 2*dedt[i];

	delete dvx;
	delete vcc;
}
