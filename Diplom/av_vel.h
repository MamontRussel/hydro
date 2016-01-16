// Function to calculate the average velocity to correct velocity
// for preventing.penetration (monaghan, 1992)
// ntotal Number of particles [in]
// mass : Particle masses [in]
// niac : Number of interaction pairs [in]
// pair_i : List of first partner of interaction pair [in]
// pair_j : List of second partner of interaction pair [in]
// w : Kernel for all interaction pairs [in]
// vx : Velocity of each particle [in]
// rho : Density of each particle [in]
// av : Average velocityof each particle [out]

void av_vel(int ntotal,double *mass,int &niac,int *pair_i,int *pair_j,double *w,double **vx,
  double *rho,double **av)
{
	int i,j;
	double *dvx = new double[dim+1];
	double epsilon;

	// epsilon a small constants chosen by experence,
	// may lead to instability,
	// for example, for the 1 dimensional shock tube problem, the E <= 0.3
	epsilon = 0.3;

	for(i=1;i<=ntotal;i++)
		for(int d=1;d<=dim;d++)
  			av[d][i] = 0;

	for(int k=1;k<=niac;k++)
	{
  		i = pair_i[k];
  		j = pair_j[k];
  		for(int d=1;d<=dim;d++)
		{
  			dvx[d] = vx[d][i] - vx[d][j];
  			av[d][i] = av[d][i] - 2*mass[j]*dvx[d]/(rho[i]+rho[j])*w[k];
  			av[d][j] = av[d][j] + 2*mass[i]*dvx[d]/(rho[i]+rho[j])*w[k];
  		}
	}

	for(i=1;i<=ntotal;i++)
		for(int d=1;d<=dim;d++)
  			av[d][i] = epsilon * av[d][i];

	delete dvx;
}
