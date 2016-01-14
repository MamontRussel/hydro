// ================================================
//  Subroutine to calculate the external forces, e.g. gravitational forces
//  The forces from the interactions with boundary virtual particles
//  are also calculated here as external forces.
//  here as the external force.
//  ntotal : Number of particles [in]
//  mass : Particle masses [in]
//  x : Coordinates of all particles [in]
//  pair_i : List of first partner of interaction pair [in]
//  pair_j : List of second partner of interaction pair [in]
//  itype : type of particles [in]
//  hsml : Smoothing Length [in]
//  dvxdt : Acceleration with respect to x, y and z [out]

void ext_force(int ntotal,double *mass,double **x,int niac,int *pair_i,int *pair_j,int *itype,
	double *hsml,double **dvxdt)
{
	int i, j;
	double *dx = new double[dim];
	double rr, f, rr0, dd, p1, p2;

	for(int i=1;i<=ntotal;i++)
		for(int d=0;d<dim;d++)
			dvxdt[d][i] = 0;

	// Consider self-gravity or not ?
	if (self_gravity)
		for(i=1;i<=ntotal;i++)
			dvxdt[dim-1][i] = -9.8;

	// Boundary particle force and penalty anti-penetration force.
	rr0 = 1.25e-5;
	dd = 1.e-2;
	p1 = 12;
	p2 = 4;
	for(int k=1;k<=niac;k++)
	{
		i = pair_i[k];
		j = pair_j[k];
		if(itype[i]>0&&itype[j]<0)
		{
			rr = 0;
			for(int d=0;d<dim;d++)
			{
				dx[d] = x[d][i] - x[d][j];
				rr = rr + dx[d]*dx[d];
			}
			rr = sqrt(rr);
			if(rr<rr0)
			{
				//непонятно с порядком действий
				//f = ((rr0/rr)**p1-(rr0/rr)**p2)/rr**2
				f = (pow(rr0/rr,p1)-(pow(rr0/rr,p2)))/rr*rr;
				for(int d=0;d<dim;d++)
					dvxdt[d][i] = dvxdt[d][i] + dd*dx[d]*f;
			}
		}
	}

	delete dx;
}
