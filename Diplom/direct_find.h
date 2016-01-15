// Subroutine to calculate the smoothing funciton for each particle and
// the interaction parameters used by the SPH algorithm. Interaction
// pairs are determined by directly comparing the particle distance
// with the corresponding smoothing length.
// itimestep Current time step [in]
// ntotal Number of particles [in]
// hsml Smoothing Length [in]
// x Coordinates of all particles [in]
// niac Number of interaction pairs [out]
// pair_i List of first partner of interaction pair [out]
// pair_j List of second partner of interaction pair [out]
// w Kernel for all interaction pairs [out]
// dwdx Derivative of kernel with respect to x, y and z [out]
// countiac Number of neighboring particles [out]

void direct_find(int itimestep, int ntotal, double *hsml, double **x,int &niac,int *pair_i,
	int *pair_j,double *w,double **dwdx,int *countiac)
{
	int i, j,sumiac, maxiac, miniac, noiac, maxp=0, minp, scale_k;
	double *dxiac = new double[dim]; 
	double driac, r, mhsml;
	double *tdwdx = new double[dim];
	if (skf==1) scale_k = 2;
	else if (skf == 2) scale_k = 3;
	else if (skf == 3) scale_k = 3;

	for(i=1;i<=ntotal;i++)
		countiac[i] = 0;

	niac = 0;
	for(i=1;i<ntotal;i++)
		for(j=i+1;j<=ntotal;j++)
		{
			dxiac[0] = x[0][i] - x[0][j];
			driac = dxiac[0]*dxiac[0];
			for(int d=0;d<dim;d++)
			{
				dxiac[d] = x[d][i] - x[d][j];
				driac = driac + dxiac[d]*dxiac[d];
			}
			mhsml = (hsml[i]+hsml[j])/2;
			if (sqrt(driac)<scale_k*mhsml)
			{
				if (niac<max_interaction)
				{
					// Neighboring pair list, and totalinteraction number and
					// the interaction number for each particle
					niac++;
					pair_i[niac] = i;
					pair_j[niac] = j;
					r = sqrt(driac);
					countiac[i]++;
					countiac[j]++;
					// Kernel and derivations of kernel
					kernel(r,dxiac,mhsml,w[niac],tdwdx);
					for(int d=0;d<dim;d++)
						dwdx[d][niac] = tdwdx[d];
				}
				else
				{
					cout<<">>> ERROR: Too many interactions\n";
					return;
				}
			}
		}
 
	// Statistics for the interaction
	sumiac = 0;
	maxiac = 0;
	miniac = 1000;
	noiac = 0;

	for(int i=1;i<=ntotal;i++)
	{
		sumiac = sumiac + countiac[i];
		if (countiac[i]>maxiac)
		{
			maxiac = countiac[i];
			maxp = i;
		}
		if (countiac[i]<miniac)
		{
			miniac = countiac[i];
			minp = i;
		}
		if (countiac[i]==0) noiac++;
		//cout << maxp<<endl;
	}

	if (itimestep % print_step==0)
	{
		if (int_stat)
		{
			cout<<"**** Statistics: interactions per particle:\n";
			cout<<"**** Particle:"<<maxp<<" maximal interactions:"<<maxiac<<endl;
			cout<<"**** Particle:"<<minp<<" minimal interactions:"<<miniac<<endl;
			cout<<"**** Average :"<<(double)(sumiac/ntotal)<<endl;
			cout<<"**** Total pairs : "<<niac<<endl;
			cout<<"**** Particles with no interactions: "<<noiac<<endl;
		}
	}

	delete dxiac;
	delete tdwdx;
}
