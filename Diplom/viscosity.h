// Subroutine to define the fluid particle viscosity
// ntotal : Number of particles
// itype : Type of particle
// x : Coordinates of all particles
// rho : Density
// eta : Dynamic viscosity

void viscosity(int ntotal,int *itype,double **x,double *rho,double *eta)
{
	//double precision 15 знаков учитывается при счете
	for (int i = 1; i <= ntotal; i++)
	{
		if(abs(itype[i]==1))eta[i]=0.;
		else if (abs(itype[i])==2)eta[i]=1.0e-3;
	}
}
