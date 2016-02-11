#ifndef VISCOSITY_H
#define VISCOSITY_H

#include <math.h>
#include <cstdlib>

// Function to define the fluid particle viscosity
// ntotal : Number of particles
// itype : Type of particle
// x : Coordinates of all particles
// rho : Density
// eta : Dynamic viscosity

void viscosity(int ntotal,int *itype,float *eta);

#endif
