#ifndef INT_FORCE_H
#define INT_FORCE_H

#include <math.h>
#include <cstdlib>
#include "modul.h"
#include "EOS.h"

// Function to calculate the internal forces on the right hand side
// of the Navier-Stokes equations, i.e. the pressure gradient and the
// gradient of the viscous stress tensor, used by the time integration.
// Moreover the entropy production due to viscous dissipation, tds/dt,
// and the change of internal energy per mass, de/dt, are calculated.
// itimestep : Current timestep number [in]
// dt : Time step [in]
// ntotal : Number of particles [in]
// hsml : Smoothing Length [in]
// mass : Particle masses [in]
// u : Particle internal energy [in]
// x : Particle coordinates [in]
// vx : Velocities of all particles [in]
// niac : Number of interaction pairs [in]
// rho :  Density [in]
// eta : Dynamic viscosity [in]
// pair_i : List of first partner of interaction pair [in]
// pair_j : List of second partner of interaction pair [in]
// dwdx : Derivative of kernel with respect to x, y and z [in]
// itype : Type of particle (material types) [in]
// itype : Particle type. [in]
// t : Particle temperature [in/out]
// c : Particle sound speed [out]
// p : Particle pressure [out]
// dvxdt : Acceleration with respect to x, y and z [out]
// tdsdt : Production of viscous entropy [out]
// dedt : Change of specific internal energy [out]

void int_force(int ntotal, float *mass,
    float **vx, int &niac, float *rho, float *eta, int *pair_i, int *pair_j, float **dwdx,
    float *u, int *itype, float *c, float *p, float **dvxdt, float *tdsdt,
    float *dedt);

#endif
