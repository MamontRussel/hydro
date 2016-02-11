#ifndef HSML_H
#define HSML_H

#include <math.h>
#include "modul.h"

// Function to evolve smoothing length
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

void h_upgrade(float &dt,int ntotal,float *mass,float **vx,float *rho,int &niac,int *pair_i,
    int *pair_j,float **dwdx,float *hsml);

#endif
