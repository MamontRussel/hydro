#ifndef ART_VISC_H
#define ART_VISC_H

#include "modul.h"

// Function to calculate the artificial viscosity
// ntotal Number of particles (including virtual particles) [in]
// hsml Smoothing Length [in]
// mass Particle masses [in]
// x Coordinates of all particles [in]
// vx Velocities of all particles [in]
// niac Number of interaction pairs [in]
// rho Density [in]
// c Temperature [in]
// pair_i List of first partner of interaction pair [in]
// pair_j List of second partner of interaction pair [in]
// w Kernel for all interaction pairs [in]
// dwdx Derivative of kernel with respect to x, y and z [in]
// dvxdt Acceleration with respect to x, y and z [out]
// dedt Change of specific internal energy [out]

void art_visc(int ntotal,float *hsml,float *mass,float **x,float **vx,int niac,float *rho,
    float *c,int *pair_i,int *pair_j, float **dwdx,float **dvxdt,float *dedt);

#endif
