#ifndef ART_HEAT_H
#define ART_HEAT_H

#include "modul.h"
#include <math.h>

// Function to calculate the artificial heat
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

void art_heat(int ntotal,float *hsml,float *mass,float **x,float **vx,int niac,float *rho,
    float *u,float *c,int *pair_i,int *pair_j,float **dwdx,float *dedt);

#endif
