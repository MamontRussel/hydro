#ifndef AV_VEL_H
#define AV_VEL_H

#include "modul.h"

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

void av_vel(int ntotal,float *mass,int &niac,int *pair_i,int *pair_j,float *w,float **vx,
  float *rho,float **av);

#endif
