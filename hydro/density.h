#ifndef DENSITY_H
#define DENSITY_H

#include "kernel.h"

// 	Function to calculate the density with SPH summation algorithm.
// 	ntotal Number of particles [in]
// 	hsml Smoothing Length [in]
// 	mass Particle masses [in]
// 	niac Number of interaction pairs [in]
// 	pair_i List of first partner of interaction pair [in]
// 	pair_j List of second partner of interaction pair [in]
// 	w Kernel for all interaction pairs [in]
//  x : Coordinates of all particles [in]
//  rho : Density [out]

void sum_density(int ntotal, float *hsml, float *mass, int niac, int *pair_i,int *pair_j,
    float *w,float *rho);

// Function to calculate 'the density with SPH continuiity approach.
// ntotal Number of particles [in]
// mass Particle masses [in]
// niac Number of interaction pairs [in]
// pair_i List of first partner of interaction pair [in]
// pair_j List of second partner of interaction pair [in]
// dwdx derivation of Kernel for all interaction pairs [in]
// vx Velocities of all particles [in]
// x Coordinates of all particles [in]
// rho Density [in]
// drhodt Density change rate of each particle [out]

void con_density(int ntotal,float *mass,int niac,int *pair_i,int *pair_j,float **dwdx,
    float **vx,float *drhodt);

#endif
