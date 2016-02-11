#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include "modul.h"

// Function for loading or generating initial particle information
// x-- coordinates of particles [out]
// vx-- velocities of particles [out]
// mass-- mass of particles [out]
// rho-- dnesities of particles [out]
// p-- pressure of particles [out]
// u-- internal energy of particles [out]
// itype-- types of particles [out]
// hsml-- smoothing lengths of particles [out]
// ntotal-- total particle number [out]
void input(float **x, float **vx, float *mass, float *rho,
    float *p, float *u, int *itype, float *hsml, int &ntotal);

// This Function is used to generate initial data for the 1d noh shock tube problem
// x-- coordinates of particles
// vx-- velocities of particles
// mass-- mass of particles
// rho-- dnesities of particles
// p-- pressure of particles
// u-- internal energy of particles
// itype-- types of particles
// =1 ideal gas
// hsml-- smoothing lengths of particles
// ntotal-- total particle number
void shock_tube(float **x, float **vx, float *mass, float *rho,
    float *p, float *u, int *itype, float *hsml, int &ntotal);

// This Function is used to generate
// 2d shear driven cavity probem with
// x-- coordinates of particles
// vx-- velocities of particles
// mass-- mass of particles
// rho-- dnesities of particles
// p-- pressure of particles
// u-- internal energy of particles
// itype-- types of particles
// 	=2 water
// h-- smoothing lengths of particles
// ntotal-- total particle number
// initial data for the Re = 1
void shear_cavity(float **x, float **vx, float *mass, float *rho,
    float *p, float *u, int *itype, float *hsml, int &ntotal);

#endif
