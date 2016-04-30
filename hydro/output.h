#ifndef OUTPUT_H
#define OUTPUT_H

#include <cstdio>
#include <QCoreApplication>
#include "modul.h"

// Function for saving particle information to external disk file
// x-- coordinates of particles
// vx-- velocities of particles
// mass-- mass of particles
// rho-- densities of particles
// p-- pressure of particles
// u-- internal energy of particles
// c-- sound velocity of particles
// itype-- types of particles
// hsml-- smoothing lengths of particles
// ntotal-- total particle number

void output(float **x,float **vx,float *mass,float *rho,float *p,float *u,int *itype,float *hsml,int ntotal);

#endif
