#ifndef DIRECT_FIND_H
#define DIRECT_FIND_H

#include <math.h>
#include <iostream>
#include "kernel.h"
#include "modul.h"

using namespace std;

// Function to calculate the smoothing funciton for each particle and
// the interaction parameters used by the SPH algorithm. Interaction
// pairs are determined by directly comparing the particle distance
// with the corresponding smoothing length.
// itimestep Current time step [in]
// ntotal Number of particles [in]
// hsml Smoothing Length [in]
// x Coordinates of all particles [in]
// niac Number of interaction pairs [out]
// pair_i List of first partner of interaction pair [out]
// pair_j List of second partner of interaction pair [out]
// w Kernel for all interaction pairs [out]
// dwdx Derivative of kernel with respect to x, y and z [out]
// countiac Number of neighboring particles [out]

void direct_find(int itimestep, int ntotal, float *hsml, float **x,int &niac,int *pair_i,
    int *pair_j,float *w,float **dwdx,int *countiac);

#endif
