#ifndef EXT_FORCE_H
#define EXT_FORCE_H

#include <math.h>
#include "modul.h"

//  Function to calculate the external forces, e.g. gravitational forces
//  The forces from the interactions with boundary virtual particles
//  are also calculated here as external forces.
//  here as the external force.
//  ntotal : Number of particles [in]
//  mass : Particle masses [in]
//  x : Coordinates of all particles [in]
//  pair_i : List of first partner of interaction pair [in]
//  pair_j : List of second partner of interaction pair [in]
//  itype : type of particles [in]
//  hsml : Smoothing Length [in]
//  dvxdt : Acceleration with respect to x, y and z [out]

void ext_force(int ntotal,float **x,int niac,int *pair_i,int *pair_j,int *itype,float **dvxdt);

#endif
