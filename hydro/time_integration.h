#ifndef TIME_INTEGRATION_H
#define TIME_INTEGRATION_H

#include "single_step.h"
#include "output.h"

using namespace std;

// x-- coordinates of particles[ input/output]
// vx-- velocities of particles [input/output]
// mass-- mass of particles [input]
// rho-- dnesities of particles [input/output]
// p-- pressure of particles [input/output]
// u-- internal energy of particles
// c-- sound velocity of particles
// s-- entropy of particles, not used here
// e-- total energy of particles
// itype-- types of particles
// 	=1 ideal gas
// 	=2 water
// 	=3 tnt
// hsml-- smoothing lengths of particles
// ntotal-- total particle number
// maxtimestep-- maximum timesteps
// dt-- timestep

void time_integration(float **x, float **vx, float *mass,
  float *rho,float *p,float *u,int *itype,float *hsml,int &ntotal,int &maxtimestep,float &dt );

#endif
