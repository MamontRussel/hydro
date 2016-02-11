#ifndef SINGLE_STEP_H
#define SINGLE_STEP_H

#include "virt_part.h"
#include "viscosity.h"
#include "art_visc.h"
#include "art_heat.h"
#include "av_vel.h"
#include "density.h"
#include "ext_force.h"
#include "int_force.h"
#include "hsml.h"
#include "direct_find.h"

// Function to determine the right hand side of a differential
// equation in a single step for performing time integration
// In this routine and its subroutines the SPH algorithms are performed
// itimestep: Current timestep number [in]
// dt Timestep [in]
// ntotal Number of particles [in]
// hsml Smoothing Length [in]
// mass Particle masses [in]
// x  Particle position [in]
// vx Particle velocity [in]
// u Particle internal energy [in]
// s Particle entropy (not used here) [in]
// rho Density [in/out]
// P Pressure [out]
// t Temperature [in/out]
// tdsdt Production of viscous entropy t*ds/dt [out]
// dx dx = vx = dx/dt [out]
// dvx dvx = dvx/dt, force per unit mass [out]
// du du = du/dt [out]
// ds ds = ds/dt [out]
// drho drho = drho/dt [out]
// itype Type of particle [in]
// av Monaghan average velocity [out]

void single_step(int &itimestep, float &dt, int &ntotal, float *hsml, float *mass,
  float **x,float **vx,float *u,float *rho,float *p,float *tdsdt,
  float **dvx,float *du,float *drho,int *itype,float **av);

#endif
