#ifndef VIRT_PART_H
#define VIRT_PART_H

#include "modul.h"
#include <iostream>

using namespace std;

// Function to determine the information of virtual particles
// Here only the Monaghan type virtual particles for the 2D shear cavity driven problem are generated.
// itimestep : Current time step [in]
// ntotal Number of particles
// nvirt Number of virtual particles
// hsml Smoothing Length
// mass Particle masses
// x Coordinates of all particles
// vx Velocities of all particles
// rho Density
// u internal energy
// itype: type of particles

void virt_part(int &itimestep,int ntotal,int &nvirt,float *hsml,float *mass,float **x,
    float **vx,float *rho,float *u,float *p,int *itype);

#endif
