#ifndef EOS_H
#define EOS_H

#include <math.h>

// 	Gamma law EOS: function to calculate the pressure and sound
// 	rho : Density [in]
//  u : Internal energy [in]
// 	p : Pressure [out]
// 	c : sound velocity [out]

void p_gas(float rho, float u,float  &p,float &c);

//  Artificial equation of s.tate for the artificial compressibility
//  rho : Density [in]
//  u : Internal energy [in]
//  p : Pressure [out]
//  c : sound velocity [out]
//  Equation of state for artificial compressibility

void p_art_water(float rho,float &p,float &c);

#endif
