#ifndef KERNEL_H
#define KERNEL_H

#include <math.h>
#include <iostream>
#include "modul.h"

using namespace std;

// Function to calculate the smoothing kernel wij and its
// derivatives dwdxij.
// if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 198S)
// 	= 2, Gauss kernel (Gingold and Monaghan 1981)
// 	= 3, Quintic kernel (Morris 1997)
// r : Distance between particles i and j [in]
// dx : x-, y- and z-distance between i and j [in]
// hsml : Smoothing length [in]
// w : Kernel for all interaction pairs [out]
// dwdx : Derivative of kernel with respect to x, y and z [out]

void kernel(float r,float* dx,float hsml,float &w,float* dwdx);

#endif
