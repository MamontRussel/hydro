#include "EOS.h"

void p_gas(float rho, float u,float  &p,float &c)
{
    float gamma = 1.4;
    // For air (idea gas)
    p = (gamma-1) * rho * u;
    c = sqrt((gamma-1) * u);

}

void p_art_water(float rho,float &p,float &c)
{
//    float gamma, rho0,b;
//    //  Artificial EOS, Form 1 (Monaghan, 1994)
//    gamma = 7;
//    rho0 = 1000;
//    b = 1.013e5;
//    p = b*(powf(rho/rho0,gamma)-1);
//    c = 1480;

    //  Artificial EOS, Form 2 (Morris, 1997)
    c = 0.01;
    p = c * c * rho;
}
