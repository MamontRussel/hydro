#include "viscosity.h"

void viscosity(int ntotal,int *itype,float *eta)
{
    //float precision 15 знаков учитывается при счете
    for (int i = 1; i <= ntotal; i++)
    {
        if(abs(itype[i])==1)eta[i]=0.;
        else if (abs(itype[i])==2)eta[i]=1.0e-3;
    }
}
