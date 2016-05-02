#include "direct_find.h"
#include <QDebug>

void direct_find(int itimestep, int ntotal, float *hsml, float **x,int &niac,int *pair_i,
    int *pair_j,float *w,float **dwdx,int *countiac)
{
    int i, j,sumiac, maxiac, miniac, noiac, maxp=0, minp, scale_k;
    float *dxiac = new float[dim+1];
    float driac, r, mhsml;
    float *tdwdx = new float[dim+1];

    if (skf==1) scale_k = 2;
    else if (skf == 2) scale_k = 3;
    else if (skf == 3) scale_k = 3;

    for(i=1;i<=ntotal;i++)
        countiac[i] = 0;

    niac = 0;

    for(i=1;i<=ntotal-1;i++)
        for(j=i+1;j<=ntotal;j++)
        {
            dxiac[1] = x[1][i] - x[1][j];
            driac = dxiac[1]*dxiac[1];
            for(int d=2;d<=dim;d++)
            {
                dxiac[d] = x[d][i] - x[d][j];
                driac = driac + dxiac[d]*dxiac[d];
            }
            mhsml = (hsml[i]+hsml[j])/2.;
            if (sqrt(driac)<scale_k*mhsml)
            {
                if (niac<max_interaction)
                {
                    // Neighboring pair list, and totalinteraction number and
                    // the interaction number for each particle
                    niac++;
                    pair_i[niac] = i;
                    pair_j[niac] = j;
                    r = sqrt(driac);
                    countiac[i]++;
                    countiac[j]++;
                    // Kernel and derivations of kernel
                    kernel(r,dxiac,mhsml,w[niac],tdwdx);
                    for(int d=1;d<=dim;d++)
                        dwdx[d][niac] = tdwdx[d];
                }
                else
                {
                    cout<<">>> ERROR: Too many interactions\n";
                    return;
                }
            }
        }

    // Statistics for the interaction
    sumiac = 0;
    maxiac = 0;
    miniac = 1000;
    noiac = 0;

    for(int i=1;i<=ntotal;i++)
    {
        sumiac = sumiac + countiac[i];
        if (countiac[i]>maxiac)
        {
            maxiac = countiac[i];
            maxp = i;
        }
        if (countiac[i]<miniac)
        {
            miniac = countiac[i];
            minp = i;
        }
        if (countiac[i]==0) noiac++;
    }

    if (itimestep % print_step==0)
    {
        if (int_stat)
        {
            cout<<"**** Statistics: interactions per particle:\n";
            cout<<"**** Particle:"<<maxp<<" maximal interactions:"<<maxiac<<endl;
            cout<<"**** Particle:"<<minp<<" minimal interactions:"<<miniac<<endl;
            cout<<"**** Average :"<<(float)(sumiac/ntotal)<<endl;
            cout<<"**** Total pairs : "<<niac<<endl;
            cout<<"**** Particles with no interactions: "<<noiac<<endl;
        }
    }

    delete[] dxiac;
    delete[] tdwdx;
}
