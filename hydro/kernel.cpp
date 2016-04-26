#include "kernel.h"

void kernel(float r,float* dx,float hsml,float &w,float* dwdx)
{
    float q,factor;
    q = r/hsml;
    w = 0.;
    for (int d = 1; d <= dim; d++)
        dwdx[d] = 0.;

    if (skf==1)
    {
        if (dim==1)factor = 1./hsml;
        else if (dim==2)factor = 15./(7.*pi*hsml*hsml);
        else if (dim==3)factor = 3./(2.*pi*hsml*hsml*hsml);
        else
        {
            cout<<">>> Error <<< : Wrong dimension: Dim ="<<dim<<endl;
            return;
        }

        if (q>=0&&q<=1)
        {
            w = factor * (2./3. - q*q + powf(q,3)/2.);
            for(int d=1;d<=dim;d++)
                dwdx[d] = factor * (-2.+3./2.*q)/powf(hsml,2) * dx[d];
        }
        else if (q>1&&q<=2)
        {
            w = factor * 1./6. * powf((2.-q),3);
            for(int d=1;d<=dim;d++)
                dwdx[d] =-factor * 1./6. * 3. *powf((2.-q),2)/hsml * (dx[d]/r);
        }
        else
        {
            w = 0.;
            for(int d=1;d<=dim;d++)
                dwdx[d] = 0.;
        }
    }
    else if (skf==2)
    {
        factor = 1./(powf(hsml,dim) * powf(pi,(dim/2.)));
        if(q>=0&&q<=3)
        {
            w = factor * exp(-q*q);
            for(int d=1;d<=dim;d++)
                dwdx[d] = w * ( -2.*dx[d]/hsml/hsml);
        }
        else
        {
            w = 0.;
            for(int d=1;d<=dim;d++)
                dwdx[d] = 0.;
        }
    }
    else if (skf==3)
    {
        if (dim==1)factor = 1./(120.*hsml);
        else if (dim==2)factor = 7./(478.*pi*hsml*hsml);
        else if (dim==3)factor = 1./(120.*pi*hsml*hsml*hsml);
        else
        {
            cout<<" >>> Error <<< : Wrong dimension: Dim ="<<dim<<endl;
            return;
        }
        if(q>=0&&q<=1)
        {
            w = factor * ( powf((3-q),5) - 6*powf((2-q),5) + 15*powf((1-q),5));
            for(int d=1;d<=dim;d++)
                dwdx[d] = factor * ( (-120 + 120*q - 50*q*q)/hsml*hsml * dx[d] );
        }
        else if(q>1&&q<=2)
        {
            w = factor * ( powf((3-q),5) - 6*powf((2-q),5));
            for(int d=1;d<=dim;d++)
                dwdx[d] = factor * (-5*powf((3-q),4) + 30*powf((2-q),4))/hsml * (dx[d]/r);
        }
        else if(q>2&&q<=3)
        {
            w = factor * powf((3-q),5);
            for(int d=1;d<=dim;d++)
                dwdx[d] = factor * (-5*powf((3-q),4)) / hsml * (dx[d]/r);
        }
        else
        {
            w = 0.;
            for(int d=1;d<=dim;d++)
                dwdx[d] = 0.;
        }
    }
}
