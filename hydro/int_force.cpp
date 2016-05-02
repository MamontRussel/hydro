#include "int_force.h"

void int_force(int ntotal, float *mass,
    float **vx, int niac, float *rho, float *eta, int *pair_i, int *pair_j, float **dwdx,
    float *u, int *itype, float *c, float *p, float **dvxdt, float *tdsdt,
    float *dedt)
{
    int i, j;
    float *dvx = new float[dim+1];
    float *txx = new float[maxn];
    float *tyy = new float[maxn];
    float *tzz = new float[maxn];
    float *txy = new float[maxn];
    float *txz = new float[maxn];
    float *tyz = new float[maxn];
    float *vcc = new float[maxn];
    float hxx, hyy, hzz, hxy, hxz, hyz, h, hvcc, he, rhoij;

    // Initialization of shear tensor, velocity divergence,
    // 	viscous energy, internal energy, acceleration
    for (i = 1; i <= ntotal; i++)
    {
        txx[i] = 0.e0;
        tyy[i] = 0.e0;
        tzz[i] = 0.e0;
        txy[i] = 0.e0;
        txz[i] = 0.e0;
        tyz[i] = 0.e0;
        vcc[i] = 0.e0;
        tdsdt[i] = 0.e0;
        dedt[i] = 0.e0;
        for (int d = 1; d <= dim; d++)
            dvxdt[d][i] = 0.e0;
    }

    // Calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab
    if (visc)
    {
        for (int k = 1; k <= niac; k++)
        {
            i = pair_i[k];
            j = pair_j[k];
            for (int d = 1; d <= dim; d++)
                dvx[d] = vx[d][j] - vx[d][i];
            //-------------------------------------
            if (dim == 1) hxx = 2 * dvx[1] * dwdx[1][k];
            else  if (dim == 2)
            {
                hxx = 2.e0 * dvx[1] * dwdx[1][k] - dvx[2] * dwdx[2][k];
                hxy = dvx[1] * dwdx[2][k] + dvx[2] * dwdx[1][k];
                hyy = 2.e0*dvx[2] * dwdx[2][k] - dvx[1] * dwdx[1][k];
            }
            else if (dim == 3)
            {
                hxx = 2.e0 * dvx[1] * dwdx[1][k] - dvx[2] * dwdx[2][k] - dvx[3] * dwdx[3][k];
                hxy = dvx[1] * dwdx[2][k] + dvx[2] * dwdx[1][k];
                hxz = dvx[1] * dwdx[3][k] + dvx[3] * dwdx[1][k];

                hyy = 2.e0 * dvx[2] * dwdx[2][k] - dvx[1] * dwdx[1][k] - dvx[3] * dwdx[3][k];
                hyz = dvx[2] * dwdx[3][k] + dvx[3] * dwdx[2][k];
                hzz = 2.e0 * dvx[3] * dwdx[3][k] - dvx[1] * dwdx[1][k] - dvx[2] * dwdx[2][k];
            }
            if (dim == 1)
            {
                hxx = 2.e0 / 3. * hxx;
                txx[i] = txx[i] + mass[j] * hxx / rho[j];
                txx[j] = txx[j] + mass[i] * hxx / rho[i];
            }
            else if (dim == 2)
            {
                hxx = 2.e0 / 3. * hxx;
                hyy = 2.e0 / 3. * hyy;
                txx[i] = txx[i] + mass[j] * hxx / rho[j];
                txx[j] = txx[j] + mass[i] * hxx / rho[i];
                txy[i] = txy[i] + mass[j] * hxy / rho[j];
                txy[j] = txy[j] + mass[i] * hxy / rho[i];
                tyy[i] = tyy[i] + mass[j] * hyy / rho[j];
                tyy[j] = tyy[j] + mass[i] * hyy / rho[i];
            }
            else if (dim == 3)
            {
                hxx = 2.e0 / 3. * hxx;
                hyy = 2.e0 / 3. * hyy;
                hzz = 2.e0 / 3. * hzz;
                txx[i] = txx[i] + mass[j] * hxx / rho[j];
                txx[j] = txx[j] + mass[i] * hxx / rho[i];
                txy[i] = txy[i] + mass[j] * hxy / rho[j];
                txy[j] = txy[j] + mass[i] * hxy / rho[i];
                txz[i] = txz[i] + mass[j] * hxz / rho[j];
                txz[j] = txz[j] + mass[i] * hxz / rho[i];
                tyy[i] = tyy[i] + mass[j] * hyy / rho[j];
                tyy[j] = tyy[j] + mass[i] * hyy / rho[i];
                tyz[i] = tyz[i] + mass[j] * hyz / rho[j];
                tyz[j] = tyz[j] + mass[i] * hyz / rho[i];
                tzz[i] = tzz[i] + mass[j] * hzz / rho[j];
                tzz[j] = tzz[j] + mass[i] * hzz / rho[i];
            }

            //	Calculate SPH sum for vc,c = dvx/dx + dvy/dy + dvz/dz:
            hvcc = 0.;
            for (int d = 1; d <= dim; d++)
                hvcc = hvcc + dvx[d] * dwdx[d][k];
            vcc[i] = vcc[i] + mass[j] * hvcc / rho[j];
            vcc[j] = vcc[j] + mass[i] * hvcc / rho[i];
        }
    }

    for (i = 1; i <= ntotal; i++)
    {
        // Viscous entropy Tds/dt = 1/2 eta/rho Tab Tab
        if (visc)
        {
            if (dim == 1)tdsdt[i] = txx[i] * txx[i];
            else if (dim == 2)
                tdsdt[i] = txx[i] * txx[i] + 2.e0 * txy[i] * txy[i] + tyy[i] * tyy[i];
            else if (dim == 3)
                tdsdt[i] = txx[i] * txx[i] + 2.e0 * txy[i] * txy[i] + 2. * txz[i] * txz[i]
                + tyy[i] * tyy[i] + 2.e0*tyz[i] * tyz[i] + tzz[i] * tzz[i];
            tdsdt[i] = 0.5e0*eta[i] / rho[i] * tdsdt[i];
        }

        // Pressure from equation of state
        if (abs(itype[i]) == 1)p_gas(rho[i], u[i], p[i], c[i]);
        else if (abs(itype[i]) == 2)p_art_water(rho[i], p[i], c[i]);
    }

    // Calculate SPH sum for pressure force -p,a/rho
    // and viscous force (eta Tab),b/rho
    // and the internal energy change de/dt due to -p/rho vc,c
    for (int k = 1; k <= niac; k++)
    {
        i = pair_i[k];
        j = pair_j[k];
        he = 0.;

        // For SPH algorithm 1
        rhoij = 1.e0 / (rho[i] * rho[j]);
        if (pa_sph == 1)
        {
            for (int d = 1; d <= dim; d++)
            {
                //Pressure part
                h = -(p[i] + p[j])*dwdx[d][k];
                he = he + (vx[d][j] - vx[d][i])*h;
                // Viscous force
                if (visc)
                {
                    if (d == 1)
                    {
                        //	x-coordinate of acceleration
                        h = h + (eta[i] * txx[i] + eta[j] * txx[j])*dwdx[1][k];
                        if (dim >= 2)
                        {
                            h = h + (eta[i] * txy[i] + eta[j] * txy[j])*dwdx[2][k];
                            if (dim == 3)
                                h = h + (eta[i] * txz[i] + eta[j] * txz[j])*dwdx[3][k];
                        }
                    }
                    else if (d == 2)
                    {
                        //	y-coordinate of acceleration
                        h = h + (eta[i] * txy[i] + eta[j] * txy[j])*dwdx[1][k] + (eta[i] * tyy[i] + eta[j] * tyy[j])*
                            dwdx[2][k];
                        if (dim == 3)
                            h = h + (eta[i] * tyz[i] + eta[j] * tyz[j])*dwdx[3][k];
                    }
                    else if (d == 3)
                    {
                        // z-coordinate of acceleration
                        h = h + (eta[i] * txz[i] + eta[j] * txz[j])*dwdx[1][k] + (eta[i] * tyz[i] + eta[j] * tyz[j])
                            *dwdx[2][k] + (eta[i] * tzz[i] + eta[j] * tzz[j])*dwdx[3][k];
                    }
                }
                h = h*rhoij;
                dvxdt[d][i] = dvxdt[d][i] + mass[j] * h;
                dvxdt[d][j] = dvxdt[d][j] - mass[i] * h;
            }
            he = he*rhoij;
            dedt[i] = dedt[i] + mass[j] * he;
            dedt[j] = dedt[j] + mass[i] * he;
            //For SPH algorithm 2 USED HERE
        }
        else if (pa_sph == 2)
        {
            for (int d = 1; d <= dim; d++)
            {
                h = -(p[i] / pow(rho[i], 2) + p[j] / pow(rho[j], 2))*dwdx[d][k];
                he = he + (vx[d][j] - vx[d][i])*h;
                //	Viscous force
                if (visc)
                {
                    if (d == 1)
                    {
                        // x-coordinate of acceleration
                        h = h + (eta[i] * txx[i] / pow(rho[i], 2) + eta[j] * txx[j] / pow(rho[j], 2))*dwdx[1][k];
                        if (dim >= 2)
                        {
                            h = h + (eta[i] * txy[i] / pow(rho[i], 2) + eta[j] * txy[j] / pow(rho[j], 2))*dwdx[2][k];
                            if (dim == 3)
                                h = h + (eta[i] * txz[i] / pow(rho[i], 2) + eta[j] * txz[j] / pow(rho[j], 2))*dwdx[3][k];
                        }
                    }
                    else if (d == 2)
                    {
                        //	y-coordinate of acceleration
                        h = h + (eta[i] * txy[i] / pow(rho[i], 2) + eta[j] * txy[j] / pow(rho[j], 2))*dwdx[1][k]
                            + (eta[i] * tyy[i] / pow(rho[i], 2) + eta[j] * tyy[j] / pow(rho[j], 2))*dwdx[2][k];
                        if (dim == 3)
                            h = h + (eta[i] * tyz[i] / pow(rho[i], 2) + eta[j] * tyz[j] / pow(rho[j], 2))*dwdx[3][k];
                    }
                    else if (d == 3)
                    {
                        // z-coordinate of acceleration
                        h = h + (eta[i] * txz[i] / pow(rho[i], 2) + eta[j] * txz[j] / pow(rho[j], 2))*dwdx[1][k]
                            + (eta[i] * tyz[i] / pow(rho[i], 2) + eta[j] * tyz[j] / pow(rho[j], 2))*dwdx[2][k]
                            + (eta[i] * tzz[i] / pow(rho[i], 2) + eta[j] * tzz[j] / pow(rho[j], 2))*dwdx[3][k];
                    }
                }
                dvxdt[d][i] = dvxdt[d][i] + mass[j] * h;
                dvxdt[d][j] = dvxdt[d][j] - mass[i] * h;
            }
            dedt[i] = dedt[i] + mass[j] * he;
            dedt[j] = dedt[j] + mass[i] * he;
        }
    }
    // Change of specific internal energy de/dt = T ds/dt - p/rho vc,c:

    for (i = 1; i <= ntotal; i++)
        dedt[i] = tdsdt[i] + 0.5e0*dedt[i];

    delete[] dvx;
    delete[] txx;
    delete[] tyy;
    delete[] tzz;
    delete[] txy;
    delete[] txz;
    delete[] tyz;
    delete[] vcc;
}
