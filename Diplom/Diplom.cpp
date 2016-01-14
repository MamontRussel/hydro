// ==============================================
// This is a three dimensional SPH code, the followings are the
// basic parameters needed in this code or calculated by this code
// mass-- mass of particles [in]
// ntotal-- total particle number ues [in]
// dt Time step used in the time integration [in]
// itype-- types of particles [in]
// x-- coordinates of particles [in/out]
// vx-- velocities of particles [in/out]
// rho-- dnesities of particles [in/out]
// p-- pressure of particles [in/out]
// u-- internal energy of particles [in/out]
// hsml-- smoothing lengths of particles [in/out]
// c-- sound velocity of particles [out]
// s-- entropy of particles [out]
// e-- total energy of particles [out]

#include "stdafx.h"
#include <iostream>
#include "modul.h"
#include "time_print.h"
#include "input.h"
#include "output.h"
#include "time_integration.h"
#include <ctime>

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	int *itype = new int[maxn];
	int ntotal, maxtimestep,yesorno=1;
	double **x = new double*[dim];
	double **vx = new double*[dim];
	double *mass = new double[maxn];
	double *rho = new double[maxn];
	double *p = new double[maxn];
	double *u = new double[maxn];
	double *c = new double[maxn];
	double *s = new double[maxn];
	double *e = new double[maxn];
	double *hsml = new double[maxn];
	double dt;

	for (int i = 0; i < dim; i++)
	{
		x[i] = new double[maxn];
		vx[i] = new double[maxn];
		for (int j = 0; j < maxn; j++)
		{
			x[i][j] = NULL;
			vx[i][j] = NULL;
			mass[j] = NULL;
			rho[j] = NULL;
			p[j] = NULL;
			u[j] = NULL;
			c[j] = NULL;
			s[j] = NULL;
			hsml[j] = NULL;
		}
	}

	//Time_Print();
	clock_t begin = clock();

	if (shocktube) dt = 0.005;
	if (shearcavity) dt = 5.e-5;

	input(x, vx, mass, rho, p, u, itype, hsml, ntotal);
	while (true)
	{
		cout << "***************************************************\n";
		cout << "Please input the maximal time steps\n";
		cout << "****************************************************\n";
		cin >> maxtimestep;

		time_integration(x, vx, mass, rho, p, u, c, s, e, itype, hsml, ntotal, maxtimestep, dt);
		//==========================================================

		output(x, vx, mass, rho, p, u, c, itype, hsml, ntotal);

		cout << "****************************************************\n";
		cout << "Are you going to run more time steps ? (0=No, 1=yes)\n";
		cout << "****************************************************\n";
		cin >> yesorno;
		if (yesorno == 0)break;
	}
	//Time_Print();
	clock_t end = clock();
	
	//Время вычислений?
	cout << "Elapsed CPU time = " << double(end - begin) / CLOCKS_PER_SEC<<endl;

	delete itype;
	delete mass;
	delete rho;
	delete p;
	delete u;
	delete c;
	delete s;
	delete e;
	delete hsml;

	return 0;
}

