// 	Gamma law EOS: function to calculate the pressure and sound
// 	rho : Density [in]
//  u : Internal energy [in]
// 	p : Pressure [out]
// 	c : sound velocity [out]

void p_gas(double rho, double u,double  &p,double &c)
{
	double gamma = 1.4;
	// For air (idea gas)
	p = (gamma-1) * rho * u;
	//cout << p <<" "<<(gamma - 1)<<" "<<rho<<" "<<u<<endl;
	c = sqrt((gamma-1) * u);
}

//  Artificial equation of s.tate for the artificial compressibility
//  rho : Density [in]
//  u : Internal energy [in]
//  p : Pressure [out]
//  c : sound velocity [out]
//  Equation of state for artificial compressibility

void p_art_water(double rho,double &p,double &c)
{
	double gamma, rho0;
	//  Artificial EOS, Form 1 (Monaghan, 1994)
	gamma = 7;
	rho0 = 1000;
	// b = 1.013e5
	//  p = b*((rho/rho0)**gamma-l)
	//  c = 148 0.
	//  Artificial EOS, Form 2 (Morris, 1997)
	c = 0.01;
	p = c*c * rho;
}
