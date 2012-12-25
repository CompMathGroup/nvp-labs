#ifndef __RIEMANN_H__
#define __RIEMANN_H__

/* Simplified ideal gas case */
/* Mostly taken from Kulikovskii/Pogorelov/Semenov */

enum Conf {
	DW_CS_DW = 1,
	SW_CS_DW = 2,
	DW_CS_SW = 4,
	SW_CS_SW = 8,
	DW_VA_DW = 16
};

class Riemann {
	double gamma;
	double chi, eta, tau;

	double f(double P, double p, double rho);
	double f_P(double P, double p, double rho);

	double r1, u1, e1, p1, r2, u2, e2, p2;	
	double c1, c2;
	double U, P;
	double W1, W1s, W2, W2s;
	double R1, E1, R2, E2;
	int iters;
public:	
	Riemann(double _gamma);
	Conf solve(double r1, double u1, double e1, double r2, double u2, double e2);
	void evaluate(double t, int np, double *x, double *r, double *u, double *e);
	void twoCS(double t, int np, double *x, double *r, double *u, double *e);
};

#endif
