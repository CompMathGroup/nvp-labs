#include "riemann.h"
#include <math.h>

Riemann::Riemann(double _gamma) {
	gamma = _gamma;
	chi = 0.5*(gamma - 1)/gamma;
	eta = 0.5*(gamma + 1)/gamma;
	tau = gamma * (gamma - 1);
}

double Riemann::f(double P, double p, double rho) {
	double pi = P/p;
	double c = sqrt(gamma * p / rho);
	if (P >= p) 
		return (P-p)/(rho * c * sqrt(eta * pi + chi));
	return 2./(gamma - 1) * c * (pow(pi, chi) - 1);
}

double Riemann::f_P(double P, double p, double rho) {
	double pi = P/p;
	double c = sqrt(gamma * p / rho);
	if (P >= p) {
		double b3 = eta * pi + chi;
		b3 *= sqrt(b3);
		return ((gamma + 1) * pi + 3 * gamma - 1) / (4 * rho * c * gamma * b3);
	}
	return c / (gamma * P) * pow(pi, chi);
}

Conf Riemann::solve(double _r1, double _u1, double _e1, double _r2, double _u2, double _e2) {
	r1 = _r1;
	u1 = _u1;
	e1 = _e1;
	r2 = _r2;
	u2 = _u2;
	e2 = _e2;
	p1 = r1 * e1 * (gamma - 1),
	p2 = r2 * e2 * (gamma - 1);
	c1 = sqrt(tau * e1),
	c2 = sqrt(tau * e2);

	Conf ret;
	double U0 = f(0, p1, r1) + f(0, p2, r2);
 
	if (U0 >= u1 - u2) {
		P = 0;
		R1 = R2 = 0;
		E1 = E2 = 0; 
		W2 = u2 + c2;
		W2s = u2 - 2 * c2 / (gamma - 1);
		W1 = u1 - c1;
		W1s = u1 + 2 * c1/ (gamma - 1);
		U = 0.5 * (W1s + W2s);
		return DW_VA_DW;
	}
	
	double dP;
	
	iters = 0;
	P = 0.5 * (p1 + p2);
	while (f(P, p1, r1) + f(P, p2, r2) - u1 + u2 > 0)
		P *= 0.5; /* Dirty way to find initial P satisfying F(P) - u1 + u2 < 0*/
	while (1) { /* Newton iterations */
		iters ++;
		dP = (f(P, p1, r1) + f(P, p2, r2) - u1 + u2) / (f_P(P, p1, r1) + f_P(P, p2, r2));
		P -= dP;
		if ((fabs(dP) < 1e-10 * fabs(P)) || (iters > 10))
			break;
	}

	U = 0.5 * (u1 + u2 + f(P, p2, r2) - f(P, p1, r1));

	if (P > p1) {
		double m1 = sqrt(0.5 * r1 * (P * (gamma + 1) + p1 * (gamma - 1)));
		R1 = r1 * m1 / (m1 - (u1 - U) * r1);
		E1 = P / (R1 * (gamma - 1));
		W1s = W1 = u1 - m1 / r1;
		ret = (Conf)(SW_CS_DW | SW_CS_SW);
	} else {
		double c1s = c1 + 0.5 * (gamma - 1) * (u1 - U);
		R1 = gamma * P / (c1s * c1s);
		E1 = P / (R1 * (gamma - 1));
		W1 = u1 - c1;
		W1s = U - c1s;
		ret = (Conf)(DW_CS_DW | DW_CS_SW);
	}

	if (P > p2) {
		double m2 = sqrt(0.5 * r2 * (P * (gamma + 1) + p2 * (gamma - 1)));
		R2 = r2 * m2 / (m2 + (u2 - U) * r2);
		E2 = P / (R2 * (gamma - 1));
		W2s = W2 = u2 + m2 / r2;
		ret = (Conf)(ret & (SW_CS_SW | DW_CS_SW));
	} else {
		double c2s = c2 - 0.5 * (gamma - 1) * (u2 - U);
		R2 = gamma * P / (c2s * c2s);
		E2 = P / (R2 * (gamma - 1));
		W2 = u2 + c2;
		W2s = U + c2s;
		ret = (Conf)(ret & (SW_CS_DW | DW_CS_DW));
	}
	return ret;
}

void Riemann::evaluate(double t, int np, double *x, double *r, double *u, double *e) {
	double x0 = .5;
	for (int i = 0; i < np; i++) {
		double xi = (x[i] - x0) / t;
		if (xi < W1) {
			r[i] = r1;
			u[i] = u1;
			e[i] = e1;
			continue;
		}
		if (xi < W1s) {
			double cs = 2/(gamma + 1) * c1 + (gamma - 1)/(gamma + 1) * (u1 - xi);
			double ps;
			r[i] = r1 * pow(cs / c1, 2/(gamma - 1));
			u[i] = xi + cs;
			ps = p1 * pow(cs / c1, 2*gamma / (gamma - 1));
			e[i] = ps / (r[i] * (gamma - 1));
			continue;
		}
		if (xi < U) {
			r[i] = R1;
			u[i] = U;
			e[i] = E1;
			continue;
		}
		if (xi < W2s) {
			r[i] = R2;
			u[i] = U;
			e[i] = E2;
			continue;
		}
		if (xi < W2) {
			double cs = 2/(gamma + 1) * c2 - (gamma - 1)/(gamma + 1) * (u2 - xi);
			double ps;
			r[i] = r2 * pow(cs / c2, 2 / (gamma - 1));
			u[i] = xi - cs;
			ps = p2 * pow(cs / c2, 2 * gamma / (gamma - 1));
			e[i] = ps / (r[i] * (gamma - 1));
			continue;
		}
		r[i] = r2;
		u[i] = u2;
		e[i] = e2;
	}
}

void Riemann::twoCS(double t, int np, double *x, double *r, double *u, double *e) {
	double v = 0.5 * (u1 + u2); /* should be v = u1 = u2*/
	double x1 = 0.4 + v * t;
	double x2 = 0.6 + v * t;
	x1 = x1 - floor(x1);
	x2 = x2 - floor(x2);

	/* also p1 = p2 <=> r1 e1 = r2 e2 */
	if (x1 < x2)
		for (int i=0; i<np; i++) {
			r[i] = (x[i] > x1 && x[i] <= x2) ? r1 : r2;
			e[i] = (x[i] > x1 && x[i] <= x2) ? e1 : e2;
		}
	else
		for (int i=0; i<np; i++) {
			r[i] = (x[i] > x2 && x[i] < x1) ? r2 : r1;
			e[i] = (x[i] > x2 && x[i] < x1) ? e2 : e1;
		}
	for (int i=0; i<np; i++) 
		u[i] = v;
}
