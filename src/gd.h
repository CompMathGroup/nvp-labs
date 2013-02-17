#ifndef __GD_H__
#define __GD_H__

#include <cstddef>
#include <cmath>
#include <algorithm>
#include "matrix.h"

#define RHO_MIN 1e-9
#define PRESS_MIN 1e-9

struct AbsoluteValue {
	double operator()(double s) const {	return std::fabs(s); }
};

struct Vars {
	double rho;
	double P;
	double E;
	double GAMMA;

	Vars() {}

	Vars(const Vector &o, double _GAMMA) {
		rho = o(0);
		P = rho * o(1);
		E = rho * (.5 * o(1) * o(1) + o(2));
		GAMMA = _GAMMA;
	}

	Vars(double _rho, double _P, double _E, double _GAMMA = 1.4) : 
		rho(_rho), P(_P), E(_E), GAMMA(_GAMMA) { }

	Vars operator +(const Vars o) const {
		return Vars(rho + o.rho, P + o.P, E + o.E, GAMMA);
	}

	Vars operator -(const Vars o) const {
		return Vars(rho - o.rho, P - o.P, E - o.E, GAMMA);
	}

	friend Vars operator *(const double m, const Vars o) {
		return Vars(m * o.rho, m * o.P, m * o.E, o.GAMMA);
	}

	Vars &operator *=(const double m) {
		rho *= m;
		P *= m;
		E *= m;
		return *this;
	}

	Vars &operator +=(const Vars o) {
		rho += o.rho;
		P += o.P;
		E += o.E;
		return *this;
	}

	Vars &operator -=(const Vars o) {
		rho -= o.rho;
		P -= o.P;
		E -= o.E;
		return *this;
	}

	void fixup() {
		if (rho < RHO_MIN) {
			rho = RHO_MIN;
			P = 0; /* WILL fail otherwise */
		}

		double v = P / rho;
		double K = .5 * P * v;
		double p =  (GAMMA - 1) * (E - K);

		if (p < PRESS_MIN) {
			p = PRESS_MIN;
			E = K + p / (GAMMA - 1);
		}
	}

	Vars f() const {
		double u = P / rho;
		double K = .5 * P * u;
		double p = (GAMMA - 1) * (E - K);
		
		return Vars(P, p + P * u, (E + p) * u);
	}
	double lmax(const Vector &lambda) const {
		double u = lambda(0);
		double c = u - lambda(1);
		return c + std::fabs(u);
	}
	double lmax() const {
		double u = P / rho;
		double K = .5 * P * u;
		double p = (GAMMA - 1) * (E - K);
		double c = sqrt(GAMMA * p / rho);
		return c + std::fabs(u);
	}
	void eigen(Matrix &W, Vector &lambda, Matrix &iW) const {
		double u = P / rho;
		double K = .5 * P * u;
		double p = (GAMMA - 1) * (E - K);
		double c = sqrt(GAMMA * p / rho);

		double u2 = u * u;
		double c2 = c * c;
		double zm = (GAMMA - 1) / (c2 + (GAMMA - 1) * u * (0.5 * u - c));
		double zp = (GAMMA - 1) / (c2 + (GAMMA - 1) * u * (0.5 * u + c));

		W(0, 0) = 1;		W(0, 1) = zm;			W(0, 2) = zp;
		W(1, 0) = u;		W(1, 1) = zm * (u - c);	W(1, 2) = zp * (u + c);
		W(2, 0) = .5 * u2;	W(2, 1) = 1;			W(2, 2) = 1;

		double zz = 0.5 / (c * (zp + zm - u2 * zm * zp));
		double _2c = 2 * c;
		double cu = c * u;

		iW(0, 0) = zz * (2 * (u * (zp - zm) + c * (zm + zp)));
		iW(0, 1) = zz * (2 * (zm - zp));
		iW(0, 2) = zz * (-4 * c * zm * zp);

		iW(1, 0) = zz * (( 2 - cu * zp - u2 * zp) * u);
		iW(1, 1) = zz * (-2 + u2 * zp);
		iW(1, 2) = zz * (_2c * zp);

		iW(2, 0) = zz * ((-2 - cu * zm + u2 * zm) * u);
		iW(2, 1) = zz * ( 2 - u2 * zm);
		iW(2, 2) = zz * (_2c * zm);

		lambda(0) = u;
		lambda(1) = (u - c);
		lambda(2) = (u + c);
	}
	template <class Functor>
	Vars lmul(const Matrix &W, const Vector &e, const Matrix &iW, const Functor &func = Functor()) const {
		double x = 0, y = 0, z = 0;
		for (int j = 0; j < 3; j++) {
			double s = func(e(j)) * (iW(j, 0) * rho + iW(j, 1) * P + iW(j, 2) * E);
			x += W(0, j) * s;
			y += W(1, j) * s;
			z += W(2, j) * s;
		}
		return Vars(x, y, z);
	}
	Vector to_nconserv() const {
		Vector v;
		v(0) = rho;
		v(1) = P / rho;
		v(2) = (E - .5 * P * v(1)) / rho;
		return v;
	}
	Vector to_vect() const {
		return Vector(rho, P, E);
	}
};

#endif
