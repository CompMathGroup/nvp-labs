#ifndef __COEFF_H__
#define __COEFF_H__

#include <iostream>
#include <cmath>

struct Rational {
	int n, m;
	double a[16];
	double b[16];

	inline double poly(double x, const double *c, int k) const {
		double y = c[k];
		switch (k) {
			case 15: y = y * x + c[14];
			case 14: y = y * x + c[13];
			case 13: y = y * x + c[12];
			case 12: y = y * x + c[11];
			case 11: y = y * x + c[10];
			case 10: y = y * x + c[ 9];
			case  9: y = y * x + c[ 8];
			case  8: y = y * x + c[ 7];
			case  7: y = y * x + c[ 6];
			case  6: y = y * x + c[ 5];
			case  5: y = y * x + c[ 4];
			case  4: y = y * x + c[ 3];
			case  3: y = y * x + c[ 2];
			case  2: y = y * x + c[ 1];
			case  1: y = y * x + c[ 0];
			case  0: break;
		}
		return y;
	}

	Rational(const double *p, const double *q) {
		for (int i = 0; i < 16; i++) {
			a[i] = p[i];
			b[i] = q[i];
		}
		n = 15;
		while ((n > 0) && (std::fabs(a[n]) > 1e-14))
			n--;
		m = 15;
		while ((m > 0) && (std::fabs(b[m]) > 1e-14))
			m--;
	}

	inline double operator()(double sigma) const {
		double p = poly(sigma, a, n);
		double q = poly(sigma, b, m);
		return p / q;
	}
};

struct Alphas {
	const Rational b1, b2, b3;
	const Rational b4, b5, b6, b7;
	const Rational g1, g2, g3;
	Alphas(
		const Rational &b1,
		const Rational &b2,
		const Rational &b3,

		const Rational &b4,
		const Rational &b5,
		const Rational &b6,
		const Rational &b7,

		const Rational &g1,
		const Rational &g2,
		const Rational &g3
	) : 
		b1(b1), b2(b2), b3(b3),
		b4(b4), b5(b5), b6(b6), b7(b7),
		g1(g1), g2(g2), g3(g3)
	{}
	bool sanity() {
		/* No checks for this type */
		return true;
	}
};

struct Gamma1 : public Alphas {
	Gamma1( const Alphas &o) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return g1(+sigma);
		else
			return g3(-sigma);
	}
};

struct Gamma2 : public Alphas {
	Gamma2( const Alphas &o) : Alphas(o) {}
	double operator()(double sigma) const {
		double s = std::fabs(sigma);
		return g2(s);
	}
};

struct Gamma3 : public Alphas {
	Gamma3( const Alphas &o) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return g3(+sigma);
		else
			return g1(-sigma);
	}
};

struct Beta1 : public Alphas {
	Beta1( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return b1(+sigma);
		else
			return b3(-sigma);
	}
};

struct Beta3 : public Alphas {
	Beta3( const Alphas &o) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return b3(+sigma);
		else
			return b1(-sigma);
	}
};

struct Beta2 : public Alphas {
	Beta2( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		double s = std::fabs(sigma);
		return b2(s);
	}
};

struct Beta4 : public Alphas {
	Beta4( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return +b4(+sigma);
		else
			return -b7(-sigma);
	}
};

struct Beta7 : public Alphas {
	Beta7( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return +b7(+sigma);
		else
			return -b4(-sigma);
	}
};

struct Beta5 : public Alphas {
	Beta5( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return +b5(+sigma);
		else
			return -b6(-sigma);
	}
};

struct Beta6 : public Alphas {
	Beta6( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return +b6(+sigma);
		else
			return -b7(-sigma);
	}
};

#endif
