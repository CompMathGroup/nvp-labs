#ifndef __COEFF_H__
#define __COEFF_H__

#include <iostream>
#include <cmath>

struct Rational {
	double a0, a1, a2, a3;
	double b0, b1, b2, b3;
	Rational(double x = 0, double y = 0, double z = 0, double w = 0)
		: a0(x), a1(y), a2(z), a3(w),
		  b0(1), b1(0), b2(0), b3(0)
		{}
	Rational(
		double a0_, double a1_, double a2_, double a3_,
		double b0_, double b1_, double b2_, double b3_)
		: a0(a0_), a1(a1_), a2(a2_), a3(a3_),
		  b0(b0_), b1(b1_), b2(b2_), b3(b3_)
		{}
	double operator()(double sigma) const {
		double nom, denom;
		nom = 0;
		nom = nom * sigma + a3;
		nom = nom * sigma + a2;
		nom = nom * sigma + a1;
		nom = nom * sigma + a0;
		denom = 0;
		denom = denom * sigma + b3;
		denom = denom * sigma + b2;
		denom = denom * sigma + b1;
		denom = denom * sigma + b0;
		return nom / denom;
	}
};

struct Alphas {
	const Rational          aup_1,   aup0,   aup1;
	const Rational a_2,     a_1,     a0,     a1,     a2;
	const Rational adown_2, adown_1, adown0, adown1, adown2;
	Alphas(
		const Rational &aup_1,
		const Rational &aup0,
		const Rational &aup1,

		const Rational &a_2,
		const Rational &a_1,
		const Rational &a0,
		const Rational &a1,
		const Rational &a2,

		const Rational &adown_2,
		const Rational &adown_1,
		const Rational &adown0,
		const Rational &adown1,
		const Rational &adown2
	) : 
		aup_1(aup_1), aup0(aup0), aup1(aup1),
		a_2(a_2), a_1(a_1), a0(a0), a1(a1), a2(a2),
		adown_2(adown_2), adown_1(adown_1), adown0(adown0), adown1(adown1), adown2(adown2)
	{}
	bool sanity() {
		double s = .25;
		double sup = aup_1(s) + aup0(s) + aup1(s);
		double smid = a_2(s) + a_1(s) + a0(s) + a1(s) + a2(s);
		double sdown = adown_2(s) + adown_1(s) + adown0(s) + adown1(s) + adown2(s);

		double msup = -aup_1(s) + aup1(s);
		double msmid = -2*a_2(s) - a_1(s) + a1(s) + 2*a2(s);
		double msdown = -2*adown_2(s) - adown_1(s) + adown1(s) + 2*adown2(s);

		bool sumzero;
		sumzero = std::fabs(sup + smid + sdown) < 1e-12;
		if (!sumzero)
			std::cerr << "Alphas not sum up to zero! sum a = " << sup + smid + sdown << std::endl;
		bool norm;
		norm = std::fabs(sup - sdown - 1) < 1e-12;
		if (!norm)
			std::cerr << "Alphas are not normalized! sum nu a = " << sup - sdown << std::endl;
		bool firstord;
		firstord = std::fabs(msup + msmid + msdown - s) < 1e-12;
		if (!firstord)
			std::cerr << "Alphas do not provide first order scheme! (sum mu a)/s = " << (msup + msdown + msmid)/s << std::endl;
		return sumzero && norm && firstord;
	}
};

struct Gamma1 : public Alphas {
	Gamma1( const Alphas &o) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return adown_2(sigma) + a_2(sigma);
		else
			return adown2(-sigma) + a2(-sigma);
	}
};

struct Gamma2 : public Alphas {
	Gamma2( const Alphas &o) : Alphas(o) {}
	double operator()(double sigma) const {
		double s = std::fabs(sigma);
		return .5 * (adown_2(s) + adown2(s) + a_2(s) + a2(s) - adown0(s) - a0(s) - aup0(s));
	}
};

struct Gamma3 : public Alphas {
	Gamma3( const Alphas &o) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return adown2(sigma) + a2(sigma);
		else
			return adown_2(-sigma) + a_2(-sigma);
	}
};

struct Beta1 : public Alphas {
	Beta1( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return aup_1(sigma);
		else
			return aup1(-sigma);
	}
};

struct Beta3 : public Alphas {
	Beta3( const Alphas &o) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return aup1(sigma);
		else
			return aup_1(-sigma);
	}
};

struct Beta2 : public Alphas {
	Beta2( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		double s = std::fabs(sigma);
		return aup0(s) - .5;
	}
};

struct Beta5 : public Alphas {
	Beta5( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return adown_2(sigma) + adown_1(sigma) - aup_1(sigma);
		else
			return adown2(-sigma) + adown1(-sigma) - aup1(-sigma);
	}
};

struct Beta6 : public Alphas {
	Beta6( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return -adown2(sigma) - adown1(sigma) + aup1(sigma);
		else
			return -adown_2(-sigma) - adown_1(-sigma) + aup_1(-sigma);
	}
};

struct Beta4 : public Alphas {
	Beta4( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return adown_2(sigma);
		else
			return adown2(-sigma);
	}
};

struct Beta7 : public Alphas {
	Beta7( const Alphas &o ) : Alphas(o) {}
	double operator()(double sigma) const {
		if (sigma > 0)
			return -adown2(sigma);
		else
			return -adown_2(-sigma);
	}
};

#endif
