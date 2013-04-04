#include "gate.h"

#include <string.h>

#include <array.h>
#include <heatsolver.h>

#include <mathlink.h>
#include <cmath>

static HeatSolver *sol = 0;

static Array<double> *u0 = 0;
static Array<double> *u1 = 0;
static Array<double> *u2 = 0;

static double *x = 0;
static double *data = 0;

static double *ua;
static double *u;

static bool firsttime;

std::vector<Alphas *> chain;

class Solution {
};

class ExpWave : public Solution {
public:	
	double operator()(double x, double t) const {
		const double c = 2.;
		return std::exp(c * (c * (t - 1.) - x));
	}
};

class Smooth : public Solution {
public:
	double operator()(double x, double t) const {
		if (t > 0)
			return 0.5 - 0.5 * erf(.5 * (x - .5) / std::sqrt(t));
		return (x < .499999) ? 1 : (x > 0.500001 ? 0 : .5);
	}
};

class FixedWave : public Solution {
	const double k, tmax;
public:	
	FixedWave(const double _k, const double _c) : k(_k), tmax(_c) {}
	double operator()(double x, double t) const {
		if (x < .5)
			return std::pow(.5 * k * (.5 - x) * (.5 - x) / (k + 2.) / (tmax - t), 1. / k);
		return 0.;
	}
};

class ConstSpeed : public Solution {
	const double k, c; 
public:	
	ConstSpeed(const double _k, const double _c) : k(_k), c(_c) {}
	double operator()(double x, double t) const {
		if (x > 1 - 1e-10)
			return 0; /* right bc */
		if (x < c * t) 
			return std::pow(c * k * (c * t - x), 1. / k);
		return 0.;
	}
};

class ConstVal : public Solution {
	const double k, c;
	double c1, c2, c3;
	double xi0;
	double sum(double eta) const {
		double z = 1 - eta;
		double z2 = z * z;
		double z3 = z2 * z;
		return 1 + c1 * z + c2 * z2 + c3 * z3;
	}
public:	
	ConstVal(const double _k, const double _c) : k(_k), c(_c) {
		c1 = -.5 / k;
		c2 = (4 * k * k + 6 * k + 3) / (24 * k * k * (2 * k * k + 3 * k + 1));
		c3 = (2 * k * k * k - 2 * k * k - 3 * k - 1) / (48 * k * k * k * (3 * k + 1) * (2 * k * k + 3 * k + 1));
		xi0 = std::pow(std::pow(.5 * k, 1 + 1. / k) * sum(0), - .5 * k / (k + 1));
	}
	double operator()(double x, double t) const {
		if (t < 1e-10)
			return 0.;
		if (x > 1 - 1e-10)
			return 0.;
		double xi = x / std::sqrt(t);
		if (xi < xi0) 
			return std::pow(.5 * k * xi0 * (xi0 - xi), 1. / k) * std::pow(sum(xi / xi0), 1. / (k + 1));
		return 0.;
	}
};

template <class T, class BC>
class Extrapolator {
	const double h;
	const BC sol;
public:
	Extrapolator(const double _h, BC _sol) : h(_h), sol(_sol) {}
	T operator()(const int idx, const double t) {
		return sol(idx * h + h, t);
	}
};

double InitializeSolver(
	const char *solver, double k, double c, const char *problem, int n, double cou) 
{
	if (sol) delete sol;
	
	if (u0) delete u0;
	if (u1)	delete u1;
	if (u2)	delete u2;

	if (x) delete[] x;
	if (data) delete[] data;

	sol = 0; 
	u0 = u1 = u2 = 0;
	x = data = 0;

	if (0 == strcasecmp(solver, "linear"))
		k = 0;

	x = new double[n+1];
	double h = 1. / n;

	for (int i = 0; i <= n; i++)
		x[i] = i * h;

	data = new double[2 * (n + 1)];
	ua = data;
	u = data + n + 1;

#define DECL(x, y, z) \
		x = new ExtrapolatingArray<double, Extrapolator<double, y> >(n - 1, 2, Extrapolator<double, y>(h, y z))

	if (0 == strcasecmp(problem, "exp")) {
		DECL(u0, ExpWave, ());
		DECL(u1, ExpWave, ());
		DECL(u2, ExpWave, ());
	} else if (0 == strcasecmp(problem, "smooth")) {
		DECL(u0, Smooth, ());
		DECL(u1, Smooth, ());
		DECL(u2, Smooth, ());
	} else if (0 == strcasecmp(problem, "fixed")) {
		DECL(u0, FixedWave, (k, c));
		DECL(u1, FixedWave, (k, c));
		DECL(u2, FixedWave, (k, c));
	} else if (0 == strcasecmp(problem, "constc")) {
		DECL(u0, ConstSpeed, (k, c));
		DECL(u1, ConstSpeed, (k, c));
		DECL(u2, ConstSpeed, (k, c));
	} else if (0 == strcasecmp(problem, "constt")) {
		DECL(u0, ConstVal, (k, c));
		DECL(u1, ConstVal, (k, c));
		DECL(u2, ConstVal, (k, c));
	}

	sol = new HeatSolver(n, k, *u0, *u1, *u2);

	u0->fillAt(0);

	for (std::vector<Alphas *>::iterator j = chain.begin(); j != chain.end(); j++)
		delete *j;
	chain.clear();
	firsttime = true;
	return sol->doFirstStep(cou);
}

int AddScheme(double *coeff, int sz) {
	if (sz != 5 * 2 * 16)
		return -1;
	Alphas *shm = new Alphas(
		Rational(&coeff[0 * 32], &coeff[0 * 32 + 16]),
		Rational(&coeff[1 * 32], &coeff[1 * 32 + 16]),
		Rational(&coeff[2 * 32], &coeff[2 * 32 + 16]),
		Rational(&coeff[3 * 32], &coeff[3 * 32 + 16]),
		Rational(&coeff[4 * 32], &coeff[4 * 32 + 16])
	);
	chain.push_back(shm);
	return 0;
}

void DoSteps(int count) {
	if (firsttime) {
		firsttime = false;
		count--;
	}
	for (;count;count--)
		sol->doStep(chain);

	double t = sol->getTime();

	for (int i = -1; i <= u1->size(); i++) {
		double w = (*u1)[i];
		u[i+1] = std::isnormal(w) ? w : 0;
		w = u1->getExtrapolated(i, t);
		ua[i+1] = std::isnormal(w) ? w : 0;
	}

	int dims[2];
	char *heads[2] = {(char *)"List", (char *)"List"};
	int depth = 2;
	
	dims[0] = 2;
	dims[1] = u1->size() + 2;
	MLPutReal64Array(stdlink, data, dims, heads, depth);
}
