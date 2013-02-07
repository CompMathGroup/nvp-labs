#include "gate.h"

#include <string.h>

#include <riemann.h>
#include <array.h>
#include <gdsolver.h>
#include <advsolver.h>

#include <mathlink.h>

#ifndef TOLERANCE
#define TOLERANCE 1e-6
#endif

bool linear;

static Riemann *rie = 0;
static GDSolver *sol = 0;
static AdvectionSolver *adv = 0;

static Array<Vars> *u0 = 0;
static Array<Vars> *u1 = 0;
static Array<Vars> *u2 = 0;

static Array<double> *w0 = 0;
static Array<double> *w1 = 0;
static Array<double> *w2 = 0;

static bool csproblem;
static double *x = 0;
static double *data = 0;

static double *ra;
static double *ua;
static double *ea;
static double *r;
static double *u;
static double *e;
static double *snum;

static bool firsttime;

std::vector<Alphas *> chain;

double step(double x) {
	if ((x < 0.4) || (x >= 0.6))
		return 0;
	return 1;
}

double cap(double x) {
	if ((x < 0.4) || (x >= 0.6))
		return 0;
	return 1 - std::pow(10*x - 5, 2);
}

double tri(double x) {
	if ((x < 0.4) || (x >= 0.6))
		return 0;
	return 1 - std::fabs(10*x - 5);
}

double wtf(double x) {
	return std::sin(100 * x);
}

double (*f)(double);

double InitializeSolver(
	const char *solver, const char *prob, int n, double C,
	double rL, double uL, double eL,
	double rR, double uR, double eR) 
{
	if (rie) delete rie;
	if (sol) delete sol;
	if (adv) delete adv;
	
	if (u0) delete u0;
	if (u1)	delete u1;
	if (u2)	delete u2;

	if (w0) delete w0;
	if (w1)	delete w1;
	if (w2)	delete w2;

	if (x) delete[] x;
	if (data) delete[] data;

	rie = 0; sol = 0; adv = 0;
	u0 = u1 = u2 = 0;
	w0 = w1 = w2 = 0;
	x = data = 0;

	linear = 0 == strcasecmp(solver, "linear");

	x = new double[n];
	double h = 1.0 / n;

	for (int i = 0; i < n; i++)
		x[i] = (i + 0.5) * h;

	if (!linear) {
		
		csproblem = 0 != strcasecmp(prob, "riemann");
		
		rie = new Riemann(GAMMA);
		rie->solve(rL, uL, eL, rR, uR, eR);

		data = new double[7 * n];
		ra = data;
		ua = data + n;
		ea = data + 2*n;
		r = data + 3*n;
		u = data + 4*n;
		e = data + 5*n;
		snum = data + 6*n;

		if (csproblem) {
			u0 = new CycledArray<Vars>(n, 2);
			u1 = new CycledArray<Vars>(n, 2);
			u2 = new CycledArray<Vars>(n, 2);
			sol = new GDSolver(n, true, *u0, *u1, *u2);

			bool isstep = 0 == strcasecmp(prob, "contacts");
			bool iscap =  0 == strcasecmp(prob, "contactc");
			bool istri =  0 == strcasecmp(prob, "contactt");
			f = isstep ? step : (iscap ? cap : (istri ? tri : wtf));

			/* Left -> inner, Right -> outer 
			* z = outer + (inner - outer) * f
			*/

			for (int i = 0; i < n; i++) {
				double rX = rR + (rL - rR) * f(x[i]);
				double uX = uR + (uL - uR) * f(x[i]);
				double eX = eR + (eL - eR) * f(x[i]);
				(*u0)[i] = Vars(Vector(rX, uX, eX));
			}
		} else {
			u0 = new ExtrapolatingArray<Vars>(n, 2);
			u1 = new ExtrapolatingArray<Vars>(n, 2);
			u2 = new ExtrapolatingArray<Vars>(n, 2);
			sol = new GDSolver(n, false, *u0, *u1, *u2);
			for (int i = 0; i < n / 2; i++)
				(*u0)[i] = Vars(Vector(rL, uL, eL));
			for (int i = n / 2; i < n; i++)
				(*u0)[i] = Vars(Vector(rR, uR, eR));
		}
		sol->tolerance = TOLERANCE;
		for (std::vector<Alphas *>::iterator j = chain.begin(); j != chain.end(); j++)
			delete *j;
		chain.clear();
		firsttime = true;
		return sol->doFirstStep(C);
	} else {
		data = new double[3 * n];
		
		ra = data;
		r = data + n;
		snum = data + 2*n;

		bool isstep = 0 == strcasecmp(prob, "step");
		bool iscap =  0 == strcasecmp(prob, "cap");
		bool istri =  0 == strcasecmp(prob, "triangle");
		f = isstep ? step : (iscap ? cap : (istri ? tri : wtf));

		w0 = new CycledArray<double>(n, 2);
		w1 = new CycledArray<double>(n, 2);
		w2 = new CycledArray<double>(n, 2);

		adv = new AdvectionSolver(n, *w0, *w1, *w2);
		adv->tolerance = TOLERANCE;
		sol = 0;

		for (int j = 0; j < n; j++) {
			double y = x[j] - C * h;
			y -= std::floor(y);
			(*w0)[j] = f(x[j]);
			(*w1)[j] = f(y);
		}
		chain.clear();
		firsttime = true;
		return adv->doFirstStep(C);
	}
}

int AddScheme(double *coeff, int sz) {
	if (sz != 10 * 2 * 16)
		return -1;
	Alphas *shm = new Alphas(
		Rational(&coeff[0 * 32], &coeff[0 * 32 + 16]),
		Rational(&coeff[1 * 32], &coeff[1 * 32 + 16]),
		Rational(&coeff[2 * 32], &coeff[2 * 32 + 16]),
		Rational(&coeff[3 * 32], &coeff[3 * 32 + 16]),
		Rational(&coeff[4 * 32], &coeff[4 * 32 + 16]),
		Rational(&coeff[5 * 32], &coeff[5 * 32 + 16]),
		Rational(&coeff[6 * 32], &coeff[6 * 32 + 16]),
		Rational(&coeff[7 * 32], &coeff[7 * 32 + 16]),
		Rational(&coeff[8 * 32], &coeff[8 * 32 + 16]),
		Rational(&coeff[9 * 32], &coeff[9 * 32 + 16])
	);
	chain.push_back(shm);
	return 0;
}

void DoSteps(int count) {
	if (firsttime) {
		firsttime = false;
		count--;
	}
	if (!linear) {
		for (;count;count--)
			sol->doStep(chain);

		double t = sol->getTime();
		if (csproblem)
			rie->twoCS(t, u1->size(), x, f, ra, ua, ea);
		else
			rie->evaluate(t, u1->size(), x, ra, ua, ea);

		const Array<int> &ss = sol->getSchemeNums();
		for (int i = 0; i < u1->size(); i++) {
			Vector w = (*u1)[i].to_nconserv();
			r[i] = std::isnormal(w(0)) ? w(0) : 0;
			u[i] = std::isnormal(w(1)) ? w(1) : 0;
			e[i] = std::isnormal(w(2)) ? w(2) : 0;
			snum[i] = ss[i];
		}

		int dims[2];
		char *heads[2] = {(char *)"List", (char *)"List"};
		int depth = 2;
		
		dims[0] = 7;
		dims[1] = u1->size();
		MLPutReal64Array(stdlink, data, dims, heads, depth);
	} else {
		for (;count;count--)
			adv->doStep(chain);

		double t = adv->getTime();
		for (int i = 0; i < w1->size(); i++) {
			r[i] = (*w1)[i];
		}

		const Array<int> &ss = adv->getSchemeNums();
		for (int j = 0; j < w1->size(); j++) {
			double y = x[j] - t;
			y -= std::floor(y);
			ra[j] = f(y);
			snum[j] = ss[j];
		}

		int dims[2];
		char *heads[2] = {(char *)"List", (char *)"List"};
		int depth = 2;
		
		dims[0] = 3;
		dims[1] = w1->size();
		MLPutReal64Array(stdlink, data, dims, heads, depth);
	}
}
