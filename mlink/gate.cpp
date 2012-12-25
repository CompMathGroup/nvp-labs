#include "gate.h"

#include <string.h>

#include <riemann.h>
#include <gd.h>
#include <array.h>
#include <solver.h>

#include <mathlink.h>

static Riemann *rie = 0;
static Solver *sol = 0;

static Array<Vars> *u0 = 0;
static Array<Vars> *u1 = 0;
static Array<Vars> *u2 = 0;

static bool csproblem;
static double *x = 0;
static double *data = 0;
static double *ra;
static double *ua;
static double *ea;
static double *r;
static double *u;
static double *e;

static bool firsttime;

std::vector<Alphas *> chain;

double InitializeSolver(const char *prob, int n, double C,
	double rL, double uL, double eL,
	double rR, double uR, double eR) {

	if (rie) delete rie;
	if (sol) delete sol;
	if (u0) delete u0;
	if (u1)	delete u1;
	if (u2)	delete u2;
	if (x) delete[] x;
	if (data) delete[] data;
	
	csproblem = 0 == strcasecmp(prob, "contact");
	
	rie = new Riemann(GAMMA);
	rie->solve(rL, uL, eL, rR, uR, eR);

	x = new double[n];
	data = new double[6 * n];
	ra = data;
	ua = data + n;
	ea = data + 2*n;
	r = data + 3*n;
	u = data + 4*n;
	e = data + 5*n;

	double h = 1.0 / n;

	for (int i = 0; i < n; i++)
		x[i] = (i + 0.5) * h;

	if (csproblem) {
		u0 = new CycledArray<Vars>(n, 2);
		u1 = new CycledArray<Vars>(n, 2);
		u2 = new CycledArray<Vars>(n, 2);
		sol = new Solver(n, true, *u0, *u1, *u2);
		for (int i = 0; i < 2 * n / 5; i++)
			(*u0)[i] = Vars(Vector(rR, uR, eR));
		for (int i = 2 * n / 5; i < 3 * n / 5; i++)
			(*u0)[i] = Vars(Vector(rL, uL, eL));
		for (int i = 3 * n / 5; i < n; i++)
			(*u0)[i] = Vars(Vector(rR, uR, eR));
	} else {
		u0 = new ExtrapolatingArray<Vars>(n, 2);
		u1 = new ExtrapolatingArray<Vars>(n, 2);
		u2 = new ExtrapolatingArray<Vars>(n, 2);
		sol = new Solver(n, true, *u0, *u1, *u2);
		for (int i = 0; i < n / 2; i++)
			(*u0)[i] = Vars(Vector(rL, uL, eL));
		for (int i = n / 2; i < n; i++)
			(*u0)[i] = Vars(Vector(rR, uR, eR));
	}
	for (std::vector<Alphas *>::iterator j = chain.begin(); j != chain.end(); j++)
		delete *j;
	chain.clear();
	firsttime = true;
	return sol->doFirstStep(C);
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
	for (;count;count--)
		sol->doStep(chain);

	double t = sol->getTime();
	if (csproblem)
		rie->twoCS(t, u1->size(), x, ra, ua, ea);
	else
		rie->evaluate(t, u1->size(), x, ra, ua, ea);
	for (int i = 0; i < u1->size(); i++) {
		Vector w = (*u1)[i].to_nconserv();
		r[i] = std::isnormal(w(0)) ? w(0) : 0;
		u[i] = std::isnormal(w(1)) ? w(1) : 0;
		e[i] = std::isnormal(w(2)) ? w(2) : 0;
	}

	int dims[2];
	char *heads[2] = {(char *)"List", (char *)"List"};
	int depth = 2;
	
	dims[0] = 6;
	dims[1] = u1->size();

	MLPutReal64Array(stdlink, data, dims, heads, depth);
}

