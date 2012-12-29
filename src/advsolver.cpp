#include "advsolver.h"
#include "linalg.h"

static HalfInteger I_2(1);

// {{{ doFirstStep
double AdvectionSolver::doFirstStep(double C) {
	/* Compute u1 using first-order scheme */
	u0.actualize();
	dt = C * h;
	cou = C;

	/* This is all, expect u1 to be already filled */
	for (int j = 0; j < N; j++)
		schemenum[j] = 0;

	t += dt;
	iter++;
	return dt;
} // }}}

// {{{ doStep
double AdvectionSolver::doStep(const std::vector<Alphas *> &schemes) {
	/* u1 is incomplete, actualize it */
	u1.actualize();

	/* Points that need update */
	std::vector<int> where(N);
	/* First pass should be done for all points */
	for (int j = 0; j < N; j++)
		where[j] = j;

	if (schemes.empty()) {
		predict();
		for (int j = 0; j < N; j++)
			schemenum[j] = 0;
	}

	int scnum = 1;
	for (std::vector<Alphas *>::const_iterator i = schemes.begin();
		i != schemes.end(); ++i) 
	{
		step(**i, where);
		for (std::vector<int>::const_iterator j = where.begin();
			j != where.end(); ++j) 
			schemenum[*j] = scnum;
		scnum ++;
		/* If any alternatives left, search for points 
		   that require recompute */
		if (i+1 != schemes.end()) 
			if (checkMono(where))
				break;
	}
	/* Cycle shift arrays */
	advance(u0, u1, u2);
	iter ++;
	t += dt;
	return cou;
} // }}}

// {{{ checkMono
bool AdvectionSolver::checkMono(std::vector<int> &where) {
	/*
	   | uC |
	---+----+
	uL | uD |
	*/
	double reltol = tolerance;
	double abstol = tolerance;
	where.clear();
	for (int j = 0; j < N; j++) {
		double wC = u2[j  ];
		double wL = u1[j-1];
		double wD = u1[j  ];

		bool need = false;

		double wmin = std::min(wL, wD);
		double wmax = std::max(wL, wD);
		double wdiff = wmax - wmin;
		wmin = wmin - reltol * wdiff - abstol;
		wmax = wmax + reltol * wdiff + abstol;
		double wmid = wC;

		if (wmid < wmin || wmid > wmax)
			need = true;

		if (need)
			where.push_back(j);
	}

	return where.empty();
}// }}}

/* only a stub for case when no scheme is specified */
void AdvectionSolver::predict() {
	for (int j = 0; j < N; j++)
		u2[j] = cou * u1[j-1] + (1 - cou) * u1[j];
	u2.actualize();
}

// {{{ step
void AdvectionSolver::step(const Alphas &scheme, const std::vector<int> &where) {
	Gamma1 _gamma1(scheme);
	Gamma2 _gamma2(scheme);
	Gamma3 _gamma3(scheme);

	Beta1 _beta1(scheme);
	Beta2 _beta2(scheme);
	Beta3 _beta3(scheme);
	Beta4 _beta4(scheme);
	Beta5 _beta5(scheme);
	Beta6 _beta6(scheme);
	Beta7 _beta7(scheme);

	double gamma1 = _gamma1(cou);
	double gamma2 = _gamma2(cou);
	double gamma3 = _gamma3(cou);

	double beta1 = _beta1(cou);
	double beta2 = _beta2(cou);
	double beta3 = _beta3(cou);
	double beta4 = _beta4(cou);
	double beta5 = _beta5(cou);
	double beta6 = _beta6(cou);
	double beta7 = _beta7(cou);

	double icou = 1. / cou;

	/* First, compute f[j + I_2] and f[j - I_2] */
	for (int j = -1; j < N; j++) {
		f[j + I_2] = .5 * (u1[j] + u1[j+1]);
		double du1(u1[j] - u1[j-1]);
		double du2(u1[j+1] - u1[j]);
		double du3(u1[j+2] - u1[j+1]);

		f[j + I_2] += icou * du1 * gamma1;
		f[j + I_2] += icou * du2 * gamma2;
		f[j + I_2] += icou * du3 * gamma3;
	}
	/* Next, compute ulower */
	for (int j = 0; j < N; j++) {
		ulower[j] = .5 * (u0[j] + u1[j]);

		ulower[j] += (u1[j-1] - u0[j-1]) * beta1;
		ulower[j] += (u1[j  ] - u0[j  ]) * beta2;
		ulower[j] += (u1[j+1] - u0[j+1]) * beta3;
		ulower[j] += (u0[j-1] - u0[j-2]) * beta4;
		ulower[j] += (u0[j  ] - u0[j-1]) * beta5;
		ulower[j] += (u0[j+1] - u0[j  ]) * beta6;
		ulower[j] += (u0[j+2] - u0[j+1]) * beta7;
	}
	/* Compute uupper (without u2 layer) */
	for (int j = 0; j < N; j++) {
		uupper[j] = .5 * u1[j];

		uupper[j] -=            u1[j-1]  * beta1;
		uupper[j] -=            u1[j  ]  * beta2;
		uupper[j] -=            u1[j+1]  * beta3;
		uupper[j] += (u1[j-1] - u1[j-2]) * beta4;
		uupper[j] += (u1[j  ] - u1[j-1]) * beta5;
		uupper[j] += (u1[j+1] - u1[j  ]) * beta6;
		uupper[j] += (u1[j+2] - u1[j+1]) * beta7;
	}
	/* build system for u2 */
	for (std::vector<int>::const_iterator i = where.begin();
		i != where.end(); ++i) 
	{
		int j = *i;
		double rh = cou * (f[j - I_2] - f[j + I_2]) + ulower[j] - uupper[j];

		B1[j] = beta1;
		B2[j] = .5 + beta2;
		B3[j] = beta3;

		rhs[j] = rh;
	}
	/* Solve system for u2 */
	EliminateCycled(B1, B2, B3, rhs, sol);
	TridiagonalSolve(B1, B2, B3, rhs, sol, true);
	/* Set u2 = sol */
	for (int j = 0; j < N; j++)
		u2[j] = sol[j];
	u2.actualize();
} // }}}
