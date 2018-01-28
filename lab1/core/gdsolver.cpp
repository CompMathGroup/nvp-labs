#include "gdsolver.h"
#include "linalg.h"

static HalfInteger I_2(1);

// {{{ doFirstStep
double GDSolver::doFirstStep(double C) {
	/* Compute u1 using first-order scheme */
	u0.actualize();
	double lmax = 0;
	for (int j = -2; j < N+2; j++) {
		F1[j] = u0[j].f();
		u0[j].eigen(W0[j], lam0[j], iW0[j]);
		double aupc = u0[j].lmax(lam0[j]);
		if (aupc > lmax)
			lmax = aupc;
	}
	dt = C * h / lmax;
	/* Note, this dt will be used for ALL steps, not only the first one */
	/* Perform some iterations to ensure dt is valid for 
	   at least two first steps */
	while (1)
	{
		cou = dt / h;
		double icou = 1. / cou;
		double i4cou = .25 * icou;
		for (int j = -2; j < N+2; j++) 
			lam0[j] *= cou;
		for (int j = -1; j < N; j++)
		{
			f[j + I_2] = .5 * (F1[j] + F1[j+1]);
			Vars udiff(u0[j+1] - u0[j]);
			f[j + I_2] -= i4cou * 
				udiff.lmul(W0[j  ], lam0[j  ], iW0[j  ], AbsoluteValue());
			f[j + I_2] -= i4cou * 
				udiff.lmul(W0[j+1], lam0[j+1], iW0[j+1], AbsoluteValue());
		}
		for (int j = 0; j < N; j++) {
			u1[j] = u0[j] + cou * (f[j - I_2] - f[j + I_2]);
			u1[j].fixup();
			schemenum[j] = 0;
		}

		lmax = 0;
		for (int j = 0; j < N; j++) {
			double aupc = u1[j].lmax(); 
			if (aupc > lmax)
				lmax = aupc;
		}
		double dtnew = C * h / lmax;
	//	printf("dt = %e, dtnew = %e\n", dt, dtnew);
		if (dtnew > 0.999 * dt)
			break;
		for (int j = -2; j < N+2; j++) 
			lam0[j] *= icou;
		dt = dtnew;
	}
	t += dt;
	iter++;
	return dt;
} // }}}

// {{{ doStep
double GDSolver::doStep(const std::vector<Alphas *> &schemes) {
	/* u1 is incomplete, actualize it */
	u1.actualize();

	/* Eigendecomposition is invalid for u1. Recompute it */
	double Citer = 0;
	for (int j = -2; j < N+2; j++) {
		u1[j].eigen(W1[j], lam1[j], iW1[j]);
		F1[j] = u1[j].f();
		lam1[j] *= cou;
		double Cv = u1[j].lmax(lam1[j]);
		if (Cv > Citer)
			Citer = Cv;
	}
	/* Control actual courant's number */
	if (Citer > Cmax)
		Cmax = Citer;
	
	/* Points that need update */
	std::vector<int> where(N);
	/* First pass should be done for all points */
	for (int j = 0; j < N; j++)
		where[j] = j;
	int scnum = 1;
	if (schemes.empty()) {
		eno();
		for (int j = 0; j < N; j++)
			schemenum[j] = 0;
	}
	for (std::vector<Alphas *>::const_iterator 
		i = schemes.begin();
		i != schemes.end(); ++i) 
	{
		/* We are going to predict u2 using first order scheme */
		predict();
		/* Now we have all needed data to correct u2 */
		correct(**i, where);
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
	/* Swap Eigendecomposition for layer 0 and 1. 
	   Decomposition for layer 0 is now valid */
	swap(W0, W1);
	swap(iW0, iW1);
	swap(lam0, lam1);
	iter ++;
	t += dt;
	return Citer;
} // }}}

// {{{ checkMono
bool GDSolver::checkMono(std::vector<int> &where) {
	/*
	   | uC |   
	---+----+---
	uL | uD | uR
	*/
	double reltol = tolerance;
	double abstol = tolerance;
	where.clear();
	for (int j = 0; j < N; j++) {

		Vector uC = iW2[j] * u2[j  ].to_vect();
		Vector uL = iW2[j] * u1[j-1].to_vect();
		Vector uD = iW2[j] * u1[j  ].to_vect();
		Vector uR = iW2[j] * u1[j+1].to_vect();

		bool need = false;
		for (int k = 0; k < 3; k++) {
			double wL = uL(k);
			double wR = uR(k);
			double wD = uD(k);
			double wC = uC(k);
			double u = lam2[j](k);
			if (u > 0) {
				double wmin = std::min(wL, wD);
				double wmax = std::max(wL, wD);
				double wdiff = wmax - wmin;
				wmin = wmin - reltol * wdiff - abstol;
				wmax = wmax + reltol * wdiff + abstol;
				double wmid = wC;

				if (wmid < wmin || wmid > wmax)
					need = true;
			}
			if (u < 0) {
				double wmin = std::min(wR, wD);
				double wmax = std::max(wR, wD);
				double wdiff = wmax - wmin;
				wmin = wmin - reltol * wdiff - abstol;
				wmax = wmax + reltol * wdiff + abstol;
				double wmid = wC;

				if (wmid < wmin || wmid > wmax)
					need = true;
			}
		}
		if (need)
			where.push_back(j);
	}

	return where.empty();
}// }}}

inline double sqr(double x) {
    return x*x;
}

void eno3(
	const Vars U[5],
	const Matrix Omega[4],
	const Matrix invOmega[4],
	Vars &Uplus,
	Vars &Uminus
	)
{
	Vector dW[4];
	Vector AdWp[4];
	Vector AdWm[4];

	Vector Up, Um;

	for (int i = 0; i < 4; i++) {
		dW[i] = Omega[i] * (U[i+1] - U[i]).to_vect();
	}

	for (int i = 0; i < 3; i++) {
		double dmm = dW[0](i);
		double dm  = dW[1](i);
		double dp  = dW[2](i);
		double dpp = dW[3](i);

		double b1 = (3 * sqr(3*dm - dmm) + 13 * sqr(dm - dmm)) / 48;
		double b2 = (3 * sqr(dm + dp) + 13 * sqr(dp - dm)) / 48;
		double b3 = (3 * sqr(3*dp - dpp) + 13 * sqr(dpp - dp)) / 48;

		int best = 0;
		double bmin = b1;
		if (b2 < bmin) {
			best = 1;
			bmin = b2;
		}
		if (b3 < bmin) {
			best = 2;
			bmin = b3;
		}

		double a[3] = {0, 0, 0};
		a[best] = 1;
		AdWp[0](i) = -a[0] / 3 * dmm;
		AdWp[1](i) = (5*a[0]+a[1]) / 6 * dm;
		AdWp[2](i) = (a[1]+2*a[2]) / 3 * dp;
		AdWp[3](i) = -a[2] / 6 * dpp;

		AdWm[0](i) = a[0] / 6 * dmm;
		AdWm[1](i) = -(2*a[0]+a[1]) / 3 * dm;
		AdWm[2](i) = -(a[1]+5*a[2]) / 6 * dp;
		AdWm[3](i) = a[2] / 3 * dpp;
	}

	Up = U[2].to_vect();
	Um = U[2].to_vect();

	for (int i = 0; i < 4; i++) {
		Up += invOmega[i] * AdWp[i];
		Um += invOmega[i] * AdWm[i];
	}

	Uplus.rho   = Up(0);
	Uplus.P     = Up(1);
	Uplus.E     = Up(2);
	Uplus.GAMMA = U[2].GAMMA;

	Uminus.rho   = Um(0);
	Uminus.P     = Um(1);
	Uminus.E     = Um(2);
	Uminus.GAMMA = U[2].GAMMA;
}

void euler_step(const int N, const Vars *U, Vars *Un) {
	std::vector<Vars> Up(N), Um(N);
	std::vector<Matrix> Omega(N+3), invOmega(N+3)
	std::vector<Vector> lambda(N+3);	

	for (int j = 0; j < N+3; j++) {
		Vars Umid = .5 * (U[j] + U[j+1]);
		Umid.eigen(invOmega[j], lambda[j], Omega[j]);
	}

	for (int j = 0; j < N; j++) {
		eno3(&U[j], &Omega[j], &invOmega[j],
			Up[j], Um[j]);
	}
}

// {{{ ENO
void GDSolver::eno() {
		
} /// }}}

// {{{ predict 
void GDSolver::predict() {
	double i4cou = .25 / cou;
	for (int j = -1; j < N; j++)
	{
		f[j + I_2] = .5 * (F1[j] + F1[j+1]);
		Vars udiff(u1[j+1] - u1[j]);
		f[j + I_2] -= i4cou * 
			udiff.lmul(W1[j  ], lam1[j  ], iW1[j  ], AbsoluteValue());
		f[j + I_2] -= i4cou * 
			udiff.lmul(W1[j+1], lam1[j+1], iW1[j+1], AbsoluteValue());
	}
	for (int j = 0; j < N; j++) {
		u2[j] = u1[j] + cou * (f[j - I_2] - f[j + I_2]);
		u2[j].fixup();
	}
	u2.actualize();
	/* Compute eigendecomposition for predicted values u2 */
	for (int j = -2; j < N+2; j++) {
		u2[j].eigen(W2[j], lam2[j], iW2[j]);
		lam2[j] *= cou;
	}
} /// }}}

// {{{ correct
void GDSolver::correct(const Alphas &scheme, const std::vector<int> &where) {
	Gamma1 gamma1(scheme);
	Gamma2 gamma2(scheme);
	Gamma3 gamma3(scheme);

	Beta1 beta1(scheme);
	Beta2 beta2(scheme);
	Beta3 beta3(scheme);
	Beta4 beta4(scheme);
	Beta5 beta5(scheme);
	Beta6 beta6(scheme);
	Beta7 beta7(scheme);

	double i2cou = .5 / cou;

	/* First, compute f[j + I_2] and f[j - I_2] */
	for (int j = -1; j < N; j++) {
		f[j + I_2] = .5 * (F1[j] + F1[j+1]);
		Vars du1(u1[j] - u1[j-1]);
		Vars du2(u1[j+1] - u1[j]);
		Vars du3(u1[j+2] - u1[j+1]);

		f[j + I_2] += i2cou * du1.lmul(W1[j-1], lam1[j-1], iW1[j-1], gamma1);
		f[j + I_2] += i2cou * du1.lmul(W1[j  ], lam1[j  ], iW1[j  ], gamma1);

		f[j + I_2] += i2cou * du2.lmul(W1[j  ], lam1[j  ], iW1[j  ], gamma2);
		f[j + I_2] += i2cou * du2.lmul(W1[j+1], lam1[j+1], iW1[j+1], gamma2);

		f[j + I_2] += i2cou * du3.lmul(W1[j+1], lam1[j+1], iW1[j+1], gamma3);
		f[j + I_2] += i2cou * du3.lmul(W1[j+2], lam1[j+2], iW1[j+2], gamma3);
	}
	/* Next, compute ulower */
	for (int j = 0; j < N; j++) {
		ulower[j] = (u0[j] + u1[j]);

		ulower[j] += (u1[j-1] - u0[j-1]).lmul(W1[j-1], lam1[j-1], iW1[j-1], beta1);
		ulower[j] += (u1[j-1] - u0[j-1]).lmul(W0[j-1], lam0[j-1], iW0[j-1], beta1);

		ulower[j] += (u1[j  ] - u0[j  ]).lmul(W1[j  ], lam1[j  ], iW1[j  ], beta2);
		ulower[j] += (u1[j  ] - u0[j  ]).lmul(W0[j  ], lam0[j  ], iW0[j  ], beta2);

		ulower[j] += (u1[j+1] - u0[j+1]).lmul(W1[j+1], lam1[j+1], iW1[j+1], beta3);
		ulower[j] += (u1[j+1] - u0[j+1]).lmul(W0[j+1], lam0[j+1], iW0[j+1], beta3);

		ulower[j] += (u0[j-1] - u0[j-2]).lmul(W0[j-2], lam0[j-2], iW0[j-2], beta4);
		ulower[j] += (u0[j-1] - u0[j-2]).lmul(W0[j-1], lam0[j-1], iW0[j-1], beta4);

		ulower[j] += (u0[j  ] - u0[j-1]).lmul(W0[j-1], lam0[j-1], iW0[j-1], beta5);
		ulower[j] += (u0[j  ] - u0[j-1]).lmul(W0[j  ], lam0[j  ], iW0[j  ], beta5);

		ulower[j] += (u0[j+1] - u0[j  ]).lmul(W0[j  ], lam0[j  ], iW0[j  ], beta6);
		ulower[j] += (u0[j+1] - u0[j  ]).lmul(W0[j+1], lam0[j+1], iW0[j+1], beta6);

		ulower[j] += (u0[j+2] - u0[j+1]).lmul(W0[j+1], lam0[j+1], iW0[j+1], beta7);
		ulower[j] += (u0[j+2] - u0[j+1]).lmul(W0[j+2], lam0[j+2], iW0[j+2], beta7);

		ulower[j] *= .5;
	}
	/* Compute uupper (without u2 layer) */
	for (int j = 0; j < N; j++) {
		uupper[j] = u1[j];

		uupper[j] -=            u1[j-1] .lmul(W2[j-1], lam2[j-1], iW2[j-1], beta1);
		uupper[j] -=            u1[j-1] .lmul(W1[j-1], lam1[j-1], iW1[j-1], beta1);

		uupper[j] -=            u1[j  ] .lmul(W2[j  ], lam2[j  ], iW2[j  ], beta2);
		uupper[j] -=            u1[j  ] .lmul(W1[j  ], lam1[j  ], iW1[j  ], beta2);

		uupper[j] -=            u1[j+1] .lmul(W2[j+1], lam2[j+1], iW2[j+1], beta3);
		uupper[j] -=            u1[j+1] .lmul(W1[j+1], lam1[j+1], iW1[j+1], beta3);

		uupper[j] += (u1[j-1] - u1[j-2]).lmul(W1[j-2], lam1[j-2], iW1[j-2], beta4);
		uupper[j] += (u1[j-1] - u1[j-2]).lmul(W1[j-1], lam1[j-1], iW1[j-1], beta4);

		uupper[j] += (u1[j  ] - u1[j-1]).lmul(W1[j-1], lam1[j-1], iW1[j-1], beta5);
		uupper[j] += (u1[j  ] - u1[j-1]).lmul(W1[j  ], lam1[j  ], iW1[j  ], beta5);

		uupper[j] += (u1[j+1] - u1[j  ]).lmul(W1[j  ], lam1[j  ], iW1[j  ], beta6);
		uupper[j] += (u1[j+1] - u1[j  ]).lmul(W1[j+1], lam1[j+1], iW1[j+1], beta6);

		uupper[j] += (u1[j+2] - u1[j+1]).lmul(W1[j+1], lam1[j+1], iW1[j+1], beta7);
		uupper[j] += (u1[j+2] - u1[j+1]).lmul(W1[j+2], lam1[j+2], iW1[j+2], beta7);

		uupper[j] *= .5;
	}
	/* build system for u2 */
	for (std::vector<int>::const_iterator i = where.begin();
		i != where.end(); ++i) 
	{
		int j = *i;
		Vars rh = cou * (f[j - I_2] - f[j + I_2]) + ulower[j] - uupper[j];

		B1[j]  = Matrix(W2[j-1], lam2[j-1], iW2[j-1], beta1);
		B1[j] += Matrix(W1[j-1], lam1[j-1], iW1[j-1], beta1);
		B1[j] *= .5;

		B3[j]  = Matrix(W2[j+1], lam2[j+1], iW2[j+1], beta3);
		B3[j] += Matrix(W1[j+1], lam1[j+1], iW1[j+1], beta3);
		B3[j] *= .5;

		B2[j]  = Matrix(W2[j  ], lam2[j  ], iW2[j  ], beta2);
		B2[j] += Matrix(W1[j  ], lam1[j  ], iW1[j  ], beta2);
		B2[j] += 1;
		B2[j] *= .5;

		rhs[j](0) = rh.rho;
		rhs[j](1) = rh.P;
		rhs[j](2) = rh.E;
	}
	/* Solve system for u2 */
	if (cycled)
		EliminateCycled(B1, B2, B3, rhs, sol);
	TridiagonalSolve(B1, B2, B3, rhs, sol, cycled);
	/* Set u2 = sol */
	for (int j = 0; j < N; j++) {
		u2[j].rho = sol[j](0);
		u2[j].P   = sol[j](1);
		u2[j].E   = sol[j](2);
		u2[j].fixup();
	}
	u2.actualize();
} // }}}
