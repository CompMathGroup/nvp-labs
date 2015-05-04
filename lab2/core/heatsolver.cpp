#include "heatsolver.h"
#include "linalg.h"

#include <cmath>

static HalfInteger I_2(1);

double HeatSolver::kappa(double u) {
	if (k < 1e-10)
		return 1;
	else {
		double v = std::fabs(u);
		if (v < tolerance)
			return std::pow(tolerance, k);
		return std::pow(v, k);
	}
}

// {{{ doFirstStep
double HeatSolver::doFirstStep(double C) {
	/* Compute u1 using analytical solution */
	u0.actualize(0);

	/* Time for woodoo */
	double T = 1;
	double kapmax = 0;
	for (int i = 0; i < u0.size(); i++) 
		kapmax = std::max(kapmax, kappa(u0.getExtrapolated(i, T)));

	dt = C * h * h / kapmax;

	u1.fillAt(dt); 
	t += dt;
	iter++;
	cou = dt / h / h;

	for (int j = -2; j <= N; j++)
		kap0[j] = kappa(u0[j]);

	return dt;
} // }}}

// {{{ doStep
double HeatSolver::doStep(const std::vector<Alphas *> &schemes) {
	/* u1 is incomplete, actualize it */
	u1.actualize(t);

	/* Kappa values are invalid for u1. Recompute them */
	double Citer = 0;
	for (int j = -2; j < N+1; j++) {
		kap1[j] = kappa(u1[j]);
		kap1[j] *= cou;
		double Cv = std::fabs(kap1[j]);
		if (Cv > Citer)
			Citer = Cv;
	}
	/* Control actual courant's number */
	if (Citer > Cmax)
		Cmax = Citer;
	
	predict();
	if (!schemes.empty())
		correct(*schemes.front());
	
	/* Cycle shift arrays */
	advance(u0, u1, u2);
	/* Swap kappa values for layer 0 and 1. 
	   Kappa values for layer 0 are now valid */
	swap(kap0, kap1);
	iter ++;
	t += dt;
	return Citer;
} // }}}

// {{{ predict 
void HeatSolver::predict() {
	double i2cou = .5 / cou;
	for (int j = -1; j < N - 1; j++)
	{
		f[j + I_2] = -i2cou * (kap1[j] + kap1[j+1]) * (u1[j+1] - u1[j]);
	}
	for (int j = 0; j < N - 1; j++) {
		u2[j] = u1[j] + cou * (f[j - I_2] - f[j + I_2]);
		if (u2[j] < tolerance && k > 1e-10)
			u2[j] = 0;
	}
	u2.actualize(t + dt);
	/* Compute kappa for predicted values u2 */
	for (int j = -2; j < N+1; j++) {
		kap2[j] = kappa(u2[j]);
		kap2[j] *= cou;
	}
} /// }}}

// {{{ correct
void HeatSolver::correct(const Alphas &scheme) {
	Gamma1 gamma1(scheme);

	Beta1 beta1(scheme);
	Beta2 beta2(scheme);
	Beta3 beta3(scheme);
	Beta4 beta4(scheme);

	double i2cou = .5 / cou;
	double icou = 1. / cou;

	/* First, compute f[j + I_2] and f[j - I_2] */
	for (int j = -1; j < N-1; j++) {
		f[j + I_2] = -i2cou * (kap1[j] + kap1[j+1]) * (u1[j+1] - u1[j]);
		double du1(u1[j  ] - u1[j-1]);
		double du2(u1[j+1] - u1[j  ]);
		double du3(u1[j+2] - u1[j+1]);

		f[j + I_2] += i2cou * gamma1(kap1[j-1]) * du1;
		f[j + I_2] += i2cou * gamma1(kap1[j  ]) * du1;

		f[j + I_2] -=  icou * gamma1(kap1[j  ]) * du2;
		f[j + I_2] -=  icou * gamma1(kap1[j+1]) * du2;

		f[j + I_2] += i2cou * gamma1(kap1[j+1]) * du3;
		f[j + I_2] += i2cou * gamma1(kap1[j+2]) * du3;

	}
	/* Next, compute ulower */
	for (int j = 0; j < N-1; j++) {
		ulower[j] = (u0[j] + u1[j]);

		ulower[j] += beta1(kap1[j-1]) * (u1[j-1] - u0[j-1]);
		ulower[j] += beta1(kap0[j-1]) * (u1[j-1] - u0[j-1]);

		ulower[j] += beta2(kap1[j  ]) * (u1[j  ] - u0[j  ]);
		ulower[j] += beta2(kap0[j  ]) * (u1[j  ] - u0[j  ]);
		
		ulower[j] += beta1(kap1[j+1]) * (u1[j+1] - u0[j+1]);
		ulower[j] += beta1(kap0[j+1]) * (u1[j+1] - u0[j+1]);
		
		ulower[j] += beta3(kap0[j-2]) * (u0[j-1] - u0[j-2]);
		ulower[j] += beta3(kap0[j-1]) * (u0[j-1] - u0[j-2]);
		
		ulower[j] += beta4(kap0[j-1]) * (u0[j  ] - u0[j-1]);
		ulower[j] += beta4(kap0[j  ]) * (u0[j  ] - u0[j-1]);

		ulower[j] += beta4(kap0[j  ]) * (u0[j  ] - u0[j+1]);
		ulower[j] += beta4(kap0[j+1]) * (u0[j  ] - u0[j+1]);
                                                          
		ulower[j] += beta3(kap0[j+1]) * (u0[j+1] - u0[j+2]);
		ulower[j] += beta3(kap0[j+2]) * (u0[j+1] - u0[j+2]);

		ulower[j] *= .5;
	}
	/* Compute uupper (without u2 layer) */
	for (int j = 0; j < N-1; j++) {
		uupper[j] = u1[j];

		uupper[j] -= beta1(kap2[j-1]) * (          u1[j-1]); //-=            u1[j-1] .lmul(W2[j-1], lam2[j-1], iW2[j-1], beta1);
		uupper[j] -= beta1(kap1[j-1]) * (          u1[j-1]); //-=            u1[j-1] .lmul(W1[j-1], lam1[j-1], iW1[j-1], beta1);

		uupper[j] -= beta2(kap2[j  ]) * (          u1[j  ]); //-=            u1[j  ] .lmul(W2[j  ], lam2[j  ], iW2[j  ], beta2);
		uupper[j] -= beta2(kap1[j  ]) * (          u1[j  ]); //-=            u1[j  ] .lmul(W1[j  ], lam1[j  ], iW1[j  ], beta2);

		uupper[j] -= beta1(kap2[j+1]) * (          u1[j+1]); //-=            u1[j+1] .lmul(W2[j+1], lam2[j+1], iW2[j+1], beta3);
		uupper[j] -= beta1(kap1[j+1]) * (          u1[j+1]); //-=            u1[j+1] .lmul(W1[j+1], lam1[j+1], iW1[j+1], beta3);

		uupper[j] += beta3(kap1[j-2]) * (u1[j-1] - u1[j-2]); //+= (u1[j-1] - u1[j-2]).lmul(W1[j-2], lam1[j-2], iW1[j-2], beta4);
		uupper[j] += beta3(kap1[j-1]) * (u1[j-1] - u1[j-2]); //+= (u1[j-1] - u1[j-2]).lmul(W1[j-1], lam1[j-1], iW1[j-1], beta4);

		uupper[j] += beta4(kap1[j-1]) * (u1[j  ] - u1[j-1]); //+= (u1[j  ] - u1[j-1]).lmul(W1[j-1], lam1[j-1], iW1[j-1], beta5);
		uupper[j] += beta4(kap1[j  ]) * (u1[j  ] - u1[j-1]); //+= (u1[j  ] - u1[j-1]).lmul(W1[j  ], lam1[j  ], iW1[j  ], beta5);

		uupper[j] += beta4(kap1[j  ]) * (u1[j  ] - u1[j+1]); //+= (u1[j+1] - u1[j  ]).lmul(W1[j  ], lam1[j  ], iW1[j  ], beta6);
		uupper[j] += beta4(kap1[j+1]) * (u1[j  ] - u1[j+1]); //+= (u1[j+1] - u1[j  ]).lmul(W1[j+1], lam1[j+1], iW1[j+1], beta6);
                                                          
		uupper[j] += beta3(kap1[j+1]) * (u1[j+1] - u1[j+2]); //+= (u1[j+2] - u1[j+1]).lmul(W1[j+1], lam1[j+1], iW1[j+1], beta7);
		uupper[j] += beta3(kap1[j+2]) * (u1[j+1] - u1[j+2]); //+= (u1[j+2] - u1[j+1]).lmul(W1[j+2], lam1[j+2], iW1[j+2], beta7);

		uupper[j] *= .5;
	}
	/* build system for u2 */
	for (int j = 0; j < N-1; j++)
	{
		double rh = cou * (f[j - I_2] - f[j + I_2]) + ulower[j] - uupper[j];

		B1[j] = .5 * (beta1(kap2[j-1]) + beta1(kap1[j-1]));

		B3[j] = .5 * (beta1(kap2[j+1]) + beta1(kap1[j+1]));

		B2[j] = .5 + .5 * (beta2(kap2[j]) + beta2(kap1[j]));

		rhs[j] = rh;
	}
	/* Solve system for u2 */
	u2.actualize(t + dt); /* for bc */

	rhs[0] -= B1[0] * u2[-1];
	rhs[N-2] -= B3[N-2] * u2[N-1];

	B1[0] = B3[N-2] = 0;

	TridiagonalSolve(B1, B2, B3, rhs, sol, false);
	/* Set u2 = sol */
	if (k > 1e-10) {
		for (int j = 0; j < N - 1; j++) {
			u2[j] = sol[j] < tolerance ? 0 : sol[j];
		}
	} else
		for (int j = 0; j < N - 1; j++)
			u2[j] = sol[j];
} // }}}
