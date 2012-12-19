#include "coeff.h"
#include "main.h"

#include <stdlib.h>
#include <stdio.h>
#include <fenv.h>
#include <cmath>
#include <iostream>

#include "schemes.h"

#define SHOCK

int main(int argc, char **argv) {

	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

	const int N = 200;
	const double C = 0.8;

#ifdef SHOCK
	ExtrapolatingArray<Vars> u0(N, 2), u1(N, 2), u2(N, 2);
#else
	CycledArray<Vars> u0(N, 2), u1(N, 2), u2(N, 2);
#endif

	selfTest();

/*
	fprintf(stderr, "s,g1,g2,g3,b1,b2,b3,b4,b5,b6,b7\n");
	for (double s = -1; s <= 1.01; s += 0.02)
		fprintf(stderr, "%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
			s, gamma1(s), gamma2(s), gamma3(s),
			beta1(s), beta2(s), beta3(s), 
			beta4(s), beta5(s), beta6(s), beta7(s));
*/
/*
	Shock tube test
	(r1,u1,p1) = (1    ,0,1)   ===> (r1,P1,E1) = (1    ,0,2.5)
	(r1,u1,p1) = (0.125,0,0.1) ===> (r2,P2,E2) = (0.125,0,0.25)
*/

/*
	Contact surface
	(r1,u1,p1) = (1  ,1,1) ===> (r1,P1,E1) = (1  ,1  ,3)
	(r2,u2,p2) = (0.1,1,1) ===> (r2,P2,E2) = (0.1,0.1,2.55)
*/

#ifdef SHOCK
	for (int i = 0; i < N / 2; i++) {
		u0[i] = Vars(1, 0, 2.5);
		u0[i].fixup();
	}
	for (int i = N / 2; i < N; i++) {
		u0[i] = Vars(.125, 0, .25);
		u0[i].fixup();
	}
#else
	for (int i = 0; i < N; i++) {
		u0[i] = Vars(.2, .2, 2.6);
		u0[i].fixup();
	}
	for (int i = .4 * N; i < .6 * N; i++) {
		u0[i] = Vars(1, 1, 3);
		u0[i].fixup();
	}
#endif

	const double t_max = 1;
	const int output = 5;

	Main m(N, 
#ifdef SHOCK
		false,
#else
		true,
#endif
		u0, u1, u2);

	std::vector<Alphas *> chain;
	chain.push_back(&rusanov);
//	chain.push_back(&test);
//	chain.push_back(&laxwendroff);
//	chain.push_back(&beamwarming);
//	chain.push_back(&upstream);

	while (m.getTime() < t_max) {
		if (!m.getIter())
			m.doFirstStep(C);
		else
			m.doStep(chain);

		if (!(m.getIter() % output))
		{
			char fn[1024];
			sprintf(fn, "sod%.5d.csv", m.getIter());
			FILE *f = fopen(fn, "w");
			fprintf(f, "rho,u,eps,scheme\n");
			const Array<int> &scnum = m.getSchemeNums();
			for (int j = 0; j < N; j++) {
				Vector v;
				v = u1[j].to_nconserv();
				fprintf(f, "%.10le,%.10le,%.10le,%d\n", v(0), v(1), v(2), scnum[j]);
			}
			fclose(f);
		}
	}

	return 0;
}
