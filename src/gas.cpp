#include "coeff.h"
#include "main.h"

#include <stdlib.h>
#include <stdio.h>
#include <fenv.h>
#include <cmath>
#include <iostream>

#define NSHOCK
#define TEST

int main(int argc, char **argv) {

	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

	const int N = 100;
	const double C = 0.5;

#ifdef SHOCK
	ExtrapolatingArray<Vars> u0(N, 2), u1(N, 2), u2(N, 2);
#else
	CycledArray<Vars> u0(N, 2), u1(N, 2), u2(N, 2);
#endif

	Rational zero(0, 0, 0, 0, 1, 0, 0, 0); /* 0 */

#ifdef UPSTREAM
	Rational aup0(1, 0, 0, 0, 1, 0, 0, 0); /* 1 */
	Rational a0(-1, 1, 0, 0, 1, 0, 0, 0); /* s - 1 */
	Rational a_2(zero);
#endif

#ifdef LAXWENDROFF
	Rational aup0(1, 0, 0, 0, 1, 0, 0, 0); /* 1 */
	Rational a0(-1, 0, 1, 0, 1, 0, 0, 0); /* s^2 - 1 */
	Rational a_2(zero);
#endif

#ifdef BEAMWORMING
	Rational aup0(1, 0, 0, 0, 1, 0, 0, 0); /* 1 */
	Rational a0(-2, 3, -1, 0, 2, 0, 0, 0); /* .5(2-s)(s-1) */
	Rational a_2(0, 1, -1, 0, 2, 0, 0, 0); /* .5(1-s)s */
#endif

#ifdef RUSANOV
	Rational aup_1(zero);
	Rational aup0(1, 0, 0, 0, 1, 0, 0, 0); /* 1 */
	Rational aup1(zero);

	Rational a_2(0, 1, 0, -1, 6, 0, 0, 0); 
	Rational a_1(0, -6, -3, 3, 6, 0, 0, 0);
	Rational a0(-6, 3, 6, -3, 6, 0, 0, 0); 
	Rational a1(0, 2, -3, 1, 6, 0, 0, 0);
	Rational a2(zero);

	Rational adown_2(zero);
	Rational adown_1(zero);
	Rational adown0(zero);
	Rational adown1(zero);
	Rational adown2(zero);
#endif

#ifdef TEST
	Rational aup_1(0, -12, 6, 0, 12, -18, 0, 0);
	Rational aup0(12, -6, -6, 0, 12, -18, 0, 0); 
	Rational aup1(zero);

	Rational a_2(0, 2, 3, 1, 12, -18, 0, 0);
	Rational a_1(zero);
	Rational a0(-12, 12, 3, -3, 12, -18, 0, 0);
	Rational a1(0, 4, -6, 2, 12, -18, 0, 0);
	Rational a2(zero);

	Rational adown_2(zero);
	Rational adown_1(zero);
	Rational adown0(zero);
	Rational adown1(zero);
	Rational adown2(zero);
#endif

	Alphas alphas(
	             aup_1,   aup0,   aup1,
		a_2,     a_1,     a0,     a1,     a2,
		adown_2, adown_1, adown0, adown1, adown2);
	if (!alphas.sanity())
		exit(1);

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
/*
	for (int i = 0; i < N / 2; i++) {
		u0[i] = Vars(1, 0, 2.5);
		u0[i].fixup();
	}
	for (int i = N / 2; i < N; i++) {
		u0[i] = Vars(.125, 0, .25);
		u0[i].fixup();
	} */
	for (int i = 0; i < N / 2; i++) {
		u0[i] = Vars(.2, .2, 2.6);
//		u0[i] = Vars(1, 1, 3);
		u0[i].fixup();
	}
	for (int i = N / 2; i < N; i++) {
//		u0[i] = Vars(.1, .1, 2.55);
		u0[i] = Vars(1, 1, 3);
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
	const int output = 1;

	Main m(N, 
#ifdef SHOCK
		false,
#else
		true,
#endif
		u0, u1, u2);

	std::vector<Alphas *> chain;
	chain.push_back(&alphas);

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
			fprintf(f, "rho,u,eps\n");
			for (int j = 0; j < N; j++) {
				Vector v;
				v = u1[j].to_nconserv();
				fprintf(f, "%.10le,%.10le,%.10le\n", v(0), v(1), v(2));
			}
			fclose(f);
		}
	}

	return 0;
}
