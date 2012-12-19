#ifndef __MAIN_H__
#define __MAIN_H__

#include "array.h"
#include "gd.h"
#include "matrix.h"
#include "coeff.h"
#include <vector>

class Main {
	Array<Vars> &u0;
	Array<Vars> &u1;
	Array<Vars> &u2;

	Array<Vars> f;

	Array<Vars> ulower;
	Array<Vars> uupper;

	Array<Vars> F1;

	Array<Vector> lam0;
	Array<Vector> lam1;
	Array<Vector> lam2;

	Array<Matrix> W0, iW0;
	Array<Matrix> W1, iW1;
	Array<Matrix> W2, iW2;

	Array<Matrix> B1;
	Array<Matrix> B2;
	Array<Matrix> B3;
	Array<Vector> rhs;
	Array<Vector> sol;

	Array<int> schemenum;

	int N;
	double cou;
	double h;
	double dt;
	int iter;
	double t;
	double Cmax;
	bool cycled;
private:
	void predict();
	void correct(const Alphas &scheme, const std::vector<int> &where);
	bool checkMono(std::vector<int> &where);
public:
	Main(int n, bool cycled, 
		Array<Vars> &u0, Array<Vars> &u1, Array<Vars> &u2) : 
		u0(u0), u1(u1), u2(u2),
		f(n-1, 1),
		ulower(n, 0), uupper(n, 0),
		F1(n, 2),
		lam0(n, 2),	lam1(n, 2),	lam2(n, 2),
		W0(n, 2), iW0(n, 2),
		W1(n, 2), iW1(n, 2),
		W2(n, 2), iW2(n, 2),
		B1(n, 0), B2(n, 0),	B3(n, 0), 
		rhs(n, 0), sol(n, 0),
		schemenum(n, 0),
		N(n),
		cycled(cycled)
	{
		h = 1.0 / N;
		t = 0;
		iter = 0;
		Cmax = 0;
	}

	double getTime() const { return t; }
	int getIter() const { return iter; }
	const Array<int> &getSchemeNums() const { return schemenum; }
	void doFirstStep(double C);
	double doStep(const std::vector<Alphas *> schemes);
};

#endif
