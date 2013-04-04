#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "array.h"
#include "coeff.h"
#include <vector>

template<class V>
struct Solver {
	double cou;
	int N;
	int iter;
	double h;
	double dt;
	double t;
	double tolerance;

	Array<V> &u0;
	Array<V> &u1;
	Array<V> &u2;

	Solver(int n, Array<V> &u0, Array<V> &u1, Array<V> &u2) : 
		N(n),
		u0(u0), u1(u1), u2(u2)
	{
		h = 1.0 / N;
		t = 0;
		iter = 0;
		tolerance = 1e-9;
	}

	double getTime() const { return t; }
	int getIter() const { return iter; }

	virtual double doFirstStep(double C) = 0;
	virtual double doStep(const std::vector<Alphas *> &schemes) = 0;
	virtual ~Solver() { }
};


#endif
