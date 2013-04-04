#ifndef __GDSOLVER_H__
#define __GDSOLVER_H__

#include "solver.h"

class HeatSolver : public Solver<double> {
	const double k;
	Array<double> f;

	Array<double> ulower;
	Array<double> uupper;

	Array<double> F1;

	Array<double> kap0;
	Array<double> kap1;
	Array<double> kap2;

	Array<double> B1;
	Array<double> B2;
	Array<double> B3;

	Array<double> rhs;
	Array<double> sol;

	double Cmax;
private:
	void predict();
	void correct(const Alphas &scheme);
	double kappa(double u);
public:
	HeatSolver(int n, double _k, 
		Array<double> &u0, Array<double> &u1, Array<double> &u2) 
	:
		Solver(n, u0, u1, u2),
		k(_k),
		f(n - 1, 1),
		ulower(n, 0), uupper(n, 0),
		F1(n, 2),
		kap0(n-1, 2), kap1(n-1, 2), kap2(n-1, 2),
		B1(n-1, 0), B2(n-1, 0), B3(n-1, 0),
		rhs(n-1, 0), sol(n-1, 0),
		Cmax(0)
	{ }

	virtual double doFirstStep(double C);
	virtual double doStep(const std::vector<Alphas *> &schemes);
	virtual ~HeatSolver() { }
};

#endif
