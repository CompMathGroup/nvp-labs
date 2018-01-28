#ifndef __ADVSOLVER_H__
#define __ADVSOLVER_H__

#include "solver.h"

class AdvectionSolver : public Solver<double> {
	Array<double> f;

	Array<double> ulower;
	Array<double> uupper;

	Array<double> B1;
	Array<double> B2;
	Array<double> B3;

	Array<double> rhs;
	Array<double> sol;
private:
	void predict();
	void eno();
	void step(const Alphas &scheme, const std::vector<int> &where);
	bool checkMono(std::vector<int> &where);
public:
	AdvectionSolver(int n, 
		Array<double> &u0, Array<double> &u1, Array<double> &u2) 
	:
		Solver(n, u0, u1, u2),
		f(n-1, 1),
		ulower(n, 0), uupper(n, 0),
		B1(n, 0), B2(n, 0),	B3(n, 0), 
		rhs(n, 0), sol(n, 0)
	{ }

	virtual double doFirstStep(double C);
	virtual double doStep(const std::vector<Alphas *> &schemes);
	virtual ~AdvectionSolver() { }
};

#endif
