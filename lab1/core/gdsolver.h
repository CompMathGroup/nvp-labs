#ifndef __GDSOLVER_H__
#define __GDSOLVER_H__

#include "gd.h"
#include "matrix.h"
#include "solver.h"

class GDSolver : public Solver<Vars> {
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

	double Cmax;
	bool cycled;
private:
	void predict();
	void eno();
	void correct(const Alphas &scheme, const std::vector<int> &where);
	bool checkMono(std::vector<int> &where);
public:
	GDSolver(int n, bool cycled, 
		Array<Vars> &u0, Array<Vars> &u1, Array<Vars> &u2) 
	:
		Solver(n, u0, u1, u2),
		f(n-1, 1),
		ulower(n, 0), uupper(n, 0),
		F1(n, 2),
		lam0(n, 2),	lam1(n, 2),	lam2(n, 2),
		W0(n, 2), iW0(n, 2),
		W1(n, 2), iW1(n, 2),
		W2(n, 2), iW2(n, 2),
		B1(n, 0), B2(n, 0),	B3(n, 0), 
		rhs(n, 0), sol(n, 0),
		Cmax(0),
		cycled(cycled)
	{ }

	virtual double doFirstStep(double C);
	virtual double doStep(const std::vector<Alphas *> &schemes);
	virtual ~GDSolver() { }
};

#endif
