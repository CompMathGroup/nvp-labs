#ifndef __GATE_H__
#define __GATE_H__

double InitializeSolver(
	const char *solv, 
	double k, double c,
	const char *prob, 
	int n, double Cou);
int AddScheme(double *coeff, int sz);
void DoSteps(int count);

#endif
