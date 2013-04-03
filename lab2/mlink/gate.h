#ifndef __GATE_H__
#define __GATE_H__

double InitializeSolver(
	const char *solv, 
	double gamma,
	const char *prob, 
	int n, double C,
	double rL, double uL, double eL,
	double rR, double uR, double eR);
int AddScheme(double *coeff, int sz);
void DoSteps(int count);

#endif
