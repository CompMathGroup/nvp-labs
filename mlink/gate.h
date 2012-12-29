#ifndef __GATE_H__
#define __GATE_H__

#ifdef __cplusplus
extern "C" {
#endif

double InitializeSolver(
	const char *solv, const char *prob, 
	int n, double C,
	double rL, double uL, double eL,
	double rR, double uR, double eR);
int AddScheme(double *coeff, int sz);
void DoSteps(int count);

#ifdef __cplusplus
}
#endif

#endif
