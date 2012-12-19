#include "solver.h"

#include <iostream>

/*
In-place solve
*/
void TridiagonalSolve(
	const Array<Matrix> &A, const Array<Matrix> &B, const Array<Matrix> &C, const Array<Vector> &d, 
	Array<Vector> &x, bool cycled) 
{
	int N = A.size() - 1;
	Matrix T;
	Array<Matrix> P(N+1, 0);
	if (cycled) {
		T = -1. * B[0].inverse();
		P[0] = T * C[0];
		x[0] = T * (A[0] * x[N] - d[0]);
	} else {
		T = -1. * (A[0] + B[0]).inverse();
		P[0] = T * C[0];
		x[0] = -1. * (T * d[0]);
	}
	for (int k = 1; k < N; k++) {
		T = -1. * (A[k] * P[k-1] + B[k]).inverse();
		P[k] = T * C[k];
		x[k] = T * (A[k] * x[k-1] - d[k]);
	}
	if (!cycled) {
		T = (A[N] * P[N-1] + B[N] + C[N]).inverse();
		x[N] = T * (d[N] - A[N] * x[N-1]);
	}
	for (int k = N-1; k >= 0; k--)
		x[k] += P[k] * x[k+1];
}

void EliminateCycled(
	const Array<Matrix> &A, const Array<Matrix> &B, const Array<Matrix> &C, const Array<Vector> &d, 
	Array<Vector> &x) 
{
	int N = A.size() - 1;
	Matrix P;
	Matrix T;
	Matrix R;
	Matrix P_N0;
	Matrix P_0N;
	Vector q_N0;
	Vector q_0N;
	Vector q;

	T = -1. * (B[0].inverse());
	P = T * C[0];
	R = T * A[0];
	q = -1. * (T * d[0]);

	for (int k = 1; k < N; k++) {
		T = -1. * (A[k] * P + B[k]).inverse();
		P = T * C[k];
		R = T * A[k] * R;
		q = T * (A[k] * q - d[k]);
	}

	T = -1. * (A[N] * (P + R) + B[N]).inverse();
	P = T * C[N];
	q = T * (A[N] * q - d[N]);

	/* x[N] = P x[0] + q */
	P_N0 = P;
	q_N0 = q;

	T = -1. * (B[N].inverse());
	P = T * A[N];
	R = T * C[N];
	q = -1. * (T * d[N]);
	for (int k = N-1; k > 0; k--) {
		T = -1. * (C[k] * P + B[k]).inverse();
		P = T * A[k];
		R = T * C[k] * R;
		q = T * (C[k] * q - d[k]);
	}
	
	T = -1. * (C[0] * (P + R) + B[0]).inverse();
	P = T * A[0];
	q = T * (C[0] * q - d[0]);

	/* x[0] = P x[N] + q */

	P_0N = P;
	q_0N = q;

	Matrix Q = P_0N * P_N0;
	Q += -1;
	Vector qq = q_0N + P_0N * q_N0;
	Vector x0 = Q.solve(-1 * qq);

	Q = P_N0 * P_0N;
	Q += -1;
	qq = q_N0 + P_N0 * q_0N;
	Vector xN = Q.solve(-1 * qq);

	x[0] = x0;
	x[N] = xN;
}
