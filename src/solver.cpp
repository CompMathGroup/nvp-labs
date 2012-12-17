#include "solver.h"

#include <iostream>

/*
In-place solve
*/
void TridiagonalSolve(
	Array<Matrix> &A, Array<Matrix> &B, Array<Matrix> &C, Array<Vector> &d) {
	int N = A.size() - 1;
	Matrix T;
	T = -1. * B[0].inverse();
	B[0] = T * C[0];
	d[0] = -1. * (T * d[0]);
	for (int k = 1; k <= N; k++) {
		T = -1. * (A[k] * B[k-1] + B[k]).inverse();
		B[k] = T * C[k];
		d[k] = T * (A[k] * d[k-1] - d[k]);
	}
	for (int k = N-1; k >= 0; k--)
		d[k] += B[k] * d[k+1];
}

/*
Debug variant
*/
void TridiagonalSolveDebug(
	Array<Matrix> &A, Array<Matrix> &B, Array<Matrix> &C, Array<Vector> &d) {
	int N = A.size() - 1;

	Array<Matrix> P(N+1, 0);
	Array<Vector> q(N+1, 0);

	Matrix T;
	T = -1. * B[0].inverse();
	P[0] = T * C[0];
	q[0] = -1. * (T * d[0]);
	for (int k = 1; k <= N; k++) {
		T = -1. * (A[k] * P[k-1] + B[k]).inverse();
		P[k] = T * C[k];
		q[k] = T * (A[k] * q[k-1] - d[k]);
	}
	for (int k = N-1; k >= 0; k--)
		q[k] += P[k] * q[k+1];

	Vector res;
	res = B[0] * q[0] + C[0] * q[1] - d[0];
	if (res.norm() > 1e-10 * d[0].norm())
		std::cerr << "problem in 0 eq " 
			<< res.norm() << " vs " << d[0].norm()
			<< std::endl;
	res = B[N] * q[N] + A[N] * q[N-1] - d[N];
	if (res.norm() > 1e-10 * d[N].norm())
		std::cerr << "problem in N eq " 
			<< res.norm() << " vs " << d[N].norm()
			<< std::endl;
	for (int i = 1; i < N; i++) {
		res = A[i] * q[i-1] + B[i]*q[i] + C[i]*q[i+1] - d[i];
		if (res.norm() > 1e-10 * d[i].norm())
			std::cerr << "problem in " << i << " eq " 
			<< res.norm() << " vs " << d[i].norm()
			<< std::endl;
	}

	for (int i = 0; i <= N; i++)
		d[i] = q[i];
}

void EliminateExtrapolated(
	Array<Matrix> &A, Array<Matrix> &B, Array<Matrix> &C, Array<Vector> &d) {
	int N = A.size() - 1;
	B[0] += A[0];
	A[0] -= A[0];

	B[N] += C[N];
	C[N] -= C[N];
}

void EliminateCycled(
	Array<Matrix> &A, Array<Matrix> &B, Array<Matrix> &C, Array<Vector> &d) {
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

	A[0] *= 0;
	B[0] *= 0;
	C[0] *= 0;
	B[0] += 1;
	d[0] = x0;

	A[N] *= 0;
	B[N] *= 0;
	C[N] *= 0;
	B[N] += 1;
	d[N] = xN;
}
