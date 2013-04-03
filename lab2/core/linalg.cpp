#ifndef __LINALG_H__
#error "This file sould be inline-compiled"
#endif

/* With this single line templates also work for <double, double> */
inline double inverse(double x) { return 1. / x; }

template <typename M, typename V>
void TridiagonalSolve(
	const Array<M> &A, const Array<M> &B, const Array<M> &C, const Array<V> &d, 
	Array<V> &x, bool cycled) 
{
	int N = A.size() - 1;
	M T;
	Array<M> P(N+1, 0);
	if (cycled) {
		T = -1. * inverse(B[0]);
		P[0] = T * C[0];
		x[0] = T * (A[0] * x[N] - d[0]);
	} else {
		T = -1. * inverse(A[0] + B[0]);
		P[0] = T * C[0];
		x[0] = -1. * (T * d[0]);
	}
	for (int k = 1; k < N; k++) {
		T = -1. * inverse(A[k] * P[k-1] + B[k]);
		P[k] = T * C[k];
		x[k] = T * (A[k] * x[k-1] - d[k]);
	}
	if (!cycled) {
		T = inverse(A[N] * P[N-1] + B[N] + C[N]);
		x[N] = T * (d[N] - A[N] * x[N-1]);
	}
	for (int k = N-1; k >= 0; k--)
		x[k] += P[k] * x[k+1];
}

template <typename M, typename V>
void EliminateCycled(
	const Array<M> &A, const Array<M> &B, const Array<M> &C, const Array<V> &d, 
	Array<V> &x) 
{
	int N = A.size() - 1;
	M P;
	M T;
	M R;
	M P_N0;
	M P_0N;
	V q_N0;
	V q_0N;
	V q;

	T = -1. * inverse(B[0]);
	P = T * C[0];
	R = T * A[0];
	q = -1. * (T * d[0]);

	for (int k = 1; k < N; k++) {
		T = -1. * inverse(A[k] * P + B[k]);
		P = T * C[k];
		R = T * A[k] * R;
		q = T * (A[k] * q - d[k]);
	}

	T = -1. * inverse(A[N] * (P + R) + B[N]);
	P = T * C[N];
	q = T * (A[N] * q - d[N]);

	/* x[N] = P x[0] + q */
	P_N0 = P;
	q_N0 = q;

	T = -1. * inverse(B[N]);
	P = T * A[N];
	R = T * C[N];
	q = -1. * (T * d[N]);
	for (int k = N-1; k > 0; k--) {
		T = -1. * inverse(C[k] * P + B[k]);
		P = T * A[k];
		R = T * C[k] * R;
		q = T * (C[k] * q - d[k]);
	}
	
	T = -1. * inverse(C[0] * (P + R) + B[0]);
	P = T * A[0];
	q = T * (C[0] * q - d[0]);

	/* x[0] = P x[N] + q */

	P_0N = P;
	q_0N = q;

	M Q = P_0N * P_N0;
	Q += -1;
	V qq = q_0N + P_0N * q_N0;
	V x0 = -1 * inverse(Q) * qq;

	Q = P_N0 * P_0N;
	Q += -1;
	qq = q_N0 + P_N0 * q_0N;
	V xN = -1 * inverse(Q) * qq;

	x[0] = x0;
	x[N] = xN;
}
