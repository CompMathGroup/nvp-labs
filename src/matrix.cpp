#include <stdio.h>

#include "matrix.h"

/* y += ax */
void saxpy(const double a, const vector_t x, vector_t y) {
	y[0] += a * x[0];
	y[1] += a * x[1];
	y[2] += a * x[2];
}

/* y *= a */
void scale(const double a, vector_t y) {
	y[0] *= a;
	y[1] *= a;
	y[2] *= a;
}

/* y = diag {u - c, u, u + c} */
void lambda(const double c, const double u, vector_t y) {
	y[0] = u;
	y[1] = u - c;
	y[2] = u + c;
}

/* W and iW */
void omega(const double c, const double u, const double gamma, matrix_t w, matrix_t iw) {
	double u2 = u * u;
	double iu = 1 / u;
	double iu2 = iu * iu;
	double c2 = c * c;
	double zm = (gamma - 1) / (c2 + (gamma - 1) * u * (0.5 * u - c));
	double zp = (gamma - 1) / (c2 + (gamma - 1) * u * (0.5 * u + c));

	w[0][0] = 2 * iu2;	w[0][1] = zm;			w[0][2] = zp;
	w[1][0] = 2 * iu;	w[1][1] = zm * (u - c);	w[1][2] = zp * (u + c);
	w[2][0] = 1;		w[2][1] = 1;			w[2][2] = 1;

	double zz = 0.5 / (c * (zp + zm - u2 * zm * zp));

	iw[0][0] = zz * (u2 * (u * (zp - zm) + c * (zm + zp)));
	iw[0][1] = zz * (u2 * (zm - zp));
	iw[0][2] = zz * (-2 * c * u2 * zm * zp);

	iw[1][0] = zz * (( 2 - c * u * zp - u2 * zp) * u);
	iw[1][1] = zz * (-2 + u2 * zp);
	iw[1][2] = zz * (2 * c * zp);

	iw[2][0] = zz * ((-2 - c * u * zm + u2 * zm) * u);
	iw[2][1] = zz * ( 2 - u2 * zm);
	iw[2][2] = zz * (2 * c * zm);
}

/* C = A*B */
void gemm(const matrix_t a, const matrix_t b, matrix_t c) {
#define MUL_ITEM(a, b, i, j, k) (a[i][k] * b[k][j])
#define MUL_SUM(a, b, c, i, j) c[i][j] = \
	MUL_ITEM(a, b, i, j, 0) + \
	MUL_ITEM(a, b, i, j, 1) + \
	MUL_ITEM(a, b, i, j, 2) 

	MUL_SUM(a, b, c, 0, 0);
	MUL_SUM(a, b, c, 0, 1);
	MUL_SUM(a, b, c, 0, 2);
	MUL_SUM(a, b, c, 1, 0);
	MUL_SUM(a, b, c, 1, 1);
	MUL_SUM(a, b, c, 1, 2);
	MUL_SUM(a, b, c, 2, 0);
	MUL_SUM(a, b, c, 2, 1);
	MUL_SUM(a, b, c, 2, 2);

#undef MUL_ITEM
#undef MUL_SUM
}

/* dump A */
void dump(const char *tag, matrix_t a) {
	printf("Matrix dump : %s\n", tag);
	for (int i = 0; i < 3; i++)
		printf("[ %10.6f %10.6f %10.6f ]\n", a[i][0], a[i][1], a[i][2]);
}

/* C = A*diag(v)*B */
void mvmm(const matrix_t a, const vector_t v, const matrix_t b, matrix_t c) {
#define MULV_ITEM(a, b, v, i, j, k) (a[i][k] * v[k] * b[k][j])
#define MULV_SUM(a, b, v, c, i, j) c[i][j] = \
	MULV_ITEM(a, b, v, i, j, 0) + \
	MULV_ITEM(a, b, v, i, j, 1) + \
	MULV_ITEM(a, b, v, i, j, 2) 

	MULV_SUM(a, b, v, c, 0, 0);
	MULV_SUM(a, b, v, c, 0, 1);
	MULV_SUM(a, b, v, c, 0, 2);
	MULV_SUM(a, b, v, c, 1, 0);
	MULV_SUM(a, b, v, c, 1, 1);
	MULV_SUM(a, b, v, c, 1, 2);
	MULV_SUM(a, b, v, c, 2, 0);
	MULV_SUM(a, b, v, c, 2, 1);
	MULV_SUM(a, b, v, c, 2, 2);

#undef MULV_ITEM
#undef MULV_SUM
}

/* Fix illegal variables */
void fix(const double gamma, vector_t u) {
	double rho = u[0];
	if (rho < RHO_MIN)
		rho = RHO_MIN;

	double P = u[1];
	double v = P / rho;
	double K = .5 * P * v;
	double E = u[2];
	double p = (gamma - 1) * (E - K);

	if (p < P_MIN)
		p = P_MIN;

	u[0] = rho;
	u[2] = p / (gamma - 1) + K;
}

/* F */
void f(const double gamma, const vector_t u, vector_t f) {
	fix(gamma, u);
	double rho = u[0];
	double P = u[1];
	double E = u[2];

	double v = P / rho;
	double K = .5 * P * u;
	double p = (gamma - 1) * (E - K);

	f[0] = P;
	f[1] = P * v + p;
	f[2] = (E + p) * v;
}
