#ifndef __MATRIX_H__
#define __MATRIX_H__

typedef double vector_t[3];
typedef double matrix_t[3][3];

/* y += ax */
void saxpy(const double a, const vector_t x, vector_t y);

/* y *= a */
void scale(const double a, vector_t y);

/* y = diag {u - c, u, u + c} */
void lambda(const double c, const double u, vector_t y);

/* W and iW */
void omega(const double c, const double u, const double gamma, matrix_t w, matrix_t iw);

/* F */
void f(const double gamma, const vector_t u, vector_t f);

/* Fix illegal values */
void fix(const double gamma, vector_t u);

/* C = A*B */
void gemm(const matrix_t a, const matrix_t b, matrix_t c);

/* C = A*diag(v)*B */
void mvmm(const matrix_t a, const vector_t v, const matrix_t b, matrix_t c);

/* dump A */
void dump(const char *tag, matrix_t a);

#endif
