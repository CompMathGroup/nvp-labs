#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "array.h"
#include "matrix.h"

void TridiagonalSolve(
	const Array<Matrix> &A, const Array<Matrix> &B, const Array<Matrix> &C, const Array<Vector> &d, 
	Array<Vector> &x, bool cycled);

void EliminateCycled(
	const Array<Matrix> &A, const Array<Matrix> &B, const Array<Matrix> &C, const Array<Vector> &d, 
	Array<Vector> &x);

#endif
