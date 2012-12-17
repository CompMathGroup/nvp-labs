#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "array.h"
#include "matrix.h"

void TridiagonalSolve(
	Array<Matrix> &A, Array<Matrix> &B, Array<Matrix> &C, Array<Vector> &d);

void TridiagonalSolveDebug(
	Array<Matrix> &A, Array<Matrix> &B, Array<Matrix> &C, Array<Vector> &d);

void EliminateExtrapolated(
	Array<Matrix> &A, Array<Matrix> &B, Array<Matrix> &C, Array<Vector> &d);

void EliminateCycled(
	Array<Matrix> &A, Array<Matrix> &B, Array<Matrix> &C, Array<Vector> &d);

#endif
