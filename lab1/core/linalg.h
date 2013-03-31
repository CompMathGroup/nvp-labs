#ifndef __LINALG_H__
#define __LINALG_H__

#include "array.h"
#include "matrix.h"

template <typename M, typename V>
void TridiagonalSolve(
	const Array<M> &A, const Array<M> &B, const Array<M> &C, const Array<V> &d, 
	Array<V> &x, bool cycled);

template <typename M, typename V>
void EliminateCycled(
	const Array<M> &A, const Array<M> &B, const Array<M> &C, const Array<V> &d, 
	Array<V> &x);

#include "linalg.cpp"

#endif
