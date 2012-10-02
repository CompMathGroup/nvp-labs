#ifndef __GRID_H__
#define __GRID_H__

#include "matrix.h"

class Grid {
	const int n;
	const int over;
	CycledArray *u, *unext, *uprev, *uswap;
	void swap();
public:
	Grid(const int _n, const int _over);
	void iterate();
	void validate();
};

#endif
