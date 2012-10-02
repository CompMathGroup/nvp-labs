#include "grid.h"

Grid::Grid(const int _n, const int _over) : n(_n), over(_over) {
	ucurr = new CycledArray(n, over);
	uprev = new CycledArray(n, over);
	unext = new CycledArray(n, over);
	upred = new CycledArray(n, over);
}

void Grid::swap() {
	uswap = uprev;
	uprev = ucurr;
	ucurr = unext;
	unext = uswap;
}

void Grid::validate() {
	uprev->validate();
	u    ->validate();
	unext->validate();
}

void Grid::iterate() {
	for (int i = 0; i < n; i++) {
		
	}
	overwrite();
}
