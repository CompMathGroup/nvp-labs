#include "cycledarray.h"

CycledArray::CycledArray(const int _n, const int _over): n(_n), over(_over) {
	back = new vector_t[2 * over + n];
	arr = back + over;
}

CycledArray::~CycledArray() {
	delete[] back;
}

CycledArray::validate() {
	for (int i = 0; i < over; i++)
		arr[n + i] = arr[ i];
	for (int i = 1; i <= over; i++)
		arr[-i] = arr[n - i];
}
