#ifndef __CYCLEDARRAY_H__
#define __CYCLEDARRAY_H__

class CycledArray {
	const int n, over;
	vector_t * const back;
public:
	vector_t * const arr;
	CycledArray(const int _n, const int _over);
	~CycledArray();
	void validate();
};

#endif
