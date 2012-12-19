#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <stdint.h>
#include <cmath>

template<class CT>
class Array {
protected:
	CT *cells;
	std::ptrdiff_t ext;
	std::ptrdiff_t n;
public:
	Array(std::ptrdiff_t _n, std::ptrdiff_t _ext) {
		n = _n;
		ext = _ext;
		**this = new CT[_n + 2 * _ext];
	}
	
	std::ptrdiff_t size() const {
		return n;
	}

	CT *& operator*() {
		return cells;
	}
	
	~Array() {
		delete[] cells;
	}

	CT &operator[] (std::ptrdiff_t i) {
		return cells[i + ext];
	}
	
	const CT &operator[] (std::ptrdiff_t i) const {
		return cells[i + ext];
	}

	virtual void actualize() { }
};

template<class CT>
class CycledArray : public Array<CT> {
public:
	CycledArray(std::ptrdiff_t _n, std::ptrdiff_t _ext) : Array<CT>(_n, _ext) {}
	void actualize() {
		std::ptrdiff_t ext = this->ext;
		std::ptrdiff_t n = this->n;
		CT *cells = this->cells + ext;
		for (std::ptrdiff_t i = 1; i <= ext; i++)
			cells[-i] = cells[n - i];
		for (std::ptrdiff_t i = 0; i < ext; i++)
			cells[n + i] = cells[ i];
	}
};

template<class CT>
class ExtrapolatingArray : public Array<CT> {
public:
	ExtrapolatingArray(std::ptrdiff_t _n, std::ptrdiff_t _ext) : Array<CT> (_n, _ext) {}
	void actualize() {
		std::ptrdiff_t ext = this->ext;
		std::ptrdiff_t n = this->n;
		CT *cells = this->cells + ext;
		for (std::ptrdiff_t i = 1; i <= ext; i++)
			cells[-i] = cells[0];
		for (std::ptrdiff_t i = 0; i < ext; i++)
			cells[n + i] = cells[n - 1];
	}
};

template <class T>
void swap(Array<T> &src, Array<T> &dst) {
	T *tmp = *src;
	*src = *dst;
	*dst = tmp;
}
template <class T>
void advance(Array<T> &u0, Array<T> &u1, Array<T> &u2) {
	T *tmp = *u2;
	*u2 = *u0;
	*u0 = *u1;
	*u1 = tmp;
}

class HalfInteger {
private:
	int v;
public:
	explicit HalfInteger(int _v) : v(_v) {}
	HalfInteger operator -() const {
		return HalfInteger(-v);
	}
	friend std::ptrdiff_t operator +(std::ptrdiff_t i, const HalfInteger &o) {
		return i + (o.v >> 1);
	}
	friend std::ptrdiff_t operator -(std::ptrdiff_t i, const HalfInteger &o) {
		return i + ((-o.v) >> 1);
	}
};

#endif
