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
	
	virtual ~Array() {
		delete[] cells;
	}

	CT &operator[] (std::ptrdiff_t i) {
		return cells[i + ext];
	}
	
	const CT &operator[] (std::ptrdiff_t i) const {
		return cells[i + ext];
	}

	virtual void actualize(double t) { }
	virtual void fillAt(double t) { }
	virtual CT getExtrapolated(const int i, double t) { return 0. / 0.; }
};

template<class CT, class E>
class ExtrapolatingArray : public Array<CT> {
	E extr;
public:
	ExtrapolatingArray(std::ptrdiff_t _n, std::ptrdiff_t _ext, E _extr) : Array<CT> (_n, _ext), extr(_extr) {}
	void actualize(double t) {
		std::ptrdiff_t ext = this->ext;
		std::ptrdiff_t n = this->n;
		CT *cells = this->cells + ext;
		for (std::ptrdiff_t i = 1; i <= ext; i++)
			cells[-i] = extr(-i, t);
		for (std::ptrdiff_t i = 0; i < ext; i++)
			cells[n + i] = extr(n + i, t);
	}
	void fillAt(double t) {
		std::ptrdiff_t ext = this->ext;
		std::ptrdiff_t n = this->n;
		CT *cells = this->cells + ext;
		for (std::ptrdiff_t i = -ext; i < n + ext; i++)
			cells[i] = extr(i, t);
	}
	CT getExtrapolated(const int i, double t) {
		return extr(i, t);
	}
	virtual ~ExtrapolatingArray() { }
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

/* This may not be portable. Valid if (-1) >> 1 === -1, 1 >> 1 === 0 */
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
