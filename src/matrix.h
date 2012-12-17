#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <cstddef>

struct Vector {
	double _v[3];
	Vector() {}
	Vector(double x, double y, double z) {
		_v[0] = x;
		_v[1] = y;
		_v[2] = z;
	}
	const double &operator ()(ptrdiff_t i) const {
		return _v[i];
	}
	double &operator ()(ptrdiff_t i) {
		return _v[i];
	}
	Vector operator +(const Vector &y) const {
		const Vector &x = *this;
		Vector z;
		for (int i = 0; i < 3; i++)
			z(i) = x(i) + y(i);
		return z;
	}
	Vector operator -(const Vector &y) const {
		const Vector &x = *this;
		Vector z;
		for (int i = 0; i < 3; i++)
			z(i) = x(i) - y(i);
		return z;
	}
	double norm() const {
		double ret = 0;
		for (int i = 0; i < 3; i++) {
			double s = std::fabs(_v[i]);
			if (s > ret)
				ret = s;
		}
		return ret;
	}
	friend Vector operator *(const double m, const Vector &y) {
		Vector z;
		for (int i = 0; i < 3; i++)
			z(i) = m * y(i);
		return z;
	}
	Vector &operator *=(const double m) {
		for (int i = 0; i < 3; i++)
			_v[i] *= m;
		return *this;
	}
	Vector &operator /=(const double m) {
		for (int i = 0; i < 3; i++)
			_v[i] /= m;
		return *this;
	}
	Vector &operator +=(const Vector &o) {
		for (int i = 0; i < 3; i++)
			_v[i] += o._v[i];
		return *this;
	}
	Vector &operator -=(const Vector &o) {
		for (int i = 0; i < 3; i++)
			_v[i] -= o._v[i];
		return *this;
	}
}; 

struct Matrix {
	double _m[3][3];
	Matrix() {
	}
	const double &operator ()(ptrdiff_t i, ptrdiff_t j) const {
		return _m[i][j];
	}
	double &operator ()(ptrdiff_t i, ptrdiff_t j) {
		return _m[i][j];
	}
	Matrix operator *(const Matrix &b) const {
		Matrix c;
		const Matrix &a = *this;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				c(i, j) = a(i, 0) * b(0, j) + a(i, 1) * b(1, j) + a(i, 2) * b(2, j);
		return c;
	}
	Matrix operator +(const Matrix &b) const {
		Matrix c;
		const Matrix &a = *this;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				c(i, j) = a(i, j) + b(i, j);
		return c;
	}
	Matrix operator -(const Matrix &b) const {
		Matrix c;
		const Matrix &a = *this;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				c(i, j) = a(i, j) - b(i, j);
		return c;
	}
	Matrix &operator +=(const double v) {
		_m[0][0] += v;
		_m[1][1] += v;
		_m[2][2] += v;
		return *this;
	}
	Matrix &operator +=(const Matrix &b) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j <3; j++)
				_m[i][j] += b._m[i][j];
		return *this;
	}
	Matrix &operator -=(const Matrix &b) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j <3; j++)
				_m[i][j] -= b._m[i][j];
		return *this;
	}
	Matrix &operator *=(const double v) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j <3; j++)
				_m[i][j] *= v;
		return *this;
	}
	Vector operator *(const Vector &x) const {
		Vector v;
		const Matrix &a = *this;
		for (int i = 0; i < 3; i++) 
			v(i) = a(i, 0) * x(0) + a(i, 1) * x(1) + a(i, 2) * x(2);
		return v;
	}
	Vector solve(const Vector &b) const {
		Matrix iA = this->inverse();
		return iA * b;
	}
	double norm() const {
		double ret = 0;
		for (int i = 0; i < 3; i++) {
			double s = 
				std::fabs(_m[i][0]) + 
				std::fabs(_m[i][1]) + 
				std::fabs(_m[i][2]);
			if (s > ret)
				ret = s;
		}
		return ret;
	}
	double cond() const {
		return this->norm() * this->inverse().norm();
	}
	template <class Functor>
	Matrix(const Matrix &W, const Vector &e, const Matrix &iW, const Functor &func = Functor()) {
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++) {
				double Aij = 0;
				for (int k = 0; k < 3; k++)
					Aij += func(e(k)) * W(i, k) * iW(k, j);
				_m[i][j] = Aij;
			}
	}
	Matrix inverse() const {
		double a = _m[0][0];
		double b = _m[0][1];
		double c = _m[0][2];

		double d = _m[1][0];
		double e = _m[1][1];
		double f = _m[1][2];

		double g = _m[2][0];
		double h = _m[2][1];
		double i = _m[2][2];

		Matrix m;
		m(0, 0) = e*i - f*h;
		m(0, 1) = c*h - b*i;
		m(0, 2) = b*f - c*e;

		m(1, 0) = f*g - d*i;
		m(1, 1) = a*i - c*g;
		m(1, 2) = c*d - a*f;

		m(2, 0) = d*h - e*g;
		m(2, 1) = b*g - a*h;
		m(2, 2) = a*e - b*d;

		double D = m(0, 0) * a + m(0, 1) * d + m(0, 2) * g;

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				m(i, j) /= D;

		return m;
	}
};

inline
Matrix operator *(double v, const Matrix &m) {
	Matrix ret;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			ret(i, j) = v * m(i, j);
	return ret;
}

#endif
