#include "Matrix.h"
#include <iostream>
#include "Exteption.h"
#include <complex>
#include <stdio.h>
/* ¬ариант 2. ћатрицы
 ласс матриц вещественных чисел произвольного размера.
 ак минимум, предоставить:
Х конструктор с параметрами : размеры матрицы и значение дл€ заполнени€;**
Х оператор() дл€ чтени€ / записи элемента матрицы по указанным индексам;*
Х операторы сложени€ и вычитани€ матриц;**
Х оператор умножени€ матриц;**
Х оператор умножени€ матрицы на скал€р(обеспечить коммутативность);**
Х оператор делени€ матрицы на скал€р;**
Х вычисление следа матрицы.**
«јƒј„ј : ѕривести заданную квадратную матрицу ј к нижнетреугольному виду. */

using namespace std;

template <typename T>
Matrix<T>::Matrix() {
	m = 0;
	n = 0;
	M = NULL;
}

template <typename T>
Matrix<T>::Matrix(int m_, int n_) {
	this->n = n_;
	this->m = m_;

	M = (T**) new T * [m];
	for (int i = 0; i < m; i++) {
		M[i] = (T*) new T[n];
	}
}

template <typename T>
Matrix<T>::Matrix(int m_, int n_, T value) {
	this->n = n_;
	this->m = m_;

	M = (T**) new T * [m];
	for (int i = 0; i < m; i++) {
		M[i] = (T*) new T[n];
	}

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++) {
			this->M[i][j] = value;
		}
}


template <typename T>
Matrix<T>::Matrix(const Matrix& M_) {
	m = M_.m;
	n = M_.n;
	M = (T**) new T * [m];
	for (int i = 0; i < m; i++) {
		M[i] = (T*) new T[n];
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			M[i][j] = M_.M[i][j];
		}
	}
}


template <typename T>
int Matrix<T>::GetM() {
	return m;
}


template <typename T>
int Matrix<T>::GetN() {
	return n;
}


template <typename T>
void Matrix<T>::Print() {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << M[i][j] << '\t';
		}
	}
}


template <typename T>
T& Matrix<T>::operator ()(int i, int j)
{
	if ((i >= 0) || (i < m) || (j >= 0) || (j < n))
		return M[i][j];
	else throw EInvalidIndex();
}


template <typename T>
Matrix<T>& Matrix<T>::operator ()(int i, int j, T value)
{
	if ((i >= 0) && (i <= m) && (j >= 0) && (j <= n))
		this->M[i][j] = value;
	else throw EInvalidIndex();
}


template <typename T>
Matrix<T> Matrix<T>::operator + (const Matrix& B) {

	if (n != B.n || m != B.m) throw EInvalidSize();
	Matrix tmp(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			tmp.M[i][j] = M[i][j] + B.M[i][j];
		}
	}
	return tmp;
}


template <typename T>
Matrix<T> Matrix<T>::operator - (const Matrix& B) {
	if (n != B.n || m != B.m) throw EInvalidSize();
	int m_ = m;
	int n_ = n;
	Matrix tmp(m_, n_);
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < n_; j++) {
			tmp.M[i][j] = M[i][j] - B.M[i][j];
		}
	}
	return tmp;
}


template <typename T>
Matrix<T> Matrix<T>::operator * (const Matrix& B) {
	if (m != B.n) throw EInvalidMull();
	int m_ = GetM();
	int n_ = B.n;
	Matrix tmp(m, n);
	for (int i = 0; i < m_; i++)
	{
		for (int j = 0; j < m_; j++)
		{
			tmp.M[i][j] = 0;
			for (int k = 0; k < n_; k++)
			{
				tmp.M[i][j] += M[i][k] * B.M[k][j];
			}
		}
	}
	return tmp;
}


template <typename T>
Matrix<T> Matrix<T>::operator * (const int a) {
	Matrix tmp(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			tmp.M[i][j] = M[i][j] * T(a);
		}
	}
	return tmp;
}


template <typename T>
Matrix<T> Matrix<T>::operator / (const int a) {
	if (T(a) == T(0)) throw EdivisionZero();
	Matrix tmp(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			tmp.M[i][j] = M[i][j] / T(a);
		}
	}
	return tmp;
}


template <typename T>
T Matrix<T>::Trace() {
	if (n != m)  throw EInvalidQuadrate();
	T sum = 0;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				sum += M[i][j];
			}
		}
	}
	return sum;
}


template <typename T>
Matrix<T> Matrix<T>::Triangular() {
	if (n != m)  throw EInvalidQuadrate();
	T r;
	for (int i = 0; i < n - 1; i++) {
		for (int j = i + 1; j < n; j++) {
			r = M[j][i] / M[i][i];
			for (int k = 0; k < n + 1; k++) {
				M[j][k] = M[j][k] - r * M[i][k];
			}
		}
	}
	return *this;
}


template <typename T>
void Matrix<T>::Transpose() {
	if (n != m) throw EInvalidQuadrate();
	T s;
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++) {
			s = M[i][j];
			M[i][j] = M[j][i];
			M[j][i] = s;
		}
}


template class Matrix <int>;
template class Matrix <float>;
template class Matrix <double>;
template class Matrix <std::complex<float>>;
template class Matrix <std::complex<double>>;


