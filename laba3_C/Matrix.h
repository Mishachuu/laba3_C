#pragma once
/* ������� 2. �������
����� ������ ������������ ����� ������������� �������.
��� �������, ������������:
� ����������� � ����������� : ������� ������� � �������� ��� ����������;**
� ��������() ��� ������ / ������ �������� ������� �� ��������� ��������;
� ��������� �������� � ��������� ������;
� �������� ��������� ������;
� �������� ��������� ������� �� ������(���������� ���������������);
� �������� ������� ������� �� ������;
� ���������� ����� �������.
������ : �������� �������� ���������� ������� � � ����������������� ����. */
#include <iostream>
#include <complex>
#include <vector>
using namespace std;

template <typename T>
class Matrix {
private:
	int n = 0;
	int m = 0;
	vector<vector<T>> M;
	auto begin() { return M.begin(); }
	auto end() { return M.end(); }
public:
	Matrix<T>();
	Matrix<T>(int m_, int n_);
	Matrix<T>(int m_, int n_, T value);
	Matrix<T>(const Matrix& M_);
	void Print();
	int GetM();
	int GetN();
	T& operator ()(int i, int j);
	Matrix& operator ()(int i, int j, T value);
	Matrix<T> operator + (const Matrix<T>& B);
	Matrix<T> operator - (const Matrix<T>& B);
	Matrix<T> operator * (const Matrix<T>& B);
	Matrix<T> operator * (const int a);
	Matrix<T> operator / (const int a);
	T Trace();
	Matrix<T> Triangular();
	void Transpose();
	auto cbegin() const { return M.cbegin(); }
	auto cend() const { return M.cend(); }

	friend std::ostream& operator << (std::ostream& s, const Matrix<T>& matrix) {
		for (auto row:matrix.M) {
			for (auto inrow:row)
				s << inrow << " ";
			s << "\n";
		}
		return s;
	}
};