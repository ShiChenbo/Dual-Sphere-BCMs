#pragma once
#include "Matrix.h"

template<typename Y>
Matrix<complex<double>> inv(Matrix<complex<Y>>& M);//������������󷵻� ��double

template<typename Y>
Matrix<double> inv(Matrix<Y>& M);//������������󷵻� ��double

template<typename Y>
Matrix<double> rref(Matrix<Y>& M, double tol = 0.000001);

template<typename Y>
Matrix<double> rref(Matrix<Y>& M, Matrix<double>& rM, double tol = 0.000001);//rM�����б任����

template<typename Y>
Matrix<complex<double>> rref(Matrix<complex<Y>>& M, double tol = 0.000001);

template<typename Y>
Matrix<complex<double>> rref(Matrix<complex<Y>>& M, Matrix<complex<double>>& rM, double tol = 0.000001);

template<typename Y>
int rank_m(Matrix<Y> M, double  tol = 0.000001);

template<typename Y>
int rank_m(Matrix<complex<Y>> M, double  tol = 0.000001);

template<typename Y>
Matrix<double> null(Matrix<Y>& M, double tol = 0.000001);

template<typename Y>
Matrix<complex<double>> null(Matrix<complex<Y>>& M, double tol = 0.000001);

template<typename Y>
Matrix<Y> les_solve(Matrix<Y>& A, Matrix<Y>& b, double tol = 1e-10);//Ax=b�Ľ⣬��˹�����

