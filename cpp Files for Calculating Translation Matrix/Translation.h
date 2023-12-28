#pragma once
#include "Matrix.h"
#include "Matrix.cpp"
#include <iostream>
#include <string>
#include <fstream>
using namespace std;
void calc_wav_TransR_sp_zaxis_simple(int degree, vector<double>& beta, double d, vector<Matrix<double>>& mat, int samplenum);//ɨƵ
void calc_wav_TransP_sp_zaxis_simple(int degree, vector<double>& beta, double d, vector<Matrix<complex<double>>>& mat, int samplenum);//ɨƵ
