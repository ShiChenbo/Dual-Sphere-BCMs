#include "MathBasic.h"

std::vector<double> MathBasic::Quad::GetLegendreRoots(int n, double precis, int maxIter)
{
	std::vector<double> xOvec(n);
	for (int i = 1; i <= n; i++)
	{
		xOvec[i - 1] = (1 - 1.0 / (8 * n * n) + 1.0 / (8 * n * n * n)) * cos(pi * (4 * i - 1) / (4 * n + 2));
	}
	int it = 0;
	double* yOvec = new double[n];//初始判断条件
	double tmpPrecis = 10000;
	while (tmpPrecis > precis && it < maxIter)
	{
		double a = 0;
		for (int i = 0; i < n; i++)
		{
			xOvec[i] -= legendre(n, xOvec[i]) / DerivativeLegendre(n, xOvec[i]);
			a = fmax(a, legendre(n, xOvec[i]));
		}
		it++;
		tmpPrecis = a;
	}
	if (it==1000)
	{
		std::cout << "Lengendre 根未收敛\n";
	}
	return std::move(xOvec);
}

double MathBasic::Quad::DerivativeLegendre(int n, double x)
{
	if (n == 0)
	{
		return 0;
	}
	return n * (x * legendre(n, x) - legendre(n - 1, x)) / (x * x - 1);
}

void MathBasic::Quad::GenerateGuassIntegInfo(int n, vector<double>& Intpoints, vector<double>& weights)
{
	Intpoints = GetLegendreRoots(n);
	weights.resize(n);
	for (int i = 0; i < n; i++)
	{
		weights[i] = 2 / ((1 - Intpoints[i] * Intpoints[i]) * pow(DerivativeLegendre(n, Intpoints[i]), 2));
	}
}

void MathBasic::vec_sph_wav::gen_vec_sph_uv1(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, double* x, double* y, double* z, std::complex<double>* res_x, std::complex<double>* res_y, std::complex<double>* res_z, int ndatasize /*= 1*/)
{
	//坐标的笛卡尔-球转换
	double* r = new double[ndatasize];
	double* theta = new double[ndatasize];
	double* phi = new double[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		r[i] = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
		theta[i] = std::acos(z[i] / r[i]);
		phi[i] = std::atan2(y[i], x[i]);
	}
	double* _1harm_r = new double[ndatasize];
	double* _1harm_theta = new double[ndatasize];
	double* _1harm_phi = new double[ndatasize];
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic(1, sigma, ndegree, norder, theta, phi, _1harm_r, _1harm_theta, _1harm_phi, ndatasize);
	std::complex<double>* uv1_r = new std::complex<double>[ndatasize];
	std::complex<double>* uv1_theta = new std::complex<double>[ndatasize];
	std::complex<double>* uv1_phi = new std::complex<double>[ndatasize];
	if (type == 'v')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double jl = MathBasic::sph_bessel(ndegree, k * r[i]);//C++提供的Bessel函数不能计算复数，后面可以使用其它的函数库
			uv1_r[i] = jl * _1harm_r[i];
			uv1_theta[i] = jl * _1harm_theta[i];
			uv1_phi[i] = jl * _1harm_phi[i];
		}
	}
	else if (type == 'u')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			std::complex<double> h2(MathBasic::sph_hankel(ndegree, 2, k * r[i]));
			uv1_r[i] = h2 * _1harm_r[i];
			uv1_theta[i] = h2 * _1harm_theta[i];
			uv1_phi[i] = h2 * _1harm_phi[i];
		}
	}
	//向量的球-笛卡尔转换
	for (int i = 0; i < ndatasize; i++)
	{
		using namespace std;
		res_x[i] = uv1_r[i] * sin(theta[i]) * cos(phi[i]) + uv1_theta[i] * cos(theta[i]) * cos(phi[i]) - uv1_phi[i] * sin(phi[i]);
		res_y[i] = uv1_r[i] * sin(theta[i]) * sin(phi[i]) + uv1_theta[i] * cos(theta[i]) * sin(phi[i]) + uv1_phi[i] * cos(phi[i]);
		res_z[i] = uv1_r[i] * cos(theta[i]) - +uv1_theta[i] * sin(theta[i]);
	}
	delete[] r;
	delete[] theta;
	delete[] phi;
	delete[] _1harm_r;
	delete[] _1harm_theta;
	delete[] _1harm_phi;
	delete[] uv1_r;
	delete[] uv1_theta;
	delete[] uv1_phi;
}

void MathBasic::vec_sph_wav::s_gen_vec_sph_uv2(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, double r, double t, double p, std::complex<double>* res_r, std::complex<double>* res_t, std::complex<double>* res_p)
{
	double _2harm_r;
	double _2harm_theta;
	double _2harm_phi;
	double _3harm_r;
	double _3harm_theta;
	double _3harm_phi;
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic(2, sigma, ndegree, norder, &t, &p, &_2harm_r, &_2harm_theta, &_2harm_phi, 1);
	gen_vec_sph_harmonic(3, sigma, ndegree, norder, &t, &p, &_3harm_r, &_3harm_theta, &_3harm_phi, 1);
	if (type == 'v')
	{
		double kr = k * r;
		double jl2 = ((ndegree + 1) * MathBasic::sph_bessel(ndegree, kr) - kr * MathBasic::sph_bessel(ndegree + 1, kr)) / kr;
		double jl3 = sqrt(ndegree * (ndegree + 1)) * MathBasic::sph_bessel(ndegree, kr) / kr;
		*res_r = jl2 * _2harm_r + jl3 * _3harm_r;
		*res_t = jl2 * _2harm_theta + jl3 * _3harm_theta;
		*res_p = jl2 * _2harm_phi + jl3 * _3harm_phi;
	}
	else if (type == 'u')
	{
		double kr = k * r;
		std::complex<double> hl2(MathBasic::sph_hankel(ndegree, 2, kr));
		std::complex<double> hl2_p1(MathBasic::sph_hankel(ndegree + 1, 2, kr));
		std::complex<double> h2 = (double(ndegree + 1) * hl2 - kr * hl2_p1) / kr;
		std::complex<double> h3 = (sqrt(ndegree * (ndegree + 1)) * hl2) / kr;
		*res_r = h2 * _2harm_r + h3 * _3harm_r;
		*res_t = h2 * _2harm_theta + h3 * _3harm_theta;
		*res_p = h2 * _2harm_phi + h3 * _3harm_phi;
	}
}

double MathBasic::vec_sph_wav::Wigner3j(int j1, int j2, int j3, int m1, int m2, int m3)
{
	// Input error checking
	int m123[] = { m1, m2, m3 };
	int j123[] = { j1, j2, j3 };
	for (int j : j123) {
		if (j < 0) {
			std::cerr << "The j must be non-negative" << std::endl;
			return 0;
		}
		if (std::fmod(j, 0.5) != 0) {
			std::cerr << "All arguments must be integers or half-integers" << std::endl;
			return 0;
		}
	}
	for (int jm : m123) {
		if (std::fmod(jm, 0.5) != 0) {
			std::cerr << "All arguments must be integers or half-integers" << std::endl;
			return 0;
		}
	}
	if (std::fmod(j1 - m1, 1) != 0 || std::fmod(j2 - m2, 1) != 0 || std::fmod(j3 - m3, 1) != 0) {
		std::cerr << "j123 and m123 do not match" << std::endl;
		return 0;
	}

	// Selection rules
	if (j3 > (j1 + j2) || j3 < std::abs(j1 - j2) || (m1 + m2 + m3 != 0) ||
		std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m3) > j3) {
		return 0;
	}

	// Simple common case
	if (m1 == 0 && m2 == 0 && m3 == 0 && std::fmod(j1 + j2 + j3, 2) != 0) {
		return 0;
	}

	// Evaluation
	int t1 = j2 - m1 - j3;
	int t2 = j1 + m2 - j3;
	int t3 = j1 + j2 - j3;
	int t4 = j1 - m1;
	int t5 = j2 + m2;
	int tmin = std::max(0, std::max(t1, t2));
	int tmax = std::min(t3, std::min(t4, t5));

	double w = 0.0;
	using namespace std;
	double expo2 = 0.5*(-lgamma(j1 + j2 + j3 + 2) +lgamma(j1 + j2 - j3 + 1) + lgamma(j1 - j2 + j3 + 1) + lgamma(-j1 + j2 + j3 + 1) +
				lgamma(j1 + m1 + 1) + lgamma(j1 - m1 + 1) + lgamma(j2 + m2 + 1) 
		+ lgamma(j2 - m2 + 1) + lgamma(j3 + m3 + 1) + lgamma(j3 - m3 + 1));
	for (int t = tmin; t <= tmax; t++) 
	{
		double expo1 = -lgamma(t + 1) - lgamma(t - t1 + 1) - lgamma(t - t2 + 1) - lgamma(t3 - t + 1) -
			lgamma(t4 - t + 1) - lgamma(t5 - t + 1);
		w += pow(-1, t) * exp(expo1 + expo2);
	}
	w *= pow(-1, j1 - j2 - m3);;
	return w;
}

double MathBasic::vec_sph_wav::Wigner3j_precalc(int j1, int j2, int j3, int m1, int m2, int m3)
{
	long long c;
	int S, L, X, B, T;
	S = -j1 + j2 + j3;
	L = j1 - j2 + j3;
	X = j1 - m1;
	B = j2 - m2;
	T = j3 + m3;
	c = L * (24 + L + L * (50 + L * (35 + L * (10 + L)))) / 120 + X * (6 + X * (11 + X * (6 + X))) / 24 
		+ T * (2 + T * (3 + T)) / 6 + B * (B + 1) / 2 + S + 1;
	return 1;
}

void MathBasic::vec_sph_wav::calc_wav_index(int degree, std::vector<std::vector<int>>& n_index)
{
	int Nsize = 2 * degree*(degree + 2);
	n_index.resize(Nsize);
	for (int i = 0; i < Nsize; i++)
	{
		n_index[i].resize(4);
	}
	//计算n index
	for (int it = 1; it <= 2; it++)
	{
		for (int id = 1; id <= degree; id++)
		{
			for (int is = 0; is <= 1; is++)
			{
				for (int io = 0; io <= id; io++)
				{
					if (io == 0 && is == 1)
					{
						continue;
					}
					int n = 2 * (id * (id + 1) + io * pow(-1, is) + -1) + it;
					n_index[n - 1][0] = id;
					n_index[n - 1][1] = io;
					n_index[n - 1][2] = is;
					n_index[n - 1][3] = it;
				}
			}
		}
	}
}

void MathBasic::vec_sph_wav::calc_wav_index_cpx(int degree, std::vector<std::vector<int>>& n_index)
{
	int Nsize = 2 * degree * (degree + 2);
	n_index.resize(Nsize);
	for (int i = 0; i < Nsize; i++)
	{
		n_index[i].resize(3);
	}
	//计算n index
	for (int it = 1; it <= 2; it++)
	{
		for (int id = 1; id <= degree; id++)
		{
			for (int io = -id; io <= id; io++)
			{
				int n = 2 * (id * (id + 1) + io + -1) + it;
				n_index[n - 1][0] = id;
				n_index[n - 1][1] = io;
				n_index[n - 1][2] = it;
			}
		}
	}
}

void MathBasic::vec_sph_wav::gen_vec_sph_uv2(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, double* x, double* y, double* z, std::complex<double>* res_x, std::complex<double>* res_y, std::complex<double>* res_z, int ndatasize /*= 1*/)
{
	//将直角坐标换算成球坐标
	double* r = new double[ndatasize];
	double* theta = new double[ndatasize];
	double* phi = new double[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		r[i] = std::sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
		theta[i] = std::acos(z[i] / r[i]);
		phi[i] = std::atan2(y[i], x[i]);
	}
	double* _2harm_r = new double[ndatasize];
	double* _2harm_theta = new double[ndatasize];
	double* _2harm_phi = new double[ndatasize];
	double* _3harm_r = new double[ndatasize];
	double* _3harm_theta = new double[ndatasize];
	double* _3harm_phi = new double[ndatasize];
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic(2, sigma, ndegree, norder, theta, phi, _2harm_r, _2harm_theta, _2harm_phi, ndatasize);
	gen_vec_sph_harmonic(3, sigma, ndegree, norder, theta, phi, _3harm_r, _3harm_theta, _3harm_phi, ndatasize);
	std::complex<double>* uv2_r = new std::complex<double>[ndatasize];
	std::complex<double>* uv2_theta = new std::complex<double>[ndatasize];
	std::complex<double>* uv2_phi = new std::complex<double>[ndatasize];
	if (type == 'v')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double kr = k * r[i];
			double jl2 = ((ndegree + 1) * MathBasic::sph_bessel(ndegree, kr) - kr * MathBasic::sph_bessel(ndegree + 1, kr)) / kr;
			double jl3 = sqrt(ndegree * (ndegree + 1)) * MathBasic::sph_bessel(ndegree, kr) / kr;
			uv2_r[i] = jl2 * _2harm_r[i] + jl3 * _3harm_r[i];
			uv2_theta[i] = jl2 * _2harm_theta[i] + jl3 * _3harm_theta[i];
			uv2_phi[i] = jl2 * _2harm_phi[i] + jl3 * _3harm_phi[i];
		}
	}
	else if (type == 'u')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double kr = k * r[i];
			std::complex<double> hl2(MathBasic::sph_hankel(ndegree, 2, kr));
			std::complex<double> hl2_p1(MathBasic::sph_hankel(ndegree + 1, 2, kr));
			std::complex<double> h2 = (double(ndegree + 1) * hl2 - kr * hl2_p1) / kr;
			std::complex<double> h3 = (sqrt(ndegree * (ndegree + 1)) * hl2) / kr;
			uv2_r[i] = h2 * _2harm_r[i] + h3 * _3harm_r[i];
			uv2_theta[i] = h2 * _2harm_theta[i] + h3 * _3harm_theta[i];
			uv2_phi[i] = h2 * _2harm_phi[i] + h3 * _3harm_phi[i];
		}
	}
	//球坐标到直角坐标转换
	for (int i = 0; i < ndatasize; i++)
	{
		using namespace std;
		res_x[i] = uv2_r[i] * sin(theta[i]) * cos(phi[i]) + uv2_theta[i] * cos(theta[i]) * cos(phi[i]) - uv2_phi[i] * sin(phi[i]);
		res_y[i] = uv2_r[i] * sin(theta[i]) * sin(phi[i]) + uv2_theta[i] * cos(theta[i]) * sin(phi[i]) + uv2_phi[i] * cos(phi[i]);
		res_z[i] = uv2_r[i] * cos(theta[i]) - +uv2_theta[i] * sin(theta[i]);
	}
	delete[] r;
	delete[] theta;
	delete[] phi;
	delete[] _2harm_r;
	delete[] _2harm_theta;
	delete[] _2harm_phi;
	delete[] _3harm_r;
	delete[] _3harm_theta;
	delete[] _3harm_phi;
	delete[] uv2_r;
	delete[] uv2_theta;
	delete[] uv2_phi;
}

void MathBasic::vec_sph_wav::gen_vec_sph_uv1_cpx(char type/*u波或v波*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize /*= 1*/)
{
	//坐标的笛卡尔-球转换
	double* r = new double[ndatasize];
	double* theta = new double[ndatasize];
	double* phi = new double[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		double& x = ((cor*)cart)[i].x;
		double& y = ((cor*)cart)[i].y;
		double& z = ((cor*)cart)[i].z;
		r[i] = std::sqrt(x * x + y * y + z * z);
		theta[i] = std::acos(z / r[i]);
		phi[i] = std::atan2(y, x);
	}
	std::complex<double>* _1harm_r = new std::complex<double>[ndatasize];
	std::complex<double>* _1harm_theta = new std::complex<double>[ndatasize];
	std::complex<double>* _1harm_phi = new std::complex<double>[ndatasize];
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic_cpx(1,ndegree, norder, theta, phi, _1harm_r, _1harm_theta, _1harm_phi, ndatasize);
	std::complex<double>* uv1_r = new std::complex<double>[ndatasize];
	std::complex<double>* uv1_theta = new std::complex<double>[ndatasize];
	std::complex<double>* uv1_phi = new std::complex<double>[ndatasize];
	if (type == 'v')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double jl = MathBasic::sph_bessel(ndegree, k * r[i]);//C++提供的Bessel函数不能计算复数，后面可以使用其它的函数库
			uv1_r[i] = jl * _1harm_r[i];
			uv1_theta[i] = jl * _1harm_theta[i];
			uv1_phi[i] = jl * _1harm_phi[i];
		}
	}
	else if (type == 'u')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			std::complex<double> h2 = MathBasic::sph_hankel(ndegree, 2, k * r[i]);
			uv1_r[i] = h2 * _1harm_r[i];
			uv1_theta[i] = h2 * _1harm_theta[i];
			uv1_phi[i] = h2 * _1harm_phi[i];
		}
	}
	//向量的球-笛卡尔转换
	for (int i = 0; i < ndatasize; i++)
	{
		using namespace std;
		((cpx_cor*)res_cart)[i].x = uv1_r[i] * sin(theta[i]) * cos(phi[i]) + uv1_theta[i] * cos(theta[i]) * cos(phi[i]) - uv1_phi[i] * sin(phi[i]);
		((cpx_cor*)res_cart)[i].y = uv1_r[i] * sin(theta[i]) * sin(phi[i]) + uv1_theta[i] * cos(theta[i]) * sin(phi[i]) + uv1_phi[i] * cos(phi[i]);
		((cpx_cor*)res_cart)[i].z = uv1_r[i] * cos(theta[i]) - uv1_theta[i] * sin(theta[i]);
	}
	delete[] r;
	delete[] theta;
	delete[] phi;
	delete[] _1harm_r;
	delete[] _1harm_theta;
	delete[] _1harm_phi;
	delete[] uv1_r;
	delete[] uv1_theta;
	delete[] uv1_phi;
}

void MathBasic::vec_sph_wav::gen_vec_sph_uv2_cpx(char type/*u波或v波*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize /*= 1*/)

{
	//将直角坐标换算成球坐标
	double* r = new double[ndatasize];
	double* theta = new double[ndatasize];
	double* phi = new double[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		double& x = ((cor*)cart)[i].x;
		double& y = ((cor*)cart)[i].y;
		double& z = ((cor*)cart)[i].z;
		r[i] = std::sqrt(x * x + y * y + z * z);
		theta[i] = std::acos(z / r[i]);
		phi[i] = std::atan2(y, x);
	}
	std::complex<double>* _2harm_r = new std::complex<double>[ndatasize];
	std::complex<double>* _2harm_theta = new std::complex<double>[ndatasize];
	std::complex<double>* _2harm_phi = new std::complex<double>[ndatasize];
	std::complex<double>* _3harm_r = new std::complex<double>[ndatasize];
	std::complex<double>* _3harm_theta = new std::complex<double>[ndatasize];
	std::complex<double>* _3harm_phi = new std::complex<double>[ndatasize];
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic_cpx(2, ndegree, norder, theta, phi, _2harm_r, _2harm_theta, _2harm_phi, ndatasize);
	gen_vec_sph_harmonic_cpx(3, ndegree, norder, theta, phi, _3harm_r, _3harm_theta, _3harm_phi, ndatasize);
	std::complex<double>* uv2_r = new std::complex<double>[ndatasize];
	std::complex<double>* uv2_theta = new std::complex<double>[ndatasize];
	std::complex<double>* uv2_phi = new std::complex<double>[ndatasize];
	if (type == 'v')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double kr = k * r[i];
			double jl2 = ((ndegree + 1) * MathBasic::sph_bessel(ndegree, kr) - kr * MathBasic::sph_bessel(ndegree + 1, kr)) / kr;
			double jl3 = sqrt(ndegree * (ndegree + 1)) * MathBasic::sph_bessel(ndegree, kr) / kr;
			uv2_r[i] = jl2 * _2harm_r[i] + jl3 * _3harm_r[i];
			uv2_theta[i] = jl2 * _2harm_theta[i] + jl3 * _3harm_theta[i];
			uv2_phi[i] = jl2 * _2harm_phi[i] + jl3 * _3harm_phi[i];
		}
	}
	else if (type == 'u')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double kr = k * r[i];
			std::complex<double> hl2(MathBasic::sph_hankel(ndegree, 2, kr));
			std::complex<double> hl2_p1(MathBasic::sph_hankel(ndegree + 1, 2, kr));
			std::complex<double> h2 = (double(ndegree + 1) * hl2 - kr * hl2_p1) / kr;
			std::complex<double> h3 = (sqrt(ndegree * (ndegree + 1)) * hl2) / kr;
			uv2_r[i] = h2 * _2harm_r[i] + h3 * _3harm_r[i];
			uv2_theta[i] = h2 * _2harm_theta[i] + h3 * _3harm_theta[i];
			uv2_phi[i] = h2 * _2harm_phi[i] + h3 * _3harm_phi[i];
		}
	}
	//球坐标到直角坐标转换
	for (int i = 0; i < ndatasize; i++)
	{
		using namespace std;
		((cpx_cor*)res_cart)[i].x = uv2_r[i] * sin(theta[i]) * cos(phi[i]) + uv2_theta[i] * cos(theta[i]) * cos(phi[i]) - uv2_phi[i] * sin(phi[i]);
		((cpx_cor*)res_cart)[i].y = uv2_r[i] * sin(theta[i]) * sin(phi[i]) + uv2_theta[i] * cos(theta[i]) * sin(phi[i]) + uv2_phi[i] * cos(phi[i]);
		((cpx_cor*)res_cart)[i].z = uv2_r[i] * cos(theta[i]) - uv2_theta[i] * sin(theta[i]);
	}
	delete[] r;
	delete[] theta;
	delete[] phi;
	delete[] _2harm_r;
	delete[] _2harm_theta;
	delete[] _2harm_phi;
	delete[] _3harm_r;
	delete[] _3harm_theta;
	delete[] _3harm_phi;
	delete[] uv2_r;
	delete[] uv2_theta;
	delete[] uv2_phi;
}

void MathBasic::vec_sph_wav::gen_far_vsw_u1(char type, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize /*= 1*/)
{
	//坐标的笛卡尔-球转换
	double* r = new double[ndatasize];
	double* theta = new double[ndatasize];
	double* phi = new double[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		double& x = ((cor*)cart)[i].x;
		double& y = ((cor*)cart)[i].y;
		double& z = ((cor*)cart)[i].z;
		r[i] = std::sqrt(x * x + y * y + z * z);
		theta[i] = std::acos(z / r[i]);
		//phi[i] = std::atan2(y, x);
		if (x == 0 && y == 0)
		{
			phi[i] = 0;
		}
		else
		{
			phi[i] = y >= 0 ? std::acos(x / std::sqrt(x * x + y * y)) : 2 * pi - std::acos(x / std::sqrt(x * x + y * y));
		}		
	}
	double* _1harm_r = new double[ndatasize];
	double* _1harm_theta = new double[ndatasize];
	double* _1harm_phi = new double[ndatasize];
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic(1, sigma, ndegree, norder, theta, phi, _1harm_r, _1harm_theta, _1harm_phi, ndatasize);
	std::complex<double>* uv1_r = new std::complex<double>[ndatasize];
	std::complex<double>* uv1_theta = new std::complex<double>[ndatasize];
	std::complex<double>* uv1_phi = new std::complex<double>[ndatasize];

	for (int i = 0; i < ndatasize; i++)
	{
		double kr = k * r[i];
		std::complex<double>jx(0, -kr + pi / 2 * (ndegree + 1));
		std::complex<double> h2 = std::exp(jx) / kr;//球hankel函数的极限
		
		uv1_r[i] = h2 * _1harm_r[i];
		uv1_theta[i] = h2 * _1harm_theta[i];
		uv1_phi[i] = h2 * _1harm_phi[i];
	}

	//向量的球-笛卡尔转换
	for (int i = 0; i < ndatasize; i++)
	{
		using namespace std;
		((cpx_cor*)res_cart)[i].x = uv1_r[i] * sin(theta[i]) * cos(phi[i]) + uv1_theta[i] * cos(theta[i]) * cos(phi[i]) - uv1_phi[i] * sin(phi[i]);
		((cpx_cor*)res_cart)[i].y = uv1_r[i] * sin(theta[i]) * sin(phi[i]) + uv1_theta[i] * cos(theta[i]) * sin(phi[i]) + uv1_phi[i] * cos(phi[i]);
		((cpx_cor*)res_cart)[i].z = uv1_r[i] * cos(theta[i]) - uv1_theta[i] * sin(theta[i]);
	}
	delete[] r;
	delete[] theta;
	delete[] phi;
	delete[] _1harm_r;
	delete[] _1harm_theta;
	delete[] _1harm_phi;
	delete[] uv1_r;
	delete[] uv1_theta;
	delete[] uv1_phi;
}

void MathBasic::vec_sph_wav::gen_far_vsw_u2(char type, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize /*= 1*/)
{
	//将直角坐标换算成球坐标
	double* r = new double[ndatasize];
	double* theta = new double[ndatasize];
	double* phi = new double[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		double& x = ((cor*)cart)[i].x;
		double& y = ((cor*)cart)[i].y;
		double& z = ((cor*)cart)[i].z;
		r[i] = std::sqrt(x * x + y * y + z * z);
		theta[i] = std::acos(z / r[i]);
		//phi[i] = std::atan2(y, x);
		if (x == 0 && y == 0)
		{
			phi[i] = 0;
		}
		else
		{
			phi[i] = y >= 0 ? std::acos(x / std::sqrt(x * x + y * y)) : 2 * pi - std::acos(x / std::sqrt(x * x + y * y));
		}
	}
	double* _2harm_r = new double[ndatasize];
	double* _2harm_theta = new double[ndatasize];
	double* _2harm_phi = new double[ndatasize];
	double* _3harm_r = new double[ndatasize];
	double* _3harm_theta = new double[ndatasize];
	double* _3harm_phi = new double[ndatasize];
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic(2, sigma, ndegree, norder, theta, phi, _2harm_r, _2harm_theta, _2harm_phi, ndatasize);
	gen_vec_sph_harmonic(3, sigma, ndegree, norder, theta, phi, _3harm_r, _3harm_theta, _3harm_phi, ndatasize);
	std::complex<double>* uv2_r = new std::complex<double>[ndatasize];
	std::complex<double>* uv2_theta = new std::complex<double>[ndatasize];
	std::complex<double>* uv2_phi = new std::complex<double>[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		double kr = k * r[i];
		std::complex<double>jx(0, -kr + pi / 2 * (ndegree + 1 + 1));//直接从不取极限的式子里近似出来。
		std::complex<double> hl2_p1 = std::exp(jx) / kr;//球hankel函数的极限
		std::complex<double> h2 = -hl2_p1;//总的效果相当于jx(0, -kr + pi / 2 * ndegree);h2 = std::exp(jx) / kr
		uv2_r[i] = h2 * _2harm_r[i];
		uv2_theta[i] = h2 * _2harm_theta[i];
		uv2_phi[i] = h2 * _2harm_phi[i];
	}
	//球坐标到直角坐标转换
	for (int i = 0; i < ndatasize; i++)
	{
		using namespace std;
		((cpx_cor*)res_cart)[i].x = uv2_r[i] * sin(theta[i]) * cos(phi[i]) + uv2_theta[i] * cos(theta[i]) * cos(phi[i]) - uv2_phi[i] * sin(phi[i]);
		((cpx_cor*)res_cart)[i].y = uv2_r[i] * sin(theta[i]) * sin(phi[i]) + uv2_theta[i] * cos(theta[i]) * sin(phi[i]) + uv2_phi[i] * cos(phi[i]);
		((cpx_cor*)res_cart)[i].z = uv2_r[i] * cos(theta[i]) - uv2_theta[i] * sin(theta[i]);
	}
	delete[] r;
	delete[] theta;
	delete[] phi;
	delete[] _2harm_r;
	delete[] _2harm_theta;
	delete[] _2harm_phi;
	delete[] _3harm_r;
	delete[] _3harm_theta;
	delete[] _3harm_phi;
	delete[] uv2_r;
	delete[] uv2_theta;
	delete[] uv2_phi;
}

void MathBasic::vec_sph_wav::s_gen_vec_sph_uv1(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, double r, double t, double p, std::complex<double>* res_r, std::complex<double>* res_t, std::complex<double>* res_p)
{
	double _1harm_r;
	double _1harm_theta;
	double _1harm_phi;
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic(1, sigma, ndegree, norder, &t, &p, &_1harm_r, &_1harm_theta, &_1harm_phi, 1);
	if (type == 'v')
	{
		double jl = MathBasic::sph_bessel(ndegree, k * r);//C++提供的Bessel函数不能计算复数，后面可以使用其它的函数库
		*res_r = jl * _1harm_r;
		*res_t = jl * _1harm_theta;
		*res_p = jl * _1harm_phi;
	}
	else if (type == 'u')
	{
		std::complex<double> h2(MathBasic::sph_hankel(ndegree, 2, k * r));
		*res_r = h2 * _1harm_r;
		*res_t = h2 * _1harm_theta;
		*res_p = h2 * _1harm_phi;
	}
}

void MathBasic::vec_sph_wav::gen_vec_sph_uv1(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize /*= 1*/)
{
	//坐标的笛卡尔-球转换
	double* r = new double[ndatasize];
	double* theta = new double[ndatasize];
	double* phi = new double[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		double& x = ((cor*)cart)[i].x;
		double& y = ((cor*)cart)[i].y;
		double& z = ((cor*)cart)[i].z;
		r[i] = std::sqrt(x * x + y * y + z * z);
		theta[i] = std::acos(z / r[i]);
		//phi[i] = std::atan2(y, x);
		if (x == 0 && y == 0)
		{
			phi[i] = 0;
		}
		else
		{
			phi[i] = y >= 0 ? std::acos(x / std::sqrt(x * x + y * y)) : 2 * pi - std::acos(x / std::sqrt(x * x + y * y));
		}
	}
	double* _1harm_r = new double[ndatasize];
	double* _1harm_theta = new double[ndatasize];
	double* _1harm_phi = new double[ndatasize];
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic(1, sigma, ndegree, norder, theta, phi, _1harm_r, _1harm_theta, _1harm_phi, ndatasize);
	std::complex<double>* uv1_r = new std::complex<double>[ndatasize];
	std::complex<double>* uv1_theta = new std::complex<double>[ndatasize];
	std::complex<double>* uv1_phi = new std::complex<double>[ndatasize];
	if (type == 'v')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double jl = MathBasic::sph_bessel(ndegree, k * r[i]);//C++提供的Bessel函数不能计算复数，后面可以使用其它的函数库
			uv1_r[i] = jl * _1harm_r[i];
			uv1_theta[i] = jl * _1harm_theta[i];
			uv1_phi[i] = jl * _1harm_phi[i];
		}
	}
	else if (type == 'u')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			std::complex<double> h2 = MathBasic::sph_hankel(ndegree, 2, k * r[i]);
			//double h2 = 1;
			uv1_r[i] = h2 * _1harm_r[i];
			uv1_theta[i] = h2 * _1harm_theta[i];
			uv1_phi[i] = h2 * _1harm_phi[i];
		}
	}
	//向量的球-笛卡尔转换
	for (int i = 0; i < ndatasize; i++)
	{
		using namespace std;
		((cpx_cor*)res_cart)[i].x = uv1_r[i] * sin(theta[i]) * cos(phi[i]) + uv1_theta[i] * cos(theta[i]) * cos(phi[i]) - uv1_phi[i] * sin(phi[i]);
		((cpx_cor*)res_cart)[i].y = uv1_r[i] * sin(theta[i]) * sin(phi[i]) + uv1_theta[i] * cos(theta[i]) * sin(phi[i]) + uv1_phi[i] * cos(phi[i]);
		((cpx_cor*)res_cart)[i].z = uv1_r[i] * cos(theta[i]) - uv1_theta[i] * sin(theta[i]);
	}
	delete[] r;
	delete[] theta;
	delete[] phi;
	delete[] _1harm_r;
	delete[] _1harm_theta;
	delete[] _1harm_phi;
	delete[] uv1_r;
	delete[] uv1_theta;
	delete[] uv1_phi;
}

void MathBasic::vec_sph_wav::gen_vec_sph_uv2(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize /*= 1*/)
{
	//将直角坐标换算成球坐标
	double* r = new double[ndatasize];
	double* theta = new double[ndatasize];
	double* phi = new double[ndatasize];
	for (int i = 0; i < ndatasize; i++)
	{
		double& x = ((cor*)cart)[i].x;
		double& y = ((cor*)cart)[i].y;
		double& z = ((cor*)cart)[i].z;
		r[i] = std::sqrt(x * x + y * y + z * z);
		theta[i] = std::acos(z / r[i]);
		//phi[i] = std::atan2(y, x);
		if (x == 0 && y == 0)
		{
			phi[i] = 0;
		}
		else
		{
			phi[i] = y >= 0 ? std::acos(x / std::sqrt(x * x + y * y)) : 2 * pi - std::acos(x / std::sqrt(x * x + y * y));
		}
	}
	double* _2harm_r = new double[ndatasize];
	double* _2harm_theta = new double[ndatasize];
	double* _2harm_phi = new double[ndatasize];
	double* _3harm_r = new double[ndatasize];
	double* _3harm_theta = new double[ndatasize];
	double* _3harm_phi = new double[ndatasize];
	//计算三种tau的矢量球谐函数
	gen_vec_sph_harmonic(2, sigma, ndegree, norder, theta, phi, _2harm_r, _2harm_theta, _2harm_phi, ndatasize);
	gen_vec_sph_harmonic(3, sigma, ndegree, norder, theta, phi, _3harm_r, _3harm_theta, _3harm_phi, ndatasize);
	std::complex<double>* uv2_r = new std::complex<double>[ndatasize];
	std::complex<double>* uv2_theta = new std::complex<double>[ndatasize];
	std::complex<double>* uv2_phi = new std::complex<double>[ndatasize];
	if (type == 'v')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double kr = k * r[i];
			double jl2 = ((ndegree + 1) * MathBasic::sph_bessel(ndegree, kr) - kr * MathBasic::sph_bessel(ndegree + 1, kr)) / kr;
			double jl3 = sqrt(ndegree * (ndegree + 1)) * MathBasic::sph_bessel(ndegree, kr) / kr;
			uv2_r[i] = jl2 * _2harm_r[i] + jl3 * _3harm_r[i];
			uv2_theta[i] = jl2 * _2harm_theta[i] + jl3 * _3harm_theta[i];
			uv2_phi[i] = jl2 * _2harm_phi[i] + jl3 * _3harm_phi[i];
		}
	}
	else if (type == 'u')
	{
		for (int i = 0; i < ndatasize; i++)
		{
			double kr = k * r[i];
			std::complex<double> hl2(MathBasic::sph_hankel(ndegree, 2, kr));
			std::complex<double> hl2_p1(MathBasic::sph_hankel(ndegree + 1, 2, kr));
			std::complex<double> h2 = (double(ndegree + 1) * hl2 - kr * hl2_p1) / kr;
			std::complex<double> h3 = (sqrt(ndegree * (ndegree + 1)) * hl2) / kr;
			uv2_r[i] = h2 * _2harm_r[i] + h3 * _3harm_r[i];
			uv2_theta[i] = h2 * _2harm_theta[i] + h3 * _3harm_theta[i];
			uv2_phi[i] = h2 * _2harm_phi[i] + h3 * _3harm_phi[i];
		}
	}
	//球坐标到直角坐标转换
	for (int i = 0; i < ndatasize; i++)
	{
		using namespace std;
		((cpx_cor*)res_cart)[i].x = uv2_r[i] * sin(theta[i]) * cos(phi[i]) + uv2_theta[i] * cos(theta[i]) * cos(phi[i]) - uv2_phi[i] * sin(phi[i]);
		((cpx_cor*)res_cart)[i].y = uv2_r[i] * sin(theta[i]) * sin(phi[i]) + uv2_theta[i] * cos(theta[i]) * sin(phi[i]) + uv2_phi[i] * cos(phi[i]);
		((cpx_cor*)res_cart)[i].z = uv2_r[i] * cos(theta[i]) - uv2_theta[i] * sin(theta[i]);
	}
	delete[] r;
	delete[] theta;
	delete[] phi;
	delete[] _2harm_r;
	delete[] _2harm_theta;
	delete[] _2harm_phi;
	delete[] _3harm_r;
	delete[] _3harm_theta;
	delete[] _3harm_phi;
	delete[] uv2_r;
	delete[] uv2_theta;
	delete[] uv2_phi;
}

double MathBasic::vec_sph_wav::Derivative_assoc_legendre(int ndegree, int norder, double theta)
{
	//该函数不能正确处理数值的奇异性（分母为0）
	double u = cos(theta);
	double res = ndegree * u * std::assoc_legendre(ndegree, norder, u) - (ndegree + norder) * std::assoc_legendre(ndegree - 1, norder, u);
	return res / sin(theta);
}

void MathBasic::vec_sph_wav::gen_vec_sph_harmonic(int tau, char sigma, int ndegree, int norder, double* theta, double* phi, double* res_r, double* res_theta, double* res_phi, int ndatasize)
{
	int em = norder == 0 ? 1 : 2;
	double jiecheng = 1;
	for (unsigned int i = ndegree + norder; i > ndegree - norder; --i)
	{
		jiecheng *= i;
	}
	double flagC = sqrt(1.0 / (4 * pi) * (2 * ndegree + 1) * em / jiecheng);
	double flag = sqrt(1.0 / (ndegree * (ndegree + 1))) * flagC;
	//assoc_legendre(3, 2, cos(pi/4)) ;//C++ 的定义中omit了前置符号位(-1)^m。
	
	for (int i = 0; i < ndatasize; i++)
	{
		//边界条件
		switch (tau)
		{
		case 1:{
			res_r[i] = 0;//径向分量为0
			//边界条件
			if (ndegree == 0 && norder == 0)
			{
				res_theta[i] = 0;
				res_phi[i] = 0;
				continue;
			}
			//theta=0/pi的极限
			if (theta[i] < 1e-15)
			{
				double sign = sqrt((2 * ndegree + 1) / 8.0 / pi);
				if (norder != 1)//delta_{m1}
				{
					res_theta[i] = 0;
					res_phi[i] = 0;
				}
				else
				{
					if (sigma == 'e')
					{
						res_theta[i] = sign * -sin(phi[i]);
						res_phi[i] = -sign * cos(phi[i]);
					}
					else if (sigma == 'o')
					{
						res_theta[i] = sign * cos(phi[i]);
						res_phi[i] = -sign * sin(phi[i]);
					}
				}
				continue;
			}
			else if(theta[i] > pi - 1e-15)
			{
				double sign0 = sqrt((2 * ndegree + 1) / 8.0 / pi);
				double sign = ndegree % 2 == 0 ? -sign0 : sign0;
				if (norder != 1)//delta_{m1}
				{
					res_theta[i] = 0;
					res_phi[i] = 0;
				}
				else
				{
					if (sigma == 'e')
					{
						res_theta[i] = sign * -sin(phi[i]);
						res_phi[i] = sign * cos(phi[i]);
					}
					else if (sigma == 'o')
					{
						res_theta[i] = sign * cos(phi[i]);
						res_phi[i] = sign * sin(phi[i]);
					}
				}
				continue;
			}
			//无数值奇异性时的计算公式
			if (sigma == 'e')
			{
				res_theta[i] = flag * std::assoc_legendre(ndegree, norder, cos(theta[i])) / sin(theta[i]) * (-norder) * sin(norder * phi[i]);
				res_phi[i] = -flag * Derivative_assoc_legendre(ndegree, norder, theta[i]) * cos(norder * phi[i]);
			}
			else if (sigma == 'o')
			{
				res_theta[i] = flag * std::assoc_legendre(ndegree, norder, cos(theta[i])) / sin(theta[i]) * (norder)*cos(norder * phi[i]);
				res_phi[i] = -flag * Derivative_assoc_legendre(ndegree, norder, theta[i]) * sin(norder * phi[i]);
			}
			break;
		}
		case 2: {
			res_r[i] = 0;//径向分量为0
			//边界条件
			if (ndegree == 0 && norder == 0)
			{
				res_theta[i] = 0;
				res_phi[i] = 0;
				continue;
			}
			//theta=0/pi的极限
			if (theta[i] < 1e-15)
			{
				double sign = sqrt((2 * ndegree + 1) / 8.0 / pi);
				if (norder != 1)//delta_{m1}
				{
					res_theta[i] = 0;
					res_phi[i] = 0;
				}
				else
				{
					if (sigma == 'e')
					{
						res_theta[i] = sign * cos(phi[i]);
						res_phi[i] = sign * -sin(phi[i]);
					}
					else if (sigma == 'o')
					{
						res_theta[i] = sign * sin(phi[i]); 
						res_phi[i] = sign * cos(phi[i]);
					}
				}
				continue;
			}
			else if (theta[i] > pi - 1e-15)
			{
				double sign0 = sqrt((2 * ndegree + 1) / 8.0 / pi);
				double sign = ndegree % 2 == 0 ? -sign0 : sign0;
				if (norder != 1)//delta_{m1}
				{
					res_theta[i] = 0;
					res_phi[i] = 0;
				}
				else
				{
					if (sigma == 'e')
					{
						res_theta[i] = -sign * cos(phi[i]); 
						res_phi[i] = sign * -sin(phi[i]);
					}
					else if (sigma == 'o')
					{
						res_theta[i] = -sign * sin(phi[i]);
						res_phi[i] = sign * cos(phi[i]);
					}
				}
				continue;
			}
			//无数值奇异性时直接计算
			if (sigma == 'e')
			{
				res_theta[i] = flag * Derivative_assoc_legendre(ndegree, norder, theta[i]) * cos(norder * phi[i]);
				res_phi[i] = flag * std::assoc_legendre(ndegree, norder, cos(theta[i])) / sin(theta[i]) * (-norder) * sin(norder * phi[i]);
			}
			else if (sigma == 'o')
			{
				res_theta[i] = flag * Derivative_assoc_legendre(ndegree, norder, theta[i]) * sin(norder * phi[i]);
				res_phi[i] = flag * std::assoc_legendre(ndegree, norder, cos(theta[i])) / sin(theta[i]) * (norder)*cos(norder * phi[i]);
			}
			break;
		}
		case 3: {
			res_theta[i] = 0;
			res_phi[i] = 0;
			//边界条件
			if (ndegree == 0 && norder == 0)
			{
				res_r[i] = sigma == 'e' ? sqrt(1.0 / 4 / pi) : 0;
				continue;
			}
			//其它情形时直接计算
			if (sigma == 'e')
			{
				res_r[i] = flagC * std::assoc_legendre(ndegree, norder, cos(theta[i])) * cos(norder * phi[i]);
			}
			else if (sigma == 'o')
			{
				res_r[i] = flagC * std::assoc_legendre(ndegree, norder, cos(theta[i])) * sin(norder * phi[i]);
			}
			break;
		}
		default:
			break;
		}
	}
}

void MathBasic::vec_sph_wav::gen_vec_sph_harmonic_cpx(int tau, int ndegree, int norder, double* theta, double* phi, std::complex<double>* res_r, std::complex<double>* res_theta, std::complex<double>* res_phi, int ndatasize)
{
	using std::complex;
	using std::exp;
	using std::abs;
	double jiecheng = 1;
	for (unsigned int i = ndegree + abs(norder); i > ndegree - abs(norder); --i)
	{
		jiecheng *= i;
	}
	double Condon_Shortley = (norder > 0 && abs(norder) % 2 != 0) ? -1 : 1;
	double flagC = Condon_Shortley * sqrt(1.0 / (4 * pi) * (2 * ndegree + 1) / jiecheng);//归一化因子
	double flag = sqrt(1.0 / (ndegree * (ndegree + 1))) * flagC;
	//assoc_legendre(3, 2, cos(pi/4)) ;//C++ 的定义中omit了前置符号位(-1)^m。

	for (int i = 0; i < ndatasize; i++)
	{
		//边界条件
		switch (tau)
		{
		case 1: {
			res_r[i] = 0;//径向分量为0
			//边界条件
			if (ndegree == 0 && norder == 0)
			{
				res_theta[i] = 0;
				res_phi[i] = 0;
				continue;
			}
			//theta=0/pi的极限
			if (theta[i] < 1e-15)
			{
				double sign = sqrt((2 * ndegree + 1) / 16.0 / pi);
				if (norder ==1)
				{
					res_theta[i] = -sign * exp(complex<double>(0, phi[i] + pi / 2));
					res_phi[i] = sign * exp(complex<double>(0, phi[i]));
				}
				else if(norder==-1)
				{
					res_theta[i] = -sign * exp(complex<double>(0, -phi[i] + pi / 2));
					res_phi[i] = -sign * exp(complex<double>(0, -phi[i]));
				}
				else//delta_{m+-1}
				{
					res_theta[i] = 0;
					res_phi[i] = 0;
				}
				continue;
			}
			else if (theta[i] > pi - 1e-15)
			{
				double sign0 = sqrt((2 * ndegree + 1) / 16.0 / pi);
				double sign = ndegree % 2 == 0 ? -sign0 : sign0;
				if (norder == 1)
				{
					res_theta[i] = -sign * exp(complex<double>(0, phi[i] + pi / 2));
					res_phi[i] = sign * exp(complex<double>(0, phi[i]));
				}
				else if (norder == -1)
				{
					res_theta[i] = -sign * exp(complex<double>(0, -phi[i] + pi / 2));
					res_phi[i] = -sign * exp(complex<double>(0, -phi[i]));
				}
				else//delta_{m+-1}
				{
					res_theta[i] = 0;
					res_phi[i] = 0;
				}
				continue;
			}
			//无数值奇异性时的计算公式
			res_theta[i] = flag * std::assoc_legendre(ndegree, abs(norder), cos(theta[i])) / sin(theta[i]) * norder * exp(complex<double>(0, norder * phi[i] + pi / 2));
			res_phi[i] = -flag * Derivative_assoc_legendre(ndegree, abs(norder), theta[i]) * exp(complex<double>(0, norder * phi[i]));
			break;
		}
		case 2: {
			res_r[i] = 0;//径向分量为0
			//边界条件
			if (ndegree == 0 && norder == 0)
			{
				res_theta[i] = 0;
				res_phi[i] = 0;
				continue;
			}
			//theta=0/pi的极限
			if (theta[i] < 1e-15)
			{
				double sign = sqrt((2 * ndegree + 1) / 16.0 / pi);
				if (norder == 1)
				{
					res_theta[i] = -sign * exp(complex<double>(0, phi[i]));
					res_phi[i] = -sign * exp(complex<double>(0, phi[i] + pi / 2));
				}
				else if (norder == -1)
				{
					res_theta[i] = sign * exp(complex<double>(0, -phi[i]));
					res_phi[i] = -sign * exp(complex<double>(0, -phi[i] + pi / 2));
				}
				else//delta_{m+-1}
				{
					res_theta[i] = 0;
					res_phi[i] = 0;
				}
				continue;
			}
			else if (theta[i] > pi - 1e-15)
			{
				double sign0 = sqrt((2 * ndegree + 1) / 16.0 / pi);
				double sign = ndegree % 2 == 0 ? -sign0 : sign0;
				if (norder == 1)
				{
					res_theta[i] = -sign * exp(complex<double>(0, phi[i]));
					res_phi[i] = -sign * exp(complex<double>(0, phi[i] + pi / 2));
				}
				else if (norder == -1)
				{
					res_theta[i] = sign * exp(complex<double>(0, -phi[i]));
					res_phi[i] = -sign * exp(complex<double>(0, -phi[i] + pi / 2));
				}
				else//delta_{m+-1}
				{
					res_theta[i] = 0;
					res_phi[i] = 0;
				}
				continue;
			}
			//无数值奇异性时直接计算
			res_theta[i] = flag * Derivative_assoc_legendre(ndegree, abs(norder), theta[i]) * exp(complex<double>(0, norder * phi[i]));
			res_phi[i] = flag * std::assoc_legendre(ndegree, abs(norder), cos(theta[i])) / sin(theta[i]) * norder * exp(complex<double>(0, norder * phi[i] + pi / 2));
			break;
		}
		case 3: {
			res_theta[i] = 0;
			res_phi[i] = 0;
			//边界条件
			if (ndegree == 0 && norder == 0)
			{
				res_r[i] = sqrt(1.0 / 4 / pi);
				continue;
			}
			//其它情形时直接计算
			res_r[i] = flagC * std::assoc_legendre(ndegree, abs(norder), cos(theta[i])) * exp(complex<double>(0, norder * phi[i]));
			break;
		}
		default:
			break;
		}
	}
}

double MathBasic::sph_bessel(int degree, double value)
{
	return std::sph_bessel(degree, value);
}

double MathBasic::sph_neumann(int degree, double value)
{
	return std::sph_neumann(degree, value);
}

std::complex<double> MathBasic::sph_hankel(int degree, int type, double value)
{
	return type == 1 ? std::complex<double>(MathBasic::sph_bessel(degree,value), MathBasic::sph_neumann(degree, value))
		: std::complex<double>(MathBasic::sph_bessel(degree, value), -MathBasic::sph_neumann(degree, value));
}
