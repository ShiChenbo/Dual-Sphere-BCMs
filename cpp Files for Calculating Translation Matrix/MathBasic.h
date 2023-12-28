#pragma once
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
namespace MathBasic
{
	const double pi = 3.141592653589793238462643383279;
	struct cor
	{
		double x, y, z;
		double norm()
		{
			return sqrt(abs(x * x + y * y + z * z));
		}
		cor operator-(cor& B)
		{
			cor C;
			C.x = x - B.x;
			C.y = y - B.y;
			C.z = z - B.z;
			return C;
		}
		void operator*=(double k)
		{
			x *= k;
			y *= k;
			z *= k;
		}
	};
	struct cpx_cor
	{
		std::complex<double> x, y, z;
		cpx_cor operator-(cpx_cor& B)
		{
			cpx_cor C;
			C.x = x - B.x;
			C.y = y - B.y;
			C.z = z - B.z;
			return C;
		}
		cpx_cor operator*(std::complex<double> k)
		{
			cpx_cor C;
			C.x = x*k;
			C.y = y*k;
			C.z = z*k;
			return C;
		}
		void operator*=(std::complex<double> k)
		{
			x *= k;
			y *= k;
			z *= k;
		}
		double norm()
		{
			return sqrt(abs(x * x + y * y + z * z));
		}
	};
	struct cpx_dyadic
	{
		std::complex<double> dya_core[3][3];
		void _2unit()
		{
			dya_core[0][0] = 1;
			dya_core[0][1] = 0;
			dya_core[0][2] = 0;
			dya_core[1][0] = 0;
			dya_core[1][1] = 1;
			dya_core[1][2] = 0;
			dya_core[2][0] = 0;
			dya_core[2][1] = 0;
			dya_core[2][2] = 1;
		}
		void make_dyadic(cpx_cor& A, cpx_cor& B)
		{
			dya_core[0][0] = A.x * B.x;
			dya_core[0][1] = A.x * B.y;
			dya_core[0][2] = A.x * B.z;
			dya_core[1][0] = A.y * B.x;
			dya_core[1][1] = A.y * B.y;
			dya_core[1][2] = A.y * B.z;
			dya_core[2][0] = A.z * B.x;
			dya_core[2][1] = A.z * B.y;
			dya_core[2][2] = A.z * B.z;
		}
		void make_dyadic(cor& A, cor& B)
		{
			dya_core[0][0] = A.x * B.x;
			dya_core[0][1] = A.x * B.y;
			dya_core[0][2] = A.x * B.z;
			dya_core[1][0] = A.y * B.x;
			dya_core[1][1] = A.y * B.y;
			dya_core[1][2] = A.y * B.z;
			dya_core[2][0] = A.z * B.x;
			dya_core[2][1] = A.z * B.y;
			dya_core[2][2] = A.z * B.z;
		}
		cpx_dyadic operator*(std::complex<double> k)
		{
			cpx_dyadic C;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					C.dya_core[i][j] = C.dya_core[i][j] * k;
				}
			}
			return C;
		}
		cpx_cor operator*(cpx_cor a)
		{
			cpx_cor C;
			C.x = dya_core[0][0] * a.x + dya_core[0][1] * a.y + dya_core[0][2] * a.z;
			C.y = dya_core[1][0] * a.x + dya_core[1][1] * a.y + dya_core[1][2] * a.z;
			C.z = dya_core[2][0] * a.x + dya_core[2][1] * a.y + dya_core[2][2] * a.z;
			return C;
		}
		void operator*=(std::complex<double> k)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					dya_core[i][j] *= k;
				}
			}
		}
		cpx_dyadic operator+(cpx_dyadic& B)
		{
			cpx_dyadic C;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					C.dya_core[i][j] = dya_core[i][j] + B.dya_core[i][j];
				}
			}
			return C;
		}
		void operator+=(cpx_dyadic& B)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					dya_core[i][j] += B.dya_core[i][j];
				}
			}
		}
		void print()
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					std::cout << dya_core[i][j]<<"\t";
				}
				std::cout << "\n";
			}
		}
	};
	inline double sph_bessel(int degree, double value);
	inline double sph_neumann(int degree, double value);
	inline std::complex<double> sph_hankel(int degree, int type, double value);
	namespace Quad
	{
		using namespace std;
		inline vector<double> GetLegendreRoots(int n, double precis = 1e-14, int maxIter = 1000);//通过迭代法找到Legendre多项式的n个根
		inline double DerivativeLegendre(int n, double x);//Lengendre函数的一阶微分
		void GenerateGuassIntegInfo(int n,vector<double>& Intpoints, vector<double>& weights);//n是积分点数
	}
	namespace vec_sph_wav
	{
		inline double Derivative_assoc_legendre(int ndegree, int norder, double theta);//连带legendre函数对theta的一阶微分
		inline void gen_vec_sph_harmonic(int tau, char sigma, int ndegree, int norder, double* theta, double* phi, double* res_r, double* res_theta, double* res_phi, int ndatasize);//生成矢量球谐
		void gen_vec_sph_harmonic_cpx(int tau, int ndegree, int norder, double* theta, double* phi, std::complex<double>* res_r, std::complex<double>* res_theta, std::complex<double>* res_phi, int ndatasize);//生成复数态矢量球谐
		//第1类u,v波（TE波）,v方向用球Bessel函数，u方向用第2类球Hankel函数
		void gen_vec_sph_uv1(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, double* x, double* y, double* z, std::complex<double>* res_x, std::complex<double>* res_y, std::complex<double>* res_z, int ndatasize = 1);
		//第2类u,v波（TM波）,v方向用球Bessel函数，u方向用第2类球Hankel函数
		void gen_vec_sph_uv2(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, double* x, double* y, double* z, std::complex<double>* res_x, std::complex<double>* res_y, std::complex<double>* res_z, int ndatasize = 1);
		//传递坐标结构的矢量波函数
		void gen_vec_sph_uv1(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize = 1);
		void gen_vec_sph_uv2(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize = 1);
		//复指数形式的矢量球面波
		void gen_vec_sph_uv1_cpx(char type/*u波或v波*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize = 1);
		void gen_vec_sph_uv2_cpx(char type/*u波或v波*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize = 1);
		//u形波的远场极限
		void gen_far_vsw_u1(char type/*形式参数*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize = 1);
		void gen_far_vsw_u2(char type/*形式参数*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, const void* cart, void* res_cart, int ndatasize = 1);
		//球坐标表示的矢量球面波
		void s_gen_vec_sph_uv1(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, double r, double t, double p, std::complex<double>* res_r, std::complex<double>* res_t, std::complex<double>* res_p);
		void s_gen_vec_sph_uv2(char type/*u波或v波*/, char sigma/*奇偶*/, double k/*波数*/, int ndegree/*度*/, int norder/*阶数*/, double r, double t, double p, std::complex<double>* res_r, std::complex<double>* res_t, std::complex<double>* res_p);
		double Wigner3j(int j1, int j2, int j3, int m1, int m2, int m3);
		double Wigner3j_precalc(int j1, int j2, int j3, int m1, int m2, int m3);
		void calc_wav_index(int degree, std::vector<std::vector<int>>& n_index);//将tslm索引映射成n索引，储存的顺序是l,m,sig,tau
		void calc_wav_index_cpx(int degree, std::vector<std::vector<int>>& n_index);//将tlm索引映射成n索引，储存的顺序是l,m,tau

	}

}


