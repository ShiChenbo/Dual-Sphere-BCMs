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
		inline vector<double> GetLegendreRoots(int n, double precis = 1e-14, int maxIter = 1000);//ͨ���������ҵ�Legendre����ʽ��n����
		inline double DerivativeLegendre(int n, double x);//Lengendre������һ��΢��
		void GenerateGuassIntegInfo(int n,vector<double>& Intpoints, vector<double>& weights);//n�ǻ��ֵ���
	}
	namespace vec_sph_wav
	{
		inline double Derivative_assoc_legendre(int ndegree, int norder, double theta);//����legendre������theta��һ��΢��
		inline void gen_vec_sph_harmonic(int tau, char sigma, int ndegree, int norder, double* theta, double* phi, double* res_r, double* res_theta, double* res_phi, int ndatasize);//����ʸ����г
		void gen_vec_sph_harmonic_cpx(int tau, int ndegree, int norder, double* theta, double* phi, std::complex<double>* res_r, std::complex<double>* res_theta, std::complex<double>* res_phi, int ndatasize);//���ɸ���̬ʸ����г
		//��1��u,v����TE����,v��������Bessel������u�����õ�2����Hankel����
		void gen_vec_sph_uv1(char type/*u����v��*/, char sigma/*��ż*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, double* x, double* y, double* z, std::complex<double>* res_x, std::complex<double>* res_y, std::complex<double>* res_z, int ndatasize = 1);
		//��2��u,v����TM����,v��������Bessel������u�����õ�2����Hankel����
		void gen_vec_sph_uv2(char type/*u����v��*/, char sigma/*��ż*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, double* x, double* y, double* z, std::complex<double>* res_x, std::complex<double>* res_y, std::complex<double>* res_z, int ndatasize = 1);
		//��������ṹ��ʸ��������
		void gen_vec_sph_uv1(char type/*u����v��*/, char sigma/*��ż*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, const void* cart, void* res_cart, int ndatasize = 1);
		void gen_vec_sph_uv2(char type/*u����v��*/, char sigma/*��ż*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, const void* cart, void* res_cart, int ndatasize = 1);
		//��ָ����ʽ��ʸ�����沨
		void gen_vec_sph_uv1_cpx(char type/*u����v��*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, const void* cart, void* res_cart, int ndatasize = 1);
		void gen_vec_sph_uv2_cpx(char type/*u����v��*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, const void* cart, void* res_cart, int ndatasize = 1);
		//u�β���Զ������
		void gen_far_vsw_u1(char type/*��ʽ����*/, char sigma/*��ż*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, const void* cart, void* res_cart, int ndatasize = 1);
		void gen_far_vsw_u2(char type/*��ʽ����*/, char sigma/*��ż*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, const void* cart, void* res_cart, int ndatasize = 1);
		//�������ʾ��ʸ�����沨
		void s_gen_vec_sph_uv1(char type/*u����v��*/, char sigma/*��ż*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, double r, double t, double p, std::complex<double>* res_r, std::complex<double>* res_t, std::complex<double>* res_p);
		void s_gen_vec_sph_uv2(char type/*u����v��*/, char sigma/*��ż*/, double k/*����*/, int ndegree/*��*/, int norder/*����*/, double r, double t, double p, std::complex<double>* res_r, std::complex<double>* res_t, std::complex<double>* res_p);
		double Wigner3j(int j1, int j2, int j3, int m1, int m2, int m3);
		double Wigner3j_precalc(int j1, int j2, int j3, int m1, int m2, int m3);
		void calc_wav_index(int degree, std::vector<std::vector<int>>& n_index);//��tslm����ӳ���n�����������˳����l,m,sig,tau
		void calc_wav_index_cpx(int degree, std::vector<std::vector<int>>& n_index);//��tlm����ӳ���n�����������˳����l,m,tau

	}

}


