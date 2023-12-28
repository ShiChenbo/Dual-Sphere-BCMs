#include "Translation.h"
#include "MathBasic.h"
void calc_wav_TransR_sp_zaxis_simple(int degree, vector<double>& beta, double d, vector<Matrix<double>>& mat, int samplenum)
{
	namespace Mvs = MathBasic::vec_sph_wav;//命名空间的别名
	//定义匿名函数，用于计算系数
	std::vector<std::vector<int>> n_index;
	Mvs::calc_wav_index(degree, n_index);
	int Nsize = n_index.size();
	vector<double> kd(samplenum);
	for (int ikd = 0; ikd < samplenum; ikd++)
	{
		mat[ikd].resize(Nsize, Nsize);
		kd[ikd] = beta[ikd] * d;
	}

	//预先计算不同频率的bessel函数
	Matrix<double>bessel_kd(NULL, 2 * degree + 1, kd.size());
	for (int ikd = 0; ikd < kd.size(); ikd++)
	{
		for (int lam = 0; lam <= 2 * degree; lam++)
		{
			bessel_kd(lam + 1, ikd + 1) = MathBasic::sph_bessel(lam, kd[ikd]);
		}
	}

	//定义匿名函数
	auto coefC = [&bessel_kd, &degree, &samplenum](vector<double>& res, double l, double m, double l_)
	{
		using namespace std;
		double em = m == 0 ? 1 : 2;
		double t00 = 0.5 * em / 2;
		std::fill(res.begin(), res.end(), 0);

		for (int lam = fabs(l - l_); lam <= l + l_; lam++)
		{
			if (std::fmod(l_ - l + lam, 2) != 0)
			{
				continue;
			}
			double t0 = Mvs::Wigner3j(l, l_, lam, 0, 0, 0) * Mvs::Wigner3j(l, l_, lam, m, -m, 0);
			if (t0 == 0)
			{
				continue;
			}
			double t1 = t00 * pow(-1, (l_ - l + lam) / 2) * (2 * lam + 1) *
				sqrt((2 * l + 1) * (2 * l_ + 1) / (l * (l + 1) * l_ * (l_ + 1)))
				* t0
				* (l * (l + 1) + l_ * (l_ + 1) - lam * (lam + 1));
			for (int ikd = 0; ikd < samplenum; ikd++)
			{
				res[ikd] += t1 * bessel_kd(lam + 1, ikd + 1);
			}
		}
	};
	auto coefD = [&bessel_kd, &degree, &samplenum, &kd](vector<double>& res, double l, double m, double l_)
	{
		using namespace std;
		double t00 = -m;
		std::fill(res.begin(), res.end(), 0);

		for (int lam = fabs(l - l_); lam <= l + l_; lam++)
		{
			if (std::fmod(l_ - l + lam, 2) != 0)
			{
				continue;
			}
			double t0 = Mvs::Wigner3j(l, l_, lam, 0, 0, 0) * Mvs::Wigner3j(l, l_, lam, m, -m, 0);
			if (t0 == 0)
			{
				continue;
			}
			double t1 = t00 * pow(-1, (l_ - l + lam + 2) / 2) * (2 * lam + 1) *
				sqrt((2 * l + 1) * (2 * l_ + 1) / (l * (l + 1) * l_ * (l_ + 1)))
				* t0;
			for (int ikd = 0; ikd < samplenum; ikd++)
			{
				res[ikd] += t1 * kd[ikd] * bessel_kd(lam + 1, ikd + 1);
			}
		}
	};

	//计算平移系数
	vector<double> coC(samplenum), coD(samplenum);
	for (int alpha = 0; alpha < Nsize; alpha += 2)//始终访问tau=1的位置
	{
		int id = n_index[alpha][0];
		int io = n_index[alpha][1];
		int is = n_index[alpha][2];
		int it = n_index[alpha][3];
		for (int alpha_ = 0; alpha_ < Nsize; alpha_++)
		{
			int id_ = n_index[alpha_][0];
			int io_ = n_index[alpha_][1];
			int is_ = n_index[alpha_][2];
			int it_ = n_index[alpha_][3];
			if (io_ != io)
			{
				continue;
			}
			double tmp = 0;
			if (it_ == 1)
			{
				if (is == is_)
				{
					coefC(coC, id, io, id_);
					for (int ikd = 0; ikd < kd.size(); ikd++)
					{
						double pm = io == 0 ? 1 + pow(-1, is) : pow(-1, io);//m=0的情形稍有一些不一样
						auto tmp = pm * coC[ikd];
						mat[ikd](alpha + 1, alpha_ + 1) = tmp;
						mat[ikd](alpha + 2, alpha_ + 2) = tmp;
					}
				}
			}
			else if (it_ == 2)
			{
				if (is != is_)
				{
					coefD(coD, id, io, id_);
					for (int ikd = 0; ikd < kd.size(); ikd++)
					{
						if (io != 0)//m=0时直接是0
						{
							double dm = pow(-1, io + is);
							auto tmp = dm * coD[ikd];
							mat[ikd](alpha + 1, alpha_ + 1) = tmp;
							mat[ikd](alpha + 2, alpha_) = tmp;
						}
					}
				}
			}
		}
	}
}

void calc_wav_TransP_sp_zaxis_simple(int degree, vector<double>& beta, double d, vector<Matrix<complex<double>>>& mat, int samplenum)
{
	namespace Mvs = MathBasic::vec_sph_wav;//命名空间的别名
	//定义匿名函数，用于计算系数
	std::vector<std::vector<int>> n_index;
	Mvs::calc_wav_index(degree, n_index);
	int Nsize = n_index.size();
	vector<double> kd(samplenum);
	for (int ikd = 0; ikd < samplenum; ikd++)
	{
		mat[ikd].resize(Nsize, Nsize);
		kd[ikd] = beta[ikd] * d;
	}

	//预先计算不同频率的bessel函数
	Matrix<complex<double>>bessel_kd(NULL, 2 * degree + 1, kd.size());
	for (int ikd = 0; ikd < kd.size(); ikd++)
	{
		for (int lam = 0; lam <= 2 * degree; lam++)
		{
			bessel_kd(lam + 1, ikd + 1) = MathBasic::sph_hankel(lam, 2, kd[ikd]);
		}
	}

	//定义匿名函数
	auto coefC = [&bessel_kd, &degree, &samplenum](vector<complex<double>>& res, double l, double m, double l_)
	{
		using namespace std;
		double em = m == 0 ? 1 : 2;
		double t00 = 0.5 * em / 2;
		std::fill(res.begin(), res.end(), 0);

		for (int lam = fabs(l - l_); lam <= l + l_; lam++)
		{
			if (std::fmod(l_ - l + lam, 2) != 0)
			{
				continue;
			}
			double t0 = Mvs::Wigner3j(l, l_, lam, 0, 0, 0) * Mvs::Wigner3j(l, l_, lam, m, -m, 0);
			if (t0 == 0)
			{
				continue;
			}
			double t1 = t00 * pow(-1, (l_ - l + lam) / 2) * (2 * lam + 1) *
				sqrt((2 * l + 1) * (2 * l_ + 1) / (l * (l + 1) * l_ * (l_ + 1)))
				* t0
				* (l * (l + 1) + l_ * (l_ + 1) - lam * (lam + 1));
			for (int ikd = 0; ikd < samplenum; ikd++)
			{
				res[ikd] += t1 * bessel_kd(lam + 1, ikd + 1);
			}
		}
	};
	auto coefD = [&bessel_kd, &degree, &samplenum, &kd](vector<complex<double>>& res, double l, double m, double l_)
	{
		using namespace std;
		double t00 = -m;
		std::fill(res.begin(), res.end(), 0);

		for (int lam = fabs(l - l_); lam <= l + l_; lam++)
		{
			if (std::fmod(l_ - l + lam, 2) != 0)
			{
				continue;
			}
			double t0 = Mvs::Wigner3j(l, l_, lam, 0, 0, 0) * Mvs::Wigner3j(l, l_, lam, m, -m, 0);
			if (t0 == 0)
			{
				continue;
			}
			double t1 = t00 * pow(-1, (l_ - l + lam + 2) / 2) * (2 * lam + 1) *
				sqrt((2 * l + 1) * (2 * l_ + 1) / (l * (l + 1) * l_ * (l_ + 1)))
				* t0;
			for (int ikd = 0; ikd < samplenum; ikd++)
			{
				res[ikd] += t1 * kd[ikd] * bessel_kd(lam + 1, ikd + 1);
			}
		}
	};

	//计算平移系数
	vector<complex<double>> coC(samplenum), coD(samplenum);
	for (int alpha = 0; alpha < Nsize; alpha += 2)//始终访问tau=1的位置
	{
		int id = n_index[alpha][0];
		int io = n_index[alpha][1];
		int is = n_index[alpha][2];
		int it = n_index[alpha][3];
		for (int alpha_ = 0; alpha_ < Nsize; alpha_++)
		{
			int id_ = n_index[alpha_][0];
			int io_ = n_index[alpha_][1];
			int is_ = n_index[alpha_][2];
			int it_ = n_index[alpha_][3];
			if (io_ != io)
			{
				continue;
			}
			double tmp = 0;
			if (it_ == 1)
			{
				if (is == is_)
				{
					coefC(coC, id, io, id_);
					for (int ikd = 0; ikd < kd.size(); ikd++)
					{
						double cm = io == 0 ? 1 + pow(-1, is) : pow(-1, io);
						auto tmp = cm * coC[ikd];
						mat[ikd](alpha + 1, alpha_ + 1) = tmp;
						mat[ikd](alpha + 2, alpha_ + 2) = tmp;
					}
				}
			}
			else if (it_ == 2)
			{
				if (is != is_)
				{
					coefD(coD, id, io, id_);
					for (int ikd = 0; ikd < kd.size(); ikd++)
					{
						if (io != 0)//io=0的情形直接是0
						{
							double dm = pow(-1, io + is);
							auto tmp = dm * coD[ikd];
							mat[ikd](alpha + 1, alpha_ + 1) = tmp;
							mat[ikd](alpha + 2, alpha_) = tmp;
						}
					}
				}
			}
		}
	}
}

