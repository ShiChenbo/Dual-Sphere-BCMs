#include "Translation.h"
void calc_Trans_sp_zaxis_publish()
{
	fstream file;
	file.open(L"TransInfo.txt", ios::in);
	if (file.fail())
	{
		cout << "fail to open TransInfo.txt\n";
		return;
	}
	double r1, r2, deg;
	double f_start, f_end, f_step;
	file >> r1 >> r2 >> deg >> f_start >> f_end >> f_step;
	file.close();

	double d = r2 - r1;

	vector<double> beta;
	for (double freq = f_start; freq <= f_end; freq += f_step)
	{
		beta.push_back(2 * pi * freq * 1e9 / c);
	}
	int nsample = beta.size();
	vector<Matrix<double>> R1(nsample);
	vector<Matrix<double>> R2(nsample);
	vector < Matrix<complex<double>>> Y(nsample);

	calc_wav_TransR_sp_zaxis_simple(deg, beta, r1, R1, nsample);
	calc_wav_TransR_sp_zaxis_simple(deg, beta, r2, R2, nsample);
	calc_wav_TransP_sp_zaxis_simple(deg, beta, d, Y, nsample);
	auto num2str = [](double num)->std::string
	{
		LONG64 num2 = static_cast<LONG64>(num + 0.00001);
		std::stringstream ss;
		ss << num2;
		return ss.str();
	};
	for (int i = 0; i < nsample; i++)
	{
		std::string filenameR1 = "TranlationMatrix\\" + num2str((f_start + i * f_step) * 1e9) + ".000000R1.dat";
		std::string filenameR2 = "TranlationMatrix\\" + num2str((f_start + i * f_step) * 1e9) + ".000000R2.dat";
		std::string filenameYr = "TranlationMatrix\\" + num2str((f_start + i * f_step) * 1e9) + ".000000Yr.dat";
		std::string filenameYi = "TranlationMatrix\\" + num2str((f_start + i * f_step) * 1e9) + ".000000Yi.dat";
		R1[i].savebin(filenameR1.c_str());
		R2[i].savebin(filenameR2.c_str());
		Y[i].savebin(filenameYr.c_str(), filenameYi.c_str());
	}

	return;
}
int main()
{
	calc_Trans_sp_zaxis_publish();
	return 0;
}