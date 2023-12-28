#include "Matrix_Method.h"

template<typename Y>
Matrix<complex<double>> inv(Matrix<complex<Y>>& M)
{
	int nRows = M.size(1);
	int nCols = M.size(2);
	Matrix<complex<double>> tM(NULL, M.size(1), M.size(2));

	int i, j, k;
	int ret = 0;

	//! �����Ƿ���
	if (nRows != nCols)
	{
		cout << "����ľ�������Ƿ���\n";
		return tM;
	}
	int *is = new int[nRows];
	int *js = new int[nCols];
	double max(0);
	complex<Y> max_c;
	tM = M;
	for (k = 0; k < nRows; k++)
	{
		//step 1, ȫѡ��Ԫ
		max = 0;
		is[k] = k;
		js[k] = k;

		for (i = k; i < nRows; i++)
		{
			for (j = k; j < nCols; j++)
			{
				complex<double> tkey = tM(i + 1, j + 1);//ǿ��ת��
				if (max < abs(tkey))
				{
					max = abs(tkey);
					max_c = tkey;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (max == 0)
		{	// �������
			cout << "���������\n";
			return tM;
		}

		//����
		if (is[k] != k)
		{
			tM.swap_row(k + 1, is[k] + 1);
		}
		if (js[k] != k)
		{
			tM.swap_col(k + 1, js[k] + 1);
		}

		tM(k + 1, k + 1) = conj(tM(k + 1, k + 1)) / pow(abs(tM(k + 1, k + 1)), 2);//1/(a+bi)

		for (j = 0; j < nCols; j++)
		{
			if (j != k)
				tM(k + 1, j + 1) *= tM(k + 1, k + 1);
		}
		for (i = 0; i < nRows; i++)
		{
			if (i != k)
			{
				for (j = 0; j < nCols; j++)
				{
					if (j != k)
						tM(i + 1, j + 1) -= tM(i + 1, k + 1) * tM(k + 1, j + 1);
				}
			}
		}
		for (i = 0; i < nRows; i++)
		{
			if (i != k)
				tM(i + 1, k + 1) *= -tM(k + 1, k + 1);
		}

	}

	//�ָ�
	//���� row <-> is[k], column <-> js[k]
	//�ָ�ʱ��row <-> js[k], column <-> is[k]
	for (k = nRows - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			tM.swap_row(k + 1, js[k] + 1);
		}
		if (is[k] != k)
		{
			tM.swap_col(k + 1, is[k] + 1);
		}
	}
	return tM;
}

template<typename Y>
Matrix<double> inv(Matrix<Y>& M)
{
	int nRows = M.size(1);
	int nCols = M.size(2);
	Matrix<double> tM(NULL, nRows, nCols);


	int i, j, k;
	int ret = 0;

	//! �����Ƿ���
	if (nRows != nCols)
	{
		cout << "����ľ�������Ƿ���\n";
		return tM;
	}
	int *is = new int[nRows];
	int *js = new int[nCols];
	double max(0);
	tM = M;
	for (k = 0; k < nRows; k++)
	{
		//step 1, ȫѡ��Ԫ
		max = 0;
		is[k] = k;
		js[k] = k;

		for (i = k; i < nRows; i++)
		{
			for (j = k; j < nCols; j++)
			{
				Y tkey = tM(i + 1, j + 1);
				if (max < abs(tkey))
				{
					max = tkey;
					is[k] = i;
					js[k] = j;
				}
			}
		}

		if (max == 0)
		{	//! �������
			cout << "���������\n";
			return tM;
		}

		//����
		if (is[k] != k)
		{
			tM.swap_row(k + 1, is[k] + 1);
		}
		if (js[k] != k)
		{
			tM.swap_col(k + 1, js[k] + 1);
		}

		tM(k + 1, k + 1) = 1 / tM(k + 1, k + 1);

		for (j = 0; j < nCols; j++)
		{
			if (j != k)
				tM(k + 1, j + 1) *= tM(k + 1, k + 1);
		}
		for (i = 0; i < nRows; i++)
		{
			if (i != k)
			{
				for (j = 0; j < nCols; j++)
				{
					if (j != k)
						tM(i + 1, j + 1) -= tM(i + 1, k + 1) * tM(k + 1, j + 1);
				}
			}
		}
		for (i = 0; i < nRows; i++)
		{
			if (i != k)
				tM(i + 1, k + 1) *= -tM(k + 1, k + 1);
		}

	}

	//�ָ�
	//���� row <-> is[k], column <-> js[k]
	//�ָ�ʱ��row <-> js[k], column <-> is[k]
	for (k = nRows - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			tM.swap_row(k + 1, js[k] + 1);
		}
		if (is[k] != k)
		{
			tM.swap_col(k + 1, is[k] + 1);
		}
	}
	return tM;
}

template<typename Y>
Matrix<double> rref(Matrix<Y>& M, double tol)
{
	Matrix<double> tM;
	int nRows = M.size(1);
	int nCols = M.size(2);
	tM = M;
	int i, j, jj, isy(0), js;
	double max(0);
	int n_exe = 0;
	for (j = 1; j <= min(nCols, nRows + n_exe); j++)
	{
		js = j - n_exe;//������ʱ��Ӧ��Ҫ����
		max = tM(js, j); isy = js;
		for (i = js + 1; i <= nRows; i++)//Ѱ����Ԫ
		{
			double tkey = tM(i, j);
			if (abs(max) < abs(tkey))
			{
				max = tkey;
				isy = i;
			}
		}
		if (max == 0)//������
		{
			n_exe++;
			continue;
		}
		//========================
		//���г�ȥmax
		for (jj = j; jj <= nCols; jj++)
		{
			tM(isy, jj) /= max;
		}
		//========================
		tM.swap_row(isy, js);//������Ԫ�У�j)�������isy
		for (i = 1; i < js; i++)//��ȥ
		{
			double tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
		}
		for (i = js + 1; i <= nRows; i++)
		{
			double tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
		}
	}
	
	return tM;
}

template<typename Y>
Matrix<double> rref(Matrix<Y>& M, Matrix<double>& rM, double tol /*= 0.000001*/)
{
	Matrix<double> tM;
	int nRows = M.size(1);
	int nCols = M.size(2);
	tM = M;
	rM.eye(nRows);
	int i, j, jj, isy(0), js;
	double max(0);
	int n_exe = 0;
	for (j = 1; j <= min(nCols, nRows + n_exe); j++)
	{
		js = j - n_exe;//������ʱ��Ӧ��Ҫ����
		max = tM(js, j); isy = js;
		for (i = js + 1; i <= nRows; i++)//Ѱ����Ԫ
		{
			double tkey = tM(i, j);
			if (abs(max) < abs(tkey))
			{
				max = tkey;
				isy = i;
			}
		}
		if (max == 0)//������
		{
			n_exe++;
			continue;
		}
		//========================
		//���г�ȥmax
		for (jj = j; jj <= nCols; jj++)
		{
			tM(isy, jj) /= max;
		}
		for (jj = 1; jj <= rM.size(2); jj++)
		{
			rM(isy, jj) /= max;
		}
		//========================
		tM.swap_row(isy, js);//������Ԫ�У�j)�������isy
		rM.swap_row(isy, js);//������Ԫ�У�j)�������isy
		for (i = 1; i < js; i++)//��ȥ
		{
			double tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
			for (jj = 1; jj <= rM.size(2); jj++)
			{
				rM(i, jj) -= tkey*rM(js, jj);
				if (abs(rM(i, jj)) < tol)
				{
					rM(i, jj) = 0;
				}
			}
		}
		for (i = js + 1; i <= nRows; i++)
		{
			double tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
			for (jj = 1; jj <= rM.size(2); jj++)
			{
				rM(i, jj) -= tkey*rM(js, jj);
				if (abs(rM(i, jj)) < tol)
				{
					rM(i, jj) = 0;
				}
			}
		}
	}

	return tM;
}


template<typename Y>
Matrix<complex<double>> rref(Matrix<complex<Y>>& M, double tol)
{
	Matrix<complex<double>> tM;
	int nRows = M.size(1);
	int nCols = M.size(2);
	tM = M;
	int i, j, jj, isy(0), js;
	double max(0);
	complex<Y> max_c;
	int n_exe = 0;
	for (j = 1; j <= min(nCols, nRows + n_exe); j++)
	{
		js = j - n_exe;//������ʱ��Ӧ��Ҫ����
		max_c = tM(js, j); max = abs(max_c); isy = js;
		for (i = js + 1; i <= nRows; i++)//Ѱ����Ԫ
		{
			complex<double> tkey = tM(i, j);
			if (abs(max) < abs(tkey))
			{
				max = abs(tkey);
				max_c = tkey;
				isy = i;
			}
		}
		if (max == 0)//������
		{
			n_exe++;
			continue;
		}
		//========================
		//���г�ȥmax
		for (jj = j; jj <= nCols; jj++)
		{
			tM(isy, jj) *= (conj(max_c) / (pow(abs(max_c), 2)));//tM(isy, jj) /= max;
		}
		//========================
		tM.swap_row(isy, js);//������Ԫ�У�j)�������isy
		
		for (i = 1; i < js; i++)//��ȥ
		{
			complex<double> tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
		}
		for (i = js + 1; i <= nRows; i++)
		{
			complex<double> tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
		}
	}

	return tM;
}

template<typename Y>
Matrix<complex<double>> rref(Matrix<complex<Y>>& M, Matrix<complex<double>>& rM, double tol /*= 0.000001*/)
{
	Matrix<complex<double>> tM;
	int nRows = M.size(1);
	int nCols = M.size(2);
	tM = M;
	rM.eye(nRows);
	int i, j, jj, isy(0), js;
	double max(0);
	complex<Y> max_c;
	int n_exe = 0;
	for (j = 1; j <= min(nCols, nRows + n_exe); j++)
	{
		js = j - n_exe;//������ʱ��Ӧ��Ҫ����
		max_c = tM(js, j); max = abs(max_c); isy = js;
		for (i = js + 1; i <= nRows; i++)//Ѱ����Ԫ
		{
			complex<double> tkey = tM(i, j);
			if (abs(max) < abs(tkey))
			{
				max = abs(tkey);
				max_c = tkey;
				isy = i;
			}
		}
		if (max == 0)//������
		{
			n_exe++;
			continue;
		}
		//========================
		//���г�ȥmax
		for (jj = j; jj <= nCols; jj++)
		{
			tM(isy, jj) *= (conj(max_c) / (pow(abs(max_c), 2)));//tM(isy, jj) /= max;
		}
		for (jj = 1; jj <= rM.size(2); jj++)
		{
			rM(isy, jj) *= (conj(max_c) / (pow(abs(max_c), 2)));
		}
		//========================
		tM.swap_row(isy, js);//������Ԫ�У�j)�������isy
		rM.swap_row(isy, js);//������Ԫ�У�j)�������isy
		for (i = 1; i < js; i++)//��ȥ
		{
			complex<double> tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
			for (jj = 1; jj <= rM.size(2); jj++)
			{
				rM(i, jj) -= tkey*rM(js, jj);
				if (abs(rM(i, jj)) < tol)
				{
					rM(i, jj) = 0;
				}
			}
		}
		for (i = js + 1; i <= nRows; i++)
		{
			complex<double> tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
			for (jj = 1; jj <= rM.size(2); jj++)
			{
				rM(i, jj) -= tkey*rM(js, jj);
				if (abs(rM(i, jj)) < tol)
				{
					rM(i, jj) = 0;
				}
			}
		}
	}

	return tM;
}

template<typename Y>
int rank_m(Matrix<Y> M, double  tol)
{
	int rank(0);
	int nRows = M.size(1);
	int nCols = M.size(2);
	int i, j, jj, isy(0), js;
	double max(0);
	int n_exe = 0;
	for (j = 1; j <= min(nCols, nRows + n_exe); j++)
	{
		js = j - n_exe;//������ʱ��Ӧ��Ҫ����
		max = M(js, j); isy = js;
		for (i = js + 1; i <= nRows; i++)//Ѱ����Ԫ
		{
			double tkey = M(i, j);
			if (abs(max) < abs(tkey))
			{
				max = tkey;
				isy = i;
			}
		}
		if (max == 0)//������
		{
			n_exe++;
			continue;
		}
		//========================
		//���г�ȥmax
		for (jj = j; jj <= nCols; jj++)
		{
			M(isy, jj) /= max;
		}
		//========================
		rank++;
		//========================
		M.swap_row(isy, js);//������Ԫ�У�j)�������isy
		for (i = 1; i < js; i++)//��ȥ
		{
			double tkey = M(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				M(i, jj) -= tkey*M(js, jj);
				if (abs(M(i, jj)) < tol)
				{
					M(i, jj) = 0;
				}
			}
		}
		for (i = js + 1; i <= nRows; i++)
		{
			double tkey = M(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				M(i, jj) -= tkey*M(js, jj);
				if (abs(M(i, jj)) < tol)
				{
					M(i, jj) = 0;
				}
			}
		}
	}

	return rank;
}

template<typename Y>
int rank_m(Matrix<complex<Y>> M, double tol /*= 0.000001*/)
{
	int rank(0);
	int nRows = M.size(1);
	int nCols = M.size(2);
	int i, j, jj, isy(0), js;
	double max(0);
	complex<Y> max_c;
	int n_exe = 0;
	for (j = 1; j <= min(nCols, nRows + n_exe); j++)
	{
		js = j - n_exe;//������ʱ��Ӧ��Ҫ����
		max_c = M(js, j); max = abs(max_c); isy = js;
		for (i = js + 1; i <= nRows; i++)//Ѱ����Ԫ
		{
			complex<double> tkey = M(i, j);
			if (abs(max) < abs(tkey))
			{
				max = abs(tkey);
				max_c = tkey;
				isy = i;
			}
		}
		if (max == 0)//������
		{
			n_exe++;
			continue;
		}
		//========================
		rank++;
		//========================
		//���г�ȥmax
		for (jj = j; jj <= nCols; jj++)
		{
			M(isy, jj) *= (conj(max_c) / (pow(abs(max_c), 2)));//tM(isy, jj) /= max;
		}
		//========================
		M.swap_row(isy, js);//������Ԫ�У�j)�������isy

		for (i = 1; i < js; i++)//��ȥ
		{
			complex<double> tkey = M(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				M(i, jj) -= tkey*M(js, jj);
				if (abs(M(i, jj)) < tol)
				{
					M(i, jj) = 0;
				}
			}
		}
		for (i = js + 1; i <= nRows; i++)
		{
			complex<double> tkey = M(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				M(i, jj) -= tkey*M(js, jj);
				if (abs(M(i, jj)) < tol)
				{
					M(i, jj) = 0;
				}
			}
		}
	}

	return rank;
}

template<typename Y>
Matrix<double> null(Matrix<Y>& M, double  tol)
{
	Matrix<double> tM;
	int nRows = M.size(1);
	int nCols = M.size(2);
	tM = M;
	int i, j, jj, isy(0), js;
	double max(0);
	int n_exe(0), rank(0);
	for (j = 1; j <= min(nCols, nRows + n_exe); j++)
	{
		js = j - n_exe;//������ʱ��Ӧ��Ҫ����
		max = tM(js, j); isy = js;
		for (i = js + 1; i <= nRows; i++)//Ѱ����Ԫ
		{
			double tkey = tM(i, j);
			if (abs(max) < abs(tkey))
			{
				max = tkey;
				isy = i;
			}
		}
		if (max == 0)//������
		{
			n_exe++;
			continue;
		}
		//========================
		rank++;
		//���г�ȥmax
		for (jj = j; jj <= nCols; jj++)
		{
			tM(isy, jj) /= max;
		}
		//========================
		tM.swap_row(isy, js);//������Ԫ�У�j)�������isy
		for (i = 1; i < js; i++)//��ȥ
		{
			double tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
		}
		for (i = js + 1; i <= nRows; i++)
		{
			double tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
		}
	}
	//���������ϵ
	int* pseries = new int[nCols + 1];//�࿪��һ���ڴ棬ʹ������1��ʼ
	for (i = 1; i <= nCols; i++)//��ʼ���к�
	{
		pseries[i] = i;
	}
	for (i = 1; i <= rank; i++)//�鿴�Խ�Ԫ�Ƿ��0
	{
		if (tM(i, pseries[i]) == 0)
		{
			for (j = pseries[i + 1]; j <= nCols; j++)
			{
				if (tM(i, j) != 0)
				{
					int tindex = pseries[i];
					pseries[i] = j;
					pseries[j] = tindex;//����������
					break;
				}
			}
		}
	}
	int dim = nCols - rank;//������ϵ�ĸ�������δ֪������n-rank
	Matrix<double> nuM(NULL, nCols, dim);
	if (dim == 0)
	{
		cout << "��η���ֻ��0��\n";
		return nuM;
	}
	for (i = 1; i <= dim; i++)
	{
		nuM(i + rank, i) = 1;//�����һ����λ��
	}
	for (i = 1; i <= rank; i++)
	{
		cout << i;
		for (j = 1; j <= dim; j++)
		{
			nuM(i, j) = -tM(i, pseries[j + rank]);
		}
	}
	j = 1;
	do 
	{
		if (pseries[j] != j)
		{
			nuM.swap_row(pseries[pseries[j]], pseries[j]);//������֮��û�����ϵ
			i = pseries[pseries[j]];
			pseries[pseries[j]] = pseries[j];
			pseries[j] = i;
		}
		else
		{
			j++;
		}
	} while (j <= rank);

	return nuM;
}

template<typename Y>
Matrix<complex<double>> null(Matrix<complex<Y>>& M, double tol)
{
	Matrix<complex<double>> tM;
	int nRows = M.size(1);
	int nCols = M.size(2);
	tM = M;
	int i, j, jj, isy(0), js;
	double max(0);
	complex<Y> max_c;
	int n_exe(0), rank(0);
	for (j = 1; j <= min(nCols, nRows + n_exe); j++)
	{
		js = j - n_exe;//������ʱ��Ӧ��Ҫ����
		max_c = tM(js, j); max = abs(max_c); isy = js;
		for (i = js + 1; i <= nRows; i++)//Ѱ����Ԫ
		{
			complex<double> tkey = tM(i, j);
			if (abs(max) < abs(tkey))
			{
				max = abs(tkey);
				max_c = tkey;
				isy = i;
			}
		}
		if (max == 0)//������
		{
			n_exe++;
			continue;
		}
		//========================
		rank++;
		//���г�ȥmax
		for (jj = j; jj <= nCols; jj++)
		{
			tM(isy, jj) *= (conj(max_c) / (pow(abs(max_c), 2)));//tM(isy, jj) /= max;
		}
		//========================
		tM.swap_row(isy, js);//������Ԫ�У�j)�������isy

		for (i = 1; i < js; i++)//��ȥ
		{
			complex<double> tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
		}
		for (i = js + 1; i <= nRows; i++)
		{
			complex<double> tkey = tM(i, j);//��¼��Ԫ����Ϊ��ȥ��ϵ����
			for (jj = j; jj <= nCols; jj++)
			{
				tM(i, jj) -= tkey*tM(js, jj);
				if (abs(tM(i, jj)) < tol)
				{
					tM(i, jj) = 0;
				}
			}
		}
	}
	//���������ϵ
	int* pseries = new int[nCols + 1];//�࿪��һ���ڴ棬ʹ������1��ʼ
	for (i = 1; i <= nCols; i++)//��ʼ���к�
	{
		pseries[i] = i;
	}
	for (i = 1; i <= rank; i++)//�鿴�Խ�Ԫ�Ƿ��0
	{
		if (tM(i, pseries[i]) == complex < double>(0, 0))
		{
			for (j = pseries[i + 1]; j <= nCols; j++)
			{
				if (tM(i, j) != complex < double>(0, 0))
				{
					int tindex = pseries[i];
					pseries[i] = j;
					pseries[j] = tindex;//����������
					break;
				}
			}
		}
	}
	int dim = nCols - rank;//������ϵ�ĸ�������δ֪������n-rank
	Matrix<complex<double>> nuM(NULL, nCols, dim);
	if (dim == 0)
	{
		cout << "��η���ֻ��0��\n";
		return nuM;
	}
	for (i = 1; i <= dim; i++)
	{
		nuM(i + rank, i) = 1;//�����һ����λ��
	}
	for (i = 1; i <= rank; i++)
	{
		cout << i;
		for (j = 1; j <= dim; j++)
		{
			nuM(i, j) = -tM(i, pseries[j + rank]);
		}
	}
	j = 1;
	do
	{
		if (pseries[j] != j)
		{
			nuM.swap_row(pseries[pseries[j]], pseries[j]);//������֮��û�����ϵ
			i = pseries[pseries[j]];
			pseries[pseries[j]] = pseries[j];
			pseries[j] = i;
		}
		else
		{
			j++;
		}
	} while (j <= rank);

	return nuM;
}

template<typename Y>
Matrix<Y> les_solve(Matrix<Y>& A, Matrix<Y>& b,double tol)
{
	Matrix<Y> Augmat = { &A, &b };//�������
	Matrix<Y> sol = rref(Augmat, tol);
	return sol.col(Augmat.size(2));
}
