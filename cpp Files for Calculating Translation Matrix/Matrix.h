#pragma once
#include <iomanip>
#include <vector>
#include <set>
#include <initializer_list>
#include <complex>
#include <algorithm>
#include <fstream>
using namespace std;

template <typename T>//ֻ��ʵ������ϡ�軯����������ϡ�軯��������
class Mat_Sparse
{
public:
	T* data;
	int* r_index;
	int* c_index;
	int m_nrows;
	int m_ncols;
	int m_nsize;
public:
	Mat_Sparse(int nrows, int ncols, int nsize)
	{
		m_nrows = nrows;
		m_ncols = ncols;
		m_nsize = nsize;
		data = nullptr;
		data = new T[m_nsize];
		r_index = new int[m_nsize];
		c_index = new int[m_nsize];
		memset(data, 0, sizeof(T)*m_nsize);
		memset(r_index, -1, sizeof(int)*m_nsize);
		memset(c_index, -1, sizeof(int)*m_nsize);
	}
	Mat_Sparse()
	{
		m_nrows = 0;
		m_ncols = 0;
		data = nullptr;		
		r_index = nullptr;
		c_index = nullptr;
	}
	~Mat_Sparse()
	{
		if (data)
		{
			delete[] data;
			data = nullptr;
		}		
		if (m_nsize != -1 && r_index)
		{
			delete[] r_index;
			r_index = nullptr;
		}
		if (m_nsize != -1 && c_index)
		{
			delete[] c_index;
			c_index = nullptr;
		}
		m_nrows = 0;
		m_ncols = 0;
		m_nsize = 0;
	}	
	void Bindindex(Mat_Sparse<T>& mat_spa_comp)//ϡ������������λ����ȫһ�£����ݲ�ͬ��
	{
		r_index = mat_spa_comp.r_index;
		c_index = mat_spa_comp.c_index;
		m_nrows = mat_spa_comp.m_nrows;
		m_ncols = mat_spa_comp.m_ncols;
		m_nsize = -1;
		if (data)
		{
			delete[] data;
		}
		data = new T[mat_spa_comp.m_nsize];
		memset(data, 0, sizeof(T)* mat_spa_comp.m_nsize);
	}
	void resize(int nrows, int ncols,int nsize)
	{
		if (m_nsize == -1)
		{
			cout << "����ϡ�����������ĳߴ�\n";
			return;
		}
		if (data)
		{
			delete[] data;
			delete[] r_index;
			delete[] c_index;
		}
		m_nrows = nrows;
		m_ncols = ncols;
		m_nsize = nsize;
		data = new T[m_nsize];
		r_index = new int[m_nsize];
		c_index = new int[m_nsize];
		memset(data, 0, sizeof(T)*m_nsize);
		memset(r_index, -1, sizeof(int)*m_nsize);
		memset(c_index, -1, sizeof(int)*m_nsize);
	}
};

template < typename T> 
class Matrix
{
private:
	template < typename T >
	class Matrix_p 
	{
	public:
		bool m_isatisfy = false;
		int m_nRows;
		int m_nCols;
		Matrix<T>* pM;//Դͷ���󣬵���Դͷthis
		int* pRelativePos;//��¼������Դͷ�����λ��

		Matrix_p(char* c, int nRows, int nCols, Matrix<T>* ptM)
		{
			m_nRows = nRows;
			m_nCols = nCols;
			pM = ptM;//�󶨾���
			if (m_nRows && m_nCols)
			{
				int length = m_nRows * m_nCols;
				pRelativePos = new int[length];
				memset(pRelativePos, 0, sizeof(int)*length);
			}
		}
		Matrix_p(const Matrix_p<T>& M)
		{
			m_nCols = M.m_nCols;
			m_nRows = M.m_nRows;
			m_isatisfy = M.m_isatisfy;
			int length = m_nRows*m_nCols;

			pRelativePos = new int[length];
			pM = M.pM;
			for (int i = 0; i < length; i++)
			{
				pRelativePos[i] = M.pRelativePos[i];
			}
		}
		~Matrix_p()
		{
			m_nCols = 0;
			m_nRows = 0;
			m_isatisfy = false;
			delete[]pRelativePos;
			pRelativePos = nullptr;
		}
		T& operator()(const int& i)
		{
			if (i > this->m_nRows*this->m_nCols || i < 1)
			{
				cout << "Matrix_p����Խ��\n";
				T tkey;
				return tkey;
			}
			return pRelativePos[i - 1];
		}
		void operator=(Matrix<T>& M)
		{
			//Ҫ���ܳ���һ��
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return;
			}
			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) = M(i + 1);
			}
		}		
		void operator=(const initializer_list<initializer_list<T>> il)
		{
			Matrix<T> M(il);
			*this = M;
		}
		void operator=(const initializer_list<T> il)
		{
			Matrix<T> M(il);
			*this = M;
		}
		void operator=(const T ele)//����һ��Ԫ�ص�ʱ�򣬱�ʾ�Ըþ���ָʾ��λ�õ�Ԫ�ض����������
		{
			int length = m_nRows*m_nCols;
			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) = ele;
			}
		}
		template<typename Y>
		Matrix operator+(Matrix<Y>& M)
		{
			Matrix<T> tM = *this;			
			if ((is_same<complex<double>, T>::value && is_same<double, Y>::value) || is_same<T, Y>::value)//����Ҫƥ�䣬���ݸ�����double
			{
				int tRows = M.size(1);
				int tCols = M.size(2);
				int length = m_nRows*m_nCols;
				if (tRows*tCols != length)
				{
					cout << "����ά�Ȳ�һ��\n";
					return tM;
				}
				for (int i = 1; i <= length; i++)
				{
					tM(i) += M(i);
				}
			}
			else
			{
				cout << "�������Ͳ�ƥ��\n";
			}
			//Ҫ���ܳ���һ��
			return tM;
		}		
		Matrix operator+(Matrix_p<T>& pMp)
		{
			//Ҫ���ܳ���һ��
			int tRows = pMp.m_nRows;
			int tCols = pMp.m_nCols;
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return tM;
			}
			for (int i = 1; i <= length; i++)
			{
				tM(i) += (*pMp.pM)(pMp(i));
			}
			return tM;
		}
		Matrix operator+(const T& tkey)
		{
			//Ҫ���ܳ���һ��
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			for (int i = 1; i <= length; i++)
			{
				tM(i) += tkey;
			}
			return tM;
		}
		Matrix_p& operator+=(Matrix_p<T>& pMp)
		{
			//Ҫ���ܳ���һ��
			int tRows = pMp.m_nRows;
			int tCols = pMp.m_nCols;
			int length = m_nRows*m_nCols;
			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return *this;
			}
			for (int i = 1; i <= length; i++)
			{
				(*(this->pM))((*this)(i)) += (*pMp.pM)(pMp(i));
			}
			return *this;
		}
		Matrix_p& operator+=(Matrix<T>& M)
		{
			//Ҫ���ܳ���һ��
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return *this;
			}
			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) += M(i + 1);
			}
			return *this;
		}
		Matrix_p& operator+=(const T& tkey)
		{
			int length = m_nRows*m_nCols;
			
			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) += tkey;
			}
			return *this;
		}
		Matrix operator-(Matrix<T>& M)
		{
			//Ҫ���ܳ���һ��
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return tM;
			}
			for (int i = 1; i <= length; i++)
			{
				tM(i) -= M(i);
			}
			return tM;
		}
		Matrix operator-(const T& tkey)
		{
			//Ҫ���ܳ���һ��
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			for (int i = 1; i <= length; i++)
			{
				tM(i) -= tkey;
			}
			return tM;
		}
		Matrix operator-(Matrix_p<T>& pMp)
		{
			//Ҫ���ܳ���һ��
			int tRows = pMp.m_nRows;
			int tCols = pMp.m_nCols;
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return tM;
			}
			for (int i = 1; i <= length; i++)
			{
				tM(i) -= (*pMp.pM)(pMp(i));
			}
			return tM;
		}
		Matrix_p& operator-=(Matrix_p<T>& pMp)
		{
			//Ҫ���ܳ���һ��
			int tRows = pMp.m_nRows;
			int tCols = pMp.m_nCols;
			int length = m_nRows*m_nCols;
			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return *this;
			}
			for (int i = 1; i <= length; i++)
			{
				(*(this->pM))((*this)(i)) -= (*pMp.pM)(pMp(i));
			}
			return *this;
		}
		Matrix_p& operator-=(Matrix<T>& M)
		{
			//Ҫ���ܳ���һ��
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return *this;
			}
			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) -= M(i + 1);
			}
			return *this;
		}
		Matrix_p& operator-=(const T& tkey)
		{
			int length = m_nRows*m_nCols;

			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) -= tkey;
			}
			return *this;
		}
		Matrix operator*(const T& tkey)
		{
			//Ҫ���ܳ���һ��
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			for (int i = 1; i <= length; i++)
			{
				tM(i) *= tkey;
			}
			return tM;
		}
		Matrix operator*(Matrix<T>& M)
		{
			//Ҫ���ܳ���һ��
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return tM;
			}
			for (int i = 1; i <= length; i++)
			{
				tM(i) *= M(i);
			}
			return tM;
		}
		Matrix_p& operator*=(Matrix<T>& M)
		{
			//Ҫ���ܳ���һ��
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return *this;
			}
			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) *= M(i + 1);
			}
			return *this;
		}
		Matrix_p& operator*=(const T& tkey)
		{
			int length = m_nRows*m_nCols;

			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) *= tkey;
			}
			return *this;
		}
		Matrix operator/(const T& tkey)
		{
			//Ҫ���ܳ���һ��
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			for (int i = 1; i <= length; i++)
			{
				tM(i) /= tkey;
			}
			return tM;
		}
		Matrix operator/(Matrix<T>& M)
		{
			//Ҫ���ܳ���һ��
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return tM;
			}
			for (int i = 1; i <= length; i++)
			{
				tM(i) /= M(i);
			}
			return tM;
		}
		Matrix_p& operator/=(const T& tkey)
		{
			int length = m_nRows*m_nCols;

			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) /= tkey;
			}
			return *this;
		}
		Matrix_p& operator/=(Matrix<T>& M)
		{
			//Ҫ���ܳ���һ��
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "����ά�Ȳ�һ��\n";
				return *this;
			}
			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) /= M(i + 1);
			}
			return *this;
		}
		friend ostream& operator<<(ostream& out, Matrix_p &pMp)
		{
			for (int i = 0; i < pMp.m_nRows; i++)
			{
				for (int j = 0; j < pMp.m_nCols; j++)
				{
					if (is_same<double, T>::value || is_same<float, T>::value)
					{
						printf("% .4lf\t", pMp.pM->operator()(pMp.pRelativePos[i + j*pMp.m_nRows]));
					}
					else
					{
						out << pMp.pM->operator()(pMp.pRelativePos[i + j*pMp.m_nRows]) << "\t";
					}
				}
				out << endl;
			}
			return out;
		}
		explicit operator bool() const//�����������Ԫ�أ���ô����false�����򷵻�true explicit��ʾֻ��ǿ��ת�������������
		{
			return m_isatisfy;
		}
		void zero()//��Ԫ������
		{
			int length = m_nRows*m_nCols;
			for (int i = 0; i < length; i++)
			{
				(*pM)(pRelativePos[i]) = 0;
			}
		}
	};

protected:
	int m_nRows;
	int m_nCols;
	T* m_pEntry = nullptr;
	bool m_blflag;//���ڱ�Ǹþ����Ƿ���һ���󶨾��󣬴���Ǳ�ڵ��ڴ��ͷŵ�����
public:
	Matrix(char c = NULL, int nRows = 0, int nCols = 0);//����1����ʽ���Ĳ���
	Matrix(const initializer_list<initializer_list<T>>& il);//���й���
	Matrix(const initializer_list<T>& il);//�й���
	Matrix(const T& ele);
	Matrix(const Matrix& M);//��������

	Matrix(const Matrix_p<T>& pMp);
	Matrix(const initializer_list<Matrix<T>*>& lM);
	Matrix(const initializer_list<initializer_list<Matrix<T>*>>& l2M);

	friend ostream& operator<<(ostream& out, Matrix &M)
	{
		for (int i = 0; i < M.m_nRows; i++)
		{
			for (int j = 0; j < M.m_nCols; j++)
			{
				if (is_same<double, T>::value || is_same<float, T>::value)
				{
					printf("% .4lf\t", M.m_pEntry[i + j*M.m_nRows]);
				}
				else if (is_same<complex<double>, T>::value || is_same<complex<float>, T>::value)
				{
					printf("%.4e + j%.4e\t", ((complex<double>*) M.m_pEntry)[i + j*M.m_nRows].real(), ((complex<double>*) M.m_pEntry)[i + j*M.m_nRows].imag());
				}
				else
				{
					out << M.m_pEntry[i + j*M.m_nRows] << "\t";
				}
			}
			out << endl;
			if (is_same<complex<double>, T>::value || is_same<complex<float>, T>::value)
			{
				out << endl;//�����һ��
			}
		}
		return out;
	}
	Matrix& operator=(const Matrix &M)
	{
		m_nRows = M.m_nRows;
		m_nCols = M.m_nCols;
		if (m_pEntry)
		{
			delete[] m_pEntry;
			m_pEntry = nullptr;
		}
		m_pEntry = new T[m_nCols*m_nRows];
		for (int j = 0; j < m_nCols*m_nRows; j++)
		{
			m_pEntry[j] = M.m_pEntry[j];
		}
		return *this;
	}
	Matrix& operator=(const initializer_list<initializer_list<T>>& il)
	{
		m_nRows = il.size();
		auto ibeg = il.begin();
		if (m_nRows)
			m_nCols = (ibeg + 0)->size();
		//����ά�Ȳ�һ�µĴ���
		for (int i = 0; i < m_nRows; i++)
		{
			if ((ibeg + i)->size() != m_nCols)
			{
				clear();
				cout << "����ά�Ȳ�һ��\n";
				break;
			}
		}

		m_pEntry = new T[m_nRows * m_nCols];

		for (int i = 0; i < m_nRows; i++)
		{
			for (int j = 0; j < m_nCols; j++)
			{
				m_pEntry[i + m_nRows*j] = *((ibeg + i)->begin() + j);
			}
		}
		return *this;
	}
	Matrix& operator=(const initializer_list<T>& il)
	{
		m_nCols = il.size();
		m_nRows = 1;

		m_pEntry = new T[m_nCols];

		for (int j = 0; j < m_nCols; j++)
		{
			m_pEntry[j] = *(il.begin() + j);
		}
		return *this;
	}
	Matrix& operator=(const T& ele)
	{
		m_nRows = 1;
		m_nCols = 1;
		if (m_pEntry)
		{
			delete[] m_pEntry;
			m_pEntry = nullptr;
		}
		m_pEntry = new T;
		*m_pEntry = ele;
		return *this;
	}
	Matrix& operator=(const Matrix_p<T>& pMp)//��ָ��ָʾ��λ�õ�Ԫ��ֵ������Ӧ����
	{
		m_nRows = pMp.m_nRows;
		m_nCols = pMp.m_nCols;
		if (m_pEntry)
		{
			delete[] m_pEntry;
			m_pEntry = nullptr;
		}
		int length = m_nCols*m_nRows;
		m_pEntry = new T[length];
		for (int i = 0; i < length; i++)
		{
			m_pEntry[i] = (*pMp.pM)(pMp.pRelativePos[i]);
		}
		return *this;
	}
	Matrix& operator=(const vector<Matrix<T>*>& vM)//��ƴ��
	{
		Matrix<T> tM;
		int nMatrixs = vM.size();//����ƴ�ӵľ�����
		if (!nMatrixs)
		{
			return tM;
		}
		int nRows = (vM[0]->size())(1);
		int nTotalCols = 0;
		//���뱣֤���������һ��
		for (int i = 0; i < nMatrixs; i++)
		{
			if (vM[i]->size(1) != nRows)
			{
				cout << "����ά�Ȳ�һ��\n";
				return tM;
			}
			nTotalCols += vM[i]->size(2);//ͳ�����е�����
		}
		T* pAddr = new T[nRows*nTotalCols];
	
		int index = 0;
		for (int k = 0; k < nMatrixs; k++)//ѭ���ڼ����������
		{
			int length = vM[k]->size(1)*vM[k]->size(2);
			for (int i = 1; i <= length; i++)
			{
				pAddr[index++] = (*vM[k])(i);
			}
		}
		if (m_pEntry)
		{
			delete[]m_pEntry;
		}
		m_nRows = nRows;
		m_nCols = nTotalCols;
		m_pEntry = pAddr;
		return *this;
	}
	Matrix& operator=(const initializer_list<Matrix<T>*>& lM)//��ƴ��
	{
		Matrix<T> tM;
		auto ibeg=lM.begin();
		int nMatrixs = lM.size();//����ƴ�ӵľ�����
		if (!nMatrixs)
		{
			return tM;
		}
		int nRows = ((*(ibeg+0))->size())(1);
		int nTotalCols = 0;
		//���뱣֤���������һ��
		for (int i = 0; i < nMatrixs; i++)
		{
			if ((*(ibeg + i))->size(1) != nRows)
			{
				cout << "����ά�Ȳ�һ��\n";
				return tM;
			}
			nTotalCols += (*(ibeg + i))->size(2);//ͳ�����е�����
		}
		
		T*paddr = new T[nRows*nTotalCols];
		
		int index = 0;
		for (int k = 0; k < nMatrixs; k++)//ѭ���ڼ����������
		{
			int length = (*(ibeg + k))->size(1)*(*(ibeg + k))->size(2);
			for (int i = 1; i <= length; i++)
			{
				paddr[index++] = (*(*(ibeg + k)))(i);
			}
		}
		m_nRows = nRows;
		m_nCols = nTotalCols;
		if (m_pEntry)
		{
			delete[]m_pEntry;
		}
		
		m_pEntry = paddr;
		return *this;
	}
	Matrix& operator=(const initializer_list<initializer_list<Matrix<T>*>>& l2M)//���ƴ��
	{
		//����ƴ�ӣ�����ƴ��
		int lRows = l2M.size();//�ֿ���������
		Matrix<T>* mM = new Matrix<T>[lRows];//�ݴ�ˮƽƴ��
		for (int i = 0; i < lRows; i++)
		{
			vector < Matrix<T>*> vM;//������ƴ��
			auto ibeg = *(l2M.begin() + i);
			for (auto it : ibeg)
			{
				vM.push_back(it);
			}
			mM[i] = vM;//�����е����о���ƴ��
		}
		//====================��ʼƴ��===========================
	
		Matrix<T> tM;
		int nMatrixs = lRows;//����ƴ�ӵľ�����
		if (!nMatrixs)
		{
			return tM;
		}
		int nCols = (mM[0].size())(2);
		int nTotalRows = 0;
		//���뱣֤���������һ��
		for (int i = 0; i < nMatrixs; i++)
		{
			if (mM[i].size(2) != nCols)
			{
				cout << "����ά�Ȳ�һ��\n";
				return tM;
			}
			nTotalRows += mM[i].size(1);//ͳ�����е�����
		}
		
		int index = 0;
		T* pAddr = new T[nTotalRows*nCols];

		for (int j = 1; j <= nCols; j++)
		{
			for (int k = 0; k < nMatrixs; k++)//ѭ���ڼ����������
			{
				int length = mM[k].size(1);

				for (int i = 1; i <= length; i++)
				{
					pAddr[index++] = mM[k](i, j);
				}
			}
		}
		if (m_pEntry)
		{
			delete[]m_pEntry;
		}
		m_nRows = nTotalRows;
		m_nCols = nCols;
		m_pEntry = pAddr;
		delete[]mM;
		return *this;
	}

	Matrix operator+(const initializer_list<initializer_list<Matrix<T>*>>& l2M)
	{
		Matrix<T> tM = l2M;
		tM += (*this);//���üӷ��Ľ�����
		return tM;
	}
	Matrix& operator+=(const initializer_list<initializer_list<Matrix<T>*>>& l2M)
	{
		Matrix<T> tM = l2M;
		(*this) += tM;
		return *this;
	}
	Matrix operator+(const initializer_list<Matrix<T>*>& lM)//��ƴ��
	{
		Matrix<T> tM = lM;
		tM += (*this);//���üӷ��Ľ�����
		return tM;
	}
	Matrix& operator+=(const initializer_list<Matrix<T>*>& lM)
	{
		Matrix<T> tM = lM;
		(*this) += tM;
		return *this;
	}
	Matrix operator+(const initializer_list<T>& il)//��ƴ��
	{
		Matrix<T> tM = il;
		tM += (*this);//���üӷ��Ľ�����
		return tM;
	}
	Matrix& operator+=(const initializer_list<T>& il)
	{
		Matrix<T> tM = il;
		(*this) += tM;
		return *this;
	}
	Matrix operator+(const initializer_list<initializer_list<T>>& il)//��ƴ��
	{
		Matrix<T> tM = il;
		tM += (*this);//���üӷ��Ľ�����
		return tM;
	}
	Matrix& operator+=(const initializer_list<initializer_list<T>>& il)
	{
		Matrix<T> tM = il;
		(*this) += tM;
		return *this;
	}
	Matrix operator-(const initializer_list<initializer_list<Matrix<T>*>>& l2M)
	{
		Matrix<T> tM = l2M;
		tM = (*this) - tM;//���üӷ��Ľ�����
		return tM;
	}
	Matrix& operator-=(const initializer_list<initializer_list<Matrix<T>*>>& l2M)
	{
		Matrix<T> tM = l2M;
		(*this) -= tM;
		return *this;
	}
	Matrix operator-(const initializer_list<Matrix<T>*>& lM)//��ƴ��
	{
		Matrix<T> tM = lM;
		tM = (*this) - tM;//���üӷ��Ľ�����
		return tM;
	}
	Matrix& operator-=(const initializer_list<Matrix<T>*>& lM)
	{
		Matrix<T> tM = lM;
		(*this) -= tM;
		return *this;
	}
	Matrix operator-(const initializer_list<T>& il)//��ƴ��
	{
		Matrix<T> tM = il;
		tM = (*this) - tM;//���üӷ��Ľ�����
		return tM;
	}
	Matrix& operator-=(const initializer_list<T>& il)
	{
		Matrix<T> tM = il;
		(*this) -= tM;
		return *this;
	}
	Matrix operator-(const initializer_list<initializer_list<T>>& il)//��ƴ��
	{
		Matrix<T> tM = il;
		tM = (*this) - tM;//���üӷ��Ľ�����
		return tM;
	}
	Matrix& operator-=(const initializer_list<initializer_list<T>>& il)
	{
		Matrix<T> tM = il;
		(*this) -= tM;
		return *this;
	}
	template<typename Y>
	Matrix operator+(Matrix<Y> &M)
	{
		Matrix<T> tM(NULL, m_nRows, m_nCols);
		if (m_nRows != M.size(1) && m_nCols != M.size(2))
		{
			cout << "���ӷ��ľ���ͬ��\n" << endl;
			return tM;
		}
		if ((is_same<complex<double>, T>::value && is_same<double, Y>::value) || is_same<T, Y>::value)//����Ҫƥ�䣬���ݸ�����double
		{
			for (int i = 1; i <= m_nRows; i++)
			{
				for (int j = 1; j <= m_nCols; j++)
				{
					tM(i, j) = (*this)(i, j) + M(i, j);
				}
			}
		}
		else
		{
			cout << "���ӷ��ľ������Ͳ�ƥ��\n";
		}
		return tM;
	}
	Matrix operator+(const T &key)
	{
		Matrix<T> tM(NULL, m_nRows, m_nCols);
		
		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) + key;
			}
		}
		return tM;
	}
	Matrix& operator+=(Matrix& M)
	{
		Matrix<T>& tM = *this;
		if (m_nRows != M.m_nRows && m_nCols != M.m_nCols)
		{
			cout << "���ӷ��ľ���ͬ��\n" << endl;
			return tM;
		}
		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) + M(i, j);
			}
		}
		return tM;
	}
	Matrix& operator+=(const T &key)
	{
		Matrix<T>& tM = *this;

		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) + key;
			}
		}
		return tM;
	}
	template<typename Y>
	Matrix operator-(Matrix<Y> &M)
	{
		Matrix<T> tM(NULL, m_nRows, m_nCols);
		if (m_nRows != M.size(1) && m_nCols != M.size(2))
		{
			cout << "�������ľ���ͬ��\n" << endl;
			return tM;
		}
		if ((is_same<complex<double>, T>::value && is_same<double, Y>::value) || is_same<T, Y>::value)//����Ҫƥ�䣬���ݸ�����double
		{
			for (int i = 1; i <= m_nRows; i++)
			{
				for (int j = 1; j <= m_nCols; j++)
				{
					tM(i, j) = (*this)(i, j) - M(i, j);
				}
			}
		}
		else
		{
			cout << "�������ľ������Ͳ�ƥ��\n";
		}
		return tM;
	}
	Matrix operator-()//���帺����
	{
		Matrix<T> tM(NULL, m_nRows, m_nCols);
		for (int i = 0; i < m_nRows*m_nCols; i++)
		{
			tM.m_pEntry[i] = -m_pEntry[i];
		}
		return tM;
	}
	Matrix operator-(const T &key)
	{
		Matrix<T> tM(NULL, m_nRows, m_nCols);

		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) - key;
			}
		}
		return tM;
	}
	Matrix& operator-=(Matrix &M)
	{
		Matrix<T>& tM = *this;
		if (m_nRows != M.m_nRows && m_nCols != M.m_nCols)
		{
			cout << "�������ľ���ͬ��\n" << endl;
			return tM;
		}
		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) - M(i, j);
			}
		}
		return tM;
	}
	Matrix& operator-=(const T &key)
	{
		Matrix<T>& tM = *this;

		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) - key;
			}
		}
		return tM;
	}
	Matrix operator*(const T& key)
	{
		Matrix<T> tM(NULL, m_nRows, m_nCols);

		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) * key;
			}
		}
		return tM;
	}
	Matrix& operator*=(const T &key)
	{
		Matrix<T>& tM = *this;

		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) * key;
			}
		}
		return tM;
	}
	Matrix operator*(Matrix &M)//���صľ���˷�
	{
		Matrix tM(NULL, m_nRows, M.m_nCols);
		if (m_nCols != M.m_nRows)
		{
			cout << "�˷�������\n";
			return tM;
		}
		int i(1), j, k;
		do{
			j = 1;
			do{
				k = 1;
				do{
					T& tkey = tM(i, j);
					tkey += (*this)(i, k)*M(k, j);
					k++;
					
				} while (k <= m_nCols);
				j++;
			} while (j <= M.m_nCols);
			i++;
		} while (i <= m_nRows);
		return tM;
	}
	Matrix operator/(const T& key)
	{
		Matrix<T> tM(NULL, m_nRows, m_nCols);

		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) / key;
			}
		}
		return tM;
	}
	Matrix& operator/=(const T &key)
	{
		Matrix<T>& tM = *this;

		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j) / key;
			}
		}
		return tM;
	}


	T& operator()(int i)//�����1��ʼ���ʵ�����
	{
		int index = i - 1;
		if (index >= m_nRows*m_nCols)
		{
			cout << "����Խ��\n";
			T t;
			return t;
		}
		return m_pEntry[index];
	}
	T& operator()(int i, int j)//�����1��ʼ���ʵ�����
	{
		if (i < 1 || j < 1)
		{
			cout << "��������i<1,j<1\n";
			T t;
			return t;
		}
		int index = (i - 1) + (j - 1)*m_nRows;
		if (index >= m_nRows*m_nCols || i > m_nRows || j > m_nCols)
		{
			cout << "����Խ��\n";
			T t;
			return t;
		}
		return m_pEntry[index];
	}
	Matrix_p<T> operator()(Matrix<int>& M)
	{
		//����һ����M�ȴ�С�ľ���ֵΪMָʾ�ľ���λ��
		Matrix_p<T> tMp(nullptr, M.m_nRows, M.m_nCols, this);
		int length = M.m_nRows*M.m_nCols;

		for (int i = 0; i < length; i++)
		{
			if (M.m_pEntry[i] < 1 || M.m_pEntry[i] > m_nRows*m_nCols)//��������Խ��
			{
				cout << "����Խ��\n";
				return tMp;
			}
		}
		for (int i = 0; i < length; i++)
		{
			tMp.pRelativePos[i] = M.m_pEntry[i];//M�����Ԫ����ΪA�����λ������
		}
		return tMp;
	}
	Matrix_p<T>& operator()(Matrix_p<T>& pMp)//������̬��Matrix_p�󶨵�����ľ�����
	{
		pMp.pM = this;
		return pMp;
	}

	Matrix_p<T> operator>(Matrix<T>& M)
	{
		//Ҫ��ȽϾ���ͬ��
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>MĬ�ϰ���A��
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "���ڱȽϵľ���ͬ��\n";
			return tMp;
		}
		int length = nRows*nCols;
		int index = 0;
		for (int i = 1; i <= length; i++)
		{
			if ((*this)(i) > M(i))
			{
				tMp.pRelativePos[index++] = i;
			}
		}
		if (index == length)//���˱�ʾ����Ԫ�ض�����
		{
			tMp.m_isatisfy = true;
		}
		tMp.m_nRows = index;
		tMp.m_nCols = 1;
		if (index == 0)
		{
			tMp.m_nCols = 0;
		}
		return tMp;
	}
	Matrix_p<T> operator>=(Matrix<T>& M)
	{
		//Ҫ��ȽϾ���ͬ��
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>MĬ�ϰ���A��
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "���ڱȽϵľ���ͬ��\n";
			return tMp;
		}
		int length = nRows*nCols;
		int index = 0;
		for (int i = 1; i <= length; i++)
		{
			if ((*this)(i) >= M(i))
			{
				tMp.pRelativePos[index++] = i;
			}
		}
		if (index == length)//���˱�ʾ����Ԫ�ض�����
		{
			tMp.m_isatisfy = true;
		}
		tMp.m_nRows = index;
		tMp.m_nCols = 1;
		if (index == 0)
		{
			tMp.m_nCols = 0;
		}
		return tMp;
	}
	Matrix_p<T> operator<(Matrix<T>& M)
	{
		//Ҫ��ȽϾ���ͬ��
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>MĬ�ϰ���A��
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "���ڱȽϵľ���ͬ��\n";
			return tMp;
		}
		int length = nRows*nCols;
		int index = 0;
		for (int i = 1; i <= length; i++)
		{
			if ((*this)(i) < M(i))
			{
				tMp.pRelativePos[index++] = i;
			}
		}
		if (index == length)//���˱�ʾ����Ԫ�ض�����
		{
			tMp.m_isatisfy = true;
		}
		tMp.m_nRows = index;
		tMp.m_nCols = 1;
		if (index == 0)
		{
			tMp.m_nCols = 0;
		}
		return tMp;
	}
	Matrix_p<T> operator<=(Matrix<T>& M)
	{
		//Ҫ��ȽϾ���ͬ��
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>MĬ�ϰ���A��
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "���ڱȽϵľ���ͬ��\n";
			return tMp;
		}
		int length = nRows*nCols;
		int index = 0;
		for (int i = 1; i <= length; i++)
		{
			if ((*this)(i) <= M(i))
			{
				tMp.pRelativePos[index++] = i;
			}
		}
		if (index == length)//���˱�ʾ����Ԫ�ض�����
		{
			tMp.m_isatisfy = true;
		}
		tMp.m_nRows = index;
		tMp.m_nCols = 1;
		if (index == 0)
		{
			tMp.m_nCols = 0;
		}
		return tMp;
	}
	Matrix_p<T> operator==(Matrix<T>& M)
	{
		//Ҫ��ȽϾ���ͬ��
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>MĬ�ϰ���A��
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "���ڱȽϵľ���ͬ��\n";
			return tMp;
		}
		int length = nRows*nCols;
		int index = 0;
		for (int i = 1; i <= length; i++)
		{
			if ((*this)(i) == M(i))
			{
				tMp.pRelativePos[index++] = i;
			}
		}
		if (index == length)//���˱�ʾ����Ԫ�ض�����
		{
			tMp.m_isatisfy = true;
		}
		tMp.m_nRows = index;
		tMp.m_nCols = 1;
		if (index == 0)
		{
			tMp.m_nCols = 0;
		}
		return tMp;
	}
	Matrix_p<T> operator!=(Matrix<T>& M)
	{
		//Ҫ��ȽϾ���ͬ��
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>MĬ�ϰ���A��
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "���ڱȽϵľ���ͬ��\n";
			return tMp;
		}
		int length = nRows*nCols;
		int index = 0;
		for (int i = 1; i <= length; i++)
		{
			if ((*this)(i) != M(i))
			{
				tMp.pRelativePos[index++] = i;
			}
		}
		if (index != 0)//ֻҪ��һ��Ԫ�ز����ھ��Ǿ������
		{
			tMp.m_isatisfy = true;
		}
		tMp.m_nRows = index;
		tMp.m_nCols = 1;
		if (index == 0)
		{
			tMp.m_nCols = 0;
		}
		return tMp;
	}
	
	template <typename Y>//���������������֮���ǿ��ת��
	operator Matrix<Y>()
	{
		Matrix<Y> tM(NULL, m_nRows, m_nCols);
		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = 1; j <= m_nCols; j++)
			{
				tM(i, j) = (*this)(i, j);
			}
		}
		return tM;
	}
	Matrix_p<T> row(int i)//ȡ����������ʼ��1
	{
		//����������ʼ��1
		Matrix_p<T> tMp(nullptr, 1, m_nCols, this);
		if (i > m_nRows || i < 1)
		{
			cout << "����������Χ\n";
			return tMp;
		}
		for (int k = 0; k < m_nCols; k++)
		{
			tMp.pRelativePos[k] = i + k*m_nRows;//��¼��Ӧλ��
		}
		
		return tMp;
	}
	Matrix_p<T> row(int i, int begin, int end)//ȡ����������ʼ��1
	{
		int length = end - begin + 1;
		Matrix_p<T> tMp(nullptr, 1, length, this);
		if (i > m_nRows || i < 1 || length < 1 || begin > m_nCols || begin < 1 || end < 1 || end > m_nCols)
		{
			cout << "����������Χ\n";
			return tMp;
		}

		for (int k = begin - 1; k < end; k++)
		{
			tMp.pRelativePos[k - begin + 1] = i + k*m_nRows;
		}
		return tMp;
	}
	Matrix_p<T> col(int j)//ȡ����������ʼ��1
	{
		Matrix_p<T> tMp(nullptr, m_nRows, 1, this);
		if (j > m_nCols || j < 1)
		{
			cout << "����������Χ\n";
			return tMp;
		}
		for (int k = 0; k < m_nRows; k++)
		{
			tMp.pRelativePos[k] = k + (j - 1)*m_nRows + 1;
		}
		return tMp;
	}
	Matrix_p<T> col(int j, int begin, int end)//ȡ����������ʼ��1
	{
		int length = end - begin + 1;
		Matrix_p<T> tMp(nullptr, length, 1, this);
		if (j > m_nCols || j < 1 || length < 1 || begin > m_nRows || begin < 1 || end < 1 || end > m_nRows)
		{
			cout << "����������Χ\n";
			return tMp;
		}
		for (int k = begin - 1; k < end; k++)
		{
			tMp.pRelativePos[k - begin + 1] = k + (j - 1)*m_nRows + 1;
		}
		return tMp;
	}
	Matrix& horzcat(vector<Matrix<T>*> vM);//ˮƽƴ�ӣ������Լ�
	Matrix& vertcat(vector<Matrix<T>*> vM);//��ֱƴ�ӣ������Լ�
	void clear();
	Matrix<int> size();
	int size(int para);
	Matrix_p<T> GetSubMatrix(int r_beg, int c_beg, int r_end, int c_end)
	{
		int nRows = r_end - r_beg + 1;
		int nCols = c_end - c_beg + 1;
		if (nRows < 1 || nCols < 1 || r_end>m_nRows || c_end>m_nCols || r_beg < 1 || c_beg < 1)
		{
			cout << "��ʼ��ֹλ�ô���\n";
			Matrix_p<T> tMp(nullptr, 0, 0, this);
			return tMp;
		}
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);
		int index = 0;
		for (int j = c_beg; j <= c_end; j++)
		{
			for (int i = r_beg; i <= r_end; i++)
			{
				tMp.pRelativePos[index++] = i + (j - 1)*m_nRows;
			}
		}
		return tMp;
	}
	Matrix GetCofactorMatrix(int r_i, int c_j);//��Ӧλ�õ�����ʽ
	void del(Matrix<int>& M);//ɾ��M��ָʾ��λ�ô���Ԫ�أ�������������
	void del(int r_i, int c_j);//ɾ����i��j�У���������ʽ����i,j=0ʱ��ʾֻɾ����Ӧ��,��
	void del(int i);//ɾ����i��λ�õ�Ԫ��
	Matrix trans();//ת��
	void trans(Matrix& mat);//ת�ã����浽����ľ�����
	Matrix repmat(int r_i, int c_j);//��A��ΪԪ�أ���չ��i*j�ķֿ����
	Matrix repmat(int r_c_i);
	Matrix elemul(Matrix<T>& M);//��ӦԪ��֮��
	Matrix elediv(Matrix<T>& M);//��ӦԪ��֮��
	Matrix& swap_row(int r_i, int r_j);//������
	Matrix& swap_col(int c_i, int c_j);//������
	Matrix& eye(int n);
	void tozero();//����Ԫ������
	void Symmetrized(bool u2d = true); 
	void resize(int nRows, int nCols);
	bool save(wchar_t* filename, wchar_t* filenameImag = nullptr);//�ڶ�������ȱʡ�������鲿
	bool savebin(const wchar_t* filename, const wchar_t* filenameImag = nullptr);//����Ϊ�������ļ����ڶ�������ȱʡ�������鲿
	bool savebin(const char* filename, const char* filenameImag = nullptr);//����Ϊ�������ļ����ڶ�������ȱʡ�������鲿
	bool load(wchar_t* filename);
	bool loadbin(wchar_t* filename);
	int split_real_imag(Matrix<double>& matReal, Matrix<double>& matimag, char C = 'A');//���ʵ�����鲿Ϊ��������,Aʵ���鲿��Ҫ��RֻҪʵ����IֻҪ�鲿
	int split_real_imag(double* dReal, double* dimag, char C = 'A');
	const T* const GetEntry();//�ǳ�����ȫ�ķ�ʽ������˼������
	void SetEntry(T* entry);//�������ݿռ��ڴ湫�ã���Ҫ����ȷ������Ŀռ�ͱ�������ڴ�һ�����ұ����ö��е��ڴ�
	void SetEntry(int nrows, int ncols, T* entry);
	void Dense2Sparse(Mat_Sparse<T>& mat_spa);//����ѹ�����棬ǰ���Ǿ�����ϡ�����
	void Dense2Sparse(Matrix<T> mat_B, Mat_Sparse<T>& mat_spa, Mat_Sparse<T>& mat_spa_comp);
	void Sparse2Dense(Mat_Sparse<T>& mat_spa);//��ȡһ��ϡ�财��ľ��󵽾����ڴ���
	template<typename Y>
	void Sparse2Dense_potential(Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Y alpha);//��Ȩ��ȡϡ������ڴ�֮��
	template<typename Y>
	void Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Y alpha);//��Ȩ��ȡϡ����󵽶�Ӧ�ֿ�֮��
	template<typename Y>
	void Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, T alpha, T beta);//��Ȩ��ȡϡ����󵽶�Ӧ�ֿ�֮��
	template<typename Y>
	void Sparse2Dense_potential_add(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, T alpha, T beta);//��Ȩ��ȡϡ����󵽶�Ӧ�ֿ�֮��
	template<typename Y>//ksg����ר��
	void Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Mat_Sparse<Y>& mat_spa_ksg, T alpha, T beta);
	void toinv();//ԭַ����
	~Matrix();//����һ��������
};


