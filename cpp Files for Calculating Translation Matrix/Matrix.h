#pragma once
#include <iomanip>
#include <vector>
#include <set>
#include <initializer_list>
#include <complex>
#include <algorithm>
#include <fstream>
using namespace std;

template <typename T>//只对实数矩阵稀疏化，复数矩阵稀疏化会有歧义
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
	void Bindindex(Mat_Sparse<T>& mat_spa_comp)//稀疏跟随矩阵（索引位置完全一致，数据不同）
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
			cout << "伴生稀疏矩阵不允许更改尺寸\n";
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
		Matrix<T>* pM;//源头矩阵，等于源头this
		int* pRelativePos;//记录数据在源头矩阵的位置

		Matrix_p(char* c, int nRows, int nCols, Matrix<T>* ptM)
		{
			m_nRows = nRows;
			m_nCols = nCols;
			pM = ptM;//绑定矩阵
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
				cout << "Matrix_p索引越界\n";
				T tkey;
				return tkey;
			}
			return pRelativePos[i - 1];
		}
		void operator=(Matrix<T>& M)
		{
			//要求总长度一致
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
		void operator=(const T ele)//传入一个元素的时候，表示对该矩阵指示的位置的元素都等于这个数
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
			if ((is_same<complex<double>, T>::value && is_same<double, Y>::value) || is_same<T, Y>::value)//类型要匹配，兼容复数和double
			{
				int tRows = M.size(1);
				int tCols = M.size(2);
				int length = m_nRows*m_nCols;
				if (tRows*tCols != length)
				{
					cout << "矩阵维度不一致\n";
					return tM;
				}
				for (int i = 1; i <= length; i++)
				{
					tM(i) += M(i);
				}
			}
			else
			{
				cout << "矩阵类型不匹配\n";
			}
			//要求总长度一致
			return tM;
		}		
		Matrix operator+(Matrix_p<T>& pMp)
		{
			//要求总长度一致
			int tRows = pMp.m_nRows;
			int tCols = pMp.m_nCols;
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
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
			//要求总长度一致
			int tRows = pMp.m_nRows;
			int tCols = pMp.m_nCols;
			int length = m_nRows*m_nCols;
			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
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
			//要求总长度一致
			int tRows = pMp.m_nRows;
			int tCols = pMp.m_nCols;
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
			int tRows = pMp.m_nRows;
			int tCols = pMp.m_nCols;
			int length = m_nRows*m_nCols;
			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
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
			//要求总长度一致
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
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
			//要求总长度一致
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;
			Matrix<T> tM = *this;
			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
			//要求总长度一致
			int tRows = M.size(1);
			int tCols = M.size(2);
			int length = m_nRows*m_nCols;

			if (tRows*tCols != length)
			{
				cout << "矩阵维度不一致\n";
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
		explicit operator bool() const//如果矩阵中无元素，那么返回false，否则返回true explicit表示只能强制转换，避免二义性
		{
			return m_isatisfy;
		}
		void zero()//将元素置零
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
	bool m_blflag;//用于标记该矩阵是否是一个绑定矩阵，存在潜在的内存释放的问题
public:
	Matrix(char c = NULL, int nRows = 0, int nCols = 0);//参数1是形式化的参数
	Matrix(const initializer_list<initializer_list<T>>& il);//行列构造
	Matrix(const initializer_list<T>& il);//行构造
	Matrix(const T& ele);
	Matrix(const Matrix& M);//拷贝构造

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
				out << endl;//多输出一个
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
		//纠正维度不一致的错误
		for (int i = 0; i < m_nRows; i++)
		{
			if ((ibeg + i)->size() != m_nCols)
			{
				clear();
				cout << "数据维度不一致\n";
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
	Matrix& operator=(const Matrix_p<T>& pMp)//将指针指示的位置的元素值传给对应矩阵
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
	Matrix& operator=(const vector<Matrix<T>*>& vM)//行拼接
	{
		Matrix<T> tM;
		int nMatrixs = vM.size();//用于拼接的矩阵数
		if (!nMatrixs)
		{
			return tM;
		}
		int nRows = (vM[0]->size())(1);
		int nTotalCols = 0;
		//必须保证矩阵的行数一致
		for (int i = 0; i < nMatrixs; i++)
		{
			if (vM[i]->size(1) != nRows)
			{
				cout << "矩阵维度不一致\n";
				return tM;
			}
			nTotalCols += vM[i]->size(2);//统计所有的列数
		}
		T* pAddr = new T[nRows*nTotalCols];
	
		int index = 0;
		for (int k = 0; k < nMatrixs; k++)//循环第几个矩阵参数
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
	Matrix& operator=(const initializer_list<Matrix<T>*>& lM)//行拼接
	{
		Matrix<T> tM;
		auto ibeg=lM.begin();
		int nMatrixs = lM.size();//用于拼接的矩阵数
		if (!nMatrixs)
		{
			return tM;
		}
		int nRows = ((*(ibeg+0))->size())(1);
		int nTotalCols = 0;
		//必须保证矩阵的行数一致
		for (int i = 0; i < nMatrixs; i++)
		{
			if ((*(ibeg + i))->size(1) != nRows)
			{
				cout << "矩阵维度不一致\n";
				return tM;
			}
			nTotalCols += (*(ibeg + i))->size(2);//统计所有的列数
		}
		
		T*paddr = new T[nRows*nTotalCols];
		
		int index = 0;
		for (int k = 0; k < nMatrixs; k++)//循环第几个矩阵参数
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
	Matrix& operator=(const initializer_list<initializer_list<Matrix<T>*>>& l2M)//混合拼接
	{
		//先行拼接，再列拼接
		int lRows = l2M.size();//分块矩阵的行数
		Matrix<T>* mM = new Matrix<T>[lRows];//暂存水平拼接
		for (int i = 0; i < lRows; i++)
		{
			vector < Matrix<T>*> vM;//用于列拼接
			auto ibeg = *(l2M.begin() + i);
			for (auto it : ibeg)
			{
				vM.push_back(it);
			}
			mM[i] = vM;//将其中的所有矩阵拼接
		}
		//====================开始拼接===========================
	
		Matrix<T> tM;
		int nMatrixs = lRows;//用于拼接的矩阵数
		if (!nMatrixs)
		{
			return tM;
		}
		int nCols = (mM[0].size())(2);
		int nTotalRows = 0;
		//必须保证矩阵的列数一致
		for (int i = 0; i < nMatrixs; i++)
		{
			if (mM[i].size(2) != nCols)
			{
				cout << "矩阵维度不一致\n";
				return tM;
			}
			nTotalRows += mM[i].size(1);//统计所有的行数
		}
		
		int index = 0;
		T* pAddr = new T[nTotalRows*nCols];

		for (int j = 1; j <= nCols; j++)
		{
			for (int k = 0; k < nMatrixs; k++)//循环第几个矩阵参数
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
		tM += (*this);//利用加法的交换律
		return tM;
	}
	Matrix& operator+=(const initializer_list<initializer_list<Matrix<T>*>>& l2M)
	{
		Matrix<T> tM = l2M;
		(*this) += tM;
		return *this;
	}
	Matrix operator+(const initializer_list<Matrix<T>*>& lM)//行拼接
	{
		Matrix<T> tM = lM;
		tM += (*this);//利用加法的交换律
		return tM;
	}
	Matrix& operator+=(const initializer_list<Matrix<T>*>& lM)
	{
		Matrix<T> tM = lM;
		(*this) += tM;
		return *this;
	}
	Matrix operator+(const initializer_list<T>& il)//行拼接
	{
		Matrix<T> tM = il;
		tM += (*this);//利用加法的交换律
		return tM;
	}
	Matrix& operator+=(const initializer_list<T>& il)
	{
		Matrix<T> tM = il;
		(*this) += tM;
		return *this;
	}
	Matrix operator+(const initializer_list<initializer_list<T>>& il)//行拼接
	{
		Matrix<T> tM = il;
		tM += (*this);//利用加法的交换律
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
		tM = (*this) - tM;//利用加法的交换律
		return tM;
	}
	Matrix& operator-=(const initializer_list<initializer_list<Matrix<T>*>>& l2M)
	{
		Matrix<T> tM = l2M;
		(*this) -= tM;
		return *this;
	}
	Matrix operator-(const initializer_list<Matrix<T>*>& lM)//行拼接
	{
		Matrix<T> tM = lM;
		tM = (*this) - tM;//利用加法的交换律
		return tM;
	}
	Matrix& operator-=(const initializer_list<Matrix<T>*>& lM)
	{
		Matrix<T> tM = lM;
		(*this) -= tM;
		return *this;
	}
	Matrix operator-(const initializer_list<T>& il)//行拼接
	{
		Matrix<T> tM = il;
		tM = (*this) - tM;//利用加法的交换律
		return tM;
	}
	Matrix& operator-=(const initializer_list<T>& il)
	{
		Matrix<T> tM = il;
		(*this) -= tM;
		return *this;
	}
	Matrix operator-(const initializer_list<initializer_list<T>>& il)//行拼接
	{
		Matrix<T> tM = il;
		tM = (*this) - tM;//利用加法的交换律
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
			cout << "作加法的矩阵不同型\n" << endl;
			return tM;
		}
		if ((is_same<complex<double>, T>::value && is_same<double, Y>::value) || is_same<T, Y>::value)//类型要匹配，兼容复数和double
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
			cout << "作加法的矩阵类型不匹配\n";
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
			cout << "作加法的矩阵不同型\n" << endl;
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
			cout << "作减法的矩阵不同型\n" << endl;
			return tM;
		}
		if ((is_same<complex<double>, T>::value && is_same<double, Y>::value) || is_same<T, Y>::value)//类型要匹配，兼容复数和double
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
			cout << "作减法的矩阵类型不匹配\n";
		}
		return tM;
	}
	Matrix operator-()//定义负矩阵
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
			cout << "作减法的矩阵不同型\n" << endl;
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
	Matrix operator*(Matrix &M)//朴素的矩阵乘法
	{
		Matrix tM(NULL, m_nRows, M.m_nCols);
		if (m_nCols != M.m_nRows)
		{
			cout << "乘法无意义\n";
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


	T& operator()(int i)//定义从1开始访问的索引
	{
		int index = i - 1;
		if (index >= m_nRows*m_nCols)
		{
			cout << "访问越界\n";
			T t;
			return t;
		}
		return m_pEntry[index];
	}
	T& operator()(int i, int j)//定义从1开始访问的索引
	{
		if (i < 1 || j < 1)
		{
			cout << "索引错误i<1,j<1\n";
			T t;
			return t;
		}
		int index = (i - 1) + (j - 1)*m_nRows;
		if (index >= m_nRows*m_nCols || i > m_nRows || j > m_nCols)
		{
			cout << "访问越界\n";
			T t;
			return t;
		}
		return m_pEntry[index];
	}
	Matrix_p<T> operator()(Matrix<int>& M)
	{
		//产生一个和M等大小的矩阵，值为M指示的矩阵位置
		Matrix_p<T> tMp(nullptr, M.m_nRows, M.m_nCols, this);
		int length = M.m_nRows*M.m_nCols;

		for (int i = 0; i < length; i++)
		{
			if (M.m_pEntry[i] < 1 || M.m_pEntry[i] > m_nRows*m_nCols)//索引不能越界
			{
				cout << "索引越界\n";
				return tMp;
			}
		}
		for (int i = 0; i < length; i++)
		{
			tMp.pRelativePos[i] = M.m_pEntry[i];//M矩阵的元素作为A矩阵的位置索引
		}
		return tMp;
	}
	Matrix_p<T>& operator()(Matrix_p<T>& pMp)//将游离态的Matrix_p绑定到具体的矩阵上
	{
		pMp.pM = this;
		return pMp;
	}

	Matrix_p<T> operator>(Matrix<T>& M)
	{
		//要求比较矩阵同型
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>M默认绑定于A上
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "用于比较的矩阵不同型\n";
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
		if (index == length)//满了表示所有元素都满足
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
		//要求比较矩阵同型
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>M默认绑定于A上
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "用于比较的矩阵不同型\n";
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
		if (index == length)//满了表示所有元素都满足
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
		//要求比较矩阵同型
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>M默认绑定于A上
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "用于比较的矩阵不同型\n";
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
		if (index == length)//满了表示所有元素都满足
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
		//要求比较矩阵同型
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>M默认绑定于A上
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "用于比较的矩阵不同型\n";
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
		if (index == length)//满了表示所有元素都满足
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
		//要求比较矩阵同型
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>M默认绑定于A上
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "用于比较的矩阵不同型\n";
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
		if (index == length)//满了表示所有元素都满足
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
		//要求比较矩阵同型
		int nRows = M.m_nRows;
		int nCols = M.m_nCols;
		Matrix_p<T> tMp(nullptr, nRows, nCols, this);//A>M默认绑定于A上
		if (nRows != m_nRows || nCols != m_nCols)
		{
			cout << "用于比较的矩阵不同型\n";
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
		if (index != 0)//只要有一个元素不等于就是矩阵不相等
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
	
	template <typename Y>//定义任意矩阵类型之间的强制转换
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
	Matrix_p<T> row(int i)//取行向量，起始于1
	{
		//控制坐标起始于1
		Matrix_p<T> tMp(nullptr, 1, m_nCols, this);
		if (i > m_nRows || i < 1)
		{
			cout << "索引超出范围\n";
			return tMp;
		}
		for (int k = 0; k < m_nCols; k++)
		{
			tMp.pRelativePos[k] = i + k*m_nRows;//记录对应位置
		}
		
		return tMp;
	}
	Matrix_p<T> row(int i, int begin, int end)//取行向量，起始于1
	{
		int length = end - begin + 1;
		Matrix_p<T> tMp(nullptr, 1, length, this);
		if (i > m_nRows || i < 1 || length < 1 || begin > m_nCols || begin < 1 || end < 1 || end > m_nCols)
		{
			cout << "索引超出范围\n";
			return tMp;
		}

		for (int k = begin - 1; k < end; k++)
		{
			tMp.pRelativePos[k - begin + 1] = i + k*m_nRows;
		}
		return tMp;
	}
	Matrix_p<T> col(int j)//取列向量，起始于1
	{
		Matrix_p<T> tMp(nullptr, m_nRows, 1, this);
		if (j > m_nCols || j < 1)
		{
			cout << "索引超出范围\n";
			return tMp;
		}
		for (int k = 0; k < m_nRows; k++)
		{
			tMp.pRelativePos[k] = k + (j - 1)*m_nRows + 1;
		}
		return tMp;
	}
	Matrix_p<T> col(int j, int begin, int end)//取行向量，起始于1
	{
		int length = end - begin + 1;
		Matrix_p<T> tMp(nullptr, length, 1, this);
		if (j > m_nCols || j < 1 || length < 1 || begin > m_nRows || begin < 1 || end < 1 || end > m_nRows)
		{
			cout << "索引超出范围\n";
			return tMp;
		}
		for (int k = begin - 1; k < end; k++)
		{
			tMp.pRelativePos[k - begin + 1] = k + (j - 1)*m_nRows + 1;
		}
		return tMp;
	}
	Matrix& horzcat(vector<Matrix<T>*> vM);//水平拼接，包括自己
	Matrix& vertcat(vector<Matrix<T>*> vM);//垂直拼接，包括自己
	void clear();
	Matrix<int> size();
	int size(int para);
	Matrix_p<T> GetSubMatrix(int r_beg, int c_beg, int r_end, int c_end)
	{
		int nRows = r_end - r_beg + 1;
		int nCols = c_end - c_beg + 1;
		if (nRows < 1 || nCols < 1 || r_end>m_nRows || c_end>m_nCols || r_beg < 1 || c_beg < 1)
		{
			cout << "起始终止位置错误\n";
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
	Matrix GetCofactorMatrix(int r_i, int c_j);//对应位置的余子式
	void del(Matrix<int>& M);//删除M所指示的位置处的元素，并产生行向量
	void del(int r_i, int c_j);//删除第i行j列，产生余子式，当i,j=0时表示只删除对应列,行
	void del(int i);//删除第i个位置的元素
	Matrix trans();//转置
	void trans(Matrix& mat);//转置，保存到传入的矩阵上
	Matrix repmat(int r_i, int c_j);//将A视为元素，扩展成i*j的分块矩阵
	Matrix repmat(int r_c_i);
	Matrix elemul(Matrix<T>& M);//对应元素之积
	Matrix elediv(Matrix<T>& M);//对应元素之商
	Matrix& swap_row(int r_i, int r_j);//交换行
	Matrix& swap_col(int c_i, int c_j);//交换列
	Matrix& eye(int n);
	void tozero();//所有元素置零
	void Symmetrized(bool u2d = true); 
	void resize(int nRows, int nCols);
	bool save(wchar_t* filename, wchar_t* filenameImag = nullptr);//第二个参数缺省，保存虚部
	bool savebin(const wchar_t* filename, const wchar_t* filenameImag = nullptr);//保存为二进制文件，第二个参数缺省，保存虚部
	bool savebin(const char* filename, const char* filenameImag = nullptr);//保存为二进制文件，第二个参数缺省，保存虚部
	bool load(wchar_t* filename);
	bool loadbin(wchar_t* filename);
	int split_real_imag(Matrix<double>& matReal, Matrix<double>& matimag, char C = 'A');//拆分实部和虚部为两个矩阵,A实部虚部都要，R只要实部，I只要虚部
	int split_real_imag(double* dReal, double* dimag, char C = 'A');
	const T* const GetEntry();//非常不安全的方式，后期思考调整
	void SetEntry(T* entry);//用于数据空间内存公用，需要自行确定传入的空间和本矩阵的内存一样大，且必须用堆中的内存
	void SetEntry(int nrows, int ncols, T* entry);
	void Dense2Sparse(Mat_Sparse<T>& mat_spa);//矩阵压缩储存，前提是矩阵是稀疏矩阵
	void Dense2Sparse(Matrix<T> mat_B, Mat_Sparse<T>& mat_spa, Mat_Sparse<T>& mat_spa_comp);
	void Sparse2Dense(Mat_Sparse<T>& mat_spa);//读取一个稀疏储存的矩阵到矩阵内存中
	template<typename Y>
	void Sparse2Dense_potential(Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Y alpha);//带权读取稀疏矩阵到内存之中
	template<typename Y>
	void Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Y alpha);//带权读取稀疏矩阵到对应分块之中
	template<typename Y>
	void Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, T alpha, T beta);//带权读取稀疏矩阵到对应分块之中
	template<typename Y>
	void Sparse2Dense_potential_add(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, T alpha, T beta);//带权读取稀疏矩阵到对应分块之中
	template<typename Y>//ksg算子专用
	void Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Mat_Sparse<Y>& mat_spa_ksg, T alpha, T beta);
	void toinv();//原址求逆
	~Matrix();//返回一个行向量
};


