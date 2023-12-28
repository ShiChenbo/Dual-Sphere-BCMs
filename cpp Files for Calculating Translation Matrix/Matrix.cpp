#include "Matrix.h"



template < typename T >
Matrix<T>::Matrix(char c ,int nRows /*= 0*/, int nCols /*= 0*/)
{
	m_nRows = nRows;
	m_nCols = nCols;
	m_blflag = false;
	if (m_nRows && m_nCols)
	{
		int length = m_nRows * m_nCols;
		m_pEntry = new T[length];
		memset(m_pEntry, 0, sizeof(T)*length);
	}
}

template <typename T>
Matrix<T>::Matrix(const initializer_list<initializer_list<T>>& il)
{
	m_blflag = false;
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
}

template <typename T>
Matrix<T>::Matrix(const initializer_list<T>& il)
{
	m_blflag = false;
	m_nCols = il.size();
	m_nRows = 1;
	
	m_pEntry = new T[m_nCols];

	for (int j = 0; j < m_nCols; j++)
	{
		m_pEntry[j] = *(il.begin()+j);
	}
}
template <typename T>
void Matrix<T>::clear()
{
	m_nRows = 0;
	m_nCols = 0;
	if (m_pEntry && m_blflag == false)//否则说明是伴生矩阵，内存手动释放
	{
		delete[] m_pEntry;
	}
	m_pEntry = nullptr;
}


template < typename T>
void Matrix<T>::tozero()
{
	int size = m_nRows * m_nCols;
	memset(m_pEntry, 0, sizeof(T) * size);
}

template <typename T>
Matrix<T>::~Matrix()
{
	clear();
}

template <typename T>
Matrix<T>::Matrix(const Matrix& M)
{
	*this = M;
}

template <typename T>
Matrix<T>::Matrix(const T& ele)
{
	*this = ele;
}

template < typename T >
Matrix<int> Matrix<T>::size()
{
	Matrix<int> tM({ m_nRows, m_nCols });
	return tM;
}


template < typename T >
Matrix<T>& Matrix<T>::vertcat(vector<Matrix<T>*> vM)
{
	vM.insert(begin(vM), this);
	int nMatrixs = vM.size();//用于拼接的矩阵数
	int nCols = m_nCols;
	int nTotalRows = 0;
	//必须保证矩阵的列数一致
	for (int i = 0; i < nMatrixs; i++)
	{
		if (vM[i]->size(2) != nCols)
		{
			Matrix<T> tM;
			cout << "矩阵维度不一致\n";
			return tM;
		}
		nTotalRows += vM[i]->size(1);//统计所有的行数
	}
	
	T* tpaddr = new  T[nTotalRows*nCols];

	int index = 0;
	for (int j = 1; j <= nCols; j++)
	{
		for (int k = 0; k < nMatrixs; k++)//循环第几个矩阵参数
		{
			int length = vM[k]->size(1);
			for (int i = 1; i <= length; i++)
			{
				tpaddr[index++] = (*vM[k])(i, j);
			}
		}
	}
	if (m_pEntry)
	{
		delete[]m_pEntry;
	}
	m_nRows = nTotalRows;
	m_nCols = nCols;
	m_pEntry = tpaddr;
	return *this;
}

template < typename T >
Matrix<T>& Matrix<T>::horzcat(vector<Matrix<T>*> vM)
{
	Matrix<T> tM;
	vM.insert(begin(vM), this);
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

	T* tpaddr = new  T[nRows*nTotalCols];
	
	int index = 0;
	for (int k = 0; k < nMatrixs; k++)//循环第几个矩阵参数
	{
		int length = vM[k]->size(1)*vM[k]->size(2);
		for (int i = 1; i <= length; i++)
		{
			tpaddr[index++] = (*vM[k])(i);
		}
	}
	if (m_pEntry)
	{
		delete[]m_pEntry;
	}
	m_pEntry = tpaddr;
	m_nRows = nRows;
	m_nCols = nTotalCols;
	return *this;
}

template < typename T >
int Matrix<T>::size(int para)
{
	if (para != 1 && para != 2)
	{
		cout << "size参数传递错误\n";
		return -1;
	}
	if (para == 1)
	{
		return m_nRows;
	}
	return m_nCols;
}


template < typename T >
Matrix<T> Matrix<T>::repmat(int r_i, int c_j)
{
	if (r_i < 1 || c_j < 1)
	{
		cout << "repmat参数传递错误\n";
		Matrix<T> tM;
		return tM;
	}
	int nRows = r_i*m_nRows;
	int nCols = c_j*m_nCols;
	Matrix<T> tM(NULL, nRows, nCols);
	int index = 1;
	for (int j = 0; j < c_j; j++)//列重复次数
	{
		for (int jj = 1; jj <= m_nCols; jj++)//遍历列
		{
			for (int i = 0; i < r_i; i++)//行重复次数
			{
				for (int ii = 1; ii <= m_nRows; ii++)//遍历行
				{
					tM(index++) = (*this)(ii, jj);
				}
			}
		}
	}
	return tM;
	
}

template < typename T >
Matrix<T> Matrix<T>::repmat(int r_c_i)
{
	int r_i = r_c_i, c_j = r_c_i;
	if (r_i < 1 || c_j < 1)
	{
		cout << "repmat参数传递错误\n";
		Matrix<T> tM;
		return tM;
	}
	int nRows = r_i*m_nRows;
	int nCols = c_j*m_nCols;
	Matrix<T> tM(NULL, nRows, nCols);
	int index = 1;
	for (int j = 0; j < c_j; j++)//列重复次数
	{
		for (int jj = 1; jj <= m_nCols; jj++)//遍历列
		{
			for (int i = 0; i < r_i; i++)//行重复次数
			{
				for (int ii = 1; ii <= m_nRows; ii++)//遍历行
				{
					tM(index++) = (*this)(ii, jj);
				}
			}
		}
	}
	return tM;
}


template < typename T >
Matrix<T> Matrix<T>::trans()
{
	int nRows = m_nCols;
	int nCols = m_nRows;
	Matrix<T> tM(NULL, nRows, nCols);
	for (int i = 1; i <= m_nRows; i++)
	{
		for (int j = 1; j <= m_nCols; j++)
		{
			tM(j, i) = (*this)(i, j);
		}
	}
	return tM;
}
template < typename T>
void Matrix<T>::trans(Matrix<T>& mat)
{
	if (mat.m_nRows!=m_nCols || mat.m_nCols != m_nRows)
	{
		cout << "存放转置的矩阵大小不一致\n";
	}
	for (int i = 1; i <= m_nRows; i++)
	{
		for (int j = 1; j <= m_nCols; j++)
		{
			if ((*this)(i, j) == 0.0)
			{
				(*this)(i, j) = mat(j, i);
			}
			else 
			{
				mat(j, i) = (*this)(i, j);
			}
		}
	}
}
template < typename T >
void Matrix<T>::del(int i)
{
	int length = m_nRows*m_nCols - 1;//删除之后的矩阵空间
	T* pEle = new T[length];
	int index = 0;
	for (int k = 0; k < i-1; k++)
	{
		pEle[index++] = m_pEntry[k];
	}
	for (int k = i; k < m_nRows*m_nCols; k++)
	{
		pEle[index++] = m_pEntry[k];
	}

	m_nRows = 1;
	m_nCols = length;
	if (m_pEntry)
	{
		delete[] m_pEntry;
	}
	m_pEntry = pEle;
}

template < typename T >
void Matrix<T>::del(int r_i, int c_j)
{
	if (r_i == 0 && c_j == 0)
	{
		cout << "删除位置指示错误\n";
		return;
	}
	if (r_i<0 || c_j<0 || r_i>m_nRows || c_j>m_nCols)
	{
		cout << "索引越界\n";
		return;
	}
	T* pEle;
	int index = 0;
	if (r_i == 0)//删除对应列，其余列左移
	{
		pEle = new T[m_nRows*(m_nCols - 1)];
		for (int j = 0; j < m_nCols; j++)
		{
			for (int i = 0; i < m_nRows; i++)
			{
				if (j + 1 == c_j)
				{
					break;
				}
				pEle[index++] = (*this)(i + 1, j + 1);
			}
		}
		m_nCols--;
	}
	else if (c_j == 0)//删除对应行，其余行上移
	{
		pEle = new T[(m_nRows - 1)*m_nCols];
		for (int j = 0; j < m_nCols; j++)
		{
			for (int i = 0; i < m_nRows; i++)
			{
				if (i + 1 == r_i)
				{
					continue;
				}
				pEle[index++] = (*this)(i + 1, j + 1);
			}
		}
		m_nRows--;
	}
	else//产生余子式
	{
		pEle = new T[(m_nRows - 1)*(m_nCols - 1)];
		for (int j = 0; j < m_nCols; j++)
		{
			for (int i = 0; i < m_nRows; i++)
			{
				if (j + 1 == c_j)
				{
					break;
				}
				if (i + 1 == r_i)
				{
					continue;
				}
				pEle[index++] = (*this)(i + 1, j + 1);
			}
		}
		m_nRows--;
		m_nCols--;
	}
	if (m_pEntry)
	{
		delete[]m_pEntry;
	}
	m_pEntry = pEle;
}

template < typename T >
void Matrix<T>::del(Matrix<int>& M)
{
	set<int> ts;
	for (int i = 1; i <= M.size(1)*M.size(2); i++)
	{
		ts.insert(M(i));
	}
	int length = m_nRows*m_nCols - ts.size();//删除之后的矩阵空间
	T* pEle = new T[length];
	int index = 0;
	for (int i = 0; i < m_nRows*m_nCols; i++)
	{
		if (ts.find(i + 1) == ts.end())
		{
			pEle[index++] = m_pEntry[i];
		}
	}
	m_nRows = 1;
	m_nCols = length;
	if (m_pEntry)
	{
		delete[] m_pEntry;
	}
	m_pEntry = pEle;
}


template < typename T >
Matrix<T> Matrix<T>::GetCofactorMatrix(int r_i, int c_j)
{
	Matrix<T> tM(NULL, m_nRows - 1, m_nCols - 1);
	if (r_i<1 || c_j<1 || r_i>m_nRows || c_j>m_nCols)
	{
		cout << "传递索引越界\n";
		return tM;
	}
	int index = 1;
	for (int j = 0; j < m_nCols; j++)
	{
		for (int i = 0; i < m_nRows; i++)
		{
			if (j + 1 == c_j)
			{
				break;
			}
			if (i + 1 == r_i)
			{
				continue;
			}
			tM(index++) = (*this)(i + 1, j + 1);
		}
	}
	return tM;
}

template < typename T >
Matrix<T>::Matrix(const initializer_list<initializer_list<Matrix<T>*>>& l2M)
{
	*this = l2M;
}

template < typename T >
Matrix<T>::Matrix(const initializer_list<Matrix<T>*>& lM)
{
	*this = lM;
}

template < typename T >
Matrix<T>::Matrix(const Matrix_p<T>& pMp)
{
	*this = pMp;
}

template < typename T >
Matrix<T> Matrix<T>::elemul(Matrix<T>& M)
{
	Matrix<T> tM(NULL, m_nRows, m_nCols);
	if (m_nRows != M.m_nRows && m_nCols != M.m_nCols)
	{
		cout << "作对应元素乘法的矩阵不同型\n" << endl;
		return tM;
	}
	for (int i = 1; i <= m_nRows; i++)
	{
		for (int j = 1; j <= m_nCols; j++)
		{
			tM(i, j) = (*this)(i, j) * M(i, j);
		}
	}
	return tM;
}

template < typename T >
Matrix<T> Matrix<T>::elediv(Matrix<T>& M)
{
	Matrix<T> tM(NULL, m_nRows, m_nCols);
	if (m_nRows != M.m_nRows && m_nCols != M.m_nCols)
	{
		cout << "作对应元素除法的矩阵不同型\n" << endl;
		return tM;
	}
	for (int i = 1; i <= m_nRows; i++)
	{
		for (int j = 1; j <= m_nCols; j++)
		{
			tM(i, j) = (*this)(i, j) / M(i, j);
		}
	}
	return tM;
}


template < typename T >
Matrix<T>& Matrix<T>::swap_col(int c_i, int c_j)
{
	if (c_i<1 || c_j<1 || c_i>m_nCols || c_j>m_nCols)
	{
		cout << "用于交换的列索引超出矩阵范围\n";
		return *this;
	}
	if (c_i == c_j)//不做交换
	{
		return *this;
	}
	for (int i = 1; i <= m_nRows; i++)
	{
		T tkey = (*this)(i, c_i);
		(*this)(i, c_i) = (*this)(i, c_j);
		(*this)(i, c_j) = tkey;
	}
	return *this;
}

template < typename T >
Matrix<T>& Matrix<T>::swap_row(int r_i, int r_j)
{
	if (r_i<1 || r_j<1 || r_i>m_nRows || r_j>m_nRows)
	{
		cout << "用于交换的行索引超出矩阵范围\n";
		return *this;
	}
	if (r_i == r_j)//不做交换
	{
		return *this;
	}
	for (int j = 1; j <= m_nCols; j++)
	{
		T tkey = (*this)(r_i, j);
		(*this)(r_i, j) = (*this)(r_j, j);
		(*this)(r_j, j) = tkey;
	}
	return *this;
}


template < typename T >
Matrix<T>& Matrix<T>::eye(int n)
{
	m_nRows = m_nCols = n;
	delete[]m_pEntry;
	int length = m_nRows*m_nCols;
	m_pEntry = new T[length];
	memset(m_pEntry, 0, sizeof(T)*length);
	for (int i = 1; i <= m_nRows; i++)
	{
		(*this)(i, i) = 1;
	}
	return *this;
}

template < typename T >
void Matrix<T>::resize(int nRows, int nCols)
{
	if (m_pEntry != nullptr)
	{
		delete[] m_pEntry;
	}
	int length = nRows*nCols;
	m_pEntry = new T[length];
	memset(m_pEntry, 0, sizeof(T)*length);
	m_nRows = nRows;
	m_nCols = nCols;
}

template < typename T >
bool Matrix<T>::save(wchar_t* filename, wchar_t* filenameImag)//
{
	fstream file, fileImag;
	file.open(filename, ios::out | ios::trunc);

	if (file.fail())
	{
		return false;
	}
	if (is_same<complex<double>, T>::value)
	{
		fileImag.open(filenameImag, ios::out | ios::trunc);
		if (fileImag.fail())
		{
			return false;
		}
		for (int i = 1; i <= m_nRows; i++)
		{
			file << "   ";
			fileImag << "   ";
			for (int j = 1; j <= m_nCols; j++)
			{
				auto& tmat = *(Matrix<complex<double>>*)this;
				
				if (abs(tmat(i, j).real()) > 1e-25)
				{
					file << tmat(i, j).real() << "   ";
				}
				else
				{
					file << 0.0000 << "   ";
				}
				if (abs(tmat(i, j).imag()) > 1e-25)
				{
					fileImag << tmat(i, j).imag() << "   ";
				}
				else
				{
					fileImag << 0.0000 << "   ";
				}
			}
			file << endl;
			fileImag << endl;
		}
		fileImag.close();
	}
	else
	{
		for (int i = 1; i <= m_nRows; i++)
		{
			file << "   ";
			for (int j = 1; j <= m_nCols; j++)
			{
				file << (*this)(i, j) << "   ";
			}
			file << endl;
		}
	}
	
	file.close();
	return true;
}


template < typename T>
bool Matrix<T>::savebin(const wchar_t* filename, const wchar_t* filenameImag /*= nullptr*/)
{
	fstream file, fileImag;
	file.open(filename, ios::out | ios::trunc | ios::binary);

	if (file.fail())
	{
		cout << "写入文件" << filename << "失败\n";
		return false;
	}
	if (is_same<complex<double>, T>::value)
	{
		fileImag.open(filenameImag, ios::out | ios::trunc | ios::binary);
		if (fileImag.fail())
		{
			cout << "写入文件" << filenameImag << "失败\n";
			return false;
		}		
		//先写出矩阵的大小，int型
		file.write(reinterpret_cast<const char*>(&this->m_nRows), sizeof(int));
		file.write(reinterpret_cast<const char*>(&this->m_nCols), sizeof(int));
		fileImag.write(reinterpret_cast<const char*>(&this->m_nRows), sizeof(int));
		fileImag.write(reinterpret_cast<const char*>(&this->m_nCols), sizeof(int));
		//再写出矩阵的元素
		for (int i = 0; i < m_nRows * m_nCols; i++)
		{
			auto pEntry = (complex<double>*)GetEntry();
			double mreal = pEntry[i].real();
			double mimag = pEntry[i].imag();
			//写出实部
			file.write(reinterpret_cast<const char*>(&mreal), sizeof(double));
			//写出虚部
			fileImag.write(reinterpret_cast<const char*>(&mimag), sizeof(double));
		}
		fileImag.close();
	}
	else
	{
		//先写出矩阵的大小，int型
		file.write(reinterpret_cast<const char*>(&m_nRows), sizeof(int));
		file.write(reinterpret_cast<const char*>(&m_nCols), sizeof(int));
		//写出矩阵元素
		file.write(reinterpret_cast<const char*>(m_pEntry), m_nRows * m_nCols * sizeof(double));
	}

	file.close();
	return true;
}



template < typename T>
bool Matrix<T>::savebin(const char* filename, const char* filenameImag /*= nullptr*/)
{
	fstream file, fileImag;
	file.open(filename, ios::out | ios::trunc | ios::binary);

	if (file.fail())
	{
		cout << "写入文件" << filename << "失败\n";
		return false;
	}
	if (is_same<complex<double>, T>::value)
	{
		fileImag.open(filenameImag, ios::out | ios::trunc | ios::binary);
		if (fileImag.fail())
		{
			cout << "写入文件" << filenameImag << "失败\n";
			return false;
		}
		//先写出矩阵的大小，int型
		file.write(reinterpret_cast<const char*>(&this->m_nRows), sizeof(int));
		file.write(reinterpret_cast<const char*>(&this->m_nCols), sizeof(int));
		fileImag.write(reinterpret_cast<const char*>(&this->m_nRows), sizeof(int));
		fileImag.write(reinterpret_cast<const char*>(&this->m_nCols), sizeof(int));
		//再写出矩阵的元素
		for (int i = 0; i < m_nRows * m_nCols; i++)
		{
			auto pEntry = (complex<double>*)GetEntry();
			double mreal = pEntry[i].real();
			double mimag = pEntry[i].imag();
			//写出实部
			file.write(reinterpret_cast<const char*>(&mreal), sizeof(double));
			//写出虚部
			fileImag.write(reinterpret_cast<const char*>(&mimag), sizeof(double));
		}
		fileImag.close();
	}
	else
	{
		//先写出矩阵的大小，int型
		file.write(reinterpret_cast<const char*>(&m_nRows), sizeof(int));
		file.write(reinterpret_cast<const char*>(&m_nCols), sizeof(int));
		//写出矩阵元素
		file.write(reinterpret_cast<const char*>(m_pEntry), m_nRows * m_nCols * sizeof(double));
	}

	file.close();
	return true;
}

template < typename T >
bool Matrix<T>::load(wchar_t* filename)
{
	//暂时只能读取向量
	fstream file;
	std::string tmp;
	file.open(filename, ios::in);
	if (file.fail())
	{
		cout << "待加载数据文件打开失败\n";
		return false;
	}
	int nlines(0);//统计行数
	while (getline(file, tmp))
	{
		nlines++;
	}
	this->resize(nlines, 1);//产生向量
	file.close();
	//重新打开文件录入数据
	file.open(filename, ios::in);
	for (int i = 1; i <= nlines; i++)
	{
		double x;
		file >> x;
		(*this)(i) = x;
	}
	file.close();
	return true;
}

template < typename T>
bool Matrix<T>::loadbin(wchar_t* filename)
{
	// 打开包含二进制数据的文件
	std::ifstream file(filename, std::ios::binary);

	if (!file) {
		std::cout << "无法打开文件" << std::endl;
		return false;
	}

	int rows, cols;
	// 读取行数和列数
	file.read(reinterpret_cast<char*>(&rows), sizeof(int));
	file.read(reinterpret_cast<char*>(&cols), sizeof(int));

	// 计算数据的总数
	int numElements = rows * cols;
	this->resize(rows, cols);
	// 读取double数据到向量中
	file.read(reinterpret_cast<char*>(m_pEntry), numElements * sizeof(double));

	// 关闭文件
	file.close();
	return true;
}

template < typename T >
void Matrix<T>::Symmetrized(bool u2d)
{
	if (u2d)
	{
		//从上对称到下
		for (int i = 1; i <= m_nRows; i++)
		{
			for (int j = i + 1; j <= m_nCols; j++)
			{
				(*this)(j, i) = (*this)(i, j);
			}
		}
	}
	else
	{
		//从下对称到上
		for (int j = 1; j <= m_nCols; j++)
		{
			for (int i = j + 1; i <= m_nRows; i++)
			{
				(*this)(j, i) = (*this)(i, j);
			}
		}
	}
	
}

template < typename T >
int Matrix<T>::split_real_imag(Matrix<double>& matReal, Matrix<double>& matImag, char C)
{
	if (matReal.size(1) != m_nRows || matReal.size(2) != m_nCols)
	{
		matReal.resize(m_nRows, m_nCols);
		matImag.resize(m_nRows, m_nCols);
	}
	int length = m_nRows*m_nCols;
	if (is_same<complex<double>, T>::value)
	{
		auto& tmat = *(Matrix<complex<double>>*)this;
		for (int i = 1; i <= length; i++)
		{
			if (C=='A')
			{
				matReal(i) = tmat(i).real();
				matImag(i) = tmat(i).imag();
			}
			else if (C == 'R')
			{
				matReal(i) = tmat(i).real();
			}
			else if (C == 'I')
			{
				matImag(i) = tmat(i).imag();
			}
		}
	}
	return length;
}

template < typename T>
int Matrix<T>::split_real_imag(double* dReal, double* dimag, char C /*= 'A'*/)
{
	int length = m_nRows*m_nCols;
	if (is_same<complex<double>, T>::value)
	{
		auto &pe = (complex<double>*)(this->m_pEntry);
		for (int i = 0; i < length; i++)
		{
			if (C == 'A')
			{
				dReal[i] = pe[i].real();
				dimag[i] = pe[i].imag();
			}
			else if (C == 'R')
			{
				dReal[i] = pe[i].real();
			}
			else if (C == 'I')
			{
				dimag[i] = pe[i].imag();
			}
		}
	}
	return length;
}


template < typename T >
const T* const Matrix<T>::GetEntry()
{
	return m_pEntry;
}


template < typename T>
void Matrix<T>::SetEntry(T* entry)
{
	if (m_pEntry && m_blflag == false)
	{
		delete[] m_pEntry;
	}
	m_pEntry = entry;
	m_blflag = true;
}

template < typename T>
void Matrix<T>::SetEntry(int nrows, int ncols, T* entry)
{
	if (m_pEntry && m_blflag == false)
	{
		delete[] m_pEntry;
	}
	m_nRows = nrows;
	m_nCols = ncols;
	m_pEntry = entry;
	m_blflag = true;
}
template < typename T>
void Matrix<T>::Dense2Sparse(Mat_Sparse<T>& mat_spa)
{
	int spa_i = 0;
	for (int j = 1; j <= m_nCols; j++)
	{
		for (int i = 1; i <= m_nRows; i++)//列优先
		{
			T& data = (*this)(i, j);
			if (data != 0)
			{
				spa_i++;//统计非0元素的个数
			}
		}
	}
	mat_spa.resize(m_nRows, m_nCols, spa_i);
	spa_i = 0;
	for (int j = 1; j <= m_nCols; j++)
	{
		for (int i = 1; i <= m_nRows; i++)//列优先
		{
			T& data = (*this)(i, j);
			if (data != 0)
			{
				mat_spa.r_index[spa_i] = i;
				mat_spa.c_index[spa_i] = j;
				mat_spa.data[spa_i++] = data;
			}
		}
	}
}

template < typename T>
void Matrix<T>::Sparse2Dense(Mat_Sparse<T>& mat_spa)
{
	this->resize(mat_spa.m_nrows, mat_spa.m_ncols);
	int spa_i = mat_spa.m_nsize;
	for (int i = 0; i < spa_i; i++)
	{
		(*this)(mat_spa.r_index[i], mat_spa.c_index[i]) = mat_spa.data[i];
	}
}

template < typename T>
void Matrix<T>::Dense2Sparse(Matrix<T> mat_B, Mat_Sparse<T>& mat_spa, Mat_Sparse<T>& mat_spa_comp)
{
	int spa_i = 0;
	for (int j = 1; j <= m_nCols; j++)
	{
		for (int i = 1; i <= m_nRows; i++)//列优先
		{
			T& data = (*this)(i, j);
			if (data != 0)
			{
				spa_i++;//统计非0元素的个数
			}
		}
	}
	mat_spa.resize(m_nRows, m_nCols, spa_i);
	mat_spa_comp.Bindindex(mat_spa);//索引绑定，共用索引空间
	spa_i = 0;
	for (int j = 1; j <= m_nCols; j++)
	{
		for (int i = 1; i <= m_nRows; i++)//列优先
		{
			T& data = (*this)(i, j);
			if (data != 0)
			{
				mat_spa.r_index[spa_i] = i;
				mat_spa.c_index[spa_i] = j;
				mat_spa.data[spa_i] = data;
				mat_spa_comp.data[spa_i++] = mat_B(i, j);
			}
		}
	}
}


template < typename T>
template < typename Y>
void Matrix<T>::Sparse2Dense_potential(Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Y alpha)
{
	if (is_same<complex<Y>, T>::value)
	{
		Matrix<complex<Y>>* p_this = (Matrix<complex<Y>>*)this;
		p_this->resize(mat_spa.m_nrows, mat_spa.m_ncols);
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*p_this)(mat_spa.r_index[i], mat_spa.c_index[i]) = mat_spa.data[i] + alpha * mat_spa_comp.data[i];
		}
	}
	else if (is_same<T, Y>::value)
	{
		this->resize(mat_spa.m_nrows, mat_spa.m_ncols);
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*this)(mat_spa.r_index[i], mat_spa.c_index[i]) = mat_spa.data[i] + alpha * mat_spa_comp.data[i];
		}
	}
	else
	{
		cout << "读取的稀疏矩阵类型不符\n";
	}
}

template < typename T>
template<typename Y>
void Matrix<T>::Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Y alpha)
{
	if (r_end - r_beg + 1 != mat_spa.m_nrows || c_end - c_beg + 1 != mat_spa.m_ncols)
	{
		cout << "读取的稀疏矩阵和本矩阵指定子块的大小不一致\n";
		return;
	}
	if (is_same<complex<Y>, T>::value)
	{
		Matrix<complex<Y>>* p_this = (Matrix<complex<Y>>*)this;
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*p_this)(mat_spa.r_index[i] + r_beg - 1, mat_spa.c_index[i] + c_beg - 1) = mat_spa.data[i] + alpha * mat_spa_comp.data[i];
		}
	}
	else if (is_same<T, Y>::value)
	{
		this->resize(mat_spa.m_nrows, mat_spa.m_ncols);
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*this)(mat_spa.r_index[i] + r_beg - 1, mat_spa.c_index[i] + c_beg - 1) = mat_spa.data[i] + alpha * mat_spa_comp.data[i];
		}
	}
	else
	{
		cout << "读取的稀疏矩阵类型不符\n";
	}
}


template < typename T>
template<typename Y>
void Matrix<T>::Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, T alpha, T beta)
{
	if (r_end - r_beg + 1 != mat_spa.m_nrows || c_end - c_beg + 1 != mat_spa.m_ncols)
	{
		cout << "读取的稀疏矩阵和本矩阵指定子块的大小不一致\n";
		return;
	}
	if (is_same<complex<Y>, T>::value)
	{
		Matrix<complex<Y>>* p_this = (Matrix<complex<Y>>*)this;
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*p_this)(mat_spa.r_index[i] + r_beg - 1, mat_spa.c_index[i] + c_beg - 1) = (mat_spa.data[i] + alpha * mat_spa_comp.data[i]) * beta;
		}
	}
	else if (is_same<T, Y>::value)
	{
		this->resize(mat_spa.m_nrows, mat_spa.m_ncols);
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*this)(mat_spa.r_index[i] + r_beg - 1, mat_spa.c_index[i] + c_beg - 1) = (mat_spa.data[i] + alpha * mat_spa_comp.data[i]) * beta;
		}
	}
	else
	{
		cout << "读取的稀疏矩阵类型不符\n";
	}
}

template < typename T>
template<typename Y>
void Matrix<T>::Sparse2Dense_potential(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, Mat_Sparse<Y>& mat_spa_ksg, T alpha, T beta)
{
	if (r_end - r_beg + 1 != mat_spa.m_nrows || c_end - c_beg + 1 != mat_spa.m_ncols)
	{
		cout << "读取的稀疏矩阵和本矩阵指定子块的大小不一致\n";
		return;
	}
	if (is_same<complex<Y>, T>::value)
	{
		Matrix<complex<Y>>* p_this = (Matrix<complex<Y>>*)this;
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*p_this)(mat_spa.r_index[i] + r_beg - 1, mat_spa.c_index[i] + c_beg - 1) = (mat_spa.data[i] + alpha * mat_spa_comp.data[i]) * beta;
		}
		spa_i = mat_spa_ksg.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*p_this)(mat_spa_ksg.r_index[i] + r_beg - 1, mat_spa_ksg.c_index[i] + c_beg - 1) -= mat_spa_ksg.data[i] * beta;
		}
	}
	else if (is_same<T, Y>::value)
	{
		this->resize(mat_spa.m_nrows, mat_spa.m_ncols);
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*this)(mat_spa.r_index[i] + r_beg - 1, mat_spa.c_index[i] + c_beg - 1) = (mat_spa.data[i] + alpha * mat_spa_comp.data[i]) * beta;
		}
		spa_i = mat_spa_ksg.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*this)(mat_spa_ksg.r_index[i] + r_beg - 1, mat_spa_ksg.c_index[i] + c_beg - 1) -= mat_spa_ksg.data[i] * beta;
		}
	}
	else
	{
		cout << "读取的稀疏矩阵类型不符\n";
	}
}


template < typename T>
template<typename Y>
void Matrix<T>::Sparse2Dense_potential_add(int r_beg, int c_beg, int r_end, int c_end, Mat_Sparse<Y>& mat_spa, Mat_Sparse<Y>& mat_spa_comp, T alpha, T beta)
{
	if (r_end - r_beg + 1 != mat_spa.m_nrows || c_end - c_beg + 1 != mat_spa.m_ncols)
	{
		cout << "读取的稀疏矩阵和本矩阵指定子块的大小不一致\n";
		return;
	}
	if (is_same<complex<Y>, T>::value)
	{
		Matrix<complex<Y>>* p_this = (Matrix<complex<Y>>*)this;
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*p_this)(mat_spa.r_index[i] + r_beg - 1, mat_spa.c_index[i] + c_beg - 1) += (mat_spa.data[i] + alpha * mat_spa_comp.data[i]) * beta;
		}
	}
	else if (is_same<T, Y>::value)
	{
		this->resize(mat_spa.m_nrows, mat_spa.m_ncols);
		int spa_i = mat_spa.m_nsize;
		for (int i = 0; i < spa_i; i++)
		{
			(*this)(mat_spa.r_index[i] + r_beg - 1, mat_spa.c_index[i] + c_beg - 1) += (mat_spa.data[i] + alpha * mat_spa_comp.data[i]) * beta;
		}
	}
	else
	{
		cout << "读取的稀疏矩阵类型不符\n";
	}
}

template < typename T>
void Matrix<T>::toinv()
{	
	if (std::is_same<complex<double>,T>::value == 0)
	{
		//该toinv方法只针对复矩阵
		return;
	}
	Matrix<complex<double>>& Mat = *(Matrix<complex<double>>*)this;
	int i, j, k;
	int ret = 0;
	//! 必须是方阵
	if (m_nRows != m_nCols)
	{
		cout << "求逆的矩阵必须是方阵\n";
		return;
	}
	int* is = new int[m_nRows];
	int* js = new int[m_nCols];
	double max(0);
	T max_c;
	for (k = 0; k < m_nRows; k++)
	{
		//step 1, 全选主元
		max = 0;
		is[k] = k;
		js[k] = k;

		for (i = k; i < m_nRows; i++)
		{
			for (j = k; j < m_nCols; j++)
			{
				complex<double> tkey = Mat(i + 1, j + 1);//强制转换
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
		{	// 无逆矩阵
			cout << "toinv：矩阵无逆矩\n";
			return;
		}

		//交换
		if (is[k] != k)
		{
			Mat.swap_row(k + 1, is[k] + 1);
		}
		if (js[k] != k)
		{
			Mat.swap_col(k + 1, js[k] + 1);
		}
		
		Mat(k + 1, k + 1) = conj(Mat(k + 1, k + 1)) / (abs(Mat(k + 1, k + 1)) * abs(Mat(k + 1, k + 1)));//1/(a+bi)

		for (j = 0; j < m_nCols; j++)
		{
			if (j != k)
				Mat(k + 1, j + 1) *= Mat(k + 1, k + 1);
		}
		for (i = 0; i < m_nRows; i++)
		{
			if (i != k)
			{
				for (j = 0; j < m_nCols; j++)
				{
					if (j != k)
						Mat(i + 1, j + 1) -= Mat(i + 1, k + 1) * Mat(k + 1, j + 1);
				}
			}
		}
		for (i = 0; i < m_nRows; i++)
		{
			if (i != k)
				Mat(i + 1, k + 1) *= -Mat(k + 1, k + 1);
		}

	}

	//恢复
	//本来 row <-> is[k], column <-> js[k]
	//恢复时：row <-> js[k], column <-> is[k]
	for (k = m_nRows - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			Mat.swap_row(k + 1, js[k] + 1);
		}
		if (is[k] != k)
		{
			Mat.swap_col(k + 1, is[k] + 1);
		}
	}
	delete[] is;
	delete[] js;
}