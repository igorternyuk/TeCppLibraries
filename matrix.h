#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <time.h>
#include <complex>
using namespace std;
template <typename T>
class Matrix
{
private:
    int rowCount, colCount;
    T **values;
    bool isQuadratic;  
public:
    Matrix();
    Matrix(const int&, bool);                    // Единичная или нулевая матрица
    Matrix(const int&, const int&, T);           // Прямоугольная матрица
    Matrix(const string&);                       // Прямоугольная матрица
    Matrix(const Matrix&);                       // Конструктор копирования
    ~Matrix();
    void fillByRandomNumbers(const int&, const int&);
    void readMatrix(const string &fileName);
    void setValueAt(int row, int col, T value);
    T getValueAt(int row, int col) const;
    T *getRow(const int&) const; // Получить строку по индексу
    T *getColumn(const int&) const; // Получить строку по индексу
    T maxElement() const;
    T minElement() const;
    int getRowCount() const; // Получить строку по индексу
    int getColumnCount() const; // Получить строку по индексу
    void swapRows(const int&, const int&);
    void swapColumns(const int&, const int&);
    void swapElements(const int&, const int&, const int&, const int&);
    Matrix appendMatrixFromRightSide(const Matrix&) const;
    void deleteRow(const int&);
    void deleteColumn(const int&);
    Matrix transpose() const;
    Matrix inverse() const;
    Matrix getMinorAt(const int&, const int&) const;
    int rang() const;
    T determinant() const;
    Matrix power(int) const;
    Matrix operator+(const Matrix&) const;
    Matrix operator+(const T&);
    Matrix operator-(const Matrix&) const;
    Matrix operator-(const T&);
    Matrix operator*(const Matrix&) const;
    Matrix operator*(const T&);
    Matrix operator/(const T&);
    Matrix &operator+=(const Matrix&);
    Matrix &operator+=(const T&);
    Matrix &operator-=(const Matrix&);
    Matrix &operator-=(const T&);
    Matrix &operator*=(const Matrix&);
    Matrix &operator*=(const T&);
    Matrix &operator/=(const T&);
    Matrix &operator=(const Matrix&);
    bool operator==(const Matrix&) const;
    bool operator!=(const Matrix& matrix) const;
    T *operator[](const int&);
};
template <typename T>
Matrix<T>::Matrix()
{
  rowCount = 0;
  colCount = 0;
  isQuadratic = false;
  values = nullptr;
}
template <typename T>
Matrix<T>::Matrix(const int& size, bool isUnitMatrix)       // Единичная или нулевая матрица
{
    if(size > 0)
       {
          values = new T*[size];
          for(int i = 0; i < size; ++i)
                values[i] = new T[size];
          for(int i = 0; i < size; ++i)
              for(int j = 0; j < size; ++j)
                values[i][j] = (isUnitMatrix && i == j) ? 1 : 0;
          this->rowCount = size;
          this->colCount = size;
       }else
       {
        values = nullptr;
        this->rowCount = 0;
        this->colCount = 0;
       }
    isQuadratic = true;
}
template <typename T>
Matrix<T>::Matrix(const int& rowCount, const int& colCount, T value) // Прямоугольная матрица
{
   if(rowCount > 0 && colCount > 0)
      {
         values = new T*[rowCount];
         for(int i = 0; i < rowCount; ++i)
               values[i] = new T[colCount];
         for(int i = 0; i < rowCount; ++i)
             for(int j = 0; j < colCount; ++j)
               values[i][j] = value;
         this->rowCount = rowCount;
         this->colCount = colCount;
         isQuadratic = rowCount == colCount ? true : false;
      }else
      {
         values = nullptr;
         this->rowCount = 0;
         this->colCount = 0;
         isQuadratic = false;
     }
}
template <typename T>
Matrix<T>::Matrix(const string &fileName)
{
  ifstream file;
  int row, col;
  file.open(fileName.c_str());
  file >> row >> col;
  if(row > 0 && col > 0)
  {
      rowCount = row;
      colCount = col;
      if(rowCount == colCount)
          isQuadratic = true;
      else
          isQuadratic = false;
      values = new T*[rowCount];
      for(int i = 0; i < rowCount; ++i)
            values[i] = new T[colCount];
      for(int i = 0; i < rowCount; ++i)
          for(int j = 0; j < colCount; ++j)
          {
              T tmp;
              file >> tmp;
              setValueAt(i, j, tmp);
          }
  }else
  {
      rowCount = 0;
      colCount = 0;
      isQuadratic = false;
      values = nullptr;
  }
  file.close();
}
template <typename T>
Matrix<T>::Matrix(const Matrix& matrix): rowCount(matrix.rowCount),
      colCount(matrix.colCount), isQuadratic(matrix.isQuadratic)    // Конструктор копирования
{
   if(this != &matrix)
   {
      values = new T*[rowCount];
      for(int j = 0; j < rowCount; ++j)
          values[j] = new T[colCount];
      for(int i = 0; i < rowCount; ++i)
          for(int j = 0; j < colCount; ++j)
              values[i][j] = matrix.values[i][j];
   }
}
template <typename T>
Matrix<T>::~Matrix()
{
    if(values != nullptr)
    {
      for(int i = 0; i < rowCount; ++i)
        if(values[i] != nullptr)
          delete[] values[i];
      delete[] values;
    }
}
template <typename T>
void Matrix<T>::fillByRandomNumbers(const int& initialValue, const int& finalValue)
{
   srand(time(NULL));
   for(int i = 0; i < rowCount; ++i)
    for(int j = 0; j < colCount; ++j)
      values[i][j] = static_cast<T>(initialValue + rand() % (finalValue - initialValue));
}
template <typename T>
void Matrix<T>::setValueAt(int row, int col, T value)
{
  if(values != nullptr && row >=0 && row < this->rowCount
     && col >=0 && col < this->colCount)
     values[row][col] = value;
}
template <typename T>
T Matrix<T>::getValueAt(int row, int col) const
{
    return (values != nullptr && row >=0 && row < rowCount &&
       col >=0 && col < colCount) ? values[row][col] : 0.0;
}
template <typename T>
T* Matrix<T>::getRow(const int& row) const // Получить строку по индексу
{
   if(values != nullptr && row >=0 && row < rowCount)
   {
      T *pRow = new T[colCount];
      for(int j = 0; j < colCount; ++j)
          pRow[j] = values[row][j];
     return pRow;
   }else
     return nullptr;
}
template <typename T>
T* Matrix<T>::getColumn(const int& col) const // Получить строку по индексу
{
    if(values != nullptr && col >=0 && col < colCount)
    {
       T *pCol = new T[rowCount];
       for(int i = 0; i < rowCount; ++i)
           pCol[i] = values[i][col];
      return pCol;
    }else
      return nullptr;
}
template <typename T>
int Matrix<T>::getRowCount() const
{
   return rowCount;
}
template <typename T>
int Matrix<T>::getColumnCount() const
{
   return colCount;
}
template <typename T>
void Matrix<T>::swapRows(const int& row1, const int& row2)
{
  if(values != nullptr && row1 >= 0 && row1 < rowCount &&
     row2 >= 0 && row2 < rowCount && row1 != row2)
     {
        T *pTmp = getRow(row1);
        for(int j = 0; j < colCount; ++j)
        {
            values[row1][j] = values[row2][j];
            values[row2][j] = pTmp[j];
        }
        delete pTmp;
     }
}
template <typename T>
void Matrix<T>::swapColumns(const int& col1, const int& col2)
{
   if(values != nullptr && col1 >= 0 && col1 < colCount &&
       col2 >= 0 && col2 < colCount && col1 != col2)
       {
          T *pTmp = getColumn(col1);
          for(int i = 0; i < rowCount; ++i)
          {
              values[i][col1] = values[i][col2];
              values[i][col2] = pTmp[i];
          }
          delete pTmp;
   }
}
template <typename T>
void Matrix<T>::swapElements(const int &row1, const int &col1,
                          const int &row2, const int &col2)
{
    if(values != nullptr && row1 >= 0 && row1 < rowCount &&
       col1 >=0 && col1 < colCount && row2 >= 0 && row2 < rowCount &&
            col2 >=0 && col2 < colCount && (row1 != row2 || col1 != col2))
    {
       T tmp = values[row1][col1];
       values[row1][col1] = values[row2][col2];
       values[row2][col2] = tmp;
    }
}
template <typename T>
T Matrix<T>::maxElement() const
{
    if(values != nullptr)
       {
          T maxValue = values[0][0];
          for(int i = 0; i < rowCount; ++i)
            for(int j = 0; j < colCount; ++j)
             {
                if(values[i][j] > maxValue)
                    maxValue = values[i][j];
             }
          return maxValue;
       } else return static_cast<T>(0);
}
template <typename T>
T Matrix<T>::minElement() const
{
    if(values != nullptr)
       {
          T minValue = values[0][0];
          for(int i = 0; i < rowCount; ++i)
            for(int j = 0; j < colCount; ++j)
             {
                if(values[i][j] < minValue)
                    minValue = values[i][j];
             }
          return minValue;
       } else return static_cast<T>(0);
}
template <typename T>
Matrix<T> Matrix<T>::appendMatrixFromRightSide(const Matrix& matrix) const
{
    if(isQuadratic && rowCount == matrix.rowCount)
    {
        Matrix mtr(rowCount, colCount + matrix.colCount, static_cast<T>(0));
        for(int i = 0; i < mtr.getRowCount(); ++i)
            for(int j = 0; j < mtr.getColumnCount(); ++j)
            {
                if(j < colCount)
                    mtr[i][j] = values[i][j];
                else if(j - i == rowCount)
                    mtr[i][j] = 1.0;
                else
                    mtr[i][j] = 0.0;
            }
        return mtr;
    }
    else
        return *this;
}
template <typename T>
void Matrix<T>::deleteRow(const int& row)
{
   if(values != nullptr && row >= 0 && row < rowCount)
   {
       for(int i = row; i < rowCount - 1; ++i)
           for(int j = 0; j < colCount; ++j)
               values[i][j] = values[i + 1][j];
       delete[] values[rowCount - 1];
       --rowCount;
   }
}
template <typename T>
void Matrix<T>::deleteColumn(const int& col)
{
    if(values != nullptr && col >= 0 && col < colCount)
    {
        for(int j = col; j < colCount - 1; ++j)
            for(int i = 0; i < rowCount; ++i)
                values[i][j] = values[i][j + 1];
        --colCount;
    }
}
template <typename T>
Matrix<T> Matrix<T>::transpose() const
{
    if(values != nullptr)
    {
        Matrix matrix(colCount, rowCount, static_cast<T>(0));
        for(int i = 0; i < matrix.getRowCount(); ++i)
            for(int j = 0; j < matrix.getColumnCount(); ++j)
                matrix[i][j] = values[j][i];
        return matrix;
    }else return *this;
}
template <typename T>
Matrix<T> Matrix<T>::inverse() const
{
    if(isQuadratic && values != nullptr)
    {
       Matrix inverseMatrix(rowCount, rowCount, static_cast<T>(0));
       T det = this->determinant();
       if (det == static_cast<T>(0)) return *this;
       for(int i = 0; i < inverseMatrix.getRowCount(); ++i)
           for(int j = 0; j < inverseMatrix.getColumnCount(); ++j)
               inverseMatrix[i][j] = pow(-1, i + j) * this->getMinorAt(i, j).determinant();
       inverseMatrix = inverseMatrix.transpose();
       inverseMatrix /= det;
       return inverseMatrix;
    } else
        return *this;
}
template <typename T>
Matrix<T> Matrix<T>::getMinorAt(const int & row, const int &col) const
{
    Matrix mtr = *this;
    mtr.deleteRow(row);
    mtr.deleteColumn(col);
    return mtr;
}
template <typename T>
int Matrix<T>::rang() const
{
    Matrix tmpMatrix = *this;
    int rang = 0;
    int limit = min(tmpMatrix.getRowCount(), tmpMatrix.getColumnCount());
    for(int i = 0; i < limit; i++)
    {
      T kk = 0;
      if(tmpMatrix.values[i][i] != static_cast<T>(0))   // если первый эл-т на главной диагонали не равен нулю делим строку на него
         kk = tmpMatrix.values[i][i];
      else
      {
          /* Если первый эл-т на главной диагонали равен нулю то перебирам в цикле
           * все элементы текущего столбца и если такой находится то складываем его строку с текущей
           * а если такого элемента нет значит все элементи вниз по столбцу нулевые и можно переходить
           * к следующей строке
           */
          bool flag = false;
          for(int ii = i; ii < limit; ++ii)
              if(tmpMatrix.values[ii][i] != static_cast<T>(0))
              {
                 kk = tmpMatrix.values[ii][i];
                 for(int col = i; col < tmpMatrix.getColumnCount(); col++)
                       tmpMatrix.values[i][col] += tmpMatrix.values[ii][col];
                 flag = true;
                 break;
              }
          if(!flag)
            continue;
      }
      for(int j = i; j < tmpMatrix.getColumnCount(); j++)  // Если достигли последней строки прекращаем цикл
           tmpMatrix[i][j] /= kk;
       //Складываем текущую строку умноженную на первые эл-ты следующих строк взятые с минусом
       // со всеми строками ниже текущей и получаем в i-м столбце нули ниже главной диагонали
      for(int row = i + 1; row < limit; row++)
      {
        T coef = -tmpMatrix.values[row][i];
        for(int col = i; col < tmpMatrix.getColumnCount(); col++)
              tmpMatrix.values[row][col] += coef * tmpMatrix.values[i][col];
      }
    }
    // Подсчитываем количество нелулевых строк полученной ступенчатой матрицы
    // Это и будет рангом матрицы
    for(int i = 0; i < limit; i++)
        for(int j = 0; j < tmpMatrix.getColumnCount(); j++)
         {
            if(tmpMatrix[i][j] != static_cast<T>(0))
            {
                rang++;
                break;
            }
          }
    return rang;
}
template <typename T>
T Matrix<T>::determinant() const
{
   if(isQuadratic && values != nullptr)
   {
       if(rowCount == 1) return values[0][0];
               else  if(rowCount == 2)
           return values[0][0] * values[1][1] -
                  values[0][1] * values[1][0];
       else
       {
          Matrix tmpMatrix = *this;
          T det = 0;
          for(int i = 0; i < rowCount; ++i)
              det += pow(-1, i) * tmpMatrix.values[0][i] * tmpMatrix.getMinorAt(0, i).determinant();
          return det;
       }
   }
   else
   return static_cast<T>(0);
}
template <typename T>
Matrix<T> Matrix<T>::power(int exponent) const
{
    if (isQuadratic && values != nullptr)
    {
        Matrix matrix = *this;
        Matrix power(matrix.getRowCount(), true);
        for(int i = 1; i <= exponent; ++i)
            power *= matrix;
        return power;
    }else
        return *this;
}
template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& matrix) const
{
    if(values == nullptr || matrix.values == nullptr ||
       rowCount != matrix.rowCount || colCount != matrix.colCount)
         return *this;
    else
    {
      Matrix result(rowCount, colCount, static_cast<T>(0));
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           result.values[i][j] = values[i][j] + matrix.values[i][j];
      return result;
    }
}
template <typename T>
Matrix<T> Matrix<T>::operator+(const T& increment)
{
    if(values == nullptr)
       return *this;
    else
    {
      Matrix result(rowCount, colCount, static_cast<T>(0));
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           result.values[i][j] = values[i][j] + increment;
      return result;
    }
}
template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix& matrix) const
{
    if(values == nullptr || matrix.values == nullptr ||
       rowCount != matrix.rowCount || colCount != matrix.colCount)
        return *this;
    else
    {
      Matrix result(rowCount, colCount, static_cast<T>(0));
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           result.values[i][j] = values[i][j] - matrix.values[i][j];
      return result;
    }
}
template <typename T>
Matrix<T> Matrix<T>::operator-(const T& decrement)
{
    if(values == nullptr)
        return *this;
    else
    {
      Matrix result(rowCount, colCount, static_cast<T>(0));
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           result.values[i][j] = values[i][j] - decrement;
      return result;
    }
}
template <typename T>
Matrix<T> Matrix<T>:: operator*(const Matrix& matrix) const
{
    if(values == nullptr || matrix.values == nullptr ||
       colCount != matrix.rowCount)
        return *this;
    else
    {
        Matrix<T> result(rowCount, matrix.colCount, static_cast<T>(0));
        T sum = static_cast<T>(0);
        for(int i = 0; i < rowCount; ++i)             // Перебор по строкам первой матрицы
            for(int j = 0; j < matrix.colCount; ++j)  // Перебор по столбцам второй матриці
            {
                sum = static_cast<T>(0);
                for(int k = 0; k < colCount; ++k)
                    sum += values[i][k] * matrix.values[k][j]; // Накопление произведений елементов i-й строки 1-й матрицы и j-го столбца 2-й матрицы
                result.values[i][j] = sum;
            }
      return result;
    }
}
template <typename T>
Matrix<T> Matrix<T>:: operator*(const T& scale)
{
    if(values == nullptr)
        return *this;
    else
    {
      Matrix result(rowCount, colCount, static_cast<T>(0));
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           result.values[i][j] = values[i][j] * scale;
      return result;
    }
}
template <typename T>
Matrix<T> Matrix<T>:: operator/(const T& scale)
{
    if(values == nullptr)
        return *this;
    else
    {
      Matrix<T> result(rowCount, colCount, static_cast<T>(0));
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           result.values[i][j] = values[i][j] / scale;
      return result;
    }
}
template <typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix& matrix)
{
    if(values == nullptr || rowCount != matrix.rowCount ||
       colCount != matrix.colCount || matrix.values == nullptr)
        return *this;
    else
    {
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           values[i][j] += matrix.values[i][j];
      return *this;
    }
}
template <typename T>
Matrix<T>& Matrix<T>::operator+=(const T& increment)
{
    if(values == nullptr)
       return *this;
   else
    {
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           values[i][j] += increment;
      return *this;
    }
}
template <typename T>
Matrix<T>& Matrix<T>:: operator-=(const Matrix& matrix)
{
    if(values == nullptr || rowCount != matrix.rowCount ||
       colCount != matrix.colCount || matrix.values == nullptr)
       return *this;
    else
    {
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           values[i][j] -= matrix.values[i][j];
      return *this;
    }
}
template <typename T>
Matrix<T>& Matrix<T>::operator-=(const T& decrement)
{
    if(values == nullptr)
        return *this;
    else
    {
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           values[i][j] -= decrement;
      return *this;
    }
}
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix& matrix)
{
    if(values == nullptr || matrix.values == nullptr || colCount != matrix.rowCount)
        return *this;
    else
    {
        Matrix result(matrix.rowCount, colCount, static_cast<T>(0));
        T sum = static_cast<T>(0);
        for(int i = 0; i < matrix.rowCount; ++i)              // Перебор по строкам первой матрицы
            for(int j = 0; j < colCount; ++j)                  // Перебор по столбцам второй матриці
            {
                sum = static_cast<T>(0);
                for(int k = 0; k < rowCount; ++k)
                    sum += matrix.values[i][k] * values[k][j]; // Накопление произведений елементов i-й строки 1-й матрицы и j-го столбца 2-й матрицы
                    result[i][j] = sum;
            }
      return *this = result;
    }
}
template <typename T>
Matrix<T>& Matrix<T>::operator*=(const T& scale)
{
    if(values == nullptr)
    {
        return *this;
    }else
    {
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           values[i][j] *= scale;
      return *this;
    }
}
template <typename T>
Matrix<T>& Matrix<T>::operator/=(const T& scale)
{
    if(values == nullptr)
    {
        return *this;
    }else
    {
      for(int i = 0; i < rowCount; ++i)
        for(int j = 0; j < colCount; ++j)
           values[i][j] /= scale;
      return *this;
    }
}
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix& matrix)
{
    if(&matrix != this)
    {
       rowCount = matrix.rowCount;
       colCount = matrix.colCount;
       isQuadratic = matrix.isQuadratic;
       if(values != nullptr)
       {
          for(int i = 0; i < rowCount; ++i)
              delete[] values[i];
          delete[] values;
       }
       values = new T*[rowCount];
       for(int j = 0; j < colCount; ++j)
           values[j] = new T[colCount];
       for(int i = 0; i < rowCount; ++i)
           for(int j = 0; j < colCount; ++j)
               values[i][j] = matrix.values[i][j];
    }
    return *this;
}
template <typename T>
Matrix<T>& operator++(Matrix<T> &mtr)
{
    return mtr += static_cast<T>(1);
}
template <typename T>
Matrix<T> operator++(Matrix<T> &mtr, int)
{
    Matrix<T> old = mtr;
    mtr += static_cast<T>(1);
    return old;
}
template <typename T>
Matrix<T>& operator--(Matrix<T> &mtr)
{
    return mtr -= static_cast<T>(1);
}
template <typename T>
Matrix<T> operator--(Matrix<T> &mtr, int)
{
    Matrix<T> old = mtr;
    mtr -= static_cast<T>(1);
    return old;
}
template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& matrix) const
{
    if(rowCount == matrix.rowCount && colCount == matrix.colCount)
    {
       for(int i = 0; i < rowCount; ++i)
           for(int j = 0; j < colCount; ++j)
           {
               if(values[i][j] != matrix.values[i][j])
                   return false;
           }
       return true;
    } else return false;
}
template <typename T>
bool Matrix<T>::operator!=(const Matrix<T>& matrix) const
{
    if(rowCount == matrix.rowCount && colCount == matrix.colCount)
    {
       for(int i = 0; i < rowCount; ++i)
           for(int j = 0; j < colCount; ++j)
           {
               if(values[i][j] != matrix.values[i][j])
                   return true;
           }
       return false;
    } else return false;
}
template <typename T>
T* Matrix<T>::operator[](const int& row)
{
    if(row >= 0 && row < rowCount)
        return values[row];
    else
        return nullptr;
}
template <typename T>
istream& operator>>(istream &stream, Matrix<T> &matrix)
{
    for(int i = 0; i < matrix.getRowCount(); ++i)
        for(int j = 0; j < matrix.getColumnCount(); ++j)
        {
            T tmp;
            stream >> tmp;
            matrix.setValueAt(i, j, tmp);
        }
    return stream;
}
template <typename T>
ostream& operator<<(ostream &stream, Matrix<T> &matrix)
{
    for(int i = 0; i < matrix.getRowCount(); ++i)
    {
        for(int j = 0; j < matrix.getColumnCount(); ++j)
           stream << matrix.getValueAt(i, j) << " ";
           stream << endl;
    }
    return stream;
}
template <typename T>
ifstream& operator>>(ifstream &stream, Matrix<T> &matrix)
{
    for(int i = 0; i < matrix.getRowCount(); ++i)
        for(int j = 0; j < matrix.getColumnCount(); ++j)
        {
            T tmp;
            stream >> tmp;
            matrix.setValueAt(i, j, tmp);
        }
    return stream;
}
template <typename T>
ofstream& operator<<(ofstream &stream, Matrix<T> &matrix)
{
    for(int i = 0; i < matrix.rowCount; ++i)
    {
        for(int j = 0; j < matrix.colCount; ++j)
            stream << matrix.getValueAt(i, j) << " ";
        stream << endl;
    }
    return stream;
}
#endif // MATRIX_H
