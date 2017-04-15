#ifndef POLYNOM_H
#define POLYNOM_H
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <time.h>
#include <complex>
#include <fstream>
#include <sstream>
#include <iomanip>
using namespace std;
template<class T>
class Polynom
{
private:
    int degree;
    vector<T> coefficients;
public:

    Polynom(int, T);                          // Парметры конструктора - степень полинома и заполнитель для коэффициентов
    Polynom(const string&fileName);           // Параметр конструктора - имя файла из которого будут считываться коэффициенты
    Polynom(std::vector<T> coefficients);     // Параметр конструктора - вектор коэффициентов полинома
    void readPolynom(const string &fileName); // Считать полином из файла
    void setValueAt(int i, const T&value);          // Установить значение коэффициента элемента степени i
    T getCoefAt(int i) const;                 // Считать значение коэффициента элемента степени i
    int getDegree() const;                       // Получить степень полинома
    T calculateValue(const T&) const;
    Polynom power(int) const;
    std::string convertToString(int prec) const;
    void showInConsole(int prec) const;
    void writeToFile(string fileName, int prec) const;
    bool operator==(const Polynom&) const;
    bool operator!=(const Polynom& Polynom) const;
    //T* operator[](const int& index);
    Polynom operator+(const Polynom&) const;
    Polynom operator+(const T&);
    Polynom operator-(const Polynom&) const;
    Polynom operator-(const T&);
    Polynom operator*(const Polynom&) const;
    Polynom operator*(const T&);
    Polynom operator/(const T&);
    Polynom &operator+=(const Polynom&);
    Polynom &operator+=(const T&);
    Polynom &operator-=(const Polynom&);
    Polynom &operator-=(const T&);
    Polynom &operator*=(const Polynom&);
    Polynom &operator*=(const T&);
    Polynom &operator/=(const T&);
    Polynom &operator=(const Polynom&);
};

template<class T>
Polynom<T>::Polynom(int degree, T value)
{
    if(degree < 0) degree = 0;
    else
    this->degree = degree;
    for(int i = 0; i <= degree; ++i)
        coefficients.push_back(value);
}

template<class T>
Polynom<T>::Polynom(const string& fileName)
{
   std::ifstream fi(fileName);
   fi >> degree;
   if(degree < 0) degree = 0;
   for(int i = 0; i <= degree; ++i)
   {
       T value;
       fi >> value;
       coefficients.push_back(value);
   }
   fi.close();
}

template<class T>
Polynom<T>::Polynom(std::vector<T> coefficients)
{
    degree = coefficients.size() - 1;
    this->coefficients = coefficients;
}

template<class T>
void Polynom<T>::readPolynom(const string &fileName)
{
    std::ifstream fi(fileName);
    fi >> degree;
    for(int i = 0; i <= degree; ++i)
    {
        T value;
        fi >> value;
        coefficients.push_back(value);
    }
}

template<class T>
void Polynom<T>::setValueAt(int i, const T &value)
{
    if(i >= 0 && i < coefficients.size())
      coefficients[i] = value;
}

template<class T>
T Polynom<T>::getCoefAt(int i) const
{
    if(i >= 0 && i < coefficients.size())
        return coefficients[i];
    else if(i >= coefficients.size())
        return static_cast<T>(0.0);
    else
        return std::numeric_limits<T>::quiet_NaN();
}

template<class T>
int Polynom<T>::getDegree() const
{
    return degree;
}

template<class T>
T Polynom<T>::calculateValue(const T &value) const
{
    T result = static_cast<T>(0);
    for(int i = 0; i <= degree; ++i)
        result += getCoefAt(i) * pow(value, i);
    return result;
}

template<class T>
Polynom<T> Polynom<T>::power(int exponent) const
{
    Polynom<T> product(0, static_cast<T>(1));
    Polynom<T> polynom = *this;
    for(int i = 1; i <= exponent; ++i)
        product *= polynom;
    return product;
}

template<class T>
void Polynom<T>::showInConsole(int prec) const
{
    std::cout << convertToString(prec) << std::endl;
}

template<class T>
std::string Polynom<T>::convertToString(int prec) const
{
    std::stringstream stream;
    stream << setprecision(prec) << getCoefAt(degree) << " * X**" << degree;
    for(int i = degree - 1; i > 0; --i)
    {
        if(getCoefAt(i) > 0)
            stream << " + " << getCoefAt(i) << " * X**" << i;
        else if(getCoefAt(i) < 0)
            stream << " - " << fabs(getCoefAt(i)) << " * X**" << i;
        else
            continue;
    }
    if(getCoefAt(0) > 0)
        stream << " + " << setprecision(prec) << getCoefAt(0);
    else if(getCoefAt(0) < 0)
        stream << " - " << setprecision(prec) << fabs(getCoefAt(0));
    return stream.str();
}

template<class T>
void Polynom<T>::writeToFile(std::string fileName, int prec) const
{
    std::ofstream of(fileName);
    of << convertToString(prec);
    of.close();
}

template<class T>
bool Polynom<T>::operator==(const Polynom &pol) const
{
    if(degree != pol.getDegree()) return false;
    bool isEqual = true;
    for(int i = 0; i < degree; ++i)
    {
        if(getCoefAt(i) != pol.getCoefAt(i))
        {
            isEqual = false;
            break;
        }
    }
    return isEqual;
}

template<class T>
bool Polynom<T>::operator!=(const Polynom &pol) const
{
    if(degree != pol.getDegree()) return true;
    bool isNotEqual = false;
    for(int i = 0; i < degree; ++i)
    {
        if(getCoefAt(i) != pol.getCoefAt(i))
        {
            isNotEqual = true;
            break;
        }
    }
    return isNotEqual;
}

/*template<class T>
T* Polynom<T>::operator[](const int& index)
{
    //if(index < 0) return std::numeric_limits<double>::quiet_NaN();
    if(index < 0 && index > getDegree()) return static_cast<T>(0.0);
    return coefficients + index;
}*/

template<class T>
Polynom<T> Polynom<T>::operator+(const Polynom& pol) const
{
    int resDegree = std::max(degree, pol.degree);
    Polynom<T> summa(resDegree, static_cast<T>(0));
    for(int i = 0; i <= resDegree; ++i)
        summa.coefficients[i] = getCoefAt(i) + pol.getCoefAt(i);
    return summa;
}

template<class T>
Polynom<T> Polynom<T>::operator+(const T& value)
{
    Polynom<T> summa(degree, static_cast<T>(0));
    for(int i = 0; i <= degree; ++i)
        summa.coefficients[i] = getCoefAt(i);
    summa.coefficients[0] += value;
    return summa;
}

template<class T>
Polynom<T> Polynom<T>::operator-(const Polynom& pol) const
{
    int resDegree = std::max(degree, pol.degree);
    Polynom<T> difference(resDegree, static_cast<T>(0));
    for(int i = 0; i <= resDegree; ++i)
        difference.coefficients[i] = getCoefAt(i) - pol.getCoefAt(i);
    return difference;
}

template<class T>
Polynom<T> Polynom<T>::operator-(const T& value)
{
    Polynom<T> difference(degree, static_cast<T>(0));
    for(int i = 0; i <= degree; ++i)
        difference.coefficients[i] = getCoefAt(i);
    difference.coefficients[0] -= value;
    return difference;
}

template<class T>
Polynom<T> Polynom<T>::operator*(const Polynom& pol) const
{
    int productDegree = degree + pol.degree;
    Polynom<T> product(productDegree, static_cast<T>(0));
    for(int i = 0; i <= degree; ++i)
        for(int j = 0; j <= pol.degree; ++j)
            product.coefficients[i + j] += this->getCoefAt(i) * pol.getCoefAt(j);
   return product;
}

template<class T>
Polynom<T> Polynom<T>::operator*(const T& value)
{
    Polynom<T> product(degree, static_cast<T>(0));
    for(int i = 0; i <= degree; ++i)
        product.coefficients[i] = getCoefAt(i) * value;
    return product;
}

template<class T>
Polynom<T> Polynom<T>::operator/(const T &value)
{
    Polynom<T> product(degree, static_cast<T>(0));
    for(int i = 0; i <= degree; ++i)
        if(value != static_cast<T>(0))
            product.coefficients[i] = getCoefAt(i) / value;
    return product;
}

template<class T>
Polynom<T>& Polynom<T>::operator+=(const Polynom &pol)
{
    if(pol.getDegree() > degree)
    {
        degree = pol.getDegree();
        coefficients.resize(degree + 1);
    }
    for(int i = 0; i <= degree; ++i)
        coefficients[i] += pol.getCoefAt(i);
    return *this;
}

template<class T>
Polynom<T>& Polynom<T>::operator+=(const T& value)
{
      coefficients[0] += value;
      return *this;
}

template<class T>
Polynom<T>& Polynom<T>::operator-=(const Polynom& pol)
{

    if(pol.getDegree() > degree) coefficients.resize(pol.getDegree());
    for(int i = 0; i <= degree; ++i)
        coefficients[i] -= pol.getCoefAt(i);
    return *this;
}

template<class T>
Polynom<T>& Polynom<T>::operator-=(const T& value)
{
    coefficients[0] -= value;
    return *this;
}

template<class T>
Polynom<T>& Polynom<T>::operator*=(const Polynom& pol)
{
    //Определяем степень результирующего полинома
    int productDegree = degree + pol.getDegree();
    Polynom<T> product(productDegree, static_cast<T>(0));
    for(int i = 0; i <= degree; ++i)
        for(int j = 0; j <= pol.getDegree(); ++j)
            product.coefficients[i + j] += getCoefAt(i) * pol.getCoefAt(j);
    degree = productDegree;
    coefficients.resize(productDegree + 1);
    for(int i = 0; i <= productDegree; ++i)
        coefficients[i] = product.getCoefAt(i);
   return *this;
}

template<class T>
Polynom<T>& Polynom<T>::operator*=(const T& value)
{
    for(int i = 0; i <= degree; ++i)
        coefficients[i] *= value;
    return *this;
}

template<class T>
Polynom<T>& Polynom<T>::operator/=(const T& value)
{
    if(value != static_cast<T>(0))
    for(int i = 0; i <= degree; ++i)
        coefficients[i] /= value;
    return *this;
}

template <class T>
Polynom<T>& Polynom<T>::operator=(const Polynom &pol)
{
   degree = pol.getDegree();
   coefficients.clear();
   coefficients.shrink_to_fit();
   for(int i = 0; i <= degree; ++i)
       coefficients.push_back(pol.coefficients[i]);
   return *this;
}
//Перегрузка операторов икрементации и декрементации
template <class T>
Polynom<T>& operator++(Polynom<T> &pol)
{
    pol.setValueAt(0, pol.getCoefAt(0) + static_cast<T>(1));
    return pol;
}

template <class T>
Polynom<T> operator++(Polynom<T> &pol, int)
{
    Polynom<T> old = pol;
    pol.setValueAt(0, pol.getCoefAt(0) + static_cast<T>(1));
    return old;
}

template <class T>
Polynom<T>& operator--(Polynom<T> &pol)
{
    pol.setValueAt(0, pol.getCoefAt(0) - static_cast<T>(1));
    return pol;
}
template <class T>
Polynom<T> operator--(Polynom<T> &pol, int)
{
    Polynom<T> old = pol;
    pol.setValueAt(0, pol.getCoefAt(0) - static_cast<T>(1));
    return old;
}
//Перегрузка операторов потокового ввода-вывода
template<class T>
istream& operator>>(istream &stream, Polynom<T> &pol)
{
    for(int i = 0; i <= pol.getDegree(); ++i)
        {
            T tmp;
            stream >> tmp;
            pol.setValueAt(i, tmp);
        }
    return stream;
}
template <typename T>
ostream& operator<<(ostream &stream, Polynom<T> &pol)
{
    for(int i = 0; i <= pol.getDegree(); ++i)
       stream << setprecision(15) << pol.getCoefAt(i) << " ";
    stream << endl;
    return stream;
}
template <typename T>
ifstream& operator>>(ifstream &stream, Polynom<T> &pol)
{
    for(int i = 0; i <= pol.getDegree(); ++i)
        {
            T tmp;
            stream >> tmp;
            pol.setValueAt(i, tmp);
        }
    return stream;
}
template <typename T>
ofstream& operator<<(ofstream &stream, Polynom<T> &pol)
{
    for(int i = 0; i <= pol.getDegree(); ++i)
       stream << setprecision(15) << pol.getCoefAt(i) << " ";
    stream << endl;
    return stream;
}
#endif // POLYNOM_H
