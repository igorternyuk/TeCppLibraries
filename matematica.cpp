#include <math.h>
#include <complex>
#include <vector>
#include <stdexcept>
#include "matrix.h"
#include "polynom.h"
#include "matematica.h"
#include <limits>
CubicSplineInterpolator::CubicSplineInterpolator()
{
   //splineTuple.reserve(20);
}
void CubicSplineInterpolator::interpolate(std::vector<std::pair<double, double>> pointsTable)
{
   // Создаем векторы для хранения трех диагоналей трехдиагональной матрицы
   nodeNumber = pointsTable.size();
   //Инициализируем сплайны
   //Вычисляем остальные коэффициенты сплайна
   for(int i = 0; i < nodeNumber; ++i)
   {
       CubicSpline cbspl;
       cbspl.x = pointsTable[i].first;
       cbspl.a = pointsTable[i].second;
       splineTuple.push_back(cbspl);
   }
   // Создаем векторы для коэффициентов прогонки
   // Размер массива с коэффициентами прогонки на единицу меньше количества узлов
   std::vector<double> alpha(nodeNumber - 1), beta(nodeNumber - 1);
   alpha[0] = 0.0;
   beta[0] = 0.0;
   double A, B, C, D, hi, hi1;
   // Вычисляем коэффициенты прогонки (прямой ход)
   for(int i = 1; i < nodeNumber - 1; ++i)
   {
       hi = pointsTable[i].first - pointsTable[i - 1].first;
       hi1 = pointsTable[i + 1].first - pointsTable[i].first;
       A = hi;
       B = 2 * (hi + hi1);
       C = hi1;
       D = 6 * ((pointsTable[i + 1].second - pointsTable[i].second) / hi1 -
               (pointsTable[i].second - pointsTable[i - 1].second) / hi);
       alpha[i] = - C / (A * alpha[i - 1] + B);
       beta[i] = (D - A * beta[i - 1]) /
                 (A * alpha[i - 1] + B);
   }
   splineTuple[0].c = 0.0;
   splineTuple[nodeNumber - 1].c = (D - A* beta[nodeNumber - 2]) /
                                      (A * alpha[nodeNumber - 2] + B);
   //Обратный ход метода прогонки
   for(int i = nodeNumber - 2; i > 0; --i)
   {
       splineTuple[i].c = alpha[i] * splineTuple[i + 1].c + beta[i];
   }
   //Вычисляем остальные коэффициенты сплайна
   for(int i = nodeNumber - 1; i > 0; --i)
   {
       double hi = pointsTable[i].first - pointsTable[i - 1].first;
       splineTuple[i].b = (pointsTable[i].second - pointsTable[i - 1].second) / hi +
                           hi * (2 * splineTuple[i].c + splineTuple[i - 1].c) / 6;
       splineTuple[i].d = (splineTuple[i].c - splineTuple[i - 1].c) / hi;
   }
}

double CubicSplineInterpolator::calculate(double x)
{
    // Если сплайны еще не построены возвращаем NaN
    if(splineTuple.size() == 0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double dx;
    int j = nodeNumber - 1;
    // Если Х меньше абсциссы нулевого узла то пользуемся сплайном нулевого узла
    if(x <= splineTuple[0].x)
    {
        dx = x - splineTuple[1].x;
        j = 1;
    } else if(x >= splineTuple[nodeNumber - 1].x) // Если Х больше абсциссы последнего узла то пользуемся сплайном последнего узла
    {
        dx = x - splineTuple[nodeNumber - 1].x;
        j = nodeNumber - 1;
    } else
    {
        // Бинарный поиск нужного сегмента
        int i = 0;
        while(i + 1 < j)
        {
           int k = i + (j - i) / 2;
           if(x <= splineTuple[k].x)
              j = k;
           else
              i = k;
        }
        dx = x - splineTuple[j].x;
    }
    return splineTuple[j].a + splineTuple[j].b * dx + splineTuple[j].c * pow(dx, 2) / 2 +
           splineTuple[j].d * pow(dx, 3) / 6;
}
std::vector<double> CubicSplineInterpolator::allCoef()
{
    std::vector<double> allCoefs;
    for(int i = 1; i < nodeNumber; ++i)
    {
        allCoefs.push_back(splineTuple[i].a);
        allCoefs.push_back(splineTuple[i].b);
        allCoefs.push_back(splineTuple[i].c);
        allCoefs.push_back(splineTuple[i].d);
    }
    return allCoefs;
}
LeastSquareMethod::LeastSquareMethod(vector<pair<double, double> > xy):
                   pointsNumber(xy.size()), xy(xy)
{}
LeastSquareMethod::LeastSquareMethod(vector<double> x, vector<double> y):
                   pointsNumber(x.size())
{
    for(int i = 0; i < pointsNumber; ++i)
    {
        std::pair<double, double> p(x[i], y[i]);
        xy.push_back(p);
    }
}
void LeastSquareMethod::linearApproximation(double &a, double &b, double &error)
{
    Matrix<double> mtrA(2, 2, 0.0);
    Matrix<double> mtrB(2, 1, 0.0);
    for(int i = 0; i < pointsNumber; ++i)
    {
        mtrA[0][0] += pow(xy[i].first, 2);
        mtrA[0][1] += xy[i].first;
        mtrB[0][0] += xy[i].first * xy[i].second;
        mtrB[1][0] += xy[i].second;
    }
    mtrA[1][0] = mtrA.getValueAt(0, 1);
    mtrA[1][1] = (double)pointsNumber;
    if(mtrA.determinant() != 0)
    {
        Matrix<double> mtrRes = mtrA.inverse() * mtrB;
        a = mtrRes.getValueAt(0, 0);
        b = mtrRes.getValueAt(1, 0);
    }else
    {
        a = 0;
        b = 0;
    }
    for(int i = 0; i < pointsNumber; ++i)
    {
        error += pow(xy[i].second - a * xy[i].first - b, 2);
        error /= pointsNumber;
        error = sqrt(error);
    }
}
void LeastSquareMethod::hyperbolicApproximation(double &a, double &b, double &error)
{
    Matrix<double> mtrA(2, 2, 0.0);
    Matrix<double> mtrB(2, 1, 0.0);
    for(int i = 0; i < pointsNumber; ++i)
    {
        mtrA[0][0] += pow(xy[i].first, -2);
        mtrA[0][1] += pow(xy[i].first, -1);
        mtrB[0][0] += xy[i].second / xy[i].first;
        mtrB[1][0] += xy[i].second;
    }
    mtrA[1][0] = mtrA.getValueAt(0, 1);
    mtrB[1][1] = (double)pointsNumber;
    if(mtrA.determinant() != 0)
    {
        Matrix<double> mtrRes = mtrA.inverse() * mtrB;
        a = mtrRes.getValueAt(0, 0);
        b = mtrRes.getValueAt(1, 0);
    }
    else
    {
        a = 0;
        b = 0;
    }
    for(int i = 0; i < pointsNumber; ++i)
    {
        error += pow(xy[i].second - a / xy[i].first - b, 2);
        error /= pointsNumber;
        error = sqrt(error);
    }
}
void LeastSquareMethod::logarithmicApproximation(double &a, double &b, double &error)
{
    Matrix<double> mtrA(2, 2, 0.0);
    Matrix<double> mtrB(2, 1, 0.0);
    for(int i = 0; i < pointsNumber; ++i)
    {
        if(xy[i].first < 0) continue;
        mtrA[0][0] += pow(log(xy[i].first), 2);
        mtrA[0][1] += log(xy[i].first);
        mtrB[0][0] += xy[i].second * log(xy[i].first);
        mtrB[1][0] += xy[i].second;
    }
    mtrA[1][0] = mtrA.getValueAt(0, 1);
    mtrB[1][1] = (double)pointsNumber;
    if(mtrA.determinant() != 0)
    {
        Matrix<double> mtrRes = mtrA.inverse() * mtrB;
        a = mtrRes.getValueAt(0, 0);
        b = mtrRes.getValueAt(1, 0);
    }
    else
    {
        a = 0;
        b = 0;
    }
    for(int i = 0; i < pointsNumber; ++i)
    {
        if(xy[i].first < 0) continue;
        error += pow(xy[i].second - a * log(xy[i].first) - b, 2);
        error /= pointsNumber;
        error = sqrt(error);
    }
}
void LeastSquareMethod::exponentialApproximation(double &a, double &b, double &error)
{
    Matrix<double> mtrA(2, 2, 0.0);
    Matrix<double> mtrB(2, 1, 0.0);
    for(int i = 0; i < pointsNumber; ++i)
    {
        if(xy[i].first < 0) continue;
        mtrA[0][0] += pow(xy[i].first, 2);
        mtrA[0][1] += xy[i].first;
        mtrB[0][0] += xy[i].first * log(xy[i].second);
        mtrB[1][0] += xy[i].second;
    }
    mtrA[1][0] = mtrA.getValueAt(0, 1);
    mtrB[1][1] = (double)pointsNumber;
    if(mtrA.determinant() != 0)
    {
        Matrix<double> mtrRes = mtrA.inverse() * mtrB;
        a = exp(mtrRes.getValueAt(0, 0));
        b = mtrRes.getValueAt(1, 0);
    }
    else
    {
        a = 0;
        b = 0;
    }
    for(int i = 0; i < pointsNumber; ++i)
    {
        error += pow(xy[i].second - a * exp(b * xy[i].first), 2);
        error /= pointsNumber;
        error = sqrt(error);
    }
}
void LeastSquareMethod::quadraticApproximation(double &a, double &b, double &c, double &error)
{
    int pointsNumber = xy.size();
    Matrix<double> mtrA(3, 3, 0.0);
    Matrix<double> mtrB(3, 1, 0.0);
    for(int i = 0; i < pointsNumber; ++i)
    {
        mtrA[0][0] += pow(xy[i].first, 4);
        mtrA[0][1] += pow(xy[i].first, 3);
        mtrA[0][2] += pow(xy[i].first, 2);
        mtrA[1][2] += xy[i].first;
        mtrB[0][0] += pow(xy[i].first,2) * xy[i].second;
        mtrB[1][0] += xy[i].first * xy[i].second;
        mtrB[2][0] += xy[i].second;
    }
    mtrA[1][0] = mtrA.getValueAt(0, 1);
    mtrA[1][1] = mtrA.getValueAt(0, 2);
    mtrA[2][0] = mtrA.getValueAt(0, 2);
    mtrA[2][1] = mtrA.getValueAt(1, 2);
    mtrA[2][2] = (double)pointsNumber;
    if(mtrA.determinant() != 0)
    {
        Matrix<double> mtrRes = mtrA.inverse() * mtrB;
        a = mtrRes.getValueAt(0, 0);
        b = mtrRes.getValueAt(1, 0);
        c = mtrRes.getValueAt(2, 0);
    }else
    {
        a = 0;
        b = 0;
        c = 0;
    }
    for(int i = 0; i < pointsNumber; ++i)
    {
        error += pow(xy[i].second - a * pow(xy[i].first, 2) - b * xy[i].first - c, 2);
        error /= pointsNumber;
        error = sqrt(error);
    }
}
std::vector<double> LeastSquareMethod::polynomialLeastSquareMethod(double polynomDegree, double &error)
{
    if(polynomDegree > pointsNumber - 1) polynomDegree = pointsNumber - 1;
    // Определяем размер матрицы СЛАУ для нахождения коэффициентов аппроксимации
    int N = polynomDegree + 1;
    Matrix<double> mtrA(N, N, 0.0);
    Matrix<double> mtrB(N, 1, 0.0);
    for(int i = 0; i < N; ++i)
    {
        for(int j = 0; j < N; ++j)
        {
            for(int k = 0; k < pointsNumber; ++k)
                mtrA[i][j] += pow(xy[k].first, i + j);
        }
        for(int k = 0; k < pointsNumber; ++k)
            mtrB[i][0] += xy[k].second * pow(xy[k].first, i);
    }
    std::vector<double> res;
    if(mtrA.determinant() != 0)
    {
        Matrix<double> mtrRes = mtrA.inverse() * mtrB;
        for(int i = 0; i < N; ++i)            
            res.push_back(mtrRes.getValueAt(i, 0));
        return res;

    }else
    {
        for(int i = 0; i < N; ++i)
            res.push_back(0.0);
        return res;
    }
    for(int i = 0; i < pointsNumber; ++i)
    {
        double polynomVal = 0.0;
        for(int k = 0; k < N; ++k)
        {
            polynomVal += res[k] * pow(xy[i].first, k); 
        }
        error += pow(xy[i].second - polynomVal, 2);
        error /= pointsNumber;
        error = sqrt(error);
    }
}
double factorial(double n)
{
    if (floor(n) == 0 || floor(n) == 1) return 1;
    else return n * factorial(n - 1);
}
int signum(double value)
{
    return value >= 0 ? 1: -1;
}
double polynomLagrange(std::vector<std::pair<double, double> > pointsTable, double value)
{
    double summa = 0, num, den, omega;                       // L = w1 * y1 + w2 * y2 + ... + wn * yn
    for(unsigned int i = 0; i < pointsTable.size(); ++i)
    {
        omega = 1;
        num = 1;
        den = 1;
        for(unsigned int j = 0; j < pointsTable.size(); ++j)
        {
            if (i != j)
            {
              num *= (value - pointsTable[j].first);
              den *= pointsTable[i].first - pointsTable[j].first;
            }
        }
        omega = num / den;
        summa += omega * pointsTable[i].second;
    }
    return summa;
}
Polynom<double> polynomLagrange(std::vector<std::pair<double, double> > pointsTable)
{
    Polynom<double> polynomSumma(0, 0);
    double omega = 1;
    double den = 1;
    for(unsigned int i = 0; i < pointsTable.size(); ++i)
    {
        omega = 1;
        den = 1;
        Polynom<double> polynomProduct(0, 1);
        for(unsigned int j = 0; j < pointsTable.size(); ++j)
        {
            if (i != j)
            {
              std::vector<double> polCoefs;
              polCoefs.push_back(-pointsTable[j].first);
              polCoefs.push_back(1);
              Polynom<double> pol(polCoefs);
              polynomProduct *= pol;
              den *= pointsTable[i].first - pointsTable[j].first;
            }
        }
        omega = pointsTable[i].second / den;
        polynomProduct *= omega;
        polynomSumma += polynomProduct;
    }
    return polynomSumma;
}
std::vector<double> polynomLagrangeVector(std::vector<std::pair<double, double> > pointsTable)
{
    Polynom<double> result = polynomLagrange(pointsTable);
    std::vector<double> resVct;
    for(int i = 0; i <= result.getDegree(); ++i)
        resVct.push_back(result.getCoefAt(i));
    return resVct;
}
std::string polynomLagrangeString(std::vector<std::pair<double, double> > pointsTable, int prec)
{
    Polynom<double> result = polynomLagrange(pointsTable);
    return result.convertToString(prec);
}

double integralTablLeftHandRect(std::vector<std::pair<double, double>> integrandTable)
{
    if(!integrandTable.empty())
    {
        int n = integrandTable.size() - 1;
        double currentStep, Integral = 0;
        for(int i = 0; i < n; ++i)
        {
           currentStep = integrandTable[i + 1].first - integrandTable[i].first;
           Integral += integrandTable[i].second * currentStep;
        }
        return Integral;
    }else return 0.0;
}
double integralTablRightHandRect(std::vector<std::pair<double, double>> integrandTable)
{
    if(!integrandTable.empty())
    {
        int n = integrandTable.size() - 1;
        double currentStep, Integral = 0;
        for(int i = 1; i <= n; ++i)
        {
           currentStep = integrandTable[i].first - integrandTable[i - 1].first;
           Integral += integrandTable[i].second * currentStep;
        }
        return Integral;
    }else return 0.0;
}
double integralTablTrapezoids(std::vector<std::pair<double, double>> integrandTable)
{
    if(!integrandTable.empty())
    {
        int n = integrandTable.size() - 1;
        double currentStep, Integral = 0;
        for(int i = 1; i <= n; ++i)
        {
           currentStep = integrandTable[i].first - integrandTable[i - 1].first;
           Integral += (integrandTable[i].second + integrandTable[i - 1].second) * currentStep;
        }
        Integral /= 2;
        return Integral;
    }else return 0.0;
}
double integralTablSimpson(std::vector<std::pair<double, double> > integrandTable)
{
    if(!integrandTable.empty())
    {
        double h1, h2; // steps
        double Integral = 0;
        int n = integrandTable.size() - 1;
        if(n % 2 != 0) throw std::runtime_error("The number of intervals must be even!");
        for(int i = 0; i <= n - 2; i += 2)
            {
                h1 = integrandTable[i + 1].first - integrandTable[i].first;
                h2 = integrandTable[i + 2].first - integrandTable[i + 1].first;
                Integral += ((h1 + h2) / (6 * h1 * h2)) *
                (h2 * (2 * h1 - h2) * integrandTable[i].second + pow(h1 + h2, 2) * integrandTable[i + 1].second +
                h1 * (2 * h2 - h1) * integrandTable[i + 2].second);
            }
        return Integral;
    }else return 0.0;

}
double integralTablPolLagr(std::vector<std::pair<double, double> > integrandTable,
                           double a, double b)
{
    if(!integrandTable.empty())
    {
        int n = integrandTable.size() - 1;
        Matrix<double> mtrX(n, n, 0.0);
        Matrix<double> mtrI(n, 1, 0.0);
        //Заполняем матрицы
        for(int i = 0; i < mtrX.getRowCount(); ++i)
        {
            mtrI.setValueAt(i, 0, (pow(b, i + 1) - pow(a, i + 1)) / (i + 1));
            for(int j = 0; j < mtrX.getColumnCount(); ++j)
            {
                mtrX.setValueAt(i, j, pow(integrandTable[j].first, i));
            }
        }
        Matrix<double> mtrInv = mtrX.inverse();
        Matrix<double> mtrA = mtrInv * mtrI;
        double Integral = 0;
        for(int i = 0; i <= n; ++i)
        {
           Integral += mtrA.getValueAt(i,0) * integrandTable[i].second;
        }
        return Integral;
    }else return 0.0;
}
std::vector<std::complex<double>> solveCubicEquation(double a, double b, double c, double d)
{
    double p = c / a - pow(b, 2) / (3 * pow(a, 2));
    double q = (2 * pow(b, 3) - 9 * a * b * c + 27 * a * a * d) / pow(3 * a, 3);
    std::complex<double> alpha, beta;
    std::complex<double> Q = pow(p / 3, 3) + pow(q / 2, 2);
    if(Q.real() >= 0)
    {
        alpha = signum(-0.5 * q + sqrt(Q.real())) * pow(fabs(-0.5 * q + sqrt(Q.real())), 1 / 3.0);
        beta = signum(-0.5 * q - sqrt(Q.real())) * pow(fabs(-0.5 * q - sqrt(Q.real())), 1 / 3.0);
    }else
    {
        std::complex<double> Zalpha = -0.5 * q + sqrt(Q);
        std::complex<double> Zbeta = -0.5 * q - sqrt(Q);
        alpha = pow(Zalpha, 1 / 3.0);
        beta = pow(Zbeta , 1 / 3.0);
    }
    std::complex<double> y1 = alpha + beta;
    std::complex<double> unit(-1,0);
    std::complex<double> j = sqrt(unit);
    std::complex<double> x1 = y1 - b / (3 * a);
    std::complex<double> x2 = -0.5 * (alpha + beta) + 0.5 * (alpha - beta) * sqrt(3.0) * j - b / (3 * a);
    std::complex<double> x3 = -0.5 * (alpha + beta) - 0.5 * (alpha - beta) * sqrt(3.0) * j - b / (3 * a);
    std::vector<std::complex<double>> roots;
    roots.push_back(x1);
    roots.push_back(x2);
    roots.push_back(x3);
    return roots;
}
std::vector<std::complex<double>> solveSquareEquation(double a, double b, double c)
{
    std::complex<double> discr(pow(b, 2) - 4 * a * c, 0.0);
    std::complex<double> x1 = (-b + sqrt(discr)) / (2 * a);
    std::complex<double> x2 = (-b - sqrt(discr)) / (2 * a);
    std::vector<std::complex<double>> roots;
    roots.push_back(x1);
    roots.push_back(x2);
    return roots;
}
double tablFuncLeftDer1(std::vector<std::pair<double, double> > pointsTable, double value)
{
    const double deltaArgValue = 0.001;
    double func = polynomLagrange(pointsTable, value);
    double funcMinusDelta = polynomLagrange(pointsTable, value - deltaArgValue);
    double funcMinusTwoDeltas = polynomLagrange(pointsTable, value - 2 * deltaArgValue);
    return (3 * func - 4 * funcMinusDelta + funcMinusTwoDeltas) / (2 * deltaArgValue);
}
double tablFuncRightDer1(std::vector<std::pair<double, double> > pointsTable, double value)
{
   const double deltaArgValue = 0.001;
   double func = polynomLagrange(pointsTable, value);
   double funcPlusDelta = polynomLagrange(pointsTable, value + deltaArgValue);
   double funcPlusTwoDeltas = polynomLagrange(pointsTable, value + 2 * deltaArgValue);
   return (-3 * func + 4 * funcPlusDelta - funcPlusTwoDeltas) / (2 * deltaArgValue);
}
double tablFuncDer1(std::vector<std::pair<double, double> > pointsTable, double value)
{
    const double deltaArgValue = 0.001;
    double funcMinusDelta = polynomLagrange(pointsTable, value - deltaArgValue);
    double funcPlusDelta = polynomLagrange(pointsTable, value + deltaArgValue);
    return (funcPlusDelta - funcMinusDelta) / (2 * deltaArgValue);
}
double tablFuncDer2(std::vector<std::pair<double, double> > pointsTable, double value)
{
    const double deltaArgValue = 0.001;
    double funcPlusDelta = polynomLagrange(pointsTable, value + deltaArgValue);
    double func= polynomLagrange(pointsTable, value);
    double funcMinusDelta = polynomLagrange(pointsTable, value - deltaArgValue);
    return (funcPlusDelta - 2 * func + funcMinusDelta) / pow(deltaArgValue, 2);
}





