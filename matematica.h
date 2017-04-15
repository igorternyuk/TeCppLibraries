#ifndef MATEMATICA_H
#define MATEMATICA_H
#include <math.h>
#include <vector>
#include <ctime>
#include <complex.h>
//#include "matrix.h"
#include "polynom.h"

class CubicSplineInterpolator
{
private:
    struct CubicSpline
    {
        double a, b, c, d, x;
    };
    std::vector<CubicSpline> splineTuple;
    unsigned int nodeNumber;
public:
    CubicSplineInterpolator();
    void interpolate(std::vector<std::pair<double, double>> pointsTable);
    double calculate(double x);
    std::vector<double> allCoef();
};

class LeastSquareMethod
{
private:
    int pointsNumber;
    std::vector<pair<double, double>> xy;
public:
    LeastSquareMethod(vector<pair<double, double>> xy);
    LeastSquareMethod(vector<double> x, vector<double> y);
    void linearApproximation(double &a, double &b, double &error);
    void hyperbolicApproximation(double &a, double &b, double &error);
    void logarithmicApproximation(double &a, double &b, double &error);
    void exponentialApproximation(double &a, double &b, double &error);
    void quadraticApproximation(double &a, double &b, double &c, double &error);
    std::vector<double> polynomialLeastSquareMethod(double polynomDegree, double &error);
};

    static void straightLineCoef(const sf::Vector2f &point1, const sf::Vector2f &point2, float &k, float &b);
    static float distanceFromPointToLine(float k, float b, sf::Vector2f point);
    static bool crossPointOfCircleAndLine(const sf::Vector2f &circle_center, const float &radius,
                const float &k, const float &b, sf::Vector2f &crossPoint1, sf::Vector2f &crossPoint2);
    static bool squareEquation(const float &A, const float &B, const float &C,
                        float &x1, float &x2);

double factorial(double n);
int signum(double value);

//Numericcal integration

double integralTablLeftHandRect(std::vector<std::pair<double, double>> integrandTable);
double integralTablRightHandRect(std::vector<std::pair<double, double>> integrandTable);
double integralTablTrapezoids(std::vector<std::pair<double, double>> integrandTable);
double integralTablSimpson(std::vector<std::pair<double, double>> integrandTable);
double integralTablPolLagr(std::vector<std::pair<double, double>> integrandTable,
                           double a, double b);

double polynomLagrange(std::vector<std::pair<double, double>> pointsTable, double value);
Polynom<double> polynomLagrange(std::vector<std::pair<double, double>> pointsTable);
std::vector<double> polynomLagrangeVector(std::vector<std::pair<double, double>> pointsTable);
std::string polynomLagrangeString(std::vector<std::pair<double, double>> pointsTable, int prec);

//Numericcal differentiation

double tablFuncLeftDer1(std::vector<std::pair<double, double>> pointsTable, double value);
double tablFuncRightDer1(std::vector<std::pair<double, double>> pointsTable, double value);
double tablFuncDer1(std::vector<std::pair<double, double>> pointsTable, double value);
double tablFuncDer2(std::vector<std::pair<double, double>> pointsTable, double value);

//Equations

std::vector<std::complex<double>> solveSquareEquation(double a, double b, double c);
std::vector<std::complex<double>> solveCubicEquation(double a, double b, double c, double d);

#endif // MATEMATICA_H
