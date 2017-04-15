#ifndef PARSER_H
#define PARSER_H
#include <string>
#include <vector>
#include <cmath>
#include <ctype.h>
#include <cstring>
#include <QString>

struct Expression {
    Expression(std::string token) : token(token) {}
    Expression(std::string token, Expression a) : token(token), args{ a } {}
    Expression(std::string token, Expression a, Expression b) : token(token), args{ a, b } {}
    std::string token;
    std::vector<Expression> args;
};

class Parser {
private:
    const char* input;
    std::string angleUnit;
    std::string inputWithVariables;
    std::vector<std::pair<char, double>> vctVariables;
    std::vector<std::string> tokens =
    { "+", "-","e", "**", "*", "/","(", ")",  "<=", ">=",
      "<", ">", "==", "!=", "&", "|","^","!", "mod",
      "abs", "factorial","sign", "inv", "cbrt","sqrt", "sqr", "cube", "gradtorad",
      "radtograd", "_exp","ln", "log2", "log8","log10","log16",
      "sin","cos","tg","ctg","secans","csecans",
      "arcsin","arccos","arctg","arcctg","arcsecans", "arccsecans",
      "sh","ch", "th", "cth","sech","csech",
      "arcsh","arcch", "arcth","arccth","arcsech","arccsech"};
public:
    Parser(const std::string &inputString, std::vector<std::pair<char,double>> vars);
    Parser(const std::string &inputString, std::vector<std::pair<char,double>> vars, std::string angleUnit);
    double calculateExpression();
private:
    std::string parseToken();
    Expression parseUnaryExpression();
    Expression parseBinaryExpression(int minPriority);
    Expression parse();
    double evaluateExpression(const Expression &e);
};
double func(std::string expr, char argName, double argValue);
double funcLeftDer1(std::string expr, char argName, double argValue);
double funcRightDer1(std::string expr, char argName, double argValue);
double funcDer1(std::string expr, char argName, double argValue);
double funcDer2(std::string expr, char argName, double argValue);
double multiVariableFunc(std::string expr, std::vector<std::pair<char, double>> varList);
double funcPartialDer1(std::string expr, std::vector<std::pair<char, double>> variables, char argName);
double integralLeftHandRect(std::string expr, char argName, double leftLim, double rightLim,
                            double eps, int &iterrNumber, std::string &iterProcess);
double integralRightHandRect(std::string expr, char argName, double leftLim, double rightLim,
                             double eps, int &iterrNumber, std::string &iterProcess);
double integralMidpoint(std::string expr, char argName, double leftLim, double rightLim,
                        double eps, int &iterrNumber, std::string &iterProcess);
double integralTrapezoids(std::string expr, char argName, double leftLim, double rightLim,
                          double eps, int &iterrNumber, std::string &iterProcess);
double integralSimpson(std::string expr, char argName, double leftLim, double rightLim,
                       double eps, int &iterrNumber, std::string &iterProcess);
double solveEquationChordsMethod(std::string expr, char argName, double a, double b, double eps, int &iterrNumber, std::string &iterProcess);
double solveEquationTangentsMethod(std::string expr, char argName, double a, double b, double eps, int &iterrNumber, std::string &iterProcess);
double solveEquationBisectionMethod(std::string expr, char argName, double a, double b, double eps, int &iterrNumber, std::string &iterProcess);
double solveEquationCombinedMethod(std::string expr, char argName, double a, double b,
                                     double eps, int &iterrNumber, std::string &iterProcess);
double solveEquationIterationMethod(std::string expr, char argName, double a, double b,
                                     double eps, double lambda, int &iterrNumber, std::string &iterProcess);
std::vector<double> solveSystemOfNonLinearEquations(std::vector<std::string> equations, std::vector<std::pair<char, double>> variables,
                                     double eps, int & iterNumber, double &err, std::string &iterProcess);
#endif // PARSER_H
