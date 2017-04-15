#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <ctype.h>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include "parser.h"
#include <QDebug>
#include "matrix.h"
#include "matematica.h"

Parser::Parser(const std::string &inputString, std::vector<std::pair<char,double>> vars):
       Parser(inputString, vars, "radian")
{}
Parser::Parser(const std::string &inputString, std::vector<std::pair<char,double>> vars, std::string angleUnit)
{
    this->angleUnit = angleUnit;
    inputWithVariables = inputString;
    int leftParentethisNumber = 0;
    int rightParentethisNumber = 0;
    for(auto &s : inputString)
    {
        if (s == '(') leftParentethisNumber++;
        if (s == ')') rightParentethisNumber++;
    }
    if(leftParentethisNumber != rightParentethisNumber)
        throw std::runtime_error("Unbalanced parentethises") ;
    //Добавляем в вектор с парами = имя переменной - значение 3 константы
    //число ПИ, число Непера и константу золотого сечения
    vctVariables = vars;
    std::pair<char, double> Pi('P', M_PI);
    vctVariables.push_back(Pi);
    std::pair<char, double> E('E', M_E);
    vctVariables.push_back(E);
    std::pair<char, double> G('G', (1 + sqrt(5)) / 2);
    vctVariables.push_back(G);

    // Во внешнем циклке перебираем все элементы вектора переменных
 //   for(std::pair<char, double> p : vctVariables) // В этом амперсанде может быть подвох!
    for(auto p : vctVariables) // В этом амперсанде может быть подвох!
    {
        // Во внутреннем циклке сравниваем каждый символ входной строки с текущим элементом в векторе переменных

        for(unsigned int i = 0; i < inputWithVariables.length(); ++i)
           {
            if(inputWithVariables.at(i) == p.first)
            {
                //Скопировать часть строки с начала до х
                std::string buf1 = inputWithVariables.substr(0,i);
                //Скопировать часть строки с х + 1 до конца
                std::string buf2 = inputWithVariables.substr(i + 1, inputWithVariables.length() - i);
                //Конвертировать значение переменной в строку обернув в скобки
                std::ostringstream sstream;
                //sstream << "(" << variable.second << ")";
                sstream << p.second;
                std::string res = sstream.str();
                // добавляем к первой части значение икса
                buf1 += res;
                //Конкатинируем строки
                inputWithVariables = buf1 + buf2;
            }
          }
    }
    input = inputWithVariables.c_str();
    qDebug() << input;
}
std::string Parser::parseToken()
{
    //std::cout << "input= " << input << std::endl;
    while (isspace(*input)) ++input;
    // Проверка является токен числом
    if (isdigit(*input))
    {
        std::string number;
        while (isdigit(*input) || *input == '.') number.push_back(*input++);
        return number;
    }
    // Проверка является токен оператором или функцией
    for (auto& t : tokens)
    {
        if (strncmp(input, t.c_str(), t.size()) == 0)
        {
            input += t.size();
            return t;
        }
    }
    // Если совпадений нет возвращаем пустую строку
    return "";
}
Expression Parser::parseUnaryExpression() {
    auto token = parseToken();
    if (token.empty()) throw std::runtime_error("Invalid input");
    if (token == "(") {
        auto result = parse();
        if (parseToken() != ")") throw std::runtime_error("Expected ')'");
        return result;
    }
    if (isdigit(token[0]))
        return Expression(token);
    return Expression(token, parseUnaryExpression());
}
int getPriority(const std::string& binaryOperation) {
    if (binaryOperation == "+" || binaryOperation == "-") return 1;
    if (binaryOperation == "*" || binaryOperation == "/" || binaryOperation == "mod" ||
        binaryOperation == "e") return 2;
    if (binaryOperation == "&" || binaryOperation == "|" || binaryOperation == "^" ) return 3;
    if (binaryOperation == "<=" || binaryOperation == ">=" || binaryOperation == "<" ||
        binaryOperation == ">" || binaryOperation == "==" || binaryOperation == "!=") return 4;
    if (binaryOperation == "**") return 5;
    return 0;
}
Expression Parser::parseBinaryExpression(int minPriority) {
    auto leftExpression = parseUnaryExpression();
    for (;;) {
        auto op = parseToken();
        auto priority = getPriority(op);
        if (priority <= minPriority) {
            input -= op.size();
            return leftExpression;
        }
        auto rightExpression = parseBinaryExpression(priority);
        leftExpression = Expression(op, leftExpression, rightExpression);
    }
}
Expression Parser::parse() {
    return parseBinaryExpression(0);
}

double Parser::evaluateExpression(const Expression &e)
{
    switch (e.args.size()) {
    case 2: {
        auto a = evaluateExpression(e.args[0]);
        auto b = evaluateExpression(e.args[1]);
        if (e.token == "+") return a + b;
        if (e.token == "-") return a - b;
        if (e.token == "e") return a * pow(10, b);
        if (e.token == "*") return a * b;
        if (e.token == "/")
        {
            if(b != 0)
            return a / b;
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "**")
        {
            return pow(a, b);
        }
        if (e.token == "mod")
        {
            if(b != 0)
                 return (int)a % (int)b;
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "<=") return a <= b ? 1 : 0;
        if (e.token == ">=") return a >= b ? 1 : 0;
        if (e.token == "<") return a < b ? 1 : 0;
        if (e.token == ">") return a > b ? 1 : 0;
        if (e.token == "==") return a == b ? 1 : 0;
        if (e.token == "!=") return a != b ? 1 : 0;
        if (e.token == "&") return (a != 0 && b != 0) ? 1 : 0;
        if (e.token == "|") return (a != 0 || b != 0) ? 1 : 0;
        if (e.token == "^") return ((a != 0 && b == 0) || (a == 0 && b != 0)) ? 1 : 0;
        throw std::runtime_error("Unknown binary operator");
    }
    case 1: {
        auto a = evaluateExpression(e.args[0]);
        if (e.token == "+") return +a;
        if (e.token == "-") return -a;
        if (e.token == "!") return (a != 0) ? 0 : 1;
        if (e.token == "factorial")
        {
            double wholePart = abs(floor(a));
            return factorial(wholePart);
        }
        if (e.token == "inv") return (a != 0) ? pow(a, -1) : 0.0;
        if (e.token == "sign") return (a >= 0) ? 1 : -1;
        if (e.token == "abs") return fabs(a);
                if (e.token == "cbrt") return pow(a,pow(3,-1));
        if (e.token == "sqr") return pow(a,2);
        if (e.token == "cube") return pow(a,3);
        if (e.token == "gradtorad") return M_PI * a / 180;
        if (e.token == "radtograd") return 180 * a / M_PI;
        if (e.token == "_exp") return exp(a);
        if (e.token == "ln")
        {
            if(a >= 0)
                 return log(a);
            else throw  std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "log8")
        {
            if(a >= 0)
                 return log10(a)/log10(8);
            else throw  std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "log10")
        {
            if(a >= 0)
                 return log10(a);
            else throw  std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "log16")
        {
            if(a >= 0)
                 return log10(a)/log10(16);
            else throw  std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "log2")
        {
            if(a >= 0)
                 return log2(a);
            else throw  std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "arcsh") return asinh(a);
        if (e.token == "arcch")
        {
            if(a >= 1)
             return acosh(a);
            else throw std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "arcth")
        {
            if(fabs(a) < 1)
                 return atanh(a);
            else throw std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "arccth")
        {
            if(a == 0)
                throw std::runtime_error("Division by zero!");
            if(fabs(a) > 1)
                return atanh(pow(a,-1));
            else throw std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "arcsech")
        {
            if(a != 0)
                 return asinh(pow(a,-1));
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "arccsech")
        {
            if(a == 0)
               throw std::runtime_error("Division by zero!");
            if(a <= 1)
             return acosh(pow(a,-1));
            else throw std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "sqrt")
        {
            if(a >= 0)
                 return sqrt(a);
            else throw std::runtime_error("Negative radicand!");
        }
        if (e.token == "sin")
        {
            if(angleUnit == "radian") return sin(a);
            if(angleUnit == "gradus") return sin(M_PI * a / 180);
            if(angleUnit == "grad") return sin(M_PI * a / 200);
            return sin(a);
        }
        if (e.token == "cos")
        {
            if(angleUnit == "radian") return cos(a);
            if(angleUnit == "gradus") return cos(M_PI * a / 180);
            if(angleUnit == "grad") return cos(M_PI * a / 200);
            return cos(a);
        }
        if (e.token == "tg")
        {
            if(angleUnit == "radian") return tan(a);
            if(angleUnit == "gradus") return tan(M_PI * a / 180);
            if(angleUnit == "grad") return tan(M_PI * a / 200);
            return tan(a);
        }
        if (e.token == "ctg")
        {
            if(tan(a) != 0)
            {
                if(angleUnit == "radian") return pow(tan(a), -1);
                if(angleUnit == "gradus") return pow(tan(M_PI * a / 180),-1);
                if(angleUnit == "grad") return pow(tan(M_PI * a / 200), -1);
                return pow(tan(a), -1);
            }
            else throw std::runtime_error("Division by zero!");
        }

        if (e.token == "csecans")
        {
            if(cos(a) != 0)
            {
                if(angleUnit == "radian") return pow(cos(a), -1);
                if(angleUnit == "gradus") return pow(cos(M_PI * a / 180),-1);
                if(angleUnit == "grad") return pow(cos(M_PI * a / 200), -1);
                return pow(cos(a), -1);
            }
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "secans")
        {
            if(sin(a) != 0)
            {
                if(angleUnit == "radian") return pow(sin(a), -1);
                if(angleUnit == "gradus") return pow(sin(M_PI * a / 180),-1);
                if(angleUnit == "grad") return pow(sin(M_PI * a / 200), -1);
            }
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "arcsin")
        {
            if(fabs(a) <= 1)
                   return asin(a);
            else throw std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "arccos")
        {
            if(fabs(a) <= 1)
                    return acos(a);
            else throw std::runtime_error("argNameument is out of the function domen");
        }
        if (e.token == "arctg") return atan(a);
        if (e.token == "arcctg")
        {
            if(a != 0)
                    return atan(1 / a);
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "arcsecans")
        {
            if(a != 0)
                    return asin(1 / a);
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "arccsecans")
        {
            if(a != 0)
                     return acos(1 / a);
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "sh") return sinh(a);
        if (e.token == "ch") return cosh(a);
        if (e.token == "th") return tanh(a);
        if (e.token == "cth")
        {
            if(tanh(a) != 0)
                      return pow(tanh(a),-1);
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "sech")
        {
            if(sinh(a) != 0)
                       return pow(sinh(a), -1);
            else throw std::runtime_error("Division by zero!");
        }
        if (e.token == "csech") return 1 / cosh(a);
        {
            if(cosh(a) != 0)
                return pow(cosh(a), -1);
            else throw std::runtime_error("Division by zero!");
        }
        throw std::runtime_error("Unknown unary operator");
    }
    case 0:
    {
        QString resString = QString::fromStdString(e.token.c_str());
        return resString.toDouble();
    }
    }
    throw std::runtime_error("Unknown expression type");
}
double Parser::calculateExpression()
{
    auto result = evaluateExpression(parse());
    return result;
}
double func(std::string expr, char argName, double argValue)
{
    std::pair<char, double> var(argName, argValue);
    std::vector<std::pair<char, double>> vctVars;
    vctVars.push_back(var);
    Parser prs(expr, vctVars);
    return prs.calculateExpression();
}
double funcLeftDer1(std::string expr, char argName, double argValue)
{
    const double deltaArgValue = 0.001;
    double funcX= func(expr, argName, argValue);
    double funcMinusTwoDeltas = func(expr, argName, argValue - 2 * deltaArgValue);
    double funcMinusDelta = func(expr, argName, argValue - deltaArgValue);
    return (3 * funcX - 4 * funcMinusDelta + funcMinusTwoDeltas) / (2 * deltaArgValue);
}
double funcRightDer1(std::string expr, char argName, double argValue)
{
    const double deltaArgValue = 0.001;
    double funcPlusDelta = func(expr, argName, argValue + deltaArgValue);
    double funcX= func(expr, argName, argValue);
    double funcPlusTwoDeltas = func(expr, argName, argValue + 2 * deltaArgValue);
    return (-3 * funcX + 4 * funcPlusDelta - funcPlusTwoDeltas) / (2 * deltaArgValue);
}
double funcDer1(std::string expr, char argName, double argValue)
{
    const double deltaArgValue = 0.001;
    double funcPlusDelta = func(expr, argName, argValue + deltaArgValue);
    double funcMinusDelta = func(expr, argName, argValue - deltaArgValue);
    double res = (funcPlusDelta - funcMinusDelta) / (2 * deltaArgValue);
    return res;
}
double funcDer2(std::string expr, char argName, double argValue)
{
    const double deltaArgValue = 0.001;
    double funcPlusDelta = func(expr, argName, argValue + deltaArgValue);
    double funcX= func(expr, argName, argValue);
    double funcMinusDelta = func(expr, argName, argValue - deltaArgValue);
    return (funcPlusDelta - 2 * funcX + funcMinusDelta) / pow(deltaArgValue, 2);
}
double multiVariableFunc(std::string expr, std::vector<std::pair<char, double>> varList)
{
    Parser prs(expr, varList);
    return prs.calculateExpression();
}
double funcPartialDer1(std::string expr, std::vector<std::pair<char, double>> variables,char argName)
{
    double delta = 0.001;

    for (auto &p : variables)
    {
        if(p.first == argName)
        {
            p.second +=delta;
            break;
        }
    }
    Parser prs_plus(expr, variables);
    double funcPlusDeltaValue =  prs_plus.calculateExpression();
    for (auto &p : variables)
    {
        if(p.first == argName)
        {
            p.second -= 2 * delta;
            break;
        }
    }
    Parser prs_minus(expr, variables);
    double funcMinusDeltaValue =  prs_minus.calculateExpression();
    return (funcPlusDeltaValue - funcMinusDeltaValue) / 2 / delta;
}

double integralLeftHandRect(std::string expr, char argName, double leftLim, double rightLim,
                                     double eps, int &iterrNumber, std::string &iterProcess)
{
    double a = leftLim, b = rightLim;
    int n = 6;
    double h = (b - a) / n;
    double Inext = 0, Iprev = 0;
    iterrNumber = 0;
    std::stringstream stream;
    stream << "Выполняется численное интегрирование методом левых прямоугольников\n";
    const int MAX_ITERRATIONS_NUMBER = 5000;
    bool flag = false;
    do
    {
        stream << "------------------------------------------------------\n" ;
        stream << "Выполняется итеррация номер : " << ++iterrNumber << "\n";
        Iprev = Inext;
        Inext = 0;
        n *= 2;
        h = (b - a) / n;

        for(int i = 0; i < n; ++i)
        {
            Inext += func(expr, argName, a + i * h);
        }
        Inext *= h;
        stream << "Количество интервалов разбиения на данной итеррации равно:\n n = " << n << "\n";
        stream << "Ширина i-го интервала на данной итеррации равна:\n hi = " << h << "\n";
        stream << "Значение интеграла на данной  итеррации равно:\n Ii = " << Inext << "\n";
        stream << "Погрешность данной  итеррации равна:\n ΔIi = Ii - Ii-1 = " << Inext << "\n";
        if(fabs(Inext - Iprev) < fabs(eps))
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ ДОСТИГНУТА!\n";
            stream << "ЗНАЧЕНИЕ ИНТЕГРАЛА РАВНО: I = " << Inext << "\n";
            flag = true;
        }
        else if (iterrNumber > MAX_ITERRATIONS_NUMBER )
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ НЕ ДОСТИГНУТА!\n";
            break;
        }
    }while(!flag);
    iterProcess = stream.str();
    return Inext;
}
double integralRightHandRect(std::string expr, char argName, double leftLim, double rightLim,
                                      double eps, int &iterrNumber, std::string &iterProcess)
{
    double a = leftLim, b = rightLim;
    int n = 6;
    double h = (b - a) / n;
    double Inext = 0, Iprev = 0;
    iterrNumber = 0;
    const int MAX_ITERRATIONS_NUMBER = 5000;
    bool flag = false;
    std::stringstream stream;
    stream << "Выполняется численное интегрирование методом правых прямоугольников\n";
    do
    {
        stream << "Выполняется итеррация номер : \n" << ++iterrNumber << "\n";
        Iprev = Inext;
        Inext = 0;
        n *= 2;
        h = (b - a) / n;
        for(int i = 1; i <= n; ++i)
        {
            Inext += func(expr, argName, a + i * h);
        }
        Inext *= h;
        stream << "Количество интервалов разбиения на данной итеррации равно:\n n = " << n << "\n";
        stream << "Ширина i-го интервала на данной итеррации равна:\n hi = " << h << "\n";
        stream << "Значение интеграла на данной  итеррации равно:\n Ii = " << Inext << "\n";
        stream << "Погрешность данной  итеррации равна:\n ΔIi = Ii - Ii-1 = " << Inext << "\n";
        if(fabs(Inext - Iprev) < fabs(eps))
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ ДОСТИГНУТА!\n";
            stream << "ЗНАЧЕНИЕ ИНТЕГРАЛА РАВНО: I = " << Inext << "\n";
            flag = true;
        }
        else if (iterrNumber > MAX_ITERRATIONS_NUMBER )
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ НЕ ДОСТИГНУТА!\n";
            break;
        }
    }while(!flag);
    iterProcess = stream.str();
    return Inext;
}
double integralMidpoint(std::string expr, char argName, double leftLim, double rightLim,
                                 double eps, int &iterrNumber, std::string &iterProcess)
{
    double a = leftLim, b = rightLim;
    int n = 6;
    double h = (b - a) / n;
    double Inext = 0, Iprev = 0;
    iterrNumber = 0;
    const int MAX_ITERRATIONS_NUMBER = 5000;
    bool flag = false;
    std::stringstream stream;
    stream << "Выполняется численное интегрирование методом средних прямоугольников\n";
    do
    {
        stream << "Выполняется итеррация номер : " << ++iterrNumber << "\n";
        Iprev = Inext;
        Inext = 0;
        n *= 2;
        h = (b - a) / n;
        for(int i = 1; i <= n; ++i)
        {
            Inext += func(expr, argName, a + i * h - 0.5 * h);
        }
        Inext *= h;
        stream << "Количество интервалов разбиения на данной итеррации равно:\n n = " << n << "\n";
        stream << "Ширина i-го интервала на данной итеррации равна:\n hi = " << h << "\n";
        stream << "Значение интеграла на данной  итеррации равно:\n Ii = " << Inext << "\n";
        stream << "Погрешность данной  итеррации равна:\n ΔIi = Ii - Ii-1 = " << Inext << "\n";
        if(fabs(Inext - Iprev) < fabs(eps))
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ ДОСТИГНУТА!\n";
            stream << "ЗНАЧЕНИЕ ИНТЕГРАЛА РАВНО: I = " << Inext << "\n";
            flag = true;
        }
        else if (iterrNumber > MAX_ITERRATIONS_NUMBER )
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ НЕ ДОСТИГНУТА!\n";
            break;
        }
    }while(!flag);
    iterProcess = stream.str();
    return Inext;
}
double integralTrapezoids(std::string expr, char argName, double leftLim, double rightLim,
                                   double eps, int &iterrNumber, std::string &iterProcess)
{
    double a = leftLim, b = rightLim;
    int n = 6;
    double h = (b - a) / n;
    double Iprev = 0, Inext = 0;
    iterrNumber = 0;
    const int MAX_ITERRATIONS_NUMBER = 5000;
    bool flag = false;
    std::stringstream stream;
    stream << "Выполняется численное интегрирование методом трапеций \n";
    do
    {
        stream << "------------------------------------------------------\n" ;
        stream << "Выполняется итеррация номер : " << ++iterrNumber << "\n";
        Iprev = Inext;
        Inext = 0;
        n *= 2;
        h = (b - a) / n;
        for(int i = 0; i < n; ++i)
        {
            Inext += 0.5 * (func(expr, argName, a + i * h) +
                     func(expr, argName, a + (i + 1) * h));
        }
        Inext *= h;
        stream << "Количество интервалов разбиения на данной итеррации равно:\n n = " << n << "\n";
        stream << "Ширина i-го интервала на данной итеррации равна:\n hi = " << h << "\n";
        stream << "Значение интеграла на данной  итеррации равно:\n Ii = " << Inext << "\n";
        stream << "Погрешность данной  итеррации равна:\n ΔIi = Ii - Ii-1 = " << Inext << "\n";
        if(fabs(Iprev - Inext) < fabs(eps))
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ ДОСТИГНУТА!\n";
            stream << "ЗНАЧЕНИЕ ИНТЕГРАЛА РАВНО: I = " << Inext << "\n";
            flag = true;
        }
        else if (iterrNumber > MAX_ITERRATIONS_NUMBER )
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ НЕ ДОСТИГНУТА!\n";
            break;
        }
    }while(!flag);
    iterProcess = stream.str();
    return Inext;
}
double integralSimpson(std::string expr, char argName, double leftLim, double rightLim,
                                double eps, int &iterrNumber, std::string &iterProcess)
{
    double a = leftLim, b = rightLim;
    int n = 6;
    double h = (b - a) / n;
    double Inext = 0, Iprev = 0;
    iterrNumber = 0;
    const int MAX_ITERRATIONS_NUMBER = 5000;
    bool flag = false;
    std::stringstream stream;
    stream << "Выполняется численное интегрирование по методу Томаса Симпсона\n";
    do
    {
        stream << "------------------------------------------------------\n" ;
        stream << "Выполняется итеррация номер : " << ++iterrNumber << "\n";
        Iprev = Inext;
        Inext = 0;
        n *= 2;
        h = (b - a) / n;
        for(int k = 1; k < n; k += 2)
        {
            Inext += func(expr, argName, a + (k - 1) * h) +
                     4 * func(expr, argName, a + k * h) +
                     func(expr, argName, a + (k + 1) * h);
        }
        Inext *= h / 3;
        stream << "Количество интервалов разбиения на данной итеррации равно:\n n = " << n << "\n";
        stream << "Ширина i-го интервала на данной итеррации равна:\n hi = " << h << "\n";
        stream << "Значение интеграла на данной  итеррации равно:\n Ii = " << Inext << "\n";
        stream << "Погрешность данной  итеррации равна:\n ΔIi = Ii - Ii-1 = " << Inext << "\n";
        if(fabs(Inext - Iprev) < fabs(eps))
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ ДОСТИГНУТА!\n";
            stream << "ЗНАЧЕНИЕ ИНТЕГРАЛА РАВНО: I = " << Inext << "\n";
            flag = true;
        }
        else if (iterrNumber > MAX_ITERRATIONS_NUMBER )
        {
            stream << "------------------------------------------------------\n" ;
            stream << "НЕОБХОДИМАЯ ТОЧНОСТЬ НЕ ДОСТИГНУТА!\n";
            break;
        }
    }while(!flag);
    iterProcess = stream.str();
    return Inext;
}
double solveEquationChordsMethod(std::string expr, char argName, double a, double b,
                                  double eps, int &iterrNumber, std::string &iterProcess)
{
    const int MAX_ITERRATIONS_NUMBER = 5000;
    std::ostringstream stream;
    if(func(expr, argName, b) * func(expr, argName, a) > 0)
    {
        stream << "На границах указанного интервала изоляции корня значения функции имеют одниковые знаки!\n";
        stream << "Укажите такие значения границ интервала изоляции, чтобы левая часть уравнения f(X)=0 меняла знак\n";
        stream << "при изменении аргумента от одной границе к другой\n";
        iterProcess = stream.str();
        iterrNumber = 0;
        return 0.0;
    }
    if(funcDer2(expr, 'X', a) < 0)
    {
        double buf = a;
        a = b;
        b = buf;
    }
    double Xnext = a;
    bool flag = false;
    do
    {
        a = Xnext;
        Xnext = a - func(expr, argName, a) * (b - a) / (func(expr, argName, b) - func(expr, argName, a));
        iterrNumber++;
        stream << "------------------------------------------------------\n" ;
        stream << "Итеррация номер: " << iterrNumber << "\n";
        stream << "Значение корня на текущей итеррации:\n Xi = " << Xnext << "\n";
        stream << "Значение погрешности на текущей итеррации:\n ΔXi+1 - Xi = " << fabs(Xnext - a) << "\n";
        if(fabs(Xnext - a) < eps)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "РЕШЕНИЕ НАЙДЕНО : Х = " << Xnext << "\n";
            flag = true;
        }else if(iterrNumber > MAX_ITERRATIONS_NUMBER)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "РЕШЕНИЕ НЕ НАЙДЕНО!ПОРОБУЙТЕ СУЗИТЬ ИЛИ ИЗМЕНИТЬ ИНТЕРВАЛ ИЗОЛЯЦИИ КОРНЯ.\n";
            break;
        }

    }while(!flag);
    iterProcess = stream.str();
    return Xnext;
}
double solveEquationTangentsMethod(std::string expr, char argName, double a, double b,
                                    double eps, int &iterrNumber, std::string &iterProcess)
{
        const int MAX_ITERRATIONS_NUMBER = 5000;
        std::ostringstream stream;
        if(func(expr, argName, b) * func(expr, argName, a) > 0)
        {
            stream << "На границах указанного интервала изоляции корня значения функции имеют одниковые знаки!\n";
            stream << "Укажите такие значения границ интервала изоляции, чтобы левая часть уравнения f(X)=0 меняла знак\n";
            stream << "при изменении аргумента от одной границе к другой\n";
            iterProcess = stream.str();
            iterrNumber = 0;
            return 0.0;
        }
        if(funcDer2(expr, 'X', b) < 0)
        {
            double buf = a;
            a = b;
            b = buf;
        }
        double Xnext = b;
        bool flag = false;
        do
        {
            b = Xnext;
            double derivative = funcDer1(expr, argName, b);
            if(derivative == 0)
            {
                stream << "------------------------------------------------------\n" ;
                stream << "Текущее значение приближения совпало с точкой экстремума функции f(X) = 0.\n";
                stream << "Попробуйте изменить границы интервала изоляции корня.\n";
                break;
            }
            Xnext -= func(expr, argName, b) / derivative;
            iterrNumber++;
            stream << "------------------------------------------------------\n" ;
            stream << "Итеррация номер: " << iterrNumber << "\n";
            stream << "Значение корня на текущей итеррации:\n Xi = " << Xnext << "\n";
            stream << "Значение погрешности на текущей итеррации:\n ΔXi+1 - Xi = " << fabs(Xnext - b)  << "\n";
            if(fabs(Xnext - b) < eps)
            {
                stream << "------------------------------------------------------\n" ;
                stream << "РЕШЕНИЕ НАЙДЕНО : Х = " << Xnext << "\n";
                flag = true;
            }else if(iterrNumber > MAX_ITERRATIONS_NUMBER)
            {
                stream << "------------------------------------------------------\n" ;
                stream << "РЕШЕНИЕ НЕ НАЙДЕНО!ПОРОБУЙТЕ СУЗИТЬ ИЛИ ИЗМЕНИТЬ ИНТЕРВАЛ ИЗОЛЯЦИИ КОРНЯ.\n";
                break;
            }
        }while(!flag);
        iterProcess = stream.str();
        return Xnext;
}
double solveEquationBisectionMethod(std::string expr, char argName, double a, double b,
                                     double eps, int &iterrNumber, std::string &iterProcess)
{
    const int MAX_ITERRATIONS_NUMBER = 5000;
    std::ostringstream stream;
    if(func(expr, argName, b) * func(expr, argName, a) > 0)
    {
        stream << "На границах указанного интервала изоляции корня значения функции имеют одниковые знаки!\n";
        stream << "Укажите такие значения границ интервала изоляции, чтобы левая часть уравнения f(X)=0 меняла знак\n";
        stream << "при изменении аргумента от одной границе к другой\n";
        iterProcess = stream.str();
        iterrNumber = 0;
        return 0.0;
    }
    double c;
    bool flag = false;
    while(!flag)
    {
        c = 0.5 * (a + b);
        if(func(expr, argName, c) * func(expr, argName, a) > 0)
            a = c;
        else
            b = c;
        iterrNumber++;
        stream << "------------------------------------------------------\n" ;
        stream << "Итеррация номер: " << iterrNumber << "\n";
        stream << "Значение корня на текущей итеррации:\n Xi = " << c << "\n";
        stream << "Значение погрешности на текущей итеррации:\n ΔXi+1 - Xi = " << fabs(a - b)  << "\n";
        if(fabs(a - b) < eps)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "РЕШЕНИЕ НАЙДЕНО : Х = " << 0.5 * (a + b) << "\n";
            flag = true;
        }else if(iterrNumber > MAX_ITERRATIONS_NUMBER)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "РЕШЕНИЕ НЕ НАЙДЕНО!ПОРОБУЙТЕ СУЗИТЬ ИЛИ ИЗМЕНИТЬ ИНТЕРВАЛ ИЗОЛЯЦИИ КОРНЯ.\n";
            break;
        }
    }    
    iterProcess = stream.str();
    return 0.5 * (a + b);
}
double solveEquationCombinedMethod(std::string expr, char argName, double a, double b,
                                     double eps, int &iterrNumber, std::string &iterProcess)
{
    const int MAX_ITERRATIONS_NUMBER = 5000;
    std::ostringstream stream;
    if(func(expr, argName, b) * func(expr, argName, a) > 0)
    {
        stream << "На границах указанного интервала изоляции корня значения функции имеют одниковые знаки!\n";
        stream << "Укажите такие значения границ интервала изоляции, чтобы левая часть уравнения f(X)=0 меняла знак\n";
        stream << "при изменении аргумента от одной границе к другой\n";
        iterProcess = stream.str();
        iterrNumber = 0;
        return 0.0;
    }
    bool flag = false;
    if(funcDer2(expr, argName, b) < 0)
    {
        double buf = a;
        a = b;
        b = buf;
    }
    stream << "Выполняется поиск решения комбинированным методом хорд и касательных.\n";
    while(!flag)
    {
        a -= func(expr, argName, a) * (b - a) / (func(expr, argName, b) - func(expr, argName, a));
        double derivative = funcDer1(expr, argName, b);
        if(derivative == 0)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "Текущее значение приближения совпало с точкой экстремума функции f(X) = 0.\n";
            stream << "Попробуйте изменить границы интервала изоляции корня.\n";
            break;
        }
        b -= func(expr, argName, b) / derivative;
        iterrNumber++;
        stream << "------------------------------------------------------\n" ;
        stream << "Итеррация номер: " << iterrNumber << "\n";
        stream << "Значение корня на текущей итеррации:\n Xi = " << 0.5 * (a + b) << "\n";
        stream << "Значение погрешности на текущей итеррации:\n ΔXi+1 - Xi = " << fabs(a - b) << "\n";
        if(fabs(a - b) < eps)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "РЕШЕНИЕ НАЙДЕНО : Х = " << 0.5 * (a + b) << "\n";
            flag = true;
        }else if(iterrNumber > MAX_ITERRATIONS_NUMBER)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "РЕШЕНИЕ НЕ НАЙДЕНО!ПОРОБУЙТЕ СУЗИТЬ ИЛИ ИЗМЕНИТЬ ИНТЕРВАЛ ИЗОЛЯЦИИ КОРНЯ.\n";
            break;
        }
    }    
    iterProcess = stream.str();
    return 0.5 * (a + b);
}
double solveEquationIterationMethod(std::string expr, char argName, double a, double b, double eps,
                                    double lambda, int &iterrNumber, std::string &iterProcess)
{
    const int MAX_ITERRATIONS_NUMBER = 5000;
    std::ostringstream stream;
    if(func(expr, argName, b) * func(expr, argName, a) > 0)
    {
        stream << "На границах указанного интервала изоляции корня значения функции имеют одниковые знаки!\n";
        stream << "Укажите такие значения границ интервала изоляции, чтобы левая часть уравнения f(X)=0 меняла знак\n";
        stream << "при изменении аргумента от одной границе к другой\n";
        iterProcess = stream.str();
        iterrNumber = 0;
        return 0.0;
    }
    double Xnext = 0.5 * (a + b), Xprev;
    bool flag = false;
    stream << "Выполняется поиск решения методом простых итерраций.\n";
    stream << "Текущее значение коэффициента λ равно: " << lambda << "\n";
    /*stream << "Для сходимости итеррационного процесса коэффициент λ должен находиться в пределах: \n";
    double lambdaLim1 = 2 / funcDer1(expr, argName, a);
    double lambdaLim2 = 2 / funcDer1(expr, argName, b);
    //stream << "0.0 < λ < " << min(lambdaLim1, lambdaLim2) << "\n";
    stream << "0.0 < λ1 < " << lambdaLim1 << "\n";
    stream << "0.0 < λ2 < " << lambdaLim2 << "\n";*/
    do
    {
        Xprev = Xnext;
        Xnext = Xprev - lambda * func(expr, argName, Xprev);
        iterrNumber++;
        stream << "------------------------------------------------------\n" ;
        stream << "Итеррация номер: " << iterrNumber << "\n";
        stream << "Значение корня на текущей итеррации:\n Xi = " << Xnext << "\n";
        stream << "Значение погрешности на текущей итеррации:\n ΔXi+1 - Xi = " << fabs(Xnext - Xprev) << "\n";
        if(fabs(Xnext - Xprev) < eps)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "РЕШЕНИЕ НАЙДЕНО : Х = " << 0.5 * (Xnext + Xprev) << "\n";
            flag = true;
        }else if(iterrNumber > MAX_ITERRATIONS_NUMBER)
        {
            stream << "------------------------------------------------------\n" ;
            stream << "РЕШЕНИЕ НЕ НАЙДЕНО!ПОРОБУЙТЕ СУЗИТЬ ИЛИ ИЗМЕНИТЬ ИНТЕРВАЛ ИЗОЛЯЦИИ КОРНЯ.\n";
            break;
        }
    }while(!flag);
    iterProcess = stream.str();
    return 0.5 * (Xnext + Xprev);
}
std::vector<double> solveSystemOfNonLinearEquations(std::vector<std::string> equations, std::vector<std::pair<char, double>> variables,
                                     double eps, int &iterNumber, double &err, std:: string &iterProcess)
{
    //Эта функция решает систему нелинейных уравнений методом Ньютона
    //Максимальное число итерраций
    const int MAX_ITERATION_NUMBER = 5000;
    std::stringstream stream;
    //ofstream fo("result_SystemOfNonLinearEquations.txt");
    int eqNumber = equations.size();
    int numberOfIterrations = 0;
    std::vector<double> result;
    //result.reserve(eqNumber + 1);
    for(int i = 0; i < eqNumber; ++i)
        result.push_back(0.0);
    // Создаем вектор значений функции F1, F2 .. Fn
    std::vector<std::pair<char, double>> varList = variables;
    do
    {
        stream << "------------------------------------------------------\n" ;
        stream << "Выполняется итеррация номер : " << ++numberOfIterrations << "\n";
        // Создаем матрицу значений функции с обратным знаком
        if(numberOfIterrations > MAX_ITERATION_NUMBER)
        {
            stream << "РЕШЕНИЕ НЕ НАЙДЕНО! ПОПРОБУЙТЕ ЗАДАТЬ ДРУГОЕ НАЧАЛЬНОЕ ПРИБЛИЖЕНИЕ.\n";
            break;
        }
        Matrix<double> mtrF(eqNumber, 1, 0.0);
        // Заполняем ее
        for(int i = 0; i < eqNumber; ++i)
        {
            mtrF.setValueAt(i, 0, -multiVariableFunc(equations[i], varList));
        }
        //Создаем матрицу Якоби
        Matrix<double> mtrJacobi(eqNumber, eqNumber, 0.0);
        //Заполняем матрицу Якоби
        stream << "Элементы матрицы Якоби на текущей итеррации:\n";
        for(int i = 0; i < mtrJacobi.getRowCount(); ++i)
        {
            for(int j = 0; j < mtrJacobi.getColumnCount(); ++j)
            {
                double tmp = funcPartialDer1(equations[i], varList, varList[j].first);
                mtrJacobi.setValueAt(i, j, tmp);

                stream << "J[" << i << "][" << j << "] = " << mtrJacobi.getValueAt(i, j) << " ";
                //qDebug() << "J[" << i << "][" << j << "] = " << mtrJacobi.getValueAt(i, j);
            }
          stream << "\n";
        }
        // Проверяем неравенство ее детерминанта нулю и находим приращения
        if(mtrJacobi.determinant() == 0)
            throw std::runtime_error("Zero Jacobian!");
        Matrix<double> mtrInvJacobi = mtrJacobi.inverse();
        //Создаем матрицу-столбец невязок
        Matrix<double> mtrResiduals = mtrInvJacobi * mtrF;
        stream << "Вектор-столбец невязок на текущей итеррации:\n";
        for(int i = 0; i < mtrResiduals.getRowCount(); ++i)
            stream << "Δ" << variables[i].first << "i = " <<  mtrResiduals.getValueAt(i, 0) << "\n";
        err = 0;
        // Рассчитываем среднеквадратическую невязку
        for(int i = 0; i < mtrResiduals.getRowCount(); ++i)
             err += pow(mtrResiduals.getValueAt(i, 0), 2);
        err /= eqNumber;
        err = sqrt(err);
        stream << "Среднеквадратическая невязка на текущей итеррации:" << err << "\n";
        // Рассчитываем новые значения переменных
        stream << "Значения переменных на текущей итеррации: \n";
        for(int i = 0; i < eqNumber; ++i)
        {
            varList[i].second += mtrResiduals.getValueAt(i, 0);
            result[i] = varList[i].second;
            stream <<  variables[i].first << "i = " << result[i] << "\n";
           // qDebug() << "result[" << i << "] = " << result[i];
        }
        // Проверяем среднеквадратическую невязку
    }while (err > eps);
    iterNumber = numberOfIterrations;
    if(numberOfIterrations < MAX_ITERATION_NUMBER)
    {
        stream << "РЕШЕНИЕ НАЙДЕНО : \n";
        for(int i = 0; i < eqNumber; ++i)
           stream << variables[i].first << " = " << result[i] << "\n";
        stream << "Среднеквадратическая невязка : " << err << "\n";
    }
    stream << "------------------------------------------------------\n" ;
    iterProcess = stream.str();
    return result;
}
