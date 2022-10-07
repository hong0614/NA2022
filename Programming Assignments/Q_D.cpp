#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include "Header.h"
using namespace std;

const double eps = numeric_limits<double>::epsilon();
const double pi = acos(-1);

class f_1 : public Function
{
public:
    double operator()(double x)
    {
        return (sin(x / 2) - 1);
    }
};

class f_2 : public Function
{
public:
    double operator()(double x)
    {
        return (exp(x) - tan(x));
    }
};

class f_3 : public Function
{
public:
    double operator()(double x)
    {
        return (x * x * x - 12 * x * x + 3 * x + 1);
    }
};

int main()
{
    cout << "Question D : Test the secant method " << endl;

    f_1 f1;
    Secant solve_1(f1, 0, pi / 2, eps, eps, 10000);
    double x_1 = solve_1.solve();
    cout << "The approximate root of function is " << x_1 << endl;
    cout << "Error : " << fabs(f1(x_1)) << endl;

    f_2 f2;
    Secant solve_2(f2, 1, 1.4, eps, eps, 100);
    double x_2 = solve_2.solve();
    cout << "The approximate root of function is " << x_2 << endl;
    cout << "Error : " << fabs(f2(x_2)) << endl;

    f_3 f3;
    Secant solve_3(f3, 0, -0.5, eps, eps, 100);
    double x_3 = solve_3.solve();
    cout << "The approximate root of function is " << x_3 << endl;
    cout << "Error : " << fabs(f3(x_3)) << endl;

    return 0;
}