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
        return (1 / x - tan(x));
    }
};

class f_2 : public Function
{
public:
    double operator()(double x)
    {
        return (1 / x - pow(2 , x));
    }
};

class f_3 : public Function
{
public:
    double operator()(double x)
    {
        return (pow(2, -x) + exp(x) + 2 * cos(x) - 6);
    }
};

class f_4 : public Function
{
public:
    double operator()(double x)
    {
        return  (x*x*x + 4*x*x + 3*x + 5) / (2*x*x*x -9*x*x +  18*x - 2);
    }
};

int main()
{
    cout << "Question B : Test the bisection method " << endl;

    f_1 f1;
    Bisection solve_1( f1, 0,  pi/2, eps, eps, 50);
    double x_1 = solve_1.solve();
    double y_1 = f1(x_1);
    cout << "The approximate root of function is " << x_1 <<  endl;
    cout << "Error : " << fabs(f1(x_1)) << endl;

    f_2 f2;
    Bisection solve_2( f2, 0, 1, eps, eps, 50);
    double x_2 = solve_2.solve();
    double y_2 = f2(x_2);
    cout << "The approximate root of function is " << x_2 << endl;
    cout << "Error : " << fabs(f2(x_2)) << endl;

    f_3 f3;
    Bisection solve_3( f3, 1, 3, eps, eps, 50);
    double x_3 = solve_3.solve();
    double y_3 = f3(x_3);
    cout << "The approximate root of function is "  << x_3 << endl;
    cout << "Error : " << fabs(f3(x_3)) << endl;
    
    f_4 f4;
    Bisection solve_4(f4, 0, 4, eps, eps, 50);
    double x_4 = solve_4.solve();
    double y_4 = f4(x_4);
    cout <<  "The approximate root of function is " << x_4 << endl;
    cout << "Error : " << fabs(f4(x_4)) << endl;
    
    return 0;
}