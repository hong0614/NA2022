#pragma once
#include <iostream> 
#include<limits>
#include <vector>
#include <cmath>
#include <map>
#include <Eigen/Dense>
#include "function.h"
#include "splines.h"
using namespace std;
using namespace Eigen;

class f1 :public Function
{
public:
    double operator()(double x)
    {
        return (1 / ( 1+ 25 * x * x));
    }
};


void A_completecubicspline(int N)
{
    f1 _f1;
    completecubicspline A1;
    A1.N = N;
    A1.m1 = _f1.rdiff(-1.0);
    A1.mN = _f1.ldiff(1.0);
    //string number = to_string(N);

    for (int i = 1; i <= N; i++)
    {
        double x_i = -1 + 2.0 * (double)(i - 1) / (N - 1);
        A1.xp.push_back(x_i);
        A1.fxp.push_back(_f1(x_i));
        // cout << x_i  << "and"  <<  _f1(x_i) << endl;
    }

    A1.solve();

    double error = 0;

    for (int i = 1; i <= N - 1; i++)
    {
        double temp = abs(A1.value((A1.xp[i] + A1.xp[i + 1]) / 2) - _f1((A1.xp[i] + A1.xp[i + 1]) / 2));
        if (error < temp)
        {
            error = temp;
        }
    }

    cout << "The error maximun norm  of cubic complete spline for " << N << " is " << error << endl;

}

void A_cubicsplinewith2nd(int N)
{
    f1 _f1;
    cubicsplinewith2nd A2;
    A2.N = N;
    A2.m_1 = _f1.drdiff(-1.0);
    A2.m_n = _f1.dlddiff(1.0);
    for (int i = 1; i <= N; i++)
    {
        double x_i = -1 + 2.0 * (double)(i - 1) / (N - 1);
        A2.xp.push_back(x_i);
        A2.fxp.push_back(_f1(x_i));
        // cout << x_i  << "and"  <<  _f1(x_i) << endl;
    }
    A2.solve();


    double error = 0;

    for (int i = 1; i <= N - 1; i++)
    {
        double temp = abs(A2.value((A2.xp[i] + A2.xp[i + 1]) / 2) - _f1((A2.xp[i] + A2.xp[i + 1]) / 2));
        if (error < temp)
        {
            error = temp;
        }
    }

    cout << "The error maximun norm  of cubic spline with specified second derivative for  " << N << " is " << error << endl;
}

void A_knotcubicspline(int N)
{
    knotcubicspline A3;
    f1 _f1;
    A3.N = N;

    for (int i = 1; i <= N; i++)
    {
        double x_i = -1 + 2.0 * (double)(i - 1) / (N - 1);
        A3.xp.push_back(x_i);
        A3.fxp.push_back(_f1(x_i));
        // cout << x_i  << "and"  <<  _f1(x_i) << endl;
    }
    A3.solve();


    double error = 0;

    for (int i = 1; i <= N - 1; i++)
    {
        double temp = abs(A3.value((A3.xp[i] + A3.xp[i + 1]) / 2) - _f1((A3.xp[i] + A3.xp[i + 1]) / 2));
        if (error < temp)
        {
            error = temp;
        }
    }

    cout << "The error maximun norm  of not a knot complete spline for " << N  << " is " << error << endl;

}

int main()
{
    A_completecubicspline(6);
    A_completecubicspline(11);
    A_completecubicspline(21);
    A_completecubicspline(41);
    A_completecubicspline(81);
    cout << endl;
    A_cubicsplinewith2nd(6);
    A_cubicsplinewith2nd(11);
    A_cubicsplinewith2nd(21);
    A_cubicsplinewith2nd(41);
    A_cubicsplinewith2nd(81);
    cout << endl;
    A_knotcubicspline(6);
    A_knotcubicspline(11);
    A_knotcubicspline(21);
    A_knotcubicspline(41);
    A_knotcubicspline(81);

    return 0;
}
