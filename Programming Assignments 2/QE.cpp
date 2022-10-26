#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include "Interpolation.h"
using namespace std;

int main()
{
    vector<double> x = { 0.0, 6.0, 10.0, 13.0, 17.0, 20.0, 28.0 };
    vector<double> sp1 = { 6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7 };
    vector<double> sp2 = { 6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89 };

    vector<double> c1 = getNewton(x, sp1);
    vector<double> c2 = getNewton(x, sp2);

    cout << "The polynomial of Sp1 is f1(x) = ";
    NewtonPolynomial(c1, x); 
    cout << endl;

    cout << "The polynomial of Sp2 is f2(x) = ";
    NewtonPolynomial(c2, x);

    double n = 43;
    cout << noshowpos << endl;

    double f1 = getvalue(c1, x, n);
    cout << "When x = 43 , f1(43) = " << f1 << endl;

    double f2 = getvalue(c2, x, n);
    cout << "When x = 43 , f2(43) = " << f2 << endl;

}