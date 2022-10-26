#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include "Interpolation.h"
using namespace std;

class f_1 : public Function
{
public:
    double operator()(double x)
    {
        return (1 / (1 + x * x));
    }
};

int main()
{
    for (int n = 2; n <= 8; n += 2)
    {
        cout << noshowpos;
        f_1 f1;
        vector<double> x;
        vector<double> f;
        for (int i = 0; i <= n; i++)
        {
            x.push_back(-5.0 + (10.0 * i / n));
            f.push_back(f1(x[i]));
        }
        vector<double> c = getNewton(x, f);
        cout << "p" << n << "(f;x ) = ";
        NewtonPolynomial(c, x);
        cout << endl;

    }
    return 0;
}