#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include "Interpolation.h"
using namespace std;
const double pi = acos(-1);

class f_2
{
public:
    double operator () (double x)
    {
        return (1.0 / (1.0 + 25.0 * x * x));
    }
};

int main()
{
    for (int n = 5; n <= 20; n += 5)
    {
        cout << noshowpos;
        f_2 f2;
        vector<double> x;
        vector<double> f;

        for (int i = 1; i <= n; i++)
        {
            x.push_back(cos((2.0 * i - 1.0) / (2.0 * n) * pi));
            f.push_back(f2(x[i - 1]));
        }
        vector<double> c = getNewton(x, f);
        cout << "p" << n << "(f;x ) = ";
        NewtonPolynomial(c, x);
        cout << endl;
    }
    return 0;
}