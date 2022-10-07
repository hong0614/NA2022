#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include "Header.h"
using namespace std;

const double eps = numeric_limits<double>::epsilon();
const double pi = acos(-1), L = 10, r = 1, V = 12.4;

class _f : public Function
{
public:
	double operator()(double h)
	{
		return (L * (0.5 * pi * r * r - r * r * asin(h / r) - h * sqrt(r * r - h * h)) - V);
	}

};

int main()
{
	// since 0 < h < r = 1
	_f f;

	cout << "Using Bisection Method" << endl;
	Bisection solve_1(f, 0, 1, 0.01, 0.0000001, 200);
	double x_1 = solve_1.solve();
	double y_1 = f(x_1);
	cout << "The approximate root  is " << x_1 << "ft" << endl;
	cout << "Volume Error : " << fabs(f(x_1)) << endl;

	cout<< endl;

	cout << "Using Newton's Method" << endl;
	Newton solve_2(f, 0, 0.01, 20);
	double x_2 = solve_2.solve();
	double y_2 = f(x_2);
	cout << "The approximate root is " << x_2 << "ft" << endl;
	cout << "Volume error : " << fabs(y_2) << endl;

	cout << endl;

	cout << "Using Secant Method" << endl;
	Secant solve_3(f, 0, 1, 0.01, 0.000001, 100);
	double x_3 = solve_3.solve();
	cout << "The approximate root of function is " << x_3 << "ft" << endl;
	cout << "Volume Error : " << fabs(f(x_3)) << endl;

	cout << endl;

	return 0;
}