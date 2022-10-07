#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include "Header.h"
using namespace std;

const double eps = numeric_limits<double>::epsilon();
const double pi = acos(-1), l = 89, h = 49, b1 = (11.5 / 180) * pi;

class _f1 : public Function
{
private:
	double D = 50;

public:
	double operator()(double a)
	{
		return (l * sin(b1) * sin(a) * cos(a) + l * cos(b1) * sin(a) * sin(a) - (((h + 0.5 * D) * sin(b1) - 0.5 * D * tan(b1)) * cos(a)) - ((h + 0.5 * D) * cos(b1) - 0.5 * D) * sin(a));
	};
};

class _f2 : public Function
{
private:
	double D = 30;

public:
	double operator()(double a)
	{
		return (l * sin(b1) * sin(a) * cos(a) + l * cos(b1) * sin(a) * sin(a) - (((h + 0.5 * D) * sin(b1) - 0.5 * D * tan(b1)) * cos(a)) - ((h + 0.5 * D) * cos(b1) - 0.5 * D) * sin(a));
	}

};

// convert radian to degree 
double Convert(double radian)
{
	return(radian * (180 / pi));
}

int main()
{
	// since 0 < a < pi/2
	_f1 f1;
	cout << "Using Bisection Method" << endl;
	Bisection solve_1(f1, 0, pi/2, eps, eps, 100);
	double x_1 = solve_1.solve();
	cout << " (a) The value of alpha is " << (x_1 * 180) / pi << endl;

	cout << endl;

	_f2 f2;
	double initialguess = (33 * pi) / 180;
	cout << "Using Newton's Method" << endl;
	Newton solve_2(f2, initialguess , eps, 10);
	double x_2 = solve_2.solve();
	cout << " (b) The value of alpha is " << (x_2 * 180) / pi << endl;

	cout << endl;

	cout << "Using Secant Method" << endl;
	double x0 = Convert(20) , x1 = Convert(50);
	Secant solve_3(f1, x0, x1, eps, eps, 100);
	double x_3 = solve_3.solve();
	Secant solve_4(f1, pi / 3, pi / 2, eps, eps, 100);
	double x_4 = solve_4.solve();
	cout << " (c) i. When chossing x0 = 20 x1 = 50 respectively, the value of alpha is " << (x_3 * 180) / pi << endl;
	cout << " (c) ii. When chossing x0 = 60 x1 = 90 respectively, the value of alpha is " << (x_4 * 180) / pi << endl;

	cout << endl;

	return 0;
}
