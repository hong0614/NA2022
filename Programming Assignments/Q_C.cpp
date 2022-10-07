#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include "Header.h"
using namespace std;

const double eps = numeric_limits<double>::epsilon();
const double pi = acos(-1);

class _f : public Function
{
public:
	double operator()(double x)
	{
		return (x - tan(x));
	}
};

int main()
{
	cout << "Question C : Test the Newton's method " << endl;
	
	_f f;
	
	Newton solve_1(f, 4.5 , eps, 10);
	double x_1 = solve_1.solve();
	double y_1= f(x_1);
	cout << "The approximate root which near 4.5 is " << x_1 << " with error : " << fabs(y_1) << endl;

	Newton solve_2(f, 7.7, eps, 10);
	double x_2 = solve_2.solve();
	double y_2 = f(x_2);
	cout << "The approximate root which near 7.7 is " << x_2 << " with error : " << fabs(y_2) << endl;

	return 0;
}