#pragma once
#include <vector>
#include <iostream>

using namespace std;

int factorial(int num)
{
	if (num == 0) //基本情况返回1;
		return 1;
	else
		return num * factorial(num - 1);
}

class Function
{
	const double _eps = 10 * numeric_limits<double>::epsilon();
public:
	virtual double operator()(double x) = 0;

	virtual double diff(double x) //calculate derivative
	{
		return ((*this)(x + _eps) - (*this)(x)) / (_eps);
	}
};

vector<double>getNewton(const vector<double>& x, vector<double>f)
{
	int n = x.size();
	vector<double> c(n), temp(n);
	c[0] = f[0];
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j < n - i; j++) temp[j] = (f[j + 1] - f[j]) / (x[j + 1] - x[j]);
		f = temp;
		c[i] = f[0];		
	}
	return c;
}

void NewtonPolynomial(vector<double>& c, const vector<double>& x)
{
	int n = x.size();
	cout << c[0];
	for (int i = 1; i < n; i++)
	{
		c[i] = c[i] / factorial(i);
		cout << showpos << c[i];
		for (int j = 0; j < i; j++) cout << "(x" << showpos << - x[j] << ")";
	}
	cout << '\n';
}

double getvalue(const vector<double>& c, const vector<double>& x, double xval)
{
	int n = c.size();
	double fx = c[0], poly = 1.0;
	for (int i = 1; i < n; i++)
	{
		poly *= (xval - x[i - 1]);
		fx += c[i] * poly;
	}
	return fx;
}