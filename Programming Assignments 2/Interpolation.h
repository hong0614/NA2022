#pragma once
#include <vector>
#include <iostream>
using namespace std;

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
		for (int j = 0; j < n - i; j++) temp[j] = (f[j + 1] - f[j]) / (x[j + i] - x[j]);
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
		cout << showpos << c[i];
		for (int j = 0; j < i; j++) cout << "(x" << showpos << -x[j] << ")";
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

void HermitePoly(vector<double> x, vector<double> y, vector<double> _x)
{
	vector<vector<double>> f;
	f.resize(x.size() + 1);
	for (int i = 0; i <= x.size(); i++)
	{
		f[i].resize(x.size() + 1);
	}
	for (int i = 0; i < x.size(); i++)
	{
		f[i][0] = y[i];
	}
	for (int i = 1; i < x.size(); i = i + 2)
	{
		f[i][1] = _x[i];
	}
	for (int i = 2; i < x.size(); i = i + 2)
	{
		f[i][1] = (f[i][0] - f[i - 1][0]) / (x[i] - x[i - 2]);
	}
	for (int i = 2; i < x.size(); i++)
	{
		for (int j = i; j < x.size(); j++)
		{
			f[j][i] = (f[j][i - 1] - f[j - 1][i - 1]) / (x[j] - x[j - i]);
		}
	}
	for (int i = 0; i < f.size(); i++)
	{
		cout << showpos << f[i][i];
		for (int j = 0; j < i; j++) cout << "(x" << showpos << -x[j] << ")";
	}
	cout << '\n';
};
