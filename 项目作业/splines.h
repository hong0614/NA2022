#pragma once
#include <iostream> 
#include<limits>
#include <vector>
#include <cmath>
#include <map>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
const double eps = numeric_limits<double>::epsilon();

class Polynomial
{
public:
	int n; //degree
	vector <double> c ; //coefficient of polynomial
	double begin;
	double end;
	double getval(double x) 
	{
		double sum = c[0];
		for (int i = 1; i <= n; i++)
		{
			sum += c[i] * pow(x - begin, i);
		}
		return sum;
	}
};


class ppf
{
public:
	int N;
	vector <double> _x, _fx;
	vector<Polynomial> poly;
	virtual double solve()
	{
		return 0;
	}	
};


class completecubicspline :public ppf
{
public:
	double m1, mN;
	vector<double> xp, fxp;

	virtual double solve()
	{
		MatrixXd A = MatrixXd::Zero(N + 1, N + 1);
		A(0, 0) = 1; A(1, 1) = 1; A(N, N) = 1;
		for (int i = 2; i < N; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				if (j == i - 1)
				{
					A(i, j) = (xp[i + 1] - xp[i]) / (xp[i + 1] - xp[i - 1]);
				}
				else if (j == i)
				{
					A(i, j) = 2;
				}
				else if (j == i + 1)
				{
					A(i, j) = (xp[i] - xp[i - 1]) / (xp[i + 1] - xp[i - 1]);
				}
			}
		}
		VectorXd b = VectorXd::Zero(N + 1);
		b(0) = 1, b(1) = m1; b(N) = mN;
		for (int i = 2; i < N; i++)
		{
			double v_1 = (xp[i] - xp[i - 1]) / (xp[i + 1] - xp[i - 1]) * (fxp[i + 1] - fxp[i]) / (xp[i + 1] - xp[i]);
			double v_2 = (xp[i + 1] - xp[i]) / (xp[i + 1] - xp[i - 1]) * (fxp[i] - fxp[i - 1]) / (xp[i] - xp[i - 1]);
			b(i) = 3 * (v_1 + v_2);
		}

		VectorXd mA = VectorXd::Zero(N + 1);
		mA = A.colPivHouseholderQr().solve(b);

		vector<double> m;
		m.push_back(0);
		for (int i = 1; i <= N; i++)
		{
			m.push_back(mA(i));
			//cout << mA(i) << endl;
		}

		Polynomial p1;
		p1.n = 3; p1.begin = 0; p1.end = 0; p1.n = 3;
		vector<double> p_c{ 0,0,0,0 };
		p1.c = p_c;
		poly.push_back(p1);

		for (int i = 1; i < N; i++)
		{
			double K_i = (fxp[i + 1] - fxp[i]) / (xp[i + 1] - xp[i]);
			Polynomial _poly;
			_poly.n = 3;
			_poly.begin = xp[i];
			_poly.end = xp[i + 1];
			_poly.c.push_back(fxp[i]);
			_poly.c.push_back(m[i]);
			_poly.c.push_back((3 * K_i - 2 * m[i] - m[i + 1]) / (xp[i + 1] - xp[i]));
			_poly.c.push_back((m[i] + m[i + 1] - 2 * K_i) / pow(xp[i + 1] - xp[i], 2));
			poly.push_back(_poly);
			//cout << "The " << i << "term " << endl;
			//cout << poly[i].n << endl;
			//cout << poly[i].c[0] << endl;
			//cout << poly[i].c[1] << endl;
			//cout << poly[i].c[2] << endl;
			//cout << poly[i].c[3] << endl;
		}
	}

	double value(double x)
	{
		for (int i = 1; i < N ; i++)
		{
			if (x >= xp[i] && x <= xp[i + 1])
			{
				return poly[i].getval(x);
			}
		}
	}

};




//completecubicspline* a = new completecubicspline();
//a->solve();

class knotcubicspline:public ppf
{
public:
	int N;
	vector<double> xp, fxp; //points x , values of f(x)
	vector<Polynomial> poly;

	virtual double solve() override
	{
		MatrixXd A = MatrixXd::Zero(N + 1, N + 1);
		int i, j;
		A(0, 0) = 1; A(1, 2) = -1; A(N, N - 1) = -1;
		A(1, 1) = (xp[3] - xp[2]) / (xp[3] - xp[1]);
		A(1, 3) = (xp[2] - xp[1]) / (xp[3] - xp[1]);
		A(N, N - 2) = (xp[N] - xp[N - 1]) / (xp[N] - xp[N - 2]);
		A(N, N) = (xp[N - 1] - xp[N - 2]) / (xp[N] - xp[N - 2]);
		for (i = 2; i < N; i++)
		{
			for (j = 1; j <= N; j++)
			{
				if (j == i - 1)
				{
					A(i, j) = (xp[i] - xp[i - 1]) / (xp[i + 1] - xp[i - 1]);
				}
				else if (j == i)
				{
					A(i, j) = 2;
				}
				else if (j == i + 1)
				{
					A(i, j) = (xp[i + 1] - xp[i]) / (xp[i + 1] - xp[i - 1]);
				}
			}
		}

		VectorXd b = VectorXd::Zero(N + 1);
		b(0) = 1; b(1) = 0; b(N) = 0;
		for (int i = 2; i < N; i++)
		{
			double miu1 = 6 * (fxp[i] - fxp[i - 1]) / (xp[i] - xp[i - 1]);
			double miu2 = 6 * (fxp[i + 1] - fxp[i]) / (xp[i + 1] - xp[i]);
			b(i) = (miu2 - miu1) / (xp[i + 1] - xp[i - 1]);

		}

		VectorXd mA = VectorXd::Zero(N + 1);
		mA = A.colPivHouseholderQr().solve(b);
		vector<double> M;
		M.push_back(0);
		for (int i = 1; i <= N; i++)
		{
			M.push_back(mA(i));
		}

		Polynomial p1;
		p1.n = 3;
		p1.begin = 0; p1.end = 0;
		vector<double> p1_c{ 0,0,0,0 };
		p1.c = p1_c;
		poly.push_back(p1);
		for (int i = 1; i < N; i++)
		{
			double K_i = (fxp[i + 1] - fxp[i]) / (xp[i + 1] - xp[i]);
			Polynomial _poly;
			_poly.n = 3;
			_poly.begin = xp[i]; _poly.end = xp[i + 1];
			_poly.c.push_back(fxp[i]);
			_poly.c.push_back((fxp[i + 1] - fxp[i]) / (xp[i + 1] - xp[i]) - (M[i + 1] + 2 * M[i]) * (xp[i + 1] - xp[i]) / 6);
			_poly.c.push_back(M[i] / 2);
			_poly.c.push_back(((M[i + 1] - M[i]) / (xp[i + 1] - xp[i])) / 6);
			poly.push_back(_poly);
			//cout << poly[i].c[0] << endl;
			//cout << poly[i].c[1] << endl;
			//cout << poly[i].c[2] << endl;
			//cout << poly[i].c[3] << endl;
		}

	}
	double value(double x)
	{
		for (int i = 1; i < N; i++)
		{
			if (x <= xp[i + 1] && x >= xp[i])
			{
				return poly[i].getval(x);
			}
		}
	}

};



class cubicsplinewith2nd :public ppf
{
public:
	int N;
	double m_1, m_n;
	vector<Polynomial>poly;
	vector<double> xp, fxp;

	virtual double solve() override
	{
		MatrixXd A = MatrixXd::Zero(N + 1, N + 1);
		A(0, 0) = 1; A(1, 1) = 1; A(N, N) = 1;
		for (int i = 2; i < N ; i++)
		{
			for (int j = 1; j <= N; j++)
			{
				if (j == i - 1)
				{
					A(i, j) = (xp[i] - xp[i - 1]) / (xp[i + 1] - xp[i - 1]);
				}
				else if (j == i)
				{
					A(i, j) = 2;
				}
				else if (j == i + 1)
				{
					A(i, j) = (xp[i + 1] - xp[i]) / (xp[i + 1] - xp[i - 1]);
				}
			}
		}

		VectorXd b = VectorXd::Zero(N + 1);
		b(0) = 1;b(1) = m_1; b(N) = m_n;
		for (int i = 2; i < N ; i++)
		{
			double temp1 = 6 * (fxp[i] - fxp[i - 1]) / (xp[i] - xp[i - 1]);
			double temp2 = 6 * (fxp[i + 1] - fxp[i]) / (xp[i + 1] - xp[i]);
			b(i) = (temp2 - temp1) / (xp[i + 1] - xp[i - 1]);
		}
		VectorXd m_A = VectorXd::Zero(N + 1);
		m_A = A.colPivHouseholderQr().solve(b);
		vector<double> M;
		M.push_back(0);
		for (int i = 1; i <= N; i++)
		{
			M.push_back(m_A(i));
		}
		Polynomial p1;
		p1.n = 3;
		p1.begin = 0; p1.end = 0;
		vector<double> p1_c{ 0,0,0,0 };
		p1.c = p1_c;
		poly.push_back(p1);
		for (int i = 1; i < N ; i++)
		{
			Polynomial _poly;
			_poly.n = 3;
			_poly.begin = xp[i]; _poly.end = xp[i + 1];
			_poly.c.push_back(fxp[i]);
			_poly.c.push_back((fxp[i + 1] - fxp[i]) / (xp[i + 1] - xp[i]) - (M[i + 1] + 2 * M[i]) * (xp[i + 1] - xp[i]) / 6);
			_poly.c.push_back(M[i] / 2);
			_poly.c.push_back(((M[i + 1] - M[i]) / (xp[i + 1] - xp[i])) / 6);
			poly.push_back(_poly);
			//cout << poly[i].c[0] << endl;
			//cout << poly[i].c[1] << endl;
			//cout << poly[i].c[2] << endl;
			//cout << poly[i].c[3] << endl;
		}
	}

	double value(double x)
	{
		for (int i = 1; i < N; i++)
		{
			if (x >= xp[i] && x <= xp[i + 1])
			{
				return poly[i].getval(x);
			}
		}
	}

};


