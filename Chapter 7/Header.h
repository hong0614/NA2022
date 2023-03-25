#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include<vector>
using namespace std;

class Function
{
public:
	virtual double operator()(double x, double y) = 0;

};

class u_1 :public Function
{
public:
	double operator()(double x, double y)
	{
		return exp(y + sin(x));
	}
};

class ux_1 :public Function
{
public:
	double operator()(double x, double y)
	{
		return exp( y + sin(x)) * (cos(x));
	}
};

class uy_1 :public Function
{
public:
	double operator()(double x, double y)
	{
		return exp(y + sin(x));
	}
};

class f_1 :public Function
{
public:
	double operator()(double  x, double y)
	{
		return -exp(y + sin(x)) * (1 + cos(x) * cos(x) - sin(x));
	}
};

void getval(int N)
{

	vector<double> x , y;
	vector<vector<double>> F;
	f_1 f1;
	double h = 1.0 / N;
	F.resize(N-1);
	for (int i = 1; i < N; i++)
	{
		F[i-1].resize(N-1);
		x.push_back( i * h);
		for(int j = 1 ;j< N;j++)
		{ 
			y.push_back( j * h);
			F[i-1][j-1] = f1(x[i-1], y[j-1]);
		}
	}

	
	for (int i = 0; i < N - 1; i++)
	{
		for (int j = 0; j < N - 1; j++)
		{
			cout << F[i][j] << " ";
		}
		cout << endl;
	}

	vector<vector<double>> A;
	A.resize(N);
	/*for (int i = 0; i <= N; i++)
	{
		A[i][] = i;
		cout << A[i][i] << endl;
	}*/
}

void getcoef(int N)
{
	vector<vector<double>> A;
	A.resize((N - 1)* (N - 1));
	for (int i = 0; i < (N - 1) * (N - 1); i++)
	{
		A[i].resize((N - 1) * (N - 1));
		for (int j = 0; j < i+1; j++)
		{
			if (i == j)
			{
				A[i][j] = 4.0;
			}

			else if (j == i - (N-1))
			{
				A[i][j] = -1.0;
			}
			
			else if (i == j + 1)
			{
				if ((j+1) % (N-1) == 0)
				{
					A[i][j] = 0;
				}
				else 
				{
					A[i][j] = -1.0;
				}
			}
			else
			{
				A[i][j] = 0;
			}
			A[j][i] = A[i][j];
		}
	}


	for (int i = 0; i < ((N - 1) * (N - 1)); i++)
	{
		for (int j = 0; j < ((N - 1) * (N - 1)); j++)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
}

class PP
{
public:
	double x, y;
};
class P0int
{
public:
	double x, y;

	vector<PP> getpoints(int n)const
	{
		vector<PP> points;
		double h = 1.0 / n;
		for (int i = 0; i < n; i++)
		{
			PP p = { i * h,i * h };
			points.push_back(p);
		}
		return points;
	}
};
//P0int X_ = {0,0};
//vector<double> points = X_.getpoints(3);
//for (int i = 0; i < N; i++)
//{
//    cout << points[i].x << "and" << points[i].y << endl;
//}


void haha()
{
	//getval(4);
	/*    std::cin >> n;
A.resize(n);
for (int i = 0; i < n; i++)
{
	A[i].resize(n);
	for (int j = 0; j < n; j++)
	{
		A[i][j] = i * j;
		cout << A[i][j] << " " ;
	}
	cout << endl;
}*/
	vector<vector<double>>getA(3);
	//getcoef(4);
}