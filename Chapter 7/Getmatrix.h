#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include<vector>
#include"Header.h"
#include"Eigen/Dense"
using namespace std;
using namespace Eigen;

class Point
{
	vector<double> x, y;
	int getpoint(int N)
	{
		double h = 1.0 / N;
		for (int i = 1; i <= N; i++)
		{
			x.push_back(i * h);
			y.push_back(i * h);
		}
	}
};

class PDE_Matrix
{
private:
	int N;
	vector<vector<double>> A, F;

public:
	vector<vector<double>>getA(int N);
	vector<vector<double>>getF(Function& f1,int N);
	MatrixXd get_eigenA(vector<vector<double>>A, int N);
	VectorXd get_eigenF(vector<vector<double>>A, int N);

};

vector<vector<double>> PDE_Matrix::getA(int N)
{
	vector<vector<double>> A;
	A.resize((N - 1) * (N - 1));
	for (int i = 0; i < (N - 1) * (N - 1); i++)
	{
		A[i].resize((N - 1) * (N - 1));
		for (int j = 0; j < i + 1; j++)
		{
			if (i == j)
			{
				A[i][j] = 4.0;
			}

			else if (j == i - (N - 1))
			{
				A[i][j] = -1.0;
			}

			else if (i == j + 1)
			{
				if ((j + 1) % (N - 1) == 0)
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

	//for (int i = 0; i < (N - 1) * (N - 1); i++)
	//{
	//	for (int j = 0; j < (N - 1) * (N - 1); j++)
	//	{
	//		cout << A[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	return A;
}

MatrixXd PDE_Matrix::get_eigenA(vector<vector<double>>A, int N)
{
	MatrixXd A2D = MatrixXd::Zero((N - 1)* (N - 1), (N - 1) * (N - 1));
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; i < N; i++)
		{
			A2D(i,j) = A[i][j];
		}
	}

	return A2D;
}


vector<vector<double>> PDE_Matrix::getF(Function& func, int N)
{
	vector<vector<double>> F;
	double h = 1.0 / N;
	F.resize(N - 1);
	vector<double> x, y;
	for (int i = 1; i < N; i++)
	{
		x.push_back(i * h);
		y.push_back(i * h);
	}
	for (int i = 1; i < N; i++)
	{
		F[i - 1].resize(N - 1);
		for (int j = 1; j < N; j++)
		{
			F[i - 1][j - 1] = h * h * func(x[i - 1], y[j - 1]);
		}
	}

	//for (int i = 0; i < N - 1; i++)
	//{
	//	for (int j = 0; j < N - 1; j++)
	//	{
	//		cout << F[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	return F;
}

VectorXd PDE_Matrix::get_eigenF(vector<vector<double>>A, int N)
{
	VectorXd F2D = VectorXd::Zero(N - 1);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; i < N; i++)
		{
			F2D((i)+(j) * (N - 1)) = F[i + 1][j + 1];
		}
	}

	return F2D;
}
class Condition :public PDE_Matrix
{
private:
	int N;
	vector<vector<double>> A, F;

public:
	//Condition(Function& func, Function& ddfunc, const string& type_name, int N);

	Condition(Function& func, Function& ddfunc, const string& type_name, int N)
	{
		PDE_Matrix p1;
		vector<vector<double>>A = p1.getA(N);
		vector<vector<double>>F = p1.getF(ddfunc,N);
		int m = N - 1;
		double h = 1.0 / N;
		vector<double> x, y;

		for (int i = 0; i <= N; i++)
		{
			x.push_back(i * h);
			y.push_back(i * h);
		}

		if (type_name == "Dirichlet")
		{
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < m; j++)
				{
					//F[i][j] = ddfunc(x[i+1], y[j+1]);

					if (i == 0)
					{
						F[i][j] += func(x[i], y[j+1]);
					}

					if (j == 0)
					{
						F[i][j] += func(x[i+1], y[j]);
					}

					if (i == m - 1)
					{
						F[i][j] += func(x[i + 2], y[j + 1]);
					}

					if (j == m - 1)
					{
						F[i][j] += func(x[i + 1], y[j + 2]);
					}
				}
			}
		}

		for (int i = 0; i < m ; i++)
		{
			for (int j = 0; j < m ; j++)
			{
				cout << F[i][j] << " ";
			}
			cout << endl;
		}

	}


};

