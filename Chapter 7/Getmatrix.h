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
	for (int i = 0; i < (N - 1) * (N - 1); i++)
	{
		for (int j = 0; j < (N - 1) * (N - 1); j++)
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

VectorXd PDE_Matrix::get_eigenF(vector<vector<double>>F, int N)
{
	VectorXd F2D = VectorXd::Zero((N - 1) * (N - 1));
	for (int i = 0; i < N-1; i++)
	{
		for (int j = 0; j < N-1; j++)
		{
			F2D((i)+(j) * (N - 1)) = F[i][j];			
		}
	}

	return F2D;
}

class True_solution
{
public:
	vector<double> true_solution(Function& func, int N)
	{
		double h = 1.0 / N;
		vector<double> u_;
		u_.resize((N - 1) * (N - 1));
		for (int i = 0; i < N - 1; i++)
		{
			for (int j = 0; j < N - 1; j++)
			{
				double t_1 = static_cast<double>(i) + 1;
				double t_2 = static_cast<double>(j) + 1;
				u_[i * (N - 1) + j] = func(t_2 * h, t_1 * h);
			}
		}

		return u_;
	}
};

class Error
{
public:
	Error(VectorXd A, vector<double>B, int norm)
	{
		int m = A.size();
		double error = 0.0;
		if (norm == 1)
		{
			for (int i = 0; i < m; i++)
			{
				error += fabs(A(i) - B[i]);
			}
			cout << "1-norm eroor is " << error << endl;
		}

		else if (norm == 2)
		{
			for (int i = 0; i < m; i++)
			{
				error += (A(i) - B[i]) * (A(i) - B[i]);
			}
			error = sqrt(error);
			cout << "2-norm eroor is " << error << endl;
		}

		else
		{
			for (int i = 0; i < m; i++)
			{
				double e = fabs(A(i) - B[i]);
				if (error < e)
				{
					error = e;
				}
			}

			cout << "max-error is " << error << endl;
		}
	}

	void convergence_rate(double error_n, double error_2n, int n)
	{
		double c = log(error_n / error_2n) / log(2);
		cout << "convergence rate is " << c << endl;
	}
};

class Equationsolver
{
public:
	VectorXd solve(MatrixXd A, VectorXd b)
	{
		MatrixXd x;
		x = A.colPivHouseholderQr().solve(b);
	}
};

class PoissonSolver:public PDE_Matrix
{
private:
	int N;
	vector<vector<double>> A, F;

public:
	//PoissonSolver(Function& func, Function& ddfunc, const string& type_name, int N);

	PoissonSolver(Function& func, Function& ddfunc, const string& type_name, int N)
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
						F[i][j] += func(x[i +2], y[j + 1]);
					}

					if (j == m - 1)
					{
						F[i][j] += func(x[i + 1], y[j + 2]);
					}
				}
			}
		}

		//for (int i = 0; i < m ; i++)
		//{
		//	for (int j = 0; j < m ; j++)
		//	{
		//		cout << F[i][j] << " ";
		//	}
		//	cout << endl;
		//}
		
		//for (int i = 0; i < (N - 1) * (N - 1); i++)
		//{
		//	for (int j = 0; j < (N - 1) * (N - 1); j++)
		//	{
		//		cout << A[i][j] << " ";
		//	}
		//	cout << endl;
		//}

		MatrixXd A2D = get_eigenA(A, N);
		//cout << A2D << endl;

		VectorXd F2D = get_eigenF(F, N);
		//cout << F2D << endl;

		VectorXd u = A2D.colPivHouseholderQr().solve(F2D);
		//cout << u << endl;

		True_solution utrue;
		vector<double> u_true = utrue.true_solution(func,N);
		//for (int i = 0; i < (N - 1) * (N - 1); i++)
		//{
		//	cout << u_true[i] << endl;
		//}

		Error err(u, u_true, 1);
		Error err1(u, u_true, 2);
		Error err2(u, u_true, 3);
	}
};






