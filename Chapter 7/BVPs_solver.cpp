#pragma once
#include"Header.h"
#include"Getmatrix.h"
#include<iostream>
#include<cmath>
#include<limits>
#include<vector>
using namespace std;
    
int main()
{

    std::cout << "Hello World!\n";
    f_1 _f1;
    u_1 _u1;
    cout << "when n=8 " << endl;
    PoissonSolver(_u1, _f1, "Dirichlet", 8);
    cout << endl;
    cout << "when n=16 " << endl;
    PoissonSolver(_u1, _f1, "Dirichlet", 16);
    cout << endl;
    cout << "when n=32 " << endl;
    PoissonSolver(_u1, _f1, "Dirichlet", 32);
    cout << endl;
    cout << "when n=64 " << endl;
    PoissonSolver(_u1, _f1, "Dirichlet", 64);
    cout << endl;
    return 0;
};

// run
    //int N = 4;
    //PDE_Matrix p1;
    //p1.getA(N);
    //f_1 _f1;
    //p1.getF(_f1, 4);

   // PDE_Matrix p1;
    //vector<vector<double>>A = p1.getA(N);
    //vector<vector<double>>F = p1.getF(_f1, N);
    //for (int i = 0; i < (N - 1) * (N - 1); i++)
    //{
       // for (int j = 0; j < (N - 1) * (N - 1); j++)
       // {
          //  cout << A[i][j] << " ";
       // }
       // cout << endl;
    //}
    //for (int i = 0; i < N - 1; i++)
    //{
    //	for (int j = 0; j < N - 1; j++)
    //	{
    //		cout << F[i][j] << " ";
    //	}
    //	cout << endl;
    //}
