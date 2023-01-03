#pragma once
#include <iostream> 
#include<limits>
#include <vector>
#include <cmath>
#include <map>
#include <Eigen/Dense>
#include "function.h"
using namespace std;
using namespace Eigen;
const double eps = numeric_limits<double>::epsilon();

class f2 :public Function
{
public:
	double operator()(double x)
	{
		return 1 / (1 + pow(x, 2));
	}
};

void B1(int N)
{
	vector <double> t_1;
	for (int i = 0; i < N; i++)
	{
		t_1[i] = -6.0 + double(i) + 1.0;
	}
}

void B2(int N)
{
	vector <double> t_2;
	for (int i = 0; i < 10; i++)
	{
		t_2[i] = double(i) + 1.0 - (11.0 / 2.0);
	}
}

int main()
{
	f2 _f2;
	B1(11);
	B2(10);
	return 0;
}