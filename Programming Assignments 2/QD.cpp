#pragma once
#include<iostream>
#include<cmath>
#include<vector>
#include"Interpolation.h"
using namespace std;

int main()
{
    vector<double> x = { 0,0,3,3,5,5,8,8,13,13 };
    vector<double> y = { 0,0,225,225,383,383,623,623,993,993 };
    vector<double> _x = { 75,75,77,77,80,80,74,74,72,72 };

    cout << " Hermite Polynomial : ";
    HermitePoly(x,y,_x);

    return 0;
}