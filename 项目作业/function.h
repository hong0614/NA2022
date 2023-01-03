#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include<math.h>
#include<vector>
#include <Eigen/Dense>
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
    
    virtual double ldiff(double x) // calculate left derivative
    {
        return ((*this)(x - _eps) - (*this)(x)) / (-_eps);
    }

    virtual double rdiff(double x) // calculate right derivative
    {
        return ((*this)(x + _eps) - (*this)(x)) / (_eps);
    }

    virtual double seconddiff(double x) //calculate second derivative
    {
        return ((*this)(x + _eps) - 2 * (*this)(x) + (*this)(x - _eps)) / (2 * _eps * _eps);
    }

    virtual double dlddiff(double x) // calculate the left second derivative
    {
        return ((*this)(x - 2 * _eps) + (*this)(x) - 2 * (*this)(x - _eps)) / (_eps * _eps);
    }
    virtual double drdiff(double x) // calculate the right second derivative
    {
        return ((*this)(x + 2 * _eps) + (*this)(x) - 2 * (*this)(x + _eps)) / (_eps * _eps);
    }

};
