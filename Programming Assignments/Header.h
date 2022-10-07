#pragma once
#include<iostream>
#include<cmath>
#include<limits>
using namespace std;

class Function
{
    const double _eps = 10 * numeric_limits<double>::epsilon();
public:
    virtual double operator()(double x) = 0;

    virtual double diff(double x) //calculate derivative
    {
        return ((*this)(x + _eps) - (*this)(x))  /( _eps);
    }
};

class EquationSolver
{
public:
    virtual double solve()
    {
        return 0;
    }
};

class Bisection : public EquationSolver
{
private:
    double a, b, delta, epsilon;
    Function & f;
    short M;

public:
    Bisection(Function& _f, double _a, double _b, double _delta, double _epsilon, short _M) :f(_f)
    {
        this->a = _a;
        this->b = _b;
        this->delta = _delta;
        this->epsilon = _epsilon;
        this->M = _M;
    }

    ~Bisection(){};

    double solve()
    {
        if (f(a) * f(b) > 0)
        {
            cout << " No real solution. " << endl;
            return 0.f;
        }
        
        double h = b - a, u = f(a), c=0, w=0;
        for (int k = 1; k <= M; k++)
        {
            h = h / 2;
            c = a + h;
            w = f(c);
            if (fabs(h) < delta || fabs(w) < epsilon)
                break;
            else if ((w > 0 && u > 0) || (w < 0 && u < 0))
            {
                a = c;
            }
        }

        return c;
    }
  
};

class Newton : public EquationSolver
{
private:
    double x0 , epsilon;
    Function& f;
    short M;

public:
    Newton(Function& _f, double x_0, double _epsilon, short _M) :f(_f)
    {
        this->x0 = x_0;
        this->epsilon = _epsilon;
        this->M = _M;
    }

    ~Newton() {};

    double solve()
    {
        double x = x0;
        double u;
        for (int k = 0; k <= M; k++)
        {
            u = f(x);
            if (fabs(u) < epsilon)
                break;
            else
            {
                x = x - u / f.diff(x);
            }              
        }

        return x;
    }
};

class Secant : public EquationSolver
{
private:
    double x0, x1, epsilon, delta;
    Function& f;
    short M;

public:
    Secant(Function& _f, double x_0, double x_1, double _delta, double _epsilon, short _M) :f(_f)
    {
        this->x0 = x_0;
        this->x1 = x_1;
        this->delta = _delta;
        this->epsilon = _epsilon;
        this->M = _M;
    }

    ~Secant() {};

    double solve()
    {
        double x_n = x1, x_m = x0;
        double u = f(x_n), v = f(x_m);
        for (int k = 2; k <= M; k++)
        {
            if (fabs(u) > fabs(v))
            {
                double t = x_n;
                x_n = x_m;
                x_m = t;

                t = u;
                u = v;
                v = t;
            }
            double s = (x_n - x_m) / (u - v);
            x_m = x_n;
            v = u;
            x_n = x_n - u * s;
            u = f(x_n);
            if (fabs(x_n - x_m) < delta || fabs(u) < epsilon)
                break;
        }

        return x_n;
    }

};
