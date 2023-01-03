#pragma once
#include <iostream> 
#include<limits>
#include <vector>
#include <cmath>
#include <map>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

//const double eps = numeric_limits<double>::epsilon();
//virtual double operator()(double x) = 0;
//
//double diff(double x) //calculate derivative
//{
//    return ((*this)(x + eps) - (*this)(x)) / (eps);
//}


class Bsplines
{

private:
    const double eps = numeric_limits<double>::epsilon();
    int n , idx; // n is n , idx: use to  
    map<int, double> tp;

public:
    Bsplines(int _n, int _idx, map<int, double>_tp)
    {
        this->n = _n;
        this->idx = _idx;
        this->tp = _tp;
    }


    double getval(double x)
    {
        if (n == 0)
        {
            if ((x > tp[idx - 1]) && (x <= tp[idx]))
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
        else
        {
            double c1 = (double)(x - tp[idx - 1]) / (double)(tp[idx + n - 1] - tp[idx - 1]);
            double c2 = (double)(tp[idx + n] - x) / (double)(tp[idx + n] - tp[idx]);
            Bsplines b1(n - 1, idx, tp);
            Bsplines b2(n - 1, idx + 1, tp);
            return (c1 * b1.getval(x) + c2 * b2.getval(x));
        }
    }

    double diff(double x)
    {
        return(getval(x + eps) - getval(x - eps)) / eps;
    }

    double seconddiff(double x)
    {
        return (getval(x + eps) - 2 * getval(x) + getval(x - eps)) / (2 * pow(eps, 3));
    }
 
    double ltdiff(double x)
    {
        return (getval(x) - 3 * getval(x - eps) + 3 * getval(x - 2 * eps) - getval(x - 3 * eps)) / pow(eps,3);
    }

    double rtdiff(double x)
    {
        return (-getval(x) + 3 * getval(x + eps) - 3 * getval(x + 2 * eps) + getval(x + 3 * eps)) / pow(eps, 3);
    }
};

class Bcompletecubicspline :public Bsplines
{
public:
    int N;
    double m1, mN;
    map<int, double> tp, ftp, c;


    virtual double solve() 
    {
        MatrixXd M = MatrixXd::Zero(N + 2, N + 2);
        Bsplines B_0(3, 0, tp);
        M(0, 1) = B_0.diff(ftp[1]);
        Bsplines B_1(3, 1, tp);
        M(0, 2) = B_1.diff(ftp[1]);
        Bsplines B1(3, -1, tp);
        M(0, 0) = B1.diff(ftp[1]);
        Bsplines B_N(3, N, tp);
        M(N + 1, N + 1) = B_N.diff(ftp[N]);        
        Bsplines B_N1(3, N - 1, tp);
        M(N + 1, N) = B_N1.diff(ftp[N]);       
        Bsplines B_N2(3, N - 2, tp);
        M(N + 1, N - 1) = B_N2.diff(ftp[N]);
        int i, j;
        for (int i = 1; i <= N; i++)
        {
            for (int j = i - 1; j <= i + 1; j++)
            {
                Bsplines B_3(3, j - 1, tp);
                M(i, j) = B_3.getval(tp[i]);
            }
        }

        VectorXd b = VectorXd::Zero(N + 2);       
        b(0) = m1; b(N + 1) = mN;        
        for (int i = 1; i <= N; i++) 
        {
            b(i) = ftp[i];
        }

        VectorXd a = VectorXd::Zero(N + 2);
        a = M.colPivHouseholderQr().solve(b);
        int k = -1;
        for (int k; k <= N; k++)
        {
            pair<int, double> c_k = make_pair(k, a(k + 1));
            c.insert(c_k);
        }
    }

    double getval(double x)
    {
        double sum = 0;
        for (int i = -1; i < N; i++)
        {
            Bsplines B_3(3, i, tp);
            sum = sum + c[i] * B_3.getval(x);
        }
    }
};

class Bknotcubicspline :public Bsplines
{
public:
    int N;
    map<int, double> tp, ftp, c;

    virtual double solve()
    {
        MatrixXd M = MatrixXd::Zero(N + 2, N + 2);
        Bsplines B_0(3, 0, tp);
        M(0, 1) = B_0.ltdiff(tp[2]) - B_0.rtdiff(tp[2]);
        Bsplines B_1(3, 1, tp);
        M(0, 2) = B_1.ltdiff(tp[2]) - B_1.rtdiff(tp[2]);
        Bsplines B_2(3, 2, tp);
        M(0, 3) = B_2.ltdiff(tp[2]) - B_2.rtdiff(tp[2]);
        Bsplines B_3(3, 3, tp);
        M(0, 4) = B_3.ltdiff(tp[2]) - B_3.rtdiff(tp[2]);
        Bsplines B1(3, -1, tp);
        M(0, 0) = B1.ltdiff(tp[2]) - B1.rtdiff(tp[2]);
        Bsplines B_N(3, N, tp);
        M(N + 1, N + 1) = B_N.ltdiff(tp[N - 1]) - B_N.rtdiff(tp[N - 1]);
        Bsplines B_N1(3, N - 1, tp);
        M(N + 1, N) = B_N1.ltdiff(tp[N - 1]) - B_N1.rtdiff(tp[N - 1]);
        Bsplines B_N2(3, N - 2, tp);
        M(N + 1, N - 1) = B_N2.ltdiff(tp[N - 1]) - B_N2.rtdiff(tp[N - 1]);
        Bsplines B_N3(3, N - 3, tp);
        M(N + 1, N - 2) = B_N3.ltdiff(tp[N - 1]) - B_N3.rtdiff(tp[N - 1]);
        Bsplines B_N4(3, N - 4, tp);
        M(N + 1, N - 3) = B_N4.ltdiff(tp[N - 1]) - B_N4.rtdiff(tp[N - 1]);
        int i, j;
        for (i = 1; i <= N; i++)
        {
            for (j = i - 1; j <= i + 1; j++)
            {
                Bsplines B_3(3, j - 1, tp);
                M(i, j) = B_3.getval(ftp[i]);
            }
        }
    
        MatrixXd b = VectorXd::Zero(N + 2);
        for (int i = 1; i <= N; i++)
        {
            b(i) = ftp[i];
        }

        VectorXd a = VectorXd::Zero(N + 2);
        a = M.colPivHouseholderQr().solve(b);
        for (int k = -1; k <= N; k++)
        {
            pair<int, double> c_k = make_pair(k, a(k + 1));
            c.insert(c_k);
        }
    } 

    double getval(double x)
    {
        double sum = 0;
        for (int i = -1; i <= N; i++)
        {
            Bsplines B_3(3, i, tp);
            sum = sum + c[i] * B_3.getval(x);
        }
        return sum;
    }
};

class Bcubicsplinewith2nd :public Bsplines
{
public:
    int N;
    double m1, mN;
    map<int, double>tp, ftp, c;

    virtual double solve()
    {
        MatrixXd M = MatrixXd::Zero(N + 2, N + 2);
        Bsplines B_0(3, 0, tp);
        M(0, 1) = B_0.seconddiff(tp[1]);
        Bsplines B_1(3, 1, tp);
        M(0, 2) = B_1.seconddiff(tp[1]);
        Bsplines B1(3, -1, tp);
        M(0, 0) = B1.seconddiff(tp[1]);
        Bsplines B_N(3, N, tp);
        M(N + 1, N + 1) = B_N.seconddiff(tp[N]);
        Bsplines B_N1(3, N - 1, tp);
        M(N + 1, N) = B_N1.seconddiff(tp[N]);
        Bsplines B_N2(3, N - 2, tp);
        M(N + 1, N - 1) = B_N2.seconddiff(tp[N]);
        int i, j;
        for (i = 1; i <= N; i++)
        {
            for (j = i - 1; j <= i + 1; j++)
            {
                Bsplines B_3(3, j - 1, tp);
                M(i, j) = B_3.getval(tp[i]);
            }
        }
        
        VectorXd b = VectorXd::Zero(N + 2);
        b(0) = m1; b(N + 1) = mN;
        for (int i = 1; i <= N; i++)
        {
            b(i) = ftp[i];
        }
        
        VectorXd a = VectorXd::Zero(N + 2);
        a = M.colPivHouseholderQr().solve(b);
        for (int k = -1; k <= N; k++)
        {
            pair<int, double> c_k = make_pair(k, a(k + 1));
            c.insert(c_k);
        }
    }


    double getval(double x)
    {
        double sum = 0;
        for (int i = -1; i <= N; i++)
        {
            Bsplines B_3(3, i, tp);
            sum = sum + c[i] * B_3.getval(x);
        }
        return sum;
    }
};