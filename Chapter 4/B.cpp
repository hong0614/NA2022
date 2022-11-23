#include<iostream>
#include<cmath>
#include<vector>
#include<algorithm>
using namespace std;

int main()
{
    char beta = 2, p = 3, L = -1, U = 1;
    double UFL = pow(beta, L);
    double OFL = pow(beta, U) * (beta - pow(beta, 1 - p));
    cout << "UFL is " << UFL << " and OFL is " << OFL << endl;

    vector<double> A;
    cout << "Normal numbers in F:" << endl;
    for (int i = L; i <= U; i++)
    {
        for (int j = pow(2, p - 1); j <= pow(2, p) - 1; j++)
        {
            double num = j;
            for (int k = 1; k <= abs(p - 1 - i); k++)
            {
                if (p - 1 - i > 0)
                {
                    num /= 2;
                }
                else
                {
                    num *= 2;
                }
            }
            A.push_back(num);
            A.push_back(-num);
        }
    }
    A.push_back(0);
    sort(A.begin(), A.end());
    cout << "x_1=[";
    for (double i : A)
    {
        cout << i << " , ";
    }
    cout << "]';" << endl;

    vector<double> B;
    cout << "All the subnormal numbers:" << endl;
    for (int i = 1; i <= pow(2, p - 1) - 1; i++)
    {
        double num = i;
        for (int j = 1; j <= p + abs(L); j++)
        {
            num /= 2;
        }
        B.push_back(num);
        B.push_back(-num);
    }
    sort(B.begin(), B.end());
    cout << "x_2=[";
    for (double i : B)
    {
        cout << i << ", ";
    }
    cout << "]';" << endl;

    return 0;
}
