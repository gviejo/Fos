#include <math.h>


#define PI 3.141592654


double max(double a, double b)
{
if (a>b)
        return a;
else
        return b;
}


double min(double a, double b)
{
if (a>b)
        return b;
else
        return a;
}
double sgn(double a)
{
if (a>0)
        return 1.;
else if (a<0)
        return -1.;
else
        return 0.;
}
double ndf(double t)
{
return 0.398942280401433*exp(-t*t/2);
}
double nc(double x)
{
double result;
if (x<-7.)
        result = ndf(x)/sqrt(1.+x*x);
else if (x>7.)
        result = 1. - nc(-x);
else
{
result = 0.2316419;
static double a[5] = {0.31938153,-0.356563782,1.781477937,-1.821255978,1.330274429};
result=1./(1+result*fabs(x));
result=1-ndf(x)*(result*(a[0]+result*(a[1]+result*(a[2]+result*(a[3]+result*a[4])))));
if (x<=0.) result=1.-result;
}
return result;
}
double fxy(double x, double y, double a, double b, double rho)
{
    double a_s;
    double b_s;
        double result;
    a_s = a / sqrt(2 * (1 - rho * rho));
    b_s = b / sqrt(2 * (1 - rho * rho));
    result = exp(a_s * (2 * x - a_s) + b_s * (2 * y - b_s) + 2 * rho * (x - a_s) * (y - b_s));
return result;
}
double Ntwo(double a, double b, double rho)
{
    static double aij[4]={0.325303,
                          0.4211071,
                          0.1334425,
                          0.006374323};
        static double bij[4]={0.1337764,
                          0.6243247,
                          1.3425378,
                          2.2626645};
    int i;
    int j;
        double result;
    result = 0;
        for(i=0;i<=3;i++) 
                {
                        for(j=0;j<=3;j++)
                        {
                                result+=aij[i] * aij[j] * fxy(bij[i], bij[j], a, b, rho); 
                        }
                }
    result = result * sqrt(1 - rho * rho) / PI;
return result;
}
double ND2(double a, double b, double rho)
{
    double rho1;
    double rho2;
    double denominator;
        double result;

    if (a * b * rho <= 0) 
            {
        if (a <= 0 && b <= 0 && rho <= 0)
            result = Ntwo(a, b, rho);
        else if (a <= 0 && b * rho >= 0)
            result = nc(a) - Ntwo(a, -b, -rho);
        else if (b <= 0 && rho >= 0)
            result = nc(b) - Ntwo(-a, b, -rho);
        else
            result = nc(a) + nc(b) - 1 + Ntwo(-a, -b, rho);
            }
    else
            {
        denominator = sqrt(a * a - 2 * rho * a * b + b * b);
        rho1 = (rho * a - b) * sgn(a) / denominator;
        rho2 = (rho * b - a) * sgn(b) / denominator;
        result = ND2(a, 0, rho1) + ND2(b, 0, rho2) - (1 - sgn(a) * sgn(b)) / 4;
    }
    if (result < 0) result = 0;
    return result;
}
