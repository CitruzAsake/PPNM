#include<math.h>
#include<assert.h>
#include<stdio.h>
#define SQR2 1.414

// We want to implement a recursive adaptive integrator that simply calculates an estimate of a given function f(x) in some interval with absolute(\delta) or relativ(r) accurary goals.

double integralrun(double f(double),double f2, double f3, double a, double b, double acc, double eps, double nrec)
{
    assert(nrec <10000);
    double f1 = f(a+1.*(b-a)/6.);
    double f4 = f(a+5.*(b-a)/6.);
    double Z = (b-a)*(2.*f1+f2+f3+2*f4)/6.;
    double z = (b-a)*(f1+f4+f2+f3)/4.;
    double toll = acc+eps*fabs(Z);
    double err = fabs(Z-z);

    if(err < toll) return Z;
    else
    {
        double Z1 = integralrun(f,f1,f2,a,(a+b)/2,acc/sqrt(2),eps,nrec+1);
        double Z2 = integralrun(f,f3,f4,(a+b)/2,b,acc/sqrt(2),eps,nrec+1);
        return Z1+Z2;
    }
}

double integrate(double f(double), double a, double b, double acc, double eps)
    {
        double f2 = f(a+2.*(b-a)/6.);
        double f3 = f(a+4.*(b-a)/6.);
        int nrec = 0;
        double value = integralrun(f,f2,f3,a,b,acc,eps,nrec);

        return value;
    }
