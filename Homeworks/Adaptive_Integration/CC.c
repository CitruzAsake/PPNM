//HERE WE WANT TO IMPLEMENT THE CLENSHAW-CURTIS
#include<math.h>
#include<assert.h>
#define SQR2 1.414
#include<stdio.h>
#include"RAI.h"
#include"CC.h"

static double A,B; // private variables caused by "static"

static double F(double f(double),double t){ // auxilliary function for CC
        return f( (A+B)/2+(A-B)/2*cos(t) )*sin(t)*(B-A)/2;
        }

double integralrun(double f(double),double f2, double f3, double a, double b, double acc, double eps, double nrec);

double CC24( double f(double),double a, double b,double acc, double eps, double f2, double f3, int nrec)
{   
    assert(nrec <10000);
    double f1=F(f,a+(b-a)/6), f4=F(f,a+5*(b-a)/6);
	double Z=(2*f1+f2+f3+2*f4)/6*(b-a), q=(f1+f4+f2+f3)/4*(b-a);
	double tolerance=acc+eps*fabs(Z), error=fabs(Z-q)/3;
	if(error < tolerance) return Z;
	else {
		double Z1=CC24(f,a,(a+b)/2,acc/SQR2,eps,f1,f2,nrec+1);
		double Z2=CC24(f,(a+b)/2,b,acc/SQR2,eps,f3,f4,nrec+1);
		return Z1+Z2; }
}

double clenshaw(double f(double),double a,double b, double acc,double eps ){
	A=a; B=b; a=0; b=M_PI;
	double f2=F(f,a+2*(b-a)/6),f3=F(f,a+4*(b-a)/6);
	int nrec=0;
	return CC24(f,a,b,2*acc,2*eps,f2,f3,nrec);
}
