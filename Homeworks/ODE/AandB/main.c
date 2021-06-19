#include<stdio.h> 
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include"../utilities.h"

// stepper function
void rkstep23(void (*f)(int n, double x, double* y, double* dydx), int n, double* yx, double h, double* yh, double* err);

// driver function
void odedriver(
	void (*f)(int n, double x, double* y, double* dydx), // right-hand-side of dy/dx = f(x, y) 
   int n,          // size of vectors 
	double  a,      // the start-point a 
	double  b,      // the end-point of the integration 
	double* ya,     // y(a) 
	double* yb,     // y(b) to be calculated 
	double  h,      // initial step-size 
	double  acc,    // absolute accuracy goal 
	double  eps,    // relative accuracy goal 
	char*   outfile // trajectory file
);

// RHS function for u'' = -u
void f(int n, double x, double* y, double* dydx) {
   dydx[0] =  y[1];
   dydx[1] = -y[0];
}

// RHS function for SIR model
double N;
double Tr;
double Tc;
void f_SIR(int n, double x, double* y, double*dydx) {
   dydx[0] = -y[0]*y[1]/(N*Tc);
   dydx[1] =  y[0]*y[1]/(N*Tc) - y[1]/Tr;
   dydx[2] =  y[1]/Tr;
}

int main() {
   {
   //
   // TEST: u'' = -u
   //
   
   // solve u'' = -u as a test (solution is a cosine)
   // the system is equivalent to
   // [ y0' ] = [  y1 ]
   // [ y1' ]   [ -y0 ]
   // with  y0 = u and y1 = u'

   // declare and set variables for odedriver
   int n = 2; 
   double a = 0;
   double b = 2*M_PI;
   double ya[n];
   double yb[n];
   double h = 0.001;
   double acc = 1e-4;
   double eps = 1e-4;
   char* outfile = "cosine.txt";

   // set initial values
   ya[0] = 1;
   ya[1] = 0;

   // call driver
   odedriver(&f, n, a, b, ya, yb, h, acc, eps, outfile);
   }

   {
   // 
   // TEST: SIR model
   //

   // declare and set variables for odedriver
   int n = 3;
   int a = 0;
   int b = 60;
   double ya[n];
   double yb[n];
   double h = 0.1;
   double acc = 1;
   double eps = 1e-6;
   char* outfile;
  
   // SCENARIO 1

   // set parameters for the SIR model
   outfile = "sir1.txt";
   N = 5.8e6;
   Tr = 7;
   Tc = 1.0;

   // set initial values
   ya[0] = N;
   ya[1] = 50;
   ya[2] = 0;
   
   // call driver
   odedriver(&f_SIR, n, a, b, ya, yb, h, acc, eps, outfile);


   // SCENARIO 2

   // set parameters for the SIR model
   outfile = "sir2.txt";
   N = 5.8e6;
   Tr = 7;
   Tc = 1.5;

   // set initial values
   ya[0] = N;
   ya[1] = 50;
   ya[2] = 0;
   
   // call driver
   odedriver(&f_SIR, n, a, b, ya, yb, h, acc, eps, outfile);



   // SCENARIO 3

   // set parameters for the SIR model
   outfile = "sir3.txt";
   N = 5.8e6;
   Tr = 7;
   Tc = 2.5;

   // set initial values
   ya[0] = N;
   ya[1] = 50;
   ya[2] = 0;
   
   // call driver
   odedriver(&f_SIR, n, a, b, ya, yb, h, acc, eps, outfile);


   // SCENARIO 4

   // set parameters for the SIR model
   outfile = "sir4.txt";
   b = 180;
   N = 5.8e6;
   Tr = 7;
   Tc = 4.0;

   // set initial values
   ya[0] = N;
   ya[1] = 50;
   ya[2] = 0;
   
   // call driver
   odedriver(&f_SIR, n, a, b, ya, yb, h, acc, eps, outfile);

   }
   

   return 0;
}

