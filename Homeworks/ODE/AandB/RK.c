#include<stdio.h> 
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include"../utilities.h"

void rkstep23(void (*f)(int n, double x, double* y, double* dydx), int n, double x, double* y_curr, double h, double* y_next, double* err) {
   
   // variables
   double k1[n];
   double k2[n];
   double k3[n];
   double k4[n];
   double y_temp[n];

   // first point
   f(n, x, y_curr, k1);
   for (int i = 0; i < n; ++i) {
      y_temp[i] = y_curr[i] + (1.0/2)*k1[i]*h;
   }

   // second point
   f(n, x + (1.0/2)*h, y_temp, k2);
   for (int i = 0; i < n; ++i) {
      y_temp[i] = y_curr[i] + (3.0/4)*k2[i]*h;
   }

   // third point and second order estimate of y(x + h)
   f(n, x + (3.0/4)*h, y_temp, k3);
   for (int i = 0; i < n; ++i) {
      y_temp[i] = y_curr[i] + ((2.0/9)*k1[i] + (1.0/3)*k2[i] + (4.0/9)*k3[i])*h;
   }

   // fourth point and third order estimate of y(x + h) and error estimate
   f(n, x + h, y_temp, k4);
   for (int i = 0; i < n; ++i) {
      y_next[i] = y_curr[i] + ((7.0/24)*k1[i] + (1./4)*k2[i] + (1.0/3)*k3[i] + (1.0/8)*k4[i])*h;
      err[i] = y_next[i] - y_temp[i];
   }
  

}

// driver function
int odedriver(
	void (*f)(int n, double x, double* y, double* dydx), 
   int n,          
	double  a,  
	double  b,  
	double* ya, 
	double* yb, 
	double  h,  
	double  acc,    // absolute accuracy goal 
	double  eps,    // relative accuracy goal 
   char*   outfile 
) { 

   // declare variables
   int k = 0;     // step counter
   double x;      // current x
   double y[n];   // current y vector
   double yh[n];  // estimate of y(x + h)
   double err[n]; // error estimate vector
   double sum;    // intermediate variable
   double tau;    // local tolerance
   FILE* file = fopen(outfile, "w");

   //first step and write it to file
   x = a;
   fprintf(file, "%20.12e ", x);
   for (int i = 0; i < n; ++i) {
      y[i] = ya[i];
      fprintf(file, "%20.12e ", y[i]);
   }
   fprintf(file, "\n");


   // loop until endpoint is reached
   while (x < b) {
      //printf("===============\n");
      //printf("Step no. %3d\n", k);
      //printf("x = %.4f\ny =\n", x);
      //cvector_print(y, n);

      // avoid overstepping
      if (x + h > b) {
         h = b - x;
      }

      // make step
      rkstep23(f, n, x, y, h, yh, err);
      
      //  local error as norm of error estimate vector
      sum = 0;
      for (int i = 0; i < n; ++i) {
         sum += err[i]*err[i];
      }  
      double e = sqrt(sum);
      //printf("Error = %.5f\n", e);

     
      sum = 0;
      for (int i = 0; i < n; ++i) {
         sum += yh[i]*yh[i];
      }
      double norm = sqrt(sum);
      //printf("Norm = %.5f\n", e);

      // calculate local tolerance
      tau = (norm*eps + acc)*sqrt(h/(b-a));

      // accept or reject step?
      if (e < tau) {
         k++;     // increment counter
         x += h;  // increment x

         // increment elements of y vector and print to file
         fprintf(file, "%20.12e ", x);
         for (int i = 0; i < n; ++i) {
            y[i] = yh[i];
            fprintf(file, "%20.12e ", y[i]);
         }
         fprintf(file, "\n");

      }
      
      // adjust step size
      if (e > 0) { h *= pow(tau/e, 0.25)*0.95; }
      else       { h *= 2; } 
     
   } /* end while */

   // set final y vector
   for (int i = 0; i < n; ++i) {
      y[i] = yh[i];
   }

   //printf("===============\n");
   //printf("Step no. %3d\n", k);
   //printf("x = %.4f\ny =\n", x);
   //cvector_print(y, n);

   // close output file
   fclose(file);

   // return number of steps
   return k;

}

