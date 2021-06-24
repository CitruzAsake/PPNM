#include<stdio.h>
#include<math.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"utilities.h"

int simplex_driver
   (  double   func(double* x)
   ,  double** simplex
   ,  int      d
   ,  double   sizegoal
   );

void simplex_initiate
   (  double   func(double* x)
   ,  double** simplex
   ,  double*  fxs
   ,  int      d
   ,  int*     hi
   ,  int*     lo
   ,  double*  centroid
   );


double rosenbrock(double* x) {
   return pow(1 - x[0], 2) + 100*pow(x[1] - x[0]*x[0], 2);
}


double himmelblau(double* x) {
   return pow(x[0]*x[0] + x[1] - 11, 2) + pow(x[0] + x[1]*x[1] - 7, 2);
}

int main() {
   

   // ROSENBROCK
   {
   printf("Searching for minimum of Rosenbrock function f(x,y) = (1 - x)**2 + 100*(y - x**2)**2\n");
   int d = 2;
   double* simplex[] = {
      (double[]){0.0, 0.0},
      (double[]){0.0, 5.0},
      (double[]){5.0, 0.0}
   };
   double sizegoal = 1e-5;
   double fxs[d + 1];
   double centroid[d];
   int hi = 0;
   int lo = 0;

   // call solver
   int iter = simplex_driver(rosenbrock, simplex, d, sizegoal);

   // print result
   printf("Converged after %d iterations!\n", iter);
   simplex_initiate(rosenbrock, simplex, fxs, d, &hi, &lo, centroid);
   printf("Minimum value y = %.4e located at x = \n", fxs[lo]);
   cvector_print(simplex[lo], d);
   } 

   // HIMMELBLAU
   {
   // set up problem
   printf("-----------------------------------------\n");
   printf("Searching for minimum of Himmelblau function g(x,y) = (x*x + y - 11)**2 + (x + y*y - 7)**2\n");
   int d = 2;
   double* simplex[] = {
      (double[]){0.0, 0.0},
      (double[]){0.0, 5.0},
      (double[]){5.0, 0.0}
   };
   double sizegoal = 1e-5;
   double fxs[d + 1];
   double centroid[d];
   int hi = 0;
   int lo = 0;

   // call solver
   int iter = simplex_driver(himmelblau, simplex, d, sizegoal);

   // print result
   printf("Converged after %d iterations!\n", iter);
   simplex_initiate(himmelblau, simplex, fxs, d, &hi, &lo, centroid);
   printf("Minimum value y = %.4e located at x = \n", fxs[lo]);
   cvector_print(simplex[lo], d);
   } 


   return 0;
}

