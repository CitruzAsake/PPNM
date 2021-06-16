#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include"../util.h"

//
// DECLARATIONS AND DEFINITIONS
//

// least squares routine
int leastsquares(gsl_vector *x, gsl_vector *y, gsl_vector *dy, gsl_vector* c, gsl_vector* dc, gsl_matrix* cov, int m, double (*F)(int, double));


// fitting functionsi (this is for linear fit)
double F(int k, double xi) {
   if      (k == 0) return 1;
   else if (k == 1) return xi;
   fprintf(stderr, "Unknown function index!\n");
   return 111; 
}

// function to evaluate fit using parameters c_k + factor*dc_k
double fit_eval(int m, double z, gsl_vector* c, double factor, gsl_vector* dc) {
   double y = 0;
   double param, delta;
   for (int k = 0; k < m; ++k) {
      param = gsl_vector_get(c, k);
      delta = gsl_vector_get(dc, k);
      y += (param + factor*delta)*F(k, z);
   }
   return y;
}


// 
// MAIN
//

int main() {

   // number of parameters (same as number of fitting functions)
   int m = 2;

   // read data from file
   int n = 9;
   gsl_vector*  x = gsl_vector_alloc(n);
   gsl_vector*  y = gsl_vector_alloc(n);
   gsl_vector* dy = gsl_vector_alloc(n);
   
   FILE* x_file  = fopen("time.txt", "r");
   FILE* y_file  = fopen("activity.txt", "r");
   FILE* dy_file = fopen("unceractivity.txt", "r");

   gsl_vector_fscanf(x_file, x);
   gsl_vector_fscanf(y_file, y);
   gsl_vector_fscanf(dy_file, dy);
   
   //vector_print(x);
   //vector_print(y);
   //vector_print(dy);

   // linearize the data by taking log of y (dyi <- dyi/yi, yi <- log(yi))
   gsl_vector_div(dy, y);
   for (int i = 0; i < y->size; ++i) {
      gsl_vector_set(y, i, log(gsl_vector_get(y, i)));
   }

   // also print linearized data to file
   FILE* logdata_file = fopen("logdata.txt", "w");
   for (int i = 0; i < n; ++i) {
      double xi  = gsl_vector_get(x, i);
      double yi  = gsl_vector_get(y, i);
      double dyi = gsl_vector_get(dy, i);
      fprintf(logdata_file, "%10.5f %10.5f %10.5f\n", xi, yi, dyi);
   }

   //vector_print(x);
   //vector_print(y);
   //vector_print(dy);
   
   // allocate stuff to hold output from least square routine
   gsl_vector* c   = gsl_vector_alloc(m);
   gsl_vector* dc  = gsl_vector_alloc(m);
   gsl_matrix* cov = gsl_matrix_alloc(m, m);

   // do least squares and write parameters with uncertainty to file
   leastsquares(x, y, dy, c, dc, cov, m, &F);
   FILE* param_file = fopen("params.txt", "w");
   for (int i = 0; i < m; ++i) {
      double ci  = gsl_vector_get(c, i);
      double dci = gsl_vector_get(dc, i);
      fprintf(param_file, "%12.6f %12.6f\n", ci, dci);
   }

   // extract lambda (ln(y) = ln(a) - lambda*t) for readability and calculate relative error
   double lambda  = fabs(gsl_vector_get(c, 1));
   double dlambda = gsl_vector_get(dc, 1);
   double eps     = dlambda/lambda;

   // calculate halflife and error (note that the relative error on halflife is the same as on lambda)
   double half  = log(2)/lambda;
   double dhalf = half*eps;
   printf("Half-life (days) from fit = %.4f +/- %.4f = [%.4f, %.4f]\n", half, dhalf, half - dhalf, half + dhalf);


   // write plot file (including upper and lower uncertainty bound)
   FILE* fit_file = fopen("fit.txt", "w");
   double min = 0;
   double max = 16;
   double step = 0.05;
   int N = (max - min)/step;
   double z, ymid, ymin, ymax;
   for (int i = 0; i < N; ++i) {
      z = i*step;
      ymid = fit_eval(m, z, c,  0,  dc);
      ymin = fit_eval(m, z, c, -1,  dc);
      ymax = fit_eval(m, z, c, +1,  dc);
      fprintf(fit_file, "%12.6f %12.6f %12.6f %12.6f\n", z, ymid, ymin, ymax);
   }

   // free memory
   gsl_vector_free(x);
   gsl_vector_free(y);
   gsl_vector_free(dy);
   gsl_vector_free(c);
   gsl_vector_free(dc);
   gsl_matrix_free(cov);

   // close files
   fclose(x_file);
   fclose(y_file);
   fclose(logdata_file);
   fclose(param_file);
   fclose(fit_file);
   
   return 0; 
}

