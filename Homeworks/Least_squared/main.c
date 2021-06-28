#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include"utilities.h"
int leastsquares(gsl_vector *x, gsl_vector *y, gsl_vector *dy, gsl_vector* c, gsl_vector* dc, gsl_matrix* cov, int m, double (*Fitfunc)(int, double));

double Fitfunc(int k, double xi) {
   if      (k == 0) return 1;
   else if (k == 1) return xi;
   fprintf(stderr, "Unknown function index!\n");
   return 123; 
}

double eval(int m, double z, gsl_vector* c, double factor, gsl_vector* dc) {
   double y = 0; double param, delta;
   for (int k = 0; k < m; ++k) {
      param = gsl_vector_get(c, k);
      delta = gsl_vector_get(dc, k);
      y += (param + factor*delta)*Fitfunc(k, z);
   }
   return y;
}

int main() {

   int m = 2; int n = 9;
   gsl_vector*  x = gsl_vector_alloc(n);
   gsl_vector*  y = gsl_vector_alloc(n);
   gsl_vector* dy = gsl_vector_alloc(n);
   FILE* x_file  = fopen("time.txt", "r");
   FILE* y_file  = fopen("activity.txt", "r");
   FILE* dy_file = fopen("unceractivity.txt", "r");
   gsl_vector_fscanf(x_file, x);
   gsl_vector_fscanf(y_file, y);
   gsl_vector_fscanf(dy_file, dy);
   gsl_vector_div(dy, y);
   for (int i = 0; i < y->size; ++i) {
      gsl_vector_set(y, i, log(gsl_vector_get(y, i)));
   }
   FILE* logdata_file = fopen("logdata.txt", "w");
   for (int i = 0; i < n; ++i) {
      double xi  = gsl_vector_get(x, i);
      double yi  = gsl_vector_get(y, i);
      double dyi = gsl_vector_get(dy, i);
      fprintf(logdata_file, "%10.5f %10.5f %10.5f\n", xi, yi, dyi);
   }
   gsl_vector* c   = gsl_vector_alloc(m);
   gsl_vector* dc  = gsl_vector_alloc(m);
   gsl_matrix* cov = gsl_matrix_alloc(m, m);

   leastsquares(x, y, dy, c, dc, cov, m, &Fitfunc);
   FILE* param_file = fopen("params.txt", "w");
   for (int i = 0; i < m; ++i) {
      double ci  = gsl_vector_get(c, i);
      double dci = gsl_vector_get(dc, i);
      fprintf(param_file, "%12.6f %12.6f\n", ci, dci);
   }

   double lambda  = fabs(gsl_vector_get(c, 1));
   double dlambda = gsl_vector_get(dc, 1);
   double eps     = dlambda/lambda;
   double half  = log(2)/lambda;
   double dhalf = half*eps;
   printf("Half-life (days) from fit = %.4f +/- %.4f = [%.4f, %.4f]\n", half, dhalf, half - dhalf, half + dhalf);

   FILE* fit_file = fopen("fit.txt", "w");
   double min = 0;double max = 25; double step = 0.08;
   int N = (max - min)/step;
   double z, ymid, ymin, ymax;
   for (int i = 0; i < N; ++i) {
      z = i*step;
      ymid = eval(m, z, c,  0,  dc);
      ymin = eval(m, z, c, -1,  dc);
      ymax = eval(m, z, c, +1,  dc);
      fprintf(fit_file, "%12.6f %12.6f %12.6f %12.6f\n", z, ymid, ymin, ymax);
   }

   gsl_vector_free(x);   gsl_vector_free(y);
   gsl_vector_free(dy);   gsl_vector_free(c);
   gsl_vector_free(dc);   gsl_matrix_free(cov);
 
   fclose(x_file);   fclose(y_file);
   fclose(logdata_file);   fclose(param_file);
   fclose(fit_file);
   return 0; 
}
