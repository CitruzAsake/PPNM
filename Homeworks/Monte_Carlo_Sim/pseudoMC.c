#include<stdio.h>
#include<complex.h>
#include<math.h>
#include<assert.h>

double corput(int n, int base) { 
   double q = 0;
   double bk = 1.0/base;
   while (n > 0) { 
      q += (n % base)*bk;
      n  /= base;
      bk /= base;
   }
   return q;
}


double Halton(int type, int n, int dim, int i) {
   if (type == 0) {
      int base[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67};
      int maxdim = sizeof(base)/sizeof(int);
      assert(dim <= maxdim);
      return corput(n, base[i]);
   }
   else {
      int base[] = {71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163};
      int maxdim = sizeof(base)/sizeof(int);
      assert(dim <= maxdim);
      return corput(n, base[i]);
   }
}

void pseudoMC(int dim, double f(int dim, double* x), double* w, double* e, double N, double* result, double* error)
{

double Q=1;
for (int i = 0; i < dim; ++i) Q*= e[i] - w[i];

double total;
double x[dim];
double yx;

total =0;
for (int i = 0; i < N; ++i) {
      for(int j = 0; j < dim; ++j) x[j] = w[j] + Halton(0, i, dim, j)*(e[j] - w[j]);
     yx = f(dim, x);
     total += yx;
   }
   double mean_w = total/N;

   total = 0;
    for (int i = 0; i < N; ++i) {
      for(int j = 0; j < dim; ++j) x[j] = w[j] + Halton(1, i, dim, j)*(e[j] - w[j]);
      yx = f(dim, x);
      total += yx;
   }
   double mean_e = total/N;

   *result = Q*(mean_w + mean_e)/2;
   *error  = fabs(Q*mean_w - Q*mean_e);
}







