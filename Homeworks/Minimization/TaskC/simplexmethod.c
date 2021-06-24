#include<stdio.h>
#include<math.h>
#include<float.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include"utilities.h"

// simplex update
void simplex_update
   (  double** simplex
   ,  double*  fxs
   ,  int      d
   ,  int*     hi
   ,  int*     lo
   ,  double*  centroid
   )
{
   //Re-initializing values
   *hi = 0;
   *lo = 0;
   double fhi = fxs[0];
   double flo = fxs[0];

   double fx;
   for (int i = 1; i < d + 1; ++i) {
      fx = fxs[i];
      if (fx > fhi) {fhi = fx; *hi = i;}
      if (fx > flo) {flo = fx; *lo = i;}
   }

   // calculate centroid
   for (int i = 0; i < d; ++i) {  
      double sum = 0;
      for (int j = 0; j < d + 1; ++j) if (j != *hi) sum += simplex[j][i];
      centroid[i] = sum/d;
   }
}


void simplex_initiate
   (  double   func(double* x)
   ,  double** simplex
   ,  double*  fxs
   ,  int      d
   ,  int*     hi
   ,  int*     lo
   ,  double*  centroid
   )
{
     for (int i = 0; i < d + 1; ++i) fxs[i] = func(simplex[i]);

     simplex_update(simplex, fxs, d, hi, lo, centroid);
}


void reflection(double* fhi, double* centroid, int d, double* reflected) {
   for (int i = 0; i < d; ++i) reflected[i] = 2*centroid[i] - fhi[i];
}

void expansion(double* fhi, double* centroid, int d, double* expanded) {
   for (int i = 0; i < d; ++i) expanded[i] = 3*centroid[i] - 2*fhi[i];
}

void contraction(double* fhi, double* centroid, int d, double* contracted) {
   for (int i = 0; i < d; ++i) contracted[i] = 0.5*centroid[i] + 0.5*fhi[i];
}

void reduction(double** simplex, int d, int lo) {
   for (int j = 0; j < d + 1; ++j) 
      if (j != lo)
         for (int i = 0; i < d; ++i) 
            simplex[j][i] = 0.5*(simplex[j][i] + simplex[lo][i]);
}

double distance(double* a, double* b, int d) {
   double sum = 0;
   for (int i = 0; i < d; ++i) sum += pow(b[i] - a[i], 2);
   return sqrt(sum);
}

double size(double** simplex, int d) { 
   double result = 0;
   double dist;
   for (int i = 0; i < d + 1; ++i) {
      dist = distance(simplex[0], simplex[i], d);
      if (dist > result) result = dist;
   }
   return result;
}

int simplex_driver
   (  double   func(double* x)
   ,  double** simplex
   ,  int      d
   ,  double   sizegoal
   )
{
   int hi = 0;
   int lo = 0;
   int k  = 0;
   double centroid[d];
   double fxs[d + 1];
   double p1[d];
   double p2[d];

   simplex_initiate(func, simplex, fxs, d, &hi, &lo, centroid);

   while (size(simplex, d) > sizegoal) {

      simplex_update(simplex, fxs, d, &hi, &lo, centroid);
      reflection(simplex[hi], centroid, d, p1);
      double f_re = func(p1);

      if (f_re < fxs[lo]) {
         expansion(simplex[hi], centroid, d, p2);
         double f_ex = func(p2);

         if (f_ex < f_re) { 
            for (int i = 0; i < d; ++i) simplex[hi][i] = p2[i]; fxs[hi] = f_ex;
         }

          else { 
            for (int i = 0; i < d; ++i) simplex[hi][i] = p1[i]; fxs[hi] = f_re;
         }
      }

       else {

         if (f_re < fxs[hi]) {
            for (int i = 0; i < d; ++i) simplex[hi][i] = p1[i]; fxs[hi] = f_re;
         } 
         
         else {
            contraction(simplex[hi], centroid, d, p1);
            double f_co = func(p1);

            if (f_co < fxs[hi]) {
               for (int i = 0; i < d; ++i) simplex[hi][i] = p1[i]; fxs[hi] = f_re;
            }

               else {
               reduction(simplex, d, lo); 
               simplex_initiate(func, simplex, fxs, d, &hi, &lo, centroid);
            }
         }

      }     
      k++;

   }
      return k;

}

