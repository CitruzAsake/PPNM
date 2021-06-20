#include<stdio.h>
#include<stdlib.h>
#include <time.h>
#include"utilities.h"
#include"GramSchmidt.h"
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

//GS decomp part.
void GS_decomp(gsl_matrix* A, gsl_matrix* R);

// back substitution part. 
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

// second part (b) inverse  
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

void GS_inverse(gsl_matrix* Q,gsl_matrix* R, gsl_matrix* B);

int main() { 
   
   // Printing relevant task.
   printf("--------------------A---------------\n");

   // random gsl matrix A, where N has to be larger than M (N>M).
   int N = 5;
   int M = 3;
   unsigned int SEED = time(NULL);
   srandom(SEED);
   gsl_matrix* A = gsl_matrix_alloc(N, M);
   for (int i = 0; i < A->size1; ++i) {
      for (int j = 0; j < A->size2; ++j) {
         gsl_matrix_set(A, i, j, 5*((double)random())/RAND_MAX);
      }
   }
   // Generation of arbitary matrix A used as original matrix.
   printf("Original A matrix:\n");
   matrix_print(A);

   // Empty gsl matrix R
   gsl_matrix* R = gsl_matrix_alloc(M, M);

   // Do QR 
   GS_decomp(A, R);

   printf("Q:\n");
   matrix_print(A);

   gsl_matrix* QtQ = gsl_matrix_alloc(M, M);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, QtQ);
   printf("Qt*Q:\n");
   matrix_print(QtQ);
   
   printf("Upper triangular part: R:\n");
   matrix_print(R);

   gsl_matrix* QR = gsl_matrix_alloc(N, M);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);
   printf(" Q*R, equivalent to A:\n");
   matrix_print(QR);

   // free memory  
   gsl_matrix_free(A);
   gsl_matrix_free(R);
   gsl_matrix_free(QtQ);
   gsl_matrix_free(QR);
 
return 0;
}


