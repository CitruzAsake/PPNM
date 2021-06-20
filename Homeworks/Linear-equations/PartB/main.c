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
   printf("--------------------B---------------\n");

   // initialize random square matrix A and empty matrix R
   int N = 4;   
   unsigned int SEED = time(NULL);
   srandom(SEED);
   gsl_matrix* A = gsl_matrix_alloc(N, N);
   for (int i = 0; i < A->size1; ++i) {
      for (int j = 0; j < A->size2; ++j) {
         gsl_matrix_set(A, i, j, 5*((double)random())/RAND_MAX);
      }
   }

   // print original A matrix
   printf("This is the original A matrix:\n");
   matrix_print(A);

   // do QR decomposition and check result
   gsl_matrix* R = gsl_matrix_alloc(N, N);
   GS_decomp(A, R);

   gsl_matrix* QR = gsl_matrix_alloc(N, N);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);
   printf("Q*R=A :\n");
   matrix_print(QR);

   // Empty matrix B and determine inverse
   gsl_matrix* B = gsl_matrix_alloc(N, N);
   GS_inverse(A, R, B);
   printf("B:\n");
   matrix_print(B);
   gsl_matrix* AB = gsl_matrix_alloc(N, N);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, QR, B, 0, AB);
   printf("A*B:\n");
   matrix_print(AB);
   gsl_matrix* BA = gsl_matrix_alloc(N, N);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, QR, 0, BA);
   printf("B*A:\n");
   matrix_print(BA);  

   // free memory
   gsl_matrix_free(A);
   gsl_matrix_free(R);
   gsl_matrix_free(QR);
   gsl_matrix_free(B);
return 0;
}


