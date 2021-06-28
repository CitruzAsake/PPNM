#include "stdio.h"
#include "math.h"

#include "utilities.h"
#include "SVD.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <time.h>


int main() {
	int n = 3;

	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* U = gsl_matrix_alloc(n,n);
	gsl_matrix* D = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);

	for (int i=0;i<n;i++) {
		for (int j=0;j<n;j++) {	gsl_matrix_set(A,i,j,(double)rand()/(double)RAND_MAX*85);
		}
	}

 	gsl_matrix* Adub = gsl_matrix_alloc(n,n);
        gsl_matrix* B1 = gsl_matrix_alloc(n,n);
        gsl_matrix* B2 = gsl_matrix_alloc(n,n);

	// RUNNING THE SVD algorithm and printing results
	gsl_matrix_memcpy(Adub,A);
	printf("Starting matrix A =\n");
	matrix_print(stdout,Adub);

	JSVD(A, V, U, D);

	printf("Final matrix A =\n");
	matrix_print(stdout,A);
	printf("Matrix V =\n");
	matrix_print(stdout,V);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,V,V,0,B1);
  printf("I from V^TV\n");
  matrix_print(stdout,B1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,V,V,0,B1);
	printf("I from VV^T\n");
	matrix_print(stdout,B1);

	printf("Matrix U =\n");
	matrix_print(stdout,U);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,U,U,0,B1);
  printf("I from U^TU\n");
	matrix_print(stdout,B1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,U,U,0,B1);
	printf("I from UU^T\n");
  matrix_print(stdout,B1);
	printf("D =\n");
	matrix_print(stdout,D);
	gsl_matrix_set_identity(B1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1,U,D,0,B1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,B1,V,0,B2);
	printf("A from UDV^T\n");
	matrix_print(stdout,B2);

	// FREE MEMORY

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(Adub);
	gsl_matrix_free(B1);
	gsl_matrix_free(B2);
	gsl_matrix_free(U);
	gsl_matrix_free(D);

	return 0;
}
