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
	printf("\n Initial Matrix A =\n");
	matrix_print(stdout,Adub);

	JSVD(A, V, U, D);

	printf("End matrix A =\n");
	matrix_print(stdout,A);
	printf(" V =\n");
	matrix_print(stdout,V);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans,1,V,V,0,B1);
  printf("I from V^TV\n");
  matrix_print(stdout,B1);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans,1,V,V,0,B1);
	printf("I from VV^T\n");
	matrix_print(stdout,B1);

	printf("U =\n");
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

// ----------------TESTING LARGE MATRIX SPEED---------------- //
	FILE* stream2 = fopen("TimeTest.txt","w");
	n = 300;
	clock_t start;
	clock_t end;
	double time_used;

	// Createing matrix and vector needed
	gsl_matrix *M = gsl_matrix_alloc(n,n);
	gsl_matrix *MV = gsl_matrix_alloc(n,n);
	gsl_matrix *MU = gsl_matrix_alloc(n,n);
	gsl_matrix *MD = gsl_matrix_alloc(n,n);

	gsl_matrix *N = gsl_matrix_alloc(n,n);
	gsl_matrix *NV = gsl_matrix_alloc(n,n);
	gsl_vector *NS = gsl_vector_alloc(n);

	gsl_matrix *BUF = gsl_matrix_alloc(n,n);

	// pseudo-Random matrix with numbers n = 0... 100)
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			gsl_matrix_set(M,i,j,(double)rand()/(double)RAND_MAX*100);
		}
	}
	gsl_matrix_memcpy(N,M); // ensure equal starting params

	// TESTING My_SVD
	start = clock();
	JSVD(M, MV, MU, MD);
	end = clock();
	time_used = ((double)(end-start)) / CLOCKS_PER_SEC;
	fprintf(stream2,"The onesided jacobi I have created is a little slower than the GSL. Though I like it very much.\n");
	fprintf(stream2, "My_jacobiSVD algorithm used %g seconds to generate the one sided jacobi-SVD of a %dx%d matrix\n",time_used,n,n);
	// TESTING GSL one-sided Jacobi SVD
	start = clock();
	gsl_linalg_SV_decomp_jacobi(N, NV, NS);
	end = clock();
	time_used = ((double)(end-start)) / CLOCKS_PER_SEC;
	fprintf(stream2, "GSL SVD algorithm used %g seconds to generate the one sided jacobi-SVD of a %dx%d matrix\n",time_used,n,n);
	

	// FREE MEMORY

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(Adub);
	gsl_matrix_free(B1);
	gsl_matrix_free(B2);
	gsl_matrix_free(U);
	gsl_matrix_free(D);
	gsl_matrix_free(M);
	gsl_matrix_free(MV);
	gsl_matrix_free(MU);
	gsl_matrix_free(MD);
	gsl_matrix_free(N);
	gsl_matrix_free(NV);
	gsl_vector_free(NS);
	gsl_matrix_free(BUF);

	return 0;
}
