#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

// QR decomposition by Gram-Schmidt algorithm (A <- Q)
void GS_decomp(gsl_matrix* A, gsl_matrix* R);

// Solve QRx = b <=> R*x = Qt*b by back-substitution
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);

// calculate inverse of A = Q*R into matrix B
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B);

// least squares routine
int leastsquares(gsl_vector *x, gsl_vector *y, gsl_vector *dy, gsl_vector* c, gsl_vector* dc, gsl_matrix* cov, int m, double (*F)(int, double)) {
   
   // number of data points 
   int n = x->size;

   // sanity checks
   assert(y->size == n && dy->size == n && "Vector x, y and dy must have same size!");
   assert(c->size == m && dc->size == m && "Vectors c and dc must have size m!");
   assert(cov->size1 == m && cov->size2 == m && "Matrix cov must be m x m!");

   // allocate 
   gsl_vector* b = gsl_vector_alloc(n);
   gsl_matrix* A = gsl_matrix_alloc(n, m);
   gsl_matrix* R = gsl_matrix_alloc(m, m);
   gsl_matrix* Rinv = gsl_matrix_alloc(m, m);

   // fill out matrix A and vector b
   double xi, yi, dyi;
   for (int i = 0; i < n; ++i) {
      xi  = gsl_vector_get(x, i);
      yi  = gsl_vector_get(y, i);
      dyi = gsl_vector_get(dy, i);
      gsl_vector_set(b, i, yi/dyi);
      for (int k = 0; k < m; ++k) {
         gsl_matrix_set(A, i, k, F(k, xi)/dyi);
      }
   }

   // do QR decomposition
   GS_decomp(A, R);

   // do back-subsitution to get parameters
   GS_solve(A, R, b, c);

   // invert R
   gsl_matrix_set_identity(cov);
   GS_inverse(cov, R, Rinv);

   // calculate covariance matrix
   gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, Rinv, Rinv, 0, cov);

   // calculate parameter uncertainties
   double covkk;
   for (int k = 0; k < m; ++k) {
      covkk = gsl_matrix_get(cov, k, k);
      gsl_vector_set(dc, k, sqrt(covkk));
   }

   return 0;
}

