#include<math.h>
#include<assert.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include"utilities.h"

// vector dot product
double dot(gsl_vector* x, gsl_vector* y) {
   double result;
   gsl_blas_ddot(x, y, &result);
   return result;
}


// vector norm
double norm(gsl_vector* x) {
   return sqrt(dot(x, x));
}


// in-place back-substitution for U*x = Qt*b (takes U and x = Qt*b and puts solution in place of x)
void backsub(gsl_matrix* R, gsl_vector* x) {
   for (int i = x->size - 1; i >= 0; --i) {
      double xi = gsl_vector_get(x, i);
      for (int k = i + 1; k < x->size; ++k) {
         xi -= gsl_matrix_get(R, i, k)*gsl_vector_get(x, k);
      }
      xi /= gsl_matrix_get(R, i, i);
      gsl_vector_set(x, i, xi);
   }
}


// QR decomposition by Gram-Schmidt algorithm (A <- Q)
void GS_decomp(gsl_matrix* A, gsl_matrix* R) {

   assert(A->size1 >= A->size2 && A->size2 == R->size2 && R->size1 == R->size2 && "A must be n x m with n >= m; R must be m x m!");
   int m = A->size2;

   // calculate Q (in place of A) and R
   for (int i = 0; i < m; ++i) {

      // calculate ai <- ai/||ai||
      gsl_vector_view view_ai = gsl_matrix_column(A, i);
      gsl_vector* ai = &view_ai.vector;
      double Rii = norm(ai);
      gsl_vector_scale(ai, 1/Rii);
      gsl_matrix_set(R, i, i, Rii);

      // calculate aj <- aj - <qi|aj>qi
      for (int j = i + 1; j < m; ++j) {
         gsl_vector_view view_aj = gsl_matrix_column(A, j);
         gsl_vector* aj = &view_aj.vector;
         double Rij = dot(ai, aj);
         gsl_blas_daxpy(-Rij, ai, aj);
         gsl_matrix_set(R, i, j, Rij);
      }
   }
}


// Solve QRx = b <=> R*x = Qt*b by back-substitution
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x) {

   // calculate x <- Qt*b
   gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);

   // do back-substitution
   backsub(R, x);

}


// calculate inverse of A = Q*R into matrix B
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B) {
   
   // check sizes
   assert(Q->size1 == Q->size2 && R->size1 == R->size2 && B->size1 == B->size2 &&
          Q->size1 == R->size1 && R->size1 == B->size1 && 
          "Q, R and B must be square matrices of same size!");
   
   int n = Q->size1;

   // allocate vector ei for solving A*bi = ei
   gsl_vector* ei = gsl_vector_alloc(n);

   // loop over columns i
   for (int i = 0; i < n; ++i) { 

      gsl_vector_view view = gsl_matrix_column(B, i);
      gsl_vector* bi = &view.vector;
      gsl_vector_set_basis(ei, i);
      
      GS_solve(Q, R, ei, bi);
   }

   // free memory
   gsl_vector_free(ei);

}

