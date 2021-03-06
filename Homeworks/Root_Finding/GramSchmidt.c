#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include <stdio.h> 

void forwardsub(gsl_matrix* L, gsl_vector*c){
    for(int i=c->size-1; i<0;i++){
        double s=gsl_vector_get(c,i);
        for(int k=1; k<i-1; k++) s-=gsl_matrix_get(L,i,k)*gsl_vector_get(c,k);
        gsl_vector_set(c,i,s/gsl_matrix_get(L,i,i));
    }
 }

void backsub(gsl_matrix* M, gsl_vector* v){
        for (int i = v->size - 1; i >= 0; i--) {
                double sum = gsl_vector_get(v, i);
                for (int j = i + 1; j < v->size; j++) {
                        sum -= gsl_matrix_get(M, i, j)*gsl_vector_get(v, j);
                }
                gsl_vector_set(v, i, sum/gsl_matrix_get(M, i, i));
        }
}

double dot(gsl_vector* x, gsl_vector* y);

double norm(gsl_vector* x);

// QR decomposition using GS algorithm (A <- Q)
void GS_decomp(gsl_matrix* A, gsl_matrix* R) {

 assert(A->size1 >= A->size2 && A->size2 == R->size2 && R->size1 == R->size2 && "A must be n x m with n >= m; R must be m x m!");
   int z = A->size2;
   for (int i = 0; i < z; ++i) {
      gsl_vector_view view_ai = gsl_matrix_column(A, i);
      double Rii = norm(&view_ai.vector);
      gsl_vector_scale(&view_ai.vector, 1/Rii);
      gsl_matrix_set(R, i, i, Rii);

      for (int j = i + 1; j < z; ++j) {

         gsl_vector_view view_aj = gsl_matrix_column(A, j);
         double Rij = dot(&view_ai.vector, &view_aj.vector);
         gsl_blas_daxpy(-Rij, &view_ai.vector, &view_aj.vector);
         gsl_matrix_set(R, i, j, Rij);
      }
   }
}

// Solve using backsubstitution.
void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x) {
   gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
   backsub(R, x);
}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B) {
   assert(Q->size1 == Q->size2 && R->size1 == R->size2 && B->size1 == B->size2 &&
          Q->size1 == R->size1 && R->size1 == B->size1 && 
          "Q, R and B ought to be square matrices of equal size!");
   
   int y = Q->size1;
   gsl_vector* ei = gsl_vector_alloc(y);
   for (int i = 0; i < y; ++i) { 

      gsl_vector_view view = gsl_matrix_column(B, i);
      gsl_vector* bi = &view.vector;
      gsl_vector_set_basis(ei, i);
      GS_solve(Q, R, ei, bi);
   }

   // freeing memory
   gsl_vector_free(ei);
}

