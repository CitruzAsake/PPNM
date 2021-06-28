#include "stdio.h"
#include "math.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include<assert.h>
#include <gsl/gsl_matrix.h>
#include "utilities.h"

double dot(gsl_vector* x, gsl_vector* y){
        double x_dot_y;
        gsl_blas_ddot(x, y, &x_dot_y); //GSL function that calculates dot product
        return x_dot_y;
}

double vector_len(gsl_vector* x){
 return sqrt(dot(x,x));
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

int matrix_print(FILE* stream, const gsl_matrix *X){
  int n = X->size1;
  int m = X->size2;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      fprintf(stream,"%0.2f ",gsl_matrix_get(X,i,j));
    }
    fprintf(stream,"\n");
  }
  fprintf(stream,"\n\n");
}



double vector_print(FILE* stream, gsl_vector* vec){
  int n = vec->size;
  for (size_t i = 0; i < n; i++) {
    fprintf(stream,"%g ",gsl_vector_get(vec,i));
  }
  fprintf(stream,"\n\n");
}
void GS_decomp(gsl_matrix* A, gsl_matrix* R){
   assert(A->size1 >= A->size2 && A->size2 == R->size2 && R->size1 == R->size2 && "A must be n x m with n >= m; R must be m x m!");
   int z = A->size2;
   for (int i = 0; i < z; ++i) {
      gsl_vector_view view_ai = gsl_matrix_column(A, i);
      double Rii = vector_len(&view_ai.vector);
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

void GS_solve(gsl_matrix *Q, gsl_matrix *R, gsl_vector *b, gsl_vector *x) {
   gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x);
   backsub(R, x);
}

void GS_inverse(gsl_matrix *A, gsl_matrix *B){
gsl_matrix *R = gsl_matrix_alloc(A->size2,A->size2);
  gsl_matrix *Acopy = gsl_matrix_alloc(A->size1,A->size2);
  gsl_matrix_memcpy(Acopy,A);
  GS_decomp(Acopy, R);

  int n = B->size1;

  gsl_matrix_set_identity(B);
  gsl_vector *e = gsl_vector_alloc(n);
  gsl_vector *buff = gsl_vector_alloc(n);

        for(int i = 0; i < n; i++){
                gsl_matrix_get_col(e, B, i);
                GS_solve(Acopy,R,e,buff);
                gsl_matrix_set_col(B,i,buff);
        }
  gsl_vector_free(e);  gsl_vector_free(buff);
  gsl_matrix_free(R);  gsl_matrix_free(Acopy);
}

