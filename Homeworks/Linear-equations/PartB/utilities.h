#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


double dot(gsl_vector* x, gsl_vector* y);

double norm(gsl_vector* x);

//Print of vector and matrix
void vector_print(gsl_vector* vec);
void matrix_print(gsl_matrix* mat);

// binary search for gsl vector and plain array
int binsearch_vector(gsl_vector* x, double x_new);
int binsearch_array(int N, double* x, double x_new);
