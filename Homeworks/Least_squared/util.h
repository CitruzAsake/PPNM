// print gsl vector
void vector_print(gsl_vector* vec);

// print gsl matrix
void matrix_print(gsl_matrix* mat);

// binary search for gsl_vector
int binsearch_vector(gsl_vector* x, double x_new); 

// binary search for plain c array
int binsearch_array(int N, double* x, double x_new); 

