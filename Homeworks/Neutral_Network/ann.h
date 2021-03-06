#include"gsl/gsl_vector.h"
#ifndef HAVE_ANN_H
#define HAVE_ANN_H
typedef struct { int n; double(*f)(double); gsl_vector* params; } ann;
ann*   ann_alloc   (int n,double(*f)(double));
void   ann_free    (ann* network);
double ann_response(ann* network,double x);
double ann_derivative(ann* network,double x);
double ann_integral(ann* network,double x);
void   ann_train(ann* network,gsl_vector* xs,gsl_vector* ys);
int quasiNewton(double F(gsl_vector* x), gsl_vector* x, double acc);
#endif
