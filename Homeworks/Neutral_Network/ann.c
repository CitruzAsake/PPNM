#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdio.h>
#include"ann.h"

int n, initialpoint;

double activation(double x){
	return x*exp(-x*x);
}

double derivative_func(double x){
	return exp(-x*x)*(1-2*pow(x,2));
}

double integral_func(double x){
	return (-1)*exp(-x*x)/2;
}

//allocation of memory
ann* ann_alloc(int n,double(*f)(double)){
	ann* network = malloc(sizeof(ann));
	network->n=n;
	network->f=f;
	network->params=gsl_vector_alloc(3*n);
	return network;
} 


//Free memory
void ann_free(ann* network){
	gsl_vector_free(network->params);
	free(network);
}


double ann_response(ann* network,double x){
	double w=0;
	for(int i=0;i<network->n;i++){
		double a=gsl_vector_get(network->params,3*i);
		double b=gsl_vector_get(network->params,3*i+1);
		double s=gsl_vector_get(network->params,3*i+2);
		w+=network->f((x-a)/b)*s;
	}
	return w;
}

//pulling derivative from ann
double ann_derivative(ann* network,double x)
{
	double total = 0;
	for(int i=0; i<network->n; i++){
		double a = gsl_vector_get(network->params,3*i);
		double b = gsl_vector_get(network->params,3*i+1);
		double s = gsl_vector_get(network->params,3*i+2);
		total += derivative_func((x-a)/b)*s/b;
		}
	return total;
}

//pulling integral from ann
double ann_feedInt(ann* network,double x)
{
	double total = 0;
	for(int i=0; i<network->n; i++){
		double a = gsl_vector_get(network->params,3*i+0);
		double b = gsl_vector_get(network->params,3*i+1);
		double s = gsl_vector_get(network->params,3*i+2);
		total += integral_func((x-a)/b)*b*s;
		total -= integral_func((initialpoint-a)/b)*b*s;
		}
	return total;
}

int quasiNewton(double F(gsl_vector* x), gsl_vector*x, double acc);

void ann_train(ann* network, gsl_vector* xs, gsl_vector* ys){

double cost(gsl_vector* p){
		gsl_vector_memcpy(network->params,p);
		double total=0;
		for(int i=0;i<xs->size;i++){
			double xi=gsl_vector_get(xs,i);
			double yi=gsl_vector_get(ys,i);
			double fi=ann_response(network,xi);
			total+=fabs(fi-yi);
		}
		return total/xs->size;
	}

	gsl_vector* p=gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(p,network->params);

	quasinewton(cost,p,1e-5);

	gsl_vector_memcpy(network->params,p);
	gsl_vector_free(p);
}
