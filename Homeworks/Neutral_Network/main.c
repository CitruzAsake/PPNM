#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<stdio.h>
#include<math.h>
#include"ann.h"


double activation(double x);

double activation2(double x){
	return exp(-x*x);
}

double activation3(double x){
	return cos(5*x)*exp(-x*x);
}

double fit(double x){
	return cos(5*x-1)*exp(-x*x);
}

int main()
{
    int n=5;
	ann* network=ann_alloc(n,activation);
	double min=-1,max=1;
	int nx=20;
	gsl_vector* vx=gsl_vector_alloc(nx);
	gsl_vector* vy=gsl_vector_alloc(nx);

	for(int i=0;i<nx;i++){
		double x=min+(max-min)*i/(nx-1);
		double f=fit(x);
		gsl_vector_set(vx,i,x);
		gsl_vector_set(vy,i,f);
	}
	
	for(int i=0;i<network->n;i++){
		gsl_vector_set(network->params,3*i,min+(max-min)*i/(network->n-1));
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,1);
	}
	ann_train(network,vx,vy);
	for(int i=0;i<vx->size;i++){
		double x=gsl_vector_get(vx,i);
		double f=gsl_vector_get(vy,i);
		printf("%g %g\n",x,f);
	}
	printf("\n\n");
	double diffz=1.0/64;
	
	for(double z=min;z<=max;z+=diffz){
	double y  = ann_response(network,z);
    	double yprime  = derivative_func(network,z);
    	double yprimeprime  = integral_func(network,z);
		printf(" %.12e %.12e %.12e %.12e\n",z,y,yprime, yprimeprime);
 	   }
	for(int i=0;i<network->n;i++){
		double ai=gsl_vector_get(network->params,3*i);
		double bi=gsl_vector_get(network->params,3*i+1);
		double si=gsl_vector_get(network->params,3*i+2);
		fprintf(stderr,"i=%i ai,bi,si = %g %g %g\n",i,ai,bi,si);
	}

gsl_vector_free(vx);
gsl_vector_free(vy);
ann_free(network);

return 0;
}

