#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>
#include <assert.h>

int binsearch(int n, double* x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
	}

double linterp(int n, double x[], double y[], double z){
		int i=binsearch(n,x,z);
                double pi=(y[i+1]-y[i])/(x[i+1]-x[i]);
                
		double s=y[i]+pi*(z-x[i]);
		return s; 
              	}


int main (void) {
	int n=9;
	double a=-2,b=2,x[n],y[n];


	for(int i=0;i<n;i++){
		x[i]=a+(b-a)*i/(n-1);
		y[i]=1/(1+pow(x[i],2));
	}

	printf("# index 0: data\n");
	for(int i=0;i<n;i++) printf("%g %g\n",x[i],y[i]);
	printf("\n\n");
	
	gsl_interp* linear     = gsl_interp_alloc(gsl_interp_linear    ,n);
	gsl_interp_init(linear    ,x,y,n);
	
	double z=1.0/5;
	printf("# index 1: interpolations\n");

	for(double i=x[0];i<=x[n-1];i+=z){
		printf("%g %g\n",i,linterp(n,x,y,i));
      	}	
	printf("\n\n");

	gsl_interp_free(linear);
	
	return 0;
}
