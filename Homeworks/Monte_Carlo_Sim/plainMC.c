#include<stdlib.h>
#include<complex.h>
#define RND (double)rand()/RAND_MAX
#include<math.h>

void plainMC(int dim,double f(int dim,double* x),double* a,double* b,int N, double* result, double* error)
{
	double Q=1;
	for(int i=0; i<dim; i++){
	Q*=b[i]-a[i]; 
		}
	double total=0, total2=0, x[dim];
	for(int i=0;i<N;i++){
		for(int j=0; j<dim;j++){
			x[i]=a[i]+RND*(b[i]-a[i]);
			}
			double yx=f(dim,x);
			total+=yx;
			total2+=yx*yx;
		}
	double mean=total/N;
	double sigma= sqrt(total2/N-mean*mean);
		*result=mean*Q;
		*error= sigma*Q/sqrt(N);
 }
