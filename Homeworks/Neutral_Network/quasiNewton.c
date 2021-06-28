#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<math.h>
#include<float.h>


#ifndef DBL_EPSILON
#define DBL_EPSILON 2.22045e-16
#endif
static const double DELTA=sqrt(DBL_EPSILON);


void numeric_gradient(double F(gsl_vector*), gsl_vector*x, gsl_vector*grad)
{
	int m = x->size;
	double fx=F(x);

	for(int i=0;i<m;i++){
		double dx;
        double xi=gsl_vector_get(x,i);
		if(fabs(xi)<sqrt(DELTA)) dx=DELTA;
		else dx=fabs(xi)*DELTA;
		gsl_vector_set(x,i,xi+dx);
		gsl_vector_set(grad,i,(F(x)-fx)/dx);
		gsl_vector_set(x,i,xi);
	}
}

int quasiNewton(double F(gsl_vector* x), gsl_vector*x, double acc) 
{
	int dim=x->size;
	int steps=0;
	int d=0;
	int g=0;
	double fz;

	gsl_matrix* B=gsl_matrix_alloc(dim,dim);
	gsl_vector* gradient=gsl_vector_alloc(dim);
	gsl_vector* Dx=gsl_vector_alloc(dim);
	gsl_vector* z=gsl_vector_alloc(dim); 
	gsl_vector* gz=gsl_vector_alloc(dim);
	gsl_vector* y=gsl_vector_alloc(dim);
	gsl_vector* u=gsl_vector_alloc(dim);
	gsl_vector* a=gsl_vector_alloc(dim);
	gsl_matrix_set_identity(B);
	numeric_gradient(F,x,gradient);
	double fx=F(x);
		
	while(steps<1000){
		steps++;
		gsl_blas_dgemv(CblasNoTrans,-1,B,gradient,0,Dx);
		if(gsl_blas_dnrm2(Dx)<DELTA*gsl_blas_dnrm2(x))
			{fprintf(stderr,"quasinewton: |Dx|<DELTA*|x|\n"); break;} 
		if(gsl_blas_dnrm2(gradient)<acc)
			{fprintf(stderr,"quasinewton: |grad|<acc\n"); break;}
		double gamma=1;
		while(1){
			gsl_vector_memcpy(z,x);
			gsl_vector_add(z,Dx);
			fz=F(z);
			double sTg; gsl_blas_ddot(Dx,gradient,&sTg);
			if(fz<fx+0.01*sTg){ g++; break; } 
			if(gamma<DELTA){
				d++;
				gsl_matrix_set_identity(B);
				break;
				}
			gamma*=0.5;
			gsl_vector_scale(Dx,0.5);
		}
		numeric_gradient(F,z,gz);
		gsl_vector_memcpy(y,gz);
		gsl_blas_daxpy(-1,gradient,y);
		gsl_vector_memcpy(u,Dx);
		gsl_blas_dgemv(CblasNoTrans,-1,B,y,1,u); 
		double sTy,uTy;
		gsl_blas_ddot(Dx,y,&sTy);
		if(fabs(sTy)>1e-12){ 
			gsl_blas_ddot(u,y,&uTy);
			double lambda=uTy/2/sTy;
			gsl_blas_daxpy(-lambda,Dx,u); 
			gsl_blas_dger(1.0/sTy,u,Dx,B);
			gsl_blas_dger(1.0/sTy,Dx,u,B);
		}
		gsl_vector_memcpy(x,z);
		gsl_vector_memcpy(gradient,gz);
		fx=fz;
	}
gsl_matrix_free(B);gsl_vector_free(gradient);
gsl_vector_free(Dx);gsl_vector_free(z);
gsl_vector_free(gz);gsl_vector_free(y);gsl_vector_free(u);gsl_vector_free(a);
fprintf(stderr,"quasiNewton: steps=%i g=%i d=%i fx=%.1e\n",steps,g,d,fx);
return steps;
}
