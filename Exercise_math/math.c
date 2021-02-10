#include<math.h>
#include<stdio.h>
#include<complex.h>

int main(){

	double complex z;
	z=csqrt(-2);
	
	printf("gamma(5) =%g \n",tgamma(5.0));
	printf("bessel(0.5) =%g \n",j1(0.5));
	printf("sqrt(-2) =%g +I%g\n",creal(z),cimag(z));
	printf("cexp(I*PI)=%g +I%g\n",creal(cexp(I*M_PI)),cimag(cexp(I*M_PI)));
	printf("cexp(I)=%g + I%g\n",creal(cexp(I)),cimag(cexp(I)));
	printf("cpow(I,M_E)=%g + I%g\n",creal(cpow(I,M_E)),cimag(cpow(I,M_E)));
	printf("cpow(I,I)=%g+ I%g\n",creal(cpow(I,I)),cimag(cpow(I,I)));



	return 0;

}


