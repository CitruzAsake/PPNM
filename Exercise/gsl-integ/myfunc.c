#include<stdio.h>
#include<math.h>
#include <gsl/gsl_integration.h>

double f (double x,void* params){
	double f = log(x)/sqrt(x);
	return f;
}

double myfunc(double z){
	gsl_function F;
	F.function=&f;
	F.params=(void*)&z;
	int limit=999;
	gsl_integration_workspace* w;
	w = gsl_integration_workspace_alloc (limit);
	double a=0,b=1,acc=1e-6,eps=1e-6,result,error;
	gsl_integration_qags(&F,a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	for(double x=1.0/8;x<=5;x+=1.0/8){}
	printf("%10gprintf("%10g\n",x,myfunc(x));
printf("%10g\n",x,myfunc(x));
 \n",myfunc(5));
return 0;
}

