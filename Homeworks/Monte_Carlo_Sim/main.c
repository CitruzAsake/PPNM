#define SQR2 1.414
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<complex.h>
#include<stdio.h>

void Halton(int n, int dim, double *a, double *b, double *x);
void pseudoMC(int dim, double f(int dim, double* x), double* a, double* b, double N, double* result2, double* error2);
void plainMC(int dim,double f(int dim,double* x),double* a,double* b,int N, double* result, double* error);

double mainfunc(int dim, double* x)
{
   return 1.0/(1 - cos(x[0])*cos(x[1])*cos(x[2]))/(M_PI*M_PI*M_PI);
}

double Func1(double x){return sqrt(x);};
double Func2(double x){return 4.*sqrt(1-(x)*(x));};
double Func3(double x){return 1/(sqrt(x));}
double Func4(double x){return log(x)/sqrt(x);}

int main(){
{
	double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        double a[dim];
        double b[dim];
	
	for (int j = 0; j < dim; ++j){
            a[j] = 0;
            b[j] = 1;
        }
	
	//Number of points
	int N = 1e3;
	double precise = 2/3.;

	plainMC(dim, Func1, a, b, N, &result, &error);
        pseudoMC(dim, Func1, a, b, N, &result2, &error2);
	
	printf("----------------------------------------\n");
        printf("∫ dx √(x) from 0 to 1 using the plain Monte Carlo algorithm \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", precise);
        printf("Diff  = %.10f \n", precise-result);
	printf("----------------------------------------\n");
                printf("∫ dx √(x) from 0 to 1 using the pseudo Monte Carlo algorithm \n");
        printf("Value = %.10f \n", result2);
        printf("Error = %.10f \n", error2);
        printf("Exact = %.10f \n", precise);
        printf("Diff  = %.10f \n", precise-result2);
    }
	
	{
        double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        double a[dim];
        double b[dim];
        for (int i = 0; i < dim; ++i){
            a[i] = 0;
            b[i] = 1;
        }
	//NUmber of points
	int N = 1e3;

	plainMC(dim, Func2, a, b, N, &result, &error);
        pseudoMC(dim, Func2, a, b, N, &result2, &error2);
	double precise = M_PI;
	printf("----------------------------------------\n");
        printf("∫ dx √(x) from 0 to 1 using the plain Monte Carlo algorithm \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", precise);
        printf("Diff  = %.10f \n", precise-result);
        printf("----------------------------------------\n");
                printf("∫ dx √(x) from 0 to 1 using the pseudo Monte Carlo algorithm \n");
        printf("Value = %.10f \n", result2);
        printf("Error = %.10f \n", error2);
        printf("Exact = %.10f \n", precise);
        printf("Diff  = %.10f \n", precise-result2);
    	printf("----------------------------------------\n");

	FILE* errorscale = fopen("errorscale.txt", "w");
        for (int N = 1e5; N <= 2e6; N += 1e5) {
        plainMC(dim, Func2, a, b, N, &result , &error);
        pseudoMC(dim, Func2, a, b, N, &result2, &error2);
        fprintf(errorscale, "%10.2e %20.15e %20.15e\n", (double)N, fabs(result - precise), fabs(result2 - precise));

	}
}

 {
        double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        double a[dim];
        double b[dim];
        for (int i = 0; i < dim; ++i){
            a[i] = 0;
            b[i] = 1;
        }
	//Number of points 
	int N = 1e3;
	plainMC(dim, Func3, a, b, N, &result, &error);
        pseudoMC(dim, Func3, a, b, N, &result2, &error2);
	double precise = 2.;
	printf("----------------------------------------------\n");
	printf("∫ dx/√(x) from 0 to 1 using plain Monte Carlo \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", precise);
        printf("Diff  = %.10f \n", precise-result);
        }

	{
        double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        double a[dim];
        double b[dim];
        for (int i = 0; i < dim; ++i){
            a[i] = 0;
            b[i] = 1;
        }
	//NUmber of points 
	int N =1e3;
	plainMC(dim, Func4, a, b, N, &result, &error);
        pseudoMC(dim, Func4, a, b, N, &result2, &error2);
        double precise = -4.;

	printf("∫ ln(x)/√(x) from 0 to 1 using the plain Monte Carlo algorithm\n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", precise);
        printf("Diff  = %.10f \n", precise-result);
            }

	{
	double result;
        double error;
        double result2;
        double error2;
        int dim = 1;
        double a[dim];
        double b[dim];
        for (int i = 0; i < dim; ++i){
            a[i] = 0;
            b[i] = 1;
        }
        //NUmber of points 
        int N =1e3;
	plainMC(dim, mainfunc, a, b, N, &result, &error);
        pseudoMC(dim, mainfunc, a, b, N, &result2, &error2);
	double precise = 1.3932;
        printf("--------------------------------------\n");
        printf("1/π³ ∫ dxdydz ([1-cos(x)cos(y)cos(z)]-1) from 0 to π using plain Monte Carlo \n");
        printf("Value = %.10f \n", result);
        printf("Error = %.10f \n", error);
        printf("Exact = %.10f \n", precise);
        printf("Diff  = %.10f \n", precise-result);
        
    }
return 0;
}



