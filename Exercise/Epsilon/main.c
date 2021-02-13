#include<stdio.h>
#include<limits.h>
#include<float.h>

int equal(double a, double b, double tau, double epsilon);

int main(){
	printf("EXERCISE 1.i\n\n");
	printf("The maximum value of INT = %d\n",INT_MAX);
	int i=1; while(i+1>i) {i++;}
	printf("my max int (WHILE-LOOP) = %i\n",i);
	int j=1;
	for(j=1; j<j+1; j++){} printf("my max int (FOR-LOOP) = %d\n",j);
	int k=1;
	do {k++;} while(k<k+1); printf("my max int (DO-WHILE-LOOP) = %d\n",k);
	
	printf("\n \nEXERCISE 1.ii\n\n");
        printf("The minimum value of INT = %d\n",INT_MIN);
        while(i-1<i) {i++;}
        printf("my min int (WHILE-LOOP) = %i\n",i);
        for(j=1; j>j-1; j++){} printf("my min int (FOR-LOOP) = %d\n",j);
        do {k++;} while(k>k-1); printf("my min int (DO-WHILE-LOOP) = %d\n",k);
	

	printf("\n\nEXERCISE 1.iii\n\n");
	printf("FLT_EPSILON value=%g\n",FLT_EPSILON);

	float x=1; while(1+x!=1){x/=2;} x*=2;
	printf("The machine FLT_epsilon for 1 (WHILE) =%g\n",x);
	float y=1;
	do{y/=2;} while(1+y!=1); y*=2;
	printf("The machine FLT_epsilon for 1 (DO-WHILE)=%g\n",y);
	float z=1;
	for(z=1;1+z!=1;z/=2){} z*=2;
	printf("The machine FLT_epsilon for 1 (FOR)=%g\n",z);

	printf("\n\nDBL_EPSILON value=%g\n",DBL_EPSILON);

        double x1=1; while(1+x1!=1){x1/=2;} x1*=2;
        printf("The machine DBL_epsilon for 1 (WHILE) =%g\n",x1);
        double y1=1;
        do{y1/=2;} while(1+y1!=1); y1*=2;
        printf("The machine DBL_epsilon for 1 (DO-WHILE)=%g\n",y1);
        double z1=1;
        for(z1=1;1+z1!=1;z1/=2){} z1*=2;
        printf("The machine DBL_epsilon for 1 (FOR)=%g\n",z1);
		

	printf("\n\nLDBL_EPSILON value=%Lg\n",LDBL_EPSILON);

        long double x2=1; while(1+x2!=1){x2/=2;} x2*=2;
        printf("The machine LDBL_epsilon for 1 (WHILE) =%Lg\n",x2);
        long double y2=1;
        do{y2/=2;} while(1+y2!=1); y2*=2;
        printf("The machine LDBL_epsilon for 1 (DO-WHILE)=%Lg\n",y2);
        long double z2=1;
        for(z2=1;1+z2!=1;z2/=2){} z2*=2;
        printf("The machine LDBL_epsilon for 1 (FOR)=%Lg\n",z2);

	printf("\n\nExercise 2\n");
	float f=INT_MAX/300;
	float sum_up_float=0;
	for(float i=1;i<f;i++){
	 sum_up_float+= 1.0*f/i;
	}
	printf("Value for sum_up_float =%g\n",sum_up_float);

	
	float sum_down_float=0;
	for(float i=1;i<f;i++){
         sum_down_float+= 1.0*f/(f-i);
        }
        printf("Value for sum_down_float =%g\n",sum_down_float);

	printf("\n\nExercise 2.ii\n");
	printf("The difference is the direction of sum");
	
	printf("\n\nExercise 2.iii\n");
	printf("Yes, just takes more compile-time");

	printf("\n\nExercise 2.iv\n");
	double d=INT_MAX/300;
        double sum_up_double=0;
        for(double i=1;i<d;i++){
         sum_up_double+= 1.0*d/i;
        }
        printf("Value for sum_up_double =%g\n",sum_up_double);


        double sum_down_double=0;
        for(double i=1;i<d;i++){
         sum_down_double+= 1.0*d/(d-i);
        }
        printf("Value for sum_down_double =%g\n",sum_down_double);	
	
	equal(2,2,2,2);
	printf("%i\n",equal(2,2,2,2));
	return 0;
}
