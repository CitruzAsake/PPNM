#include <stdio.h>
#include <math.h>

int main(){
	FILE *IF=fopen("input.txt","r");
        FILE *OF=fopen("outputfiles.txt","w");

	double x; 
	int items;
	do{
		items=fscanf(IF,"%lg",&x);
		fprintf(OF,"x=%g cos(x)=%g\n",x,cos(x));
	}while(items!=EOF);

fclose(IF);
fclose(OF);
	return 0;
}
