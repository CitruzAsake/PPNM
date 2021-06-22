
#include<math.h>
#include<gsl/gsl_matrix.h>

//Multiplies arbitary matrix A with Jacobi matrix
void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
		}
}
//Dots jacobi matrix and arbitary A chosen
void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
        for(int j=0;j<A->size2;j++){
                double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
                double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
                gsl_matrix_set(A,p,j,new_apj);
                gsl_matrix_set(A,q,j,new_aqj);
                }

}
//Diagonolizes using the rotation of Jacobi
int Jacobi(gsl_matrix* A, gsl_matrix* V){
	gsl_matrix_set_identity(V);
	int n=A->size1,changed,sweeps=0;
	do{
		changed=0;
		sweeps++;
		for(int p=0;p<n-1;p++)
		for(int q=p+1;q<n;q++)
		{
			double apq=gsl_matrix_get(A,p,q);
			double app=gsl_matrix_get(A,p,p);
			double aqq=gsl_matrix_get(A,q,q);
			double theta=0.5*atan2(2*apq,aqq-app);
			double c=cos(theta),s=sin(theta);
			double new_app=c*c*app-2*s*c*apq+s*s*aqq;
			double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
			if(new_app!=app || new_aqq!=aqq)
				{
				changed=1;
				timesJ(A,p,q,theta);
				Jtimes(A,p,q,-theta);
				timesJ(V,p,q,theta);
				}
		}
	}while(changed!=0);
	return sweeps;
}

