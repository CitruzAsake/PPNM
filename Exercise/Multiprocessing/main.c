#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#define N 100

void* monte_c(void* arg1){
  FILE* datafile=fopen("datafile.txt","a");
  double j,k,l;
  int *count= (int*)arg1;
  unsigned int SEED;
  for (int i=0; i<=N; i++){
    j= (double)rand_r(&SEED)/RAND_MAX;
    k= (double)rand_r(&SEED)/RAND_MAX;
    l=pow((j),2)+pow((k),2);
    fprintf(datafile,"%10g %10g \n",j,k);
    if (l<=1){ *count=*count+1;}
  }
  fclose(datafile);
  return NULL;
  
}

int main() {
  double pi; int c1=0, c2=0;
  pthread_t thread1;
  FILE* data1=fopen("datafile.txt","w");

  pthread_create(&thread1, NULL, monte_c, (void*) &c1);

  monte_c((void*)&c2);

  void* returnval = NULL;

  pthread_join(thread1,returnval);

  int tot_count=c1+c2;
  
  pi=(double)(tot_count)/(N*2) *4;
  printf("Number of interation is equal to %i, estimate of pi is given as %g\n",N*2,pi);
  printf("For larger N we get a better estimation of pi");
  FILE* circle=fopen("circledatafile.txt","w");
  double zk=1e-3;
  for(double j=0;j<=M_PI/2;j=j+zk){
    fprintf(circle,"%10g %10g\n",cos(j),sin(j));
  }

  fclose(data1);
  fclose(circle);
  return 0;
}
