#include <stdio.h>
#include <stdlib.h>
void sub(int,int,double*,double*);
int main(int argc,char** argv) {
  int m=3,n=2;
  int i,j;
  double *matrix=(double*) malloc(m*n*sizeof(double));
  double *vector=(double*) malloc(n*sizeof(double));

  printf("in main\n");
  for (j=0;j<n;j++) {
    for (i=0;i<m;i++) {
      matrix[i+j*m]=1./(double)(i+j+1);
    }
    vector[j]=1.;
  }
  sub(m,n,vector,matrix);
  printf("back in main\n");

  free (matrix); matrix=0;
  free (vector); vector=0;
  return EXIT_SUCCESS;
}
