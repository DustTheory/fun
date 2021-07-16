#include <stdlib.h> // drand48
#include <stdio.h> // printf
#include <string.h> // memcpy
int main(int argc,char **argv) {
  int i,j,inc,info;
  int n=10;
  int ipiv[n];
  double A[n][n],LU[n][n],b[n],x[n],alpha,beta;
  char trans='N';

  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) A[j][i]=drand48();
  }
  memcpy(LU[0],A[0],n*n*sizeof(double));
  dgetrf_(&n,&n,LU[0],&n,ipiv,&info);
  printf("info = %d\n",info);
  printf("ipiv = ");
  for (i=0;i<n-1;i++) printf("%d ",ipiv[i]);
  printf("\n");

  for (i=0;i<n;i++) b[i]=drand48();
  memcpy(x,b,n*sizeof(double));
  inc=1;
  dgetrs_(&trans,&n,&inc,LU[0],&n,ipiv,x,&n,&info);
  printf("x = ");
  for (j=0;j<n;j++) printf("%.8g ",x[j]);
  printf("\n");

  alpha=-1.;
  beta=1.;
  dgemv_(&trans,&n,&n,&alpha,A[0],&n,x,&inc,&beta,b,&inc);
  printf("error = ");
  for (i=0;i<n;i++) printf("%.8g ",b[i]);
  printf("\n");
}
