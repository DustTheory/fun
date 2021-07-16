#include <stdlib.h> // drand48
#include <stdio.h> // printf
#include <string.h> // memcpy
int main(int argc,char **argv) {
  int i,j,inc,info;
  int n=10;
  int ipiv[n],jpiv[n];
  double A[n][n],LU[n][n],b[n],x[n],s,alpha,beta;
  char trans='N';

  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) A[j][i]=drand48();
  }
  memcpy(LU[0],A[0],n*n*sizeof(double));
  dgetc2_(&n,LU[0],&n,ipiv,jpiv,&info);
  printf("info = %d\n",info);
  printf("ipiv = ");
  for (i=0;i<n-1;i++) printf("%d ",ipiv[i]);
  printf("\n");
  printf("jpiv = ");
  for (j=0;j<n-1;j++) printf("%d ",jpiv[j]);
  printf("\n");

  for (i=0;i<n;i++) b[i]=drand48();
  memcpy(x,b,n*sizeof(double));
  dgesc2_(&n,LU[0],&n,x,ipiv,jpiv,&s);
  printf("scale = %.8g\n",s);
  printf("x = ");
  for (j=0;j<n;j++) printf("%.8g ",x[j]);
  printf("\n");

  alpha=-1.;
  beta=1.;
  inc=1;
  dgemv_(&trans,&n,&n,&alpha,A[0],&n,x,&inc,&beta,b,&inc);
  printf("error = ");
  for (i=0;i<n;i++) printf("%.8g ",b[i]);
  printf("\n");
}
