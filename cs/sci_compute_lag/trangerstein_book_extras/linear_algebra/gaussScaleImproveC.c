#include <stdlib.h> // drand48
#include <stdio.h> // printf
#include <string.h> // memcpy
int main(int argc,char **argv) {
  int i,j,info;
  double rcond,ferr[1],berr[1];
  int n=1024;
  int *ipiv,*iwork;
  double **A,**LU,*b,*x,*r,*c,*work,*space;
  char fact='E';
  char trans='N';
  char equed;
  int nrhs=1;

/*arrays are too big for static memory allocation*/
  A=(double**) malloc(n*sizeof(double*));
  space=(double*) malloc(n*n*sizeof(double));
  for (j=0;j<n;j++) {
    A[j]=space;
    space+=n;
  }
  LU=(double**) malloc(n*sizeof(double*));
  space=(double*) malloc(n*n*sizeof(double));
  for (j=0;j<n;j++) {
    LU[j]=space;
    space+=n;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) A[j][i]=drand48();
  }

  b=(double*) malloc(n*sizeof(double));
  x=(double*) malloc(n*sizeof(double));
  r=(double*) malloc(n*sizeof(double));
  c=(double*) malloc(n*sizeof(double));
  work=(double*) malloc(4*n*sizeof(double));
  ipiv=(int*) malloc(n*sizeof(int));
  iwork=(int*) malloc(n*sizeof(int));

  for (i=0;i<n;i++) b[i]=drand48();
  dgesvx_(&fact,&trans,&n,&nrhs,A[0],&n,LU[0],&n,ipiv,&equed,r,c,b,&n,x,&n,
    &rcond,ferr,berr,work,iwork,&info);
  printf("info,equed,rcond,ferr,berr = %d %c %.8g %.8g %.8g\n",
    info,equed,rcond,ferr[0],berr[0]);

  free(iwork); iwork=0;
  free(ipiv); ipiv=0;
  free(work); work=0;
  free(c); c=0;
  free(r); r=0;
  free(x); x=0;
  free(b); b=0;
  free(LU[0]); LU[0]=0;
  free(LU); LU=0;
  free(A[0]); A[0]=0;
  free(A); A=0;
}
