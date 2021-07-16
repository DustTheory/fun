#include <stdlib.h> // drand48
#include <stdio.h> // printf
#include <string.h> // memcpy
int main(int argc,char **argv) {
  int i,j,inc;
  int n=10;
  float R[n][n],y[n],x[n],Rx[n];
  char uplo='U';
  char trans='N';
  char diag='N';
  float a=-1.;

  for (j=0;j<n;j++) {
    for (i=0;i<=j;i++) R[j][i]=drand48();
  }

  for (i=0;i<n;i++) y[i]=drand48();
  printf("y = ");
  for (i=0;i<n;i++) printf("%.8g ",y[i]);
  printf("\n");

  memcpy(x,y,n*sizeof(float));
  inc=1;
  strsv_(&uplo,&trans,&diag,&n,R[0],&n,x,&inc);
  printf("x = ");
  for (i=0;i<n;i++) printf("%.8g ",x[i]);
  printf("\n");

  memcpy(Rx,x,n*sizeof(float));
  strmv_(&uplo,&trans,&diag,&n,R[0],&n,Rx,&inc);
  saxpy_(&n,&a,Rx,&inc,&y,&inc);
  printf("error = ");
  for (i=0;i<n;i++) printf("%.8g ",y[i]);
  printf("\n");
}
