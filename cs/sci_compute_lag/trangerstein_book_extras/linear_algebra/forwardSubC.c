#include <stdlib.h> // drand48
#include <stdio.h> // printf
#include <string.h> // memcpy
int main(int argc,char **argv) {
  int i,j,inc;
  int m=10;
  float L[m][m],b[m],y[m],Ly[m];
  char uplo='L';
  char trans='N';
  char diag='U';
  float a=-1.;

  for (j=0;j<m;j++) {
    L[j][j]=1.;
    for (i=j+1;i<m;i++) L[j][i]=drand48();
  }

  for (i=0;i<m;i++) b[i]=drand48();
  printf("b = ");
  for (i=0;i<m;i++) printf("%.8g ",b[i]);
  printf("\n");

  memcpy(y,b,m*sizeof(float));
  inc=1;
  strsv_(&uplo,&trans,&diag,&m,L[0],&m,y,&inc);
  printf("y = ");
  for (i=0;i<m;i++) printf("%.8g ",y[i]);
  printf("\n");

  memcpy(Ly,y,m*sizeof(float));
  strmv_(&uplo,&trans,&diag,&m,L[0],&m,Ly,&inc);
  saxpy_(&m,&a,Ly,&inc,&b,&inc);
  printf("error = ");
  for (i=0;i<m;i++) printf("%.8g ",b[i]);
  printf("\n");
}
