#include <stdio.h>
void sub(int m,int n,double *vector,double *matrix) {
  int i,j;
  printf("in sub\n");
  printf("vector = ");
  for (j=0;j<n;j++) printf("%.4g ",vector[j]);
  printf("\n");
  for (i=0;i<m;i++) {
    printf("matrix = ");
    for (j=0;j<n;j++) printf("%.4g ",matrix[i+j*m]);
    printf("\n");
  }
}
