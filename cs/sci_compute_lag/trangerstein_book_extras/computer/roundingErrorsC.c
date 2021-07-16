#include <math.h> /* for log */
#include <stdio.h> /* for printf */
#include <stdlib.h> /* for EXIT_SUCCESS */
#include <values.h> /* for DBL_EPSILON */

int main(int argc,char **argv) {
  double eps,lambda,rho;
  int exponent,i;

  eps=DBL_EPSILON;
  printf("find smallest number eps such that 1 + eps > 1:\n");
  for (i=0;i<2;i++) {
    lambda=1.+eps;
    rho=(lambda-1.)/eps;
    exponent=(int) (log(eps)/log(2.)+0.5);
    printf("eps = 2^{%d} = %lX , ((1.+eps)-1.)/eps = %1.0f\n",exponent,
      *(unsigned long*) &eps,rho);
    eps*=0.5;
  }
  
  printf("\nfind largest number eps such that 1 + eps = 1:\n");
  eps=DBL_EPSILON;
  for (i=0;i<2;i++) {
    lambda=1.+eps;
    rho=(lambda-1.)/eps;
    exponent=(int) (log(eps)/log(2.)+0.5);
    printf("eps = 2^{%d} = %lX , ((1.+eps)-1.)/eps = %1.0f\n",exponent,
      *(unsigned long*) &eps,rho);
    eps*=0.5;
  }
  
  printf("\nfind largest number eps such that ( 1 / eps - 1 ) * eps = 1:\n");
  eps=DBL_EPSILON;
  for (i=0;i<3;i++) {
    lambda=1./eps;
    rho=(lambda-1.)*eps;
    exponent=(int) (log(eps)/log(2.)+0.5);
    printf("eps = 2^{%d} = %lX , (1./eps-1.)*eps = %18.16g\n",exponent,
      *(unsigned long*) &eps,rho);
    eps*=0.5;
  }

  printf("\nfind largest number eps such that ( 1 - 1 / eps ) + 1 / eps = 0:\n");
  eps=DBL_EPSILON;
  for (i=0;i<3;i++) {
    lambda=1./eps;
    rho=(1.-lambda)+lambda;
    exponent=(int) (log(eps)/log(2.)+0.5);
    printf("eps = 2^{%d} = %lX , (1.-1./eps)+1./eps = %1.0f\n",exponent,
      *(unsigned long*) &eps,rho);
    eps*=0.5;
  }
  
  return EXIT_SUCCESS;
}
