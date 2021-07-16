#include <math.h> // cos, sin
#include <stdio.h> /* for FILE, scanf, printf, fprintf */
#include <stdlib.h> /* for EXIT_SUCCESS */

double f(double x) { return sin(x); }
double fprime(double x) { return cos(x); }

int main(int argc,char** argv) {
  double x,h,exact,diff;
  int i,nlevels;
  FILE *out;

  x=2.;
  printf("enter x\n");
  scanf("%lf",&x); /*must use lf to scanf a double; lf same as f for printf*/
  h=1.;
  nlevels=64;

  exact=fprime(x);
  out=fopen("c_derivative_output","w");
  for (i=0;i<nlevels;i++,h*=0.5) {
    diff=( f(x+h) - f(x) )/h;
    fprintf(out,"%lf %lf\n",-log10(h),log10(fabs(diff-exact)));
  }
  fclose(out);
  return EXIT_SUCCESS;
}
