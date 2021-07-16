#include <math.h> /* for abs, ceil, max, pow, sqrt, HUGE_VAL */
#include <stdio.h> /* for scanf, printf */
#include <stdlib.h> /* for EXIT_SUCCESS */

double max(double x,double y) { return ( x > y ? x : y ); }
int main(int argc,char **argv) {
  double a,b,c,m,disc,large,small;
  a=HUGE_VAL;
  b=HUGE_VAL;
  c=HUGE_VAL;
  printf("enter a,b,c for quadratic a * x^2 + b * x + c = 0\n");
  scanf("%lf %lf %lf",&a,&b,&c);
  printf("a,b,c = %24.16lg %24.16lg %24.16lg\n",a,b,c);

  m=max(abs(a),max(abs(b),abs(c)));
  if (m<=0.) {
    printf("all coefficients zero, so roots are arbitrary\n");
  } else {
    m = pow(2.,ceil(log(m)/log(2.)));
    a /= m;
    b /= m;
    c /= m;
    if (abs(a)>0.) {
      b /= -2. * a;
      c /= a;
      disc = b * b - c;
      if (disc>=0.) {
        disc=sqrt(disc);
        large = b + ( b > 0. ? disc : -disc );
        small = ( abs( large ) > 0. ? c / large : 0.);
        printf("two real roots = %24.16lg %24.16lg\n",large,small);
      } else {
        printf("complex roots: real part = %24.16lg imaginary part = %24.16lg\n",b,sqrt(-disc));
      }
    } else {
      if (abs(b)>0.) {
        printf("one real root = %24.16lg\n",-c/b);
      } else {
        printf("no roots: a = 0 = b");
      }
    }
  }
  return EXIT_SUCCESS;
}
