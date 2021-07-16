#include <math.h> // for HUGE_VALF
#include <stdio.h> /* for printf */
#include <stdlib.h> /* for EXIT_SUCCESS */
#include <values.h> /* for FLT_MAX, FLT_MIN, FLT_EPSILON */

typedef union {
  long double ld;
  char bytes[sizeof(long double)];
} converter;
void print_converter(converter c) {
  int i;
  printf("0x");
  for (i=0;i<sizeof(long double);i++) printf("%x",c.bytes[i]);
}

int main(int argc,char **argv) {
  float f;
  double d;
  converter c;
  float aold,anew,a;
  printf("sizeof( float ) = %d\n",sizeof(float));
  printf("sizeof( double ) = %d\n",sizeof(double));
  printf("sizeof( long double ) = %d\n",sizeof(long double));
#if (__BYTE_ORDER == __LITTLE_ENDIAN)
  printf("__BYTE_ORDER == __LITTLE_ENDIAN\n");
#else
  printf("__BYTE_ORDER == __BIG_ENDIAN\n");
#endif

  f=HUGE_VALF;
  printf("HUGE_VALF = %32.24g %#x\n",f,*(unsigned int*) &f);
  f=FLT_MAX;
  printf("FLT_MAX = %32.24g %#x\n",f,*(unsigned int*) &f);
  f=FLT_MIN;
  printf("FLT_MIN = %32.24g %#x\n",f,*(unsigned int*) &f);
  f=FLT_EPSILON;
  printf("FLT_EPSILON = %32.24g %#x\n\n",f,*(unsigned int*) &f);

  d=HUGE_VAL;
  printf("HUGE_VAL = %56.48g 0x%lX\n",d,*(unsigned long*) &d);
  d=DBL_MAX;
  printf("DBL_MAX = %56.48g 0x%lX\n",d,*(unsigned long*) &d);
  d=DBL_MIN;
  printf("DBL_MIN = %56.48g 0x%lX\n",d,*(unsigned long*) &d);
  d=DBL_EPSILON;
  printf("DBL_EPSILON = %56.48g 0x%lX\n\n",d,*(unsigned long*) &d);

  c.ld=HUGE_VALL;
  printf("HUGE_VALL = %120.112Lg ",c.ld);
  print_converter(c);
  printf("\n");
  c.ld=LDBL_MAX;
  printf("LDBL_MAX = %120.112Lg ",c.ld);
  print_converter(c);
  printf("\n");
  c.ld=LDBL_MIN;
  printf("LDBL_MIN = %120.112Lg ",c.ld);
  print_converter(c);
  printf("\n");
  c.ld=LDBL_EPSILON;
  printf("LDBL_EPSILON = %120.112Lg ",c.ld);
  print_converter(c);
  printf("\n");

//rounding error occurs after mantissa bits become all 1's
  printf("compute a_n = 0.5*a_{n-1} + 1., a_1 = 1 until rounding error\n");
  aold=1.;
  anew=1.+0.5*aold;
  while (anew < 2.) {
    aold=anew;
    printf("%32.24g %#x\n",anew,*(unsigned int*) &anew);
    anew=1.+0.5*aold;
  }
  printf("%32.24g %#x\n\n",anew,*(unsigned int*) &anew);

//overflow
//bits in mantissa ordered most significant first
//rounding error occurs when number becomes even
//  the rounding error turns number into a power of 2
  printf("compute a_n = 2.*a_{n-1} + 1., a_1 = 1 until overflow\n");
  a=1.;
  while (isfinite(a)) {
    printf("%32.24g %#x\n",a,*(unsigned int*) &a);
    a=2.*a+1.;
  }
  printf("%32.24g %#x\n\n",a,*(unsigned int*) &a);

//underflow
//de-normal numbers occur when biased exponent becomes zero, and
//  then the mantissa bits are no longer all ones
  printf("compute a_n = 0.5*a_{n-1}, a_1 = 2-epsilon while positive\n");
  while (aold > 0.) {
    printf("%32.24g %#x\n",aold,*(unsigned int*) &aold);
    aold=0.5*aold;
  }
  printf("%32.24g %#x\n",aold,*(unsigned int*) &aold);
  return EXIT_SUCCESS;
}
