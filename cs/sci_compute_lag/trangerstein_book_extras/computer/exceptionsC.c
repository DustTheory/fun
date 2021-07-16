#include <math.h>
#include <stdio.h> /* for printf */
#include <stdlib.h> /* for EXIT_SUCCESS */
#include <values.h> // for FLT_MAX, FLT_MIN, FLT_EPSILON

#define __USE_GNU 1
#include <fenv.h>

/*
 *  from bits/fenv.h:
 *  FE_INVALID = 0x1
 *  __FE_DENORM = 0x2
 *  FE_DIVBYZERO = 0x4
 *  FE_OVERFLOW = 0x8
 *  FE_UNDERFLOW = 0x10
 *  FE_INEXACT = 0x20
*/

int main(int argc,char **argv) {
  double d,z;
#if (_TRAP_FPE_==1)
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
  printf("exception codes:\n\t invalid = %d\n\t denorm = %d\n\t divide by zero = %d\n",FE_INVALID,__FE_DENORM,FE_DIVBYZERO);
  printf("\t overflow = %d\n\t underflow = %d\n\t inexact = %d\n\n",
    FE_OVERFLOW,FE_UNDERFLOW,FE_INEXACT);

  d=HUGE_VAL;
  printf("HUGE_VAL = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#if (_TRAP_FPE_!=1)
  d=1.;
  d /= HUGE_VAL;
  printf("1./HUGE_VAL = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  d=DBL_MAX;
  d *= 2.;
  printf("DBL_MAX * 2. = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  d=HUGE_VAL;
  d *= HUGE_VAL;
  printf("infinity * infinity = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#if (_TRAP_FPE_!=1)
  d=1.;
  z=0.;
  d /= z;
  printf("1. / 0. = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  d=0.;
  z=0.;
  d /= z;
  printf("0. / 0. = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  d=HUGE_VAL;
  d += HUGE_VAL;
  printf("infinity + infinity = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

#if (_TRAP_FPE_!=1)
  d=HUGE_VAL;
  d -= HUGE_VAL;
  printf("infinity - infinity = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  d=HUGE_VAL;
  d /= HUGE_VAL;
  printf("infinity / infinity = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  d=HUGE_VAL;
  d *= 0.;
  printf("infinity * 0. = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  d=log(-1.);
  printf("log(-1.) = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  d=log(0.);
  printf("log(0.) = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  d=sqrt(-1.);
  printf("sqrt(-1.) = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

#if (_TRAP_FPE_!=1)
  d=acos(2.);
  printf("acos(2.) = %16.4g 0x%lX %d\n",d,*(unsigned long*) &d,
    fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
  feclearexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  return EXIT_SUCCESS;
}
