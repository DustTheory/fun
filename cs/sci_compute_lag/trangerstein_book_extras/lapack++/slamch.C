#include <assert.h>
#include <float.h>
#include <iostream>
#include "arch.H"

// slamch_ determines single precision machine parameters.   
extern "C" {
  static bool first=true;
  float F77NAME(slamch)(char *cmach) {
    static float rmax, sfmin, prec;

    if (first) {
      first=false;
      prec = FLT_EPSILON * 2.;
      rmax  = FLT_MAX * (1.f-FLT_EPSILON);
      sfmin = FLT_MIN;
      if (1.f >= rmax * FLT_MIN) sfmin = (FLT_EPSILON+1.f) / rmax;
    }

    switch (*cmach) {
      case 'E' : // eps = relative machine precision
      case 'e' :
	return FLT_EPSILON;
      case 'S' : // sfmin = safe minimum: 1/sfmin does not overflow
      case 's' :
	return sfmin;
      case 'B' : // base
      case 'b' :
	return 2.;
      case 'P' : // eps*base
      case 'p' :
	return prec;
      case 'N' : // number of (base) digits in the mantissa
      case 'n' :
	return FLT_MANT_DIG;
      case 'R' : // 1.0 when rounding occurs in addition, 0.0 otherwise
      case 'r' :
	return 1.;
      case 'M' : // minimum exponent before (gradual) underflow
      case 'm' :
	return FLT_MIN_EXP;
      case 'U' : // underflow threshold = base**(emin-1)
      case 'u' :
	return FLT_MIN;
      case 'L' : // largest exponent before overflow
      case 'l' :
	return FLT_MAX_EXP;
      case 'O' : // overflow threshold  = (base**emax)*(1-eps)
      case 'o' :
	return rmax;
      default:
	assert(0);
	return 0.;
    }
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   
*/
  }
}
