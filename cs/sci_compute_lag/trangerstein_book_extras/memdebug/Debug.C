//**********************************************************************
// Copyright 2006 John A. Trangenstein
//
// This software is made available for research and instructional use 
// only. 
// You may copy and use this software without charge for these 
// non-commercial purposes, provided that the copyright notice and 
// associated text is reproduced on all copies.  
// For all other uses (including distribution of modified versions), 
// please contact the author at
//   John A. Trangenstein
//   Department of Mathematics
//   Duke University
//   Durham, NC 27708-0320
//   USA
// or
//   johnt@math.duke.edu
// 
// This software is made available "as is" without any assurance that it
// is completely correct, or that it will work for your purposes.  
// Use the software at your own risk.
//**********************************************************************
// "$Header:$"
#include "Debug.H"
#include <complex>
#include <iostream>
#include <sys/times.h>
#include "Errors.H"
extern "C" {
  void F77_NAME(print_loc)(const void *p) {
    cout << p << endl;
    cout << "\tdouble = " << *static_cast<const double*>(p) << endl;
    cout << "\tfloat = " << *static_cast<const float*>(p) << endl;
    cout << "\tint = " << *static_cast<const int*>(p) << endl;
  }
  void F77_NAME(print_double)(const double *p) {
    cout << p << "\tdouble = " << *static_cast<const double*>(p) << endl;
  }
  void F77_NAME(print_float)(const float *p) {
    cout << p << "\tfloat = " << *static_cast<const float*>(p) << endl;
  }
  void F77_NAME(print_int)(const int *p) {
    cout << p << "\tint = " << *static_cast<const int*>(p) << endl;
  }
  void F77_NAME(check_double_size)(const double *p,const double *p1) {
    CHECK_SAME(p1,p+1)
  }
  void F77_NAME(check_single_size)(const float *p,const float *p1) {
    CHECK_SAME(p1,p+1)
  }
  void F77_NAME(check_complex_size)(const complex<double> *p,
  const complex<double> *p1) {
    CHECK_SAME(p1,p+1)
  }
  void F77_NAME(check_integer_size)(const int *p,const int *p1){
    CHECK_SAME(p1,p+1)
   }
  void F77_NAME(check_logical_size)(const bool *p,const bool *p1) {
    CHECK_SAME(p1,p+1)
  }
}
