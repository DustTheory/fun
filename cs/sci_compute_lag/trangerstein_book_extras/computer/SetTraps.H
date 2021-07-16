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
//taken from
// http://gcc.gnu.org/onlinedocs/gcc-3.3.6/g77/Floating_002dpoint-Exception-Handling.html
#ifndef setTraps_H
#define setTraps_H "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/SetTraps.H,v 1.1 2009/08/20 17:33:33 johnt Exp $"
#ifdef DEBUG
#define __USE_GNU 1
#include <fenv.h>
#include <iostream>
static void __attribute__((constructor)) trapfpe () {
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
}
#endif
#endif
