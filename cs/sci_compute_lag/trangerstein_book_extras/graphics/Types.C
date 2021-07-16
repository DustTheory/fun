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
// "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/Types.C,v 1.1 2009/08/20 17:33:33 johnt Exp $"
#include "Types.H"
#include <iostream>
#include <stdlib.h>
#include "Arch.H"

extern "C" {
//type-specific and not otherwise available in Fortran
  double F77NAME(drand48)() { return drand48(); }
}

ostream& operator<<(ostream &os,HAND hand) {
  switch (hand) {
    case LEFT_HAND:
      os << "LEFT_HAND";
      break;
    case RIGHT_HAND:
      os << "RIGHT_HAND";
      break;
    default:
      os << "unknown";
      break;
  } 
  return os;
}

/*
ostream& operator<<(ostream &os,READ_PLOT_RECORD rpr) {
  switch (rpr) {
    case read_first_marks:
      os << "read_first_marks";
      break;
    case read_until_bracket:
      os << "read_until_bracket";
      break;
    default:
      os << "unknown";
      break;
  } 
  return os;
}
*/
