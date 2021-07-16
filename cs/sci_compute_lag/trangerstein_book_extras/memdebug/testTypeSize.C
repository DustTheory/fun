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
#include <stdlib.h>
#include "Arch.H"
#include "MemoryDebugger.H"
#include "Tracer.H"

extern "C" {
  void F77_NAME(test_types)();
}

int main( int /*argc*/, char** /*argv*/) {
  cout << boolalpha;
  { TRACER_CALL(t,"test_types");
    F77_NAME(test_types)();
  }
  { MemoryDebugger md(1);
//  constructor calls test_types
  }
  return 1;
}
