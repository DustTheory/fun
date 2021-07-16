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
// "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/testTimedObject.C,v 1.1 2009/08/20 17:33:33 johnt Exp $"
#include <iostream>
#include <stdlib.h>
#include "MemoryDebugger.H"
#include "TimedObject.H"

int main(int /*argc*/,char** /*argv[]*/) {
  cout << boolalpha;
  {
    MemoryDebugger md(1);
    {
      TimedObject *z=OPERATOR_NEW TimedObject("TimedObject");
      { Timer y(z);
        int j=0;
        for (int i=0;i<INT_MAX;i++) j+=i;
      }
      z->printOn(cout);
    }
  }
  return EXIT_SUCCESS;
}
