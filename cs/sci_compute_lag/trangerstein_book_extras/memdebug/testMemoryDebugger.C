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
// "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/testMemoryDebugger.C,v 1.1 2009/08/20 17:33:33 johnt Exp $"
#include <iostream>
#include <stdlib.h>
#include "Arch.H"
#include "MemoryDebugger.H"
#include "Tracer.H"
extern "C" {
  void F77NAME(fort)(int &,double* );
}

class IntAlloc {
  private:
    int *data;
  public:
    IntAlloc() { data=OPERATOR_NEW int; }
    ~IntAlloc() { delete data; }
};

class IntArray {
  private:
    int *array;
  public:
    IntArray(int n) { array=OPERATOR_NEW_BRACKET(int,n); }
    ~IntArray() { delete [] array; }
};

int main(int /*argc*/,char** /*argv[]*/) {
  cout << boolalpha;
  cout << "entering main" << endl;

//array a is allocated outside the scope of the MemoryDebugger,
//  so we cannot catch errors in its use.
  int *a=OPERATOR_NEW_BRACKET(int,10);

  size_t buffer_size=1;
  {
    MemoryDebugger md(buffer_size);

//  the OPERATOR_NEW macro will store the file name and line number for
//    the next memory allocation, and print those with error messages
    int *b=OPERATOR_NEW int;

//  If the following 3 lines are uncommented, the MemoryDebugger will
//    generate an error message and abort
//  Both Tracer::Tracer and Tracer::~Tracer will call mem_check
//  { TRACER_CALL(t,"underrunning integer array");
//    b[-1]=0; 
//  }

//  the following checks a specific pointer for out-of_bounds errors
//  mem_check_ptr(b);

//  the following checks all pointers for out-of_bounds errors
//  mem_check();

//  the following may not generate an error, because sizeof(int) may be 
//    less than the alignment size
//  MemoryDebugger can write special bit patterns into memory only at
//    aligned locations
    { TRACER_CALL(t,"overrunning integer array");
      b[1]=0;
    }

//  the following will generate an error
//  { TRACER_CALL(t,"overrunning integer array");
//     b[2]=0;
//  }

//  the following will not generate an error, because a was allocated 
//    before the MemoryDebugger was constructed
    { TRACER_CALL(t,"overrunning array allocated outside MemoryDebugger scope");
       a[1]=0;
    }

//  the following will generate an error, because the underrun 
//    overwrites memory used by malloc
//  { TRACER_CALL(t,"underrunning array allocated outside MemoryDebugger scope");
//     a[-1]=0;
//  }

    double *c=OPERATOR_NEW_BRACKET(double,5);

//  If the following 3 lines are uncommented, the MemoryDebugger will
//    generate an error message and abort
//  { TRACER_CALL(t,"underrunning double array");
//    c[-1]=1.;
//  }

//  mem_check_ptr(c);
//  mem_check();

//  If the following 3 lines are uncommented, the MemoryDebugger will
//    generate an error message and abort
//  { TRACER_CALL(t,"overrunning double array");
//    c[5]=-1.;
//  }

//  mem_check_ptr(c);
//  mem_check();

//  if b was underrun or overrun, the following will generate an error 
//    message
//  { TRACER_CALL(t,"deleting integer array");
      delete b;
//  }

//  if c was underrun or overrun, the following will generate an error 
//    message
//  { TRACER_CALL(t,"deleting double array");
      delete [] c;
//  }
//  mem_check();
    
  
//  if not allocated using OPERATOR_NEW or OPERATOR_NEW_BRACKET,
//  MemoryDebugger will not know file name and line number for call to new
    double *d=new double [3];
//  { TRACER_CALL(t,"calling Fortran with correct array dimension");
      int n=3;
      fort_(n,d);
//  }

//  uncommenting the following 3 lines will cause a write out of bounds
//    and an error message
//  { TRACER_CALL(t,"calling Fortran with too-large array dimension");
//      int m=4;
//      fort_(m,d);
//  }

//  the following 3 lines will not cause MemoryDebugger to abort
//    because the write out of bounds occurs beyond the buffer size
//  however, the code may get a segmentation fault as a result of the 
//    write
//  { TRACER_CALL(t,"calling Fortran to write beyond buffer size");
//    int m=4+buffer_size;
//    fort_(m,d);
//  }

//  commenting any of the deletes inside the scope of the MemoryDebugger
//    will generate an error message
    delete [] d;

//  check that operator new and operator delete call 
//    matching constructors and destructors
    IntAlloc *ial=OPERATOR_NEW IntAlloc();
    delete ial;

    IntArray *iar=OPERATOR_NEW IntArray(2);
    delete iar;
  }
  delete [] a;
  cout << "leaving main" << endl;
  return EXIT_SUCCESS;
}
