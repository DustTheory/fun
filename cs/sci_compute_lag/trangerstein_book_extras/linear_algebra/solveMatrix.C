#include <iostream>
#include <limits>
#include <stdlib.h>
#include "SetTraps.H"
#include "MemoryDebugger.H"
#include "matrix.H"
#include "vector.H"
int main(int /*argc*/,char** /*argv*/) {
  { MemoryDebugger md(1);
    int n=2;

    { Matrix A(n,n);
      Vector b(n);

      A(0,0)=numeric_limits<float>::epsilon()*.25; A(1,0)=A(0,1)=A(1,1)=1.;
      b[0]=1.; b[1]=0.;

      Vector *x=A.solve(b);
      cout << "solution = " << *x << endl;

      delete x;
    } // A and bx go out of scope here
  }
}
