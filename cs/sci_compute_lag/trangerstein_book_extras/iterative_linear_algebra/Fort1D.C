#include <iostream>
#include "Fort1D.H"

extern "C" {
  void F77NAME(print_loc)(void *p) {
    cout << p << endl;
    cout << "\tdouble = " << *(double*) p << endl;
    cout << "\tfloat = " << *(float*) p << endl;
    cout << "\tint = " << *(int*) p << endl;
  }
}
