#include <iostream>
#include <limits> // for numeric_limits
#include <stdlib.h>
#include "SetTraps.H"
#include "MemoryDebugger.H"
extern "C" {
  void forward_sub__(int&,float*,float*);
  void back_sub__(int&,float*,float*);
  void gauss_elim__(int&,float*);
}
int main(int /*argc*/,char** /*argv*/) {
  { MemoryDebugger md(1);
    int n=2;

    float *A=new float[n*n];
    float *bx=new float[n];

//  1.-FLT_EPSILON < 1. and 1.-FLT_EPSILON*0.5 < 1., so...
    A[0]=numeric_limits<double>::epsilon()*0.25;
    A[1]=A[2]=A[3]=1.;
    bx[0]=1.;
    bx[1]=0.;

    gauss_elim__(n,A);
    forward_sub__(n,A,bx);
    back_sub__(n,A,bx);

//  float one=1.;
//  float onema=one-A[0];
//  cout << "1-(1-epsilon) = " << one-onema << endl;
    cout << "solution = " << bx[0] << " " << bx[1] << endl;
    delete A; A=0;
    delete bx; bx=0;
  }
}
