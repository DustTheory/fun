#include <iostream>
#include <limits>
#include <stdlib.h>
#include "SetTraps.H"
#include "MemoryDebugger.H"
extern "C" {
  void sgetrf_simple__(const int&,float*,int*);
  void sgetrs_simple__(const int&,const float*,const int*,float*);
}
int main(int /*argc*/,char** /*argv*/) {
  { MemoryDebugger md(1);
    int n=2;

    float *A=new float[n*n];
    float *bx=new float[n];
    int *ipiv=new int[n];

    A[0]=numeric_limits<float>::epsilon()*0.25; A[1]=A[2]=A[3]=1.;
    bx[0]=1.; bx[1]=0.;

    sgetrf_simple__(n,A,ipiv);
    sgetrs_simple__(n,A,ipiv,bx);

    cout << "solution = " << bx[0] << " " << bx[1] << endl;

    delete ipiv; ipiv=0;
    delete bx; bx=0;
    delete A; A=0;
  }
}
