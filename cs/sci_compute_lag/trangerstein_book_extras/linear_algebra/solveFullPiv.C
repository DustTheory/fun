#include <iostream>
#include <limits>
#include <stdlib.h>
#include "SetTraps.H"
#include "MemoryDebugger.H"
extern "C" {
  void forward_sub__(int&,float*,float*);
  void back_sub__(int&,float*,float*);
  void gauss_elim__(int&,float*,int*,int*,int&);
}
int main(int /*argc*/,char** /*argv*/) {
  { MemoryDebugger md(1);
    int n=2;

    float *A=new float[n*n];
    float *bx=new float[n];
    int *ipiv=new int[n];
    int *jpiv=new int[n];

    A[0]=numeric_limits<float>::epsilon()*0.25;
    A[1]=A[2]=A[3]=1.;
    bx[0]=1.;
    bx[1]=0.;

//  factor matrix with pivoting
    int npivots=0;
    gauss_elim__(n,A,ipiv,jpiv,npivots);
//  reorder equations
    for (int i=0;i<n;i++) {
      int ipiv1=ipiv[i]-1;
      if (i != ipiv1) { // Fortran indices start at 1
        float temp=bx[i]; bx[i]=bx[ipiv1]; bx[ipiv1]=temp;
      }
    }
//  forward- and back-solve
    forward_sub__(n,A,bx);
    back_sub__(n,A,bx);
//  reorder unknowns
    for (int i=n-1;i>=0;i--) {
      int jpiv1=jpiv[i]-1;
      if (i != jpiv1) { // Fortran indices start at 1
        float temp=bx[i]; bx[i]=bx[jpiv1]; bx[jpiv1]=temp;
      }
    }

    cout << "solution = " << bx[0] << " " << bx[1] << endl;

    delete jpiv; jpiv=0;
    delete ipiv; ipiv=0;
    delete bx; bx=0;
    delete A; A=0;
  }
}
