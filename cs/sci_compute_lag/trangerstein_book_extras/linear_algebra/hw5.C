#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "CholeskyFactorization.H"
#include "GaussianFactorization.H"
#include "MDMtFactorization.H"
#include "MemoryDebugger.H"
#include "SymmetricMatrix.H"
#include "Tracer.H"

using namespace std;

int main(int /*argc*/,char** /*argv*/) {
  MemoryDebugger md(1);
  {
    int n=6;
    SquareMatrix<float,float> A(n,0.);
    Vector<float,float> b(n);
    for (int j=0;j<n;j++) {
      for (int i=0;i<=j;i++) {
        A(i,j)=1./static_cast<float>(i+j+1);
        if (i<j) A(j,i)=A(i,j);
      }
      b[j]=1./static_cast<float>(j+n+1);
    }
    GaussianFactorization<float,float,SquareMatrix<float,float> >
      gf(A,Factorization::PIVOT_ROWS_AND_COLUMNS);
    float rcond_gf=gf.reciprocalConditionNumber(); // 1-norm by default
    Vector<float,float> x_gf(n);
    gf.solve(b,x_gf);
    Vector<float,float> r_gf(n);
    r_gf.copy(b);
    A.gemv(1.,x_gf,-1.,r_gf); // r=A*x-b
    float relative_error_in_x_gf=r_gf.asum()/(b.asum()*rcond_gf);// 1-norms
    cout << "\t|| r ||, || b ||, kappa = " << r_gf.asum() << " "
         << b.asum() << " " << 1./rcond_gf << endl;
    cout << "\trelative error in x_gf = " << relative_error_in_x_gf << endl;
    cout << "\tx_gf = " << endl;
    x_gf.printOn(cout);

    SymmetricMatrix<float,float> S(n);
    S.copy(A);
    MDMtFactorization<float,float> mdmt(S);
    float rcond_mdmt=mdmt.reciprocalConditionNumber();
    Vector<float,float> x_mdmt(n);
    mdmt.solve(b,x_mdmt);
    Vector<float,float> r_mdmt(n);
    r_mdmt.copy(b);
    A.gemv(1.,x_mdmt,-1.,r_mdmt);
    float relative_error_in_x_mdmt
      =r_mdmt.asum()/(b.asum()*rcond_mdmt);
    cout << "\n\t|| r ||, || b ||, kappa = " << r_mdmt.asum() << " "
         << b.asum() << " " << 1./rcond_mdmt << endl;
    cout << "\trelative error in x_mdmt = " << relative_error_in_x_mdmt
         << endl;
    cout << "\tx_mdmt = " << endl;
    x_mdmt.printOn(cout);

    SymmetricPositiveMatrix<float,float> P(n);
    P.copy(A);
    CholeskyFactorization<float,float,
      SymmetricPositiveMatrix<float,float> > chol(P);
    float rcond_chol=chol.reciprocalConditionNumber();
    Vector<float,float> x_chol(n);
    chol.solve(b,x_chol);
    Vector<float,float> r_chol(n);
    r_chol.copy(b);
    A.gemv(1.,x_chol,-1.,r_chol);
    float relative_error_in_x_chol
      =r_chol.asum()/(b.asum()*rcond_chol);
    cout << "\n\t|| r ||, || b ||, kappa = " << r_chol.asum() << " "
         << b.asum() << " " << 1./rcond_chol << endl;
    cout << "\trelative error in x_chol = " << relative_error_in_x_chol
         << endl;
    cout << "\tx_chol = " << endl;
    x_chol.printOn(cout);

//  float rcondP=P.reciprocalConditionNumber();
//  cout << "\n\trcond for SPM = " << rcondP << endl;
  }
}
