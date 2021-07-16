#ifndef _GAUSSIAN_FACTORIZATION_H
#define _GAUSSIAN_FACTORIZATION_H
#include <iostream>
//#include "LaEnum.H"
#include "Factorization.H"
#include "SquareMatrix.H"
#include "TrapezoidalMatrix.H"

template<typename F,typename Z,class M> class GaussianFactorization {
  private:
//  if piv_op==pivot_rows, then ipiv is set in Fortran (1 -- m)
//  otherwise it is set in C++ (0 -- m-1)
    const M *A_original;
    M *LU;
    Vector<F,F> *r,*c;
    int *ipiv,*jpiv;
    Factorization::PIVOT_OPTION piv_op;
    Factorization::EQUILIBRATE_OPTION equ_op;
    F anormi,anormm,anormo;
      // infinity-norm, max entry norm, one-norm of equilibrated matrix
    F colcnd,rowcnd;
    F rpvgrw;
    char equed;
  public:
    GaussianFactorization(const M &A,
      Factorization::PIVOT_OPTION po=Factorization::PIVOT_ROWS,
      Factorization::EQUILIBRATE_OPTION=Factorization::NO_EQUILIBRATION);
      // dgetc2,dgetf2,dgetrf, dlatdf
    ~GaussianFactorization() {
      if (LU) delete LU; LU=0; 
      if (r) delete r; r=0;
      if (c) delete c; c=0;
      if (ipiv) delete [] ipiv;ipiv=0; 
      if (jpiv) delete [] jpiv;jpiv=0;
    }

//  const int* rowPivots() const { return ipiv; }
//  const int* columnPivots() const { return jpiv; }

    bool isFactored() const { return LU != 0; }
    double reciprocalPivotGrowth() const { return rpvgrw; } // dla_gerpvgrw

    void solve(const Vector<F,Z> &b,Vector<F,Z> &x,
      Factorization::TRANSPOSE_OPTION to=Factorization::NO_TRANSPOSE);
    void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,
      Factorization::TRANSPOSE_OPTION to=Factorization::NO_TRANSPOSE,
      Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE);
    F reciprocalConditionNumber(
      Factorization::CONDITION_NUMBER_NORM cnn=Factorization::ONE_NORM);
    void improve(const Vector<F,Z> &b,Vector<F,Z> &x,F &berr,F &ferr,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE);
    void improve(const Matrix<F,Z> &B,Matrix<F,Z> &X,Vector<F,F> &berr,
      Vector<F,F> &ferr,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE,
      Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE);
    void printOn(ostream& = cerr) const ;
};
template<typename F,typename Z> void testGaussianFactorization(F,Z);

// M = SquareMatrix==> UnitLowerTriangularMatrix, UpperTriangularMatrix
// M = UpperHessenbergMatrix==> UnitLowerBandMatrix, UpperTriangularMatrix
// M = TridiagonalMatrix==> BandMatrix,BandMatrix
// M = BandMatrix==> BandMatrix,BandMatrix

/*
template<typename F,typename Z> const UnitLowerTrapezoidalMatrix<F,Z>&
GaussianFactorization<F,Z,SquareMatrix<F,Z> >::lower() const {
  return *reinterpret_cast<const
    const UnitLowerTrapezoidalMatrix<F,Z>*>(LU);
}
template<typename F,typename Z,SquareMatrix<F,Z> >
const UpperTrapezoidalMatrix<F,Z>&
GaussianFactorization<F,Z,SquareMatrix<F,Z> >::upper() const {
  return *reinterpret_cast<const
    const UpperTrapezoidalMatrix<F,Z>*>(LU);
}
*/

// Modified from ltgmd.H by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
#endif 
