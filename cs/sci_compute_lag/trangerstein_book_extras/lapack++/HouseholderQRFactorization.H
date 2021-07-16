#ifndef HOUSEHOLDER_FACTORIZATION_H
#define HOUSEHOLDER_FACTORIZATION_H
#include "Factorization.H"
#include "OrthogonalMatrix.H"
#include "TrapezoidalMatrix.H"

template<typename F,typename Z> class HouseholderQRFactorization {
  private:
    const Matrix<F,Z> *A_original;
    Matrix<F,Z> *QR;
    Vector<F,Z> *tau;
    int *jpvt;
    Factorization::PIVOT_OPTION piv_op;
    int iascl;
    F ascl,anrm;
  public:
    HouseholderQRFactorization(const Matrix<F,Z>& A,
      Factorization::PIVOT_OPTION po=Factorization::NO_PIVOTING);
    ~HouseholderQRFactorization() {
      delete QR; QR=0; 
      delete tau; tau=0;
      if (jpvt!=0) delete [] jpvt; jpvt=0;
    }
    TrapezoidalMatrix<F,Z>* rightTrapezoidalPart() const {
      int m=QR->size(0),n=QR->size(1);
      if (m>=n) {
        UpperTrapezoidalMatrix<F,Z> *R=
          OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n);
        R->copyFrom(n,n,*QR);
        return R;
      }
      LowerTrapezoidalMatrix<F,Z> *L=
        OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(m,m);
      L->copyFrom(m,m,*QR);
      return L;
    }
    OrthogonalMatrix<F,Z>* orthogonalPart() const;
    F solve(const Vector<F,Z> &b,Vector<F,Z> &x,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,
      Vector<F,Z> &residual_norm,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void improve(const Vector<F,Z> &b,Vector<F,Z> &x,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void improve(const Matrix<F,Z> &B,Matrix<F,Z> &X,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void printOn(ostream& = cerr) const ;
};
template<typename F,typename Z> void testHouseholderQRFactorization(F,Z);

// Modified from ltgmd.H by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
#endif 
