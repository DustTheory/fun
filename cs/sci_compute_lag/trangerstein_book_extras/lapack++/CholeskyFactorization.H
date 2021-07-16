#ifndef _CHOLESKY_FACTORIZATION_H
#define _CHOLESKY_FACTORIZATION_H
#include <iostream>
#include "Factorization.H"
#include "Matrix.H"

// SymmetricPositiveMatrix: zpotrf ==> Cholesky
// SymmetricPositiveTridiagonalMatrix: zpttrf ==> L D L^T
// SymmetricPositiveBandMatrix : zpbtrf ==> Cholesky
template<typename F,typename Z,class M> class CholeskyFactorization {
  private:
    const M *A_original;
    M *L;
    Vector<F,F> *s;
    Factorization::EQUILIBRATE_OPTION equ_op;
    F anormi,anormm,anormo;
      // infinity-norm, max entry norm, one-norm of equilibrated matrix
    F scond;
    F rpvgrw;
    char equed;
  public:
    CholeskyFactorization(const M &A,
      Factorization::EQUILIBRATE_OPTION=Factorization::NO_EQUILIBRATION);
    ~CholeskyFactorization() {
      if (L) delete L; L=0; 
      if (s) delete s; s=0;
    }

    bool isFactored() const { return L != 0; }
    double reciprocalPivotGrowth() const { return rpvgrw; }

    void solve(const Vector<F,Z> &b,Vector<F,Z> &x);
    void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,
      Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE);
    F reciprocalConditionNumber();
    void improve(const Vector<F,Z> &b,Vector<F,Z> &x,F &berr,F &ferr);
    void improve(const Matrix<F,Z> &B,Matrix<F,Z> &X,Vector<F,F> &berr,
      Vector<F,F> &ferr,
      Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE);
    void printOn(ostream& = cerr) const ;
};
template<typename F,typename Z> void testCholeskyFactorization(F,Z);

#endif 
