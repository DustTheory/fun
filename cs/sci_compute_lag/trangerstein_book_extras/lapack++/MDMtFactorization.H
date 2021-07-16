#ifndef _MDMt_FACTORIZATION_H
#define _MDMt_FACTORIZATION_H
#include <iostream>
#include "Factorization.H"
#include "SymmetricMatrix.H"

// SymmetricMatrix: dsytrf
// condition number: dsycon
// improve: dsyrfs, dsyrfsx
// solve: dsysv, dsysvx, dsysvxx, dsytrs, dsytrs2
// factor: dsytf2, dsytrf
template<typename F,typename Z> class MDMtFactorization {
  private:
    const SymmetricMatrix<F,Z> *A_original;
    SymmetricMatrix<F,Z> *MD;
    int *ipiv;
    Vector<F,F> *s;
    Factorization::EQUILIBRATE_OPTION equ_op;
    F anormi,anormm,anormo;
      // infinity-norm, max entry norm, one-norm of equilibrated matrix
    F scond;
    F rpvgrw;
    char equed;
  public:
    MDMtFactorization(const SymmetricMatrix<F,Z> &A,
      Factorization::EQUILIBRATE_OPTION=Factorization::NO_EQUILIBRATION);
    ~MDMtFactorization() {
      if (MD) delete MD; MD=0; 
      if (s) delete s; s=0;
      if (ipiv) delete [] ipiv; ipiv=0;
    }

    bool isFactored() const { return MD != 0; }
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
template<typename F,typename Z> void testMDMtFactorization(F,Z);

#endif 
