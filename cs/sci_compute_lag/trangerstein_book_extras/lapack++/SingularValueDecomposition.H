#ifndef SINGULAR_VALUE_DECOMPOSITION_H
#define SINGULAR_VALUE_DECOMPOSITION_H
#include "Factorization.H"
#include "OrthogonalMatrix.H"
#include "TrapezoidalMatrix.H"

template<typename F,typename Z> class SingularValueDecomposition {
  private:
    const Matrix<F,Z> *A_original;
    OrthogonalMatrix<F,Z> *U;
    OrthogonalMatrix<F,Z> *Vtranspose;
    Vector<F,F> *s;
  public:
    SingularValueDecomposition(const Matrix<F,Z>& A);
    ~SingularValueDecomposition() {
      delete U; U=0; 
      delete Vtranspose; Vtranspose=0; 
      delete s; s=0;
    }
    OrthogonalMatrix<F,Z>* leftSingularVectors() const { return U; }
    OrthogonalMatrix<F,Z>* tranposedRightSingularVectors() const {
      return Vtranspose;
    }
    Vector<F,F>* singularValues() const { return s; }
    void solve(const Vector<F,Z> &b,Vector<F,Z> &x,F rcond,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,F rcond,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void regularize(const Vector<F,Z> &b,Vector<F,Z> &x,F ridge,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void regularize(const Matrix<F,Z> &B,Matrix<F,Z> &X,F ridge,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void printOn(ostream& = cerr) const ;
};
template<typename F,typename Z>
  void testSingularValueDecomposition(F,Z);

// Modified from ltgmd.H by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
#endif 
