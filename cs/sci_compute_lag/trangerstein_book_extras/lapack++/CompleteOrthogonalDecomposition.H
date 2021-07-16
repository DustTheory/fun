#ifndef COMPLETE_ORTHOGONAL_DECOMPOSITION_H
#define COMPLETE_ORTHOGONAL_DECOMPOSITION_H
#include "Factorization.H"
#include "OrthogonalMatrix.H"
#include "TrapezoidalMatrix.H"

template<typename F,typename Z> class CompleteOrthogonalDecomposition {
  private:
    const Matrix<F,Z> *A_original;
    Matrix<F,Z> *URV;
    Vector<F,Z> *utau;
    Vector<F,Z> *vtau;
    int *jpvt;
    int iascl,the_rank;
    F ascl,anrm;
  public:
    CompleteOrthogonalDecomposition(const Matrix<F,Z>& A,F rcond);
    ~CompleteOrthogonalDecomposition() {
      delete URV; URV=0; 
      delete utau; utau=0;
      delete vtau; vtau=0;
      delete jpvt; jpvt=0;
    }
    int rank() const { return the_rank; }
    UpperTriangularMatrix<F,Z>* upperTriangularPart() const {
      UpperTriangularMatrix<F,Z> *R=
        OPERATOR_NEW UpperTriangularMatrix<F,Z>(the_rank,the_rank);
      R->copyFrom(the_rank,the_rank,*URV);
      return R;
    }
    OrthogonalMatrix<F,Z>* leftOrthogonalPart() const;
    OrthogonalMatrix<F,Z>* rightOrthogonalPartTransposed() const;
    F solve(const Vector<F,Z> &b,Vector<F,Z> &x,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,
      Vector<F,F> &residual_norm,
      Factorization::TRANSPOSE_OPTION=Factorization::NO_TRANSPOSE) const;
    void printOn(ostream& = cerr) const ;
};
template<typename F,typename Z>
  void testCompleteOrthogonalDecomposition(F,Z);

// Modified from ltgmd.H by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
#endif 
