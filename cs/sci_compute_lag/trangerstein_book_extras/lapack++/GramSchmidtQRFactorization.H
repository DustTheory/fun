#ifndef GRAM_SCHMIDT_FACTORIZATION_H
#define GRAM_SCHMIDT_FACTORIZATION_H
//#include "LaEnum.H"
#include "Factorization.H"
#include "OrthogonalMatrix.H"
//#include "SymmetricMatrix.H"
#include "TrapezoidalMatrix.H"

template<typename F,typename Z> class GramSchmidtQRFactorization {
  private:
    static F tol;
    const Matrix<F,Z> *A_original;
    OrthogonalMatrix<F,Z> *Q;
    UpperTrapezoidalMatrix<F,Z> *R;
  protected:
    void reorthogonalize(int);
  public:
    GramSchmidtQRFactorization(const Matrix<F,Z>& A);
    GramSchmidtQRFactorization(const GramSchmidtQRFactorization &gsf) :
    A_original(gsf.A_original),Q(0),R(0) {
      int m=gsf.Q->size(0),n=gsf.Q->size(1);
      Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(m,n);
      Q->copyFrom('A',m,n,*gsf.Q);
      m=gsf.R->size(0),n=gsf.R->size(1);
      R=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(m,n);
      R->copyFrom(m,n,*gsf.R);
    }
    inline ~GramSchmidtQRFactorization() { 
      if (Q) delete Q; Q = 0;
      if (R) delete R; R=0;
    }

    inline bool isFactored() const { return Q != 0; }
    inline const OrthogonalMatrix<F,Z>& orthogonalPart() const {
      return *Q;
    }
    inline const UnitUpperTrapezoidalMatrix<F,Z>* rightTrapezoidalPart()
    const {
      UnitUpperTrapezoidalMatrix<F,Z> *Runit=OPERATOR_NEW
        UnitUpperTrapezoidalMatrix<F,Z>(R->size(0),R->size(1));
      Runit->copy(*R);
      return Runit;
    }

//  x = R^{-1} Q^T b:
    void solveOverdetermined(const Vector<F,Z> &b,Vector<F,Z> &x,
      Vector<F,Z> &residual) const;
    void solveOverdetermined(const Matrix<F,Z> &B,Matrix<F,Z> &X,
      Matrix<F,Z> &Residual) const;
//  x = Q R^{-T} b:
    void solveUnderdetermined(const Vector<F,Z> &b,Vector<F,Z> &x) const;
    void solveUnderdetermined(const Matrix<F,Z> &B,Matrix<F,Z> &X) const;

    void improveOverdetermined(const Vector<F,Z> &b,Vector<F,Z> &x,
      Vector<F,Z> &residual) const;
    void improveOverdetermined(const Matrix<F,Z> &B,Matrix<F,Z> &X,
      Matrix<F,Z> &Residual) const;
    void improveUnderdetermined(const Vector<F,Z> &b,Vector<F,Z> &x) const;
    void improveUnderdetermined(const Matrix<F,Z> &B,Matrix<F,Z> &X) const;

    void exchangeColumn(int j,int jAin,int jAout,const Matrix<F,Z> &A);
    void exchangeRow(int i,int iAin,int iAout,const Matrix<F,Z> &A);
    void dropColumn(int);
    void addColumn(int,const Matrix<F,Z>&);
    void dropRow(int);
    void dropRowAndColumn(int,int);
    void addRow(int,const Matrix<F,Z>&);
    void printOn(ostream& = cerr) const ;
};
template<typename F,typename Z> void testGramSchmidtQRFactorization(F,Z);

// Modified from ltgmd.H by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
#endif 
