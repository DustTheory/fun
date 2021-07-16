#ifndef _SYM_MAT_H
#define _SYM_MAT_H
#include <iostream>
#include <string.h>
#include "OrthogonalMatrix.H"
#include "SquareMatrix.H"
#include "TrapezoidalMatrix.H"

//for complex arrays, SymmetricMatrix is Hermetian
//only the lower triangular part is guaranteed to be stored correctly
template<typename F,typename Z> class SymmetricMatrix :
public SquareMatrix<F,Z> {
  private:
//  prevent some Matrix functions
    void interchangeColumns(int,int) { assert(0); }
    void interchangeRows(int,int) { assert(0); }
    void gemv(char trans,Z alpha,const Vector<F,Z> &x,Z beta,
    Vector<F,Z> &y) const {
      symv(alpha,x,beta,y);
    }
    void ger(Z alpha,const Vector<F,Z> &x,const Vector<F,Z> &y) {
      syr(alpha,x);
    }
    void gerc(Z alpha,const Vector<F,Z> &x,const Vector<F,Z> &y) {
      syr(alpha,x);
    }
    void geru(Z alpha,const Vector<F,Z> &x,const Vector<F,Z> &y) {
      syr(alpha,x);
    }
    void gemm(char transa,char transb,Z alpha,const Matrix<F,Z> &A,
    const Matrix<F,Z> &B,Z beta) {
      OBSOLETE("does not apply to this class");
    }
    virtual void copyFrom(char uplo,int m,int n,const Matrix<F,Z> &A) {
      SquareMatrix<F,Z>::copyFrom('L',n,n,A);
    }
    int lastNonzeroColumn() const { assert(0); return 0;}
    int lastNonzeroRow() const { assert(0); return 0;}
    F equilibrate(Vector<F,F> &r,Vector<F,F> &c,F &rowcnd,F &colcnd)
    const {
      return equilibrate(r,rowcnd);
    }
    void scale(char type,int kl,int ku,F cfrom,F cto) {
      scale(cfrom,cto);
    }
    void set(char uplo,Z diag,Z offdiag) { set(diag,offdiag); }

//  SquareMatrix<F,Z>* transpose() const { assert(0); return 0; }
//  SquareMatrix<F,Z>* conjugateTranspose() const { assert(0); return 0; }
//  prevent copy constructor
    explicit SymmetricMatrix(const SymmetricMatrix &S) :
      SquareMatrix<F,Z>(S.size(0)) { copy(S); }
//  use copy instead
    SymmetricMatrix& operator=(const SymmetricMatrix &S) {
      SquareMatrix<F,Z>::resize(S.size(0)); copy(S); return *this;
    }
  protected:
    virtual void fillWith(Z d);
  public:
    SymmetricMatrix() : SquareMatrix<F,Z>() {;}
    explicit SymmetricMatrix(int n) : SquareMatrix<F,Z>(n) {;}
    explicit SymmetricMatrix(int n,Z d) : SquareMatrix<F,Z>(n) {
      fillWith(d);
    }
    virtual ~SymmetricMatrix() {;}
    virtual SquareMatrix<F,Z>* makeMatrix() const;

    virtual SymmetricMatrix& operator=(Z scalar) {
      fillWith(scalar); return *this;
    }
    virtual Z* addr() { return SquareMatrix<F,Z>::addr(); }
    virtual const Z* addr() const { return SquareMatrix<F,Z>::addr(); }
    virtual Z* addr(int i,int j) {
      CHECK_TEST(j<=i)
      return SquareMatrix<F,Z>::addr(i,j);
    }
    virtual const Z* addr(int i,int j) const {
      CHECK_TEST(j<=i)
      return SquareMatrix<F,Z>::addr(i,j);
    }
    virtual Z operator()(int i,int j) const; // may return complex conjugate
//  require j<=i: for complex matrices, cannot return reference to conjugate
    virtual Z& operator()(int i,int j) { return *(this->addr(i,j)); }
    virtual void copy(const Matrix<F,Z> &M) {
      int n=this->size(0);
      CHECK_SAME(n,M.size(0));
      CHECK_SAME(n,M.size(1));
      Matrix<F,Z>::copyFrom('L',n,n,M);
    }

    SymmetricMatrix& operator+=(const SymmetricMatrix&);
    SymmetricMatrix& operator-=(const SymmetricMatrix&);
//  if Z is complex, scale by real part only:
    virtual SymmetricMatrix& operator*=(Z);
    virtual SymmetricMatrix& operator/=(Z);
    virtual SymmetricMatrix* operator+(const SymmetricMatrix &A) const {
      SymmetricMatrix *S=OPERATOR_NEW SymmetricMatrix(*this);
      *S+=A; return S;
    }
    virtual SquareMatrix<F,Z>* operator+(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator+(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator+(const SquareMatrix<F,Z> &A) const;
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &A) const;

    virtual SymmetricMatrix* operator-(const SymmetricMatrix &A) const {
      SymmetricMatrix *S=OPERATOR_NEW SymmetricMatrix(*this);
      *S-=A; return S;
    }
    virtual SquareMatrix<F,Z>* operator-(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator-(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator-(const SquareMatrix<F,Z> &A) const;
    virtual SquareMatrix<F,Z>* operator-(const Matrix<F,Z> &A) const;

    virtual SquareMatrix<F,Z>* operator*(Z scalar) const;
    virtual SquareMatrix<F,Z>* operator/(Z scalar) const;
    virtual SquareMatrix<F,Z>* operator*(const SymmetricMatrix<F,Z>&) const;
    virtual Matrix<F,Z>* operator*(const UpperTrapezoidalMatrix<F,Z>&)
      const;
    virtual Matrix<F,Z>* operator*(const LowerTrapezoidalMatrix<F,Z>&)
      const;
    virtual SquareMatrix<F,Z>* operator*(const SquareMatrix<F,Z>&) const;
    virtual Matrix<F,Z>* operator*(const Matrix<F,Z>&) const;
    virtual Vector<F,Z>* operator*(const Vector<F,Z>&) const;
//  Blas2 routines
    void symv(Z alpha,const Vector<F,Z> &x,Z beta,Vector<F,Z> &y) const;
    void syr(Z alpha,const Vector<F,Z> &x);
    void syr2(Z alpha,const Vector<F,Z> &x,const Vector<F,Z> &y);
//  Blas3 routines
    void symm(Z alpha,const Matrix<F,Z> &B,Z beta,Matrix<F,Z> &C,
      char side='L') const; // C=A*alpha*B+C*beta or C=B*alpha*A+C*beta
    void syrk(Z alpha,const Matrix<F,Z> &A,Z beta,char transa='N');
    void syr2k(Z alpha,const Matrix<F,Z> &A,const Matrix<F,Z> &B,Z beta,
      char transab='N');
//  LaPack routines
//  void syamv(F alpha,const Vector<F,Z> &x,F beta,Vector<F,F> &y) const;
    F equilibrate(Vector<F,F> &s,F &scond) const;
    void scale(F cfrom,F cto) {
      Matrix<F,Z>::scale('L',1,1,cfrom,cto);
    }
    void set(Z diag,Z offdiag) { Matrix<F,Z>::set('L',diag,offdiag); }
    virtual F normFrobenius() const;
    virtual F normInfinity() const;
    virtual F normMaxEntry() const;
    virtual F normOne() const;

    F reciprocalConditionNumber() const;
//  SymmetricMatrix* inverse() const;
    Vector<F,F>* eigenvalues(OrthogonalMatrix<F,Z> *&Q) const;

//  the final char argument is ignored:
    virtual void solve(const Matrix<F,Z>&,Matrix<F,Z>&,char side='L',
      char='N') const;
    virtual void solve(const Vector<F,Z>&,Vector<F,Z>&,char='N') const;

    virtual void printOn(ostream& = cerr) const;
};

template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const UnitUpperTrapezoidalMatrix<F,Z> &U,const SymmetricMatrix<F,Z> &S) {
  return S+U;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const UpperTrapezoidalMatrix<F,Z> &U,const SymmetricMatrix<F,Z> &S) {
  return S+U;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const UnitLowerTrapezoidalMatrix<F,Z> &L,const SymmetricMatrix<F,Z> &S) {
  return S+L;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const LowerTrapezoidalMatrix<F,Z> &L,const SymmetricMatrix<F,Z> &S) {
  return S+L;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const SquareMatrix<F,Z> &SQ,const SymmetricMatrix<F,Z> &S) {
  return S+SQ;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const Matrix<F,Z> &M,const SymmetricMatrix<F,Z> &S) {
  return S+M;
}

template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const UnitUpperTrapezoidalMatrix<F,Z>&,const SymmetricMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const UpperTrapezoidalMatrix<F,Z>&,const SymmetricMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const UnitLowerTrapezoidalMatrix<F,Z>&,const SymmetricMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const LowerTrapezoidalMatrix<F,Z>&,const SymmetricMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const SquareMatrix<F,Z>&,const SymmetricMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const Matrix<F,Z>&,const SymmetricMatrix<F,Z>&);

template<typename F,typename Z> Matrix<F,Z>* operator*(
const UnitUpperTrapezoidalMatrix<F,Z> &U,const SymmetricMatrix<F,Z> &S) {
  return static_cast<const UpperTrapezoidalMatrix<F,Z>&>(U) * S;
}
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const UpperTrapezoidalMatrix<F,Z>&,const SymmetricMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
const UnitLowerTrapezoidalMatrix<F,Z> &L,const SymmetricMatrix<F,Z> &S) {
  return static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L) * S;
}
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const LowerTrapezoidalMatrix<F,Z>&,const SymmetricMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SquareMatrix<F,Z>&,const SymmetricMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const Matrix<F,Z>&,const SymmetricMatrix<F,Z>&);

template<typename F,typename Z> void testSymmetricMatrix(F,Z);


template<typename F,typename Z> class SymmetricPositiveMatrix : 
public SymmetricMatrix<F,Z> {
  private:
    explicit SymmetricPositiveMatrix(const SymmetricPositiveMatrix &S) : 
      SymmetricMatrix<F,Z>(S.size(0)) { SymmetricMatrix<F,Z>::copy(S); }
    SymmetricPositiveMatrix& operator=(const SymmetricPositiveMatrix &S) {
      SquareMatrix<F,Z>::resize(S.size(0));
      SymmetricMatrix<F,Z>::copy(S); return *this;
    }

//  prevent some SymmetricMatrix functions
    virtual SymmetricPositiveMatrix& operator-=(
    const SymmetricMatrix<F,Z> &M) {
      OBSOLETE("inappropriate for this class"); return *this;
    }
    virtual SymmetricPositiveMatrix& operator-=(const Matrix<F,Z> &M) {
      OBSOLETE("inappropriate for this class"); return *this;
    }

  public:
    SymmetricPositiveMatrix() : SymmetricMatrix<F,Z>() {;}
    explicit SymmetricPositiveMatrix(int n): SymmetricMatrix<F,Z>(n) {;}
    explicit SymmetricPositiveMatrix(int n,Z d) :
      SymmetricMatrix<F,Z>(n) { SymmetricMatrix<F,Z>::fillWith(d); }
    virtual ~SymmetricPositiveMatrix() {;}

    virtual SymmetricPositiveMatrix& operator=(Z scalar) { 
      SymmetricMatrix<F,Z>::fillWith(scalar); return *this; 
    }
    SymmetricPositiveMatrix& operator+=(const SymmetricPositiveMatrix &S) {
      SymmetricMatrix<F,Z>::operator+=(S); return *this;
    }
//  scale by abs(Z):
    virtual SymmetricPositiveMatrix& operator*=(Z);
    virtual SymmetricPositiveMatrix& operator/=(Z);
    virtual SymmetricPositiveMatrix* operator+(
    const SymmetricPositiveMatrix &A) const {
      SymmetricPositiveMatrix *S=
        OPERATOR_NEW SymmetricPositiveMatrix(*this);
      *S+=A; return S;
    }
    virtual SymmetricMatrix<F,Z>* operator+(const SymmetricMatrix<F,Z> &S)
    const {
      return SymmetricMatrix<F,Z>::operator+(S);
    }
    virtual SquareMatrix<F,Z>* operator+(
    const UpperTrapezoidalMatrix<F,Z> &U) const {
      return SymmetricMatrix<F,Z>::operator+(U);
    }
    virtual SquareMatrix<F,Z>* operator+(
    const LowerTrapezoidalMatrix<F,Z> &L) const {
      return SymmetricMatrix<F,Z>::operator+(L);
    }
    virtual SquareMatrix<F,Z>* operator+(const SquareMatrix<F,Z> &S) const {
      return SymmetricMatrix<F,Z>::operator+(S);
    }
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &M) const {
      return SymmetricMatrix<F,Z>::operator+(M);
    }

//  LaPack routines:
    F equilibrate(Vector<F,F> &s,F &scond) const;
    F reciprocalConditionNumber() const;
//  SymmetricPositiveMatrix* inverse() const;
    virtual void solve(const Matrix<F,Z>&,Matrix<F,Z>&,char side='L',
      char='N') const;
    virtual void solve(const Vector<F,Z>&,Vector<F,Z>&,char='N') const;
    void printOn(ostream& = cerr) const;
};
template<typename F,typename Z> void testSymmetricPositiveMatrix(F,Z);

// Modified from gmd.H by John Trangenstein, 11/9/96 and 5/25/12
// gmd.H etc was unnecessarily repetitive; ie, not object-oriented.
#endif
