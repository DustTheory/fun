#ifndef _BAND_MATRIX_H
#define _BAND_MATRIX_H
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "OrthogonalMatrix.H"
#include "SquareMatrix.H"
#include "SymmetricMatrix.H"
#include "TrapezoidalMatrix.H"

template<typename F,typename Z> class BandMatrix;
template<typename F,typename Z> class UpperHessenbergMatrix;

//we cannot derive this from either a Matrix or a BandMatrix because
//  LaPack uses very different storage schemes for each
template<typename F,typename Z> class TridiagonalMatrix {
  private:
    Vector<F,Z> *L,*D,*U;
    int dim;

//  prevent public use of copy constructor
    explicit TridiagonalMatrix(const TridiagonalMatrix &T) : dim(T.dim) {
      L=OPERATOR_NEW Vector<F,Z>(dim-1); L->copy(*T.L);
      D=OPERATOR_NEW Vector<F,Z>(dim);   D->copy(*T.D);
      U=OPERATOR_NEW Vector<F,Z>(dim-1); U->copy(*T.U);
    }
// Use copy instead
    TridiagonalMatrix& operator=(const TridiagonalMatrix& T) {
      resize(T); copy(T); return *this;
    }
  protected:
    void fillWith(Z scalar) { *L=scalar; *D=scalar; *U=scalar; }
  public:
    static const Z outofbounds_;
    static const Z undefined_;
    static Z safety_;

    TridiagonalMatrix() : L(0),D(0),U(0),dim(0) {;}

    explicit TridiagonalMatrix(int n): dim(n) { 
      CHECK_TEST(n>=2);
      L=OPERATOR_NEW Vector<F,Z>(n-1);
      D=OPERATOR_NEW Vector<F,Z>(n);
      U=OPERATOR_NEW Vector<F,Z>(n-1);
    }
    explicit TridiagonalMatrix(int n,Z d): dim(n) {
      CHECK_TEST(n>=2);
      L=OPERATOR_NEW Vector<F,Z>(n-1,d);
      D=OPERATOR_NEW Vector<F,Z>(n,d);
      U=OPERATOR_NEW Vector<F,Z>(n-1,d);
    }
    virtual ~TridiagonalMatrix() {
      delete L; L=0; delete D; D=0; delete U; U=0;
    }
    virtual SquareMatrix<F,Z>* makeMatrix() const;
    virtual TridiagonalMatrix& operator=(Z scalar) {
      *L=scalar; *D=scalar; *U=scalar; return *this;
    }

    int size(int) const { return dim; }
    virtual Z* lowerDiagonalAddr() const { return L->addr(); }
    virtual Z* diagonalAddr() const { return D->addr(); }
    virtual Z* upperDiagonalAddr() const { return U->addr(); }
    virtual Z* addr(int i,int j) {
      if (i==j) return D->addr(i);
      if (i==j-1) return U->addr(i);
      if (i==j+1) return L->addr(j);
      OBSOLETE("invalid array indices");
      return 0;
    }
    virtual const Z* addr(int i,int j) const {
      if (i==j) return D->addr(i);
      if (i==j-1) return U->addr(i);
      if (i==j+1) return L->addr(j);
      OBSOLETE("invalid array indices");
      return 0;
    }
    virtual Z& operator()(int i,int j) { 
      if (i==j) return (*D)[i];
      if (i==j-1) return (*U)[i];
      if (i==j+1) return (*L)[j];
      return safety_=outofbounds_;
    }
    virtual Z operator()(int i,int j) const { 
      if (i==j) return (*D)[i];
      if (i==j-1) return (*U)[i];
      if (i==j+1) return (*L)[j];
      return outofbounds_;
    }
    virtual void resize(int n) { 
      if (n==dim) return;
      CHECK_TEST(n>=2);
      if (L) { delete L; L=0; delete D; D=0; delete U; U=0; }
      L=OPERATOR_NEW Vector<F,Z>(n-1);
      D=OPERATOR_NEW Vector<F,Z>(n);
      U=OPERATOR_NEW Vector<F,Z>(n-1);
      dim=n;
    }
    virtual void resize(const TridiagonalMatrix<F,Z>& T) { resize(T.dim); }
    virtual void copy(const TridiagonalMatrix<F,Z>& T) {
      CHECK_SAME(dim,T.dim);
      L->copy(*T.L); D->copy(*T.D); U->copy(*T.U);
    }
    TridiagonalMatrix& operator+=(const TridiagonalMatrix<F,Z> &T) {
      CHECK_SAME(dim,T.dim)
      (*L) += (*T.L); (*D) += (*T.D); (*U) += (*T.U);
      return *this;
    }
    TridiagonalMatrix& operator-=(const TridiagonalMatrix<F,Z> &T) {
      CHECK_SAME(dim,T.dim)
      (*L) -= (*T.L); (*D) -= (*T.D); (*U) -= (*T.U); return *this;
    }
    virtual TridiagonalMatrix& operator*=(Z d) {
      (*L) *= d; (*D) *= d; (*U) *= d; return *this;
    }
    virtual TridiagonalMatrix& operator/=(Z d) {
      (*L) /= d; (*D) /= d; (*U) /= d; return *this;
    }
    virtual TridiagonalMatrix* operator+(const TridiagonalMatrix &T) const
    {
      TridiagonalMatrix *sum=OPERATOR_NEW TridiagonalMatrix(dim);
      sum->copy(*this); (*sum)+=T; return sum;
    }
    virtual SquareMatrix<F,Z>* operator+(
      const SymmetricMatrix<F,Z> &S) const;
    virtual SquareMatrix<F,Z>* operator+(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual UpperHessenbergMatrix<F,Z>* operator+(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &M) const;

    virtual TridiagonalMatrix* operator-(const TridiagonalMatrix &T) const
    {
      TridiagonalMatrix *dif=OPERATOR_NEW TridiagonalMatrix(dim);
      dif->copy(*this); (*dif)-=T; return dif;
    }
    virtual SquareMatrix<F,Z>* operator-(const SymmetricMatrix<F,Z> &M)
      const;
    virtual SquareMatrix<F,Z>* operator-(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual UpperHessenbergMatrix<F,Z>* operator-(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator-(const Matrix<F,Z> &M) const;
    virtual TridiagonalMatrix* operator*(Z d) const {
      TridiagonalMatrix *prod=OPERATOR_NEW TridiagonalMatrix(dim);
      prod->copy(*this); (*prod) *= d; return prod;
    }
    virtual TridiagonalMatrix* operator/(Z d) const {
      TridiagonalMatrix *prod=OPERATOR_NEW TridiagonalMatrix(dim);
      prod->copy(*this); (*prod) /= d; return prod;
    }

    virtual BandMatrix<F,Z>* operator*(const TridiagonalMatrix &T) const;
    virtual SquareMatrix<F,Z>* operator*(const SymmetricMatrix<F,Z> &S)
      const;
    virtual Matrix<F,Z>* operator*(const LowerTrapezoidalMatrix<F,Z> &L)
      const;
    virtual Matrix<F,Z>* operator*(const UpperTrapezoidalMatrix<F,Z> &U)
      const;
    virtual SquareMatrix<F,Z>* operator*(const SquareMatrix<F,Z> &M) const;
    virtual Matrix<F,Z>* operator*(const Matrix<F,Z> &M) const;
    virtual Vector<F,Z>* operator*(const Vector<F,Z> &v) const;

// LaPack:
    F normFrobenius() const;
    F normInfinity() const;
    F normMaxEntry() const;
    F normOne() const;
    F reciprocalConditionNumber(char norm) const; // dgtcon
    void gtmv(Z alpha,const Vector<F,Z> &x,Z beta,Vector<F,Z> &b,
      char trans='N') const; // b=T*x*alpha+b*beta or b=T^T*x*alpha+b*beta
    void gtmm(Z alpha,const Matrix<F,Z> &X,Z beta,Matrix<F,Z> &B,
      char side='L',char trans='N') const;
      // B=T*X*alpha+B*beta or B=T^T*X*alpha+B*beta or
      // B=X*T*alpha+B*beta or B=X*T^T*alpha+B*beta
//  to find eigenvalues and eigenvectors, construct UpperHessenbergMatrix

    virtual void solve(const Vector<F,Z> &b,Vector<F,Z> &x,char trans='N')
      const; // dgtts2, dgttrs, dgtsv, dgtsvx, dgttrf
    virtual void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,char side='L',
      char trans='N') const;
    virtual void printOn(ostream& =cerr) const ;

};
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const SymmetricMatrix<F,Z> &S,const TridiagonalMatrix<F,Z> &T) {
  return T+S;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const LowerTrapezoidalMatrix<F,Z> &L,const TridiagonalMatrix<F,Z> &T) {
  return T+L;
}
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator+(
const UpperTrapezoidalMatrix<F,Z> &U,const TridiagonalMatrix<F,Z> &T) {
  return T+U;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const Matrix<F,Z> &M,const TridiagonalMatrix<F,Z> &T) {
  return T+M;
}

template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const SymmetricMatrix<F,Z>&,const TridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const LowerTrapezoidalMatrix<F,Z>&,const TridiagonalMatrix<F,Z>&);
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator-(
  const UpperTrapezoidalMatrix<F,Z>&,const TridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const Matrix<F,Z>&,const TridiagonalMatrix<F,Z>&);

template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SymmetricMatrix<F,Z>&,const TridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const LowerTrapezoidalMatrix<F,Z>&,const TridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const UpperTrapezoidalMatrix<F,Z>&,const TridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SquareMatrix<F,Z>&,const TridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const Matrix<F,Z>&,const TridiagonalMatrix<F,Z>&);

template<typename F,typename Z> void testTridiagonalMatrix(F,Z);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//we cannot derive this from a BandMatrix or a TridiagonalMatrix
//  because LaPack uses very different storage schemes for these three

template<typename F,typename Z>
class SymmetricTridiagonalMatrix {
  private:
//  prevent public use of copy constructor
    explicit SymmetricTridiagonalMatrix(
    const SymmetricTridiagonalMatrix &T) : dim(T.dim) {
      L=OPERATOR_NEW Vector<F,Z>(dim-1); L->copy(*T.L);
      D=OPERATOR_NEW Vector<F,F>(dim);   D->copy(*T.D);
    }
// Use copy instead
    SymmetricTridiagonalMatrix& operator=(
    const SymmetricTridiagonalMatrix& T) {
      resize(T); copy(T); return *this;
    }
  protected:
    Vector<F,F> *D;
    Vector<F,Z> *L;
    int dim;

    void fillWith(Z scalar);
  public:
    static const Z outofbounds_;
    static const Z undefined_;
    static Z safety_;

    SymmetricTridiagonalMatrix() : L(0),D(0),dim(0) {;}
    explicit SymmetricTridiagonalMatrix(int n): dim(n) { 
      L=OPERATOR_NEW Vector<F,Z>(n-1);
      D=OPERATOR_NEW Vector<F,F>(n);
    }
    explicit SymmetricTridiagonalMatrix(int n,Z d);
    virtual ~SymmetricTridiagonalMatrix() { 
      delete L; L=0; delete D; D=0;
    }
    virtual SquareMatrix<F,Z>* makeMatrix() const;
    virtual SymmetricTridiagonalMatrix& operator=(Z scalar);
    virtual void resize(int n) { 
      if (n==dim) return;
      CHECK_TEST(n>=2);
      if (L) { delete L; L=0; delete D; D=0; }
      L=OPERATOR_NEW Vector<F,Z>(n-1);
      D=OPERATOR_NEW Vector<F,F>(n);
      dim=n;
    }
    virtual void resize(const SymmetricTridiagonalMatrix<F,Z>& T) { 
      resize(T.size(0)); 
    }
    virtual void copy(const SymmetricTridiagonalMatrix<F,Z>& T) {
      CHECK_SAME(dim,T.dim);
      L->copy(*T.L); D->copy(*T.D);
    }

    int size(int ) const { return dim; }
    virtual const Vector<F,F>* diagonal() const { return D; }
    virtual const Vector<F,Z>* lowerDiagonal() const { return L; }
    virtual Z* lowerDiagonalAddr(int i) const { return L->addr(i); }
    virtual F* diagonalAddr(int i) const { return D->addr(i); }
    virtual Z lowerDiagonalValue(int i) const { return (*L)[i]; }
    virtual Z& lowerDiagonalValue(int i) { return (*L)[i]; }
    virtual F diagonalValue(int i) const { return (*D)[i]; }
    virtual F& diagonalValue(int i) { return (*D)[i]; }
    virtual Z upperDiagonalValue(int i) const;
    virtual Z operator()(int i,int j) const {
      if (i==j) return diagonalValue(i);
      if (i==j-1) return upperDiagonalValue(i);
      if (j==i-1) return lowerDiagonalValue(j);
      return outofbounds_;
    }

    SymmetricTridiagonalMatrix& operator+=(
    const SymmetricTridiagonalMatrix &T) {
      CHECK_SAME(dim,T.dim)
      (*L) += (*T.L); (*D) += (*T.D); return *this;
    }
    SymmetricTridiagonalMatrix& operator-=(
    const SymmetricTridiagonalMatrix &T) {
      CHECK_SAME(dim,T.dim)
      (*L) -= (*T.L); (*D) -= (*T.D); return *this;
    }
    virtual SymmetricTridiagonalMatrix& operator*=(Z d);
    virtual SymmetricTridiagonalMatrix& operator/=(Z d);
    virtual SymmetricTridiagonalMatrix* operator+(
    const SymmetricTridiagonalMatrix &T) const {
      SymmetricTridiagonalMatrix *S=
        OPERATOR_NEW SymmetricTridiagonalMatrix(dim);
      S->copy(*this); (*S)+=T; return S;
    }
    virtual TridiagonalMatrix<F,Z>* operator+(
      const TridiagonalMatrix<F,Z> &T) const;
    virtual SymmetricMatrix<F,Z>* operator+(
      const SymmetricMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator+(
      const LowerTrapezoidalMatrix<F,Z> &T) const;
    virtual UpperHessenbergMatrix<F,Z>* operator+(
      const UpperTrapezoidalMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &T) const;

    virtual SymmetricTridiagonalMatrix<F,Z>* operator-(
      const SymmetricTridiagonalMatrix &S) const;
    virtual TridiagonalMatrix<F,Z>* operator-(
      const TridiagonalMatrix<F,Z> &T) const;
    virtual SymmetricMatrix<F,Z>* operator-(
      const SymmetricMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator-(
      const LowerTrapezoidalMatrix<F,Z> &T) const;
    virtual UpperHessenbergMatrix<F,Z>* operator-(
      const UpperTrapezoidalMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator-(const Matrix<F,Z> &T) const;

    virtual TridiagonalMatrix<F,Z>* operator*(Z d) const;
    virtual TridiagonalMatrix<F,Z>* operator/(Z d) const;

    virtual BandMatrix<F,Z>* operator*(
      const SymmetricTridiagonalMatrix &S) const;
    virtual BandMatrix<F,Z>* operator*(
      const TridiagonalMatrix<F,Z> &S) const;
    virtual SquareMatrix<F,Z>* operator*(
      const SymmetricMatrix<F,Z> &S) const;
    virtual Matrix<F,Z>* operator*(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual Matrix<F,Z>* operator*(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator*(const SquareMatrix<F,Z> &S) const;
    virtual Matrix<F,Z>* operator*(const Matrix<F,Z> &M) const;
    virtual Vector<F,Z>* operator*(const Vector<F,Z> &v) const;

// LaPack:
    F normFrobenius() const;
    F normInfinity() const;
    F normMaxEntry() const;
    F normOne() const;
    F reciprocalConditionNumber(char norm) const;
    void stmv(Z alpha,const Vector<F,Z> &x,Z beta,Vector<F,Z> &b) const;
      // y=T*x*alpha+y*beta
    void stmm(Z alpha,const Matrix<F,Z> &X,Z beta,Matrix<F,Z> &B,
      char side='L') const; // C=T*B*alpha+C*beta or C=B*T*alpha+C*beta

    virtual void solve(const Vector<F,Z> &b,Vector<F,Z> &x) const;
    virtual void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,char side='L')
      const;

    void printOn(ostream& =cerr) const ;

};

// to compute eigenvalues of complex SymmetricTridiagonalMatrix,
//   construct a SymmetricBandMatrix
// zsteqr works only on real SymmetricTridiagonalMatrix
template<typename F> Vector<F,F>* eigenvalues(
  const SymmetricTridiagonalMatrix<F,F> &T,OrthogonalMatrix<F,F> *&Q);

template<typename F,typename Z> TridiagonalMatrix<F,Z>* operator+(
const TridiagonalMatrix<F,Z> &T,
const SymmetricTridiagonalMatrix<F,Z> &S) {
  return S+T;
}
template<typename F,typename Z> SymmetricMatrix<F,Z>* operator+(
const SymmetricMatrix<F,Z> &M,
const SymmetricTridiagonalMatrix<F,Z> &T) {
  return T+M;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const LowerTrapezoidalMatrix<F,Z> &L,
const SymmetricTridiagonalMatrix<F,Z> &T) {
  return T+L;
}
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator+(
const UpperTrapezoidalMatrix<F,Z> &U,
const SymmetricTridiagonalMatrix<F,Z> &T) {
  return T+U;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const Matrix<F,Z> &M,const SymmetricTridiagonalMatrix<F,Z> &T) {
  return T+M;
}

template<typename F,typename Z> TridiagonalMatrix<F,Z>* operator-(
  const TridiagonalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SymmetricMatrix<F,Z>* operator-(
  const SymmetricMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const LowerTrapezoidalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator-(
  const UpperTrapezoidalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const Matrix<F,Z>&,const SymmetricTridiagonalMatrix<F,Z>&);

template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const TridiagonalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SymmetricMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const LowerTrapezoidalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const UpperTrapezoidalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const Matrix<F,Z>&,const SymmetricTridiagonalMatrix<F,Z>&);

template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const TridiagonalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SymmetricMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const LowerTrapezoidalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const UpperTrapezoidalMatrix<F,Z>&,
  const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SquareMatrix<F,Z>&,const SymmetricTridiagonalMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const Matrix<F,Z>&,const SymmetricTridiagonalMatrix<F,Z>&);

template<typename F,typename Z> void
  testSymmetricTridiagonalMatrix(F,Z);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z>
class SymmetricPositiveTridiagonalMatrix :
public SymmetricTridiagonalMatrix<F,Z> {
  private:
//  prevent public use of copy constructor
    explicit SymmetricPositiveTridiagonalMatrix(
    const SymmetricPositiveTridiagonalMatrix &T) :
    SymmetricTridiagonalMatrix<F,Z>() {
      OBSOLETE(0);
    }
// Use copy instead
    SymmetricPositiveTridiagonalMatrix& operator=(
    const SymmetricTridiagonalMatrix<F,Z>& T) {
      SymmetricTridiagonalMatrix<F,Z>::resize(T); copy(T); return *this;
    }
  public:
    SymmetricPositiveTridiagonalMatrix() :
      SymmetricTridiagonalMatrix<F,Z>() {}
    explicit SymmetricPositiveTridiagonalMatrix(int n): 
      SymmetricTridiagonalMatrix<F,Z>(n) {}
    explicit SymmetricPositiveTridiagonalMatrix(int n,Z d): 
      SymmetricTridiagonalMatrix<F,Z>(n,d) {}
    virtual ~SymmetricPositiveTridiagonalMatrix() {} 

    virtual SymmetricPositiveTridiagonalMatrix& operator=(Z scalar) {
      SymmetricTridiagonalMatrix<F,Z>::operator=(scalar); return *this;
    }

    virtual void copy(const SymmetricTridiagonalMatrix<F,Z>& T) {
      CHECK_NONZERO(dynamic_cast<const
        SymmetricPositiveTridiagonalMatrix*>(&T));
      SymmetricTridiagonalMatrix<F,Z>::copy(T);
    }
    SymmetricPositiveTridiagonalMatrix& operator+=(
    const SymmetricPositiveTridiagonalMatrix &T) {
      SymmetricTridiagonalMatrix<F,Z>::operator+=(T); return *this;
    }

    virtual SymmetricPositiveTridiagonalMatrix* operator+(
    const SymmetricPositiveTridiagonalMatrix &T) const {
      SymmetricPositiveTridiagonalMatrix *S=OPERATOR_NEW
        SymmetricPositiveTridiagonalMatrix(this->size(0));
      S->copy(*this); (*S)+=T; return S;
    }
    virtual SymmetricTridiagonalMatrix<F,Z>* operator+(
    const SymmetricTridiagonalMatrix<F,Z> &T) const {
      SymmetricTridiagonalMatrix<F,Z> *S=
        OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(this->size(0));
      S->copy(*this); (*S)+=T; return S;
    }
    virtual TridiagonalMatrix<F,Z>* operator+(
    const TridiagonalMatrix<F,Z> &T) const {
      return SymmetricTridiagonalMatrix<F,Z>::operator+(T);
    }
    virtual SymmetricPositiveMatrix<F,Z>* operator+(
      const SymmetricPositiveMatrix<F,Z> &S) const;
    virtual SymmetricMatrix<F,Z>* operator+(
    const SymmetricMatrix<F,Z> &S) const {
      return SymmetricTridiagonalMatrix<F,Z>::operator+(S);
    }
    virtual SquareMatrix<F,Z>* operator+(
    const LowerTrapezoidalMatrix<F,Z> &L) const {
      return SymmetricTridiagonalMatrix<F,Z>::operator+(L);
    }
    virtual UpperHessenbergMatrix<F,Z>* operator+(
    const UpperTrapezoidalMatrix<F,Z> &U) const {
      return SymmetricTridiagonalMatrix<F,Z>::operator+(U);
    }
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &M) const {
      return SymmetricTridiagonalMatrix<F,Z>::operator+(M);
    }

// LaPack:
    F reciprocalConditionNumber(char norm) const; // dptcon

    virtual void solve(const Vector<F,Z> &b,Vector<F,Z> &x) const;
    virtual void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,char side='L')
      const; // dptts2, dpttrf, dpttrs, dptsv, dptsvx

    void printOn(ostream& =cerr) const ;

};
template<typename F,typename Z> void
  testSymmetricPositiveTridiagonalMatrix(F,Z);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> class DiagonalMatrix {
  private:
    Vector<F,Z> *D;
    int dim;

//  prevent public use of copy constructor
    explicit DiagonalMatrix(const DiagonalMatrix &M) : dim(M.dim) {
      D=OPERATOR_NEW Vector<F,Z>(dim);   D->copy(*M.D);
    }
// Use copy instead
    DiagonalMatrix& operator=(const DiagonalMatrix &M) {
      resize(M); copy(M); return *this;
    }
  protected:
    virtual void fillWith(Z scalar) { *D=scalar; }
  public:
    DiagonalMatrix() : D(0),dim(0) {;}
    explicit DiagonalMatrix(int n) : dim(n) {
      CHECK_TEST(n>=1);
      D=OPERATOR_NEW Vector<F,Z>(n);
    }
    explicit DiagonalMatrix(int n,Z d) : dim(n) {
      CHECK_TEST(n>=1);
      D=OPERATOR_NEW Vector<F,Z>(n,d);
    }
    virtual ~DiagonalMatrix() {
      delete D; D=0;
    }
    virtual SquareMatrix<F,Z>* makeMatrix() const;
    virtual DiagonalMatrix& operator=(Z scalar) {
      *D=scalar; return *this;
    }

    int size(int) const { return dim; }
    virtual Z* addr(int i=0) { return D->addr(i); }
    virtual const Z* addr(int i=0) const { return D->addr(i); }
    virtual Z& operator[](int i) { return (*D)[i]; }
    virtual Z operator[](int i) const { return (*D)[i]; }
    virtual void resize(int n) {
      if (n==dim) return;
      CHECK_TEST(n>=1);
      if (D) { delete D; D=0; }
      D=OPERATOR_NEW Vector<F,Z>(n);
      dim=n;
    }
    virtual void resize( const DiagonalMatrix &M) { resize(M.dim); }
    virtual void copy(const DiagonalMatrix &M) {
      CHECK_SAME(dim,M.dim);
      D->copy(*M.D);
    }
    DiagonalMatrix& operator+=(const DiagonalMatrix &M) {
      CHECK_SAME(dim,M.dim); (*D) += (*M.D); return *this;
    }
    DiagonalMatrix& operator-=(const DiagonalMatrix &M) {
      CHECK_SAME(dim,M.dim); (*D) -= (*M.D); return *this;
    }
    virtual DiagonalMatrix& operator*=(Z d) {
      (*D) *= d; return *this;
    }
    virtual DiagonalMatrix& operator/=(Z d) {
      (*D) /= d; return *this;
    }

    virtual DiagonalMatrix* operator+(const DiagonalMatrix &M) const {
      DiagonalMatrix *S=OPERATOR_NEW DiagonalMatrix(dim);
      S->copy(*this); *S+=M; return S;
    }
    virtual TridiagonalMatrix<F,Z>* operator+(
      const TridiagonalMatrix<F,Z> &T) const;
    virtual UpperTriangularMatrix<F,Z>* operator+(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual LowerTriangularMatrix<F,Z>* operator+(
      const LowerTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &U) const;

    virtual DiagonalMatrix* operator-(const DiagonalMatrix &M) const {
      DiagonalMatrix *S=OPERATOR_NEW DiagonalMatrix(dim);
      S->copy(*this); *S-=M; return S;
    }
    virtual TridiagonalMatrix<F,Z>* operator-(
      const TridiagonalMatrix<F,Z> &M) const;
    virtual UpperTriangularMatrix<F,Z>* operator-(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual LowerTriangularMatrix<F,Z>* operator-(
      const LowerTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator-(const Matrix<F,Z> &U) const;

    virtual DiagonalMatrix* operator*(Z d) const {
      DiagonalMatrix *P=OPERATOR_NEW DiagonalMatrix(dim);
      P->copy(*this); *P*=d; return P;
    }
    virtual DiagonalMatrix* operator/(Z d) const {
      DiagonalMatrix *P=OPERATOR_NEW DiagonalMatrix(dim);
      P->copy(*this); *P/=d; return P;
    }

    virtual DiagonalMatrix* operator*(const DiagonalMatrix &M) const;
    TridiagonalMatrix<F,Z>* operator*(
      const SymmetricTridiagonalMatrix<F,Z> &T) const;
    virtual TridiagonalMatrix<F,Z>* operator*(
      const TridiagonalMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator*(const SymmetricMatrix<F,Z> &S)
      const;
    virtual UpperTrapezoidalMatrix<F,Z>* operator*(
      const UpperTrapezoidalMatrix<F,Z> &L) const;
    virtual LowerTrapezoidalMatrix<F,Z>* operator*(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator*(const SquareMatrix<F,Z> &S)
      const;
    virtual Matrix<F,Z>* operator*(const Matrix<F,Z> &M) const;
    virtual Vector<F,Z>* operator*(const Vector<F,Z> &v) const;

    void dmv(Z alpha,const Vector<F,Z> &x,Z beta,Vector<F,Z> &b) const;
     // b=D*x*alpha+b*beta
    void dmm(Z alpha,const Matrix<F,Z> &X,Z beta,Matrix<F,Z> &B,
      char side='L') const; // B=D*X*alpha+B*beta or B=X*D*alpha+B*beta

    virtual void solve(const Vector<F,Z> &b,Vector<F,Z> &x) const;
    virtual void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,char side='L')
      const;
    virtual void printOn(ostream& =cerr) const ;
};

template<typename F> SymmetricTridiagonalMatrix<F,F>* operator+(
  const DiagonalMatrix<F,F> &A,const SymmetricTridiagonalMatrix<F,F> &T);
template<typename F> SymmetricTridiagonalMatrix<F,F>* operator+(
const SymmetricTridiagonalMatrix<F,F> &T,const DiagonalMatrix<F,F> &A) {
  return A+T;
}
template<typename F> SymmetricTridiagonalMatrix<F,complex<F> >*
  operator+(const DiagonalMatrix<F,F> &A,
  const SymmetricTridiagonalMatrix<F,complex<F> > &T);
template<typename F> SymmetricTridiagonalMatrix<F,complex<F> >*
operator+(const SymmetricTridiagonalMatrix<F,complex<F> > &T,
const DiagonalMatrix<F,F> &A) {
  return A+T;
}
template<typename F> TridiagonalMatrix<F,complex<F> >* operator+(
  const DiagonalMatrix<F,complex<F> > &A,
  const SymmetricTridiagonalMatrix<F,complex<F> > &T);
template<typename F> TridiagonalMatrix<F,complex<F> >* operator+(
const SymmetricTridiagonalMatrix<F,complex<F> > &T,
const DiagonalMatrix<F,complex<F> > &A) {
  return A+T;
}
template<typename F,typename Z>
TridiagonalMatrix<F,Z>* operator+(const TridiagonalMatrix<F,Z> &T,
const DiagonalMatrix<F,Z> &A) {
  return A+T;
}
template<typename F> SymmetricMatrix<F,F>* operator+(
  const DiagonalMatrix<F,F> &A,const SymmetricMatrix<F,F> &T);
template<typename F> SymmetricMatrix<F,F>* operator+(
const SymmetricMatrix<F,F> &T,const DiagonalMatrix<F,F> &A) {
  return A+T;
}
template<typename F> SymmetricMatrix<F,complex<F> >* operator+(
  const DiagonalMatrix<F,F> &A,const SymmetricMatrix<F,complex<F> > &T);
template<typename F> SymmetricMatrix<F,complex<F> >* operator+(
const SymmetricMatrix<F,complex<F> > &T,
const DiagonalMatrix<F,F> &A) {
  return A+T;
}
template<typename F> SquareMatrix<F,complex<F> >* operator+(
  const DiagonalMatrix<F,complex<F> > &A,
  const SymmetricMatrix<F,complex<F> > &T);
template<typename F> SquareMatrix<F,complex<F> >* operator+(
const SymmetricMatrix<F,complex<F> > &T,
const DiagonalMatrix<F,complex<F> > &A) {
  return A+T;
}
template<typename F,typename Z> UpperTriangularMatrix<F,Z>* operator+(
const UpperTrapezoidalMatrix<F,Z> &U,const DiagonalMatrix<F,Z> &A) {
  return A+U;
}
template<typename F,typename Z> LowerTriangularMatrix<F,Z>* operator+(
const LowerTrapezoidalMatrix<F,Z> &L,const DiagonalMatrix<F,Z> &A) {
  return A+L;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const Matrix<F,Z> &M,const DiagonalMatrix<F,Z> &A) {
  return A+M;
}

template<typename F> SymmetricTridiagonalMatrix<F,F>* operator-(
  const DiagonalMatrix<F,F> &A,const SymmetricTridiagonalMatrix<F,F> &T);
template<typename F> SymmetricTridiagonalMatrix<F,F>* operator-(
  const SymmetricTridiagonalMatrix<F,F> &T,const DiagonalMatrix<F,F> &A);
template<typename F> SymmetricTridiagonalMatrix<F,complex<F> >*
  operator-(const DiagonalMatrix<F,F> &A,
  const SymmetricTridiagonalMatrix<F,complex<F> > &T);
template<typename F> SymmetricTridiagonalMatrix<F,complex<F> >*
  operator-(const SymmetricTridiagonalMatrix<F,complex<F> > &T,
  const DiagonalMatrix<F,F> &A);
template<typename F> TridiagonalMatrix<F,complex<F> >* operator-(
  const DiagonalMatrix<F,complex<F> > &A,
  const SymmetricTridiagonalMatrix<F,complex<F> > &T);
template<typename F> TridiagonalMatrix<F,complex<F> >* operator-(
  const SymmetricTridiagonalMatrix<F,complex<F> > &T,
  const DiagonalMatrix<F,complex<F> > &A);
template<typename F,typename Z>
  TridiagonalMatrix<F,Z>* operator-(const TridiagonalMatrix<F,Z> &T,
  const DiagonalMatrix<F,Z> &A);
template<typename F> SymmetricMatrix<F,F>* operator-(
  const DiagonalMatrix<F,F> &A,const SymmetricMatrix<F,F> &T);
template<typename F> SymmetricMatrix<F,F>* operator-(
  const SymmetricMatrix<F,F> &T,const DiagonalMatrix<F,F> &A);
template<typename F> SymmetricMatrix<F,complex<F> >* operator-(
  const DiagonalMatrix<F,F> &A,const SymmetricMatrix<F,complex<F> > &T);
template<typename F> SymmetricMatrix<F,complex<F> >* operator-(
  const SymmetricMatrix<F,complex<F> > &T,
  const DiagonalMatrix<F,F> &A);
template<typename F> SquareMatrix<F,complex<F> >* operator-(
  const DiagonalMatrix<F,complex<F> > &A,
  const SymmetricMatrix<F,complex<F> > &T);
template<typename F> SquareMatrix<F,complex<F> >* operator-(
  const SymmetricMatrix<F,complex<F> > &T,
  const DiagonalMatrix<F,complex<F> > &A);
template<typename F,typename Z> UpperTriangularMatrix<F,Z>* operator-(
  const UpperTrapezoidalMatrix<F,Z> &U,const DiagonalMatrix<F,Z> &A);
template<typename F,typename Z> LowerTriangularMatrix<F,Z>* operator-(
  const LowerTrapezoidalMatrix<F,Z> &U,const DiagonalMatrix<F,Z> &A);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const Matrix<F,Z> &U,const DiagonalMatrix<F,Z> &A);

template<typename F,typename Z> TridiagonalMatrix<F,Z>* operator*(
  const SymmetricTridiagonalMatrix<F,Z> &T,const DiagonalMatrix<F,Z> &A);
template<typename F,typename Z> TridiagonalMatrix<F,Z>* operator*(
  const TridiagonalMatrix<F,Z> &T,const DiagonalMatrix<F,Z> &A);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const DiagonalMatrix<F,Z> &A,const SymmetricMatrix<F,Z> &M);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SymmetricMatrix<F,Z> &M,const DiagonalMatrix<F,Z> &A);
template<typename F,typename Z> UpperTrapezoidalMatrix<F,Z>* operator*(
  const UpperTrapezoidalMatrix<F,Z> &U,const DiagonalMatrix<F,Z> &A);
template<typename F,typename Z> LowerTrapezoidalMatrix<F,Z>* operator*(
  const LowerTrapezoidalMatrix<F,Z> &U,const DiagonalMatrix<F,Z> &A);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SquareMatrix<F,Z> &U,const DiagonalMatrix<F,Z> &A);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const Matrix<F,Z> &U,const DiagonalMatrix<F,Z> &A);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> class UpperHessenbergMatrix :
public SquareMatrix<F,Z> {
  private:
//  prevent copy constructor
    explicit UpperHessenbergMatrix(const UpperHessenbergMatrix &H) : 
      SquareMatrix<F,Z>(H.size(0)) { copy(H); }
//  use copy instead of operator=
    UpperHessenbergMatrix& operator=(
    const UpperHessenbergMatrix& H) {
      SquareMatrix<F,Z>::resize(H.size(0)); copy(H); return *this;
    }
//  prevent SquareMatrix functions
    virtual void copyFrom(char uplo,int m,int n,const Matrix<F,Z> &A) {
      return SquareMatrix<F,Z>::copyFrom(uplo,m,n,A);
    }

    int factor(int *ipiv);
  protected:
    virtual void fillWith(Z d);
  public:
    static const Z outofbounds_;
    static const Z undefined_;
    static Z safety_;

    UpperHessenbergMatrix() : SquareMatrix<F,Z>() { }
    explicit UpperHessenbergMatrix(int n): SquareMatrix<F,Z>(n,n) {;}
    explicit UpperHessenbergMatrix(int n,Z d): SquareMatrix<F,Z>(n,n) {
      fillWith(d);
    }
//  this constructor allows us to find eigenvalues for Tridiagonal's
    explicit UpperHessenbergMatrix(const TridiagonalMatrix<F,Z> &T);
    virtual ~UpperHessenbergMatrix() {;}
    virtual SquareMatrix<F,Z>* makeMatrix() const;

    virtual UpperHessenbergMatrix& operator=(Z scalar) { 
      fillWith(scalar); return *this; 
    }
    virtual Z* addr() { return SquareMatrix<F,Z>::addr(); }
    virtual const Z* addr() const { return SquareMatrix<F,Z>::addr(); }
    virtual Z* addr(int i,int j) {
      CHECK_TEST(i>=0 && i<=j+1);
      return SquareMatrix<F,Z>::addr(i,j);
    }
    virtual const Z* addr(int i,int j) const {
      CHECK_TEST(i>=0 && i<=j+1);
      return SquareMatrix<F,Z>::addr(i,j);
    }
    virtual Z& operator()(int i,int j) {
      if (i<=j+1) return SquareMatrix<F,Z>::operator()(i,j);
      return safety_=outofbounds_;
    }
    virtual Z operator()(int i,int j) const {
      if (i<=j+1) return SquareMatrix<F,Z>::operator()(i,j);
      return outofbounds_;
    }
    virtual void copy(const Matrix<F,Z> &S);
    UpperHessenbergMatrix& operator+=(const UpperHessenbergMatrix &H);
    UpperHessenbergMatrix& operator-=(const UpperHessenbergMatrix &H);
    virtual UpperHessenbergMatrix& operator*=(Z scalar);
    virtual UpperHessenbergMatrix& operator/=(Z scalar);

    virtual UpperHessenbergMatrix* operator+(
    const UpperHessenbergMatrix &H) const {
      UpperHessenbergMatrix *S=OPERATOR_NEW UpperHessenbergMatrix(*this);
      *S+=H; return S;
    }
    virtual UpperHessenbergMatrix* operator+(
      const DiagonalMatrix<F,Z> &T) const;
    virtual UpperHessenbergMatrix* operator+(
      const SymmetricTridiagonalMatrix<F,Z> &T) const;
    virtual UpperHessenbergMatrix* operator+(
    const SymmetricPositiveTridiagonalMatrix<F,Z> &T) const {
      return
        operator+(static_cast<const SymmetricTridiagonalMatrix<F,Z>&>(T));
    }
    virtual UpperHessenbergMatrix* operator+(
      const TridiagonalMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator+(const SymmetricMatrix<F,Z> &S)
      const;
    virtual SquareMatrix<F,Z>* operator+(
    const SymmetricPositiveMatrix<F,Z> &S) const {
      return operator+(static_cast<const SymmetricMatrix<F,Z>&>(S));
    }
    virtual UpperHessenbergMatrix* operator+(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual UpperHessenbergMatrix* operator+(
    const UnitUpperTriangularMatrix<F,Z> &U) const {
      return operator+(static_cast<const UpperTrapezoidalMatrix<F,Z>&>(U));
    }
    virtual SquareMatrix<F,Z>* operator+(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator+(
    const UnitLowerTriangularMatrix<F,Z> &L) const {
      return operator+(static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L));
    }
    virtual SquareMatrix<F,Z>* operator+(
    const UnitLowerTrapezoidalMatrix<F,Z> &L) const {
      return operator+(static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L));
    }
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &M) const;
    virtual SquareMatrix<F,Z>* operator+(const SquareMatrix<F,Z> &S)
    const {
      return operator+(static_cast<const Matrix<F,Z>&>(S));
    }
    virtual SquareMatrix<F,Z>* operator+(const OrthogonalMatrix<F,Z> &Q)
    const {
      return operator+(static_cast<const Matrix<F,Z>&>(Q));
    }
    virtual UpperHessenbergMatrix* operator-(
    const UpperHessenbergMatrix &H) const {
      UpperHessenbergMatrix *S=OPERATOR_NEW UpperHessenbergMatrix(*this);
      *S-=H; return S;
    }
    virtual UpperHessenbergMatrix* operator-(const DiagonalMatrix<F,Z> &T)
      const;
    virtual UpperHessenbergMatrix* operator-(
      const SymmetricTridiagonalMatrix<F,Z> &T) const;
    virtual UpperHessenbergMatrix* operator-(
      const TridiagonalMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator-(const SymmetricMatrix<F,Z> &S)
      const;
    virtual UpperHessenbergMatrix* operator-(
      const UpperTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator-(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator-(const Matrix<F,Z> &M) const;
    virtual SquareMatrix<F,Z>* operator-(const SquareMatrix<F,Z> &S)
    const {
      return operator-(static_cast<const Matrix<F,Z>&>(S));
    }
    virtual SquareMatrix<F,Z>* operator-(const OrthogonalMatrix<F,Z> &Q)
    const {
      return operator-(static_cast<const Matrix<F,Z>&>(Q));
    }

    virtual UpperHessenbergMatrix* operator*(Z scalar) const {
      UpperHessenbergMatrix *S=OPERATOR_NEW UpperHessenbergMatrix(*this);
      *S*=scalar; return S;
    }
    virtual UpperHessenbergMatrix* operator/(Z scalar) const {
      UpperHessenbergMatrix *S=OPERATOR_NEW UpperHessenbergMatrix(*this);
      *S/=scalar; return S;
    }

    virtual SquareMatrix<F,Z>* operator*(const UpperHessenbergMatrix<F,Z>&)
      const;
    virtual UpperHessenbergMatrix* operator*(const DiagonalMatrix<F,Z>&)
      const;
    virtual SquareMatrix<F,Z>* operator*(
      const SymmetricTridiagonalMatrix<F,Z>&) const;
    virtual SquareMatrix<F,Z>* operator*(const TridiagonalMatrix<F,Z>&)
      const;
    virtual SquareMatrix<F,Z>* operator*(const SymmetricMatrix<F,Z>&)
      const;
    virtual Matrix<F,Z>* operator*(const UpperTrapezoidalMatrix<F,Z>&)
      const;
    virtual Matrix<F,Z>* operator*(const LowerTrapezoidalMatrix<F,Z>&)
      const;
    virtual SquareMatrix<F,Z>* operator*(const SquareMatrix<F,Z>&) const;
    virtual Matrix<F,Z>* operator*(const Matrix<F,Z>&) const;
    virtual Vector<F,Z>* operator*(const Vector<F,Z>&) const;

    void uhmv(Z alpha,const Vector<F,Z> &x,Z beta,Vector<F,Z> &b,
      char trans='N') const; // b=H*x*alpha+b*beta or b=H^T*x*alpha+b*beta
    void uhmm(Z alpha,const Matrix<F,Z> &X,Z beta,Matrix<F,Z> &B,
      char side='L',char trans='N') const;
      // B=H*X*alpha+B*beta or B=X*H*alpha+B*beta or
      // B=H^T*X*alpha+B*beta or B=X*H^T*alpha+B*beta

//  SquareMatrix<F,Z>* transpose() const ;
//  LaPack:
    virtual void copyFrom(int m,const SquareMatrix<F,Z> &M);
    F normFrobenius() const;
    F normInfinity() const;
    F normMaxEntry() const;
    F normOne() const;
    Vector<F,complex<F> >* eigenvalues(SquareMatrix<F,complex<F> > *&V,
      SquareMatrix<F,complex<F> > *&U) const; // _geev
    virtual void solve(const Vector<F,Z>&,Vector<F,Z>&,char trans='N')
      const;
    virtual void solve(const Matrix<F,Z>&,Matrix<F,Z>&,char side='L',
      char trans='N') const;
    void printOn(ostream& = cerr) const;
};

template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator+(
const DiagonalMatrix<F,Z> &T,const UpperHessenbergMatrix<F,Z> &H) {
  return H+T;
}
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator+(
const SymmetricTridiagonalMatrix<F,Z> &T,
const UpperHessenbergMatrix<F,Z> &H) {
  return H+T;
}
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator+(
const SymmetricPositiveTridiagonalMatrix<F,Z> &T,
const UpperHessenbergMatrix<F,Z> &H) {
  return H+T;
}
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator+(
const TridiagonalMatrix<F,Z> &T,
const UpperHessenbergMatrix<F,Z> &H) {
  return H+T;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
SymmetricMatrix<F,Z> &S,const UpperHessenbergMatrix<F,Z> &H) {
  return H+S;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
SymmetricPositiveMatrix<F,Z> &S,const UpperHessenbergMatrix<F,Z> &H) {
  return H+S;
}
template<typename F,typename Z>  UpperHessenbergMatrix<F,Z>* operator+(
const UnitUpperTriangularMatrix<F,Z> &U,
const UpperHessenbergMatrix<F,Z> &H) {
  return H+U;
}
template<typename F,typename Z>  UpperHessenbergMatrix<F,Z>* operator+(
const UpperTrapezoidalMatrix<F,Z> &U,
const UpperHessenbergMatrix<F,Z> &H) {
  return H+U;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator+(
const UnitLowerTriangularMatrix<F,Z> &L,
const UpperHessenbergMatrix<F,Z> &H) {
  return H+L;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator+(
const UnitLowerTrapezoidalMatrix<F,Z> &L,
const UpperHessenbergMatrix<F,Z> &H) {
  return H+L;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator+(
const LowerTrapezoidalMatrix<F,Z> &L,
const UpperHessenbergMatrix<F,Z> &H) {
  return H+L;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator+(
const OrthogonalMatrix<F,Z> &Q,const UpperHessenbergMatrix<F,Z> &H) {
  return H+Q;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator+(
const SquareMatrix<F,Z> &S,const UpperHessenbergMatrix<F,Z> &H) {
  return H+S;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator+(
const Matrix<F,Z> &M,const UpperHessenbergMatrix<F,Z> &H) {
  return H+M;
}

template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator-(
  const DiagonalMatrix<F,Z> &T,const UpperHessenbergMatrix<F,Z> &H);
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator-(
  const SymmetricTridiagonalMatrix<F,Z> &T,
  const UpperHessenbergMatrix<F,Z> &H);
template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator-(
  const TridiagonalMatrix<F,Z> &T,const UpperHessenbergMatrix<F,Z> &H);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const SymmetricMatrix<F,Z> &S,const UpperHessenbergMatrix<F,Z> &H);
template<typename F,typename Z>  UpperHessenbergMatrix<F,Z>* operator-(
  const UpperTrapezoidalMatrix<F,Z> &U,
  const UpperHessenbergMatrix<F,Z> &H);
template<typename F,typename Z>  UpperHessenbergMatrix<F,Z>* operator-(
const UnitUpperTriangularMatrix<F,Z> &U,
const UpperHessenbergMatrix<F,Z> &H) {
  return static_cast<const UpperTrapezoidalMatrix<F,Z>&>(U)-H;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator-(
  const LowerTrapezoidalMatrix<F,Z> &L,
  const UpperHessenbergMatrix<F,Z> &H);
template<typename F,typename Z>  SquareMatrix<F,Z>* operator-(
const UnitLowerTriangularMatrix<F,Z> &L,
const UpperHessenbergMatrix<F,Z> &H) {
  return static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L)-H;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator-(
const UnitLowerTrapezoidalMatrix<F,Z> &L,
const UpperHessenbergMatrix<F,Z> &H) {
  return static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L)-H;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator-(
  const Matrix<F,Z> &M,const UpperHessenbergMatrix<F,Z> &H);
template<typename F,typename Z>  SquareMatrix<F,Z>* operator-(
const SquareMatrix<F,Z> &S,const UpperHessenbergMatrix<F,Z> &H) {
  return static_cast<const Matrix<F,Z>&>(S)-H;
}
template<typename F,typename Z>  SquareMatrix<F,Z>* operator-(
const OrthogonalMatrix<F,Z> &Q,const UpperHessenbergMatrix<F,Z> &H) {
  return static_cast<const Matrix<F,Z>&>(Q)-H;
}

template<typename F,typename Z> UpperHessenbergMatrix<F,Z>* operator*(
  const DiagonalMatrix<F,Z>&,const UpperHessenbergMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SymmetricTridiagonalMatrix<F,Z>&,
  const UpperHessenbergMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const TridiagonalMatrix<F,Z>&,const UpperHessenbergMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SymmetricMatrix<F,Z>&,const UpperHessenbergMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SquareMatrix<F,Z>&,const UpperHessenbergMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const UpperTrapezoidalMatrix<F,Z>&,const UpperHessenbergMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
const UnitUpperTriangularMatrix<F,Z> &U,
const UpperHessenbergMatrix<F,Z> &H) {
  return static_cast<const UpperTrapezoidalMatrix<F,Z>&>(U)*H;
}
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const LowerTrapezoidalMatrix<F,Z>&,const UpperHessenbergMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
const UnitLowerTriangularMatrix<F,Z> &L,
const UpperHessenbergMatrix<F,Z> &H) {
  return static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L)*H;
}
template<typename F,typename Z> Matrix<F,Z>* operator*(
const UnitLowerTrapezoidalMatrix<F,Z> &L,
const UpperHessenbergMatrix<F,Z> &H) {
  return static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L)*H;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SquareMatrix<F,Z>&,const UpperHessenbergMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const Matrix<F,Z>&,const UpperHessenbergMatrix<F,Z>&);

template<typename F,typename Z> void testUpperHessenbergMatrix(F,Z);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> class BandMatrix {
  private:
    Matrix<F,Z> *AB;
    int dim,nsub,nsup,nt;
//  note that we should have nt < dim, 
//    else SquareMatrix uses less storage than BandMatrix

//  prevent public use of copy constructor
    explicit BandMatrix(const BandMatrix &B) : dim(B.dim), nsub(B.nsub),
    nsup(B.nsup),nt(B.nt) {
      AB=OPERATOR_NEW Matrix<F,Z>(B.AB->size(0),B.AB->size(1));
      AB->copy(*B.AB);
    }
// Use copy instead
    BandMatrix& operator=(const BandMatrix &B) { 
      resize(B.dim,B.nsub,B.nsup); AB->copy(*B.AB); return *this;
    }
  protected:
    virtual void fillWith(Z scalar);
  public:
    static const Z outofbounds_;
    static Z safety_;

    BandMatrix() : AB(0),dim(0),nsub(0),nsup(0),nt(1) {;}
    explicit BandMatrix(int n,int nb,int np): dim(n) {
//    nsub=min(nb,n-1);
//    nsup=min(np,n-1);
      nsub=nb; // LaPack routines do not check for more bands than rows
      nsup=np;
      nt=nsub+nsup+1;
      AB=OPERATOR_NEW Matrix<F,Z>(nt,n);
    }
    explicit BandMatrix(int n,int nb,int np,Z d): dim(n) {
//    nsub=min(nb,n-1);
//    nsup=min(np,n-1);
      nsub=nb; // LaPack routines do not check for more bands than rows
      nsup=np;
      nt=nsub+nsup+1;
      AB=OPERATOR_NEW Matrix<F,Z>(nt,n);
      fillWith(d);
    }
    virtual ~BandMatrix() { delete AB; AB=0; } 
    virtual SquareMatrix<F,Z>* makeMatrix() const;
    virtual BandMatrix& operator=(Z scalar) {
      fillWith(scalar); return *this;
    }

    int size(int ) const { return dim; }
    int subDiags() const { return nsub; }
    int supDiags() const { return nsup; }
    int bands() const { return nt; }
    virtual Z* addr() { return AB->addr(); }
    virtual const Z* addr() const { return AB->addr(); }
    virtual Z* addr(int i, int j) {
      CHECK_BOUNDS(i,0,dim)
      CHECK_BOUNDS(j,0,dim)
      CHECK_TEST(i-j<=nsub && j-i<=nsup);
      return AB->addr(nsup+i-j,j);
    }
    virtual const Z* addr(int i, int j) const {
      CHECK_BOUNDS(i,0,dim)
      CHECK_BOUNDS(j,0,dim)
      CHECK_TEST(i-j<=nsub && j-i<=nsup);
      return AB->addr(nsup+i-j,j);
    }
    virtual Z& operator()(int i,int j) { 
      if (i-j<=nsub && j-i<=nsup) return *AB->addr(nsup+i-j,j);
      else return safety_=outofbounds_;
    }
    virtual Z operator()(int i,int j) const { 
      if (i-j<=nsub && j-i<=nsup) return *AB->addr(nsup+i-j,j);
      else return outofbounds_;
    }
    virtual void resize(int n,int nb,int np) { 
      if (dim==n && nsub==nb && nsup==np) return;
      dim=n; nsub=nb; nsup=np; nt=nsub+nsup+1;
      if (AB) delete AB;
      AB=OPERATOR_NEW Matrix<F,Z>(nt,n);
    }
    virtual void resize(const BandMatrix &B) {
      resize(B.size(0),B.nsub,B.nsup);
    }
    virtual void copy(const BandMatrix &B) {
      CHECK_SAME(dim,B.dim);
      AB->copy(*B.AB);
    }

    BandMatrix& operator+=(const BandMatrix &B);
    BandMatrix& operator-=(const BandMatrix &B);
    virtual BandMatrix& operator*=(Z d);
    virtual BandMatrix& operator/=(Z d);

    virtual BandMatrix* operator+(const BandMatrix &B) const;
    virtual SquareMatrix<F,Z>* operator+(
      const UpperHessenbergMatrix<F,Z> &H) const;
    virtual BandMatrix* operator+(const DiagonalMatrix<F,Z> &D) const;
    virtual BandMatrix* operator+(
      const SymmetricTridiagonalMatrix<F,Z> &T) const;
    virtual BandMatrix* operator+(
    const SymmetricPositiveTridiagonalMatrix<F,Z> &T) const {
      return
        operator+(static_cast<const SymmetricTridiagonalMatrix<F,Z>&>(T));
    }
    virtual BandMatrix* operator+(const TridiagonalMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator+(const SymmetricMatrix<F,Z> &H)
      const;
    virtual SquareMatrix<F,Z>* operator+(
    const SymmetricPositiveMatrix<F,Z> &H) const {
      return operator+(static_cast<const SymmetricMatrix<F,Z>&>(H));
    }
    virtual SquareMatrix<F,Z>* operator+(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator+(
    const UnitUpperTriangularMatrix<F,Z> &U) const {
      return operator+(static_cast<const UpperTrapezoidalMatrix<F,Z>&>(U));
    }
    virtual SquareMatrix<F,Z>* operator+(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator+(
    const UnitLowerTriangularMatrix<F,Z> &L) const {
      return operator+(static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L));
    }
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &M) const;
    virtual SquareMatrix<F,Z>* operator+(const SquareMatrix<F,Z> &S)
    const {
      return operator+(static_cast<const Matrix<F,Z>&>(S));
    }
    virtual SquareMatrix<F,Z>* operator+(const OrthogonalMatrix<F,Z> &Q)
    const {
      return operator+(static_cast<const Matrix<F,Z>&>(Q));
    }

    virtual BandMatrix* operator-(const BandMatrix &B) const;
    virtual SquareMatrix<F,Z>* operator-(
      const UpperHessenbergMatrix<F,Z> &H) const;
    virtual BandMatrix* operator-(const DiagonalMatrix<F,Z> &D) const;
    virtual BandMatrix* operator-(
      const SymmetricTridiagonalMatrix<F,Z> &T) const;
    virtual BandMatrix* operator-(
    const SymmetricPositiveTridiagonalMatrix<F,Z> &T) const {
      return
        operator-(static_cast<const SymmetricTridiagonalMatrix<F,Z>&>(T));
    }
    virtual BandMatrix* operator-(const TridiagonalMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator-(const SymmetricMatrix<F,Z> &H)
      const;
    virtual SquareMatrix<F,Z>* operator-(
    const SymmetricPositiveMatrix<F,Z> &H) const {
      return operator-(static_cast<const SymmetricMatrix<F,Z>&>(H));
    }
    virtual SquareMatrix<F,Z>* operator-(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator-(
    const UnitUpperTriangularMatrix<F,Z> &U) const {
      return operator-(static_cast<const UpperTrapezoidalMatrix<F,Z>&>(U));
    }
    virtual SquareMatrix<F,Z>* operator-(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator-(
    const UnitLowerTriangularMatrix<F,Z> &L) const {
      return operator-(static_cast<const LowerTrapezoidalMatrix<F,Z>&>(L));
    }
    virtual SquareMatrix<F,Z>* operator-(const Matrix<F,Z> &M) const;
    virtual SquareMatrix<F,Z>* operator-(const SquareMatrix<F,Z> &S)
    const {
      return operator-(static_cast<const Matrix<F,Z>&>(S));
    }
    virtual SquareMatrix<F,Z>* operator-(const OrthogonalMatrix<F,Z> &Q)
    const {
      return operator-(static_cast<const Matrix<F,Z>&>(Q));
    }

    virtual BandMatrix* operator*(Z d) const {
      BandMatrix *prod=OPERATOR_NEW BandMatrix(*this);
      (*prod) *= d; return prod;
    }
    virtual BandMatrix* operator/(Z d) const {
      BandMatrix *quot=OPERATOR_NEW BandMatrix(*this);
      (*quot) /= d; return quot;
    }

    virtual BandMatrix* operator*(const BandMatrix &B) const;
    virtual SquareMatrix<F,Z>* operator*(
      const UpperHessenbergMatrix<F,Z>&) const;
    virtual BandMatrix* operator*(const DiagonalMatrix<F,Z> &D) const;
    virtual BandMatrix* operator*(
      const SymmetricTridiagonalMatrix<F,Z> &T) const;
    virtual BandMatrix* operator*(const TridiagonalMatrix<F,Z> &T) const;
    virtual SquareMatrix<F,Z>* operator*(const SymmetricMatrix<F,Z> &S)
      const;
    virtual Matrix<F,Z>* operator*(const UpperTrapezoidalMatrix<F,Z> &U)
      const;
    virtual Matrix<F,Z>* operator*(const LowerTrapezoidalMatrix<F,Z> &U)
      const;
    virtual SquareMatrix<F,Z>* operator*(const SquareMatrix<F,Z> &S) const;
    virtual Matrix<F,Z>* operator*(const Matrix<F,Z> &S) const;
    virtual Vector<F,Z>* operator*(const Vector<F,Z> &S) const;

//  Blas:
    void gbmv(Z alpha,const Vector<F,Z> &x,Z beta,Vector<F,Z> &b,
      char trans='N') const; // b=A*x*alpha+b*beta or b=A^T*x*alpha+b*beta
    void gbmm(Z alpha,const Matrix<F,Z> &X,Z beta,Matrix<F,Z> &B,
      char side='L',char trans='N') const;
      // B=A*X*alpha+B*beta or B=A^T*X*alpha+B*beta or
      // B=X*A*alpha+B*beta or B=X*A^T*alpha+B*beta
//  LaPack:
    F equilibrate(Vector<F,F> &r,Vector<F,F> &c,F &rowcnd,F &colcnd)
      const;
//  void gbamv(F alpha,const Vector<F,Z> &x,F beta,Vector<F,F> &y)
//    const; // dla_gbamv
    F normFrobenius() const;
    F normInfinity() const;
    F normMaxEntry() const;
    F normOne() const;
    F reciprocalConditionNumber(char norm) const;

//  BandMatrix* transpose() const ;
    virtual void solve(const Vector<F,Z> &b,Vector<F,Z> &x,char trans='N')
      const;
    virtual void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,char side='L',
      char trans='N') const;
    virtual void printOn(ostream& =cerr) const ;
};

template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const UpperHessenbergMatrix<F,Z> &H,const BandMatrix<F,Z> &B) {
  return B+H;
}
template<typename F,typename Z> BandMatrix<F,Z>* operator+(
const DiagonalMatrix<F,Z> &D,const BandMatrix<F,Z> &B) {
  return B+D;
}
template<typename F,typename Z> BandMatrix<F,Z>* operator+(
const SymmetricTridiagonalMatrix<F,Z> &T,const BandMatrix<F,Z> &B) {
  return B+T;
}
template<typename F,typename Z> BandMatrix<F,Z>* operator+(
const TridiagonalMatrix<F,Z> &T,const BandMatrix<F,Z> &B) {
  return B+T;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const SymmetricMatrix<F,Z> &H,const BandMatrix<F,Z> &B) {
  return B+H;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const UpperTrapezoidalMatrix<F,Z> &U,const BandMatrix<F,Z> &B) {
  return B+U;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const LowerTrapezoidalMatrix<F,Z> &L,const BandMatrix<F,Z> &B) {
  return B+L;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const Matrix<F,Z> &M,const BandMatrix<F,Z> &B) {
  return B+M;
}

template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const UpperHessenbergMatrix<F,Z> &H,const BandMatrix<F,Z> &B);
template<typename F,typename Z> BandMatrix<F,Z>* operator-(
  const DiagonalMatrix<F,Z> &D,const BandMatrix<F,Z> &B);
template<typename F,typename Z> BandMatrix<F,Z>* operator-(
  const SymmetricTridiagonalMatrix<F,Z> &,const BandMatrix<F,Z> &B);
template<typename F,typename Z> BandMatrix<F,Z>* operator-(
  const TridiagonalMatrix<F,Z> &T,const BandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const SymmetricMatrix<F,Z> &H,const BandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const UpperTrapezoidalMatrix<F,Z> &U,const BandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const LowerTrapezoidalMatrix<F,Z> &L,const BandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const Matrix<F,Z> &M,const BandMatrix<F,Z> &B);

template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const UpperHessenbergMatrix<F,Z>&,const BandMatrix<F,Z>&);
template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const DiagonalMatrix<F,Z>&,const BandMatrix<F,Z>&);
template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const SymmetricTridiagonalMatrix<F,Z>&,const BandMatrix<F,Z>&);
template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const TridiagonalMatrix<F,Z>&,const BandMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SymmetricMatrix<F,Z>&,const BandMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const UpperTrapezoidalMatrix<F,Z>&,const BandMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const LowerTrapezoidalMatrix<F,Z>&,const BandMatrix<F,Z>&);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SquareMatrix<F,Z>&,const BandMatrix<F,Z>&);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const Matrix<F,Z>&,const BandMatrix<F,Z>&);

template<typename F,typename Z> void testBandMatrix(F,Z);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//need class SymmetricBandMatrix and
//  constructor from SymmetricTridiagonalMatrix
template<typename F,typename Z> class SymmetricBandMatrix {
  private:
    Matrix<F,Z> *AB;
    int dim,nsub,nt;
//  note that we should have nsub < (dim-1)/2, else
//    packed SymmetricMatrix uses less storage than SymmetricBandMatrix

//  prevent public use of copy constructor
    explicit SymmetricBandMatrix(const SymmetricBandMatrix &B) :
    dim(B.dim), nsub(B.nsub), nt(B.nt) {
      AB=OPERATOR_NEW Matrix<F,Z>(B.AB->size(0),B.AB->size(1));
      AB->copy(*B.AB);
    }
// Use copy instead
    SymmetricBandMatrix& operator=(const SymmetricBandMatrix &B) {
      resize(B.dim,B.nsub); AB->copy(*B.AB); return *this;
    }
  protected:
    virtual void fillWith(Z scalar);
  public:
    static const Z outofbounds_;
    static Z safety_;

    SymmetricBandMatrix() : AB(0),dim(0),nsub(0),nt(1) {;}
    explicit SymmetricBandMatrix(int n,int nb) : dim(n) {
//    nsub=min(nb,n-1);
      nsub=nb;
      nt=nsub+1;
      AB=OPERATOR_NEW Matrix<F,Z>(nt,n);
    }
    explicit SymmetricBandMatrix(int n,int nb,Z d): dim(n) {
//    nsub=min(nb,n-1);
      nsub=nb;
      nt=nsub+1;
      AB=OPERATOR_NEW Matrix<F,Z>(nt,n);
      fillWith(d);
    }
    explicit SymmetricBandMatrix(const SymmetricTridiagonalMatrix<F,Z> &T);
    virtual ~SymmetricBandMatrix() { delete AB; AB=0; }
    virtual SquareMatrix<F,Z>* makeMatrix() const;
    virtual SymmetricBandMatrix& operator=(Z scalar) {
      fillWith(scalar); return *this;
    }

    SymmetricMatrix<F,Z>* makeSymmetricMatrix() const;
//  to solve a linear system or compute a condition number
    BandMatrix<F,Z>* makeBandMatrix() const;

    int size(int) const { return dim; }
    int subDiags() const { return nsub; }
    int bands() const { return nt; }
    virtual Z* addr() { return AB->addr(); }
    virtual const Z* addr() const { return AB->addr(); }
    virtual Z* addr(int i, int j) {
      CHECK_BOUNDS(i,0,dim)
      CHECK_BOUNDS(j,0,dim)
      CHECK_BOUNDS(i-j,0,nt)
      return AB->addr(i-j,j);
    }
    virtual const Z* addr(int i, int j) const {
      CHECK_BOUNDS(i,0,dim)
      CHECK_BOUNDS(j,0,dim)
      CHECK_BOUNDS(i-j,0,nt)
      return AB->addr(i-j,j);
    }
    virtual Z& operator()(int i,int j) {
      if (i-j>=0 && i-j<=nsub) return *AB->addr(i-j,j);
      else return safety_=outofbounds_;
    }
//  requires conjugation for Z complex
    virtual Z operator()(int i,int j) const;
    virtual void resize(int n,int nb) {
      if (dim==n && nsub==nb) return;
      dim=n; nsub=nb; nt=nsub+1;
      if (AB) delete AB;
      AB=OPERATOR_NEW Matrix<F,Z>(nt,n);
    }
    virtual void resize(const SymmetricBandMatrix &B) {
      resize(B.size(0),B.nsub);
    }
    virtual void copy(const SymmetricBandMatrix &B) {
      CHECK_SAME(dim,B.dim);
      AB->copy(*B.AB);
    }

    SymmetricBandMatrix& operator+=(const SymmetricBandMatrix &B);
    SymmetricBandMatrix& operator-=(const SymmetricBandMatrix &B);
//  if Z is complex, scale by real part only:
    virtual SymmetricBandMatrix& operator*=(Z d);
    virtual SymmetricBandMatrix& operator/=(Z d);

    virtual SymmetricBandMatrix* operator+(const SymmetricBandMatrix &B)
      const;
    virtual BandMatrix<F,Z>* operator+(const BandMatrix<F,Z> &B) const;
    virtual SquareMatrix<F,Z>* operator+(
      const UpperHessenbergMatrix<F,Z> &H) const;
    virtual SymmetricBandMatrix* operator+(const DiagonalMatrix<F,Z> &D)
      const;
    virtual SymmetricBandMatrix* operator+(
      const SymmetricTridiagonalMatrix<F,Z>&) const;
    virtual BandMatrix<F,Z>* operator+(const TridiagonalMatrix<F,Z> &T)
      const;
    virtual SymmetricMatrix<F,Z>* operator+(const SymmetricMatrix<F,Z> &H)
      const;
    virtual SquareMatrix<F,Z>* operator+(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator+(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &M) const;

    virtual SymmetricBandMatrix* operator-(const SymmetricBandMatrix &B)
      const;
    virtual BandMatrix<F,Z>* operator-(const BandMatrix<F,Z> &B) const;
    virtual SquareMatrix<F,Z>* operator-(
      const UpperHessenbergMatrix<F,Z> &H) const;
    virtual SymmetricBandMatrix* operator-(const DiagonalMatrix<F,Z> &D)
      const;
    virtual SymmetricBandMatrix* operator-(
      const SymmetricTridiagonalMatrix<F,Z>&) const;
    virtual BandMatrix<F,Z>* operator-(const TridiagonalMatrix<F,Z> &T)
      const;
    virtual SymmetricMatrix<F,Z>* operator-(const SymmetricMatrix<F,Z> &H)
      const;
    virtual SquareMatrix<F,Z>* operator-(
      const UpperTrapezoidalMatrix<F,Z> &U) const;
    virtual SquareMatrix<F,Z>* operator-(
      const LowerTrapezoidalMatrix<F,Z> &L) const;
    virtual SquareMatrix<F,Z>* operator-(const Matrix<F,Z> &M) const;
    
    virtual BandMatrix<F,Z>* operator*(Z scalar) const {
      BandMatrix<F,Z> *P=makeBandMatrix(); (*P)*=scalar;
      return P;
    }
    virtual BandMatrix<F,Z>* operator/(Z scalar) const {
      BandMatrix<F,Z> *P=makeBandMatrix(); (*P)/=scalar;
      return P;
    }
    virtual BandMatrix<F,Z>* operator*(const SymmetricBandMatrix &B) const;
    virtual BandMatrix<F,Z>* operator*(const BandMatrix<F,Z> &B) const;
    virtual SquareMatrix<F,Z>* operator*(
      const UpperHessenbergMatrix<F,Z>&) const;
    virtual BandMatrix<F,Z>* operator*(const DiagonalMatrix<F,Z> &D) const;
    virtual BandMatrix<F,Z>* operator*(
      const SymmetricTridiagonalMatrix<F,Z> &T) const;
    virtual BandMatrix<F,Z>* operator*(const TridiagonalMatrix<F,Z> &T)
      const;
    virtual SquareMatrix<F,Z>* operator*(const SymmetricMatrix<F,Z> &S)
      const;
    virtual Matrix<F,Z>* operator*(const UpperTrapezoidalMatrix<F,Z> &U)
      const;
    virtual Matrix<F,Z>* operator*(const LowerTrapezoidalMatrix<F,Z> &U)
      const;
    virtual SquareMatrix<F,Z>* operator*(const SquareMatrix<F,Z> &S) const;
    virtual Matrix<F,Z>* operator*(const Matrix<F,Z> &S) const;
    virtual Vector<F,Z>* operator*(const Vector<F,Z> &S) const;

//  Blas:
    void sbmv(Z alpha,const Vector<F,Z> &x,Z beta,Vector<F,Z> &b) const;
      // b=A*x*alpha+b*beta
    void sbmm(Z alpha,const Matrix<F,Z> &X,Z beta,Matrix<F,Z> &B,
      char side='L') const; // B=A*X*alpha+B*beta or B = X*A*alpha+B*beta
//  LaPack:
    F normFrobenius() const; // dlansb
    F normInfinity() const;
    F normMaxEntry() const;
    F normOne() const;

    Vector<F,F>* eigenvalues(OrthogonalMatrix<F,Z> *&Q) const;

    virtual void printOn(ostream& =cerr) const ;
};

template<typename F,typename Z> BandMatrix<F,Z>* operator+(
const BandMatrix<F,Z> &B,const SymmetricBandMatrix<F,Z> &S) {
  return S+B;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const UpperHessenbergMatrix<F,Z> &H,const SymmetricBandMatrix<F,Z> &B) {
  return B+H;
}
template<typename F,typename Z> SymmetricBandMatrix<F,Z>* operator+(
const DiagonalMatrix<F,Z> &D,const SymmetricBandMatrix<F,Z> &B) {
  return B+D;
}
template<typename F,typename Z> SymmetricBandMatrix<F,Z>* operator+(
const SymmetricTridiagonalMatrix<F,Z> &S,const SymmetricBandMatrix<F,Z> &B)
{
  return B+S;
}
template<typename F,typename Z> BandMatrix<F,Z>* operator+(
const TridiagonalMatrix<F,Z> &T,const SymmetricBandMatrix<F,Z> &B) {
  return B+T;
}
template<typename F,typename Z> SymmetricMatrix<F,Z>* operator+(
const SymmetricMatrix<F,Z> &H,const SymmetricBandMatrix<F,Z> &B) {
  return B+H;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const UpperTrapezoidalMatrix<F,Z> &U,const SymmetricBandMatrix<F,Z> &B) {
  return B+U;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const LowerTrapezoidalMatrix<F,Z> &L,const SymmetricBandMatrix<F,Z> &B) {
  return B+L;
}
template<typename F,typename Z> SquareMatrix<F,Z>* operator+(
const Matrix<F,Z> &M,const SymmetricBandMatrix<F,Z> &B) {
  return B+M;
}

template<typename F,typename Z> BandMatrix<F,Z>* operator-(
  const BandMatrix<F,Z> &B,const SymmetricBandMatrix<F,Z> &S);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const UpperHessenbergMatrix<F,Z> &H,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> SymmetricBandMatrix<F,Z>* operator-(
  const DiagonalMatrix<F,Z> &D,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> SymmetricBandMatrix<F,Z>* operator-(
  const SymmetricTridiagonalMatrix<F,Z> &S,
  const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> BandMatrix<F,Z>* operator-(
  const TridiagonalMatrix<F,Z> &T,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> SymmetricMatrix<F,Z>* operator-(
  const SymmetricMatrix<F,Z> &H,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const UpperTrapezoidalMatrix<F,Z> &U,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const LowerTrapezoidalMatrix<F,Z> &L,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator-(
  const Matrix<F,Z> &M,const SymmetricBandMatrix<F,Z> &B);

template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const BandMatrix<F,Z> &B,const SymmetricBandMatrix<F,Z> &S);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const UpperHessenbergMatrix<F,Z>&,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const DiagonalMatrix<F,Z> &D,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const SymmetricTridiagonalMatrix<F,Z> &T,
  const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> BandMatrix<F,Z>* operator*(
  const TridiagonalMatrix<F,Z> &T,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SymmetricMatrix<F,Z> &S,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const UpperTrapezoidalMatrix<F,Z> &U,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const LowerTrapezoidalMatrix<F,Z> &U,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
  const SquareMatrix<F,Z> &S,const SymmetricBandMatrix<F,Z> &B);
template<typename F,typename Z> Matrix<F,Z>* operator*(
  const Matrix<F,Z> &S,const SymmetricBandMatrix<F,Z> &B);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> class SymmetricPositiveBandMatrix :
public SymmetricBandMatrix<F,Z> {
  private:
//  prevent public use of copy constructor
    explicit SymmetricPositiveBandMatrix(
    const SymmetricPositiveBandMatrix &B) :
    SymmetricBandMatrix<F,Z>(B.size(0),B.subDiags()) {
      SymmetricBandMatrix<F,Z>::copy(B);
    }
// Use copy instead
    SymmetricPositiveBandMatrix& operator=(
    const SymmetricPositiveBandMatrix& A) {
      SymmetricBandMatrix<F,Z>::copy(A); return *this;
    }
  public:
    SymmetricPositiveBandMatrix() : SymmetricBandMatrix<F,Z>() {;}
    explicit SymmetricPositiveBandMatrix(int n,int nb):
      SymmetricBandMatrix<F,Z>(n,nb) {}
    explicit SymmetricPositiveBandMatrix(int n,int nb,Z d):
      SymmetricBandMatrix<F,Z>(n,nb,d) {}
    virtual ~SymmetricPositiveBandMatrix() {;} 

    SymmetricPositiveBandMatrix& operator=(Z scalar) {
      SymmetricBandMatrix<F,Z>::fillWith(scalar); return *this;
    }

    SymmetricPositiveBandMatrix& operator+=(
    const SymmetricPositiveBandMatrix &B) {
      SymmetricBandMatrix<F,Z>::operator+=(B); return *this;
    }

    virtual SymmetricPositiveBandMatrix* operator+(
      const SymmetricPositiveBandMatrix &B) const;
    virtual SymmetricBandMatrix<F,Z>* operator+(
    const SymmetricBandMatrix<F,Z> &B) const {
      return SymmetricBandMatrix<F,Z>::operator+(B);
    }
    virtual SymmetricPositiveBandMatrix* operator+(
      const SymmetricPositiveTridiagonalMatrix<F,Z> &T) const;
    virtual SymmetricBandMatrix<F,Z>* operator+(
      const SymmetricTridiagonalMatrix<F,Z> &T) const;
    virtual SymmetricPositiveMatrix<F,Z>* operator+(
      const SymmetricPositiveMatrix<F,Z> &S) const;
    virtual SymmetricMatrix<F,Z>* operator+(const SymmetricMatrix<F,Z> &S)
      const;
    virtual BandMatrix<F,Z>* operator+(const BandMatrix<F,Z> &B) const {
      return SymmetricBandMatrix<F,Z>::operator+(B);
    }
    virtual SquareMatrix<F,Z>* operator+(
    const UpperHessenbergMatrix<F,Z> &H) const {
      return SymmetricBandMatrix<F,Z>::operator+(H);
    }
    virtual SymmetricBandMatrix<F,Z>* operator+(
    const DiagonalMatrix<F,Z> &D) const {
      return SymmetricBandMatrix<F,Z>::operator+(D);
    }
    virtual BandMatrix<F,Z>* operator+(
    const TridiagonalMatrix<F,Z> &T) const {
      return SymmetricBandMatrix<F,Z>::operator+(T);
    }
    virtual SquareMatrix<F,Z>* operator+(
    const UpperTrapezoidalMatrix<F,Z> &U) const {
      return SymmetricBandMatrix<F,Z>::operator+(U);
    }
    virtual SquareMatrix<F,Z>* operator+(
    const LowerTrapezoidalMatrix<F,Z> &L) const {
      return SymmetricBandMatrix<F,Z>::operator+(L);
    }
    virtual SquareMatrix<F,Z>* operator+(const Matrix<F,Z> &M) const {
      return SymmetricBandMatrix<F,Z>::operator+(M);
    }

//  LaPack:
    F equilibrate(Vector<F,F> &s,F &scond) const;
    F reciprocalConditionNumber() const; // dpbcon

    virtual void solve(const Vector<F,Z> &b,Vector<F,Z> &x,char trans='N')
      const;
    virtual void solve(const Matrix<F,Z> &b,Matrix<F,Z> &x,char side='L',
      char trans='N') const;
    virtual void printOn(ostream& =cerr) const ;
};

template<typename F,typename Z> SymmetricPositiveBandMatrix<F,Z>*
operator+(const SymmetricPositiveTridiagonalMatrix<F,Z> &T,
SymmetricPositiveBandMatrix<F,Z> &B) {
  return B+T;
}

template<typename F,typename Z> void testSymmetricPositiveBandMatrix(F,Z);

// Modified from gmd.H by John Trangenstein, 11/9/96
// gmd.H etc was unnecessarily repetitive; ie, not object-oriented.
#endif
