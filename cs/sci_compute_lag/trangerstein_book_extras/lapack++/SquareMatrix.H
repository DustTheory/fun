#ifndef _SQUARE_MATRIX_H
#define _SQUARE_MATRIX_H
#include <complex>
#include <iostream>
//#include <string.h>
#include "Matrix.H"

template<typename R,typename C> class SquareMatrix : public Matrix<R,C> {
  private:
//  prevent some Matrix functions
    void setDimensions(int m,int n) { // SquareMatrix requires n == m
      CHECK_SAME(m,n)
      Matrix<R,C>::setDimensions(m,n);
    }
    virtual void resize(int m,int n) { // SquareMatrix requires n == m
      CHECK_SAME(m,n)
      Matrix<R,C>::resize(m,n); 
    }
    virtual void resize(const Matrix<R,C> &A) { // A must be SquareMatrix 
      int m=A.size(0),n=A.size(1);
      CHECK_SAME(m,n)
      Matrix<R,C>::resize(m,n);
    }
//  prevent copy constructor
    explicit SquareMatrix(const SquareMatrix &S) : 
      Matrix<R,C>(S.size(0),S.size(1)) { Matrix<R,C>::copy(S);}
    SquareMatrix& operator=(const SquareMatrix &A) { // do not make virtual
      resize(A.size(0)); Matrix<R,C>::copy(A); return *this;
    }
  protected:
    void setDimension(int n) { setDimensions(n,n); }
  public:
    SquareMatrix() : Matrix<R,C>() { }
    explicit SquareMatrix(int n): Matrix<R,C>(n,n) {;}
    explicit SquareMatrix(int n,C d): Matrix<R,C>(n,n,d) {;}
    virtual ~SquareMatrix() {;}

    virtual SquareMatrix& operator=(C scalar) {
      Matrix<R,C>::fillWith(scalar); return *this; 
    }
    virtual void resize(int n) { Matrix<R,C>::resize(n,n); }
    virtual void resize(const SquareMatrix &A) { Matrix<R,C>::resize(A); }
    SquareMatrix& operator+=(const Matrix<R,C> &A) {
      Matrix<R,C>::operator+=(A); return *this;
    }
    SquareMatrix& operator-=(const Matrix<R,C> &A) {
      Matrix<R,C>::operator-=(A); return *this;
    }
    virtual SquareMatrix& operator*=(C scalar) {
      Matrix<R,C>::operator*=(scalar);
      return *this;
    }
    virtual SquareMatrix& operator/=(C scalar) {
      Matrix<R,C>::operator/=(scalar);
      return *this;
    }
//  these return SquareMatrix* etc. to remind us that the return value
//    must be deleted by the calling procedure, to avoid memory leaks
    virtual SquareMatrix<R,C>* operator+(const SquareMatrix<R,C> &M)
    const {
      SquareMatrix *S=OPERATOR_NEW SquareMatrix(*this); *S+=M; return S;
    }
    virtual SquareMatrix<R,C>* operator+(const Matrix<R,C> &M) const {
      SquareMatrix *S=OPERATOR_NEW SquareMatrix(*this); *S+=M; return S;
    }
    virtual SquareMatrix<R,C>* operator-(const SquareMatrix<R,C> &M)
    const {
      SquareMatrix *S=OPERATOR_NEW SquareMatrix(*this); *S-=M; return S;
    }
    virtual SquareMatrix<R,C>* operator-(const Matrix<R,C> &M) const {
      SquareMatrix *S=OPERATOR_NEW SquareMatrix(*this); *S-=M; return S;
    }
    virtual SquareMatrix* operator*(C scalar) const {
      SquareMatrix *S=OPERATOR_NEW SquareMatrix(*this);
      *S*=scalar; return S;
    }
    virtual SquareMatrix* operator/(C scalar) const {
      SquareMatrix *S=OPERATOR_NEW SquareMatrix(*this);
      *S/=scalar; return S;
    }
    virtual SquareMatrix<R,C>* operator*(const SquareMatrix<R,C>&) const;
    virtual Matrix<R,C>* operator*(const Matrix<R,C> &M) const;
    virtual Vector<R,C>* operator*(const Vector<R,C> &v) const {
      return Matrix<R,C>::operator*(v);
    }

//  use Matrix::gemv or Matrix::gemm instead of transpose:
//  virtual SquareMatrix* transpose() const ;
//  virtual SquareMatrix* conjugateTranspose() const ;
  
//  LaPack:
    R reciprocalConditionNumber(char norm) const;
//  SquareMatrix* inverse() const;
    Vector<R,complex<R> >* eigenvalues(SquareMatrix<R,complex<R> > *&V,
      SquareMatrix<R,complex<R> > *&U) const; // V=left eigenvector,U=right
    virtual void printOn(ostream&) const;
};
//if you need a determinant, get it from an LUFactor

//  these return SquareMatrix* etc. to remind us that the return value
//    must be deleted by the calling procedure, to avoid memory leaks
template<typename R,typename C> SquareMatrix<R,C>* operator+(
const Matrix<R,C> &M,const SquareMatrix<R,C> &S) {
  SquareMatrix<R,C> *T=OPERATOR_NEW SquareMatrix<R,C>(S.size(0));
  T->copy(S); *T+=M; return T;
}
template<typename R,typename C> SquareMatrix<R,C>* operator-(
const Matrix<R,C> &M,const SquareMatrix<R,C> &S) {
  SquareMatrix<R,C> *T=OPERATOR_NEW SquareMatrix<R,C>(S.size(0));
  T->copy(M); *T-=S; return T;
}
template<typename R,typename C> Matrix<R,C>* operator*(
  const Matrix<R,C> &M,const SquareMatrix<R,C> &S);

template<typename R,typename C> void testSquareMatrix(R,C);
#endif
