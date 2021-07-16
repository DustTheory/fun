#ifndef _MATRIX_H
#define _MATRIX_H
#include <iostream>
#include "Vector.H"

//compare to TNT class Array2D<T> in file tnt_array2d.h
//
//we use two template parameters:
//  C for Vector entries: C can be real or complex
//  R for modulus of C: R must be real
//the possibilities are <R,C> =
//    <double,double>
//    <double,complex<double> >
//    <float,float>
//    <float,complex<float> >
template<typename R,typename C> class Matrix {
  private:
    Vector<R,C> *data;
    int dim[2];

//  prevent public use of copy constructor
    explicit Matrix(const Matrix &M) {
      setDimensions(M.dim[0],M.dim[1]);
      data=OPERATOR_NEW Vector<R,C>(dim[0]*dim[1]);
      data->copy(*M.data);
    }
//  use copy instead of operator=
    Matrix& operator=(const Matrix& A) { // do not make this virtual
      resize(A.size(0),A.size(1));
      data->resize(dim[0]*dim[1]);
      data->copy(*A.data);
      return *this;
    }
  protected:
    void setDimensions(int m,int n) { dim[0]=m; dim[1]=n; }
    virtual void fillWith(C scalar) { *data=scalar; }
  public:
#ifdef DEBUG
    static const C undefined_;
#endif
    static const C mone_,one_,zero_;

    Matrix() : data(0) { setDimensions(0,0); }
    explicit Matrix(int m,int n) {
      data=OPERATOR_NEW Vector<R,C>(m*n);
      setDimensions(m,n);
    }
    explicit Matrix(int m,int n,C d) {
      data=OPERATOR_NEW Vector<R,C>(m*n,d);
      setDimensions(m,n);
    }
    virtual ~Matrix() { if (data!=0) delete data; data=0; }

    virtual Matrix& operator=(C scalar) {
      *data=scalar; return *this;
    }

    int size(int d) const { return dim[d]; }
    C* addr() { return data->addr(); }
    const C* addr() const { return data->addr(); }
    virtual C* addr(int i, int j) {
      CHECK_BOUNDS(i,0,dim[0])
      CHECK_BOUNDS(j,0,dim[1])
      return data->addr(i+j*dim[0]);
    }
    virtual const C* addr(int i, int j) const {
      CHECK_BOUNDS(i,0,dim[0])
      CHECK_BOUNDS(j,0,dim[1])
      return data->addr(i+j*dim[0]);
    }
    virtual C& operator()(int i,int j) { return *addr(i,j); }
    virtual C operator()(int i,int j) const { return *addr(i,j);}
    virtual void resize(int m,int n) { 
      setDimensions(m,n);
      if (data!=0) data->resize(m*n);
      else data=OPERATOR_NEW Vector<R,C>(m*n);
    }
    virtual void resize(const Matrix& A) { resize(A.size(0),A.size(1)); }
    virtual void copy(const Matrix& A) {
      CHECK_SAME(size(0),A.size(0));
      CHECK_SAME(size(1),A.size(1));
      data->copy(*A.data);
    }

    Matrix& operator+=(const Matrix &M) { // do not make this virtual
      CHECK_SAME(size(0),M.size(0))
      CHECK_SAME(size(1),M.size(1))
      data->operator+=(*M.data);
      return *this;
    }
    Matrix& operator-=(const Matrix &M) { // do not make this virtual
      CHECK_SAME(size(0),M.size(0))
      CHECK_SAME(size(1),M.size(1))
      data->operator-=(*M.data);
      return *this;
    }
    virtual Matrix& operator*=(C d) {
      data->operator*=(d); return *this;
    }
    virtual Matrix& operator/=(C d) {
      data->operator/=(d); return *this;
    }
//  these return Matrix* to remind us that the return value
//    must be deleted by the calling procedure, to avoid memory leaks
    virtual Matrix* operator+(const Matrix &M) const {
      CHECK_SAME(size(0),M.size(0))
      CHECK_SAME(size(1),M.size(1))
      Matrix *S=OPERATOR_NEW Matrix(*this);
      S->data->operator+=(*M.data); return S;
    }
    virtual Matrix* operator-(const Matrix &M) const {
      CHECK_SAME(size(0),M.size(0))
      CHECK_SAME(size(1),M.size(1))
      Matrix *D=OPERATOR_NEW Matrix(*this);
      D->data->operator-=(*M.data); return D;
    }
    virtual Matrix* operator*(C d) const {
      Matrix *P=OPERATOR_NEW Matrix(*this);
      P->data->operator*=(d); return P;
    }
    virtual Matrix* operator/(C d) const {
      Matrix *Q=OPERATOR_NEW Matrix(*this);
      Q->data->operator/=(d); return Q;
    }
    virtual Matrix<R,C>* operator*(const Matrix<R,C>&) const;
    virtual Vector<R,C>* operator*(const Vector<R,C>&) const;

//  use gemv or gemm instead of transpose:
//  virtual Matrix* transpose() const ;
//  virtual Matrix* conjugateTranspose() const ;
    void interchangeColumns(int,int);
    void interchangeRows(int,int);

//  BLAS routines
    void gemv(C alpha,const Vector<R,C> &x,C beta,
      Vector<R,C> &y,char trans='N') const; // y=op(THIS)*x*alpha+y*beta
    void ger(C alpha,const Vector<R,C> &x,const Vector<R,C> &y);
      // for real C: THIS+=x*alpha*y^T
    void gerc(C alpha,const Vector<R,C> &x,const Vector<R,C> &y);
      // for complex C: THIS+=x*alpha*y^H
    void geru(C alpha,const Vector<R,C> &x,const Vector<R,C> &y);
      // for complex C: THIS+=x*alpha*y^T

//  Blas3 routines
    void gemm(C alpha,const Matrix<R,C> &A,const Matrix<R,C> &B,C beta,
      char transa='N',char transb='N'); // THIS=op(A)*alpha*op(B)+THIS*beta

//  LaPack routines
    R equilibrate(Vector<R,R> &r,Vector<R,R> &c,R &rowcnd,R &colcnd) const;
    virtual void copyFrom(char uplo,int m,int n,const Matrix &A);//copy from A
    void scale(char type,int kl,int ku,R denominator,R numerator);
    void set(char uplo,C offdiag,C diag); //set diagonal, off-diag entries
    int lastNonzeroColumn() const;
    int lastNonzeroRow() const;
    virtual R normFrobenius() const;
    virtual R normInfinity() const;
    virtual R normMaxEntry() const;
    virtual R normOne() const;
    virtual void solve(const Matrix&,Matrix&,char side='L',char trans='N')
      const;
    virtual void solve(const Vector<R,C>&,Vector<R,C>&,char trans='N')
      const;

    virtual void printOn(ostream&) const ;
};
template<typename R,typename C> void testMatrix(R f,C scalar);
#endif
