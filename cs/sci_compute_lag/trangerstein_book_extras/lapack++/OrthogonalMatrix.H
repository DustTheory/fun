#ifndef _ORTHOGONAL_MATRIX_H_
#define _ORTHOGONAL_MATRIX_H_
#include "TrapezoidalMatrix.H"

template<typename F,typename Z> class OrthogonalMatrix : public Matrix<F,Z> {
  private:
    explicit OrthogonalMatrix(const OrthogonalMatrix &Q) : 
    Matrix<F,Z>(Q.size(0),Q.size(1)) {
      Matrix<F,Z>::copy(Q);
    }
    OrthogonalMatrix& operator=(const OrthogonalMatrix &A) {
      Matrix<F,Z>::copy(A); return *this;
    }
//  prevent public use of some Matrix functions
    virtual void fillWith(Z s) { Matrix<F,Z>::fillWith(s); }
    virtual OrthogonalMatrix& operator=(Z scalar) {
      OBSOLETE("inappropriate for this class"); return *this;
    }
    F equilibrate(Vector<F,F>&,Vector<F,F>&,F&,F&) const {
      OBSOLETE("inappropriate for this class");
      return 0.;
    }
    void scale(char,int,int,F,F) {
      OBSOLETE("inappropriate for this class");
    }
    void set(char,Z,Z) {
      OBSOLETE("inappropriate for this class");
    }
    virtual Matrix<F,Z>* operator+(const Matrix<F,Z>&) const {
      OBSOLETE("inappropriate for this class");
      return 0;
    }
    virtual Matrix<F,Z>* operator-(const Matrix<F,Z>&) const {
      OBSOLETE("inappropriate for this class");
      return 0;
    }
  protected:
  public:
    OrthogonalMatrix() : Matrix<F,Z>() {;}
    explicit OrthogonalMatrix(int m,int n);
    virtual ~OrthogonalMatrix() {;}

    virtual void resize(int,int);
    virtual void resize(const Matrix<F,Z> &M) {
      resize(M.size(0),M.size(1));
    }
    virtual void copy(const Matrix<F,Z> &M) {
      Matrix<F,Z>::copy(M);
    }

    virtual OrthogonalMatrix<F,Z>& operator*=(Z d) { 
      CHECK_TEST(abs(d)==1.);
      Matrix<F,Z>::operator*=(d);
      return *this;
    }
    virtual OrthogonalMatrix<F,Z>& operator/=(Z d) {
      CHECK_TEST(abs(d)==1.);
      Matrix<F,Z>::operator/=(d);
      return *this;
    }
//  use gemv or gemm instead of transpose
//  OrthogonalMatrix* transpose() const;
//  OrthogonalMatrix* conjugateTranspose() const;

//  in the solve routines, trans!='N' is assumed to mean conjugate transpose
    virtual void solve(const Matrix<F,Z> &B,Matrix<F,Z> &X,char side='L',
      char trans='N') const;
    virtual void solve(const Vector<F,Z> &b,Vector<F,Z> &x,char trans='N')
      const;
    void printOn(ostream& = cerr) const ;
};

template<typename F,typename Z> void testOrthogonalMatrix(F,Z);
#endif 
