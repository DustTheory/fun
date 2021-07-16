#ifndef _MATRIX_H
#define _MATRIX_H
#include <iostream>
#include <string.h>
#include "Errors.H"

class Vector;
class Matrix {
  private:
    int sz,sizes[2];
    float *data;

    void setSizes(int m,int n) { 
      if (data) delete [] data; sz=m*n; data=new float[sz];
      sizes[0]=m; sizes[1]=n; 
    }
  protected:
  public:
    Matrix() : sz(0),data(0) { sizes[0]=0; sizes[1]=0; }
    Matrix(const int &m,const int &n): data(0) {
      ASSERT(m>0); ASSERT(n>0); setSizes(m,n);
    }
    Matrix(const Matrix &G) : data(0) { copy(G); }
    virtual ~Matrix() {if (data) delete [] data; data=0;}

    int size(int d) const { ASSERT(d==0 || d==1); return sizes[d]; }
    int sameSizeAs(const Matrix &G) const {
      return sizes[0]==G.sizes[0] && sizes[1]==G.sizes[1];
    }
    virtual float& operator()(int i,int j) {
      ASSERT(i>=0); ASSERT(i<sizes[0]);
      ASSERT(j>=0); ASSERT(j<sizes[1]);
      return data[i+j*sizes[0]];
    }
    virtual float operator()(int i,int j) const {
      ASSERT(i>=0); ASSERT(i<sizes[0]);
      ASSERT(j>=0); ASSERT(j<sizes[1]);
      return data[i+j*sizes[0]];
    }
    Matrix& operator=(const Matrix& A) { copy(A); return *this; }

    virtual void fill(float);
    virtual void copy(const Matrix&);
    virtual void copyInto(Matrix&) const ;
    Matrix& operator+=(const Matrix&);
    Matrix& operator-=(const Matrix&);
    Matrix& operator*=(float);
    Matrix& operator/=(float);

    virtual void printOn(ostream& =cerr) const ;

//these allocate memory for the result
//  to avoid memory leaks, the address of the result must be deleted
    Matrix* operator+(const Matrix&) const ;
    Matrix* operator-(const Matrix&) const ;
    Matrix* operator*(float d) const {
      Matrix *prod=new Matrix(*this); (*prod)*=d; return prod;
    }
    Matrix* operator/(float d) const {
      ASSERT(d!=0.);
      Matrix *quot=new Matrix(*this); (*quot)*=(1./d); return quot;
    }
    Vector* operator*(const Vector&) const;
    Matrix* operator*(const Matrix&) const;
    Matrix* transpose() const ;
    Vector* solve(const Vector&) const;

};
inline ostream& operator<<(ostream &s,const Matrix &G) {
  G.printOn(s); return s;
}
#endif
