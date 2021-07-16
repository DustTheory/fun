#ifndef _VECTOR_H_
#define _VECTOR_H_
#include <iostream>
#include <math.h>
#include "Errors.H"

//This class is not polymorphic. It is specific to type float.
//The template class Vector<F,Z> is polymorphic; see
//  lapack++/Vector.H
class Vector {
  private:
    float *data;
    int sz;

    explicit Vector(const Vector&);
    void fillWith(float);               // x.fill(0)
  public:
    Vector() : data(0), sz(0) {;}
    explicit Vector(int m) : sz(m) { ASSERT(m>0); data=new float[m]; }
    explicit Vector(int m,float v) : sz(m) {
      ASSERT(m>0); data=new float[m]; fillWith(v);
    }
    virtual ~Vector() { delete data; data=0; sz=0; }

    Vector& operator=(float v) { fillWith(v); return *this; }
    int size() const { return sz; }

//  operator[] is inline: ASSERT may trace back to another line number
    float operator[](int i) const {
      ASSERT(i>=0 && i<sz); return data[i];
    }
    float& operator[](int i) { ASSERT(i>=0 && i<sz); return data[i]; }

//real linear space operations
    Vector& operator+=(const Vector&);   // x+=y
    Vector& operator-=(const Vector&);   // x-=y
    Vector& operator*=(float);           // x*=alpha
    Vector& operator/=(float);           // x/=alpha

//these allocate memory for the result
//  to avoid memory leaks, the address of the result must be deleted
    Vector& operator+(const Vector&) const; // x+y
    Vector& operator-(const Vector&) const; // x-y
    Vector& operator*(float) const;         // x*y
    Vector& operator/(float) const;         // x/y

//blas1 routines
    int amax() const;               // i=x.imax()
    float asum() const;             // alpha=x.asum()
    void axpy(float,const Vector&); // y.axpy(alpha,x)
    void copy(const Vector&);       // y.copy(x)
    void copyFrom(int,Vector&) const;   // y.copyFrom(n,x)
    float dot(const Vector&) const; // alpha=y.dot(x)
    float nrm2() const;             // alpha=x.nrm2()
    void swap(Vector&);             // y.swap(x)

    void printOn(ostream& = cout) const;
};
inline ostream& operator<<(ostream &os,const Vector &x) {
  x.printOn(os); return os;
}
#endif
