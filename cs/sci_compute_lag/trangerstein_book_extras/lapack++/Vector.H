#ifndef _VECTOR_H_
#define _VECTOR_H_    
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <string.h>
#include "Errors.H"
/*
using namespace std;
# ifndef ASSERT
# define ASSERT(s) { if (!(s)) { cerr << "assertion " << __STRING(s) << " failed in file " << __FILE__ << " at line " << __LINE__ << endl; abort(); } }
# endif
# ifndef OBSOLETE
# define OBSOLETE(s) { cerr << s << " in file " << __FILE__ << " at line " << __LINE__ << endl; abort(); }
# endif
# ifdef DEBUG
#   define CHECK_TEST(s) { ASSERT(s) }
#   define CHECK_POINTER(p) { ASSERT(p!=0) }
#   define CHECK_NONNEGATIVE(r) { ASSERT(r>=0) }
#   define CHECK_POSITIVE(r) { ASSERT(r>0) }
#   define CHECK_BOUNDS(v,l,u) { ASSERT(v>=l && v<u) }
#   define CHECK_NONZERO(q) { ASSERT(q!=0) }
#   define CHECK_SAME(v1,v2) { ASSERT((v1)==(v2)) }
# else
#   define CHECK_TEST(s) ;
#   define CHECK_POINTER(p) ;
#   define CHECK_NONNEGATIVE(r) ;
#   define CHECK_POSITIVE(r) ;
#   define CHECK_BOUNDS(v,l,u) ;
#   define CHECK_NONZERO(q) ;
#   define CHECK_SAME(v1,v2) ;
# endif
*/

#include "MemoryDebugger.H"
/*
# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif
*/

//#include "Arch.H"
#ifndef F77NAME
#define F77NAME(x) x ## _
#endif
#ifndef F77_NAME
#define F77_NAME(x) x ## __
#endif

//compare to TNT class Vector<T> in file tnt_vec.h
//
//we use two template parameters:
//  C for Vector entries: C can be real or complex
//  R for modulus of C: R must be real
//the possibilities are <R,C> =
//    <double,double>
//    <double,complex<double> >
//    <float,float>
//    <float,complex<float> >
template<typename R,typename C> class Vector {
  private:
    C *data;
    int sz;

//  prevent public use of copy constructor
    explicit Vector(const Vector &m) : data(0),sz(m.size()) {
      if (sz>0) {
        data=OPERATOR_NEW_BRACKET(C,sz);
#ifdef DEBUG
        fillWith(undefined_);
#endif
        copy(m);
      }
    }
//  use Vector::copy instead of Vector::operator=
    Vector& operator=(const Vector &m) { 
      resize(m.size());
      copy(m); 
      return *this;
    }
  protected:
    virtual void fillWith(C);
  public:                                                            
#ifdef DEBUG
    static const C undefined_;
#endif
    static const C mone_,one_,zero_;

    Vector() : data(0),sz(0) {;}
    explicit Vector(int n) : sz(n) { 
      assert(sz>0);
      data = OPERATOR_NEW_BRACKET(C,sz);
#ifdef DEBUG
      fillWith(undefined_);
#endif
    }
    explicit Vector(int n,C scalar) : data(0),sz(n) { 
      assert(sz>0);
      data=OPERATOR_NEW_BRACKET(C,n);
      fillWith(scalar); 
    }
    virtual ~Vector() { if (data) delete data; data=0; sz=0; }

    void resize(int n) {
      if (n==sz) return;
      assert(n >= 0);
      if (data) delete data; 
      data=0;
      sz=n;
      if (sz>0) {
        data=OPERATOR_NEW_BRACKET(C,sz);
#ifdef DEBUG
        fillWith(undefined_);
#endif
      }
    }
    virtual Vector& operator=(C scalar) { fillWith(scalar); return *this; }

    int size() const { return sz; }
    C* addr() const { return data; }
    C* addr(int i) { 
      CHECK_BOUNDS(i,0,sz)
      return data+i; 
    }
    const C* addr(int i) const { 
      CHECK_BOUNDS(i,0,sz)
      return data+i; 
    }
    C& operator[](int i) { 
      CHECK_BOUNDS(i,0,sz)
      return data[i]; 
    }
    C operator[](int i) const { 
      CHECK_BOUNDS(i,0,sz)
      return data[i]; 
    }

    virtual Vector& operator+=(const Vector &y) { // x+=y
      CHECK_SAME(sz,y.sz);
      axpy(one_,y); return *this;
    }
    virtual Vector& operator-=(const Vector &y) { // x-=y
      CHECK_SAME(sz,y.sz);
      axpy(mone_,y); return *this;
    }
    virtual Vector& operator*=(C alpha) {         // x*=alpha
      scal(alpha); return *this;
    }
    virtual Vector& operator/=(C alpha) {         // x/=alpha
      scal(one_/alpha); return *this;
    }
//  these allocate memory for the result
//  to avoid memory leaks, when no longer needed delete return val
    virtual Vector* operator+(const Vector &b) const {        // x+y
      Vector *v=OPERATOR_NEW Vector(*this); *v += b; return v;
    }
    virtual Vector* operator-(const Vector &b) const {        // x-y
      Vector *v=OPERATOR_NEW Vector(*this); *v -= b; return v;
    }
    virtual Vector* operator*(C d) const {                    // x*alpha
      Vector *v=OPERATOR_NEW Vector(*this); *v *= d; return v;
    }
    virtual Vector* operator/(C d) const {                    // x/alpha
      Vector *v=OPERATOR_NEW Vector(*this); *v /= d; return v;
    }

//  BLAS routines
    int amax() const;                   // index of first entry of max modulus
    R asum() const;                     // sum of moduli of entries
    void axpy(C alpha,const Vector &x); // this += x * alpha
    void copy(const Vector &x) {        // copy x into this
      CHECK_SAME(sz,x.sz);
      memcpy(data,x.data,sz*sizeof(C)); //faster than Blas copy
    }
    void copyFrom(int n,const Vector &x) {
      n=min(n,min(sz,x.sz));
      memcpy(data,x.data,n*sizeof(C));
    }
    C dot(const Vector &x) const;       // C real: x . this
    C dotc(const Vector &x) const;      // C complex: hermitian(x) . this
    C dotu(const Vector &x) const;      // C complex: transpose(x) . this
    R nrm2() const;                     // || this ||_2
    void rot(Vector &x,R c,R s);        // apply plane rotation to [ x \\ y ]
    void scal(C s);                     // this *= s
    void swap(Vector &x);               // swap entries between x and this

//  LaPack routines
    void largv(Vector &y,Vector<R,R> &c); // generate plane rotation
    void lartv(Vector &y,const Vector<R,R> &c,const Vector &s);
      //apply rotation; similar to BLAS rot, but numerically better

    virtual void printOn(ostream&) const;
};                                                                     

// _VECTOR_H_
#endif 
