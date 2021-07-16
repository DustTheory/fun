#ifndef _VECTOR_H_
#define _VECTOR_H_ "$Header:$"
//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.
//
//      Lapack++ "Shared" Vector Double Class
//
//----------------------------  vector.h  ---------------------------
//    $Id: vector.h,v 1.42 2003/05/31 20:43:24 wolf Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  vector.h  ---------------------------
//
//modified from LAPACK++ file vd.h and deal.II/lac/include/lac/vector.h
//  by John Trangenstein, August 2009
//**********************************************************************
// Copyright 2009 John A. Trangenstein
//
// This software is made available for research and instructional use 
// only. 
// You may copy and use this software without charge for these 
// non-commercial purposes, provided that the copyright notice and 
// associated text is reproduced on all copies.  
// For all other uses (including distribution of modified versions), 
// please contact the author at
//   John A. Trangenstein
//   Department of Mathematics
//   Duke University
//   Durham, NC 27708-0320
//   USA
// or
//   johnt@math.duke.edu
// 
// This software is made available "as is" without any assurance that it
// is completely correct, or that it will work for your purposes.  
// Use the software at your own risk.
//**********************************************************************

#include <iostream>
#include <assert.h>
#include "Errors.H"
#include "MemoryDebugger.H"
#include "Tracer.H"

//this class used to be templated; we only used double for the template
//not in DEAL
class LAVector {
  protected:
#ifdef DEBUG
    static double undefined_;
#endif
    static double zero_,one_;
    double *data;
    int sz;
    bool owns_data;
  public:                                                            
    virtual void fillWith(double);
    void copy(const LAVector&); // calls memcpy

    LAVector() : data(0),sz(0),owns_data(true) {;}
    LAVector(int n) : data(0),sz(n),owns_data(true) { 
      if (sz>0) {
        data=OPERATOR_NEW_BRACKET(double,sz);
#ifdef DEBUG
        fillWith(undefined_);
#endif
      }
    }
    LAVector(int n,double scalar) : data(0),sz(n),owns_data(true) { 
      if (sz>0) {
        data=OPERATOR_NEW_BRACKET(double,n); 
        fillWith(scalar); 
      }
    }
    LAVector(const LAVector &m) : data(0),sz(0),owns_data(true) { 
      copy(m); 
    }
    LAVector(double *d,int s) : data(d),sz(s),owns_data(false) {}
    virtual ~LAVector() { 
      if (owns_data && data!=0) delete data; data=0; 
    }

    int resize(int n) {
      if (n==sz) return sz;
      CHECK_TEST(owns_data);
      if (data!=0) delete data; 
      sz=n; data=0;
      if (sz>0) {
        data=OPERATOR_NEW_BRACKET(double,sz);
        owns_data=true;
      }
      return sz;
    }

    double* addr() const { return data; }
    double* addr(int i) { 
      CHECK_BOUNDS(i,0,sz)
      return data+i; 
    }
    const double* addr(int i) const { 
      CHECK_BOUNDS(i,0,sz)
      return data+i; 
    }
    double& operator[](int i) { 
      CHECK_BOUNDS(i,0,sz)
      return data[i]; 
    }
    double operator[](int i) const { 
      CHECK_BOUNDS(i,0,sz)
      return data[i]; 
    }
    int size() const { return sz; }
    bool ownsData() const { return owns_data; }

    virtual LAVector& operator=(const LAVector &m) { 
      copy(m); return *this; 
    }
    LAVector& operator=(double s) { fillWith(s); return *this; } 
    LAVector& operator+=(const LAVector&); // calls lapack daxpy
    LAVector& operator-=(const LAVector&); // call lapack daxpy
    LAVector& operator*=(double); // calls lapack dscal
    LAVector& operator/=(double); // calls lapack dscal
    LAVector* operator+(const LAVector &B) const {
      LAVector *v=OPERATOR_NEW LAVector(*this); *v += B; return v;
    }
    LAVector* operator-(const LAVector &B) const {
      LAVector *v=OPERATOR_NEW LAVector(*this); *v -= B; return v;
    }
    LAVector* operator*(double d) const {
      LAVector *v=OPERATOR_NEW LAVector(*this); *v *= d; return v;
    }
    LAVector* operator/(double d) const {
      LAVector *v=OPERATOR_NEW LAVector(*this); *v /= d; return v;
    }
    double l2Norm() const; // calls lapack dnrm2

//  this function allows operator<< to access class members 
//  without being a friend
    virtual void printOn(ostream&) const;
};                                                                     
inline ostream& operator<<(ostream &s,const LAVector &m) {
  m.printOn(s); return s;
}

//compare to DEAL Vector, which does not use LaPack
class Vector {
  protected:
    LAVector d;
  public:
    Vector() { }
    explicit Vector(int n) : d(n) { }
    Vector(int n,double scalar) : d(n,scalar) { } 
    Vector(const Vector &m) : d(m.d) { }
    Vector(double *r,int s) : d(r,s) {}
    virtual ~Vector() { d.LAVector::~LAVector(); }
    Vector* clone() const { return OPERATOR_NEW Vector(d.size()); }
    Vector* cloneEmpty() const { return OPERATOR_NEW Vector; }
    Vector* subVector(int offset,int sz) {
      CHECK_TEST(offset+sz<=d.size());
      return OPERATOR_NEW Vector(d.addr()+offset,sz);
    }
    int size() const { return d.size(); }
    bool ownsData() const { return d.ownsData(); }
    void reinit(int n,bool fast=false) {
      d.resize(n);
      if (!fast) d.fillWith(0.);
    }
    double* addr() { return d.addr(); }
    const double* addr() const { return d.addr(); }
    double& operator[](int i) { return d[i]; }
    double operator[](int i) const { return d[i]; }
    Vector& operator=(const Vector &v) {
      CHECK_SAME(d.size(),v.d.size());
      d.copy(v.d);
      return *this;
    }
    Vector& operator=(double s) { d=s; return *this; }
    Vector& operator+=(const Vector &v) { 
      d+=v.d;
      return *this;
    }
    Vector& operator-=(const Vector &v) { 
      d-=v.d;
      return *this;
    }
    void scale(double a) { d*=a; }
    void equ(double,const Vector&); // this=r*vector; calls lapack
    double operator*(const Vector&) const; // dot product; call lapack
    double linftyNorm() const;
    void add(double,const Vector&); //this+=r*vector; call lapack
    void swapWith(Vector&);
    double l2Norm() const { return d.l2Norm(); }
    bool allZero() const { return linftyNorm()<=0.; }
    void printOn(ostream &os) const {
      os << "Vector:" << endl;
      d.printOn(os);
    }
};
#endif 
