#ifndef _MATRIX_H
#define _MATRIX_H "$Header:$"
//
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge 
// National Laboratory) nor the Authors make any representations about 
// the suitability of this software for any purpose.  This software is 
// provided ``as is'' without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.

//----------------------------  full_matrix.h  ---------------------------
//    $Id: full_matrix.h,v 1.60 2003/06/02 22:10:09 wolf Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  full_matrix.h  ---------------------------

//Modified from LaPack++ gmc.h and deal.II/lac/include/lac/full_matrix.h
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
#include <string.h>
#include "SolverArrays.H"
#include "Vector.H"

class LAMatrix : public LAVector {
  protected:
    int dim[2];

// disable public use of LAVector functions not appropriate for LAMatrix
    double& operator[](int i) { return LAVector::operator[](i); }
    double operator[](int i) const { return LAVector::operator[](i); }
    int resize(int n) { return LAVector::resize(n); }
//  can't resize from a LAVector
    void copy(const LAVector&) { assert(0); }
//  Use copy instead, to avoid confusion with copy constructor
    LAMatrix& operator=(const LAVector& A) {
      OBSOLETE("disallowed");
      return *this;
    }
    void setDimensions(int m,int n) { dim[0]=m; dim[1]=n; }
  public:
    void copy(const LAMatrix& A) {
      LAVector::copy(A); setDimensions(A.size(0),A.size(1));
    }

    LAMatrix() : LAVector() { dim[0]=0; dim[1]=0; }
    LAMatrix(int m,int n): LAVector(m*n) { dim[0]=m; dim[1]=n; }
    LAMatrix(int m,int n,double d): LAVector(m*n,d) { dim[0]=m; dim[1]=n; }
    LAMatrix(const LAMatrix &G) : LAVector(G) {
      setDimensions(G.dim[0],G.dim[1]);
    }
    virtual ~LAMatrix() {;}

    int size(int d) const { return dim[d]; }
    double* addr() { return LAVector::addr(); }
    const double* addr() const { return LAVector::addr(); }
    virtual double* addr(int i, int j) {
      CHECK_BOUNDS(i,0,dim[0])
      CHECK_BOUNDS(j,0,dim[1])
      return LAVector::addr(i+j*dim[0]);
    }
    virtual const double* addr(int i, int j) const {
      CHECK_BOUNDS(i,0,dim[0])
      CHECK_BOUNDS(j,0,dim[1])
      return LAVector::addr(i+j*dim[0]);
    }
    virtual double& operator()(int i,int j) { return *addr(i,j); }
    virtual double operator()(int i,int j) const { return *addr(i,j);}
    void resize(int m,int n) { resize(m*n); setDimensions(m,n); }
    void resize(const LAMatrix& A) { resize(A.size(0),A.size(1)); }

    LAMatrix& operator+=(const LAMatrix &M) {
      CHECK_SAME(size(0),M.size(0))
      CHECK_SAME(size(1),M.size(1))
      *static_cast<LAVector*>(this) += M; return *this;
    }
//  this has a return type of LAMatrix* to remind us that the return
//  pointer must be deleted by the calling procedure
    LAMatrix* operator+(const LAMatrix &M) const {
      CHECK_SAME(size(0),M.size(0))
      CHECK_SAME(size(1),M.size(1))
      LAMatrix *sum=new LAMatrix(*this);
      *static_cast<LAVector*>(sum)+=M; return sum;
    }
    LAMatrix& operator-=(const LAMatrix &M) {
      CHECK_SAME(size(0),M.size(0))
      CHECK_SAME(size(1),M.size(1))
      *static_cast<LAVector*>(this)-=M; return *this;
    }
    LAMatrix* operator-(
    const LAMatrix &M) const {
      CHECK_SAME(size(0),M.size(0))
      CHECK_SAME(size(1),M.size(1))
      LAMatrix *dif=new LAMatrix(*this);
      *static_cast<LAVector*>(dif)-=M; return dif;
    }
    LAMatrix& operator*=(double d) {
      *static_cast<LAVector*>(this)*=d; return *this;
    }
    LAMatrix* operator*(double d) const {
      LAMatrix *prod=new LAMatrix(*this);
      *static_cast<LAVector*>(prod)*=d; return prod;
    }
    LAMatrix& operator/=(double d) {
      *static_cast<LAVector*>(this)/=d; return *this;
    }
    LAMatrix* operator/(double d) const {
      LAMatrix *quot=new LAMatrix(*this);
      *static_cast<LAVector*>(quot)/=d; return quot;
    }

    void copyInto(LAMatrix&) const ;
    LAMatrix* transpose() const ;
    void printOn(ostream& =cerr) const ;

};
inline ostream& operator<<(ostream &s,const LAMatrix &G) {
  G.printOn(s); return s;
}

int LaEnvBlockSize(const char*,const LAMatrix&);

LAMatrix* productOf(const LAMatrix&,const LAMatrix&);
LAMatrix* factorNSolve(const LAMatrix &,
  const LAMatrix &);
void switchColumns(LAMatrix &,int,int);
void switchRows(LAMatrix &,int,int);

//compare to DEAL FullMatrix, which does not use LaPack
class Matrix : public SolverMatrix {
  protected:
    LAMatrix M;
  public:
    Matrix() { }
    Matrix(int m,int n) : M(m,n) { }
    Matrix(int m,int n,double d): M(m,n,d) { }
    Matrix(const Matrix &G) : M(G.M) { }
    void copy(const LAMatrix &A) { M.copy(A); }
    virtual ~Matrix() {}
    virtual int size(int d) const { M.size(d); }
    double* addr() { return M.addr(); }
    const double* addr() const { return M.addr(); }
    Matrix& operator=(double s) { M.fillWith(s); return *this; }
    virtual double& operator()(int i,int j) { return M(i,j); }
    virtual double operator()(int i,int j) const { return M(i,j);}
    void resize(int m,int n) { M.resize(m,n); }
    void copyInto(Matrix &A) const { M.copyInto(A.M); }
    void vmult(Vector&,const Vector&) const;
    void transposeVmult(Vector&,const Vector&) const;
    void vmultAdd(Vector&,const Vector&) const;
    void transposeVmultAdd(Vector&,const Vector&) const;
    double residual(Vector &r,const Vector &x,const Vector &b) const;
    void preconditionJacobi(Vector&,const Vector&,double=1.) 
      const;
    void PSOR(Vector&,const NumPtr<int>&,const NumPtr<int>&,double=1.) 
      const;
    void TPSOR(Vector&,const NumPtr<int>&,const NumPtr<int>&,double) 
      const;
    void preconditionSOR(Vector&,const Vector&,double=1.) const;
    void preconditionSSOR(Vector&,const Vector&,double=1.) const;
    void preconditionTSOR(Vector&,const Vector&,double) const;
    double matrixNormSquare(const Vector &v) const {
      Vector *Mv=v.clone();
      vmult(*Mv,v);
      double prod=(*Mv)*v;
      delete Mv;
      return prod;
    }
    Matrix* clone() const { 
      return OPERATOR_NEW Matrix(size(0),size(1)); 
    }
    Matrix* cloneEmpty() const { return OPERATOR_NEW Matrix; }
    void printOn(ostream &os) const {
      os << "Matrix:" << endl;
      M.printOn(os);
    }
};

//replaces LaPack gmc.H gmd.H gmf.H gmi.H gmli.H
// Modified from gmd.H by John Trangenstein, 11/9/96
#endif
