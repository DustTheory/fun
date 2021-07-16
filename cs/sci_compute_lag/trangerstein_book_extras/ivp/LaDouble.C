// "$Header:$"
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
// Neither the Institutions (University of Tennessee, and Oak Ridge National
// Laboratory) nor the Authors make any representations about the suitability 
// of this software for any purpose.  This software is provided ``as is'' 
// without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.

//----------------------------  vector.cc  ---------------------------
// $Id: vector.cc,v 1.17 2003/01/08 17:58:17 wolf Exp $
// Version: $Name: Version-4-0-0 $
//
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
// This file is subject to QPL and may not be  distributed
// without copyright and license information. Please refer
// to the file deal.II/doc/license.html for the  text  and
// further information on this license.
//
//----------------------------  vector.cc  ---------------------------

//----------------------------  full_matrix.double.cc  ---------------------
//  $Id: full_matrix.double.cc,v 1.25 2003/05/06 23:16:30 wolf Exp $
//  Version: $Name: Version-4-0-0 $
//
//  Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//  This file is subject to QPL and may not be  distributed
//  without copyright and license information. Please refer
//  to the file deal.II/doc/license.html for the  text  and
//  further information on this license.
//
//----------------------------  full_matrix.double.cc  ---------------------

//modified from LAPACK++ file vd.cc and deal.II/lac/source/vector.cc
//  by John Trangenstein, August 2009

#include <float.h>
#include <math.h>
//#include "Arch.H"
//#include "Matrix.H"
//#include "MemoryDebugger.H"
//#include "MyInline.H"
#include "Vector.H"
//double double_zero_=0.;
//double double_one_=1.;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  int F77NAME(idamax)(const int &n,const double *sx,const int &incx);
  void F77NAME(daxpy)(const int &n,const double &a,double *x,
    const int &incx,double *y,const int &incy);
  void F77NAME(dcopy)(const int &n,const double *sx,const int &incx,
    double *sy,const  int &incy);
  double F77NAME(ddot)(const int &n,const double *sx,const int &incx,
                       const double *sy,const int &incy);
  void F77NAME(dgels)(const char *trans,const int &m,const int &n,
                      const int &nrhs,double *a,const int &lda,
                      double *b,const int &ldb,double *work,
                      const int &lwork,int &info);
  void F77NAME(dgemv)(char *trans,const int &m,const int &n,
                      const double &alpha,const double *a,
                      const int &lda,const double *x,const int &incx,
                      const double &beta,double *y,const int &incy);
  void F77NAME(dgesv)(const int &n,const int &nrhs,double *a,
                      const int &lda,const int *ipiv,double *b,
                      const int &ldb,int &info);
  double F77NAME(dnrm2)(const int &n,double *x,const int &incx);
  void F77NAME(dscal)(const int &n,const double &a,double *x,
                      const int &incx);
  void F77NAME(dswap)(const int &n,double *sx,const int &incx,
                      double *sy,const int &incy);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
LAVector& LAVector::operator+=(const LAVector &m) {
  F77NAME(daxpy)(sz,1.,m.data,1,data,1);
  return *this;
}

LAVector& LAVector::operator-=(const LAVector &m) {
  F77NAME(daxpy)(sz,-1.,m.data,1,data,1);
  return *this;
}

LAVector& LAVector::operator*=(double d) {
  F77NAME(dscal)(sz,d,data,1);
  return *this;
}

LAVector& LAVector::operator/=(double d) {
  assert(abs(d)>zero_);
  F77NAME(dscal)(sz,1./d,data,1);
  return *this;
}

double LAVector::l2Norm() const { return F77NAME(dnrm2)(sz,data,1); }
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void Vector::equ(double a,const Vector &v) {
  int inc=1;
  int n=d.size();
  F77NAME(dcopy)(n,v.d.addr(),inc,d.addr(),inc);
  if (a!=1.) F77NAME(dscal)(n,a,d.addr(),inc);
}

double Vector::operator*(const Vector &v) const {
  int inc=1;
  int n=d.size();
  return F77NAME(ddot)(n,d.addr(),inc,v.d.addr(),inc);
}

double Vector::linftyNorm() const {
  int inc=1;
  int n=d.size();
  int i=F77NAME(idamax)(n,d.addr(),inc);
  return abs(d[i-1]);
}

void Vector::add(double a,const Vector &y) {
  int inc=1;
  int n=d.size();
  F77NAME(daxpy)(n,a,y.d.addr(),inc,d.addr(),inc); 
}

void Vector::swapWith(Vector &v) {
  int inc=1;
  int n=d.size();
  F77NAME(dswap)(n,d.addr(),inc,v.d.addr(),inc); 
}
/*
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int LaEnvBlockSize(const char*,const LAMatrix&);
ostream& operator<<(ostream&,const LAMatrix&);

LAMatrix* productOf(const LAMatrix &A,const LAMatrix &X){
  int m=A.size(0),n=A.size(1),L=X.size(1);
  CHECK_SAME(n,X.size(0))
  LAMatrix *Y=OPERATOR_NEW LAMatrix(m,L);

  char trans= 'n';
  int incx=1,incy=1;
  double alpha=double_one_, beta=double_zero_;
  for (int j=0;j<L;j++) {
    F77NAME(dgemv)(&trans,m,n,alpha,A.addr(),m,X.addr(0,j),incx,
                   beta,&((*Y)(0,j)),incy);
  }
  return Y;
}

LAMatrix* factorNSolve(const LAMatrix &A,const LAMatrix &B) {
  int m=A.size(0), n=A.size(1), k=B.size(1);
  CHECK_SAME(m,B.size(0))
  LAMatrix Af(A);
  LAMatrix *X=OPERATOR_NEW LAMatrix(n,k);

  int info;
  if (m==n) {
    X->copy(B);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(dgesv)(m,k,Af.addr(),m,ipiv,X->addr(),n,info);
    CHECK_SAME(info,0)
    assert(info==0);
    delete [] ipiv;
  } else {
    char trans='n';
    int nb = LaEnvBlockSize("DGELSV", Af);
    if (m>n) {
      LAMatrix *Xtmp=OPERATOR_NEW LAMatrix(B);

      int W=n + nb * (m>k ? m : k);
      double *work=OPERATOR_NEW_BRACKET(double,W);
      F77NAME(dgels)(&trans,m,n,k,Af.addr(),m,Xtmp->addr(),m,work,W,
        info);
      CHECK_SAME(info,0)

      Xtmp->copyInto(*X);
      delete [] work;
      delete Xtmp;
    } else {
      B.copyInto(*X);
      int W=m + nb * (n>k ? n : k);
      double *work=OPERATOR_NEW_BRACKET(double,W);
      F77NAME(dgels)(&trans,m,n,k,Af.addr(),m,X->addr(),n,work,W,info);
      CHECK_SAME(info,0)
      delete [] work;
    }
  }
  return X;
}

void switchColumns(LAMatrix &G,int i,int j) {
  int m=G.size(0),n=G.size(1),inc0=1;
  CHECK_BOUNDS(i,0,n)
  CHECK_BOUNDS(j,0,n)
  F77NAME(dswap)(m,&(G(0,i)),inc0,&(G(0,j)),inc0);
}

void switchRows(LAMatrix &G,int i,int j) {
  int m=G.size(0),n=G.size(1),inc1=m;
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,m)
  F77NAME(dswap)(n,&(G(i,0)),inc1,&(G(j,0)),inc1);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Matrix::vmult(Vector &Mv,const Vector &v) const {
  int m=Mv.size();
  int n=v.size();
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),n);
  char trans='n';
  double alpha=1.;
  double beta=0.;
  int inc=1;
  F77NAME(dgemv)(&trans,m,n,alpha,M.addr(),m,v.addr(),inc,beta,
    Mv.addr(),inc);
}

void Matrix::transposeVmult(Vector &Mv,const Vector &v) const {
  int n=Mv.size();
  int m=v.size();
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),n);
  char trans='t';
  double alpha=1.;
  double beta=0.;
  int inc=1;
  F77NAME(dgemv)(&trans,m,n,alpha,M.addr(),m,v.addr(),inc,beta,
    Mv.addr(),inc);
}

void Matrix::vmultAdd(Vector &Mv,const Vector &v) const {
  int m=Mv.size();
  int n=v.size();
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),n);
  char trans='n';
  double alpha=1.;
  double beta=1.;
  int inc=1;
  F77NAME(dgemv)(&trans,m,n,alpha,M.addr(),m,v.addr(),inc,beta,
    Mv.addr(),inc);
}

void Matrix::transposeVmultAdd(Vector &Mv,const Vector &v) const {
  int n=Mv.size();
  int m=v.size();
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),n);
  char trans='t';
  double alpha=1.;
  double beta=1.;
  int inc=1;
  F77NAME(dgemv)(&trans,m,n,alpha,M.addr(),m,v.addr(),inc,beta,
    Mv.addr(),inc);
}

double Matrix::residual(Vector &r,const Vector &x,const Vector &b) const {
  int m=b.size();
  int n=x.size();
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),n);
  int inc=1;
  F77NAME(dcopy)(m,b.addr(),inc,r.addr(),inc);
  char trans='n';
  double alpha=-1.;
  double beta=1.;
  F77NAME(dgemv)(&trans,m,n,alpha,M.addr(),m,x.addr(),inc,beta,r.addr(),
    inc);
  return F77NAME(dnrm2)(m,r.addr(),inc);
}

void Matrix::preconditionJacobi(Vector &dst,const Vector &src,double om) 
const {
  int m=dst.size();
  int n=src.size();
  CHECK_SAME(m,n);
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),n);
  for (int i=0;i<m;i++) dst[i]=om*src[i]/M(i,i);
}

void Matrix::PSOR(Vector &dst,
const NumPtr<int> &permutation,
const NumPtr<int> &inverse_permutation,double om) const {
  int m=dst.size();
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),m);
  CHECK_SAME(permutation.getNumber(),m);
  CHECK_SAME(inverse_permutation.getNumber(),m);
  for (int i=0;i<m;i++) {
    int row=permutation[i];
    double s=dst[row];
    for (int j=0;j<m;j++) {
      if (inverse_permutation[j]<i) s-=M(row,j)*dst[j];
    }
    dst[row]=s*om/M(row,row);
  }
}

void Matrix::TPSOR(Vector &dst,const NumPtr<int> &permutation,
const NumPtr<int> &inverse_permutation,double om) const {
  int m=dst.size();
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),m);
  CHECK_SAME(permutation.getNumber(),m);
  CHECK_SAME(inverse_permutation.getNumber(),m);
  for (int i=m-1;i>=0;i--) {
    int row=permutation[i];
    double s=dst[row];
    for (int j=0;j<m;j++) {
      if (inverse_permutation[j]>i) s-=M(row,j)*dst[j];
    }
    dst[row]=s*om/M(row,row);
  } 
}

void Matrix::preconditionSOR(Vector &dst,const Vector &src,double om) 
const {
  int m=dst.size();
  CHECK_SAME(src.size(),m);
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),m);
  for (int row=0;row<m;row++) {
    double s=src[row];
    for (int col=0;col<row;col++) s-=M(row,col)*dst[col]; 
    dst[row]=s*om/M(row,row);
  }
}

void Matrix::preconditionTSOR(Vector &dst,const Vector &src,double om) 
const {
  int m=dst.size();
  CHECK_SAME(src.size(),m);
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),m);
  for (int row=m-1;row>=0;row--) {
    double s=src[row];
    for (int col=row+1;col<m;col++) s-=M(row,col)*dst[col]; 
    dst[row]=s*om/M(row,row);
  }
}

void Matrix::preconditionSSOR(Vector &dst,const Vector &src,double om) 
const {
  int m=dst.size();
  CHECK_SAME(src.size(),m);
  CHECK_SAME(size(0),m);
  CHECK_SAME(size(1),m);
  for (int i=0;i<m;i++) {
    double s=0.;
    for (int j=0;j<i;j++) s+=M(i,j)*dst[j];
    dst[i]-=s*om/M(i,i);
  }
  for (int i=m-1;i>=0;i--) {
    double s=0.;
    for (int j=i+1;j<m;j++) s+=M(i,j)*dst[j];
    dst[i]-=s*om/M(i,i);
  }
}
*/
