// "$Header:$"
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
#include <ConditionNumber.H>
#include <MemoryDebugger.H>
#include <Vector.H>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  void F77NAME(dgttrf)(const int &n,double *dl,double *d,double *du,
    double *du2,int *ipiv,int &info);
  void F77NAME(dgtcon)(const char&,const int &n,double *dl,double *d,
    double *du,double *du2,int *ipiv,const REAL &anorm,REAL &rcond,
    REAL *work,int *iwork,int &info);
  REAL F77NAME(dasum)(const int &n,const REAL *dx,const int &incx);
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
REAL symmetricConditionNumber(const SolverMatrix &A,int mn) {
  int n=min(mn,A.size(0));

//symmetric Lanczos tridiagonalization process:
//see Golub & vanLoan p. 472
  Vector *qold=OPERATOR_NEW Vector(n);
  for (int i=0;i<n;i++) (*qold)[i]=drand48();
  qold->scale(1./qold->l2Norm());

  Vector Aq(n);
  REAL *D=OPERATOR_NEW_BRACKET(REAL,n);
  A.vmult(Aq,*qold);
  D[0]=(*qold)*Aq;

  Vector r(n);
  REAL *DL=OPERATOR_NEW_BRACKET(REAL,n-1);
  REAL *DU=OPERATOR_NEW_BRACKET(REAL,n-1);
  r=Aq;
  r.add(-D[0],*qold);

  Vector *q=OPERATOR_NEW Vector(n);
  int m=0;
  for (int k=1;k<n;k++) {
    DL[k-1]=DU[k-1]=r.l2Norm();
    if (DL[k-1]<=0.) break;

    m=k;
    q->equ(1./DL[k-1],r);
    A.vmult(Aq,*q);
    D[k]=(*q)*Aq;
    r=Aq;
    r.add(-D[k],*q);
    r.add(-DL[k-1],*qold);
//  in bits/stl_algobase.h:
    swap(q,qold);
  }

//tridiagonal matrix norm
  REAL anorm=0.;
  for (int j=0;j<m;j++) {
    REAL colmax=abs(D[j]);
    if (j<m-1) colmax=max(colmax,abs(DL[j]));
    if (j>0) colmax=max(colmax,abs(DU[j-1]));
    anorm=max(anorm,colmax);
  }

//tridiagonal factorization
  int info=0;
  REAL *DU2=OPERATOR_NEW_BRACKET(REAL,m-2);
  int *ipiv=OPERATOR_NEW_BRACKET(int,m);
  F77NAME(dgttrf)(m,DL,D,DU,DU2,ipiv,info);
  REAL rcond=numeric_limits<double>::infinity();
  REAL *work=OPERATOR_NEW_BRACKET(REAL,2*m);
  int *iwork=OPERATOR_NEW_BRACKET(int,m);
  info=0;
  F77NAME(dgtcon)('1',m,DL,D,DU,DU2,ipiv,anorm,rcond,work,iwork,info);

  delete iwork;
  delete work;
  delete ipiv;
  delete q;
  delete DU2;
  delete DU;
  delete DL;
  delete D;
  delete qold;
  return 1./rcond;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
REAL unSymmetricConditionNumber(const SolverMatrix &A,int mn) {
  int n=min(mn,A.size(0));

//symmetric Lanczos tridiagonalization process applied to A^T * A
//see Golub & vanLoan p. 472
  Vector *qold=OPERATOR_NEW Vector(n);
  for (int i=0;i<n;i++) (*qold)[i]=drand48();
  qold->scale(1./qold->l2Norm());

  Vector Aq(n);
  Vector AtAq(n);
  REAL *D=OPERATOR_NEW_BRACKET(REAL,n);
  A.vmult(Aq,*qold);
  D[0]=Aq*Aq;
  A.transposeVmult(AtAq,Aq);

  Vector r(n);
  REAL *DL=OPERATOR_NEW_BRACKET(REAL,n-1);
  REAL *DU=OPERATOR_NEW_BRACKET(REAL,n-1);
  r=AtAq;
  r.add(-D[0],*qold);
  DL[0]=DU[0]=r.l2Norm();

  Vector *q=OPERATOR_NEW Vector(n);
  int m=0;
  for (int k=1;k<n;k++) {
    m=k;

    q->equ(1./DL[k-1],r);
    A.vmult(Aq,*q);
    A.transposeVmult(AtAq,Aq);
    D[k]=Aq*Aq;
    r=AtAq;
    r.add(-D[k],*q);
    r.add(-DU[k-1],*qold);

    if (k<n-1) {
      DL[k]=DU[k]=r.l2Norm();
      if (DL[k]<=0.) break;
//    in bits/stl_algobase.h:
      swap(q,qold);
    }
  }

//tridiagonal matrix norm
  REAL anorm=0.;
  for (int j=0;j<m;j++) {
    REAL colmax=abs(D[j]);
    if (j<m-1) colmax=max(colmax,abs(DL[j]));
    if (j>0) colmax=max(colmax,abs(DU[j-1]));
    anorm=max(anorm,colmax);
  }

//tridiagonal factorization
  int info=0;
  REAL *DU2=OPERATOR_NEW_BRACKET(REAL,m-2);
  int *ipiv=OPERATOR_NEW_BRACKET(int,m);
  F77NAME(dgttrf)(m,DL,D,DU,DU2,ipiv,info);
  REAL rcond=numeric_limits<double>::infinity();
  REAL *work=OPERATOR_NEW_BRACKET(REAL,2*m);
  int *iwork=OPERATOR_NEW_BRACKET(int,m);
  info=0;
  F77NAME(dgtcon)('1',m,DL,D,DU,DU2,ipiv,anorm,rcond,work,iwork,info);

  delete iwork;
  delete work;
  delete ipiv;
  delete q;
  delete DU2;
  delete DU;
  delete DL;
  delete D;
  delete qold;
  return 1./sqrt(rcond);
}
