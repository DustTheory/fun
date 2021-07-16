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
// Neither the Institutions (University of Tennessee, and Oak Ridge 
// National Laboratory) nor the Authors make any representations about 
// the suitability of this software for any purpose.  This software is 
// provided ``as is'' without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.

// Modified from LaPack++ gmc.C by John Trangenstein, 11/7/96
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

#include "Matrix.H"

extern "C" {
  int F77NAME(ilaenv)(const int &i,const char *n,const char &opts,
    const int &n1,const int &n2,const int &n3,const int &n4);
}

void LAMatrix::copyInto(LAMatrix &G) const {
  int m=(size(0)<G.size(0) ? size(0) : G.size(0));
  int n=(size(1)<G.size(1) ? size(1) : G.size(1));
  for (int j=0;j<n;j++) {
    memcpy(G.addr(0,j),addr(0,j),m*sizeof(double));
  }
}

LAMatrix* LAMatrix::transpose() const {
  int m=size(0),n=size(1);
  LAMatrix *X=new LAMatrix(n,m);
  for (int j=0;j<n;j++) {
    for (int i=0;i<m;i++) (*X)(j,i)=operator()(i,j);
  }
  return X;
}

int LaEnvBlockSize(const char *fname,
const LAMatrix& A) {
  char opts='U';
  int m=A.size(0),n=A.size(1),one=1,junk=-1;
  return F77NAME(ilaenv)(one,fname,opts,m,n,junk,junk);
}

void LAMatrix::printOn(ostream& s) const {
  s << "LAMatrix(" << size(0) << " x " << size(1) << ")\n" ;
  for (int i=0; i< size(0); i++) {
    for (int j=0; j< size(1); j++) s << operator()(i,j) << "  ";
    s << "\n";
  }
}

