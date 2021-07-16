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
// the suitability // of this software for any purpose.  This software 
// is provided ``as is'' without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.
  
// Modified from vc.C by John Trangenstein, 11/7/96
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

#include "Vector.H"
#include <limits>
#include <string.h>

#ifdef DEBUG
double LAVector::undefined_ = numeric_limits<double>::infinity();
#endif
double LAVector::zero_ = 0.;
double LAVector::one_ = 1.;

void LAVector::fillWith(double scalar) {
  for (int i=0; i<sz; i++) data[i] = scalar;
}

void LAVector::copy(const LAVector &m) {
  if (sz!=m.sz) resize(m.sz);
  memcpy(data,m.data,sz*sizeof(double));
}

void LAVector::printOn(ostream& s) const {
  s << "LAVector(" << sz << "):\n" ;
  for (int i=0;i<sz;i++) {
    s << "\tv[" << i << "] = " << data[i] << endl;
  }
}

void testVector(double scalar) {
  LAVector *a1=new LAVector;
  cout << "\nafter LAVector()" << endl;
  cout << "a1:" << endl;
  cout << *a1 << endl;

  cout << "enter n" << endl;
  int n;
  cin >> n;
  cout << "n = " << n << endl;
  LAVector *a2=new LAVector(n);;
  cout << "\nafter LAVector(int)" << endl;
  cout << "a2:" << endl;
  cout << *a2 << endl;

  a2->fillWith(scalar);
  cout << "after a2->fillWith(double)" << endl;
  cout << "a2:" << endl;
  cout << *a2 << endl;

  a1->copy(*a2);
  cout << "after a1->copy(a2)" << endl;
  cout << "a1:" << endl;
  cout << *a1 << endl;

  LAVector *a3=new LAVector(n,scalar);
  cout << "after a3=LAVector(n,scalar)" << endl;
  cout << "a3:" << endl;
  cout << *a3 << endl;
  delete a3;

  a3=new LAVector(*a2);
  cout << "after a3=LAVector(a2)" << endl;
  cout << "a3:" << endl;
  cout << *a3 << endl;

  LAVector *a4=new LAVector;
  a4->resize(n);
  cout << "after a4->resize(n)" << endl;
  cout << "a4:" << endl;
  cout << *a4 << endl;
  delete a4;

  delete a3;
  cout << "after delete a3" << endl;
  cout << "a2:" << endl;
  cout << *a2 << endl;

  a3=new LAVector();
  a3->resize(a2->size());
  cout << "after a3->resize(a2->size)" << endl;
  cout << "a3:" << endl;
  cout << *a3 << endl;
  delete a3;

  LAVector *temp;
  temp=(*a1)+(*a2);
  cout << "a1+a2:" << endl;
  cout << *temp << endl;
  delete temp;

  temp=(*a1)-(*a2);
  cout << "a1-a2:" << endl;
  cout << *temp << endl;
  delete temp;
  delete a2;

  temp=(*a1)*(3.);
  cout << "a1*3.:" << endl;
  cout << *temp << endl;
  delete temp;

  temp=(*a1)/(4.);
  cout << "a1/4.:" << endl;
  cout << *temp << endl;
  delete temp;
  delete a1;
}
