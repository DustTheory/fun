// "$Header:$"
//----------------------------  compressed_sparsity_pattern.cc  ------------
//   $Id: compressed_sparsity_pattern.cc,v 1.8 2003/01/08 17:58:16 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  compressed_sparsity_pattern.cc  ------------
//
//modified from deal.II/lac/source/compressed_sparsity_pattern.cc
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

#include "CompressedSparsityPattern.H"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "NumPtr.C"
INSTANTIATE_NUMPTR(int)

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparsityPatternBase::printOn(ostream &os) const {
  os << "SparsityPatternBase:"
     << "\n\trows = " << rows
     << "\n\tcols = " << cols 
     << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int CompressedSparsityPattern::maxEntriesPerRow() const {
  int m=0;
  for (int i=1;i<nRows();i++) {
    m=max(m,static_cast<int>(column_indices[i].size()));
  }
  return m;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CompressedSparsityPattern::symmetrize() {
  CHECK_SAME(nRows(),nCols());
  for (int row=0;row<nRows();row++) {
    for (set<int>::const_iterator i=
    column_indices[row].begin();i!=column_indices[row].end();i++) {
      if (row != *i) add(*i,row);
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int CompressedSparsityPattern::columnNumber(int row,int index) const {
  CHECK_TEST(index<column_indices[row].size());
  set<int>::const_iterator p=column_indices[row].begin();
  advance(p,index);
  return *p;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int CompressedSparsityPattern::bandwidth() const {
  int b=0;
  for (int row=0;row<nRows();row++) {
    for (set<int>::const_iterator i=
    column_indices[row].begin();i!=column_indices[row].end();i++) {
      if (abs(row-*i)>b) b=abs(row-*i);
    }
  }
  return b;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int CompressedSparsityPattern::nNonzeroElements() const {
  int n=0;
  for (int i=0;i<nRows();i++) n+=column_indices[i].size();
  return n;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void CompressedSparsityPattern::printOn(ostream &os) const {
  os << "CompressedSparsityPattern:"
     << "\n\tcolumn_indices.getNumber = " << column_indices.getNumber()
     << endl;
  for (int row=0;row<column_indices.getNumber();row++) {
    os << "\tcolumn_indices[" << row << "] = ";
    for (set<int>::const_iterator i=
    column_indices[row].begin();i!=column_indices[row].end();i++) {
      os << *i << " ";
    }
    os << endl;
  }
  SparsityPatternBase::printOn(os);
}
