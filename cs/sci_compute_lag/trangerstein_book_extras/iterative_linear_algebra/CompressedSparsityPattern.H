#ifndef COMPRESSED_SPARSITY_PATTERN_H
#define COMPRESSED_SPARSITY_PATTERN_H "$Header:$"
//----------------------------  compressed_sparsity_pattern.h  -------------
//   $Id: compressed_sparsity_pattern.h,v 1.6 2003/01/08 17:58:13 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  compressed_sparsity_pattern.h  -------------
//
//modified fromdeal.II/lac/include/lac/compressed_sparsity_pattern.h
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

#include <set>
#include "NumPtr.H"

//not in DEAL, except as a template argument
class SparsityPatternBase {
  private:
    int rows;
    int cols;
  protected:
    explicit SparsityPatternBase(const SparsityPatternBase &spb) : 
    rows(spb.rows),cols(spb.cols) { 
      CHECK_SAME(spb.rows,0);
      CHECK_SAME(spb.cols,0);
    }
  public:
    SparsityPatternBase() : rows(0),cols(0) { }
    SparsityPatternBase(int m,int n) : rows(m),cols(n)
      { }
    SparsityPatternBase(int n) : rows(n),cols(n) { }
    virtual ~SparsityPatternBase() {}
    int nRows() const { return rows; }
    int nCols() const { return cols; }
    virtual bool empty() const { return (rows==0)&&(cols==0); }
    void reinit(int m,int n) { rows=m; cols=n; }
    virtual void add(int,int)=0;
    virtual bool exists(int,int) const=0;
    virtual void symmetrize()=0;
    virtual int columnNumber(int,int) const=0;
    virtual int bandwidth() const=0;
    virtual int nNonzeroElements() const=0;
    virtual int rowLength(int) const=0;
    virtual void printOn(ostream&) const;
};

//compare to DEAL CompressedSparsityPattern
class CompressedSparsityPattern : public SparsityPatternBase {
  private:
    NumPtr<set<int> > column_indices;

    CompressedSparsityPattern(const CompressedSparsityPattern &s) :
      SparsityPatternBase(s) { }
    CompressedSparsityPattern& operator=(
    const CompressedSparsityPattern &s) {
      CHECK_SAME(s.nRows(),0);
      CHECK_SAME(s.nCols(),0);
      CHECK_SAME(nRows(),0);
      CHECK_SAME(nCols(),0);
      return *this;
    }
  public:
    void reinit(int m,int n) {
      SparsityPatternBase::reinit(m,n);
      column_indices.cleanup();
      column_indices.allocate(nRows());
    }
    CompressedSparsityPattern() : SparsityPatternBase() { }
    CompressedSparsityPattern(int m,int n) : 
      SparsityPatternBase(m,n) { reinit (m,n); }
    CompressedSparsityPattern(int n) : SparsityPatternBase(n) {
      reinit(n,n);
    }
    ~CompressedSparsityPattern() {
      column_indices.cleanup();
    }
    int maxEntriesPerRow() const;
    void add(int i,int j) {
      CHECK_TEST(j<nCols());
      column_indices[i].insert(j);
    }
    bool exists(int i,int j) const {
      CHECK_TEST(j<nCols());
      return column_indices[i].find(j)!=column_indices[i].end();
    }
    void compress() {}

    void symmetrize();
    int columnNumber(int,int) const;
    int bandwidth() const;
    int nNonzeroElements() const;

    int rowLength(int row) const {
      return column_indices[row].size();
    }
    virtual void printOn(ostream&) const;
};

#endif
