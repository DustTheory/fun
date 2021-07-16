#ifndef SPARSITY_PATTERN_H
#define SPARSITY_PATTERN_H "$Header:$"
//----------------------------  sparsity_pattern.h  ------------------------
//   $Id: sparsity_pattern.h,v 1.29 2003/01/08 17:58:15 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  sparsity_pattern.h  ------------------------
//
//modified from deal.II/lac/include/lac/sparsity_pattern.h
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

#include "CompressedSparsityPattern.H"
#include "Matrix.H"
#include "NumPtr.H"

/**
 * Structure representing the sparsity pattern of a sparse matrix.
 * 
 * The following picture will illustrate the relation between the
 * @p{SparsityPattern} an the @p{SparseMatrix}.
 *
 * @begin{verbatim}
 *  SparsityPattern:                               \
 *                                                 |
 *              _________________________          |
 *  rowstart   |0 | 4| 8|11|13|....                |
 *             |__|__|__|__|__|__________          | 
 *              |   \  \                           |
 *              |    \  \__                        | 
 *              |     \    \                       |
 *              |      \    \__                    |
 *             \ /      \      \                   |
 *              |        \      \__                |     
 *              |         \        \               |
 *              |          \        \__            |   
 *              |           \               \           |            
 *              0___________4___________8____       \ Position         
 *  colnums    |3 | 2| 9|17| 1| 4| 6| 8| 4|..       /         
 *             /__|/_|__|__|__|__|__|__|__|__      |         
 *            /   /                                |        
 *           / \______  _____/ \_____  _____/      |         
 *          /         \/             \/            |                 
 *         /   /  row = 0        row = 1           |    
 *        /   /                                    |
 *       /   /                                     |
 *      /   /                                      | 
 *     /   /___colnums[1]                          |  
 *    /                                            |
 *   /_________colnums[0]                          |
 *                                                 |                    
 *                                                /                    
 * @end{verbatim}
 *
 * @begin{verbatim}
 * For row = 0
 *   
 * it exists: (0| 3) = colnums[0]
 *            (0| 2) = colnums[1]
 *            (0| 9) = colnums[2]
 *            (0|17) = colnums[3]
 *
 * For row = 1
 *   
 * it exists: (1| 1) = colnums[4]
 *            (1| 4) = colnums[5]
 * ....
 *
 * @end{verbatim}
 *
 * @begin{verbatim}
 * SparseMatrix:                                  \
 *                                                 |
 *              _____________________________      |
 *  val        |  |  |  |  |  |  |  |  | 3|..       \ Value
 *             |__|__|__|__|__|__|__|__|__|__       /
 *                                                 |
 *                                                 |
 *                                                /
 * @end{verbatim}
 *
 * If you want to get the @p{3} you need to get its position in the
 * table above and its value by returning the value of the element on
 * which the pointer shows, using @p{*val}.  For example @p{val[8]=3}. Its
 * position is @p{colnums[8]} so @p{row=2}. In other words, if you want to get
 * the element @p{a_{24}} you know that @p{row=2}. To get the element, a
 * search of @p{4} form @p{colnums[rowstart[2]]} to @p{colnums[rowstart[3]]} is
 * needed. Then @p{a_{24}=val[number of the found element] = 3}.
 *
 *
 * @author Wolfgang Bangerth and others
 */

//compare to DEAL SparsityPattern
class SparsityPattern : public SparsityPatternBase {
  private:
    int max_row_length;
    NumPtr<int> rowstart;
    NumPtr<int> colnums;
    bool compressed;
    bool insert_diagonal_if_square;
  protected:
    void reinit(int,int,const NumPtr<int>&);
    void reinit(int m,int n,int max_per_row)
    {
      const NumPtr<int> row_lengths(m);
      row_lengths.initialize(max_per_row);
      reinit(m,n,row_lengths);
    }
  public:
//Define a value which is used to indicate that a certain value in the 
//@p{colnums} array is unused, i.e. does not represent a certain column
//number index.
//Indices with this invalid value are used to insert new entries to the
//sparsity pattern using the @p{add} member function, and are removed 
//when calling @p{compress}.
//You should not assume that the variable declared here has a certain 
//value. The initialization is given here only to enable the compiler to
//perform some optimizations, but the actual value of the variable may 
//change over time.
    static const int invalid_entry=INT_MAX; 

    explicit SparsityPattern(const SparsityPattern &s) : 
    SparsityPatternBase(s), /* cols(0), */ max_row_length(0),rowstart(0),
    colnums(0),insert_diagonal_if_square(true) {
      CHECK_SAME(s.rowstart.getNumber(),0);
      CHECK_SAME(s.colnums.getNumber(),0);
      reinit(0,0,0);
    }
    SparsityPattern(int m,int n,int max_per_row,bool idis=true) :
    SparsityPatternBase(m,n), /* cols(0), */ max_row_length(0),rowstart(0),
    colnums(0),insert_diagonal_if_square(idis) {
      reinit(m,n,max_per_row);
    }
    ~SparsityPattern () {
      rowstart.cleanup();
      colnums.cleanup();
    }
    SparsityPattern& operator=(const SparsityPattern &s) {
      CHECK_SAME(s.rowstart.getNumber(),0);
      CHECK_SAME(s.colnums.getNumber(),0);
      CHECK_SAME(rowstart.getNumber(),0);
      CHECK_SAME(colnums.getNumber(),0);
      return *this;
    }
    static const int* optimizedLowerBound(const int*,int,int);
    void compress();
    void copyFrom(const CompressedSparsityPattern&);
    void copyFrom(const SparsityPattern&);
    void copyFrom(const Matrix&);
    bool empty() const {
      if ((rowstart.getData()==0) || (nRows()==0) || (nCols()==0)) {
        CHECK_SAME(rowstart.getNumber(),0);
        CHECK_SAME(colnums.getNumber(),0);
        return true;
      }
      return false;
    }
    int maxEntriesPerRow() const;
    int diagonalEntry(int) const;
    int operator()(int,int) const;
    void add(int,int);
    bool exists(int,int) const;
    pair<int,int> matrixPosition(int) const;
    void symmetrize();
    int bandwidth() const;
    int nRows() const { 
      int rs=rowstart.getNumber();
      return (rs>0 ? rs-1 : 0); 
    }
    bool isCompressed() const { return compressed; }
    bool storesDiagonalFirst() const {
      return insert_diagonal_if_square && (nCols()==nRows());
    }
    int getRowstartIndex(int row) const { return rowstart[row]; }
    const NumPtr<int>& getColumnNumbers() const { return colnums; }
    int getColumnNumber(int col) const { return colnums[col]; }
    int& getColumnNumber(int col) { return colnums[col]; }
    int rowLength(int row) const { return rowstart[row+1]-rowstart[row]; }
    int columnNumber(int row,int index) const {
      return colnums[rowstart[row]+index];
    }
    int nNonzeroElements() const {
      CHECK_TEST((rowstart.getData()) && (colnums.getData()));
      CHECK_TEST(compressed);
      return rowstart[rowstart.getNumber()-1]-rowstart[0];
    }
    virtual void printOn(ostream&) const;
};

#endif
