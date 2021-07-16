// "$Header:$"
//------------------------  sparsity_pattern.cc  ---------------------------
//   $Id: sparsity_pattern.cc,v 1.39 2003/01/08 17:58:17 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  sparsity_pattern.cc  -----------------------
//
//modified from deal.II/lac/source/sparsity_pattern.cc
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
#include "SparsityPattern.H"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

//#include "CompressedSparsityPattern.H"
#include "Matrix.H"
#include "MemoryDebugger.H"

#ifndef PRINT_BOOLEAN
#define PRINT_BOOLEAN(b) ((b)?"TRUE":"FALSE")
#endif

const int SparsityPattern::invalid_entry;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
const int* SparsityPattern::optimizedLowerBound(const int *first,
int len,int val) {
  if (len==0) return first;
  while (true) {
    if (len<8) {
      for (int loop=len;loop>=1;loop--) {
        if (*first>=val) return first;
        first++;
      }
      return first;
    }
    int half=len>>1;
    const int *middle=first+half;
    if (*middle<val) {
      first=middle+1;
      len-=half+1;
    } else len=half;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//if m==n, will assume that matrix has a nonzero diagonal; store first
void SparsityPattern::reinit(int m,int n,const NumPtr<int> &row_lengths) {
  CHECK_SAME(row_lengths.getNumber(),m);
  SparsityPatternBase::reinit(m,n);
  if (m==0 || n==0) {
    SparsityPatternBase::reinit(0,0);
    rowstart.cleanup();
    colnums.cleanup();
    max_row_length=0;
    compressed=false;
    return;
  }
  bool matrix_is_square=(m==n && insert_diagonal_if_square);
  int vec_len=0;
  for (int i=0;i<m;i++) {
    vec_len+=
      min((matrix_is_square?max(row_lengths[i],1):row_lengths[i]),n);
  }
  if (vec_len==0) vec_len=1;
  int rlm=0;
  for (int i=0;i<row_lengths.getNumber();i++) {
    rlm=max(rlm,row_lengths[i]);
  }
  max_row_length=min(rlm,n);
  if (matrix_is_square&&(max_row_length==0)&&(m!=0)) max_row_length=1;
  int rs=rowstart.getNumber();
  if (m+1!=rs) {
    rowstart.cleanup();
    rowstart.allocate(m+1);
  }
  if (vec_len>colnums.getNumber()) {
    colnums.cleanup();
    colnums.allocate(vec_len);
  }
  rowstart[0]=0;
  for (int i=1;i<=m;i++) {
    rowstart[i]=rowstart[i-1]
      +(matrix_is_square ? max(min(row_lengths[i-1],n),1) 
                         : min(row_lengths[i-1],n));
  }
  CHECK_TEST((rowstart[m]==vec_len) ||
             ((vec_len==1) && (rowstart[m]==0)));
  colnums.initialize(invalid_entry);
  if (matrix_is_square) {
    for (int i=0;i<m;i++) colnums[rowstart[i]] = i;
  }
  compressed=false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparsityPattern::compress() {
  CHECK_TEST((rowstart.getData()!=0) && (colnums.getData()!=0));
  if (compressed) return;
  int next_free_entry=0,next_row_start=0;
  int nonzero_elements=0;
  for (int i=rowstart[0];i<rowstart[nRows()];i++) {
    nonzero_elements+=(colnums[i]!=invalid_entry);
  }
  NumPtr<int> new_colnums(nonzero_elements); 
  vector<int> tmp_entries(max_row_length);
  bool is_square=(nRows()==nCols() && insert_diagonal_if_square);
  for (int line=0;line<nRows();line++) {
    int row_length=0;
    for (int j=rowstart[line];j<rowstart[line+1]; 
    j++,row_length++) {
      if (colnums[j] != invalid_entry) {
        tmp_entries[row_length] = colnums[j];
      } else break;
    }
    if (row_length>1) {
      sort( is_square ? tmp_entries.begin()+1 
        : tmp_entries.begin(), tmp_entries.begin()+row_length);
    }
    for (int j=0;j<row_length;j++) {
      new_colnums[next_free_entry++]=tmp_entries[j];
    }
    rowstart[line]=next_row_start;
    next_row_start=next_free_entry;
#ifdef DEBUG
    CHECK_TEST(!is_square || (new_colnums[rowstart[line]] == line));
    if (rowstart[line]<next_row_start) {
      int col=next_row_start;
      int nc=new_colnums[rowstart[line]];
      for (int c=rowstart[line]+1;c<next_row_start;c++) {
        CHECK_TEST(new_colnums[c]!=nc);
        if (c+1<next_row_start) {
          CHECK_TEST(new_colnums[c]!=new_colnums[c+1]);
        }
      }
      CHECK_TEST(col==next_row_start);
    }
#endif
  }
  CHECK_SAME(next_free_entry,nonzero_elements);
  rowstart[nRows()]=next_row_start;
  colnums.copy(new_colnums);
  compressed=true;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparsityPattern::copyFrom(const CompressedSparsityPattern &csp) {
  int m=csp.nRows(),n=csp.nCols();
  bool is_square=(m==n);
  insert_diagonal_if_square=true;
  NumPtr<int> row_lengths(m);
  for (int i=0;i<m;i++) {
    row_lengths[i]=csp.rowLength(i);
  }
  reinit(m,n,row_lengths);
  for (int row=0;row<m;row++) {
    int *row_cols=&colnums[rowstart[row]]+is_square;
    for (int j=0;j<csp.rowLength(row);j++) {
      int col=csp.columnNumber(row,j);
      CHECK_TEST(col<n);
      if ((col!=row) || !is_square) *row_cols++=col;
    }
  }
  compress();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparsityPattern::copyFrom(const SparsityPattern &sp) {
  int m=sp.nRows(),n=sp.nCols();
  SparsityPatternBase::reinit(m,n);
  max_row_length=sp.max_row_length;
  compressed=sp.compressed;
  insert_diagonal_if_square=sp.insert_diagonal_if_square;
  rowstart.cleanup();
  rowstart.allocate(sp.rowstart.getNumber());
  rowstart.copyFrom(sp.rowstart);
  colnums.cleanup();
  colnums.allocate(sp.colnums.getNumber());
  colnums.copyFrom(sp.colnums);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparsityPattern::copyFrom(const Matrix &matrix) {
  NumPtr<int> entries_per_row(matrix.size(0));
  entries_per_row.initialize(0);
  for (int row=0;row<matrix.size(0);row++) {
    for (int col=0; col<matrix.size(1);col++) {
      if (matrix(row,col)!=0.) ++entries_per_row[row];
    }
  }
  insert_diagonal_if_square=true;
  reinit(matrix.size(0),matrix.size(1),entries_per_row);
  for (int row=0;row<matrix.size(0);row++) {
    for (int col=0;col<matrix.size(1);col++) {
      if (matrix(row,col)!=0.) add(row,col);
    }
  }
  compress();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int SparsityPattern::maxEntriesPerRow() const {
  if (!compressed) return max_row_length;
  int m=0;
  for (int i=1;i<nRows();i++) {
    m=max(m,rowstart[i]-rowstart[i-1]);
  }
  return m;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int SparsityPattern::diagonalEntry(int row) const {
  CHECK_SAME(nRows(),nCols());
  if (insert_diagonal_if_square) return rowstart[row];
  const int *sorted_region_start=&colnums[rowstart[row]];
  int len=rowstart[row+1]-rowstart[row];
  const int*const p=optimizedLowerBound(sorted_region_start,len,row);
  if (p<colnums.getData()+rowstart[row+1] && *p==row) {
    return p-colnums.getData();
  } else return invalid_entry;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int SparsityPattern::operator()(int i,int j) const {
  CHECK_TEST((rowstart.getData()!=0) && (colnums.getData()!=0));
  CHECK_TEST(j<nCols());
  CHECK_TEST(compressed);
  bool matrix_is_square=(nRows()==nCols() && insert_diagonal_if_square);
  if (rowstart[i]==rowstart[i+1]) return invalid_entry;
  if ((i==j)&&matrix_is_square) return rowstart[i];
  const int *sorted_region_start=0;
  int len=0;
  if (matrix_is_square) {
    sorted_region_start=&colnums[rowstart[i]+1];
    len=rowstart[i+1]-rowstart[i]-1;
  } else {
    sorted_region_start=&colnums[rowstart[i]];
    len=rowstart[i+1]-rowstart[i];
  }
  const int*const p=
    optimizedLowerBound(sorted_region_start,len,j);
  if (p<colnums.getData()+rowstart[i+1] && *p==j) {
    return p-colnums.getData();
  } else return invalid_entry;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparsityPattern::add(int i,int j) {
  CHECK_TEST((rowstart.getData()!=0) && (colnums.getData()!=0));
  CHECK_TEST(j<nCols());
  CHECK_TEST(!compressed);
  for (int k=rowstart[i];k<rowstart[i+1];k++) {
    if (colnums[k]==j) return;
    if (colnums[k]==invalid_entry) {
      colnums[k]=j;
      return;
    }
  }
  OBSOLETE("not enough space in SparsityPattern::add");
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool SparsityPattern::exists(int i,int j) const {
  CHECK_TEST((rowstart.getData()!=0) && (colnums.getData()!=0));
  CHECK_TEST(j<nCols());
  for (int k=rowstart[i];k<rowstart[i+1];k++) {
    if (colnums[k]==j) return true;
  }
  return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
pair<int,int> SparsityPattern::matrixPosition(
int global_index) const {
  CHECK_TEST(compressed);
  CHECK_TEST(global_index<nNonzeroElements());
  int row=
    upper_bound(&rowstart[0],&rowstart[nRows()],global_index)
    -&rowstart[0]-1;
  int col=colnums[global_index];
  return make_pair(row,col);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparsityPattern::symmetrize () {
  CHECK_TEST((rowstart.getData()!=0) && (colnums.getData()!=0));
  CHECK_TEST(!compressed);
  CHECK_SAME(nRows(),nCols());
  for (int row=0;row<nRows();row++) {
    for (int k=rowstart[row];k<rowstart[row+1];k++) {
      if (colnums[k]==invalid_entry) break;
      if (colnums[k]!=row) add(colnums[k],row);
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int SparsityPattern::bandwidth () const {
  CHECK_TEST((rowstart.getData()!=0) && (colnums.getData()!=0));
  int b=0;
  for (int i=0;i<nRows();i++) {
    for (int j=rowstart[i];j<rowstart[i+1];j++) {
      if (colnums[j]!=invalid_entry) {
        int ij=abs(static_cast<signed int>(i-colnums[j]));
        if (ij>b) b=ij;
      } else break;
    }
  }
  return b;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparsityPattern::printOn(ostream &os) const {
  os << "SparsityPattern:" << endl;
  os << "\tmax_row_length = " << max_row_length
     << "\n\trowstart.getNumber = " << rowstart.getNumber()
     << "\n\tcolnums.getNumber = " << colnums.getNumber()
     << "\n\tcompressed = " << PRINT_BOOLEAN(compressed)
     << endl;
  for (int i=0;i+1<rowstart.getNumber();i++) {
    cout << "\tcolumns in row " << i << " = ";
    for (int j=rowstart[i];j<rowstart[i+1];j++) {
      if (colnums[j]!=invalid_entry) cout << " " << colnums[j];
    }
    cout << endl;
  }
  SparsityPatternBase::printOn(os);
}

template int const* lower_bound<int const*, int>(
  int const*, int const*, int const&);
template int const* upper_bound<int const*, int>(
  int const*, int const*, int const&);
template int* adjacent_find<int*>(int*, int*);

template void __introsort_loop<
  __gnu_cxx::__normal_iterator<int*, vector<int, allocator<int> > >,
  int>
(__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,
__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,int);

template void __final_insertion_sort<
  __gnu_cxx::__normal_iterator<int*,vector<int,allocator<int> > > >
(__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,
__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >);

template void partial_sort<
  __gnu_cxx::__normal_iterator<int*,vector<int,allocator<int> > > >
(__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,
__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,
__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >);

template __gnu_cxx::__normal_iterator<
  int*, 
  vector<int, allocator<int> > >
__unguarded_partition<
  __gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,
  int>
(__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,
__gnu_cxx::__normal_iterator<int*, 
                             vector<int, allocator<int> > >,int);

template void __insertion_sort<
  __gnu_cxx::__normal_iterator<int*,vector<int,allocator<int> > > >
(__gnu_cxx::__normal_iterator<int*,vector<int,allocator<int> > >,
__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >);

template void make_heap<
  __gnu_cxx::__normal_iterator<int*,vector<int,allocator<int> > > >
(__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,
__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >);

template void sort_heap<
  __gnu_cxx::__normal_iterator<int*,vector<int,allocator<int> > > >
(__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,
__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >);

template void __unguarded_linear_insert<
  __gnu_cxx::__normal_iterator<int*,vector<int,allocator<int> > >,int>
(__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,int);

template void __adjust_heap<
  __gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,int,int>
(__gnu_cxx::__normal_iterator<int*,vector<int, allocator<int> > >,int,int,
  int);

template void __push_heap<
  __gnu_cxx::__normal_iterator<int*, 
                               vector<int,allocator<int> > >,
  int, 
  int>
(__gnu_cxx::__normal_iterator<int*, 
                              vector<int, allocator<int> > >,
int, int, int);

template int* fill_n<int*, int, int>(int*, int, int const&);

template void __heap_select<
  __gnu_cxx::__normal_iterator<int*, vector<int, allocator<int> > > 
>(
  __gnu_cxx::__normal_iterator<int*, vector<int, allocator<int> > >, 
  __gnu_cxx::__normal_iterator<int*, vector<int, allocator<int> > >, 
  __gnu_cxx::__normal_iterator<int*, vector<int, allocator<int> > >
);
