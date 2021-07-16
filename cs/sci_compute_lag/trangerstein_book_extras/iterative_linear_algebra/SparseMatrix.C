// "$Header:$"
//----------------------------  sparse_matrix.templates.h  -----------------
// $Id: sparse_matrix.templates.h,v 1.72.2.1 2003/06/05 12:11:22 guido Exp $
// Version: $Name: Version-4-0-0 $
//
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
// This file is subject to QPL and may not be  distributed
// without copyright and license information. Please refer
// to the file deal.II/doc/license.html for the  text  and
// further information on this license.
//
//----------------------------  sparse_matrix.templates.h  -----------------
//
//modified from deal.II/lac/include/lac/sparse_matrix.templates.h
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

#include "SparseMatrix.H"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <numeric>
#include <ostream>
#include "MemoryDebugger.H"
#include "Tracer.H"

#include "NumPtr.C"
INSTANTIATE_NUMPTR(double)

double SparseMatrix::bogus=numeric_limits<double>::infinity();
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::reinit() {
  CHECK_POINTER(cols);
  CHECK_TEST(cols->isCompressed() || cols->empty());
  if (cols->empty()) {
    val.cleanup();
    return;
  }
  int N=cols->nNonzeroElements();
  int max_len=val.getNumber();
  if (N>max_len || max_len==0) {
    val.cleanup();
    val.allocate(N);
  }
  val.initialize(0.);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int SparseMatrix::nActuallyNonzeroElements() const {
  CHECK_POINTER(cols);
  int nnz=0;
  for (int i=0;i<nNonzeroElements();i++) {
    nnz+=(abs(val[i])!=0.);
  }
  return nnz;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::symmetrize() {
  CHECK_POINTER(cols);
  CHECK_SAME(cols->nRows(),cols->nCols());
  int n_rows=cols->nRows();
  for (int row=0;row<n_rows;row++) {
    int first_offdiag_col=cols->getRowstartIndex(row)+1;
    NumPtr<double> val_ptr(val,first_offdiag_col);
    NumPtr<int> 
      colnum_ptr(cols->getColumnNumbers(),first_offdiag_col);      
    for (int index=0;index<cols->rowLength(row);index++) {
      double mean_value=
        (val_ptr[index]+val[(*cols)(colnum_ptr[index],row)])*0.5;
      val_ptr[index]=mean_value;
      set(colnum_ptr[index],row,mean_value);
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::copyFrom(const Matrix &matrix) {
  reinit();
  for (int row=0;row<matrix.size(0);row++) {
    for (int col=0;col<matrix.size(1);col++) {
      double mij=matrix(row,col);
      if (abs(mij)!=0.) set(row,col,mij);
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::addScaled(double factor,const SparseMatrix &matrix) {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  if (cols==matrix.cols) {
    for (int i=0;i<val.getNumber();i++) {
      val[i]+=factor*matrix.val[i];
    }
  } else {
//  the following will fail if matrix.cols not contained in this->cols
    for (int i=0;i<matrix.cols->nRows();i++) {
      for (int j=matrix.cols->getRowstartIndex(i);
      j<matrix.cols->getRowstartIndex(i+1);j++) {
        add(i,matrix.cols->getColumnNumber(j),factor*matrix.val[j]);
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::vmult(Vector &dst,const Vector &src) const {
//TRACER_CALL(t,"SparseMatrix::vmult");
//cols->printOn(cout);
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  CHECK_SAME(n_rows,dst.size());
  CHECK_SAME(cols->nCols(),src.size());
  for (int row=0;row<n_rows;row++) {
    double s=0.;
    int row_start=cols->getRowstartIndex(row);
    NumPtr<double> val_ptr(val,row_start);
    NumPtr<int> colnum_ptr(cols->getColumnNumbers(),row_start);
    for (int index=0;index<cols->rowLength(row);index++) {
      s+=val_ptr[index]*src[colnum_ptr[index]];
    }
    dst[row]=s;
  };
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::transposeVmult(Vector &dst,const Vector &src) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  CHECK_SAME(n_rows,src.size());
  CHECK_SAME(cols->nCols(),dst.size());
  dst=0.;
  for (int i=0;i<n_rows;i++) {
    for (int j=cols->getRowstartIndex(i);
    j<cols->getRowstartIndex(i+1);j++) {
      int p=cols->getColumnNumber(j);
      dst[p]+=val[j]*src[i];
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::vmultAdd(Vector& dst, const Vector& src) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  CHECK_SAME(n_rows,dst.size());
  CHECK_SAME(cols->nCols(),src.size());
  for (int row=0;row<n_rows;row++) {
    double s=0.;
    int row_start=cols->getRowstartIndex(row);
    NumPtr<double> val_ptr(val,row_start);
    NumPtr<int> colnum_ptr(cols->getColumnNumbers(),row_start);
    for (int index=0;index<cols->rowLength(row);index++) {
      s+=val_ptr[index]*src[colnum_ptr[index]];
    }
    dst[row]+=s;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::transposeVmultAdd(Vector& dst,const Vector& src) 
const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  CHECK_SAME(n_rows,src.size());
  CHECK_SAME(cols->nCols(),dst.size());
  for (int i=0;i<n_rows;i++) {
    for (int j=cols->getRowstartIndex(i);
    j<cols->getRowstartIndex(i+1);j++) {
      int p=cols->getColumnNumber(j);
      dst[p]+=val[j]*src[i];
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double SparseMatrix::matrixNormSquare(const Vector &v) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  CHECK_SAME(n_rows,v.size());
  CHECK_SAME(cols->nCols(),v.size());
  double sum=0.;
  for (int row=0;row<n_rows;row++) {
    double s=0.;
    int row_start=cols->getRowstartIndex(row);
    NumPtr<double> val_ptr(val,row_start);
    NumPtr<int> colnum_ptr(cols->getColumnNumbers(),row_start);
    for (int index=0;index<cols->rowLength(row);index++) {
      s+=val_ptr[index]*v[colnum_ptr[index]];
    }
    sum+=s*v[row];
  }
  return sum;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double SparseMatrix::matrixScalarProduct(const Vector &u,const Vector &v) 
const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  CHECK_SAME(n_rows,u.size());
  CHECK_SAME(cols->nCols(),v.size());
  double sum=0.;
  for (int row=0;row<n_rows;row++) {
    double s=0.;
    int row_start=cols->getRowstartIndex(row);
    NumPtr<double> val_ptr(val,row_start);
    NumPtr<int> colnum_ptr(cols->getColumnNumbers(),row_start);
    for (int index=0;index<cols->rowLength(row);index++) {
      s+=val_ptr[index]*v[colnum_ptr[index]];
    }
    sum+=s*u[row];
  }
  return sum;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double SparseMatrix::l1Norm() const {
  CHECK_POINTER(cols);
  Vector column_sums(cols->nCols());
  column_sums=0.;
  int n_rows=cols->nRows();
  for (int row=0;row<n_rows;row++) {
    int row_start=cols->getRowstartIndex(row);
    NumPtr<double> val_ptr(val,row_start);
    NumPtr<int> colnum_ptr(cols->getColumnNumbers(),row_start);
    for (int index=0;index<cols->rowLength(row);index++) {
      column_sums[colnum_ptr[index]]+=abs(val_ptr[index]);
    }
  }
  return column_sums.linftyNorm();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double SparseMatrix::linftyNorm() const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  double max=0;
  int n_rows=cols->nRows();
  for (int row=0;row<n_rows;row++) {
    double sum=0;
    NumPtr<double> val_ptr(val,cols->getRowstartIndex(row));
    for (int index=0;index<cols->rowLength(row);index++) {
      sum+=abs(val_ptr[index]);
    }
    if (sum>max) max=sum;
  }
  return max;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double SparseMatrix::residual(Vector &dst,const Vector &u,const Vector &b) 
const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  CHECK_SAME(n_rows,dst.size());
  CHECK_SAME(n_rows,b.size());
  CHECK_SAME(cols->nCols(),u.size());
  double norm=0.;   
  for (int i=0;i<n_rows;i++) {
    double s=b[i];
    for (int j=cols->getRowstartIndex(i);
    j<cols->getRowstartIndex(i+1);j++) {
      int p=cols->getColumnNumber(j);
      s-=val[j]*u[p];
    }
    dst[i]=s;
    norm+=dst[i]*dst[i];
  }
  return sqrt(norm);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::preconditionJacobi(Vector &dst,const Vector &src,
double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(dst.size(),n_cols);
  CHECK_SAME(src.size(),n_cols);
  for (int i=0;i<n_cols;i++) {
    dst[i]=om*src[i]/val[cols->getRowstartIndex(i)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::preconditionSSOR(Vector &dst,const Vector &src,
double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(dst.size(),n_rows);
  CHECK_SAME(src.size(),n_cols);
  NumPtr<int> colnum_ptr(cols->getColumnNumbers());
  for (int row=0;row<n_cols;row++) {
    dst[row]=src[row];
    int row_start=cols->getRowstartIndex(row);
    int first_right_of_diagonal_index = 
      (SparsityPattern::optimizedLowerBound(&colnum_ptr[row_start+1],
      cols->getRowstartIndex(row+1)-row_start-1,row)-&colnum_ptr[0]);
    for (int j=row_start+1;j<first_right_of_diagonal_index;j++)
    {
      dst[row]-=om*val[j]*dst[colnum_ptr[j]];
    }
    dst[row]/=val[row_start];
  }
  for (int row=0;row<n_cols;row++) {
    dst[row]*=(2.-om)*val[cols->getRowstartIndex(row)];
  }
  for (int row=n_cols-1; row>=0; --row) {
    int next_row_start=cols->getRowstartIndex(row+1);
    const int first_right_of_diagonal_index = 
      (SparsityPattern::optimizedLowerBound(
      &colnum_ptr[cols->getRowstartIndex(row)+1],
      next_row_start-cols->getRowstartIndex(row)-1,
      static_cast<int>(row))
      -&colnum_ptr[0]);
    for (int j=first_right_of_diagonal_index;j<next_row_start;
    j++) {
      if (colnum_ptr[j]>static_cast<int>(row)) {
        dst[row]-=om*val[j]*dst[colnum_ptr[j]];
      }
    }
    dst[row]/=val[cols->getRowstartIndex(row)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::SOR(Vector &dst,double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(n_rows,dst.size());
  for (int row=0;row<n_rows;row++) {
    double s=dst[row];
    for (int j=cols->getRowstartIndex(row);
    j<cols->getRowstartIndex(row+1);j++) {
      int col=cols->getColumnNumber(j);
      if (col<row) s-=val[j]*dst[col];
    }
    dst[row]=s*om/val[cols->getRowstartIndex(row)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::TSOR(Vector &dst,double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(n_rows,dst.size());
  for (int row=n_rows;row!=0;) {
    --row;
    double s=dst[row];
    for (int j=cols->getRowstartIndex(row);
    j<cols->getRowstartIndex(row+1);j++) {
      if (cols->getColumnNumber(j)>row) {
        s-=val[j]*dst[cols->getColumnNumber(j)];
      }
    }
    dst[row]=s*om/val[cols->getRowstartIndex(row)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::PSOR(Vector &dst,const NumPtr<int> &permutation,
const NumPtr<int> &inverse_permutation,double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(n_rows,dst.size());
  CHECK_SAME(n_rows,permutation.getNumber());
  CHECK_SAME(n_rows,inverse_permutation.getNumber());
  for (int urow=0;urow<n_rows;urow++) {
    int row=permutation[urow];
    double s=dst[row];
    for (int j=cols->getRowstartIndex(row);
    j<cols->getRowstartIndex(row+1);j++) {
      int col=cols->getColumnNumber(j);
      if (inverse_permutation[col]<urow) s-=val[j]*dst[col];
    }
    dst[row]=s*om/val[cols->getRowstartIndex(row)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::TPSOR(Vector &dst,const NumPtr<int> &permutation,
const NumPtr<int> &inverse_permutation,double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(n_rows,dst.size());
  CHECK_SAME(n_rows,permutation.getNumber());
  CHECK_SAME(n_rows,inverse_permutation.getNumber());
  for (int urow=n_rows;urow!=0;) {
    --urow;
    int row=permutation[urow];
    double s=dst[row];
    for (int j=cols->getRowstartIndex(row);
    j<cols->getRowstartIndex(row+1);j++) {
      int col=cols->getColumnNumber(j);
      if (inverse_permutation[col]>urow) s-=val[j]*dst[col];
    }
    dst[row]=s*om/val[cols->getRowstartIndex(row)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::SORStep(Vector &v,const Vector &b,double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(n_rows,v.size());
  CHECK_SAME(n_rows,b.size());
  for (int row=0;row<n_rows;row++) {
    double s=b[row];
    for (int j=cols->getRowstartIndex(row);
    j<cols->getRowstartIndex(row+1);j++) {
      s-=val[j]*v[cols->getColumnNumber(j)];
    }
    v[row]+=s*om/val[cols->getRowstartIndex(row)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::TSORStep(Vector &v,const Vector &b,double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(n_rows,v.size());
  CHECK_SAME(n_rows,b.size());
  for (int row=n_rows-1;row>=0;row--) {
    double s=b[row];
    for (int j=cols->getRowstartIndex(row);
    j<cols->getRowstartIndex(row+1);j++) {
      s-=val[j]*v[cols->getColumnNumber(j)];
    }
    v[row]+=s*om/val[cols->getRowstartIndex(row)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::SSOR(Vector &dst,double om) const {
  CHECK_POINTER(cols);
  CHECK_POINTER(val.getData());
  int n_rows=cols->nRows();
  int n_cols=cols->nCols();
  CHECK_SAME(n_rows,n_cols);
  CHECK_SAME(n_rows,dst.size());
  for (int i=0; i<n_rows; i++) {
    double s=0.;
    for (int j=cols->getRowstartIndex(i);
    j<cols->getRowstartIndex(i+1);j++) {
      int p=cols->getColumnNumber(j);
      if (p!=SparsityPattern::invalid_entry) {
        if (i>j) s+=val[j]*dst[p];
      }
    }
    dst[i]-=s*om;
    dst[i]/=val[cols->getRowstartIndex(i)];
  }
  for (int i=n_rows-1;i>=0;i--) {
    double s=0.;
    for (int j=cols->getRowstartIndex(i);
    j<cols->getRowstartIndex(i+1);j++) {
      int p=cols->getColumnNumber(j);
      if (p!=SparsityPattern::invalid_entry) {
        if (static_cast<int>(i)<j) s+=val[j]*dst[p];
      }
    }
    dst[i]-=s*om/val[cols->getRowstartIndex(i)];
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SparseMatrix::printOn(ostream &os) const {
  CHECK_POINTER(cols);
  CHECK_TEST(val.getNumber()>0);
  os << "\tSparseMatrix:" << endl;
  int n_rows=cols->nRows();
  os << "\tcols->nRows() = " << n_rows << endl;
  for (int i=0;i<=n_rows;i++) {
    os << "\tgetRowstartIndex(" << i << ") = " 
       << cols->getRowstartIndex(i) << endl; 
  }
  os << "\tval.getNumber = " << val.getNumber() << endl;
  for (int i=0;i<n_rows;i++) {
    for (int j=cols->getRowstartIndex(i);
    j<cols->getRowstartIndex(i+1);j++) {
      os << "(" << i << "," << cols->getColumnNumber(j) << ") " 
         << val[j] << endl;
    }
  }
}
