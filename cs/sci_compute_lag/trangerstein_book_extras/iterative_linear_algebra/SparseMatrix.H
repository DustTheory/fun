#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H "$Header:$"
//----------------------------  sparse_matrix.h  ---------------------------
//   $Id: sparse_matrix.h,v 1.97 2003/06/03 09:16:36 guido Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  sparse_matrix.h  ---------------------------
//
//modified from deal.II/lac/include/lac/sparse_matrix.h
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

#include "Solver.H"
#include "SolverArrays.H"
#include "SparsityPattern.H"

//compare to DEAL SparseMatrix
class SparseMatrix : public SolverMatrix {
  private:
    const SparsityPattern *cols;
    NumPtr<double> val;

    static double bogus;
  public:
    virtual void reinit();
    virtual void reinit(const SparsityPattern *sparsity) {
      cols=sparsity; reinit();
    }
    SparseMatrix() : cols(0),val(0) {}
    SparseMatrix(const SparseMatrix &m) : cols(0),val(0) {
      CHECK_SAME(m.cols,0);
      CHECK_SAME(m.val.getData(),0);
    }
    explicit SparseMatrix(const SparsityPattern *c) : cols(c),val(0) { 
      reinit(); 
    }
    virtual ~SparseMatrix() { 
      val.cleanup(); 
    }
    SparseMatrix& operator=(const SparseMatrix &m) {
      CHECK_SAME(m.cols,0);
      CHECK_SAME(m.val.getData(),0);
      return *this;
    }
    SparseMatrix* clone() const {
      CHECK_POINTER(cols);
      return OPERATOR_NEW SparseMatrix(cols);
    }
    SparseMatrix* cloneEmpty() const {
      return OPERATOR_NEW SparseMatrix;
    }
    virtual void clear() { cols=0; val.cleanup(); }
    bool empty() const {
      if (cols==0) return true;
      return cols->empty();
    }
    int size(int i) const {
      CHECK_BOUNDS(i,0,2);
      CHECK_POINTER(cols);
      return (i==0 ? cols->nRows() : cols->nCols());
    }
    int nNonzeroElements() const {
      CHECK_POINTER(cols);
      return cols->nNonzeroElements();
    }
    void set(int i,int j,double value) {
      CHECK_POINTER(cols);
      int index=cols->operator()(i,j);
      if (index!=SparsityPattern::invalid_entry) {
        val[index]=value;
      }
#ifdef DEBUG
      else CHECK_TEST(value==val[index]);
#endif
    }
    double rawEntry(int row,int index) const {
      CHECK_TEST(row<cols->nRows());
      CHECK_TEST(index<cols->rowLength(row));
      return val[cols->getRowstartIndex(row)+index];
    }
    const SparsityPattern& getSparsityPattern() const { 
      CHECK_POINTER(cols);
      return *cols;
    }
    NumPtr<double>& valueArray() { return val; }
    const NumPtr<double>& valueArray() const { return val; }
    SparseMatrix& operator*=(double factor) {
      CHECK_POINTER(cols);
      CHECK_POINTER(val.getData());
      for (int i=0;i<val.getNumber();i++) val[i]*=factor;
      return *this;
    }
    SparseMatrix& operator/=(double factor) {
      CHECK_POINTER(cols);
      CHECK_POINTER(val.getData());
      CHECK_TEST(abs(factor)>0.);
      double factor_inv=1./factor;
      for (int i=0;i<val.getNumber();i++) val[i]*=factor_inv;
      return *this;
    }
    void add(int i,int j,double value) {
      CHECK_POINTER(cols);
      int index=cols->operator()(i,j);
      CHECK_TEST((index!=SparsityPattern::invalid_entry) || 
        abs(value)==0.);
      if (abs(value)>0.) val[index]+=value;
    }
    void symmetrize();
    SparseMatrix& copyFrom(const SparseMatrix &matrix) {
      CHECK_POINTER(cols);
      CHECK_POINTER(val.getData());
      CHECK_SAME(cols,matrix.cols);
      val.copyFrom(matrix.val);
      return *this;
    }
    double operator()(int i,int j) const {
      CHECK_POINTER(cols);
      int index=cols->operator()(i,j);
      CHECK_TEST(index!=SparsityPattern::invalid_entry);
      return val[index];
    }
    double& operator()(int i,int j) {
      CHECK_POINTER(cols);
      int index=cols->operator()(i,j);
      CHECK_TEST(index!=SparsityPattern::invalid_entry);
      return val[index];
    }
    double el(int i,int j) const {
      CHECK_POINTER(cols);
      int index=cols->operator()(i,j);
      if (index!=SparsityPattern::invalid_entry) return val[index];
      return 0.;
    }
    double diagElement(int row) const {
      CHECK_POINTER(cols);
      int entry=cols->diagonalEntry(row);
      return (entry==SparsityPattern::invalid_entry ? bogus : val[entry]);
    }
    double& diagElement(int row) {
      CHECK_POINTER(cols);
      int entry=cols->diagonalEntry(row);
      return (entry==SparsityPattern::invalid_entry ? bogus : val[entry]);
    }
    double globalEntry(int j) const { return val[j]; }
    double& globalEntry(int j) { return val[j]; }
    int nActuallyNonzeroElements() const;
    int bandwidth() const { return cols->bandwidth(); }
    void copyFrom(const Matrix&);
    void addScaled(double,const SparseMatrix&);
    void vmult(Vector&,const Vector&) const;
    void transposeVmult(Vector&,const Vector&) const;
    void vmultAdd(Vector&,const Vector&) const;
    void transposeVmultAdd(Vector&,const Vector&) const;
    double matrixNormSquare(const Vector&) const;
    double matrixScalarProduct(const Vector&,const Vector&) const;
    double l1Norm() const;
    double linftyNorm() const;
    double residual(Vector&,const Vector&,const Vector&) const;
    void preconditionJacobi(Vector&,const Vector&,double=1.) const;
    void preconditionSSOR(Vector&,const Vector&,double=1.) const;
    void SOR(Vector&,double=1.) const;
    void preconditionSOR(Vector &dst,const Vector &src,double om=1.) 
    const {
      CHECK_POINTER(cols);
      CHECK_POINTER(val.getData());
      CHECK_SAME(cols->nRows(),cols->nCols());
      dst=src;
      SOR(dst,om);
    }
    void TSOR(Vector&,double=1.) const;
    void preconditionTSOR(Vector &dst,const Vector &src,double om=1.) 
    const {
      CHECK_POINTER(cols);
      CHECK_POINTER(val.getData());
      CHECK_SAME(cols->nRows(),cols->nCols());
      dst=src;
      TSOR(dst,om);
    }
    void PSOR(Vector&,const NumPtr<int>&,const NumPtr<int>&,double=1.)
      const;
    void TPSOR(Vector&,const NumPtr<int>&,const NumPtr<int>&,double=1.)
      const;
    void SORStep(Vector&,const Vector&,double=1.) const;
    void TSORStep(Vector&,const Vector&,double=1.) const;
    void SSORStep(Vector &v,const Vector &b,double om=1.) const {
      SORStep(v,b,om);
      TSORStep(v,b,om);
    }
    void SSOR(Vector&,double=1.) const;
    void printOn(ostream&) const;
};

#endif
