#ifndef SHIFTED_MATRIX_H
#define SHIFTED_MATRIX_H "$Header:$"
//----------------------------  shifted_matrix.h  --------------------------
//    $Id: shifted_matrix.h,v 1.7 2003/01/08 17:58:14 wolf Exp $
//    Version: $Name: Version-4-0-0 $
//
//    Copyright (C) 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  shifted_matrix.h  --------------------------
//
//modified from deal.II/lac/include/lac/shifted_matrix.h
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

#include <Solver.H>

//compare to DEAL ShiftedMatrix
class ShiftedMatrix : public SolverMatrix {
  private:
    const SolverMatrix *A;
    REAL sigma;
  public:
    ShiftedMatrix(const SolverMatrix &rA,REAL s) : A(&rA), sigma(s) {}
    void setShift(REAL s) { sigma = s; }
    REAL getShift() const { return sigma; }
    int size(int i) const { return A->size(i); }
    virtual REAL operator()(int i,int j) const { 
      return (i==j ? (*A)(i,j)+sigma : (*A)(i,j));
    }
    virtual REAL& operator()(int i,int j) { 
      OBSOLETE("should not call this");
      return sigma;
    }
    REAL matrixNormSquare(const Vector &x) const {
      REAL r=A->matrixNormSquare(x);
      return r+(x*x)*sigma;
    }
    void vmult(Vector &dst,const Vector &src) const {
      A->vmult(dst,src);
      if (abs(sigma)!=0.) dst.add(sigma,src);
    }
    virtual void transposeVmult(Vector &dst,const Vector &src) const {
      A->transposeVmult(dst,src);
      if (abs(sigma)!=0.) dst.add(sigma,src);
    }
    void vmultAdd(Vector &dst,const Vector &src) const {
      A->vmultAdd(dst,src);
      if (abs(sigma)!=0.) dst.add(sigma,src);
    }
    void transposeVmultAdd(Vector &dst,const Vector &src) const {
      A->transposeVmultAdd(dst,src);
      if (abs(sigma)!=0.) dst.add(sigma,src);
    }
    void PSOR(Vector &dst,const NumPtr<int> &permutation,
    const NumPtr<int> &inverse_permutation,REAL om) const {
      OBSOLETE("not programmed");
    }
    void TPSOR(Vector &dst,
    const NumPtr<int> &permutation,
    const NumPtr<int> &inverse_permutation,REAL om) const {
      OBSOLETE("not programmed");
    }
    void preconditionJacobi(Vector &dst,const Vector &src,REAL om) 
    const {
      OBSOLETE("not programmed");
    }
    void preconditionSOR(Vector &dst,const Vector &src,REAL om) const {
      OBSOLETE("not programmed");
    }
    void preconditionSSOR(Vector &dst,const Vector &src,REAL om) 
    const {
      OBSOLETE("not programmed");
    }
    void preconditionTSOR(Vector &dst,const Vector &src,REAL om) 
    const {
      OBSOLETE("not programmed");
    }
    REAL residual(Vector &dst,const Vector& src,const Vector &rhs) 
    const {
      A->vmult(dst,src);
      if (abs(sigma)!=0.) dst.add(sigma,src);
      dst.scale(-1.);
      dst.add(1.,rhs);
      return dst.l2Norm();
    }
    ShiftedMatrix* clone() const {
      return OPERATOR_NEW ShiftedMatrix(*A,sigma);
    }
    ShiftedMatrix* cloneEmpty() const {
      OBSOLETE("not programmed");
    }
    void printOn(ostream &os) const {
      os << "ShiftedMatrix: sigma = " << sigma << endl; 
      A->printOn(os);
    }
};
    
//compare to DEAL ShiftedMatrixGeneralized
class ShiftedMatrixGeneralized : public SolverMatrix {
  private:
    const SolverMatrix *A;
    const SolverMatrix *M;
    REAL sigma;
  public:
    ShiftedMatrixGeneralized(const SolverMatrix &rA,
      const SolverMatrix &rM,REAL s) : A(&rA), M(&rM), sigma(s) {}
    void setShift(REAL s) { sigma = s; }
    REAL getShift() const { return sigma; }
    void vmult(Vector &dst,const Vector &src) const {
      A->vmult(dst,src);
      if (abs(sigma)!=0.) {
        Vector *aux=dst.clone();
        M->vmult(*aux,src);
        dst.add(sigma,*aux);
        delete aux;
      }
    }
    REAL residual(Vector &dst,const Vector &src,const Vector &rhs) 
    const {
      A->vmult(dst,src);
      if (abs(sigma)!=0.) {
        Vector *aux=dst.clone();
        M->vmult(*aux,src);
        dst.add(sigma,*aux);
        delete aux;
      }
      dst.scale(-1.);
      dst.add(1.,rhs);
      return dst.l2Norm();
    }
};
    
#endif
