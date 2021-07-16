#ifndef SOLVER_ARRAYS_H
#define SOLVER_ARRAYS_H "$Header:$"

//**********************************************************************
// Copyright 2006 John A. Trangenstein
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

#include "NumPtr.H"

class Vector;
class SolverMatrix {
  public:
    SolverMatrix() {}
    virtual ~SolverMatrix() {}
    virtual int size(int) const=0;
    virtual double operator()(int,int) const=0;
    virtual double& operator()(int,int)=0;
    virtual double matrixNormSquare(const Vector&) const=0;
    virtual void vmult(Vector&,const Vector&) const=0;
    virtual void transposeVmult(Vector&,const Vector&) const=0;
    virtual void vmultAdd(Vector&,const Vector&) const=0;
    virtual void transposeVmultAdd(Vector&,const Vector&) const=0;
    virtual double residual(Vector&,const Vector&,const Vector&) 
      const=0;
    virtual void preconditionJacobi(Vector&,const Vector&, double) 
      const=0;
    virtual void PSOR(Vector&,const NumPtr<int>&,const NumPtr<int>&,
      double) const=0;
    virtual void preconditionSOR(Vector&,const Vector&,double) const=0;
    virtual void preconditionSSOR(Vector&,const Vector&,double) const=0;
    virtual void TPSOR(Vector&,const NumPtr<int>&,const NumPtr<int>&,
      double) const=0;
    virtual void preconditionTSOR(Vector&,const Vector&,double) const=0;
    virtual SolverMatrix* clone() const=0;
    virtual SolverMatrix* cloneEmpty() const=0;
    virtual void printOn(ostream&) const=0;
};

#endif
