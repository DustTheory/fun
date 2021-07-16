#ifndef PRECONDITION_H
#define PRECONDITION_H "$Header:$"
//----------------------------  precondition.h  ---------------------------
//   $Id: precondition.h,v 1.37 2003/04/25 12:52:32 guido Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  precondition.h  ---------------------------
//
//modified from deal.II/lac/include/lac/precondition.h
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

#include "SolverArrays.H"
#include "Vector.H"

class Solver;
class SolverMatrix;
class Vector;
//not in DEAL
class SolverPreconditioner {
  public:
    SolverPreconditioner() {}
    virtual ~SolverPreconditioner() {}
    virtual SolverPreconditioner* clone() const {
      OBSOLETE("not programmed by any derived class");
    }
    virtual SolverPreconditioner* cloneEmpty() const {
      OBSOLETE("not programmed by any derived class");
    }
    virtual void setMatrix(const SolverMatrix&) {
      OBSOLETE("not programmed by any derived class");
    }
    virtual void vmult(Vector&,const Vector&) const {
      OBSOLETE("not programmed by any derived class");
    }
    virtual void transposeVmult(Vector&,const Vector&) const {
      OBSOLETE("not programmed by any derived class");
    }
    virtual void printOn(ostream&) const;
};

//compare to DEAL PreconditionIdentity
class PreconditionIdentity : public SolverPreconditioner {
  public:
    PreconditionIdentity() {}
    virtual ~PreconditionIdentity() {}
    virtual PreconditionIdentity* clone() const { 
      return OPERATOR_NEW PreconditionIdentity; 
    }
    virtual PreconditionIdentity* cloneEmpty() const { 
      return OPERATOR_NEW PreconditionIdentity; 
    }
    virtual void setMatrix(const SolverMatrix&) {;}
    virtual void vmult(Vector &dst,const Vector &src) const { dst=src; }
    virtual void transposeVmult(Vector &dst,const Vector &src) const { 
      dst=src; 
    }
    virtual void printOn(ostream &os) const;
};

//compare to DEAL PreconditionUseMatrix
class PreconditionUseMatrix : public SolverPreconditioner {
  public:
    typedef void (*function_ptr)(Vector&,const Vector&);
  private:
    const SolverMatrix *matrix;
    const function_ptr precondition;
  public:
    PreconditionUseMatrix(const SolverMatrix &M,
      const function_ptr method) : matrix(&M),precondition(method) {}
    PreconditionUseMatrix(const function_ptr method) : matrix(0),
      precondition(method) {}
    virtual ~PreconditionUseMatrix() {}
    PreconditionUseMatrix* clone() const { 
      return OPERATOR_NEW PreconditionUseMatrix(precondition); 
    }
    PreconditionUseMatrix* cloneEmpty() const { 
      return OPERATOR_NEW PreconditionUseMatrix(precondition); 
    }
    void setMatrix(const SolverMatrix &m) { matrix=&m; }
    void vmult(Vector &dst,const Vector &src) const {
      (*precondition)(dst,src);
    }
    virtual void printOn(ostream &os) const;
};

//compare to DEAL PreconditionRelaxation
class PreconditionRelaxation : public SolverPreconditioner {
  protected:
    const SolverMatrix *A;
    double relaxation;
  public:
    PreconditionRelaxation(const SolverMatrix &rA,double r=1.) : A(&rA),
      relaxation(r) { }
    PreconditionRelaxation(double r=1.) : A(0),relaxation(r) { }
    virtual ~PreconditionRelaxation() {}
    void setMatrix(const SolverMatrix &m) { A=&m; }
    virtual void printOn(ostream &os) const;
};

// this leads to non-convergence of multigrid with r=1
// instead, take r = 1 / max number nonzero entries in a row of the matrix
//compare to DEAL PreconditionJacobi
class PreconditionJacobi : public PreconditionRelaxation {
  public:
    PreconditionJacobi(const SolverMatrix &rA,double r=1.) :
      PreconditionRelaxation(rA,r) {}
    PreconditionJacobi(double r=1.) : PreconditionRelaxation(r) {}
    virtual ~PreconditionJacobi() {}
    PreconditionJacobi* clone() const {
      return OPERATOR_NEW PreconditionJacobi(relaxation);
    }
    PreconditionJacobi* cloneEmpty() const {
      return OPERATOR_NEW PreconditionJacobi(relaxation);
    }
    void vmult(Vector &dst,const Vector &src) const {
      CHECK_POINTER(A);
      A->preconditionJacobi(dst,src,relaxation);
    }
    void transposeVmult(Vector &dst,const Vector &src) const {
      CHECK_POINTER(A);
      A->preconditionJacobi(dst,src,relaxation);
    }
    virtual void printOn(ostream &os) const;
};

//compare to DEAL PreconditionSOR
class PreconditionSOR : public PreconditionRelaxation {
  public:
    PreconditionSOR(const SolverMatrix &rA,double r=1.) :
      PreconditionRelaxation(rA,r) {}
    PreconditionSOR(double r=1.) : PreconditionRelaxation(r) {}
    virtual ~PreconditionSOR() {}
    PreconditionSOR* clone() const {
      return OPERATOR_NEW PreconditionSOR(relaxation);
    }
    PreconditionSOR* cloneEmpty() const {
      return OPERATOR_NEW PreconditionSOR(relaxation);
    }
    void vmult(Vector &dst,const Vector &src) const {
      CHECK_POINTER(A);
      A->preconditionSOR(dst,src,relaxation);
    }
    void transposeVmult(Vector &dst,const Vector &src) const
    {
      CHECK_POINTER(A);
      A->preconditionTSOR(dst,src,relaxation);
    }
    virtual void printOn(ostream &os) const;
};

//compare to DEAL PreconditionSSOR
class PreconditionSSOR : public PreconditionRelaxation {
  public:
    PreconditionSSOR(const SolverMatrix &rA,double r=1.) :
      PreconditionRelaxation(rA,r) {}
    PreconditionSSOR(double r=1.) : PreconditionRelaxation(r) {}
    virtual ~PreconditionSSOR() {}
    PreconditionSSOR* clone() const {
      return OPERATOR_NEW PreconditionSSOR(relaxation);
    }
    PreconditionSSOR* cloneEmpty() const {
      return OPERATOR_NEW PreconditionSSOR(relaxation);
    }
    void vmult(Vector &dst,const Vector &src) const {
      CHECK_POINTER(A);
      A->preconditionSSOR(dst,src,relaxation);
    }
    void transposeVmult(Vector &dst,const Vector &src) const {
      CHECK_POINTER(A);
      A->preconditionSSOR(dst,src,relaxation);
    }
    virtual void printOn(ostream &os) const;
};

//compare to DEAL PreconditionPSOR
class PreconditionPSOR : public PreconditionRelaxation {
  private:
    const NumPtr<int> &permutation;
    const NumPtr<int> &inverse_permutation;
  public:
    PreconditionPSOR(const SolverMatrix &rA,
    const NumPtr<int> &p,const NumPtr<int> &ip,double r=1.) : 
      PreconditionRelaxation(rA,r),permutation(p),
      inverse_permutation(ip) {}
    virtual ~PreconditionPSOR() {}
    void vmult(Vector &dst,const Vector &src) const {
      CHECK_POINTER(A);
      dst=src;
      A->PSOR(dst,permutation,inverse_permutation,relaxation);
    }
    void transposeVmult(Vector &dst,const Vector &src) const {
      CHECK_POINTER(A);
      dst=src;
      A->TPSOR(dst,permutation,inverse_permutation,relaxation);
    }
    virtual void printOn(ostream &os) const;
};

//compare to DEAL PreconditionLACSolver
class PreconditionLACSolver {
  private:
    Solver *solver;
    const SolverMatrix *matrix;
    const SolverPreconditioner *precondition;
  public:
    PreconditionLACSolver(Solver &s,const SolverMatrix &m,
    const SolverPreconditioner &p) : solver(&s),matrix(&m),
      precondition(&p) {}
    virtual ~PreconditionLACSolver() { }
    void vmult(Vector&,const Vector&) const;
    virtual void printOn(ostream &os) const;
};

//compare to DEAL PreconditionedMatrix
class PreconditionedMatrix {
  private:
    const SolverMatrix *A;
    const SolverPreconditioner *P;
  public:
    PreconditionedMatrix(const SolverMatrix *rA,
      const SolverPreconditioner *rP) : A(rA),P(rP) {}
    void vmult(Vector &dst,const Vector &src) const {
      Vector *h=src.clone();
      A->vmult(*h,src);
      P->vmult(dst,*h);
      delete h; h=0;
    }
    void transposeVmult(Vector &dst,const Vector &src) const {
      Vector *h=src.clone();
      A->transposeVmult(*h,src);
      P->transposeVmult(dst,*h);
      delete h; h=0;
    }
    double residual(Vector &dst,const Vector &src,const Vector &rhs) 
    const {
      Vector *h=src.clone();
      A->vmult(*h,src);
      P->vmult(dst,*h);
      delete h; h=0;
      dst.scale(-1.);
      dst.add(1.,rhs);
      return dst.l2Norm();
    }
};

#endif
