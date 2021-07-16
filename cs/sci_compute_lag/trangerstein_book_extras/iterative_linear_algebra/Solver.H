#ifndef SOLVER_H
#define SOLVER_H "$Header:$"
//---------------------------  solver_control.h  ---------------------------
//   $Id: solver_control.h,v 1.29 2003/01/08 17:58:14 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//---------------------------  solver_control.h  ---------------------------
//----------------------------  solver.h  ---------------------------
//   $Id: solver.h,v 1.29 2003/02/10 09:10:08 guido Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  solver.h  ---------------------------
//--------------------------  solver_bicgstab.h  ---------------------------
//   $Id: solver_bicgstab.h,v 1.47 2003/01/08 17:58:14 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//--------------------------  solver_bicgstab.h  ---------------------------
//----------------------------  solver_cg.h  ---------------------------
//   $Id: solver_cg.h,v 1.48 2003/02/14 16:35:45 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  solver_cg.h  ---------------------------
//----------------------------  solver_gmres.h  ---------------------------
//   $Id: solver_gmres.h,v 1.64 2003/05/18 02:26:47 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  solver_gmres.h  ---------------------------
//----------------------------  solver_minres.h  ---------------------------
//   $Id: solver_minres.h,v 1.22 2003/04/26 16:10:31 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  solver_minres.h  -----------------------
//----------------------------  solver_qmrs.h  ---------------------------
//   $Id: solver_qmrs.h,v 1.31 2003/04/30 23:07:00 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//----------------------------  solver_qmrs.h  ---------------------------
//------------------------  solver_richardson.h  ---------------------------
//   $Id: solver_richardson.h,v 1.35 2003/04/26 16:11:11 wolf Exp $
//   Version: $Name: Version-4-0-0 $
//
//   Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//   This file is subject to QPL and may not be  distributed
//   without copyright and license information. Please refer
//   to the file deal.II/doc/license.html for the  text  and
//   further information on this license.
//
//-------------------------  solver_richardson.h  --------------------------
//
//modified from deal.II/lac/include/lac/solver_control, solver.h,
//  solver_bicgstab.h, solver_cg.h, solver_gmres.h, solver_minres.h,
//  solver_qmrs.h and solver_richardson.h
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

#include <float.h>
#include <limits>
#include <math.h>

#include "Matrix.H"
#include "MemoryDebugger.H"
#include "NumPtr.H"
#include "SolverArrays.H"

enum SOLVER_STATE { SOLVER_ITERATE=0,SOLVER_SUCCESS,SOLVER_FAILURE };
ostream& operator<<(ostream&,SOLVER_STATE);

//SOLVER_SUCCESS <==> check_value <= tol
//compare to DEAL SolverControl
class SolverControl {
  protected:
    int maxsteps;
    double tol;
    bool check_failure;
    double relative_failure_residual;
    double failure_residual;
    bool m_log_history;
    int m_log_frequency;
    bool m_log_result;
  public:
    SolverControl(int m=1,double t=DBL_EPSILON,bool lh=false,bool lr=true) : 
      maxsteps(m),tol(t),check_failure(false),
      relative_failure_residual(0.),failure_residual(0.),m_log_history(lh),
      m_log_frequency(1),m_log_result(lr) { }
    virtual ~SolverControl() {}
    virtual SOLVER_STATE check(int,double);
    int maxSteps() const { return maxsteps; }
    int setMaxSteps(int newval) {
      int old=maxsteps; maxsteps=newval; return old;
    }
    void setFailureCriterion(double rel_failure_residual) {
      relative_failure_residual=rel_failure_residual;
      check_failure=true;
    }
    void clearFailureCriterion() {
      relative_failure_residual=0; 
      failure_residual=0;
      check_failure=false;
    }
    double tolerance() const { return tol; }
    double setTolerance(double t) { double old=tol; tol=t; return old; }
    void logHistory(bool newval) { m_log_history = newval; }
    bool logHistory() const { return m_log_history; }
    int logFrequency(int f) {
      if (f==0) f=1; int old=m_log_frequency;
      m_log_frequency=f; return old;
    }
    void logResult(bool newval) { m_log_result=newval; }
    bool logResult() const { return m_log_result; }
    virtual void printOn(ostream&) const;
};

//SOLVER_SUCCESS <==> check_value <= tol*initial_val
//compare to DEAL ReductionControl
class ReductionControl : public SolverControl {
  protected:
    double initial_val;
  public:
    ReductionControl(int n=1,double red=1.e-2,bool m_log_history=false,
    bool m_log_result=true) :
      SolverControl(n,red,m_log_history,m_log_result) {}
    virtual ~ReductionControl() {}
    virtual SOLVER_STATE check(int,double);
    double setReduction(double t) { return setTolerance(t); }
    virtual void printOn(ostream&) const;
};

//SOLVER_SUCCESS <==> check_value <= tol * rhs_norm
//not in DEAL
class RelativeErrorControl : public SolverControl {
  protected:
    double rhs_norm;
  public:
    RelativeErrorControl(int n=1,double tol=DBL_EPSILON,
      bool m_log_history=false,bool m_log_result=true) :
      SolverControl(n,tol,m_log_history,m_log_result),rhs_norm(0.) {}
    virtual ~RelativeErrorControl() {}
    virtual SOLVER_STATE check(int,double);
    double setRhsNorm(double r) {
      double old=rhs_norm;
      rhs_norm=r;
      return old;
    }
    virtual void printOn(ostream&) const;
};

class SolverPreconditioner;

//compare to DEAL Solver
class Solver { 
  protected:
    SolverControl &cntrl;
  public:
    Solver(SolverControl &cn) : cntrl(cn) {}
    virtual ~Solver() {}
    SolverControl& control() const { return cntrl; }
    virtual void solve(const SolverMatrix&,Vector&,const Vector&,
      const SolverPreconditioner&)=0;
    virtual void printOn(ostream&) const;
};

//compare to DEAL SolverBicgstab
class BicgstabSolver : public Solver {
  private:
    double alpha;
    double beta;
    double omega;
    double rho;
    double rhobar;
    double res;
    bool exact_residual;
    double breakdown;

    double criterion(const SolverMatrix &A,const Vector &x,
    const Vector &b,Vector &t) {
      A.vmult(t,x); 
      t.scale(-1.);
      t.add(1.,b);
      res=t.linftyNorm();
      return res;
    }
  public:
    BicgstabSolver(SolverControl &cn,
    bool er=true,double bd=1./sqrt(DBL_EPSILON)) : Solver(cn /* ,mem */ ),
      alpha(numeric_limits<double>::infinity()),beta(numeric_limits<double>::infinity()),omega(numeric_limits<double>::infinity()),
      rho(numeric_limits<double>::infinity()),rhobar(numeric_limits<double>::infinity()),res(numeric_limits<double>::infinity()),exact_residual(er),
      breakdown(bd) {}
    virtual ~BicgstabSolver() {}
    void solve(const SolverMatrix&,Vector&,const Vector&,
      const SolverPreconditioner&);
    virtual void printOn(ostream&) const;
};

//compare to DEAL SolverCG
class CGSolver : public Solver {
  private:
    double res2;
  public:
    CGSolver(SolverControl &cn) : Solver(cn) {}
    virtual ~CGSolver() {}
    void solve(const SolverMatrix&,Vector&,const Vector&,
      const SolverPreconditioner&);
    virtual void printOn(ostream&) const;
};

//compare to DEAL SolverGMRES
class GMRESSolver : public Solver {
  private:
    int max_n_tmp_vectors;
    bool right_preconditioning;
    bool use_default_residual;
    virtual double criterion() {
      OBSOLETE("criterion meaningless for GMRES");
      return 0;
    }
    void givensRotation(Vector&,Vector&,Vector&,Vector&,int) const;
    GMRESSolver(const GMRESSolver&);
  public:
    GMRESSolver(SolverControl &cn,int mntv=30,bool rp=false,bool udr=true) :
      Solver(cn),max_n_tmp_vectors(mntv),
      right_preconditioning(rp),use_default_residual(udr) {}
    virtual ~GMRESSolver() {}
    template<class T> SOLVER_STATE templateSolve(const T&,
      void (T::*)(Vector&,const Vector&) const,
      Vector&,const Vector&,const SolverPreconditioner&);
    template<class T> SOLVER_STATE daskrSolve(const T&,
      void (T::*)(Vector&,const Vector&) const,
      Vector&,const Vector&,const SolverPreconditioner&);
    void solve(const SolverMatrix &A,Vector &x,const Vector &b,
    const SolverPreconditioner &p) {
      SOLVER_STATE iteration_state=
        templateSolve(A,&SolverMatrix::vmult,x,b,p);
//    CHECK_TEST(iteration_state==SOLVER_SUCCESS);
    }
    virtual void printOn(ostream&) const;
};

//compare to DEAL SolverMinRes
class MinResSolver : public Solver {
  private:
    double res2;

    virtual double criterion() { return res2; }
  public:
    MinResSolver(SolverControl &cn) : Solver(cn) {}
    virtual ~MinResSolver() {}
    void solve(const SolverMatrix&,Vector&,const Vector&,
      const SolverPreconditioner&);
    virtual void printOn(ostream&) const;
};

//compare to DEAL SolverQMRS
class QMRSSolver : public Solver {
  private:
    double res2;
    bool exact_residual;
    double breakdown;
    virtual double criterion() { return res2; }
  public:
    QMRSSolver(SolverControl &cn,bool er=false,double bd=DBL_EPSILON): 
      Solver(cn),res2(numeric_limits<double>::infinity()),exact_residual(er),breakdown(bd) {}
    void solve(const SolverMatrix&,Vector&,const Vector&,
      const SolverPreconditioner&);
    virtual void printOn(ostream&) const;
};

//compare to DEAL SolverRichardson
class RichardsonSolver : public Solver {
  private:
//  must have omega < 2/spectral radius(A)
    double omega;
    bool use_preconditioned_residual;
  public:
    RichardsonSolver(SolverControl &cn,double om,bool upr) : 
      Solver(cn),omega(om),use_preconditioned_residual(upr) {}
    virtual ~RichardsonSolver() {}
    void solve(const SolverMatrix&,Vector&,const Vector&,
      const SolverPreconditioner&);
    void transposeSolve(const SolverMatrix&,Vector&,const Vector&,
      const SolverPreconditioner&);
    void setOmega(double om=1.) { omega=om; }
    virtual void printOn(ostream&) const;
};

//not in DEAL
class HdivCGSolver : public Solver {
  private:
    double res2;
    void solve(const SolverMatrix&,Vector&,const Vector&,
    const SolverPreconditioner&) {
      OBSOLETE("not appropriate for this class");
    }
  public:
    HdivCGSolver(SolverControl &cn) : Solver(cn) {}
    virtual ~HdivCGSolver() {}
    void solve(const SolverMatrix&,const SolverMatrix&,Vector&,
      const Vector&,const SolverPreconditioner&);
    virtual void printOn(ostream&) const;
};

#endif
