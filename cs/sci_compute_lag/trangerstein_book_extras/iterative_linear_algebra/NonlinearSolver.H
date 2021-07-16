#ifndef NONLINEAR_SOLVER_H
#define NONLINEAR_SOLVER_H "$Header:$"
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

#include "Solver.H"
#include "Vector.H"

template<class T> class NonlinearSolver {
  private:
    T &nonlinear_system;
    bool (T::*function)(const Vector &x,Vector &f) const;
    Vector &x;
    Vector f;
    double x_norm;
    double f_norm;
  public:
    NonlinearSolver(T &ns,
    bool (T::*function_in)(const Vector&,Vector&) const,Vector &x_in) : 
    nonlinear_system(ns),function(function_in),x(x_in),
    f(x.size()),x_norm(numeric_limits<double>::infinity()),
    f_norm(numeric_limits<double>::infinity()) {
      (nonlinear_system.*function)(x,f);
      x_norm=x.l2Norm();
      f_norm=f.l2Norm();
    }
    double fNorm() const { return f.l2Norm(); }
    const Vector& solution() const { return x; }
//  df/dx(x) * v approx [ f(x+v*c)-f(x) ]/c
    void jacobianTimesVector(Vector &jv,const Vector &v) const;
//  changes x and f:
    bool lineSearch(double step_tolerance,const Vector &dx);
    void krylovSubspaceIteration(SolverControl&,SolverControl&);
    void minpackHybridIteration(SolverControl&);
};

#endif
