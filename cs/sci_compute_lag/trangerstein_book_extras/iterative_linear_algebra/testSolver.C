// "$Header:$"
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
#include "MemoryDebugger.H"
#include "Precondition.H"
#include "SetTraps.H"
#include "Solver.H"
#include "Tracer.H"

int size=13;
REAL omega=1.;

void debugSolver(Solver *solver) {
  Matrix A(size,size,0.);
  Matrix M(size,size,0.);
  Vector b(size,0.);
  Vector x(size,numeric_limits<double>::infinity());
  Vector r(size,numeric_limits<double>::infinity());
  for (int j=0;j<size;j++) {
    A(j,j)=16.;
    if (j>0) A(j-1,j)=-8.;
    if (j<size-1) A(j+1,j)=-8.;
    M(j,j)=1./A(j,j);
//  x[j]=drand48();
    x[j]=0.;
  }
//b[0]=1.;
  b[size-1]=0.01;

  PreconditionIdentity p;
//PreconditionUseMatrix p(M,M.solve);
//PreconditionJacobi p(M,omega);
//PreconditionSOR p(M,omega);
//PreconditionSSOR p(M,omega);
//PreconditionLACSolver p(s,M,p);

  solver->solve(A,x,b,p);
  REAL rnorm=A.residual(r,x,b);
  cout << "\trnorm = " << rnorm << endl;
}

int main(int argc,char* argv[]) {
#ifdef DEBUG
//setTraps();
#endif
  {
#ifdef MEM_DEBUG
    MemoryDebugger md(1);
#endif
    SolverControl sc(size,DBL_EPSILON*128);
//
    {
      TRACER_CALL(t,"CGSolver");
      CGSolver cgs(sc);
      debugSolver(&cgs);
    }
//
/*
    {
      TRACER_CALL(t,"BicgstabSolver");
      BicgstabSolver bs(sc);
      debugSolver(&bs);
    }
*/
/*
    {
      TRACER_CALL(t,"GMRESSolver");
//    GMRESSolver gms(sc);
      GMRESSolver gms(sc,5);
      debugSolver(&gms);
    }
*/
/*
    {
      TRACER_CALL(t,"MinResSolver");
      MinResSolver mrs(sc);
      debugSolver(&mrs);
    }
*/
/*
    {
      TRACER_CALL(t,"QMRSSolver");
      QMRSSolver qmrs(sc);
      debugSolver(&qmrs);
    }
*/
/*
    {
//    must have omega < 2/spectral radius(A)
      TRACER_CALL(t,"RichardsonSolver");
      RichardsonSolver rs(sc,0.25,false);
      debugSolver(&rs);
    }
*/
  }
  return EXIT_SUCCESS;
}
