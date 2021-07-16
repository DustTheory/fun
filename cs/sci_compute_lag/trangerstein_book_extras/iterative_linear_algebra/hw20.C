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
#include "SparsityPattern.H"
#include "SparseMatrix.H"
//#include "Tracer.H"

double omega=1.;

void debugSolver(Solver *solver) {
  ifstream in_file;
  in_file.open("bcsstk14.mtx",ios::in);
  int m,n,nonzero;
  in_file >> m >> n >> nonzero;
  SparsityPattern sp(m,n,48);
  for (int k=0;k<nonzero;k++) {
    int i,j;
    double aij;
    in_file >> i >> j >> aij;
    sp.add(i-1,j-1);
  }
  in_file.close();
  sp.symmetrize();
  sp.compress();

  SparseMatrix A(&sp);
  in_file.open("bcsstk14.mtx",ios::in);
  in_file >> m >> n >> nonzero;
  for (int k=0;k<nonzero;k++) {
    int i,j;
    double aij;
    in_file >> i >> j >> aij;
    A.set(i-1,j-1,aij);
    if (i!=j) A.set(j-1,i-1,aij);
  }
  in_file.close();

/*
  int m=20,n=20;
  SparsityPattern sp(m,n,3);
  for (int k=0;k<m;k++) {
    sp.add(k,k);
    if (k>0) sp.add(k-1,k);
  }
  sp.symmetrize();
  sp.compress();
  SparseMatrix A(&sp);
  for (int k=0;k<m;k++) {
    A.set(k,k,2.);
    if (k>0) {
      A.set(k-1,k,-1.);
      A.set(k,k-1,-1.);
    }
  }
*/

  solver->control().setMaxSteps(2*m);
  Vector b(m,0.);
  Vector x(n,1.);
  A.vmult(b,x);
  Vector r(m,numeric_limits<double>::infinity());
  for (int j=0;j<n;j++) x[j]=drand48();
//x.printOn(cout);

  PreconditionIdentity p;
//PreconditionUseMatrix p(M,M.solve);
//PreconditionJacobi p(M,omega);
//PreconditionSOR p(M,omega);
//PreconditionSSOR p(M,omega);
//PreconditionLACSolver p(s,M,p);

  solver->solve(A,x,b,p);
  double rnorm=A.residual(r,x,b);
  cout << "\trnorm = " << rnorm << endl;
  for (int j=0;j<n;j++) x[j]-=1.;
  cout << "\tnorm x error = " << x.linftyNorm() << endl;
}

int main(int argc,char* argv[]) {
#ifdef DEBUG
//setTraps();
#endif
  {
#ifdef MEM_DEBUG
    MemoryDebugger md(1);
#endif
    ReductionControl sc(1,1.e-6);
//
    {
      TRACER_CALL(t,"CGSolver");
      CGSolver cgs(sc);
      debugSolver(&cgs);
    }
//
//
    {
      TRACER_CALL(t,"BicgstabSolver");
      BicgstabSolver bs(sc);
      debugSolver(&bs);
    }
//
//
    {
      TRACER_CALL(t,"GMRESSolver");
//    GMRESSolver gms(sc);
      GMRESSolver gms(sc,5);
      debugSolver(&gms);
    }
//
//
    {
      TRACER_CALL(t,"MinResSolver");
      MinResSolver mrs(sc);
      debugSolver(&mrs);
    }
//
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
