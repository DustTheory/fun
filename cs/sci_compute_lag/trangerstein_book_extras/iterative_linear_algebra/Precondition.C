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
#include "Precondition.H"
#include "Solver.H"

void PreconditionIdentity::printOn(ostream &os) const {
  os << "PreconditionIdentity:" << endl;
  SolverPreconditioner::printOn(os);
}

void PreconditionUseMatrix::printOn(ostream &os) const {
  os << "PreconditionUseMatrix:" << endl;
  os << "\tmatrix = " << matrix << endl;
  SolverPreconditioner::printOn(os);
}

void PreconditionRelaxation::printOn(ostream &os) const {
  os << "PreconditionRelaxation:" << endl;
  os << "\trelaxation = " << relaxation << endl;
  os << "\tA = " << A << endl;
  A->printOn(os);
  SolverPreconditioner::printOn(os);
}

void PreconditionJacobi::printOn(ostream &os) const {
  os << "PreconditionJacobi:" << endl;
  PreconditionRelaxation::printOn(os);
}

void PreconditionSOR::printOn(ostream &os) const {
  os << "PreconditionSOR:" << endl;
  PreconditionRelaxation::printOn(os);
}

void PreconditionSSOR::printOn(ostream &os) const {
  os << "PreconditionSSOR:" << endl;
  PreconditionRelaxation::printOn(os);
}

void PreconditionPSOR::printOn(ostream &os) const {
  os << "PreconditionPSOR:" << endl;
  PreconditionRelaxation::printOn(os);
}

void PreconditionLACSolver::vmult(Vector &dst,const Vector &src) const {
  CHECK_POINTER(solver);
  CHECK_POINTER(matrix);
  CHECK_POINTER(precondition);
  solver->solve(*matrix,dst,src,*precondition);
}
void PreconditionLACSolver::printOn(ostream &os) const {
  os << "PreconditionLACSolver:" << endl;
}
