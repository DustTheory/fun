#ifndef NONLINEAR_SYSTEM_SOLVER_H
#define NONLINEAR_SYSTEM_SOLVER_H "$Header:$"
#include "Matrix.H"

class NonlinearSystemSolver {
  public:
    NonlinearSystemSolver() {}
    virtual ~NonlinearSystemSolver() {}
    virtual void solve(const double &tol,Vector<double,double> &x,
      Vector<double,double>* (*f)(const Vector<double,double> &x),
      Matrix<double,double>* (*j)(const Vector<double,double> &x) = 0) { }
};

void minpack_f(const int &n,const double *x,double *fvec,
  const int &iflag);
void minpack_fj(const int &n,const double *x,double *fvec,
  double *fjac,const int &ldfjac,const int &iflag);

class MinpackSolver : public NonlinearSystemSolver {
  public:
    MinpackSolver() : NonlinearSystemSolver() {}
    ~MinpackSolver() {}
    virtual void solve(const double &tol,Vector<double,double> &x,
      Vector<double,double>* (*f)(const Vector<double,double> &x),
      Matrix<double,double>* (*j)(const Vector<double,double> &x) = 0);
};
#endif
