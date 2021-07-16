#include "NonlinearSystemSolver.H"

Vector<double,double>* 
  (*nonlinear_system_fcn)(const Vector<double,double> &x);
Matrix<double,double>* 
  (*nonlinear_system_jac)(const Vector<double,double> &x);

extern "C" {
  void F77NAME(hybrd1)(
    void (*f)(const int &n,const double *x,double *fvec,const int &iflag),
    const int &n,double *x,double *fvec,const double &tol,int &info,
    double *wa,const int &lwa);
  void F77NAME(hybrj1)(
    void (*f)(const int &n,const double *x,double *fvec,double *fjac,
      const int &ldfjac,const int &iflag),
    const int &n,double *x,double *fvec,const double *fjac,
    const int &ldfjac,const double &tol,int &info,
    double *wa,const int &lwa);
}

void minpack_f(const int &n,const double *x,double *fvec,
const int &iflag) {
  Vector<double,double> xx(n);
  memcpy(xx.addr(),x,n*sizeof(double));
  Vector<double,double> *fx=nonlinear_system_fcn(xx);
  memcpy(fvec,fx->addr(),n*sizeof(double));
  delete fx; fx=0;
}

void minpack_fj(const int &n,const double *x,double *fvec,
double *fjac,const int &ldfjac,const int &iflag) {
  if (nonlinear_system_jac == 0 || iflag == 1 ) {
    Vector<double,double> xx(n);
    memcpy(xx.addr(),x,n*sizeof(double));
    Vector<double,double> *fx=nonlinear_system_fcn(xx);
    memcpy(fvec,fx->addr(),n*sizeof(double));
    delete fx; fx=0;
  } else if (iflag == 2) {
    CHECK_POINTER(nonlinear_system_jac)
    Vector<double,double> xx(n);
    memcpy(xx.addr(),x,n*sizeof(double));
    Matrix<double,double> *jx=nonlinear_system_jac(xx);
    memcpy(fjac,jx->addr(),n*n*sizeof(double));
    delete jx; jx=0;
  }
}


void MinpackSolver::solve(const double &tol,Vector<double,double> &x,
Vector<double,double>* (*f)(const Vector<double,double> &x),
Matrix<double,double>* (*j)(const Vector<double,double> &x)) {
  nonlinear_system_fcn=f;
  nonlinear_system_jac=j;
  int n=x.size();
  int info=0;
  Vector<double,double> fvec(n);
  if (nonlinear_system_jac == 0) {
    int lwa = (n*(3*n+13))/2;
    double *wa=OPERATOR_NEW_BRACKET(double,lwa);
    F77NAME(hybrd1)(minpack_f,n,x.addr(),fvec.addr(),tol,info,wa,lwa);
    delete wa; wa=0;
  } else {
    Matrix<double,double> fjac(n,n);
    int lwa = (n*(n+13))/2;
    double *wa=OPERATOR_NEW_BRACKET(double,lwa);
    F77NAME(hybrj1)(minpack_fj,n,x.addr(),fvec.addr(),fjac.addr(),n,
      tol,info,wa,lwa);
    delete wa; wa=0;
  }
}
