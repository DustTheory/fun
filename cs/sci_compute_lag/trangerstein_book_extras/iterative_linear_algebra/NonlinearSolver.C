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

#include "NonlinearSolver.H"
#include "Precondition.H"
#include "Solver.H"
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  void F77NAME(daxpy)(const int &n,const double &da,const double *dx,
    const int &incx,double *dy,const int &incy);
  void F77NAME(dgesvd)(const char &jobu,const char &jobvt,const int &m,
    const int &n,double *a,const int &lda,double *s,double *u,
    const int &ldu,double *vt,const int &ldvt,double *work,
    const int &lwork,int &info);
  double F77NAME(ddot)(const int &n,const double *dx,const int &incx,
    const double *dy,const int &incy);
  void F77NAME(dogleg)(const int &n,double *r,const int &lr,
    const double *diag,const double *qtb,const double &delta,
    double *x,double *wa1,double *wa2);
  void F77NAME(qrfac)(const int &m,const int &n,double *a,
    const int &lda,const bool &pivot,int *ipvt,const int &lipvt,
    double *rdiag,double *acnorm,double *wa);
  void F77NAME(qform)(const int &m,const int &n,double *q,
    const int &ldq,double *wa);
  void F77NAME(r1mpyq)(const int &m,const int &n,double *a,
    const int &lda,const double *v,const double *w);
  void F77NAME(r1updt)(const int &m,const int &n,double *s,
    const int &ls,const double *u,double *v,double *w,bool &sing);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//  df/dx(x) * v approx [ f(x+v*c)-f(x) ]/c
template<class T> void NonlinearSolver<T>::jacobianTimesVector(
Vector &jv,const Vector &v) const {
  REAL c=x_norm;
  if (c<=0.) c=1.;
  c*=sqrt(DBL_EPSILON);
  Vector *xnew=x.clone();
  (*xnew)=v;
  xnew->scale(c);
  (*xnew)+=x;
  (nonlinear_system.*function)(*xnew,jv);
  jv-=f;
  jv.scale(1./c);
  delete xnew;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> bool NonlinearSolver<T>::lineSearch(
REAL step_tolerance,const Vector &dx) {
  REAL dx_norm=dx.l2Norm();
  bool success=false;
  if (dx_norm>=x_norm*step_tolerance) {
    REAL alpha=1.e-2;
    Vector *xnew=x.clone();
    Vector *fnew=f.clone();
    for (REAL rl=1.;rl*dx_norm>=x_norm*step_tolerance;rl*=0.5) {
      (*xnew)=dx;
      xnew->scale(-rl);
      (*xnew)+=x;
      bool f_ok=(nonlinear_system.*function)(*xnew,*fnew);
      if (f_ok) {
        REAL fnew_norm=fnew->l2Norm();
        if (fnew_norm*fnew_norm<=f_norm*f_norm*(1.-alpha*rl)) {
          success=true;
          x=(*xnew);
          f=(*fnew);
          x_norm=x.l2Norm();
          f_norm=fnew_norm;
          break;
        }
      }
    }
    delete xnew;
    delete fnew;
  }
  return success;
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//adapted from subroutine DNSK of DDASKR
//  ***AUTHORS   Linda R. Petzold, Peter N. Brown, Alan C. Hindmarsh, and
//               Clement W. Ulrich
//  Center for Computational Sciences & Engineering, L-316
//  Lawrence Livermore National Laboratory
//  P.O. Box 808,
//  Livermore, CA 94551
template<class T> void NonlinearSolver<T>::krylovSubspaceIteration(
SolverControl &linear_system_control,
SolverControl &nonlinear_system_control) {
  int neq=x.size();
  SOLVER_STATE solver_state=
    nonlinear_system_control.check(0,f_norm);
  if (solver_state==SOLVER_ITERATE) {
    Vector *dx=x.clone();
//  int kmp=min(5,neq);
    int kmp=neq;
    GMRESSolver gmres_solver(linear_system_control,kmp);
    PreconditionIdentity pi;
    int it=0;
    REAL step_tolerance=pow(DBL_EPSILON,2./3.);
    while (solver_state==SOLVER_ITERATE) {
      it++;
      (*dx)=0.;
//    gmres_solver.daskrSolve(*this,
//      &NonlinearSolver<T>::jacobianTimesVector,*dx,f,pi);
      gmres_solver.templateSolve<NonlinearSolver<T> >(*this,
        &NonlinearSolver<T>::jacobianTimesVector,*dx,f,pi);
      bool line_search_success=lineSearch(step_tolerance,*dx);
      solver_state=nonlinear_system_control.check(it,f_norm);
    }
    delete dx;
  }
  CHECK_TEST(solver_state==SOLVER_SUCCESS);
}
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//adapted from subroutine hybrd:
//argonne national laboratory. minpack project. march 1980.
//burton s. garbow, kenneth e. hillstrom, jorge j. more
template<class T> void NonlinearSolver<T>::minpackHybridIteration(
SolverControl &nonlinear_system_control) {
  double xtol=nonlinear_system_control.tolerance();
  int maxfev=nonlinear_system_control.maxSteps();
  int max_slow1=500;
  int max_slow2=20;
  int max_fail=20;
  int n=x.size();
  int mode=1;
  double factor=100.;
  int nprint=0;
  double epsfcn=DBL_EPSILON;
  int ldfjac=n;
  Matrix fjac(n,n);
  int lr=(n*(n+1))/2;
  Vector diag(n);
  Vector r(lr);
  Vector qtf(n);
  Vector wa1(n),wa2(n),wa3(n),wa4(n);
  Vector xnew(n);
//evaluate the function at the starting point
//and calculate its norm.
  bool f_ok=(nonlinear_system.*function)(x,f);
  CHECK_TEST(f_ok);
  int nfev=1;
  f_norm=f.l2Norm();
//number of fcn calls needed to compute the jacobian matrix.
//initialize iteration counter and monitors.
  int iter=1;
  int ncsuc=0;
  int ncfail=0;
  int nslow1=0;
  int nslow2=0;
//beginning of the outer loop.
  bool converged=false;
  bool failed=false;
  double delta=0.;
  do {
    bool jeval=true;
//  calculate the jacobian matrix.
    double eps=sqrt(max(epsfcn,DBL_EPSILON));
    xnew=x;
    for (int j=0;j<n;j++) {
      double temp=x[j];
      double h=eps*abs(temp);
      if (h<=0.) h=eps;
      xnew[j]=temp+h;
      (nonlinear_system.*function)(xnew,wa1);
      xnew[j]=temp;
      for (int i=0;i<n;i++) fjac(i,j)=(wa1[i]-f[i])/h;
    }
    nfev+=n;
//  compute the qr factorization of the jacobian.
    int *iwa=0;
    bool pivot=false;
    F77NAME(qrfac)(n,n,fjac.addr(),ldfjac,pivot,iwa,1,
      wa1.addr(),wa2.addr(),wa3.addr());
//  on the first iteration and if mode is 1, scale according
//  to the norms of the columns of the initial jacobian.
    if (iter==1) {
      if (mode!=2) {
        for (int j=0;j<n;j++) {
          diag[j]=wa2[j];
          if (wa2[j]<=0.) diag[j]=1.;
        }
      }
//    calculate norm of scaled x and initialize step bound delta.
      for (int j=0;j<n;j++) wa3[j]=diag[j]*x[j];
      x_norm=wa3.l2Norm();
      delta=factor*x_norm;
      if (delta<=0.) delta=factor;
      if (x_norm<=0.) {
        double max_r_diag=0.;
        for (int i=0;i<n;i++) {
          max_r_diag=max(max_r_diag,abs(wa1[i]));
        }
        x_norm=f_norm/max_r_diag;
      }
    }
//  form (q transpose)*f and store in qtf.
    qtf=f;
    for (int j=0;j<n;j++) {
      if (abs(fjac(j,j))>0.) {
        double sum=F77NAME(ddot)(n-j,&fjac(j,j),1,&qtf[j],1);
        double temp=-sum/fjac(j,j);
        F77NAME(daxpy)(n-j,temp,&fjac(j,j),1,&qtf[j],1);
      }
    }
//  copy the triangular factor of the qr factorization into r.
    bool sing=false;
    for (int j=0;j<n;j++) {
      int l=j;
      for (int i=0;i<j;i++,l+=n-i) r[l]=fjac(i,j);
      r[l]=wa1[j];
      if (wa1[j]<=0.) sing=true;
    }
//  accumulate the orthogonal factor in fjac.
    F77NAME(qform)(n,n,fjac.addr(),ldfjac,wa1.addr());
//  rescale if necessary.
    if (mode!=2) {
      for (int j=0;j<n;j++) diag[j]=max(diag[j],wa2[j]);
    }
//  beginning of the inner loop.
    do {
//    determine the direction p.
      F77NAME(dogleg)(n,r.addr(),lr,diag.addr(),qtf.addr(),delta,
        wa1.addr(),wa2.addr(),wa3.addr());
//    store the direction p and x + p. calculate the norm of p.
      for (int j=0;j<n;j++) {
        wa1[j]=-wa1[j];
        wa2[j]=x[j]+wa1[j];
        wa3[j]=diag[j]*wa1[j];
      }
      double pnorm=wa3.l2Norm();
//    on the first iteration, adjust the initial step bound.
      if (iter==1) delta=min(delta,pnorm);
//    evaluate the function at x + p and calculate its norm.
      bool f_ok=(nonlinear_system.*function)(wa2,wa4);
      nfev++;
      double fnorm1=wa4.l2Norm();
//    compute the scaled actual reduction.
      double actred=(fnorm1<f_norm ? 1.-pow(fnorm1/f_norm,2) : -1.);
//    compute the scaled predicted reduction.
      int l=0;
      for (int i=0;i<n;i++) {
        double sum=0.;
        for (int j=i;j<n;j++,l++) sum+=r[l]*wa1[j];
        wa3[i]=qtf[i]+sum;
      }
      double temp=wa3.l2Norm();
      double prered=0.;
      if (temp<f_norm) prered=1.-pow(temp/f_norm,2);
//    compute the ratio of the actual to the predicted reduction.
      double ratio=(prered>0. ? actred/prered : 0.);
//    update the step bound.
      if (!f_ok || ratio<0.1) {
        ncsuc=0;
        ncfail++;
        delta=0.5*delta;
      } else {
        ncfail=0;
        ncsuc++;
        if (ratio>=0.5 || ncsuc>1) delta=max(delta,2.*pnorm);
        if (abs(ratio-1.)<=0.1) delta=2.*pnorm;
      }
//    test for successful iteration.
      if (f_ok && ratio>=1.e-4) {
//      successful iteration. update x, f, and their norms.
        for (int j=0;j<n;j++) {
          x[j]=wa2[j];
          wa2[j]=diag[j]*x[j];
          f[j]=wa4[j];
        }
        x_norm=wa2.l2Norm();
        f_norm=fnorm1;
        iter++;
      }
//    determine the progress of the iteration.
      nslow1++;
      if (f_ok && actred>=1.e-3) nslow1=0;
      if (jeval) nslow2++;
      if (f_ok && actred>=0.1) nslow2=0;
//    test for convergence.
      if (f_ok && (delta<=xtol*x_norm || f_norm<=0.)) { 
        converged=true;
        break;
      }
//    tests for termination and stringent tolerances.
      if (nfev>=maxfev || 0.1*max(0.1*delta,pnorm)<=DBL_EPSILON*x_norm ||
      nslow2==max_slow2 || nslow1==max_slow1) {
        failed=true;
        break;
      }
//    criterion for recalculating jacobian approximation
//    by forward differences.
      if (ncfail==max_fail) break;
//    calculate the rank one modification to the jacobian
//    and update qtf if necessary.
      for (int j=0;j<n;j++) {
        double sum=F77NAME(ddot)(n,&fjac(0,j),1,wa4.addr(),1);
        wa2[j]=(sum-wa3[j])/pnorm;
        wa1[j]=diag[j]*((diag[j]*wa1[j])/pnorm);
        if (ratio>=1.e-4) qtf[j]=sum;
      }
//    compute the qr factorization of the updated jacobian.
      F77NAME(r1updt)(n,n,r.addr(),lr,wa1.addr(),wa2.addr(),wa3.addr(),
        sing);
      F77NAME(r1mpyq)(n,n,fjac.addr(),ldfjac,wa2.addr(),wa3.addr());
      F77NAME(r1mpyq)(1,n,qtf.addr(),1,wa2.addr(),wa3.addr());
//    end of the inner loop.
      jeval=false;
    } while (true);
//  end of the outer loop.
  } while (!converged && !failed);
//termination, either normal or user imposed.
  CHECK_TEST(converged);
}
