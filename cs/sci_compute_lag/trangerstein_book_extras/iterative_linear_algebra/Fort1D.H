#ifndef FORT1D_H
#define FORT1D_H "$Header:$"
#include <iostream>
#include "Arch.H"

extern "C" {
  void F77_NAME(algebraic_multigrid_prolongation)(const int &fic,
    const int &lac,const int &fif,const int &laf,const int &ifirstc,
    const int &ilastc,const double *matrixf,const double *solutionc,
    double *solutionf);
  void F77_NAME(algebraic_multigrid_restriction)(const int &fic,
    const int &lac,const int &fif,const int &laf,const int &ifirstc,
    const int &ilastc,const double *matrixf,const double *residualf,
    double *residualc);
  void F77_NAME(compute_residual)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double *a,const double *b,
    const double *x,double *ax_b);
  double F77NAME(d1mach)(const int &i);
  void F77NAME(daxpy)(const int &n,const double &da,const double *dx,
    const int &incx,double *dy,const int &incy);
  void F77NAME(dcopy)(const int &n,const double *dx,const int &incx,
    double *dy,const int &incy);
  double F77NAME(ddot)(const int &n,const double *dx,const int &incx,
    const double *dy,const int &incy);
  void F77NAME(dgbtrf)(const int &m,const int &n,const int &kl,
    const int &ku,double *ab,const int &ldab,int *ipiv,int &info);
  void F77NAME(dgbtrs)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const double *ab,const int &ldab,
    const int *ipiv,double *b,const int &ldb,int &info);
  void F77NAME(dgmres)(const int &n,const double *b,double *x,
    const int &nelt,const int *ia,const int *ja,const double *a,
    const int &isym,
    void (*matvec)(const int &n,const double *x,double *y,const int &nelt,const int *ia,const int *ja,const double *a,const int &isym), 
    void (*msolve)(const int &n,const double *r,double *z,const int &nelt,const int *ia,const int *ja,const double *a,const int &isym,double *rwork,int *iwork), 
    const int &itol,const double &tol,const int &itmax,
    int &iter,double &err,int &ierr,const int &iunit,double *sb,
    double *sx,double *rgwk,const int &lrgw,int *igwk,const int &ligw,
    double *rwork,int *iwork);
  void F77NAME(dgtsv)(const int &n,const int &nrhs,double *dl,
    double *d,double *du,double *b,const int &ldb,int &info);
  void F77NAME(dhels)(const double *a,const int &lda,const int &n,
    const double *q,double *b);
  void F77NAME(dheqr)(double *a,const int &lda,const int &n,double *q,
    int &info,int &ijob);
  double F77NAME(dnrm2)(const int &n,const double *dx,const int &incx);
  void F77NAME(dorth)(double *vnew,const double *v,double *hes,
    const int &n,const int &ll,const int &ldhes,const int &kmp,
    double &snormw);
  void F77NAME(dptsv)(const int &n,const int &nrhs,double *d,double *e,
    double *b,const int &ldb,int &info);
  void F77NAME(dscal)(const int &n,const double &da,const double *dx,
    const int &incx);
  void F77_NAME(finite_element_prolongation)(const int &fic,
    const int &lac,const int &fif,const int &laf,const int &ifirstc,
    const int &ilastc,const double *solutionc,double *solutionf);
  void F77_NAME(finite_element_restriction)(const int &fic,
    const int &lac,const int &fif,const int &laf,const int &ifirstc,
    const int &ilastc,const double *residualf,double *restrictionc);
  void F77_NAME(gauss_seidel)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double *matrix,
    const double *rhs, double *solution);
  void F77_NAME(gauss_seidel_red_black)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const int &iteration,
    const double *matrix,const double *rhs, double *solution);
  void F77_NAME(gauss_seidel_reverse)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double *matrix,
    const double *rhs, double *solution);
  void F77_NAME(gauss_seidel_to_fro)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const int &iteration,
    const double *matrix,const double *rhs, double *solution);
//void F77NAME(jacobi)(const int &fi,const int &la,
//  const int &ifirst,const int &ilast,
//  double *matrix,double *rhs, double *solution, double *residual);
  int F77NAME(idamax)(const int &n,const double *dx,const int &incx);
  void F77_NAME(jacobi_omega)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double &omega,
    const double *matrix,const double *rhs, double *solution, 
    double *residual);
  void F77_NAME(jacobi_preconditioner)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double *matrix,
    const double *rhs, double *solution);
  void F77_NAME(matrix_multiply)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double *a,const double *x,
    double *ax);
  void F77NAME(domn)(const int &n,const double *b,double *x,
    const int &nelt,const int *ia,const int *ja,const double *a,
    const int &isym,
    void (*matvec)(const int &n,const double *x,double *y,const int &nelt,const int *ia,const int *ja,const double *a,const int &isym), 
    void (*msolve)(const int &n,const double *r,double *z,const int &nelt,const int *ia,const int *ja,const double *a,const int &isym,double *rwork,int *iwork), 
    const int &nsave,const int &itol,const double &tol,const int &itmax,
    int &iter,double &err,int &ierr,const int &iunit,double *r,
    double *z,double *p,double *ap,double *emap,double *dz,
    double *csav,double *rwork,int *iwork);
  void F77NAME(richardson)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double &mu,
    const double *matrix,const double *rhs, 
    double *solution, double *residual);
  void F77NAME(sor)(const int &fi,const int &la,
    const int &ifirst,const int &ilast, const double &omega,
    const double *matrix,const double *rhs, double *solution);
  void F77NAME(print_loc)(void *p);
}
#endif
