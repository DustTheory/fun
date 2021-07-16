#ifndef FORT2D_H
#define FORT2D_H "$Header:$"
#include "Arch.H"

extern "C" {
  void F77_NAME(block_jacobi_preconditioner)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double *matrix,const double *cg_residual, 
    double *diagonal_block,double *cg_preconditioner_solve);
  void F77_NAME(compute_residual)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    double &residual_max,const double *matrix,const double *rhs,
    const double *solution,double *residual);
  void F77_NAME(copy_dlap_solution)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double *dlap_, double *x);
  double F77NAME(d1mach)(const int &i);
  void F77NAME(dgesvd)(const char &jobu,const char &jobvt,const int &m,
    const int &n,double *a,const int &lda,double *s,double *u,
    const int &ldu,double *vt,const int &ldvt,double *work,
    const int &lwork,int &info);
//void F77NAME(dicf)(const int &n,const int &nnz,double *a,double *diag,
//  int *col_ptr,int *row_ind,const int &nbands,int &info,int *indr,
//  int *ndf,int *ist,double *w);
  void F77_NAME(dicf_preconditioner)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const int &n,const int &nnz, const double *a,const double *diag,
    const int *col_ptr,const double *cg_residual,const int *row_ind, 
    double *cg_preconditioner_solve,double *residual);
  double F77NAME(ddot)(const int &n,const double *dx,const int &incx,
    const double *dy,const int &incy);
//void F77NAME(dpotrf)(const char *uplo,const int &n,double *a,
//  const int &n,int &info);
  void F77NAME(dgmres)(const int &n,const double *b,double *x,
    const int &nelt,const int *ia,const int *ja,const double *a,
    const int &isym,
    void (*matvec)(const int &n,const double *x,double *y,const int &nelt,
      const int *ia,const int *ja,const double *a,const int &isym), 
    void (*msolve)(const int &n,const double *r,double *z,const int &nelt,
      const int *ia,const int *ja,const double *a,const int &isym,
      double *rwork,int *iwork), 
    const int &itol,const double &tol,const int &itmax,
    int &iter,double &err,int &ierr,const int &iunit,double *sb,
    double *sx,double *rgwk,const int &lrgw,int *igwk,const int &ligw,
    double *rwork,int *iwork);
  void F77NAME(domn)(const int &n,const double *b,double *x,
    const int &nelt,const int *ia,const int *ja,const double *a,
    const int &isym,
    void (*matvec)(const int &n,const double *x,double *y,const int &nelt,
      const int *ia,const int *ja,const double *a,const int &isym), 
    void (*msolve)(const int &n,const double *r,double *z,const int &nelt,
      const int *ia,const int *ja,const double *a,const int &isym,
      double *rwork,int *iwork), 
    const int &nsave,const int &itol,const double &tol,const int &itmax,
    int &iter,double &err,int &ierr,const int &iunit,double *r,
    double *z,double *p,double *ap,double *emap,double *dz,
    double *csav,double *rwork,int *iwork);
  void F77_NAME(gauss_seidel_red_black)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const int &iteration,double &residual_max,const double *matrix,
    const double *rhs,double *solution);
  void F77_NAME(gauss_seidel_to_fro)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const int &iteration,double &residual_max,const double *matrix,
    const double *rhs,double *solution);
  int F77NAME(istdic)(const int &n,double *diag,double *a,
    const int *ia,const int *ja,double *ta,int *ifirst,int *list);
  void F77_NAME(iterative_improvement_dicf)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const int &n,const int &nnz,double &residual_max,const double *a,
    const double *diag,const int *col_ptr,const double *matrix,
    const double *rhs,const int *row_ind,double *solution,
    double *residual); 
  void F77_NAME(iterative_improvement_dicf_setup)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double *matrix, 
    int &n,int &nnz,double *a,double *diag,int *col_ptr,int *row_ind);
  void F77_NAME(jacobi_omega)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double &omega,
    double &residual_max,const double *matrix,const double *rhs, 
    double *solution, double *residual);
  void F77_NAME(jacobi_preconditioner)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double *matrix,
    const double *cg_residual, double *cg_preconditioner_solve);
  void F77_NAME(matrix_multiply)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double *a,const double *x,
    double *ax);
  void F77_NAME(matrix_setup)(
    const int &fc_0,const int &fc_1,const int &lc_0,const int &lc_1,
    const int &start_0,const int &start_1,
    const double *first_residual,double *matrix);
  void F77_NAME(prolong_am)(
    const int &fc_0,const int &fc_1,const int &lc_0,const int &lc_1,
    const int &ff_0,const int &ff_1,const int &lf_0,const int &lf_1,
    const double *prolongcell,
    const double *prolongside0,const double *prolongside1,
    const double *vectorc,double *vectorf);
  void F77_NAME(prolong_fem)(
    const int &fc_0,const int &fc_1,const int &lc_0,const int &lc_1,
    const int &ff_0,const int &ff_1,const int &lf_0,const int &lf_1,
    const double *vectorc,double *vectorf);
  void F77NAME(rbugcell)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double *array);
  void F77NAME(rcmpcell)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double *array,
    double &arraymax,double &arraymin);
  void F77_NAME(restrict_am)(
    const int &fc_0,const int &fc_1,const int &lc_0,const int &lc_1,
    const int &ff_0,const int &ff_1,const int &lf_0,const int &lf_1,
    const double *prolongcell,
    const double *prolongside0,const double *prolongside1,
    const double *vectorf,double *vectorc);
  void F77_NAME(restrict_fem)(
    const int &fc_0,const int &fc_1,const int &lc_0,const int &lc_1,
    const int &ff_0,const int &ff_1,const int &lf_0,const int &lf_1,
    const double *vectorf,double *vectorc);
  void F77NAME(richardson)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double &mu,
    const double &residual_max,const double *matrix,const double *rhs, 
    double *solution, double *residual);
  void F77_NAME(setup_am)(
    const int &fc_0,const int &fc_1,const int &lc_0,const int &lc_1,
    const int &ff_0,const int &ff_1,const int &lf_0,const int &lf_1,
    const double *matrix,
    double *prolongcell,double *prolongside0,double *prolongside1);
  void F77_NAME(setup_convection_diffusion)(
    const int &fc_0,const int &fc_1,const int &lc_0,const int &lc_1,
    const int &ff_0,const int &ff_1,const int &lf_0,const int &lf_1,
    const double &decay,
    const double &peclet, const double *true_solution, 
    double *matrix,double *solution);
  void F77_NAME(setup_dlap_matrix)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    int &nelt,
    const double *matrix,double *dlap_a,int *dlap_ia,int *dlap_ja);
  void F77_NAME(setup_dlap_system)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
    const double *x, double *dlap_x);
  void F77_NAME(setup_heat)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double &decay,
    const double *diffusion,const double *true_solution,
    double *matrix,double *solution);
  void F77_NAME(setup_laplace)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double *diffusion,const double *true_solution,
    double *matrix,double *solution);
  void F77NAME(sor)(
    const int &fi_0,const int &fi_1,const int &la_0,const int &la_1,
    const int &ifirst_0,const int &ifirst_1,const int &ilast_0,
      const int &ilast_1,
    const double &omega,
    double &residual_max,const double *matrix,const double *rhs, 
    double *solution);
  void F77_NAME(vector_setup)(
    const int &fc_0,const int &fc_1,const int &lc_0,const int &lc_1,
    const int &start_0,const int &start_1,
    double *solution_increment);
}
#endif
