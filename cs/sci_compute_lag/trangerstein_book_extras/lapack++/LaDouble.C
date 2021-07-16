#include <complex>
#include <float.h>
#include <limits>
#include <math.h>
#include <string.h>
#include "Vector.H"
const double double_half_=0.5;
const double double_mone_=-1.;
const double double_one_=1.;
const double double_undefined_=numeric_limits<double>::infinity();
const double double_zero_=0.;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  int F77NAME(idamax)(const int &n,const double *x,const int &incx);
  int F77NAME(idmax)(const int &n,const double *a,const int &inca);
  int F77NAME(idmin)(const int &n,const double *a,const int &inca);
  int F77NAME(idsumn)(const int &n,const double *a,const int &inca);
  int F77NAME(idsump)(const int &n,const double *a,const int &inca);
  double F77NAME(dasum)(const int&,const double*,const int&);
  void F77NAME(daxpy)(const int &n,const double &a,const double *x,
    const int &incx,double *y,const int &incy);
  void F77NAME(dcopy)(const int &n,const double *sx,const int &incx,
    double *sy,const int &incy);
  double F77NAME(ddot)(const int &n,const double *sx,const int &incx,
    const double *sy,const int &incy);
  void F77NAME(dgbcon)(const char &norm,const int &n,const int &kl,
    const int &ku,const double *AB,const int &ldab,const int *ipiv,
    const double &anorm,double &rcond,double *work,int *iwork,int &info);
  void F77NAME(dgbconnp)(const char &norm,const int &n,const int &kl,
    const int &ku,const double *AB,const int &ldab,const double &anorm,
    double &rcond,double *work,int *iwork,int &info);
  void F77NAME(dgbequ)(const int &m,const int &n,const int &kl,
    const int &ku,const double *AB,const int &ldab,double *r,double *c,
    double &rowcnd,double &colcnd,double &amax,int &info);
  void F77NAME(dgbmv)(const char &trans,const int &m,const int &n,
    const int &kl,const int &ku,const double &alpha,const double *A,
    const int &lda,const double *x,const int &incx,const double &beta,
    double *y,const int &incy);
  void F77NAME(dgbtf2)(const int &m,const int &n,const int &kl,
    const int &ku,double *AB,const int &ldab,int *ipiv,int &info);
  void F77NAME(dgbtf2np)(const int &m,const int &n,const int &kl,
    const int &ku,double *AB,const int &ldab,int &info);
  void F77NAME(dgbtrf)(const int &m,const int &n,const int &kl,
    const int &ku,double *AB,const int &ldab,int *ipiv,int &info);
  void F77NAME(dgbtrs)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const double *AB,const int &ldab,
    const int *ipiv,double *B,const int &ldb,int &info);
  void F77NAME(dgbtrsnp)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const double *AB,const int &ldab,
    double *B,const int &ldb,int &info);
  void F77NAME(dgecon)(const char &trans,const int &n,const double *a,
    const int &lda,const double &anorm,
    double &rcond,double *work,int *iwork,int &info);
  double F77NAME(dgeequ)(const int &m,const int &n,const double *A,
    const int &lda,double *r,double *c,double &rowcnd,double &rolcnd,
    double &amax,int &info);
  void F77NAME(dgeev)(const char &jobvl,const char &jobvr,const int &n,
    double *A,const int &lda,double *wr,double *wi,double *vl,
    const int &ldvl,double *vr,const int &ldvr,double *work,
    const int &lwork,int &info);
  void F77NAME(dgels)(const char &trans,const int &m,const int &n,
    const int &nrhs,double *a,const int &lda,double *b,const int &ldb,
    double *work,const int &lwork,const int &info);
  void F77NAME(dgemm)(const char &transa,const char &transb,
    const int &m,const int &n,const int &k, 
    const double &alpha,const double *A,
    const int &lda,const double *B,const int &ldb,
    const double &beta,double *C,const int &ldc);
  void F77NAME(dgemv)(const char &trans,const int &m,const int &n,
    const double &alpha,const double *a,const int &lda,const double *x,
    const int &incx,const double &beta,double *y,const int &incy);
  void F77NAME(dgelqf)(const int &m,const int &n,double *A,const int &lda,
    double *tau,double *work,int &lwork,int &info);
  void F77NAME(dgeqp3)(const int &m,const int &n,double *A,const int &lda,
    int *jpvt,double *tau,double *work,const int &lwork,int &info);
  void F77NAME(dgeqrf)(const int &m,const int &n,double *A,const int &lda,
    double *tau,double *work,int &lwork,int &info);
  void F77NAME(dger)(const int&,const int&,const double&,const double*,
    const int&,const double*,const int&,double*,const int&);
  void F77NAME(dgesv)(int &n,int &nrhs,double *a,int &lda,int *ipiv, 
    double *b, int &ldb, int &info);
  void F77NAME(dgesvd)(const char &jobu,const char &jobvt,const int &m, 
    const int &n,const double *A,const int &lda,
    double *S,double *U,const int &ldu,double *Vt,
    const int &ldvt,double *work,const int &lwork, 
    int &info);
  void F77NAME(dgesc2)(const int &n,const double *A,const int &lda,
    double *rhs,const int *ipiv,const int *jpiv,double &scale);
  void F77NAME(dgetc2)(const int &n,double *A,const int &lda,
    int *ipiv,int *jpiv,int &info);
  void F77NAME(dgetrf)(int &m,int &n,double *a,int &lda,int *ipiv,
    int &info);
  void F77NAME(dgetrs)(const char &trans,const int &n,const int &nrhs,
    double *a,const int &lda,int *ipiv,double *b,const int &ldb,
    int &info);
  void F77NAME(dgetri)(int &n,double *a,int &lda,int *ipiv,double *work,
    int &lwork,int &info);
  void F77NAME(dgtcon)(const char &norm,const int &n,const double *dl,
    const double *d,const double *du,const double *du2,const int *ipiv,
    const double &anorm,double &rcond,double *work,int *iwork,int &info);
  void F77NAME(dgtconnp)(const char &norm,const int &n,const double *dl,
    const double *d,const double *du,const double &anorm,double &rcond,
    double *work,int *iwork,int &info);
  void F77NAME(dgtmv)(const int &n,const double &alpha,const double *dl,
    const double *d,const double *du,const double *x,const int &incx,
    const double &beta,double *y,const int &incy);
  void F77NAME(dgtrfs)(const char &trans,const int &n,const int &nrhs,
    const double *dl,const double *d,const double *du,const double *dlf,
    const double *df,const double *duf,const double *du2,const int *ipiv,
    const double *B,const int &ldb,double *X,const int &ldx,double *ferr,
    double *berr,double *work,int *iwork,int &info);
  void F77NAME(dgtrfsnp)(const char &trans,const int &n,const int &nrhs,
    const double *dl,const double *d,const double *du,const double *dlf,
    const double *df,const double *duf,const double *B,const int &ldb,
    double *X,const int &ldx,double *ferr,double *berr,double *work,
    int *iwork,int &info);
  void F77NAME(dgtrfsr)(const char &trans,const int &n,const int &nrhs,
    const double *dl,const double *d,const double *du,const double *dlf,
    const double *df,const double *duf,const double *du2,const int *ipiv,
    const double *B,const int &ldb,double *X,const int &ldx,double *ferr,
    double *berr,double *work,int *iwork,int &info);
  void F77NAME(dgtrfsrnp)(const char &trans,const int &n,const int &nrhs,
    const double *dl,const double *d,const double *du,const double *dlf,
    const double *df,const double *duf,const double *B,const int &ldb,
    double *X,const int &ldx,double *ferr,double *berr,double *work,
    int *iwork,int &info);
  void F77NAME(dgtsv)(const int &n,const int &nrhs,const double *dl,
    const double *d,const double *du,double *b,const int &ldb,int &info);
  void F77NAME(dgtsvnp)(const int &n,const int &nrhs,const double *dl,
    const double *d,const double *du,double *b,const int &ldb,int &info);
  void F77NAME(dgtsvr)(const int &n,const int &nrhs,const double *dl,
    const double *d,const double *du,double *b,const int &ldb,int &info);
  void F77NAME(dgtsvrnp)(const int &n,const int &nrhs,const double *dl,
    const double *d,const double *du,double *b,const int &ldb,int &info);
  void F77NAME(dgttrf)(const int &n,double *dl,double *d,double *du,
    double *du2,int *ipiv,int &info);
  void F77NAME(dgttrfnp)(const int &n,double *dl,double *d,double *du,
    int &info);
  void F77NAME(dgttrs)(const char &trans,const int &n,const int &nrhs,
    double *dl,double *d,double *du,double *du2,int *ipiv,
    double *B,const int &ldb,int &info);
  void F77NAME(dhseqr)(const char &job,const char &compz,const int &n,
    const int &ilo,const int &ihi,double *H,const int &ldh,double *wr,
    double *wi,double *Z,const int &ldz,double *work,const int &lwork,
    int &info);
  void F77NAME(dlabad)(double &small,double &large);
  void F77NAME(dlacn2)(const int &n,double *v,double *x,int *isgn,
    double &est,int &kase,int *isave);
  void F77NAME(dlacpy)(const char &uplo,const int &m,const int &n,
    const double *A,const int &lda,double *B,const int &ldb);
  void F77NAME(dladiv)(const double&,const double&,const double&,
    const double&,double&,double&);
  void F77NAME(dgbamv)(const char &trans,const int &m,const int &n,
    const int &kl,const int &ku,const double &alpha,const double *A,
    const int &lda,const double *x,const int &incx,const double &beta,
    double *y,const int &incy);
  void F77NAME(dlaic1)(const int &job,const int &j,const double *x,
    const double &sest,const double *w,const double &gamma,
    double &sestpr,double &s,double &c);
  double F77NAME(dlamch)(const char&);
  double F77NAME(dlangb)(const char &norm,const int &n,const int &kl,
    const int &lu,const double *AB,const int &ldab,double *work);
  double F77NAME(dlange)(const char &norm,const int &m,const int &n,
    const double *A,const int &lda,double *work);
  double F77NAME(dlangt)(const char &norm,const int &n,const double *dl,
    const double *d,const double *du);
  double F77NAME(dlanhs)(const char &norm,const int &n,const double *A,
    const int &lda,double *work);
  double F77NAME(dlansb)(const char &norm,const char &uplo,const int &n,
    const int &k,const double *AB,const int &ldab,double *work);
  double F77NAME(dlanst)(const char &norm,const int &n,const double *d,
    const double *e);
  void F77NAME(dlapmt)(const bool &forwrd,const int &m,const int &n,
    double *X,const int &ldx,int *k);
  void F77NAME(dlaqgb)(const int &m,const int &n,const int &kl,
    const int &ku,double *AB,const int &ldab,const double *r,
    const double *c,const double &rowcnd,const double &colcnd,
    const double &amax,char &equed);
  void F77NAME(dlaqsb)(const char &uplo,const int &n,const int &kd,
    double *AB,const int &ldab,const double *s,const double &scond,
    const double &amax,char &equed);
  void F77NAME(dlaqsy)(const char &uplo,const int &n,double *A,
    const int &lda,const double *s,const double &scond,
    const double &amax,char &equed);
  void F77NAME(dlaswp)(const int &n,double *A,const int &lda,
    const int &k1,const int &k2,const int *ipiv,const int &incx);
  double F77NAME(dlansy)(const char &norm,const char &uplo,const int &n,
    const double *A,const int &lda,double *work);
  double F77NAME(dlantb)(const char &norm,const char &uplo,
    const char &diag,const int &n,const int &k,const double *AB,
    const int &ldab,double *work);
  double F77NAME(dlantr)(const char &norm,const char &uplo,
    const char &diag,const int &m,const int &n,const double *A,
    const int &lda,double *work);
  void F77NAME(dlaqge)(const int &M,const int &n,double *A,const int &lda,
    const double *r,const double *c,const double &rowcnd,
    const double &colcnd,const double &amax,char &equed);
  void F77NAME(dlargv)(const int &n,double *x,const int &incx,
    double *y,const int &incy,double *c,const int &incc);
  void F77NAME(dlartv)(const int &n,double *x,const int &incx,
    double *y,const int &incy,const double *c,const double *s,
    const int &incc);
  void F77NAME(dlascl)(const char &type,const int &kl,const int &ku,
    const double &cfrom,const double &cto,const int &m,const int &n,
    double *A,const int &lda,int &info);
  void F77NAME(dlaset)(const char &uplo,const int &m,const int &n,
    const double &alpha,const double &beta,double *A,const int &lda);
  void F77_NAME(dla_geamv)(const int &trans,const int &m,const int &n,
    const double &alpha,const double *a,const int &lda,const double *x,
    const int &incx,const double &beta,double *y,const int &incy);
  double F77_NAME(dla_porpvgrw)(const char &uplo,const int &ncols,
    const double *A,const int &lda,const double *AF,const int &ldaf,
    double *work);
  void F77_NAME(dla_syamv)(const int &uplo,const int &n,
    const double &alpha,const double *A,const int &lda,
    const double *x,const int &ldx,const double &beta,
    double *y,const int &incy);
  double F77_NAME(dla_syrpvgrw)(const char &uplo,const int &n,
    const int &info,const double *A,const int &lda,const double *AF,
    const int &ldaf,const int *ipiv,double *work);
  double F77NAME(dnrm2)(const int &n,double *x,const int &incx);
  void F77NAME(dorglq)(int &m,int &n,int &k,double *a,int &lda,
    const double *tau,double *work,int &lwork,
    int &info);
  void F77NAME(dorgqr)(const int &m,const int &n,const int &k,double *A,
    const int &lda,const double *tau,double *work,const int &lwork,
    int &info);
  void F77NAME(dormlq)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,double *A,const int &lda,double *tau,
    double *C,const int &ldc,double *work,int &lwork,int &info);
  void F77NAME(dormqr)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,double *A,const int &lda,double *tau,
    double *C,const int &ldc,double *work,int &lwork,int &info);
  void F77NAME(dormrz)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,const int &l,const double *A,const int &lda,
    const double *tau,double *C,const int &ldc,double *work,
    const int &lwork,int &info);
  void F77NAME(dpbcon)(const char &uplo,const int &n,const int &kd,
    double *AB,const int &ldab,const double &anorm,double &rcond,
    double *work,int *iwork,int &info);
  void F77NAME(dpbequ)(const char &uplo,const int &n,const int &kd,
    const double *AB,const int &ldab,double *s,double &scond,
    double &amax,int &info);
  void F77NAME(dpbsv)(const char &uplo,const int &n,const int &kd,
    const int &nrhs,double *AB,const int &ldab,double *B,const int &ldb,
    int &info);
  void F77NAME(dpbtrf)(const char &uplo,const int &n,const int &kd,
    double *AB,const int &ldab,int &info);
  void F77NAME(dpbtrs)(const char &uplo,const int &n,const int &kd,
    const int &nrhs,const double *AB,const int &ldab,double *B,
    const int &ldb,int &info);
  void F77NAME(dpocon)(const char &uplo,const int &n,const double *a,
    const int &lda,const double &anorm,double &rcond,double *work,
    int *iwork,int &info);
//void F77NAME(dpoequ)(const int &n,const double *A,const int &lda,
//  double *s,double &scond,double &amax,int &info);
  void F77NAME(dpoequb)(const int &n,const double *A,const int &lda,
    double *s,double &scond,double &amax,int &info);
  void F77NAME(dpotrf)(const char &uplo,const int &n,double *A,
    const int &lda,int &info);
  void F77NAME(dpotri)(const char &uplo,const int &n,double *A,
    const int &lda,int &info);
  void F77NAME(dpotrs)(const char &uplo,const int &n,const int &nrhs,
    const double *A,const int &lda,double *B,const int &ldb,int &info);
  void F77NAME(dptcon)(const int &n,const double *d,const double *e,
    const double &anorm,double &rcond,double *work,int &info);
  void F77NAME(dptsv)(const int &n,const int &nrhs,double *d,double *e,
    double *b,const int &ldb,int &info);
  void F77NAME(dptrfs)(const int &n,const int &nrhs,const double *d,
    const double *e,const double *df,const double *ef,const double *B,
    const int &ldb,double *X,const int &ldx,double *ferr,double *berr,
    double *work,int &info);
  void F77NAME(dptrfsr)(const int &n,const int &nrhs,const double *d,
    const double *e,const double *df,const double *ef,const double *B,
    const int &ldb,double *X,const int &ldx,double *ferr,double *berr,
    double *work,int &info);
  void F77NAME(dpttrf)(const int &n,double *d,double *e,int &info);
  void F77NAME(dpttrs)(const int &n,const int &nrhs,const double *d,
    const double *e, double *B,const int &ldb,int &info);
  void F77NAME(drot)(const int &n,double *sx,const int &incx,double *sy,
    const int &incy,const double &c,const double &s);
  void F77NAME(drotg)(double &sa,double &sb,double &c,double &s);
  void F77NAME(drscl)(const int &n,const double &sa,double *sx,
    const int &incx);
  void F77NAME(dsbev)(const char &jobz,const char &uplo,const int &n,
    const int &kd,double *AB,const int &ldab,double *w,double *Z,
    const int &ldz,double *work,int &info);
  void F77NAME(dsbamv)(const char &uplo,const int &n,const int &k,
    const double &alpha,const double *A,const int &lda,const double *x,
    const int &incx,const double &beta,double *y,const int &incy);
  void F77NAME(dsbmv)(const char &uplo,const int &n,const int &k,
    const double &alpha,const double *A,const int &lda,const double *x,
    const int &incx,const double &beta,double *y,const int &incy);
  void F77NAME(dscal)(const int &n,const double &a,double *x,
    const int &incx);
  void F77NAME(dsteqr)(const char &jobz,const int &n,double *d,double *e,
    double *Z,const int &ldz,double *work,int &info);
  void F77NAME(dstmv)(const int &n,const double &alpha,const double *dl,
    const double *d,const double *x,const int &incx,const double &beta,
    double *y,const int &incy);
  void F77NAME(dswap)(const int &n,double *sx,const int &incx,double *sy,
    const int &incy);
  void F77NAME(dsycon)(const char &uplo,const int &n,const double *a,
    const int &lda,const int *ipiv,const double &anorm,double &rcond,
    double *work,int *iwork,int &info);
  void F77NAME(dsyequb)(const char &uplo,const int &n,const double *A,
    const int &lda,double *s,double &scond,double &amax,double *work,
    int &info);
  void F77NAME(dsyev)(const char &jobz,const char &uplo,const int &n,
    double *A,const int &lda,double *W,double *work,int &lwork,
    int &info);
  void F77NAME(dsymm)(const char &side,const char &uplo,const int &m,
    const int &n,const double &alpha,const double *A,const int &lda,
    const double *B,const int &ldb,const double &beta,double *C,
    const int &ldc);
  void F77NAME(dsymv)(const char &uplo,const int &n,const double &alpha,
    const double *A,const int &lda,const double *x,const int &incx,
    const double &beta,double *y,const int &incy);
  void F77NAME(dsyr)(const char &uplo,const int &n,const double &alpha,
    const double *x,const int &incx,double *A,const int &lda);
  void F77NAME(dsyrk)(const char &uplo,const char &trans,const int &n,
    const int &k,const double &alpha,const double *A,const int &lda,
    const double &beta,double *C,const int &ldc);
  void F77NAME(dsyr2k)(const char &uplo,const char &trans,const int &n,
    const int &k,const double &alpha,const double *A,const int &lda,
    const double *B,const int &ldb,const double &beta,double *C,
    const int &ldc);
  void F77NAME(dsyr2)(const char &uplo,const int &n,const double &alpha,
    const double *x,const int &incx,const double *y,const int &incy,
    double *A,const int &lda);
  void F77NAME(dsytrf)(const char &uplo,const int &n,double *a,
    const int &lda,int *ipiv,double *work,const int &lwork,int &info);
  void F77NAME(dsytri)(const char &uplo,const int &n,double *a,
    const int &lda,int *ipiv,double *work,int &info);
  void F77NAME(dsytrs)(const char &uplo,const int &n,const int &nrhs,
    double *a,const int &lda,int *ipiv,double *b,const int &ldb,
    int &info);
  void F77NAME(dtbsv)(const char &uplo,const char &trans,const char &diag,
    const int &n,const int &k,const double *A,const int &lda,double *x,
    const int &incx);
  void F77NAME(dtrcon)(const char &norm,const char &uplo,
    const char &diag,const int &n,const double *A,const int &lda,
    double &rcond,double *work,int *iwork,int &info);
  void F77NAME(dtrevc)(const char &side,const char &howmny,bool *select,
    const int &n,const double *T,const int &ldt,double *Vl,
    const int &ldvl,double *Vr,const int &ldvr,const int &mm,
    int &m,double *work,int &info);
  void F77NAME(dtrmv)(const char &uplo,const char &trans,
    const char &diag,const int &n,const double *A,const int &lda,
    double *x,const int &incx);
  void F77NAME(dtrsv)(const char&,const char&,const char&,const int &n,
    const double *A,const int &lda,double *x,const int &incx);
  void F77NAME(dtrmm)(const char &side,const char &uplo,
    const char &transa,const char &diag,const int &m,const int &n,
    const double &alpha,const double *A,const int &lda,double *B,
    const int &ldb);
  void F77NAME(dtrsm)(const char &side,const char &uplo,
    const char &transa,const char &diag,const int &m,const int &n,
    const double &alpha,const double *a,const int &lda,double *b,
    const int &ldb);
  void F77NAME(dtrtri)(const char &uplo,const char &diag,const int &n,
    double *A,const int &lda,int &info);
  void F77NAME(dtzrzf)(const int &m,const int &n,double *A,const int &lda,
    double *tau,double *work,const int &lwork,int &info);
  double F77NAME(dznrm2)(const int &n,const complex<double> *x,
    const int &incx);
  int F77NAME(iladlc)(const int &m,const int &n,const double *A,
    const int &lda);
  int F77NAME(iladlr)(const int &m,const int &n,const double *A,
    const int &lda);
  int F77NAME(ilaenv)(const int &ispec,const char *name,const char *opts,
    const int &n1,const int &n2,const int &n3,const int &n4);
  int F77NAME(ilatrans)(const char &trans);
  int F77NAME(ilauplo)(const char &trans);
  void F77NAME(zdscal)(const int &n,const double &da,
    complex<double> *zx,const int &incx);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Vector.C"

#ifdef DEBUG
template<> const double Vector<double,double>::undefined_ = 
  numeric_limits<double>::infinity();
#endif
template<> const double Vector<double,double>::zero_ = double_zero_;
template<> const double Vector<double,double>::one_ = double_one_;
template<> const double Vector<double,double>::mone_ = double_mone_;

template<> int Vector<double,double>::amax() const { 
  return F77NAME(idamax)(sz,data,1)-1;
}
template<> double Vector<double,double>::asum() const {
  return F77NAME(dasum)(sz,data,1);
}
template<> void Vector<double,double>::axpy(double a,
const Vector<double,double> &x) {
  F77NAME(daxpy)(min(sz,x.sz),a,x.data,1,data,1);
}
template<> double Vector<double,double>::dot(
const Vector<double,double> &x) const {
  return F77NAME(ddot)(min(sz,x.sz),x.data,1,data,1);
}
template<> double Vector<double,double>::dotc(
const Vector<double,double> &x) const {
  OBSOLETE(0);
}
template<> double Vector<double,double>::dotu(
const Vector<double,double> &x) const {
  OBSOLETE(0);
}
template<> double Vector<double,double>::nrm2() const {
  return F77NAME(dnrm2)(sz,data,1);
}
template<> void Vector<double,double>::rot(Vector<double,double> &x,
double c,double s) {
  F77NAME(drot)(min(sz,x.sz),x.data,1,data,1,c,s);
}
template<> void Vector<double,double>::scal(double a) {
  F77NAME(dscal)(sz,a,data,1);
}
template<> void Vector<double,double>::swap(Vector<double,double> &x) {
  F77NAME(dswap)(min(sz,x.sz),x.data,1,data,1);
}
template<> void Vector<double,double>::largv(Vector<double,double> &y,
Vector<double,double> &c) {
  int n=size();
  CHECK_SAME(n,y.size());
  CHECK_SAME(n,c.size());
  F77NAME(dlargv)(n,data,1,y.data,1,c.data,1);
}
template<> void Vector<double,double>::lartv(Vector<double,double> &y,
const Vector<double,double> &c,const Vector<double,double> &s) {
  int n=size();
  CHECK_SAME(n,y.size());
  CHECK_SAME(n,c.size());
  CHECK_SAME(n,s.size());
  F77NAME(dlartv)(n,addr(),1,y.addr(),1,c.addr(),s.addr(),1);
}

template class Vector<double,double>; 
template void testVector(const double&,const double&);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Matrix.C"

#ifdef DEBUG
template<> const double Matrix<double,double>::undefined_ = 
  numeric_limits<double>::infinity();
#endif
template<> const double Matrix<double,double>::zero_ = double_zero_;
template<> const double Matrix<double,double>::one_ = double_one_;
template<> const double Matrix<double,double>::mone_ = double_mone_;

/*
template<> Matrix<double,double>* Matrix<double,double>::transpose()
const {
  int m=size(0),n=size(1);
  Matrix<double,double> *X=OPERATOR_NEW Matrix<double,double>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(m,addr(0,j),1,X->addr(j,0),n);
  }
  return X;
}

template<> Matrix<double,double>*
Matrix<double,double>::conjugateTranspose() const {
  return transpose();
}
*/

template<> void Matrix<double,double>::interchangeColumns(int i,int j) {
  int m=size(0),n=size(1);
  CHECK_BOUNDS(i,0,n)
  CHECK_BOUNDS(j,0,n)
  F77NAME(dswap)(m,addr(0,i),1,addr(0,j),1);
}

template<> void Matrix<double,double>::interchangeRows(int i,int j) {
  int m=size(0),n=size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,m)
  F77NAME(dswap)(n,addr(i,0),m,addr(j,0),m);
}

template<> void Matrix<double,double>::gemv(double alpha,
const Vector<double,double> &x,double beta,Vector<double,double> &y,
char trans) const {
//TRACER_CALL(t,"Matrix::gemv");
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(x.size(),m);
    CHECK_SAME(y.size(),n);
    F77NAME(dgemv)('T',m,n,alpha,addr(),m,x.addr(),1,beta,y.addr(),1);
  } else {
    CHECK_SAME(x.size(),n);
    CHECK_SAME(y.size(),m);
    F77NAME(dgemv)('N',m,n,alpha,addr(),m,x.addr(),1,beta,y.addr(),1);
  }
}

template<> void Matrix<double,double>::ger(double alpha,
const Vector<double,double> &x,const Vector<double,double> &y) {
  int m=size(0),n=size(1);
  CHECK_SAME(x.size(),m);
  CHECK_SAME(y.size(),n);
  F77NAME(dger)(m,n,alpha,x.addr(),1,y.addr(),1,addr(),m);
}

template<> void Matrix<double,double>::gerc(double alpha,
const Vector<double,double> &x,const Vector<double,double> &y) {
  OBSOLETE("inappropriate for this class");
}

template<> void Matrix<double,double>::geru(double alpha,
const Vector<double,double> &x,const Vector<double,double> &y) {
  OBSOLETE("inappropriate for this class");
}

template<> void Matrix<double,double>::gemm(double alpha,
const Matrix<double,double> &A,const Matrix<double,double> &B,
double beta,char transa,char transb) {
  int m=size(0),n=size(1);
  int k=0;
  if (transa=='N' || transa=='n') {
    CHECK_SAME(A.size(0),m);
    k=A.size(1);
    if (transb=='N' || transb=='n') {
      CHECK_SAME(B.size(1),n);
      CHECK_SAME(B.size(0),k);
    } else {
      CHECK_SAME(B.size(0),n);
      CHECK_SAME(B.size(1),k);
    }
  } else {
    CHECK_SAME(A.size(1),m);
    k=A.size(0);
    if (transb=='N' || transb=='n') {
      CHECK_SAME(B.size(1),n);
      CHECK_SAME(B.size(0),k);
    } else {
      CHECK_SAME(B.size(0),n);
      CHECK_SAME(B.size(1),k);
    }
  }
  F77NAME(dgemm)(transa,transb,m,n,k,alpha,A.addr(),A.size(0),
    B.addr(),B.size(0),beta,addr(),m);
}

template<> double Matrix<double,double>::equilibrate(
Vector<double,double> &r,Vector<double,double> &c,double &rowcnd,
double &colcnd) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,r.size());
  CHECK_SAME(n,c.size());
  double amax=double_undefined_;
  int info;
  F77NAME(dgeequ)(m,n,addr(),m,r.addr(),c.addr(),rowcnd,colcnd,
    amax,info);
  CHECK_SAME(info,0);
  return amax;
}

template<> void Matrix<double,double>::copyFrom(char uplo,int m,int n,
const Matrix<double,double> &A) {
  int s0=size(0),as0=A.size(0);
  m=min(m,min(s0,as0));
  n=min(n,min(size(1),A.size(1)));
  if (uplo=='A' || uplo=='a') {
    F77NAME(dlacpy)(uplo,m,n,A.addr(),as0,addr(),s0);
//  otherwise, dlacpy only copies a triangular part:
  } else if (uplo=='L' || uplo=='l') {
    for (int j=0;j<n;j++) {
      if (j<m) {
        F77NAME(dcopy)(m-j,A.addr(j,j),1,addr(j,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+1,m),A.addr(0,j),1,addr(0,j),1);
    }
  }
}

template<> void Matrix<double,double>::scale(char type,int kl,int ku,
double denominator,double numerator) {
  int m=size(0),info;
  F77NAME(dlascl)(type,kl,ku,denominator,numerator,m,size(1),addr(),
    m,info);
  CHECK_SAME(info,0);
}

template<> void Matrix<double,double>::set(char uplo,double offdiag,
double diag) {
  int m=size(0);
  F77NAME(dlaset)(uplo,m,size(1),offdiag,diag,addr(),m);
}

template<> int Matrix<double,double>::lastNonzeroColumn() const {
  int m=size(0);
  return F77NAME(iladlc)(m,size(1),addr(),m);
}

template<> int Matrix<double,double>::lastNonzeroRow() const {
  int m=size(0);
  return F77NAME(iladlr)(m,size(1),addr(),m);
}

template<> double Matrix<double,double>::normFrobenius() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlange)('F',m,size(1),addr(),m,work);
}

template<> double Matrix<double,double>::normInfinity() const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double val=F77NAME(dlange)('I',m,size(1),addr(),m,work);
  delete [] work; work=0;
  return val;
}

template<> double Matrix<double,double>::normMaxEntry() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlange)('M',m,size(1),addr(),m,work);
}

template<> double Matrix<double,double>::normOne() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlange)('1',m,size(1),addr(),m,work);
}

template<> Matrix<double,double>* Matrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(M.size(0),k);
  Matrix<double,double> *P=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dgemm)('N','N',m,n,k,double_one_,addr(),m,M.addr(),k,
    double_zero_,P->addr(),m);
  return P;
}

template<> Vector<double,double>* Matrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(v.size(),n);
  Vector<double,double> *P=OPERATOR_NEW Vector<double,double>(m);
  F77NAME(dgemv)('N',m,n,double_one_,addr(),m,v.addr(),1,double_zero_,
    P->addr(),1);
  return P;
}

template<> void Matrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  int m=size(0), n=size(1);
  Matrix<double,double> AF(*this);
  int info;
  if (m==n) {
    CHECK_SAME(m,b.size())
    CHECK_SAME(n,x.size())
    x.copy(b);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(dgetrf)(m,m,AF.addr(),m,ipiv,info);
    if (info==0) {
      F77NAME(dgetrs)(trans,m,1,AF.addr(),m,ipiv,x.addr(),m,info);
    }
    CHECK_SAME(info,0)
    delete [] ipiv;
  } else {
    double w=numeric_limits<double>::infinity();
    int lwork=-1;
    bool transposed=(trans!='N' && trans!='n');
    if (m>n) {
      if (transposed) {
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())
        x.copyFrom(n,b);
        F77NAME(dgels)(trans,m,n,1,AF.addr(),m,x.addr(),m,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dgels)(trans,m,n,1,AF.addr(),m,x.addr(),m,work,
          lwork,info);
        CHECK_SAME(info,0)
        delete [] work;
      } else {
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        Vector<double,double> xtmp(m);
        xtmp.copy(b);
        F77NAME(dgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),m,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),m,work,
          lwork,info);
        CHECK_SAME(info,0)

        x.copyFrom(n,xtmp);
        delete [] work;
      }
    } else {
      if (transposed) {
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())
        Vector<double,double> xtmp(n);
        xtmp.copy(b);
        F77NAME(dgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),n,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),n,work,
          lwork,info);
        CHECK_SAME(info,0)

        x.copyFrom(m,xtmp);
        delete [] work;
      } else {
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        x.copyFrom(m,b);
        F77NAME(dgels)(trans,m,n,1,AF.addr(),m,x.addr(),n,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dgels)(trans,m,n,1,AF.addr(),m,x.addr(),n,work,
          lwork,info);
        CHECK_SAME(info,0)
        delete [] work;
      }
    }
  }
}

template<> void Matrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  int m=size(0), n=size(1);
  Matrix<double,double> AF(*this);
  int info;
  bool left_side=(side=='L' || side=='l');
  bool transposed=(trans!='N' && trans!='n');
  if (m==n) {
    X.copy(B);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(dgetrf)(m,m,AF.addr(),m,ipiv,info);
    if (info==0) {
      if (left_side) {
        int k=B.size(1);
        CHECK_SAME(k,X.size(1))
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        F77NAME(dgetrs)(trans,m,k,AF.addr(),m,ipiv,X.addr(),m,info);
      } else {
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        CHECK_SAME(m,B.size(1))
        CHECK_SAME(n,X.size(1))
        if (transposed) {
          for (int j=0;j<m-1;j++) {
            if (ipiv[j]-1!=j) {
              F77NAME(dswap)(k,X.addr(0,j),1,X.addr(0,ipiv[j]-1),1);
            }
          }
          F77NAME(dtrsm)('R','L','T','U',k,m,double_one_,AF.addr(),m,
            X.addr(),k);
          F77NAME(dtrsm)('R','U','T','N',k,m,double_one_,AF.addr(),m,
            X.addr(),k);
        } else {
          F77NAME(dtrsm)('R','U','N','N',k,m,double_one_,AF.addr(),m,
            X.addr(),k);
          F77NAME(dtrsm)('R','L','N','U',k,m,double_one_,AF.addr(),m,
            X.addr(),k);
          for (int j=m-2;j>=0;j--) {
            if (ipiv[j]-1!=j) {
              F77NAME(dswap)(k,X.addr(0,j),1,X.addr(0,ipiv[j]-1),1);
            }
          }
        }
      }
    }
    CHECK_SAME(info,0)
    delete [] ipiv; ipiv=0;
  } else {
    double w=numeric_limits<double>::infinity();
    int lwork=-1;
    if (m>n) {
      if (left_side) {
        int k=B.size(1);
        CHECK_SAME(k,X.size(1))
        if (transposed) {
          CHECK_SAME(n,B.size(0))
          CHECK_SAME(m,X.size(0))
          X=zero_;
          X.copyFrom('A',n,k,B);
          F77NAME(dgels)(trans,m,n,k,AF.addr(),m,X.addr(),m,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w);
          double *work=OPERATOR_NEW_BRACKET(double,lwork);
          F77NAME(dgels)(trans,m,n,k,AF.addr(),m,X.addr(),m,work,
            lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        } else {
          CHECK_SAME(m,B.size(0))
          CHECK_SAME(n,X.size(0))
          Matrix<double,double> Xtmp(B);
          F77NAME(dgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),m,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w);
          double *work=OPERATOR_NEW_BRACKET(double,lwork);
          F77NAME(dgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),m,work,
            lwork,info);
          CHECK_SAME(info,0)

          X.copyFrom('A',n,k,Xtmp);
          delete [] work; work=0;
        }
      } else { // m>n, right side
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        double *tau=OPERATOR_NEW_BRACKET(double,n);
        F77NAME(dgeqrf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dgeqrf)(m,n,AF.addr(),m,tau,work,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        if (transposed) {
          CHECK_SAME(m,B.size(1))
          CHECK_SAME(n,X.size(1))
          Matrix<double,double> Xtmp(B);

          lwork=-1;
          F77NAME(dormqr)('R','N',k,m,n,AF.addr(),m,tau,Xtmp.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w);
          work=OPERATOR_NEW_BRACKET(double,lwork);
          F77NAME(dormqr)('R','N',k,m,n,AF.addr(),m,tau,Xtmp.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
          X.copyFrom('A',k,n,Xtmp);

          F77NAME(dtrsm)('R','U','T','N',k,n,double_one_,AF.addr(),m,
            X.addr(),k);
        } else {
          CHECK_SAME(n,B.size(1))
          CHECK_SAME(m,X.size(1))
          X=zero_;
          X.copyFrom('A',k,n,B);
          F77NAME(dtrsm)('R','U','N','N',k,n,double_one_,AF.addr(),m,
            X.addr(),k);
          lwork=-1;
          F77NAME(dormqr)('R','T',k,m,n,AF.addr(),m,tau,X.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w);
          work=OPERATOR_NEW_BRACKET(double,lwork);
          F77NAME(dormqr)('R','T',k,m,n,AF.addr(),m,tau,X.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        }
        delete [] tau; tau=0;
      }
    } else { // m < n
      if (left_side) {
        int k=B.size(1);
        CHECK_SAME(k,X.size(1))
        if (transposed) {
          CHECK_SAME(n,B.size(0))
          CHECK_SAME(m,X.size(0))
          Matrix<double,double> Xtmp(B);
          F77NAME(dgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),n,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w);
          double *work=OPERATOR_NEW_BRACKET(double,lwork);
          F77NAME(dgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),n,work,
            lwork,info);
          CHECK_SAME(info,0)

          X.copyFrom('A',m,k,Xtmp);
          delete [] work; work=0;
        } else {
          CHECK_SAME(m,B.size(0))
          CHECK_SAME(n,X.size(0))
          X.copyFrom('A',m,k,B);
          F77NAME(dgels)(trans,m,n,k,AF.addr(),m,X.addr(),n,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w);
          double *work=OPERATOR_NEW_BRACKET(double,lwork);
          F77NAME(dgels)(trans,m,n,k,AF.addr(),m,X.addr(),n,work,
            lwork,info);
          CHECK_SAME(info,0)
          delete [] work;
        }
      } else { // m<n, right side
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        double *tau=OPERATOR_NEW_BRACKET(double,n);
        F77NAME(dgelqf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)
        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dgelqf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        if (transposed) {
          CHECK_SAME(m,B.size(1))
          CHECK_SAME(n,X.size(1))
          X=zero_;
          X.copyFrom('A',k,m,B);
          F77NAME(dtrsm)('R','L','T','N',k,m,double_one_,AF.addr(),m,
            X.addr(),k);
          lwork=-1;
          F77NAME(dormlq)('R','N',k,n,m,AF.addr(),m,tau,X.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w);
          work=OPERATOR_NEW_BRACKET(double,lwork);
          F77NAME(dormlq)('R','N',k,n,m,AF.addr(),m,tau,X.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        } else {
          CHECK_SAME(n,B.size(1))
          CHECK_SAME(m,X.size(1))
          Matrix<double,double> Xtmp(B);
          lwork=-1;
          F77NAME(dormlq)('R','T',k,n,m,AF.addr(),m,tau,Xtmp.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w);
          work=OPERATOR_NEW_BRACKET(double,lwork);
          F77NAME(dormlq)('R','T',k,n,m,AF.addr(),m,tau,Xtmp.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
          X.copyFrom('A',k,m,Xtmp);

          F77NAME(dtrsm)('R','L','N','N',k,m,double_one_,AF.addr(),m,
            X.addr(),k);
        }
        delete [] tau; tau=0;
      }
    }
  }
}

template class Matrix<double,double>; 
template void testMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SquareMatrix.C"

template<> SquareMatrix<double,double>*
SquareMatrix<double,double>::operator*(
const SquareMatrix<double,double> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  F77NAME(dgemm)('N','N',n,n,n,double_one_,addr(),n,S.addr(),n,
    double_zero_,P->addr(),n);
  return P;
}

template<> Matrix<double,double>*
SquareMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *P=(m==n ?
    OPERATOR_NEW SquareMatrix<double,double>(n) :
    OPERATOR_NEW Matrix<double,double>(m,n));
  F77NAME(dgemm)('N','N',m,n,m,double_one_,addr(),m,M.addr(),m,
    double_zero_,P->addr(),m);
  return P;
}

/*
template<> SquareMatrix<double,double>*
SquareMatrix<double,double>::transpose() const {
  int n=size(0);
  SquareMatrix<double,double> *X=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n,addr(0,j),1,X->addr(j,0),n);
  }
  return X;
}

template<> SquareMatrix<double,double>*
SquareMatrix<double,double>::conjugateTranspose() const {
  return transpose();
}
*/

template<> double 
SquareMatrix<double,double>::reciprocalConditionNumber(char norm) const {
  int n=size(0),info;
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,4*n);
  double anorm=F77NAME(dlange)(norm,n,n,addr(),n,work);

  SquareMatrix<double,double> *AF=
    OPERATOR_NEW SquareMatrix<double,double>(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(dgetrf)(n,n,AF->addr(0,0),n,ipiv,info);
  CHECK_SAME(info,0)
  delete [] ipiv; ipiv=0;

  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(dgecon)(norm,n,AF->addr(),n,anorm,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete[] work;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SquareMatrix<double,double>*
SquareMatrix<double,double>::inverse() const {
  int n=size(0);
  SquareMatrix<double,double> *Ainv=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  *Ainv = *this;
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(dgetrf)(n,n,Ainv->addr(0,0),n,ipiv,info);
  CHECK_SAME(info,0)

  double w=numeric_limits<double>::infinity();
  int lwork=-1;
  F77NAME(dgetri)(n,Ainv->addr(0,0),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dgetri)(n,Ainv->addr(0,0),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)

  delete [] ipiv;
  delete [] work;
  return Ainv;
}
*/

template<> Vector<double,complex<double> >*
SquareMatrix<double,double>::eigenvalues(
SquareMatrix<double,complex<double> > *&V,
SquareMatrix<double,complex<double> > *&U) const {
  int n=size(0);
  if (V!=0) CHECK_SAME(n,V->size(0));
  if (U!=0) CHECK_SAME(n,U->size(0));
  char jobvl=(V==0 ? 'N' : 'V');
  char jobvr=(U==0 ? 'N' : 'V');
  double *wr=OPERATOR_NEW_BRACKET(double,n);
  double *wi=OPERATOR_NEW_BRACKET(double,n);
  double *vl=(V==0 ? 0 : OPERATOR_NEW_BRACKET(double,n*n));
  double *vr=(U==0 ? 0 : OPERATOR_NEW_BRACKET(double,n*n));
  SquareMatrix<double,double> AF(n);
  AF.copy(*this);
  double w;
  int lwork=-1,info;
  F77NAME(dgeev)(jobvl,jobvr,n,AF.addr(),n,wr,wi,vl,n,vr,n,&w,
    lwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dgeev)(jobvl,jobvr,n,AF.addr(),n,wr,wi,vl,n,vr,n,work,
    lwork,info);
  CHECK_TEST(info==0);

  Vector<double,complex<double> > *lambda =
    OPERATOR_NEW Vector<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    complex<double> &ev=(*lambda)[j];
    ev.real()=wr[j];
    ev.imag()=wi[j];
  }
  for (int j=0;j<n;) {
    int jn=j*n;
    if (wi[j]>zero_) {
      int jp1n=jn+n;
      if (V!=0) {
        for (int i=0;i<n;i++) {
          complex<double> &z=V->operator()(i,j);
          z.real()=vl[i+jn];
          z.imag()=vl[i+jp1n];
          complex<double> &zz=V->operator()(i,j+1);
          zz.real()=vl[i+jn];
          zz.imag()=-vl[i+jp1n];
        }
      }
      if (U!=0) {
        for (int i=0;i<n;i++) {
          complex<double> &z=U->operator()(i,j);
          z.real()=vr[i+jn];
          z.imag()=vr[i+jp1n];
          complex<double> &zz=U->operator()(i,j+1);
          zz.real()=vr[i+jn];
          zz.imag()=-vr[i+jp1n];
        }
      }
      j+=2;
    } else {
      if (V!=0) {
        for (int i=0;i<n;i++) {
          complex<double> &z=V->operator()(i,j);
          z.real()=vl[i+jn];
          z.imag()=zero_;
        }
      }
      if (U!=0) {
        for (int i=0;i<n;i++) {
          complex<double> &z=U->operator()(i,j);
          z.real()=vr[i+jn];
          z.imag()=zero_;
        }
      }
      j++;
    }
  }
  delete [] work;
  if (vl!=0) delete [] vr;
  if (vr!=0) delete [] vl;
  delete [] wr;
  delete [] wi;
  return lambda;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,const SquareMatrix<double,double> &S) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<double,double> *P=(m==n ?
    OPERATOR_NEW SquareMatrix<double,double>(n) :
    OPERATOR_NEW Matrix<double,double>(m,n));
  F77NAME(dgemm)('N','N',m,n,n,double_one_,M.addr(),m,S.addr(),n,
    double_zero_,P->addr(),m);
  return P;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "TrapezoidalMatrix.C"

template<> double TrapezoidalMatrix<double,double>::outofbounds_ =
  double_zero_;
template<> double TrapezoidalMatrix<double,double>::safety_ =
  double_zero_;

template class TrapezoidalMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

template<> Matrix<double,double>*
LowerTrapezoidalMatrix<double,double>::makeMatrix() const {
  int m=size(0),n=size(1);
  Matrix<double,double> *M=OPERATOR_NEW
    Matrix<double,double>(m,n,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(this)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(m-j,addr(j,j),1,M->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*M)(j,j)=double_one_;
      if (j+1<m) F77NAME(dcopy)(m-j-1,addr(j+1,j),1,M->addr(j+1,j),1);
    }
  }
  return M;
}

template<> LowerTrapezoidalMatrix<double,double>& 
LowerTrapezoidalMatrix<double,double>::operator+=(
const LowerTrapezoidalMatrix<double,double> &L) {
  bool L_non_unit=
    (dynamic_cast<const UnitLowerTrapezoidalMatrix<double,double>*>(&L)
    ==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0))
  CHECK_SAME(n,L.size(1))
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(m-j,one_,L.addr(j,j),1,addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(daxpy)(m-j-1,one_,L.addr(j+1,j),1,addr(j+1,j),1);
      }
      (*this)(j,j)+=double_one_;
    }
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<double,double>&
LowerTrapezoidalMatrix<double,double>::operator-=(
const LowerTrapezoidalMatrix<double,double> &L) {
  bool L_non_unit=
    (dynamic_cast<const UnitLowerTrapezoidalMatrix<double,double>*>(&L)
    ==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0))
  CHECK_SAME(n,L.size(1))
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(m-j,mone_,L.addr(j,j),1,addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(daxpy)(m-j-1,mone_,L.addr(j+1,j),1,addr(j+1,j),1);
      }
      (*this)(j,j)-=double_one_;
    }
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<double,double>& 
LowerTrapezoidalMatrix<double,double>::operator*=(double d) {
  int m=size(0),n=size(1);
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(m-j,d,addr(j,j),1);
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<double,double>&
LowerTrapezoidalMatrix<double,double>::operator/=(double d) {
  int m=size(0),n=size(1);
  double dinv=double_one_/d;
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(m-j,dinv,addr(j,j),1);
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<double,double>*
LowerTrapezoidalMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<double,double> *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(m,n);
  S->copy(*this);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(m-j,double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(m-j-1,double_one_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> Matrix<double,double>*
LowerTrapezoidalMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *S=OPERATOR_NEW Matrix<double,double>(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(m-j,double_one_,addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> LowerTrapezoidalMatrix<double,double>*
LowerTrapezoidalMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  bool L_non_unit=
   (dynamic_cast<const UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<double,double> *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(m,n);
  S->copy(*this);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(m-j,double_mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(m-j-1,double_mone_,L.addr(j+1,j),1,
        S->addr(j+1,j),1);
      (*S)(j,j)-=double_one_;
    }
  }
  return S;
}

template<> Matrix<double,double>*
LowerTrapezoidalMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *D=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dlaset)('U',m,n,double_zero_,double_zero_,D->addr(),m);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(m-j,addr(j,j),1,D->addr(j,j),1);
    F77NAME(daxpy)(m,double_mone_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> LowerTrapezoidalMatrix<double,double>* 
LowerTrapezoidalMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &L) const {
// to compute the jth column of the product
//   [ L_11   0  ] [ 0 ] = [    0   ]
//   [ L_21 L_22 ] [ m ] = [ L_22 m ]
//   [ L_31 L_32 ]       = [ L_32 m ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  LowerTrapezoidalMatrix<double,double> *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(m,n);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(k-j,L.addr(j,j),1,P->addr(j,j),1);
      F77NAME(dtrmv)('L','N','N',k-j,addr(j,j),m,P->addr(j,j),1);
      if (m>k) {
        F77NAME(dgemv)('N',m-k,k-j,double_one_,addr(k,j),m,L.addr(j,j),1,
          double_zero_,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*P)(j,j)=(*this)(j,j);
      if (j<k-1) {
        F77NAME(dcopy)(k-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
        F77NAME(dtrmv)('L','N','N',k-j-1,addr(j+1,j+1),m,
          P->addr(j+1,j),1);
        F77NAME(daxpy)(k-j-1,double_one_,addr(j+1,j),1,P->addr(j+1,j),1);
      }
      if (m>k) {
        if (j<k-1) {
          F77NAME(dgemv)('N',m-k,k-j-1,double_one_,addr(k,j+1),m,
            L.addr(j+1,j),1,double_zero_,P->addr(k,j),1);
        }
        F77NAME(daxpy)(m-k,double_one_,addr(k,j),1,P->addr(k,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<double,double>*
LowerTrapezoidalMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<double,double> *P=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dlacpy)('A',k,n,M.addr(),k,P->addr(),m);
  F77NAME(dtrmm)('L','L','N','N',k,n,double_one_,addr(),m,
    P->addr(),m);
  if (m>k) {
    F77NAME(dgemm)('N','N',m-k,n,k,double_one_,addr(k,0),m,
      M.addr(),k,double_zero_,P->addr(k,0),m);
  }
  return P;
}

template<> Vector<double,double>*
LowerTrapezoidalMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<double,double> *p=OPERATOR_NEW Vector<double,double>(m);
  F77NAME(dcopy)(n,v.addr(),1,p->addr(),1);
  F77NAME(dtrmv)('L','N','N',n,addr(),m,p->addr(),1);
  if (m>n) {
    F77NAME(dgemv)('N',m-n,n,double_one_,addr(n,0),m,v.addr(),1,
      double_zero_,p->addr(n),1);
  }
  return p;
}

/*
template<> UpperTrapezoidalMatrix<double,double>*
LowerTrapezoidalMatrix<double,double>::transpose() const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<double,double> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(m-j,addr(j,j),1,U->addr(j,j),n);
  }
  return U;
}

template<> UpperTrapezoidalMatrix<double,double>*
LowerTrapezoidalMatrix<double,double>::conjugateTranspose() const {
  return transpose();
}
*/

template<> Vector<double,double>*
LowerTrapezoidalMatrix<double,double>::trmv(
const Vector<double,double> &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<double,double> *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<double,double>(n);
    F77NAME(dcopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(dtrmv)('L','T','N',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(dgemv)('T',m-n,n,one_,addr(n,0),m,x.addr(n),1,one_,
        p->addr(),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<double,double>(m);
    F77NAME(dcopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(dtrmv)('L','N','N',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(dgemv)('N',m-n,n,one_,addr(n,0),m,x.addr(),1,zero_,
        p->addr(n),1);
    }
  }
  return p;
}

template<> Matrix<double,double>*
LowerTrapezoidalMatrix<double,double>::trmm(
const Matrix<double,double> &X,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<double,double> *P=0;
  if (side=='L' || side=='l') {
    int k=X.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,X.size(0));
      P=OPERATOR_NEW Matrix<double,double>(n,k);
      F77NAME(dlacpy)('A',n,k,X.addr(),m,P->addr(),n);
      F77NAME(dtrmm)('L','L','T','N',n,k,one_,addr(),m,P->addr(),n);
      if (m>n) {
        F77NAME(dgemm)('T','N',n,k,m-n,one_,addr(n,0),m,X.addr(n,0),m,
          one_,P->addr(),n);
      }
    } else {
      CHECK_SAME(n,X.size(0));
      P=OPERATOR_NEW Matrix<double,double>(m,k);
      F77NAME(dlacpy)('A',n,k,X.addr(),n,P->addr(),m);
      F77NAME(dtrmm)('L','L','N','N',n,k,one_,addr(),m,P->addr(),m);
      if (m>n) {
        F77NAME(dgemm)('N','N',m-n,k,n,one_,addr(n,0),m,X.addr(),n,
          zero_,P->addr(n,0),m);
      }
    }
  } else {
    int k=X.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,X.size(1));
      P=OPERATOR_NEW Matrix<double,double>(k,m);
      F77NAME(dlacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(dtrmm)('R','L','T','N',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(dgemm)('N','T',k,m-n,n,one_,X.addr(),k,addr(n,0),m,
          zero_,P->addr(0,n),k);
      }
    } else {
      CHECK_SAME(m,X.size(1));
      P=OPERATOR_NEW Matrix<double,double>(k,n);
      F77NAME(dlacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(dtrmm)('R','L','N','N',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(dgemm)('N','N',k,n,m-n,one_,X.addr(0,n),k,addr(n,0),m,
          one_,P->addr(),k);
      }
    }
  }
  return P;
}

template<> void LowerTrapezoidalMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(dcopy)(n,b.addr(),1,x.addr(),1);
    if (m>n) { // use trailing entries of x as free variables
      F77NAME(dgemv)('T',m-n,n,mone_,addr(n,0),m,x.addr(n),1,one_,
        x.addr(),1);
    }
    F77NAME(dtrsv)('L','T','N',n,addr(),m,x.addr(),1);
  } else { // calling routine will have to check consistency conditions
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(dcopy)(n,b.addr(),1,x.addr(),1);
    F77NAME(dtrsv)('L','N','N',n,addr(),m,x.addr(),1);
  }
}

template<> void LowerTrapezoidalMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(dlacpy)('A',n,k,B.addr(),n,X.addr(),n);
      if (m>n) { // use these entries of X as free variables
        F77NAME(dgemm)('T','N',n,k,m-n,mone_,addr(n,0),m,X.addr(n,0),m,
          one_,X.addr(),m);
      }
      F77NAME(dtrsm)('L','L','T','N',n,k,one_,addr(),m,X.addr(),m);
    } else { // calling routine will have to check consistency conditions
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(dlacpy)('A',n,k,B.addr(),m,X.addr(),n);
      F77NAME(dtrsm)('L','L','N','N',n,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(dlacpy)('A',k,n,B.addr(),k,X.addr(),k);
      F77NAME(dtrsm)('R','L','T','N',k,n,one_,addr(),m,X.addr(),k);
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(dlacpy)('A',k,n,B.addr(),k,X.addr(),k);
      if (m>n) { // use these entries of X as free variables
        F77NAME(dgemm)('N','N',k,n,m-n,mone_,X.addr(0,n),n,addr(n,0),m,
          one_,X.addr(),k);
      }
      F77NAME(dtrsm)('R','L','N','N',k,n,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> double LowerTrapezoidalMatrix<double,double>::normFrobenius()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('F','L','N',m,size(1),addr(),m,work);
}

template<> double LowerTrapezoidalMatrix<double,double>::normInfinity()
const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double result=F77NAME(dlantr)('I','L','N',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> double LowerTrapezoidalMatrix<double,double>::normMaxEntry()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('M','L','N',m,size(1),addr(),m,work);
}

template<> double LowerTrapezoidalMatrix<double,double>::normOne()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('O','L','N',m,size(1),addr(),m,work);
}

template<> Matrix<double,double>* operator+(
const Matrix<double,double> &M,
const LowerTrapezoidalMatrix<double,double> &L) {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *S=OPERATOR_NEW Matrix<double,double>(m,n);
  S->copy(M);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(m-j,double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(daxpy)(m-j-1,double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> Matrix<double,double>* operator-(
const Matrix<double,double> &M,
const LowerTrapezoidalMatrix<double,double> &L) {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *D=OPERATOR_NEW Matrix<double,double>(m,n);
  D->copy(M);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(m-j,double_mone_,L.addr(j,j),1,D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(daxpy)(m-j-1,double_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
      (*D)(j,j)-=double_one_;
    }
  }
  return D;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,
const LowerTrapezoidalMatrix<double,double> &L) {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  char diagL=
    (dynamic_cast<const UnitLowerTrapezoidalMatrix<double,double>*>(&L)
    ==0 ? 'N' : 'U');
  int m=M.size(0),k=L.size(0),n=L.size(1);
  CHECK_SAME(k,M.size(1));
  Matrix<double,double> *P=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dlacpy)('A',m,n,M.addr(),m,P->addr(),m);
  F77NAME(dtrmm)('R','L','N',diagL,m,n,double_one_,L.addr(),k,
    P->addr(),m);
  if (k>n) {
    F77NAME(dgemm)('N','N',m,n,k-n,double_one_,M.addr(0,n),m,
      L.addr(n,0),k,double_one_,P->addr(),m);
  }
  return P;
}

template class LowerTrapezoidalMatrix<double,double>;
template void testLowerTrapezoidalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double
LowerTriangularMatrix<double,double>::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  F77NAME(dtrcon)(norm,'L','N',n,addr(),n,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete [] work;
  return rcond;
}

/*
template<> LowerTriangularMatrix<double,double>*
LowerTriangularMatrix<double,double>::inverse() const {
  LowerTriangularMatrix<double,double> *L=
    OPERATOR_NEW LowerTriangularMatrix<double,double>(*this);
  int n=size(0);
  int info;
  F77NAME(dtrtri)('L','N',n,L->addr(),n,info);
  CHECK_TEST(info==0);
  return L;
}
*/

template class LowerTriangularMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void
UnitLowerTrapezoidalMatrix<double,double>::copyFrom(int m,int n,
const Matrix<double,double> &L) {
  m=min(m,min(size(0),L.size(0)));
  n=min(n,min(size(1),L.size(1)));
  for (int j=0;j<n;j++) {
    if (j+1<m) {
      F77NAME(dcopy)(m-j-1,L.addr(j+1,j),1,addr(j+1,j),1);
    }
  }
}

/*
template<> UnitUpperTrapezoidalMatrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::transpose() const {
  int m=size(0),n=size(1);
  UnitUpperTrapezoidalMatrix<double,double> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<double,double>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(m-j-1,addr(j+1,j),1,U->addr(j,j+1),n);
  }
  return U;
}

template<> UnitUpperTrapezoidalMatrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::conjugateTranspose() const {
  return transpose();
}
*/

template<> LowerTrapezoidalMatrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<double,double> *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(m,n);
  if (L_non_unit) {
    S->copy(L);
    for (int j=0;j<n;j++) {
      (*S)(j,j)+=double_one_;
      F77NAME(daxpy)(m-j-1,double_one_,addr(j+1,j),1,S->addr(j+1,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=2.;
      if (j<m-1) {
        F77NAME(dcopy)(m-j-1,addr(j+1,j),1,S->addr(j+1,j),1);
        F77NAME(daxpy)(m-j-1,double_one_,addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> Matrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *S=OPERATOR_NEW Matrix<double,double>(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    (*S)(j,j)+=double_one_;
    F77NAME(daxpy)(m-j-1,double_one_,addr(j+1,j),1,S->addr(j+1,j),1);
  }
  return S;
}

template<> LowerTrapezoidalMatrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<double,double> *dif=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(m,n);
  for (int j=0;j<n;j++) {
    (*dif)(j,j)=(L_non_unit ? double_one_-L(j,j) : double_zero_);
    if (j<m-1) {
      F77NAME(dcopy)(m-j-1,addr(j+1,j),1,dif->addr(j+1,j),1);
      F77NAME(daxpy)(m-j-1,double_mone_,L.addr(j+1,j),1,
        dif->addr(j+1,j),1);
    }
  }
  return dif;
}

template<> Matrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *dif=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dlaset)('U',m,n,double_zero_,double_zero_,dif->addr(),m);
  for (int j=0;j<n;j++) {
    (*dif)(j,j)=double_one_-M(j,j);
    F77NAME(dcopy)(m-j-1,addr(j+1,j),1,dif->addr(j+1,j),1);
    F77NAME(daxpy)(m-j-1,double_mone_,M.addr(j+1,j),1,dif->addr(j+1,j),1);
  }
  return dif;
}

template<> LowerTrapezoidalMatrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator*(double d) const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<double,double> *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(m,n);
  for (int j=0;j<n;j++) {
    (*P)(j,j)=d;
    if (j<m-1) {
      F77NAME(dcopy)(m-j-1,addr(j+1,j),1,P->addr(j+1,j),1);
      F77NAME(dscal)(m-j-1,d,P->addr(j+1,j),1);
    }
  }
  return P;
}

template<> LowerTrapezoidalMatrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator/(double d) const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<double,double> *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(m,n);
  double dinv=double_one_/d;
  for (int j=0;j<n;j++) {
    (*P)(j,j)=dinv;
    if (j<m-1) {
      F77NAME(dcopy)(m-j-1,addr(j+1,j),1,P->addr(j+1,j),1);
      F77NAME(dscal)(m-j-1,dinv,P->addr(j+1,j),1);
    }
  }
  return P;
}

template<> LowerTrapezoidalMatrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &L) const {
// to compute the jth column of the product
//   [ L_11   0  ] [ 0 ] = [    0   ]
//   [ L_21 L_22 ] [ m ] = [ L_22 m ]
//   [ L_31 L_32 ]       = [ L_32 m ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  LowerTrapezoidalMatrix<double,double> *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(m,n);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(k-j,L.addr(j,j),1,P->addr(j,j),1);
      F77NAME(dtrmv)('L','N','U',k-j,addr(j,j),m,P->addr(j,j),1);
      if (m>k) {
        F77NAME(dgemv)('N',m-k,k-j,double_one_,addr(k,j),m,L.addr(j,j),1,
          double_zero_,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*P)(j,j)=double_one_;
      if (j<k-1) {
        F77NAME(dcopy)(k-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
        F77NAME(dtrmv)('L','N','U',k-j-1,addr(j+1,j+1),m,
          P->addr(j+1,j),1);
        F77NAME(daxpy)(k-j-1,double_one_,addr(j+1,j),1,P->addr(j+1,j),1);
      }
      if (m>k) {
        if (j<k-1) {
          F77NAME(dgemv)('N',m-k,k-j-1,double_one_,addr(k,j+1),m,
            L.addr(j+1,j),1,double_zero_,P->addr(k,j),1);
        }
        F77NAME(daxpy)(m-k,double_one_,addr(k,j),1,P->addr(k,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  F77NAME(dlacpy)('A',k,n,M.addr(),k,P->addr(),m);
  F77NAME(dtrmm)('L','L','N','U',k,n,double_one_,addr(),m,
    P->addr(),m);
  if (m>k) {
    F77NAME(dgemm)('N','N',m-k,n,k,double_one_,addr(k,0),m,
      M.addr(),k,double_one_,P->addr(k,0),m);
  }
  return P;
}

template<> Vector<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<double,double> *p=
    OPERATOR_NEW Vector<double,double>(m,double_zero_);
  F77NAME(dcopy)(n,v.addr(),1,p->addr(),1);
  F77NAME(dtrmv)('L','N','U',n,addr(),m,p->addr(),1);
  if (m>n) {
    F77NAME(dgemv)('N',m-n,n,double_one_,addr(n,0),m,
      v.addr(),1,double_one_,p->addr(n),m);
  }
  return p;
}

template<> Vector<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::trmv(
const Vector<double,double> &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<double,double> *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<double,double>(n);
    F77NAME(dcopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(dtrmv)('L','T','U',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(dgemv)('T',m-n,n,one_,addr(n,0),m,x.addr(n),1,one_,
        p->addr(),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<double,double>(m);
    F77NAME(dcopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(dtrmv)('L','N','U',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(dgemv)('N',m-n,n,one_,addr(n,0),m,x.addr(),1,zero_,
        p->addr(n),1);
    }
  }
  return p;
}

template<> Matrix<double,double>*
UnitLowerTrapezoidalMatrix<double,double>::trmm(
const Matrix<double,double> &X,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<double,double> *P=0;
  if (side=='L' || side=='l') {
    int k=X.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,X.size(0));
      P=OPERATOR_NEW Matrix<double,double>(n,k);
      F77NAME(dlacpy)('A',n,k,X.addr(),m,P->addr(),n);
      F77NAME(dtrmm)('L','L','T','U',n,k,one_,addr(),m,P->addr(),n);
      if (m>n) {
        F77NAME(dgemm)('T','N',n,k,m-n,one_,addr(n,0),m,X.addr(n,0),m,
          one_,P->addr(),n);
      }
    } else {
      CHECK_SAME(n,X.size(0));
      P=OPERATOR_NEW Matrix<double,double>(m,k);
      F77NAME(dlacpy)('A',n,k,X.addr(),n,P->addr(),m);
      F77NAME(dtrmm)('L','L','N','U',n,k,one_,addr(),m,P->addr(),m);
      if (m>n) {
        F77NAME(dgemm)('N','N',m-n,k,n,one_,addr(n,0),m,X.addr(),n,
          zero_,P->addr(n,0),m);
      }
    }
  } else {
    int k=X.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,X.size(1));
      P=OPERATOR_NEW Matrix<double,double>(k,m);
      F77NAME(dlacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(dtrmm)('R','L','T','U',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(dgemm)('N','T',k,m-n,n,one_,X.addr(),k,addr(n,0),m,
          zero_,P->addr(0,n),k);
      }
    } else {
      CHECK_SAME(m,X.size(1));
      P=OPERATOR_NEW Matrix<double,double>(k,n);
      F77NAME(dlacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(dtrmm)('R','L','N','U',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(dgemm)('N','N',k,n,m-n,one_,X.addr(0,n),k,addr(n,0),m,
          one_,P->addr(),k);
      }
    }
  }
  return P;
}

template<> void UnitLowerTrapezoidalMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(dcopy)(n,b.addr(),1,x.addr(),1);
    if (m>n) { // use trailing entries of x as free variables
      F77NAME(dgemv)('T',m-n,n,mone_,addr(n,0),m,x.addr(n),1,one_,
        x.addr(),1);
    }
    F77NAME(dtrsv)('L','T','U',n,addr(),m,x.addr(),1);
  } else { // calling routine will have to check consistency conditions
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(dcopy)(n,b.addr(),1,x.addr(),1);
    F77NAME(dtrsv)('L','N','U',n,addr(),m,x.addr(),1);
  }
}

template<> void UnitLowerTrapezoidalMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(dlacpy)('A',n,k,B.addr(),n,X.addr(),n);
      if (m>n) { // use these entries of X as free variables
        F77NAME(dgemm)('T','N',n,k,m-n,mone_,addr(n,0),m,X.addr(n,0),m,
          one_,X.addr(),m);
      }
      F77NAME(dtrsm)('L','L','T','U',n,k,one_,addr(),m,X.addr(),m);
    } else { // calling routine will have to check consistency conditions
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(dlacpy)('A',n,k,B.addr(),m,X.addr(),n);
      F77NAME(dtrsm)('L','L','N','U',n,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(dlacpy)('A',k,n,B.addr(),k,X.addr(),k);
      F77NAME(dtrsm)('R','L','T','U',k,n,one_,addr(),m,X.addr(),k);
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(dlacpy)('A',k,n,B.addr(),k,X.addr(),k);
      if (m>n) { // use these entries of X as free variables
        F77NAME(dgemm)('N','N',k,n,m-n,mone_,X.addr(0,n),n,addr(n,0),m,
          one_,X.addr(),k);
      }
      F77NAME(dtrsm)('R','L','N','U',k,n,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> double
UnitLowerTrapezoidalMatrix<double,double>::normFrobenius() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('F','L','U',m,size(1),addr(),m,work);
}

template<> double
UnitLowerTrapezoidalMatrix<double,double>::normInfinity() const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double result=F77NAME(dlantr)('I','L','U',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> double
UnitLowerTrapezoidalMatrix<double,double>::normMaxEntry() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('M','L','U',m,size(1),addr(),m,work);
}

template<> double
UnitLowerTrapezoidalMatrix<double,double>::normOne() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('O','L','U',m,size(1),addr(),m,work);
}

template class UnitLowerTrapezoidalMatrix<double,double>;
template void testUnitLowerTrapezoidalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double
UnitLowerTriangularMatrix<double,double>::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  F77NAME(dtrcon)(norm,'L','U',n,addr(),n,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete [] work;
  return rcond;
}

/*
template<> UnitLowerTriangularMatrix<double,double>*
UnitLowerTriangularMatrix<double,double>::inverse() const {
  UnitLowerTriangularMatrix<double,double> *result=
    OPERATOR_NEW UnitLowerTriangularMatrix<double,double>(*this);
  int n=size(0);
  int info;
  F77NAME(dtrtri)('L','U',n,result->addr(),n,info);
  CHECK_TEST(info==0);
  return result;
}
*/

template class UnitLowerTriangularMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> Matrix<double,double>*
UpperTrapezoidalMatrix<double,double>::makeMatrix() const {
  int m=size(0),n=size(1);
  Matrix<double,double> *M=OPERATOR_NEW
    Matrix<double,double>(m,n,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(this)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+1,m),addr(0,j),1,M->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j,m),addr(0,j),1,M->addr(0,j),1);
      if (j<m) (*M)(j,j)=double_one_;
    }
  }
  return M;
}

template<> UpperTrapezoidalMatrix<double,double>& 
UpperTrapezoidalMatrix<double,double>::operator+=(
const UpperTrapezoidalMatrix<double,double> &M) {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0))
  CHECK_SAME(n,M.size(1))
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(j+1,m),one_,M.addr(0,j),1,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<double,double>&
UpperTrapezoidalMatrix<double,double>::operator-=(
const UpperTrapezoidalMatrix<double,double> &M) {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0))
  CHECK_SAME(n,M.size(1))
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(j+1,m),mone_,M.addr(0,j),1,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<double,double>& 
UpperTrapezoidalMatrix<double,double>::operator*=(double d) {
  int m=size(0),n=size(1);
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(min(j+1,m),d,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<double,double>&
UpperTrapezoidalMatrix<double,double>::operator/=(double d) {
  int m=size(0),n=size(1);
  double dinv=double_one_/d;
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(min(j+1,m),dinv,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<double,double>*
UpperTrapezoidalMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<double,double> *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(m,n);
  if (U_non_unit) {
    S->copy(U);
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(min(j+1,m),double_one_,addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      if (j>0) {
        F77NAME(daxpy)(min(j,m),double_one_,addr(0,j),1,S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
UpperTrapezoidalMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(m,zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(j+1,addr(0,j),1,S->addr(0,j),1);
      F77NAME(daxpy)(m-j,double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(j+1,addr(0,j),1,S->addr(0,j),1);
      if (j<m-1) {
        F77NAME(daxpy)(m-j-1,double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> Matrix<double,double>*
UpperTrapezoidalMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *S=OPERATOR_NEW Matrix<double,double>(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(m,j+1),double_one_,addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> UpperTrapezoidalMatrix<double,double>*
UpperTrapezoidalMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<double,double> *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(m,n);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      F77NAME(daxpy)(min(j+1,m),double_mone_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      if (j>0) {
        F77NAME(daxpy)(min(j,m),double_mone_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)-=double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
UpperTrapezoidalMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *D=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(j+1,addr(0,j),1,D->addr(0,j),1);
      F77NAME(daxpy)(m-j,double_mone_,L.addr(j,j),1,D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(j+1,addr(0,j),1,D->addr(0,j),1);
      if (j<m-1) {
        F77NAME(daxpy)(m-j-1,double_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
      (*D)(j,j)-=double_one_;
    }
  }
  return D;
}

template<> Matrix<double,double>*
UpperTrapezoidalMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *D=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(j+1,addr(0,j),1,D->addr(0,j),1);
    F77NAME(daxpy)(m,double_mone_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> UpperTrapezoidalMatrix<double,double>* 
UpperTrapezoidalMatrix<double,double>::operator*(
const UpperTrapezoidalMatrix<double,double> &U) const {
// [ U_1 , U_2 ] [ V_11 , V_12 , V_13 ]
//               [   0  , V_22 , V_23 ]
//   = [ U_1 V_11 , U_1 V_12 + U_2 V_22 , U_1 V_13 + U_2 V_23 ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  int m=size(0),k=size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  UpperTrapezoidalMatrix<double,double> *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(m,n,double_zero_);
  if (U_non_unit) {
    for (int j=0;j<m;j++) {
      F77NAME(dcopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(dtrmv)('U','N','N',j+1,addr(),m,P->addr(0,j),1);
    }
    for (int j=m;j<k;j++) {
      F77NAME(dcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(dtrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
      F77NAME(dgemv)('N',m,j-m+1,double_one_,addr(0,m),m,U.addr(m,j),1,
        double_one_,P->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(dcopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(dtrmv)('U','N','N',j,addr(),m,P->addr(0,j),1);
        F77NAME(daxpy)(j,double_one_,addr(0,j),1,P->addr(0,j),1);
      }
      (*P)(j,j)=(*this)(j,j);
    }
    for (int j=m;j<k;j++) {
      F77NAME(dcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(dtrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
      if (j>m) {
        F77NAME(dgemv)('N',m,j-m,double_one_,addr(0,m),m,U.addr(m,j),1,
          double_one_,P->addr(0,j),1);
      }
      F77NAME(daxpy)(m,double_one_,addr(0,j),1,P->addr(0,j),1);
    }
  }
  for (int j=k;j<n;j++) {
    F77NAME(dcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
    F77NAME(dtrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
    F77NAME(dgemv)('N',m,k-m,double_one_,addr(0,m),m,U.addr(m,j),1,
      double_one_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<double,double>*
UpperTrapezoidalMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  if (m>=n) {
//  note that
//  [ U_11 U_12 U_13 ] [ L_1 ] = [ U_11 L_1 + U_12 L_2 + U_13 L_3 ]
//  [      U_22 U_23 ] [ L_2 ] = [            U_22 L_2 + U_23 L_3 ]
//                     [ L_3 ]
    if (L_non_unit) { // U_11 L_1:
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(dgemv)('N',j,n-j,double_one_,addr(0,j),m,L.addr(j,j),1,
            double_zero_,P->addr(0,j),1);
        }
        F77NAME(dcopy)(n-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(dtrmv)('U','N','N',n-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(daxpy)(j,double_one_,addr(0,j),1,P->addr(0,j),1);
          if (j<n-1) {
            F77NAME(dgemv)('N',j,n-j-1,double_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,double_zero_,P->addr(0,j),1);
          }
        }
        if (j<n-1) {
          F77NAME(dcopy)(n-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(dtrmv)('U','N','N',n-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(ddot)(n-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (n<k) { // U_12 L_2 + U_13 L_3:
      F77NAME(dgemm)('N','N',n,n,k-n,double_one_,addr(0,n),m,
        L.addr(n,0),k,double_one_,P->addr(0,0),m);
    }
    if (n<m) {
      F77NAME(dlacpy)('A',m-n,n,L.addr(n,0),k,P->addr(n,0),m);
      F77NAME(dtrmm)('L','U','N','N',m-n,n,double_one_,addr(n,n),m,
        P->addr(n,0),m); // U_22 L_2
      if (m<k) {
        F77NAME(dgemm)('N','N',m-n,n,k-m,double_one_,addr(n,m),m,
          L.addr(m,0),k,double_one_,P->addr(n,0),m); // U_23 L_3
      }
    }
  } else {
//  [ U_1 U_2 U_3 ] [ L_11   0  ]
//                  [ L_21 L_22 ]
//                  [ L_31 L_32 ]
//  = [ U_1 L_11 + U_2 L_21 _ U_3 L_31 , U_2 L_22 , U_2 L_22 + U_3 L_32 ]
    if (L_non_unit) { // U_1 L_11:
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(dgemv)('N',j,m-j,double_one_,addr(0,j),m,L.addr(j,j),1,
            double_zero_,P->addr(0,j),1);
        }
        F77NAME(dcopy)(m-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(dtrmv)('U','N','N',m-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(daxpy)(j,double_one_,addr(0,j),1,P->addr(0,j),1);
          if (j<m-1) {
            F77NAME(dgemv)('N',j,m-j-1,double_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,double_zero_,P->addr(0,j),1);
          }
        }
        if (j<m-1) {
          F77NAME(dcopy)(m-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(dtrmv)('U','N','N',m-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(ddot)(m-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (m<k) { // U_2 L_21 + U_3 L_31 
      F77NAME(dgemm)('N','N',m,m,k-m,double_one_,addr(0,m),m,
        L.addr(m,0),k,double_one_,P->addr(0,0),m);
    }
    char diagL=(L_non_unit ? 'N' : 'U');
    F77NAME(dlacpy)('A',m,n-m,addr(0,m),m,P->addr(0,m),m);
    F77NAME(dtrmm)('R','L','N',diagL,m,n-m,double_one_,L.addr(m,m),k,
      P->addr(0,m),m); // U_2 L_22
    if (n<k) {
      F77NAME(dgemm)('N','N',m,n-m,k-n,double_one_,addr(0,n),m,
        L.addr(n,m),k,double_one_,P->addr(0,m),m); // U_3 L_32
    }
  }
  return P;
}

template<> Matrix<double,double>*
UpperTrapezoidalMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
// [ U_1 U_2 ] [ M_1 ] = U_1 M_1 + U_2 M_2
//             [ M_2 ]
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<double,double> *P=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dlacpy)('A',m,n,M.addr(),k,P->addr(),m);
  F77NAME(dtrmm)('L','U','N','N',m,n,double_one_,addr(),m,
    P->addr(),m); // U_1 M_1
  if (k>m) {
    F77NAME(dgemm)('N','N',m,n,k-m,double_one_,addr(0,m),m,
      M.addr(m,0),k,double_one_,P->addr(),m); // U_2 M_2
  }
  return P;
}

template<> Vector<double,double>*
UpperTrapezoidalMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<double,double> *p=OPERATOR_NEW Vector<double,double>(m);
  F77NAME(dcopy)(m,v.addr(),1,p->addr(),1);
  F77NAME(dtrmv)('U','N','N',m,addr(),m,p->addr(),1);
  F77NAME(dgemv)('N',m,n-m,double_one_,addr(0,m),m,v.addr(m),1,
    double_one_,p->addr(),1);
  return p;
}

/*
template<> LowerTrapezoidalMatrix<double,double>*
UpperTrapezoidalMatrix<double,double>::transpose() const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<double,double> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(j+1,m),addr(0,j),1,L->addr(j,0),n);
  }
  return L;
}

template<> LowerTrapezoidalMatrix<double,double>*
UpperTrapezoidalMatrix<double,double>::conjugateTranspose() const {
  return transpose();
}
*/

template<> Vector<double,double>*
UpperTrapezoidalMatrix<double,double>::trmv(
const Vector<double,double> &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<double,double> *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<double,double>(n);
    F77NAME(dcopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(dtrmv)('U','T','N',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(dgemv)('T',m,n-m,one_,addr(0,m),m,x.addr(),1,zero_,
        p->addr(m),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<double,double>(m);
    F77NAME(dcopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(dtrmv)('U','N','N',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(dgemv)('N',m,n-m,one_,addr(0,m),m,x.addr(m),1,one_,
        p->addr(),1);
    }
  }
  return p;
}

template<> Matrix<double,double>*
UpperTrapezoidalMatrix<double,double>::trmm(
const Matrix<double,double> &M,char side,char trans) const {
//TRACER_CALL(tr,"UpperTrapezoidalMatrix::trmm");
  int m=size(0),n=size(1);
  Matrix<double,double> *P=0;
  if (side=='L' || side=='l') {
    int k=M.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,M.size(0));
      P=OPERATOR_NEW Matrix<double,double>(n,k);
      F77NAME(dlacpy)('A',m,k,M.addr(),m,P->addr(),n);
      F77NAME(dtrmm)('L','U','T','N',m,k,one_,addr(),m,P->addr(),n);
      if (n>m) {
        F77NAME(dgemm)('T','N',n-m,k,m,one_,addr(0,m),m,M.addr(),m,
          zero_,P->addr(m,0),n);
      }
    } else {
      CHECK_SAME(n,M.size(0));
      P=OPERATOR_NEW Matrix<double,double>(m,k);
      F77NAME(dlacpy)('A',m,k,M.addr(),n,P->addr(),m);
      F77NAME(dtrmm)('L','U','N','N',m,k,one_,addr(),m,P->addr(),m);
      if (n>m) {
        F77NAME(dgemm)('N','N',m,k,n-m,one_,addr(0,m),m,M.addr(m,0),n,
          one_,P->addr(),m);
      }
    }
  } else {
    int k=M.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,M.size(1));
      P=OPERATOR_NEW Matrix<double,double>(k,m);
      F77NAME(dlacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(dtrmm)('R','U','T','N',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(dgemm)('N','T',k,m,n-m,one_,M.addr(0,m),k,addr(0,m),m,
          one_,P->addr(),k);
      }
    } else {
      CHECK_SAME(m,M.size(1));
      P=OPERATOR_NEW Matrix<double,double>(k,n);
      F77NAME(dlacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(dtrmm)('R','U','N','N',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(dgemm)('N','N',k,n-m,m,one_,M.addr(),k,addr(0,m),m,
          zero_,P->addr(0,m),k);
      }
    }
  }
  return P;
}

template<> double UpperTrapezoidalMatrix<double,double>::normFrobenius()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('F','U','N',m,size(1),addr(),m,work);
}

template<> double UpperTrapezoidalMatrix<double,double>::normInfinity()
const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double result=F77NAME(dlantr)('I','U','N',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> double UpperTrapezoidalMatrix<double,double>::normMaxEntry()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('M','U','N',m,size(1),addr(),m,work);
}

template<> double UpperTrapezoidalMatrix<double,double>::normOne()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('O','U','N',m,size(1),addr(),m,work);
}

template<> void UpperTrapezoidalMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    // calling routine will have to check consistency conditions
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(dcopy)(m,b.addr(),1,x.addr(),1);
    F77NAME(dtrsv)('U','T','N',m,addr(),m,x.addr(),1);
  } else {
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(dcopy)(m,b.addr(),1,x.addr(),1);
    if (n>m) { // use trailing entries of x as free variables
      F77NAME(dgemv)('N',m,n-m,mone_,addr(0,m),m,x.addr(m),1,one_,
        x.addr(),1);
    }
    F77NAME(dtrsv)('U','N','N',m,addr(),m,x.addr(),1);
  }
}

template<> void UpperTrapezoidalMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(dlacpy)('A',m,k,B.addr(),n,X.addr(),m);
      F77NAME(dtrsm)('L','U','T','N',m,k,one_,addr(),m,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(dlacpy)('A',m,k,B.addr(),m,X.addr(),n);
      if (n>m) { // use these entries of X as free variables
        F77NAME(dgemm)('N','N',m,k,n-m,mone_,addr(0,m),m,X.addr(m,0),n,
          one_,X.addr(),n);
      }
      F77NAME(dtrsm)('L','U','N','N',m,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(dlacpy)('A',k,m,B.addr(),k,X.addr(),k);
      if (n>m) { // use these entries of X as free variables
        F77NAME(dgemm)('N','T',k,m,n-m,mone_,X.addr(0,m),k,addr(0,m),m,
          one_,X.addr(),k);
      }
      F77NAME(dtrsm)('R','U','T','N',k,m,one_,addr(),m,X.addr(),k);
    } else { // calling routine will have to check consistency
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(dlacpy)('A',k,m,B.addr(),k,X.addr(),k);
      F77NAME(dtrsm)('R','U','N','N',k,m,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> SquareMatrix<double,double>* operator+(
const UnitLowerTrapezoidalMatrix<double,double> &L,
const UpperTrapezoidalMatrix<double,double> &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    (*S)(j,j)=double_one_;
    if (j<m-1) {
      F77NAME(dcopy)(m-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
    }
    if (U_non_unit) {
      F77NAME(daxpy)(j+1,double_one_,U.addr(0,j),1,S->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(daxpy)(j,double_one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator+(
const LowerTrapezoidalMatrix<double,double> &L,
const UpperTrapezoidalMatrix<double,double> &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    if (L_non_unit) {
      F77NAME(dcopy)(m-j,L.addr(j,j),1,S->addr(j,j),1);
    } else {
      (*S)(j,j)=double_one_;
      if (j<m-1) {
        F77NAME(dcopy)(m-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
    }
    if (U_non_unit) {
      F77NAME(daxpy)(j+1,double_one_,U.addr(0,j),1,S->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(daxpy)(j,double_one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> Matrix<double,double>* operator+(
const Matrix<double,double> &M,
const UpperTrapezoidalMatrix<double,double> &U) {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(not_trapezoidal);
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  Matrix<double,double> *S=OPERATOR_NEW Matrix<double,double>(m,n);
  S->copy(M);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(min(j+1,m),double_one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j,m),double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const UnitLowerTrapezoidalMatrix<double,double> &L,
const UpperTrapezoidalMatrix<double,double> &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  SquareMatrix<double,double> *D=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    (*D)(j,j)=double_one_;
    if (j<m-1) {
      F77NAME(dcopy)(m-j-1,L.addr(j+1,j),1,D->addr(j+1,j),1);
    }
    if (U_non_unit) {
      F77NAME(daxpy)(j+1,double_mone_,U.addr(0,j),1,D->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(daxpy)(j,double_mone_,U.addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)-=double_one_;
    }
  }
  return D;
}

template<> SquareMatrix<double,double>* operator-(
const LowerTrapezoidalMatrix<double,double> &L,
const UpperTrapezoidalMatrix<double,double> &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  SquareMatrix<double,double> *D=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    if (L_non_unit) {
      F77NAME(dcopy)(m-j,L.addr(j,j),1,D->addr(j,j),1);
    } else {
      (*D)(j,j)=double_one_;
      if (j<m-1) {
        F77NAME(dcopy)(m-j-1,L.addr(j+1,j),1,D->addr(j+1,j),1);
      }
    }
    if (U_non_unit) {
      F77NAME(daxpy)(j+1,double_mone_,U.addr(0,j),1,D->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(daxpy)(j,double_mone_,U.addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)-=double_one_;
    }
  }
  return D;
}

template<> Matrix<double,double>* operator-(
const Matrix<double,double> &M,
const UpperTrapezoidalMatrix<double,double> &U) {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(not_trapezoidal);
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  Matrix<double,double> *D=OPERATOR_NEW Matrix<double,double>(m,n);
  D->copy(M);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(min(j+1,m),double_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j,m),double_mone_,U.addr(0,j),1,
          D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)-=double_one_;
    }
  }
  return D;
}

template<> Matrix<double,double>* operator*(
const UnitLowerTrapezoidalMatrix<double,double> &L,
const UpperTrapezoidalMatrix<double,double> &U) {
//  Note that
//  [ L_1 ] [ U_1 U_2 ] = [ L_1 U_1 , L_1 U_2 ]
//  [ L_2 ]             = [ L_2 U_1 , L_2 U_2 ]
//  and that
//  [ L_11      ] [ u ] = [ L_11 u ]
//  [ L_21 L_22 ] [   ] = [ L_21 u ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U) ==0);
  int m=L.size(0),k=L.size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<k;j++) { // L_1 U_1
    if (U_non_unit) {
      F77NAME(dcopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(dtrmv)('L','N','U',j+1,L.addr(),m,P->addr(0,j),1);
      if (j+1<k) {
        F77NAME(dgemv)('N',k-j-1,j+1,double_one_,L.addr(j+1,0),m,
          U.addr(0,j),1,double_zero_,P->addr(j+1,j),1);
      }
    } else {
      if (j>0) {
        F77NAME(dcopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(dtrmv)('L','N','U',j,L.addr(),m,P->addr(0,j),1);
        (*P)(j,j)=F77NAME(ddot)(j,L.addr(j,0),m,U.addr(0,j),1);
        if (j+1<k) {
          F77NAME(dgemv)('N',k-j-1,j,double_one_,L.addr(j+1,0),m,
            U.addr(0,j),1,double_zero_,P->addr(j+1,j),1);
        }
      }
      (*P)(j,j)+=double_one_;
      if (j+1<k) {
        F77NAME(daxpy)(k-j-1,double_one_,L.addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
    }
  }
  if (n>k) { // L_1 U_2
    F77NAME(dlacpy)('A',k,n-k,U.addr(0,k),k,P->addr(0,k),m);
    F77NAME(dtrmm)('L','L','N','U',k,n-k,double_one_,L.addr(),m,
      P->addr(0,k),m);
  }
  if (m>k) {
    char diagU=(U_non_unit ? 'N' : 'U');
    F77NAME(dlacpy)('A',m-k,k,L.addr(k,0),m,P->addr(k,0),m);
    F77NAME(dtrmm)('R','U','N',diagU,m-k,k,double_one_,U.addr(),k,
      P->addr(k,0),m); // L_2 U_1
    if (n>k) { // L_2 U_2
      F77NAME(dgemm)('N','N',m-k,n-k,k,double_one_,L.addr(k,0),m,
        U.addr(0,k),k,double_zero_,P->addr(k,k),m);
    }
  }
  return P;
}

template<> Matrix<double,double>* operator*(
const LowerTrapezoidalMatrix<double,double> &L,
const UpperTrapezoidalMatrix<double,double> &U) {
//  Note that
//  [ L_1 ] [ U_1 U_2 ] = [ L_1 U_1 , L_1 U_2 ]
//  [ L_2 ]             = [ L_2 U_1 , L_2 U_2 ]
//  and that
//  [ L_11      ] [ u ] = [ L_11 u ]
//  [ L_21 L_22 ] [   ] = [ L_21 u ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  char diagL=(L_non_unit ? 'N' : 'U');
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  int m=L.size(0),k=L.size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<k;j++) { // L_1 U_1
    if (U_non_unit) {
      F77NAME(dcopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(dtrmv)('L','N',diagL,j+1,L.addr(),m,P->addr(0,j),1);
      if (j+1<k) {
        F77NAME(dgemv)('N',k-j-1,j+1,double_one_,L.addr(j+1,0),m,
          U.addr(0,j),1,double_zero_,P->addr(j+1,0),1);
      }
    } else {
      if (j>0) {
        F77NAME(dcopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(dtrmv)('L','N',diagL,j,L.addr(),m,P->addr(0,j),1);
        (*P)(j,j)=F77NAME(ddot)(j,L.addr(j,0),m,U.addr(0,j),1);
        if (j+1<k) {
          F77NAME(dgemv)('N',k-j-1,j,double_one_,L.addr(j+1,0),m,
            U.addr(0,j),1,double_zero_,P->addr(j+1,j),1);
        }
      }
      (*P)(j,j)+=(L_non_unit ? L(j,j) : double_one_);
      if (j+1<k) {
        F77NAME(daxpy)(k-j-1,double_one_,L.addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
    }
  }
  if (n>k) { // L_1 U_2
    F77NAME(dlacpy)('A',k,n-k,U.addr(0,k),k,P->addr(0,k),m);
    F77NAME(dtrmm)('L','L','N',diagL,k,n-k,double_one_,L.addr(),m,
      P->addr(0,k),m);
  }
  if (m>k) {
    char diagU=(U_non_unit ? 'N' : 'U');
    F77NAME(dlacpy)('A',m-k,k,L.addr(k,0),m,P->addr(k,0),m);
    F77NAME(dtrmm)('R','U','N',diagU,m-k,k,double_one_,U.addr(),k,
      P->addr(k,0),m); // L_2 U_1
    if (n>k) { // L_2 U_2
      F77NAME(dgemm)('N','N',m-k,n-k,k,double_one_,L.addr(k,0),m,
        U.addr(0,k),k,double_zero_,P->addr(k,k),m); 
    }
  }
  return P;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,
const UpperTrapezoidalMatrix<double,double> &U) {
// Note that
// M [ U_1 , U_2 ] = [ M U_1 , M U_2 ]
  char diagU=
    (dynamic_cast<const UnitUpperTrapezoidalMatrix<double,double>*>(&U)
    ==0 ? 'N' : 'U');
  int m=M.size(0),k=U.size(0),n=U.size(1);
  CHECK_SAME(k,M.size(1));
  Matrix<double,double> *P=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dlacpy)('A',m,k,M.addr(),m,P->addr(),m);
  F77NAME(dtrmm)('R','U','N',diagU,m,k,double_one_,U.addr(),k,
    P->addr(),m); // M U_1
  if (n>k) { // M U_2
    F77NAME(dgemm)('N','N',m,n-k,k,double_one_,M.addr(),m,U.addr(0,k),k,
      double_zero_,P->addr(0,k),m);
  }
  return P;
}

template class UpperTrapezoidalMatrix<double,double>;
template void testUpperTrapezoidalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double
UpperTriangularMatrix<double,double>::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  F77NAME(dtrcon)(norm,'U','N',n,addr(),n,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete [] work;
  return rcond;
}

/*
template<> UpperTriangularMatrix<double,double>*
UpperTriangularMatrix<double,double>::inverse() const {
  UpperTriangularMatrix<double,double> *I=
    OPERATOR_NEW UpperTriangularMatrix<double,double>(*this);
  int n=size(0);
  int info;
  F77NAME(dtrtri)('U','N',n,I->addr(),n,info);
  CHECK_TEST(info==0);
  return I;
}
*/

template class UpperTriangularMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void
UnitUpperTrapezoidalMatrix<double,double>::copyFrom(int m,int n,
const Matrix<double,double> &U) {
  m=min(m,min(size(0),U.size(0)));
  n=min(n,min(size(1),U.size(1)));
  for (int j=1;j<n;j++) {
    F77NAME(dcopy)(min(j,m),U.addr(0,j),1,addr(0,j),1);
  }
}

template<> UpperTrapezoidalMatrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<double,double> *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(m,n);
  if (U_non_unit) {
    S->copy(U);
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j,m),double_one_,addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=double_one_;
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(min(j,m),addr(0,j),1,S->addr(0,j),1);
        F77NAME(daxpy)(min(j,m),double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)=2.;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(m,double_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(j,addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)=double_one_;
      F77NAME(daxpy)(m-j,double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(j,addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)=2.;
      if (m<j-1) {
        F77NAME(daxpy)(m-j-1,double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> Matrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *S=OPERATOR_NEW Matrix<double,double>(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(daxpy)(min(j,m),double_one_,addr(0,j),1,S->addr(0,j),1);
    }
    if (j<m) (*S)(j,j)+=double_one_;
  }
  return S;
}

template<> UpperTrapezoidalMatrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<double,double> *D=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(m,n);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)=double_one_;
      F77NAME(daxpy)(min(j+1,m),double_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
        F77NAME(daxpy)(min(j,m),double_mone_,addr(0,j),1,D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)=double_zero_;
    }
  }
  return D;
}

template<> SquareMatrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *D=
    OPERATOR_NEW SquareMatrix<double,double>(m,double_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(j,addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)=double_one_;
      F77NAME(daxpy)(m-j,double_mone_,L.addr(j,j),1,D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(j,addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)=double_zero_;
      if (m<j-1) {
        F77NAME(daxpy)(m-j-1,double_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
    }
  }
  return D;
}

template<> Matrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *D=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(dcopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
    }
    if (j<m) (*D)(j,j)=double_one_;
    F77NAME(daxpy)(m,double_one_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> UpperTrapezoidalMatrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator*(double d) const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<double,double> *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(m,n);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(dcopy)(min(j,m),addr(0,j),1,P->addr(0,j),1);
      F77NAME(dscal)(min(j,m),d,P->addr(0,j),1);
    }
    if (j<m) (*P)(j,j)=d;
  }
  return P;
}

template<> UpperTrapezoidalMatrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator/(double d) const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<double,double> *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(m,n);
  double dinv=double_one_/d;
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(dcopy)(min(j,m),addr(0,j),1,P->addr(0,j),1);
      F77NAME(dscal)(min(j,m),dinv,P->addr(0,j),1);
    }
    if (j<m) (*P)(j,j)=dinv;
  }
  return P;
}

template<> UpperTrapezoidalMatrix<double,double>* 
UnitUpperTrapezoidalMatrix<double,double>::operator*(
const UpperTrapezoidalMatrix<double,double> &U) const {
// [ U_1 , U_2 ] [ V_11 , V_12 , V_13 ]
//               [   0  , V_22 , V_23 ]
//   = [ U_1 V_11 , U_1 V_12 + U_2 V_22 , U_1 V_13 + U_2 V_23 ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  int m=size(0),k=size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  UpperTrapezoidalMatrix<double,double> *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(m,n,double_zero_);
  if (U_non_unit) {
    for (int j=0;j<m;j++) { // U_1 V_11
      F77NAME(dcopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(dtrmv)('U','N','U',j+1,addr(),m,P->addr(0,j),1);
    }
    for (int j=m;j<k;j++) {
//    the next two lines could have been dlacpy & dtrmm, outside j loop
      F77NAME(dcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(dtrmv)('U','N','U',m,addr(),m,P->addr(0,j),1); // U_1 V_12
//    the next two lines could have been dtrmm('R',...)
//    but we could not replace both xxxmv with xxxmm
      F77NAME(dgemv)('N',m,j-m+1,double_one_,addr(0,m),m,U.addr(m,j),1,
        double_one_,P->addr(0,j),1); // U_2 V_22
    }
  } else {
    for (int j=0;j<m;j++) { // U_1 V_11
      if (j>0) {
        F77NAME(dcopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(dtrmv)('U','N','U',j,addr(),m,P->addr(0,j),1);
        F77NAME(daxpy)(j,double_one_,addr(0,j),1,P->addr(0,j),1);
      }
      (*P)(j,j)=double_one_;
    }
    for (int j=m;j<k;j++) {
      F77NAME(dcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(dtrmv)('U','N','U',m,addr(),m,P->addr(0,j),1); // U_1 V_12
      if (j>m) { // U_2 V_22
        F77NAME(dgemv)('N',m,j-m,double_one_,addr(0,m),m,U.addr(m,j),1,
          double_one_,P->addr(0,j),1);
      }
      F77NAME(daxpy)(m,double_one_,addr(0,j),1,P->addr(0,j),1);
    }
  }
  for (int j=k;j<n;j++) {
    F77NAME(dcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
    F77NAME(dtrmv)('U','N','U',m,addr(),m,P->addr(0,j),1); // U_1 V_13
    F77NAME(dgemv)('N',m,k-m,double_one_,addr(0,m),m,U.addr(m,j),1,
      double_one_,P->addr(0,j),1); // U_2 V_23
  }
  return P;
}

template<> Matrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  if (m>=n) {
//  note that
//  [ U_11 U_12 U_13 ] [ L_1 ] = [ U_11 L_1 + U_12 L_2 + U_13 L_3 ]
//  [      U_22 U_23 ] [ L_2 ] = [            U_22 L_2 + U_23 L_3 ]
//                     [ L_3 ]
    if (L_non_unit) { // U_11 L_1:
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(dgemv)('N',j,n-j,double_one_,addr(0,j),m,L.addr(j,j),1,
            double_zero_,P->addr(0,j),1);
        }
        F77NAME(dcopy)(n-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(dtrmv)('U','N','U',n-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(daxpy)(j,double_one_,addr(0,j),1,P->addr(0,j),1);
          if (j<n-1) {
            F77NAME(dgemv)('N',j,n-j-1,double_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,double_zero_,P->addr(0,j),1);
          }
        }
        if (j<n-1) {
          F77NAME(dcopy)(n-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(dtrmv)('U','N','U',n-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(ddot)(n-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=double_one_;
      }
    }
    if (n<k) { // U_12 L_2 + U_13 L_3:
      F77NAME(dgemm)('N','N',n,n,k-n,double_one_,addr(0,n),m,
        L.addr(n,0),k,double_one_,P->addr(0,0),m);
    }
    if (n<m) {
      F77NAME(dlacpy)('A',m-n,n,L.addr(n,0),k,P->addr(n,0),m);
      F77NAME(dtrmm)('L','U','N','U',m-n,n,double_one_,addr(n,n),m,
        P->addr(n,0),m); // U_22 L_2
      if (m<k) {
        F77NAME(dgemm)('N','N',m-n,n,k-m,double_one_,addr(n,m),m,
          L.addr(m,0),k,double_one_,P->addr(n,0),m); // U_23 L_3
      }
    }
  } else {
//  [ U_1 U_2 U_3 ] [ L_11   0  ]
//                  [ L_21 L_22 ]
//                  [ L_31 L_32 ]
//  = [ U_1 L_11 + U_2 L_21 _ U_3 L_31 , U_2 L_22 , U_2 L_22 + U_3 L_32 ]
    if (L_non_unit) { // U_1 L_11:
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(dgemv)('N',j,m-j,double_one_,addr(0,j),m,L.addr(j,j),1,
            double_zero_,P->addr(0,j),1);
        }
        F77NAME(dcopy)(m-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(dtrmv)('U','N','U',m-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(daxpy)(j,double_one_,addr(0,j),1,P->addr(0,j),1);
          if (j<m-1) {
            F77NAME(dgemv)('N',j,m-j-1,double_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,double_zero_,P->addr(0,j),1);
          }
        }
        if (j<m-1) {
          F77NAME(dcopy)(m-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(dtrmv)('U','N','U',m-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(ddot)(m-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=double_one_;
      }
    }
    if (m<k) { // U_2 L_21 + U_3 L_31 
      F77NAME(dgemm)('N','N',m,m,k-m,double_one_,addr(0,m),m,
        L.addr(m,0),k,double_one_,P->addr(0,0),m);
    }
    char diagL=(L_non_unit ? 'N' : 'U');
    F77NAME(dlacpy)('A',m,n-m,addr(0,m),m,P->addr(0,m),m);
    F77NAME(dtrmm)('R','L','U',diagL,m,n-m,double_one_,L.addr(m,m),k,
      P->addr(0,m),m); // U_2 L_22
    if (n<k) {
      F77NAME(dgemm)('N','N',m,n-m,k-n,double_one_,addr(0,n),m,
        L.addr(n,m),k,double_one_,P->addr(0,m),m); // U_3 L_32
    }
  }
  return P;
}

template<> Matrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
// [ U_1 U_2 ] [ M_1 ] = U_1 M_1 + U_2 M_2
//             [ M_2 ]
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<double,double> *P=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dlacpy)('A',m,n,M.addr(),k,P->addr(),m);
  F77NAME(dtrmm)('L','U','N','U',m,n,double_one_,addr(),m,
    P->addr(),m);
  if (k>m) {
    F77NAME(dgemm)('N','N',m,n,k-m,double_one_,addr(0,m),m,
      M.addr(m,0),k,double_one_,P->addr(),m);
  }
  return P;
}

template<> Vector<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<double,double> *p=OPERATOR_NEW Vector<double,double>(m);
  F77NAME(dcopy)(m,v.addr(),1,p->addr(),1);
  F77NAME(dtrmv)('U','N','U',m,addr(),m,p->addr(),1);
  F77NAME(dgemv)('N',m,n-m,double_one_,addr(0,m),m,v.addr(m),1,
    double_one_,p->addr(),1);
  return p;
}

/*
template<> UnitLowerTrapezoidalMatrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::transpose() const {
  int m=size(0),n=size(1);
  UnitLowerTrapezoidalMatrix<double,double> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<double,double>(n,m);
  for (int i=0;i<m;i++) {
    F77NAME(dcopy)(n-i-1,addr(i,i+1),m,L->addr(i+1,i),1);
  }
  return L;
}

template<> UnitLowerTrapezoidalMatrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::conjugateTranspose() const {
  return transpose();
}
*/

template<> Vector<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::trmv(
const Vector<double,double> &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<double,double> *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<double,double>(n);
    F77NAME(dcopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(dtrmv)('U','T','U',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(dgemv)('T',m,n-m,one_,addr(0,m),m,x.addr(),1,zero_,
        p->addr(m),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<double,double>(m);
    F77NAME(dcopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(dtrmv)('U','N','U',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(dgemv)('N',m,n-m,one_,addr(0,m),m,x.addr(m),1,one_,
        p->addr(),1);
    }
  }
  return p;
}

template<> Matrix<double,double>*
UnitUpperTrapezoidalMatrix<double,double>::trmm(
const Matrix<double,double> &M,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<double,double> *P=0;
  if (side=='L' || side=='l') {
    int k=M.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,M.size(0));
      P=OPERATOR_NEW Matrix<double,double>(n,k);
      F77NAME(dlacpy)('A',m,k,M.addr(),m,P->addr(),n);
      F77NAME(dtrmm)('L','U','T','U',m,k,one_,addr(),m,P->addr(),n);
      if (n>m) {
        F77NAME(dgemm)('T','N',n-m,k,m,one_,addr(0,m),m,M.addr(),m,
          zero_,P->addr(m,0),n);
      }
    } else {
      CHECK_SAME(n,M.size(0));
      P=OPERATOR_NEW Matrix<double,double>(m,k);
      F77NAME(dlacpy)('A',m,k,M.addr(),n,P->addr(),m);
      F77NAME(dtrmm)('L','U','N','U',m,k,one_,addr(),m,P->addr(),m);
      if (n>m) {
        F77NAME(dgemm)('N','N',m,k,n-m,one_,addr(0,m),m,M.addr(m,0),n,
          one_,P->addr(),m);
      }
    }
  } else {
    int k=M.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,M.size(1));
      P=OPERATOR_NEW Matrix<double,double>(k,m);
      F77NAME(dlacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(dtrmm)('R','U','T','U',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(dgemm)('N','T',k,m,n-m,one_,M.addr(0,m),k,addr(0,m),m,
          one_,P->addr(),k);
      }
    } else {
      CHECK_SAME(m,M.size(1));
      P=OPERATOR_NEW Matrix<double,double>(k,n);
      F77NAME(dlacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(dtrmm)('R','U','N','U',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(dgemm)('N','N',k,n-m,m,one_,M.addr(),k,addr(0,m),m,
          zero_,P->addr(0,m),k);
      }
    }
  }
  return P;
}

template<> double
UnitUpperTrapezoidalMatrix<double,double>::normFrobenius() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('F','U','U',m,size(1),addr(),m,work);
}

template<> double
UnitUpperTrapezoidalMatrix<double,double>::normInfinity() const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double result=F77NAME(dlantr)('I','U','U',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> double
UnitUpperTrapezoidalMatrix<double,double>::normMaxEntry() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('M','U','U',m,size(1),addr(),m,work);
}

template<> double
UnitUpperTrapezoidalMatrix<double,double>::normOne() const {
  int m=size(0);
  double *work=0;
  return F77NAME(dlantr)('O','U','U',m,size(1),addr(),m,work);
}

template<> void UnitUpperTrapezoidalMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    // calling routine will have to check consistency conditions
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(dcopy)(m,b.addr(),1,x.addr(),1);
    F77NAME(dtrsv)('U','T','U',m,addr(),m,x.addr(),1);
  } else {
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(dcopy)(m,b.addr(),1,x.addr(),1);
    if (n>m) { // use trailing entries of x as free variables
      F77NAME(dgemv)('N',m,n-m,mone_,addr(0,m),m,x.addr(m),1,one_,
        x.addr(),1);
    }
    F77NAME(dtrsv)('U','N','U',m,addr(),m,x.addr(),1);
  }
}

template<> void UnitUpperTrapezoidalMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<double,double>*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(dlacpy)('A',m,k,B.addr(),n,X.addr(),m);
      F77NAME(dtrsm)('L','U','T','U',m,k,one_,addr(),m,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(dlacpy)('A',m,k,B.addr(),m,X.addr(),n);
      if (n>m) { // use these entries of X as free variables
        F77NAME(dgemm)('N','N',m,k,n-m,mone_,addr(0,m),m,X.addr(m,0),n,
          one_,X.addr(),n);
      }
      F77NAME(dtrsm)('L','U','N','U',m,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(dlacpy)('A',k,m,B.addr(),k,X.addr(),k);
      if (n>m) { // use these entries of X as free variables
        F77NAME(dgemm)('N','T',k,m,n-m,mone_,X.addr(0,m),k,addr(0,m),m,
          one_,X.addr(),k);
      }
      F77NAME(dtrsm)('R','U','T','U',k,m,one_,addr(),m,X.addr(),k);
    } else { // calling routine will have to check consistency
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(dlacpy)('A',k,m,B.addr(),k,X.addr(),k);
      F77NAME(dtrsm)('R','U','N','U',k,m,one_,addr(),m,X.addr(),k);
    }
  }
}

template class UnitUpperTrapezoidalMatrix<double,double>;
template void testUnitUpperTrapezoidalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double
UnitUpperTriangularMatrix<double,double>::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  F77NAME(dtrcon)(norm,'U','U',n,addr(),n,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete [] work;
  return rcond;
}

/*
template<> UnitUpperTriangularMatrix<double,double>*
UnitUpperTriangularMatrix<double,double>::inverse() const {
  UnitUpperTriangularMatrix<double,double> *I=
    OPERATOR_NEW UnitUpperTriangularMatrix<double,double>(*this);
  int n=size(0);
  int info;
  F77NAME(dtrtri)('U','U',n,I->addr(),n,info);
  CHECK_TEST(info==0);
  return I;
}
*/

template class UnitUpperTriangularMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "OrthogonalMatrix.C"

/*
template<> OrthogonalMatrix<double,double>*
OrthogonalMatrix<double,double>::transpose() const {
  int m=size(0),n=size(1);
  OrthogonalMatrix<double,double> *Q=
    OPERATOR_NEW OrthogonalMatrix<double,double>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(m,addr(0,j),1,Q->addr(j,0),n);
  }
  return Q;
}

template<> OrthogonalMatrix<double,double>*
OrthogonalMatrix<double,double>::conjugateTranspose() const {
  return transpose();
}
*/

/*
template<> void OrthogonalMatrix<double,double>::solve(
const UpperTrapezoidalMatrix<double,double> &B,
Matrix<double,double> &X,char side,char trans) const {
  int m=size(0), n=size(1);
  bool B_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&B)==0);
  X=double_zero_;
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<=min(j,n-1);i++) {
            F77NAME(daxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<min(j,n);i++) {
            F77NAME(daxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
          if (j<n) {
            F77NAME(daxpy)(m,double_one_,addr(0,j),1,X.addr(0,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(ddot)(min(j+1,m),addr(0,i),1,B.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(ddot)(min(j,m),addr(0,i),1,B.addr(0,j),1);
            if (j<m) X(i,j)+=(*this)(j,i);
          }
        }
      }
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      if (B_non_unit) {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=F77NAME(ddot)(m-i,B.addr(i,i),k,addr(i,j),1);
          }
        }
      } else {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=(*this)(i,j)
              +F77NAME(ddot)(m-i-1,B.addr(i,i+1),k,addr(i+1,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      if (B_non_unit) {
        for (int i=0;i<k;i++) {
          for (int j=i;j<n;j++) {
            F77NAME(daxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      } else {
        for (int i=0;i<k;i++) {
          F77NAME(daxpy)(m,double_one_,addr(0,i),1,X.addr(i,0),k);
          for (int j=i+1;j<n;j++) {
            F77NAME(daxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      }
    }
  }
}
*/

/*
template<> void OrthogonalMatrix<double,double>::solve(
const LowerTrapezoidalMatrix<double,double> &B,
Matrix<double,double> &X,char side,char trans) const {
  int m=size(0), n=size(1);
  bool B_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&B)==0);
  X=double_zero_;
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=j;i<n;i++) {
            F77NAME(daxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          F77NAME(daxpy)(m,double_one_,addr(0,j),1,X.addr(0,j),1);
          for (int i=j+1;i<n;i++) {
            F77NAME(daxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(ddot)(m-j,addr(j,i),1,B.addr(j,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=(*this)(j,i)
              +F77NAME(ddot)(m-j-1,addr(j+1,i),1,B.addr(j+1,j),1);
          }
        }
      }
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      if (B_non_unit) {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=F77NAME(ddot)(min(i+1,m),B.addr(i,0),k,addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=F77NAME(ddot)(min(i,m),B.addr(i,0),k,addr(0,j),1);
            if (i<m) X(i,j)+=(*this)(i,j);
          }
        }
      }
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      if (B_non_unit) {
        for (int i=0;i<k;i++) {
          for (int j=0;j<=min(i,n-1);j++) {
            F77NAME(daxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      } else {
        for (int i=0;i<k;i++) {
          for (int j=0;j<min(i,n);j++) {
            F77NAME(daxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
          if (i<n) {
            F77NAME(daxpy)(m,double_one_,addr(0,i),1,X.addr(i,0),k);
          }
        }
      }
    }
  }
}
*/

template<> void OrthogonalMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  int m=size(0), n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(dgemm)('N','N',m,k,n,double_one_,addr(),m,B.addr(),n,
        double_zero_,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
//    F77NAME(dgemm)('T','N',n,k,m,double_one_,addr(),m,B.addr(),m,
//      double_zero_,X.addr(),n);
      Vector<double,double> *r=OPERATOR_NEW Vector<double,double>(m);
      for (int l=0;l<k;l++) {
        F77NAME(dcopy)(m,B.addr(0,l),1,r->addr(),1);
        for (int j=0;j<n;j++) {
          double &Xjl=X(j,l);
          Xjl=F77NAME(ddot)(m,r->addr(),1,addr(0,j),1); 
          F77NAME(daxpy)(m,-Xjl,addr(0,j),1,r->addr(),1);
        }
      }
      delete r; r=0;
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
//    F77NAME(dgemm)('N','N',k,n,m,double_one_,B.addr(),k,addr(),m,
//      double_zero_,X.addr(),k);
      Vector<double,double> *r=OPERATOR_NEW Vector<double,double>(m);
      for (int l=0;l<k;l++) {
        F77NAME(dcopy)(m,B.addr(l,0),k,r->addr(),1);
        for (int j=0;j<n;j++) {
          double &Xlj=X(l,j);
          Xlj=F77NAME(ddot)(m,r->addr(),1,addr(0,j),1);
          F77NAME(daxpy)(m,-Xlj,addr(0,j),1,r->addr(),1);
        }
      }
      delete r; r=0;
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(dgemm)('N','T',k,m,n,double_one_,B.addr(),k,addr(),m,
        double_zero_,X.addr(),k);
    }
  }
}

template<> void OrthogonalMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  int m=size(0), n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(dgemv)('N',m,n,double_one_,addr(),m,b.addr(),1,double_zero_,
      x.addr(),1);
  } else {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
//  F77NAME(dgemv)('T',m,n,double_one_,addr(),m,b.addr(),1,double_zero_,
//    x.addr(),1);
    Vector<double,double> *r=OPERATOR_NEW Vector<double,double>(m);
    r->copy(b);
    for (int j=0;j<n;j++) {
      x[j]=F77NAME(ddot)(m,r->addr(),1,addr(0,j),1); 
      F77NAME(daxpy)(m,-x[j],addr(0,j),1,r->addr(),1);
    }
    delete r; r=0;
  }
}

template class OrthogonalMatrix<double,double>;
template void testOrthogonalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SymmetricMatrix.C"

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,double> *M=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,addr(j,j),1,M->addr(j,j),1);
    if (j+1<n) {
      F77NAME(dcopy)(n-j-1,addr(j+1,j),1,M->addr(j,j+1),n);
    }
  }
  return M;
}

template<> void SymmetricMatrix<double,double>::fillWith(double d) {
  Matrix<double,double>::set('L',d,d);
}

template<> double
SymmetricMatrix<double,double>::operator()(int i,int j) const {
  return (j<=i ? *(this->addr(i,j)) : *(this->addr(j,i)) );
}

template<> SymmetricMatrix<double,double>&
SymmetricMatrix<double,double>::operator+=(
const SymmetricMatrix<double,double> &S) {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0))
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_one_,S.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<double,double>&
SymmetricMatrix<double,double>::operator-=(
const SymmetricMatrix<double,double> &S) {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0))
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,S.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<double,double>&
SymmetricMatrix<double,double>::operator*=(double scalar) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(n-j,scalar,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<double,double>&
SymmetricMatrix<double,double>::operator/=(double scalar) {
  int n=this->size(0);
  CHECK_NONZERO(scalar)
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(n-j,double_one_/scalar,addr(j,j),1);
  }
  return *this;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &U) const {
  int n=this->size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(dcopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(min(j+1,n),one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j,n),one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      if (j<n) (*S)(j,j)+=one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int n=this->size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(dcopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(n-j,one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(daxpy)(n-j-1,one_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      (*S)(j,j)+=one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator+(
const SquareMatrix<double,double> &S) const {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  T->copy(S);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,one_,addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n-j-1,one_,addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  int n=this->size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,one_,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n-j-1,one_,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &U) const {
  int n=this->size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(dcopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(min(j+1,n),mone_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j,n),mone_,U.addr(0,j),1,S->addr(0,j),1);
      }
      if (j<n) (*S)(j,j)+=mone_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const UnitUpperTrapezoidalMatrix<double,double> &U,
const SymmetricMatrix<double,double> &S) {
  int n=S.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(dcopy)(min(j,n),U.addr(0,j),1,T->addr(0,j),1);
    }
    if (j<n) (*T)(j,j)=double_one_;
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n-j-1,double_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<double,double>* operator-(
const UpperTrapezoidalMatrix<double,double> &U,
const SymmetricMatrix<double,double> &S) {
  int n=S.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+1,n),U.addr(0,j),1,T->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(min(j,n),U.addr(0,j),1,T->addr(0,j),1);
      }
      if (j<n) (*T)(j,j)=double_one_;
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n-j-1,double_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int n=this->size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(dcopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(n-j,mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(daxpy)(n-j-1,mone_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      (*S)(j,j)+=mone_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const UnitLowerTrapezoidalMatrix<double,double> &L,
const SymmetricMatrix<double,double> &S) {
  int n=S.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    if (j<n-1) {
      F77NAME(dcopy)(n-j-1,L.addr(j+1,j),1,T->addr(j+1,j),1);
    }
    (*T)(j,j)=double_one_;
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n-j-1,double_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<double,double>* operator-(
const LowerTrapezoidalMatrix<double,double> &L,
const SymmetricMatrix<double,double> &S) {
  int n=S.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(n-j,L.addr(j,j),1,T->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(dcopy)(n-j-1,L.addr(j+1,j),1,T->addr(j+1,j),1);
      }
      (*T)(j,j)=double_one_;
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n-j-1,double_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator-(
const SquareMatrix<double,double> &S) const {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(dcopy)(n-j-1,addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n,mone_,S.addr(0,j),1,T->addr(0,j),1);
  }
  return T;
}

template<> SquareMatrix<double,double>* operator-(
const SquareMatrix<double,double> &S,
const SymmetricMatrix<double,double> &SS) {
  int n=SS.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  T->copy(S);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,SS.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n-j-1,double_mone_,SS.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  int n=this->size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(dcopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n,mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const Matrix<double,double> &M,
const SymmetricMatrix<double,double> &S) {
  int n=S.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  T->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n-j-1,double_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator*(double d) const {
  int n=size(0);
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) { 
    const double *col_j=addr(j,j);
    double *S_col_j=S->addr(j,j);
    *S_col_j=(*col_j)*d;
    if (j<n-1) {
      S_col_j++;
      col_j++;
      double *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_col_j++,S_row_j+=n) {
        *S_row_j=*S_col_j=(*col_j)*d;
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator/(double d) const {
  CHECK_NONZERO(d);
  return operator*(double_one_/d);
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator*(
const SymmetricMatrix<double,double> &S) const {
// compute by bordering: note that
// [ sigma s^T ] [ tau t^T ] = [ sigma tau + s^T t , sigma t^T + s^T T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^T +  S  T ]
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int k=n-1;k>=0;k--) {
    if (k<n-1) {
      F77NAME(dsymv)('L',n-k-1,double_one_,addr(k+1,k+1),n,
        S.addr(k+1,k),1,double_zero_,T->addr(k+1,k),1); // S t
      F77NAME(dsymv)('L',n-k-1,double_one_,S.addr(k+1,k+1),n,
        addr(k+1,k),n,double_zero_,T->addr(k,k+1),n); // s^T T
      (*T)(k,k)=F77NAME(ddot)(n-k-1,addr(k+1,k),1,S.addr(k+1,k),1);//s^T t
    }
    F77NAME(dger)(n-k,n-k,double_one_,addr(k,k),1,S.addr(k,k),1,
      T->addr(k,k),n);
  }
  return T;
}

template<> Matrix<double,double>*
SymmetricMatrix<double,double>::operator*(
const UpperTrapezoidalMatrix<double,double> &U) const {
// compute by bordering: note that
// S [ U_1 , U_2 ] = [ S U_1 , S U_2 ]
// and that
// [ sigma s^T ] [ upsilon u^T ] = [ sigma upsilon , sigma u^T + s^T U ]
// [   s    S  ] [          U  ] [ [   s   upsilon ,   s   u^T +  S  U ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  char diag=(
    dynamic_cast<const UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0
    ? 'N' : 'U');
  for (int k=m-1;k>=0;k--) { // S U_1
    if (k<m-1) { // s^T U
      F77NAME(dcopy)(m-k-1,addr(k+1,k),1,M->addr(k,k+1),m);
      F77NAME(dtrmv)('U','T',diag,m-k-1,U.addr(k+1,k+1),m,
        M->addr(k,k+1),m);
    }
    if (diag=='N') {
      F77NAME(dger)(m-k,m-k,double_one_,addr(k,k),1,U.addr(k,k),m,
        M->addr(k,k),m);
    } else {
      F77NAME(daxpy)(m-k,double_one_,addr(k,k),1,M->addr(k,k),1);
      F77NAME(dger)(m-k,m-k-1,double_one_,addr(k,k),1,U.addr(k,k+1),m,
        M->addr(k,k+1),m);
    }
  }
  if (n>m) { // S U_2
    F77NAME(dsymm)('L','L',m,n-m,double_one_,addr(),m,U.addr(0,m),m,
      double_zero_,M->addr(0,m),m);
  }
  return M;
// compute by columns: note that
// S [ U_1 , U_2 ] = [ S U_1 , S U_2 ]
// and that
// [ S_11 S_21^T ] [ u ] = [ S_11 u ]
// [ S_21  S_22  ] [   ] = [ S_21 u ]
//bool U_non_unit=(dynamic_cast<const
//  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
//if (U_non_unit) { // S U_1
//  for (int j=0;j<m;j++) {
//    F77NAME(dsymv)('L',j+1,double_one_,addr(),m,U.addr(0,j),1,
//      double_zero_,M->addr(0,j),1); // S_11 u
//    if (j+1<m) {
//      F77NAME(dgemv)('N',m-j-1,j+1,double_one_,addr(j+1,0),m,
//        U.addr(0,j),1,double_zero_,M->addr(j+1,j),1); // S_21 u
//    }
//  }
//} else {
//  for (int j=0;j<m;j++) {
//    if (j>0) {
//      F77NAME(dsymv)('L',j,double_one_,addr(),m,U.addr(0,j),1,
//        double_zero_,M->addr(0,j),1);
//      F77NAME(daxpy)(j,double_one_,addr(j,0),m,M->addr(0,j),1);
//      (*M)(j,j)=F77NAME(ddot)(j,addr(j,0),m,U.addr(0,j),1);
//    }
//    (*M)(j,j)+=(*this)(j,j);
//    if (j+1<m) {
//      if (j>0) {
//        F77NAME(dgemv)('N',m-j-1,j,double_one_,addr(j+1,0),m,
//          U.addr(0,j),1,double_zero_,M->addr(j+1,j),1);
//      }
//      F77NAME(daxpy)(m-j-1,double_one_,addr(j+1,j),1,M->addr(j+1,j),1);
//    }
//  }
//}
//if (n>m) { // S U_2
//  F77NAME(dsymm)('L','L',m,n-m,double_one_,addr(),m,U.addr(0,m),m,
//    double_zero_,M->addr(0,m),m);
//}
}

template<> Matrix<double,double>* operator*(
const UpperTrapezoidalMatrix<double,double> &U,
const SymmetricMatrix<double,double> &S) {
//compute by bordering: note that
// [ U_1 U_2 ] [ S_11 S_21^T ]
//             [ S_21  S_22  ]
//   = [ U_1 S_11 + U_2 S_21 , U_1 S_21^T + U_2 S_22 ]
// and that
// [ U    u    ] [  S    s   ] = [ U S +    u    s^T , U s + u    sigma ]
// [   upsilon ] [ s^T sigma ] = [       upsilon s^T ,    upsilon sigma ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  char diag=(
    dynamic_cast<const UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0
    ? 'N' : 'U');
  for (int k=0;k<m;k++) { // U_1 S_11
    if (k>0) {
      F77NAME(dcopy)(k,S.addr(k,0),n,M->addr(0,k),1);
      F77NAME(dtrmv)('U','N',diag,k,U.addr(),m,M->addr(0,k),1); // U s
    }
    if (diag=='N') {
      F77NAME(dger)(k+1,k+1,double_one_,U.addr(0,k),1,S.addr(k,0),n,
        M->addr(),m);
    } else {
      F77NAME(dger)(k,k+1,double_one_,U.addr(0,k),1,S.addr(k,0),n,
        M->addr(),m);
      F77NAME(daxpy)(k+1,double_one_,S.addr(k,0),n,M->addr(k,0),m);
    }
  }
  if (n>m) {
    F77NAME(dgemm)('N','N',m,m,n-m,double_one_,U.addr(0,m),m,
      S.addr(m,0),n,double_one_,M->addr(),m); // U_2 S_21
    for (int j=m;j<n;j++) {
      F77NAME(dcopy)(m,S.addr(j,0),n,M->addr(0,j),1);
    }
    F77NAME(dtrmm)('L','U','N',diag,m,n-m,double_one_,U.addr(),m,
      M->addr(0,m),m); // U_1 S_21^T
    F77NAME(dsymm)('R','L',m,n-m,double_one_,S.addr(m,m),n,
      U.addr(0,m),m,double_one_,M->addr(0,m),m); // U_2 S_22
  }
  return M;
}

template<> Matrix<double,double>*
SymmetricMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &L) const {
// compute by bordering: note that
// [ S_11 S_21^T ] [ L_1 ] = [ S_11 L_1 + S_21^T L_2 ]
// [ S_21  S_22  ] [ L_2 ] = [ S_21 L_1 +  S_22  L_2 ]
// and that
// [  S    s   ] [   L          ] = [  S  L +   s   ell^T ,   s   lambda ]
// [ s^T sigma ] [ ell^T lambda ] = [ s^T L + sigma ell^T , sigma lambda ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(1));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0 ? 'N' : 'U');
  for (int k=0;k<n;k++) { // S_11 L_1
    if (k>0) {
      F77NAME(dcopy)(k,addr(k,0),m,M->addr(k,0),m);
      F77NAME(dtrmv)('L','T',diag,k,L.addr(),m,M->addr(k,0),m); // s^T L
    }
    if (diag=='N') {
      F77NAME(dger)(k+1,k+1,double_one_,addr(k,0),m,L.addr(k,0),m,
        M->addr(),m);
    } else {
      F77NAME(dger)(k+1,k,double_one_,addr(k,0),m,L.addr(k,0),m,
        M->addr(),m);
      F77NAME(daxpy)(k+1,double_one_,addr(k,0),m,M->addr(0,k),1);
    }
  }
  if (m>n) {
    F77NAME(dgemm)('T','N',n,n,m-n,double_one_,addr(n,0),m,
      L.addr(n,0),m,double_one_,M->addr(),m); // S_21^T L_2
    F77NAME(dlacpy)('A',m-n,n,addr(n,0),m,M->addr(n,0),m);
    F77NAME(dtrmm)('R','L','N',diag,m-n,n,double_one_,L.addr(),m,
      M->addr(n,0),m); // S_21 L_1
    F77NAME(dsymm)('L','L',m-n,n,double_one_,addr(n,n),m,L.addr(n,0),m,
      double_one_,M->addr(n,0),m); // S_22 L_2
  }
  return M;
// compute by columns: note that
// [ S_11 S_21^T ] [     ] = [ S_21^T ell ]
// [ S_21  S_22  ] [ ell ] = [  S_22  ell ]
//bool L_non_unit=(dynamic_cast<const
//  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
//if (L_non_unit) {
//  for (int j=0;j<n;j++) {
//    if (j>0) {
//      F77NAME(dgemv)('T',m-j,j,double_one_,addr(j,0),m,L.addr(j,j),1,
//        double_zero_,M->addr(0,j),1);
//    }
//    F77NAME(dsymv)('L',m-j,double_one_,addr(j,j),m,L.addr(j,j),1,
//      double_zero_,M->addr(j,j),1);
//  }
//} else {
//  for (int j=0;j<n;j++) {
//    if (j>0) {
//      if (j+1<m) {
//        F77NAME(dgemv)('T',m-j-1,j,double_one_,addr(j+1,0),m,
//          L.addr(j+1,j),1,double_zero_,M->addr(0,j),1);
//      }
//      F77NAME(daxpy)(j,double_one_,addr(j,0),m,M->addr(0,j),1);
//    }
//    (*M)(j,j)=F77NAME(ddot)(m-j-1,addr(j+1,j),1,L.addr(j+1,j),1);
//    F77NAME(dsymv)('L',m-j-1,double_one_,addr(j+1,j+1),m,
//      L.addr(j+1,j),1,double_zero_,M->addr(j+1,j),1);
//    (*M)(j,j)+=(*this)(j,j);
//    F77NAME(daxpy)(m-j-1,double_one_,addr(j+1,j),1,M->addr(j+1,j),1);
//  }
//}
//return M;
}

template<> Matrix<double,double>* operator*(
const LowerTrapezoidalMatrix<double,double> &L,
const SymmetricMatrix<double,double> &S) {
// compute by bordering: note that
// [ L_1 ] S = [ L_1 S ]
// [ L_2 ]   = [ L_2 S ]
// and that
// [ lambda   ] [ sigma s^T ] = [ lambda sigma       , lambda s^T       ]
// [  ell   L ] [   s    S  ] = [  ell   sigma + L s ,   ell  s^T + L S ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0 ? 'N' : 'U');
  for (int k=n-1;k>=0;k--) { // L_1 S
    if (k<n-1) {
      F77NAME(dcopy)(n-k-1,S.addr(k+1,k),1,M->addr(k+1,k),1);
      F77NAME(dtrmv)('L','N',diag,n-k-1,L.addr(k+1,k+1),m,
        M->addr(k+1,k),1); // L s
    }
    if (diag=='N') {
      F77NAME(dger)(n-k,n-k,double_one_,L.addr(k,k),1,S.addr(k,k),1,
        M->addr(k,k),m);
    } else {
      F77NAME(daxpy)(n-k,double_one_,S.addr(k,k),1,M->addr(k,k),m);
      F77NAME(dger)(n-k-1,n-k,double_one_,L.addr(k+1,k),1,S.addr(k,k),1,
        M->addr(k+1,k),m);
    }
  }
  if (m>n) { // L_2 S
    F77NAME(dsymm)('R','L',m-n,n,double_one_,S.addr(),n,L.addr(n,0),m,
      double_zero_,M->addr(n,0),m);
  }
  return M;
}

template<> SquareMatrix<double,double>*
SymmetricMatrix<double,double>::operator*(
const SquareMatrix<double,double> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  F77NAME(dsymm)('L','L',n,n,double_one_,addr(),n,S.addr(),n,
    double_zero_,T->addr(),n);
  return T;
}

template<> SquareMatrix<double,double>* operator*(
const SquareMatrix<double,double> &A,
const SymmetricMatrix<double,double> &B) {
// compute by bordering: note that
// [  M   v ] [  S    s   ] = [  M  S +  v s^T ,  M  s +  v sigma ]
// [ m^T mu ] [ s^T sigma ] = [ m^T S + mu s^T , m^T s + mu sigma ]
  int n=A.size(0);
  CHECK_SAME(n,B.size(0));
  SquareMatrix<double,double> *T=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
//for (int k=0;k<n;k++) {
//  if (k>0) {
//    F77NAME(dgemv)('N',k,k,double_one_,A.addr(),n,B.addr(k,0),n,
//      double_zero_,T->addr(0,k),1); // M s
//    F77NAME(dsymv)('L',k,double_one_,B.addr(),n,A.addr(k,0),n,
//      double_zero_,T->addr(k,0),n); // m^T S
//    (*T)(k,k)=F77NAME(ddot)(k,A.addr(k,0),n,B.addr(k,0),n); // m^T s
//  }
//  F77NAME(dger)(k+1,k+1,double_one_,A.addr(0,k),1,B.addr(k,0),n,
//    T->addr(),n);
//}
//return T;
// compute by bordering: note that
// [ mu v^T ] [ sigma s^T ] = [ mu sigma + v^T s , mu s^T + v^T S ]
// [  m  M  ] [   s    S  ] = [  m sigma +  M  s ,  m s^T +  M  S ]
  for (int k=n-1;k>=0;k--) {
    if (k+1<n) {
      F77NAME(dgemv)('N',n-k-1,n-k-1,double_one_,A.addr(k+1,k+1),n,
        B.addr(k+1,k),1,double_zero_,T->addr(k+1,k),1); // M s
      F77NAME(dsymv)('L',n-k-1,double_one_,B.addr(k+1,k+1),n,
        A.addr(k,k+1),n,double_zero_,T->addr(k,k+1),n);//v^T S = (S v)^T
      (*T)(k,k)=
        F77NAME(ddot)(n-k-1,A.addr(k,k+1),n,B.addr(k+1,k),1); // v^T s
    }
    F77NAME(dger)(n-k,n-k,double_one_,A.addr(k,k),1,B.addr(k,k),1,
      T->addr(k,k),n);
  }
  return T;
}

template<> Matrix<double,double>*
SymmetricMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  int m=size(0),n=M.size(1);
  CHECK_SAME(m,M.size(0));
  Matrix<double,double> *T=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dsymm)('L','L',m,n,double_one_,addr(),m,M.addr(),m,
    double_zero_,T->addr(),m);
  return T;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,const SymmetricMatrix<double,double> &S) {
  int m=M.size(0),n=S.size(0);
  CHECK_SAME(n,M.size(1));
  Matrix<double,double> *T=OPERATOR_NEW Matrix<double,double>(m,n);
  F77NAME(dsymm)('R','L',m,n,double_one_,S.addr(),n,M.addr(),m,
    double_zero_,T->addr(),m);
  return T;
}

template<> Vector<double,double>*
SymmetricMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<double,double> *w=OPERATOR_NEW Vector<double,double>(n);
  F77NAME(dsymv)('L',n,double_one_,addr(),n,v.addr(),1,double_zero_,
    w->addr(),1);
  return w;
}

// y := A * x * alpha + y * beta
template<> void SymmetricMatrix<double,double>::symv(double alpha,
const Vector<double,double> &x,double beta,Vector<double,double> &y)
const {
  int n=this->size(0);
  F77NAME(dsymv)('L',n,alpha,addr(),n,x.addr(),1,beta,y.addr(),1);
}

// A += x * alpha * x^T
template<> void SymmetricMatrix<double,double>::syr(double alpha,
const Vector<double,double> &x) {
  int n=this->size(0);
  F77NAME(dsyr)('L',n,alpha,x.addr(),1,addr(),n);
}

// A += x * alpha * y^T + y * alpha * x^T
template<> void SymmetricMatrix<double,double>::syr2(double alpha,
const Vector<double,double> &x,const Vector<double,double> &y) {
  int n=this->size(0);
  F77NAME(dsyr2)('L',n,alpha,x.addr(),1,y.addr(),1,addr(),n);
}

// C := A * alpha * B + C * beta
template<> void SymmetricMatrix<double,double>::symm(double alpha,
const Matrix<double,double> &B,double beta,Matrix<double,double> &C,
char side) const {
  int m=C.size(0),n=C.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,size(0));
    F77NAME(dsymm)('L','L',m,n,alpha,addr(),m,B.addr(),m,beta,
      C.addr(),m);
  } else {
    CHECK_SAME(n,size(0));
    F77NAME(dsymm)('R','L',m,n,alpha,addr(),n,B.addr(),m,beta,
      C.addr(),m);
  }
}

// C := A * alpha * A^T + C * beta
template<> void SymmetricMatrix<double,double>::syrk(double alpha,
const Matrix<double,double> &A,double beta,char transa) {
  int m=A.size(0),n=A.size(1);
  if (transa!='N' && transa!='n') {
    CHECK_SAME(n,size(0));
    F77NAME(dsyrk)('L',transa,n,m,alpha,A.addr(),m,beta,addr(),n);
  } else {
    CHECK_SAME(m,size(0));
    F77NAME(dsyrk)('L',transa,m,n,alpha,A.addr(),m,beta,addr(),m);
  }
}

// C := A * alpha * B^T + B * alpha + A^T + C * beta
template<> void SymmetricMatrix<double,double>::syr2k(double alpha,
const Matrix<double,double> &A,const Matrix<double,double> &B,
double beta,char transab) {
  int m=A.size(0),n=A.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (transab!='N' && transab!='n') {
    CHECK_SAME(n,size(0));
    F77NAME(dsyr2k)('L',transab,n,m,alpha,A.addr(),m,B.addr(),m,beta,
      addr(),n);
  } else {
    CHECK_SAME(m,size(0));
    F77NAME(dsyr2k)('L',transab,m,n,alpha,A.addr(),m,B.addr(),m,beta,
      addr(),m);
  }
}

/*
// y := abs(A) * abs(x) * alpha + abs(y) * beta
template<> void SymmetricMatrix<double,double>::syamv(double alpha,
const Vector<double,double> &x,double beta,Vector<double,double> &y)
const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,y.size());
  F77_NAME(dla_syamv)('L',n,alpha,addr(),n,x.addr(),1,beta,y.addr(),1);
}
*/

template<> double SymmetricMatrix<double,double>::equilibrate(
Vector<double,double> &s,double &scond) const {
  int n=size(0);
  CHECK_SAME(n,s.size());
  double amax=double_undefined_;
  int info;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  F77NAME(dsyequb)('L',n,addr(),n,s.addr(),scond,amax,work,info);
  CHECK_SAME(info,0);
  delete [] work;
  return amax;
}

template<> double SymmetricMatrix<double,double>::normFrobenius() const {
  int n=size(0);
  double *work=0;
  return F77NAME(dlansy)('F','L',n,addr(),n,work);
}

template<> double SymmetricMatrix<double,double>::normInfinity() const {
  int n=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,n);
  double val=F77NAME(dlansy)('I','L',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> double SymmetricMatrix<double,double>::normMaxEntry() const {
  int n=size(0);
  double *work=0;
  return F77NAME(dlansy)('M','L',n,addr(),n,work);
}

template<> double SymmetricMatrix<double,double>::normOne() const {
  int n=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,n);
  double val=F77NAME(dlansy)('O','L',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> double
SymmetricMatrix<double,double>::reciprocalConditionNumber() const {
  int n=size(0);

  SymmetricMatrix<double,double> *AF=
    OPERATOR_NEW SymmetricMatrix<double,double>(*this);
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  double w=numeric_limits<double>::infinity();
  F77NAME(dsytrf)('L',n,AF->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0);

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dsytrf)('L',n,AF->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;

  work=OPERATOR_NEW_BRACKET(double,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  double rcond;
  double anorm=F77NAME(dlansy)('O','L',n,addr(),n,work);
  F77NAME(dsycon)('L',n,AF->addr(),n,ipiv,anorm,rcond,work,iwork,info);
  delete [] ipiv; ipiv=0;
  delete [] work; work=0;
  delete [] iwork; iwork=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SymmetricMatrix<double,double>*
SymmetricMatrix<double,double>::inverse() const {
  int n=size(0);
  SymmetricMatrix<double,double> *Ainv=
    OPERATOR_NEW SymmetricMatrix<double,double>(*this);
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  double w=numeric_limits<double>::infinity();
  F77NAME(dsytrf)('L',n,Ainv->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0);

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dsytrf)('L',n,Ainv->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;

  work=OPERATOR_NEW_BRACKET(double,n);
  F77NAME(dsytri)('L',n,Ainv->addr(),n,ipiv,work,info);
  CHECK_SAME(info,0)

  delete [] ipiv;
  delete [] work;
  return Ainv;
}
*/

template<> Vector<double,double>*
SymmetricMatrix<double,double>::eigenvalues(
OrthogonalMatrix<double,double> *&Q) const {
  int n=size(0);
  if (Q!=0) CHECK_SAME(n,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Q->Matrix<double,double>::copyFrom('L',n,n,*this);
  Vector<double,double> *lambda =OPERATOR_NEW Vector<double,double>(n);
  double w;
  int lwork=-1,info;
  F77NAME(dsyev)(jobz,'L',n,Q->addr(),n,lambda->addr(),&w,lwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dsyev)(jobz,'L',n,Q->addr(),n,lambda->addr(),work,lwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;

  return lambda;
}

template<> void SymmetricMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char) const {
  int n=size(0);
  SymmetricMatrix<double,double> AF(*this);
  int info;
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  x.copy(b);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  double w=numeric_limits<double>::infinity();
  F77NAME(dsytrf)('L',n,AF.addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dsytrf)('L',n,AF.addr(),n,ipiv,work,lwork,info);
  delete [] work; work=0;

  if (info==0) {
    F77NAME(dsytrs)('L',n,1,AF.addr(),n,ipiv,x.addr(),n,info);
  }
  CHECK_SAME(info,0)
  delete [] ipiv; ipiv=0;
}

template<> void SymmetricMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,char)
const {
  int n=size(0);
  SymmetricMatrix<double,double> AF(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  double w=numeric_limits<double>::infinity();
  int info;
  F77NAME(dsytrf)('L',n,AF.addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dsytrf)('L',n,AF.addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;

  bool left_side=(side=='L' || side=='l');
  if (left_side) {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1))
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    X.copy(B);
    F77NAME(dsytrs)('L',n,1,AF.addr(),n,ipiv,X.addr(),n,info);
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0))
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    double *t=OPERATOR_NEW_BRACKET(double,n);
    for (int i=0;i<k;i++) {
      F77NAME(dcopy)(n,B.addr(i,0),k,t,1);
      F77NAME(dsytrs)('L',n,1,AF.addr(),n,ipiv,t,n,info);
      F77NAME(dcopy)(n,t,1,X.addr(i,0),k);
      CHECK_SAME(info,0)
    }
    delete [] t; t=0;
  }
  delete [] ipiv; ipiv=0;
}

template class SymmetricMatrix<double,double>;
template void testSymmetricMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveMatrix<double,double>&
SymmetricPositiveMatrix<double,double>::operator*=(double scalar) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(n-j,abs(scalar),addr(j,j),1);
  }
  return *this;
}

template<> SymmetricPositiveMatrix<double,double>&
SymmetricPositiveMatrix<double,double>::operator/=(double scalar) {
  int n=this->size(0);
  CHECK_NONZERO(scalar)
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(n-j,double_one_/abs(scalar),addr(j,j),1);
  }
  return *this;
}

template<> double SymmetricPositiveMatrix<double,double>::equilibrate(
Vector<double,double> &s,double &scond) const {
  int n=size(0);
  CHECK_SAME(n,s.size());
  double amax=double_undefined_;
  int info;
  F77NAME(dpoequb)(n,addr(),n,s.addr(),scond,amax,info);
  CHECK_SAME(info,0);
  return amax;
}

template<> double
SymmetricPositiveMatrix<double,double>::reciprocalConditionNumber(
) const {
  int n=size(0);
  SymmetricPositiveMatrix<double,double> *AF=
    OPERATOR_NEW SymmetricPositiveMatrix<double,double>(*this);
  int info;
  F77NAME(dpotrf)('L',n,AF->addr(),n,info);
  CHECK_SAME(info,0)

  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  double rcond;
  double anorm=F77NAME(dlansy)('O','L',n,addr(),n,work);
  F77NAME(dpocon)('L',n,AF->addr(),n,anorm,rcond,work,iwork,info);
  delete [] work; work=0;
  delete [] iwork; iwork=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SymmetricPositiveMatrix<double,double>*
SymmetricPositiveMatrix<double,double>::inverse() const {
  int n=size(0);
  SymmetricPositiveMatrix<double,double> *Ainv=
    OPERATOR_NEW SymmetricPositiveMatrix<double,double>(*this);
  int info;
  F77NAME(dpotrf)('L',n,Ainv->addr(),n,info);
  CHECK_SAME(info,0);

  F77NAME(dpotri)('L',n,Ainv->addr(),n,info);
  CHECK_SAME(info,0)
  return Ainv;
}
*/

template<> void SymmetricPositiveMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char) const {
  int n=size(0);
  SymmetricPositiveMatrix<double,double> AF(*this);
  int info;
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  x.copy(b);
  F77NAME(dpotrf)('L',n,AF.addr(),n,info);
  CHECK_SAME(info,0)

  F77NAME(dpotrs)('L',n,1,AF.addr(),n,x.addr(),n,info);
  CHECK_SAME(info,0)
}

template<> void SymmetricPositiveMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,char)
const {
  int n=size(0);
  SymmetricPositiveMatrix<double,double> AF(*this);
  int info;
  F77NAME(dpotrf)('L',n,AF.addr(),n,info);
  CHECK_SAME(info,0)

  bool left_side=(side=='L' || side=='l');
  if (left_side) {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1))
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    X.copy(B);
    F77NAME(dpotrs)('L',n,1,AF.addr(),n,X.addr(),n,info);
    CHECK_SAME(info,0)
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0))
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    double *t=OPERATOR_NEW_BRACKET(double,n);
    for (int i=0;i<k;i++) {
      F77NAME(dcopy)(n,B.addr(i,0),k,t,1);
      F77NAME(dpotrs)('L',n,1,AF.addr(),n,t,n,info);
      F77NAME(dcopy)(n,t,1,X.addr(i,0),k);
      CHECK_SAME(info,0)
    }
    delete [] t; t=0;
  }
}

template class SymmetricPositiveMatrix<double,double>;
template void testSymmetricPositiveMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "BandMatrix.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double TridiagonalMatrix<double,double>::safety_ =
  double_zero_;
template<> const double TridiagonalMatrix<double,double>::outofbounds_ =
  double_zero_;
template<> const double TridiagonalMatrix<double,double>::undefined_ =
  numeric_limits<double>::infinity();

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,double> *M=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  F77NAME(dcopy)(n,D->addr(),1,M->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,M->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,U->addr(),1,M->addr(0,1),n+1);
  return M;
}

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::operator+(
const SymmetricMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int k=0;k<n;k++) {
    F77NAME(dcopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) F77NAME(dcopy)(n-k-1,M.addr(k+1,k),1,S->addr(k,k+1),n);
  }
  F77NAME(daxpy)(n,double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(dcopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)=double_one_;
      if (k<n-1) F77NAME(dcopy)(n-k-1,M.addr(k+1,k),1,S->addr(k+1,k),1);
    }
  }
  F77NAME(daxpy)(n,double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<double,double>*
TridiagonalMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<double,double> *S=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(dcopy)(k+1,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(dcopy)(k,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)=double_one_;
    }
  }
  F77NAME(daxpy)(n,double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  S->copy(M);
  F77NAME(daxpy)(n,double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::operator-(
const SymmetricMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    F77NAME(daxpy)(n-k,double_mone_,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) {
      F77NAME(daxpy)(n-k-1,double_mone_,M.addr(k+1,k),1,S->addr(k,k+1),n);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const SymmetricMatrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) { 
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int k=0;k<n;k++) {
    F77NAME(dcopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) F77NAME(dcopy)(n-k-1,M.addr(k+1,k),1,S->addr(k,k+1),n);
  }
  F77NAME(daxpy)(n,double_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0);
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(daxpy)(n-k,double_mone_,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)-=double_one_;
      if (k<n-1) {
        F77NAME(daxpy)(n-k-1,double_mone_,M.addr(k+1,k),1,
          S->addr(k+1,k),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const LowerTrapezoidalMatrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) { 
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(dcopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)=double_one_;
      if (k<n-1) {
        F77NAME(dcopy)(n-k-1,M.addr(k+1,k),1,S->addr(k+1,k),1);
      }
    }
  }
  F77NAME(daxpy)(n,double_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<double,double>*
TridiagonalMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<double,double> *S=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&M)==0);
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(daxpy)(k+1,double_mone_,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(daxpy)(k,double_mone_,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)-=double_one_;
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<double,double>* operator-(
const UpperTrapezoidalMatrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<double,double> *S=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(dcopy)(k+1,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(dcopy)(k,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)=double_one_;
    }
  }
  F77NAME(daxpy)(n,double_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  *S-=M;
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const Matrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  S->copy(M);
  F77NAME(daxpy)(n,double_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::operator*(
const SymmetricMatrix<double,double> &M) const {
// compute by bordering: note that
// [    tau     upsilon e_0^T ] [ sigma s^T ]
// [ e_0 lambda       T       ] [   s    S  ]
//   = [ tau sigma + upsilon e_0^T s , tau s^T + upsilon e_0^T S ]
//   = [ e_0 lambda sigma +      T s , e_0 lambda s^T + T S      ]
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  (*S)(n-2,n-2)=(*D)[n-2]*M(n-2,n-2)+(*U)[n-2]*M(n-1,n-2);
  (*S)(n-1,n-2)=(*L)[n-2]*M(n-2,n-2)+(*D)[n-1]*M(n-1,n-2);
  (*S)(n-2,n-1)=(*D)[n-2]*M(n-1,n-2)+(*U)[n-2]*M(n-1,n-1);
  (*S)(n-1,n-1)=(*L)[n-2]*M(n-1,n-2)+(*D)[n-1]*M(n-1,n-1);
  for (int k=n-3;k>=0;k--) {
    F77NAME(dgtmv)(n-k-1,double_one_,L->addr(k+1),D->addr(k+1),
      U->addr(k+1),M.addr(k+1,k),1,double_zero_,S->addr(k+1,k),1); // T s
    F77NAME(daxpy)(n-k-1,(*U)[k],M.addr(k+1,k+1),1,S->addr(k,k+1),n);
      // upsilon e_0^T S
    (*S)(k,k)=(*U)[k]*M(k+1,k); // upsilon e_0^T s
    F77NAME(daxpy)(n-k,(*D)[k],M.addr(k,k),1,S->addr(k,k),n);
      // tau [ sigma , s^T ]
    F77NAME(daxpy)(n-k,(*L)[k],M.addr(k,k),1,S->addr(k+1,k),n);
      // e_0 lambda [ sigma , s^T ]
  }
  return S;
}

template<> SquareMatrix<double,double>* operator*(
const SymmetricMatrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) {
// compute by bordering: note that
// [ sigma s^T ] [    tau     upsilon e_0^T ]
// [   s    S  ] [ e_0 lambda       T       ]
//   = [ sigma tau + s^T e_0 lambda , sigma upsilon e_0^T + s^T T ]
//   = [     s tau +   S e_0 lambda ,     s upsilon e_0^T +   S T ]
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  (*S)(n-2,n-2)=M(n-2,n-2)*T(n-2,n-2)+M(n-1,n-2)*T(n-1,n-2);
  (*S)(n-1,n-2)=M(n-1,n-2)*T(n-2,n-2)+M(n-1,n-1)*T(n-1,n-2);
  (*S)(n-2,n-1)=M(n-2,n-2)*T(n-2,n-1)+M(n-1,n-2)*T(n-1,n-1);
  (*S)(n-1,n-1)=M(n-1,n-2)*T(n-2,n-1)+M(n-1,n-1)*T(n-1,n-1);
  for (int k=n-3;k>=0;k--) {
    F77NAME(dgtmv)(n-k-1,double_one_,T.addr(k+1,k+2),T.addr(k+1,k+1),
      T.addr(k+2,k+1),M.addr(k+1,k),1,double_zero_,S->addr(k,k+1),n);
    F77NAME(daxpy)(n-k-1,T(k+1,k),M.addr(k+1,k+1),1,S->addr(k+1,k),1);
      // S e_0 lambda
    (*S)(k,k)=M(k+1,k)*T(k+1,k); // s^T e_0 lambda
    F77NAME(daxpy)(n-k,T(k,k),M.addr(k,k),1,S->addr(k,k),1);
      // [ sigma ] tau
      // [   s   ]
    F77NAME(daxpy)(n-k,T(k,k+1),M.addr(k,k),1,S->addr(k,k+1),1);
      // [ sigma ] upsilon
      // [   s   ]
  }
  return S;
}

template<> Matrix<double,double>*
TridiagonalMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &M) const {
// compute by columns: note that
// [      T_11       ,e_j upsilon e_0^T ][ 0 ] = [ e_j upsilon e_0^T ell ]
// [ e_0 lambda e_j^T,      T_22        ][ell] = [             T_22  ell ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0);
  if (M_non_unit) {
    F77NAME(dgtmv)(m,double_one_,L->addr(),D->addr(),U->addr(),
      M.addr(),1,double_zero_,S->addr(),1);
    for (int j=1;j<n;j++) {
      (*S)(j-1,j)=(*U)[j-1]*M(j,j); // e_j upsilon e_0^T ell
      if (j<m-1) { // T_22 ell
        F77NAME(dgtmv)(m-j,double_one_,L->addr(j),D->addr(j),
          U->addr(j),M.addr(j,j),1,double_one_,S->addr(j,j),1);
      } else (*S)(j,j)=(*D)[j]*M(j,j);
    }
  } else {
// note that
// [     tau    , upsilon e_0^T ] [ 1 ] = [ tau + upsilon e_0^T ell ]
// [ e_0 lambda ,        T      ] [ell] = [ e_0 lambda      + T ell ]
    (*S)(0,0)=(*D)[0]+(*U)[0]*M(1,0);
    (*S)(1,0)=(*L)[0];
    F77NAME(dgtmv)(m-1,double_one_,L->addr(1),D->addr(1),
      U->addr(1),M.addr(1,0),1,double_one_,S->addr(1,0),1);
    for (int j=1;j<n;j++) {
      (*S)(j-1,j)=(*U)[j-1]; // e_j upsilon e_0^T ell
      (*S)(j,j)=(*D)[j]; // tau
      if (j<m-1) {
        (*S)(j,j)+=(*U)[j]*M(j+1,j); // upsilon e_0^T ell
        (*S)(j+1,j)=(*L)[j]; // e_0 lambda
      }
      if (j<m-2) {
        F77NAME(dgtmv)(m-j-1,double_one_,L->addr(j+1),
          D->addr(j+1),U->addr(j+1),M.addr(j+1,j),1,double_one_,
          S->addr(j+1,j),1);
      } else if (j<m-1) (*S)(j+1,j)+=(*D)[j+1]*M(j+1,j);
    }
  }
  return S;
}

template<> Matrix<double,double>* operator*(
const LowerTrapezoidalMatrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) {
// compute by columns: note that
// L T e_j = L e_{j-1} T_{j-1,j} + L e_j T_{j,j} + L e_{j+1{ T){j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(m-j+1,T(j-1,j),M.addr(j-1,j-1),1,
          S->addr(j-1,j),1);
      }
      F77NAME(daxpy)(m-j,T(j,j),M.addr(j,j),1,S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(m-j-1,T(j+1,j),M.addr(j+1,j+1),1,S->addr(j+1,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        (*S)(j-1,j)=T(j-1,j);
        F77NAME(daxpy)(m-j,T(j-1,j),M.addr(j,j-1),1,S->addr(j,j),1);
      }
      (*S)(j,j)+=T(j,j);
      if (j<m-1) {
        F77NAME(daxpy)(m-j-1,T(j,j),M.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      if (j<n-1) {
        (*S)(j+1,j)+=T(j+1,j);
        if (j<m-2) {
          F77NAME(daxpy)(m-j-2,T(j+1,j),M.addr(j+2,j+1),1,
            S->addr(j+2,j),1);
        }
      }
    }
  }
  return S;
}

template<> Matrix<double,double>*
TridiagonalMatrix<double,double>::operator*(
const UpperTrapezoidalMatrix<double,double> &M) const {
// compute by columns: note that
// T [ U_1 , U_2 ] = [ T U_1 , T U_2 ]
// and that
// [      T_11       ,e_j upsilon e_0^T ][ u ] = [             T_11 u ]
// [ e_0 lambda e_j^T,      T_22        ][ 0 ] = [ e_0 lambda e_j^T u ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&M)==0);
  if (M_non_unit) {
    (*S)(0,0)=(*D)[0]*M(0,0);
    (*S)(1,0)=(*L)[0]*M(0,0);
    for (int j=1;j<m-1;j++) { // T U_1
      F77NAME(dgtmv)(j+1,double_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,double_zero_,S->addr(0,j),1); // T_11 u
      (*S)(j+1,j)=(*L)[j]*M(j,j); // e_0 lambda e_j^T u
    }
    for (int j=m-1;j<n;j++) { // T U_2
      F77NAME(dgtmv)(m,double_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,double_zero_,S->addr(0,j),1);
    }
  } else {
    (*S)(0,0)=(*D)[0];
    (*S)(1,0)=(*L)[0];
    for (int j=1;j<m;j++) {
      F77NAME(dgtmv)(j,double_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,double_zero_,S->addr(0,j),1);
      (*S)(j-1,j)+=(*U)[j-1];
      (*S)(j,j)=(*L)[j-1]*M(j-1,j)+(*D)[j];
      if (j<m-1) (*S)(j+1,j)=(*L)[j];
    }
    for (int j=m;j<n;j++) { // T U_2
      F77NAME(dgtmv)(m,double_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,double_zero_,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<double,double>* operator*(
const UpperTrapezoidalMatrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) {
// compute by columns: note that
// U T e_j = U e_{j-1} T_{j-1,j} + U e_j T_{j,j} + U e_{j+1{ T){j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j,m),T(j-1,j),M.addr(0,j-1),1,S->addr(0,j),1);
      }
      F77NAME(daxpy)(min(j+1,m),T(j,j),M.addr(0,j),1,S->addr(0,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(min(j+2,m),T(j+1,j),M.addr(0,j+1),1,
          S->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j-1,m),T(j-1,j),M.addr(0,j-1),1,
          S->addr(0,j),1);
        if (j<=m) (*S)(j-1,j)=T(j-1,j);
      }
      F77NAME(daxpy)(min(j,m),T(j,j),M.addr(0,j),1,S->addr(0,j),1);
      if (j<m) (*S)(j,j)+=T(j,j);
      if (j<n-1) {
        F77NAME(daxpy)(min(j+1,m),T(j+1,j),M.addr(0,j+1),1,
          S->addr(0,j),1);
        if (j<m-1) (*S)(j+1,j)+=T(j+1,j);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
TridiagonalMatrix<double,double>::operator*(
const SquareMatrix<double,double> &M) const {
  int m=M.size(0);
  CHECK_SAME(m,size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(m,double_zero_);
  for (int j=0;j<m;j++) {
    F77NAME(dgtmv)(m,double_one_,L->addr(),D->addr(),U->addr(),
      M.addr(0,j),1,double_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator*(
const SquareMatrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) {
  int m=M.size(0);
  CHECK_SAME(m,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(m,double_zero_);
  for (int j=0;j<m;j++) {
    for (int k=max(0,j-1);k<=min(m-1,j+1);k++) {
      F77NAME(daxpy)(m,T(k,j),M.addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<double,double>*
TridiagonalMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dgtmv)(m,double_one_,L->addr(),D->addr(),U->addr(),
      M.addr(0,j),1,double_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,
const TridiagonalMatrix<double,double> &T) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
      F77NAME(daxpy)(m,T(k,j),M.addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Vector<double,double>*
TridiagonalMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<double,double> *w=
    OPERATOR_NEW Vector<double,double>(n,double_zero_);
  F77NAME(dgtmv)(n,double_one_,L->addr(),D->addr(),U->addr(),
    v.addr(),1,double_zero_,w->addr(),1);
  return w;
}

template<> double TridiagonalMatrix<double,double>::normFrobenius()
const {
  return F77NAME(dlangt)('F',dim,L->addr(),D->addr(),U->addr());
}

template<> double TridiagonalMatrix<double,double>::normInfinity()
const {
  return F77NAME(dlangt)('I',dim,L->addr(),D->addr(),U->addr());
}

template<> double TridiagonalMatrix<double,double>::normMaxEntry()
const {
  return F77NAME(dlangt)('M',dim,L->addr(),D->addr(),U->addr());
}

template<> double TridiagonalMatrix<double,double>::normOne() const {
  return F77NAME(dlangt)('O',dim,L->addr(),D->addr(),U->addr());
}

template<> double
TridiagonalMatrix<double,double>::reciprocalConditionNumber(char norm)
const {
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  Vector<double,double> *UF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *UF2=
    OPERATOR_NEW Vector<double,double>(max(0,dim-2));
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(dgttrf)(dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),ipiv,
    info);
  double rcond=undefined_;
  double *work=OPERATOR_NEW_BRACKET(double,2*dim);
  int *iwork=OPERATOR_NEW_BRACKET(int,dim);
  double anorm=F77NAME(dlangt)(norm,dim,L->addr(),D->addr(),U->addr());
  F77NAME(dgtcon)(norm,dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
    ipiv,anorm,rcond,work,iwork,info);
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  delete [] ipiv; ipiv=0;
  delete UF2; UF2=0;
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
  return rcond;
}

template<> void TridiagonalMatrix<double,double>::gtmv(double alpha,
const Vector<double,double> &x,double beta,Vector<double,double> &b,
char trans) const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  if (trans=='N' || trans=='n') {
    F77NAME(dgtmv)(n,alpha,L->addr(),D->addr(),U->addr(),
      x.addr(),1,beta,b.addr(),1);
  } else {
    F77NAME(dgtmv)(n,alpha,U->addr(),D->addr(),L->addr(),
      x.addr(),1,beta,b.addr(),1);
  }
}

template<> void TridiagonalMatrix<double,double>::gtmm(double alpha,
const Matrix<double,double> &X,double beta,Matrix<double,double> &B,
char side,char trans) const {
  if (side=='L' || side=='l') {
    int m=size(0),n=X.size(1);
    CHECK_SAME(n,B.size(1));
    CHECK_SAME(m,X.size(0));
    CHECK_SAME(m,B.size(0));
    if (trans=='N' || trans=='n') {
      for (int j=0;j<n;j++) {
        F77NAME(dgtmv)(m,alpha,L->addr(),D->addr(),U->addr(),
          X.addr(0,j),1,beta,B.addr(0,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        F77NAME(dgtmv)(m,alpha,U->addr(),D->addr(),L->addr(),
          X.addr(0,j),1,beta,B.addr(0,j),1);
      }
    }
  } else {
    int m=X.size(0),n=X.size(1);
    CHECK_SAME(n,size(0));
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,B.size(1));
    if (abs(beta)==double_zero_) B=double_zero_;
    else if (beta!=double_one_) B*=beta;
    if (trans=='N' || trans=='n') {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(daxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(daxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

template<> void TridiagonalMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  Vector<double,double> *UF=OPERATOR_NEW Vector<double,double>(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  x.copy(b);
  int info;
  if (trans!='N' && trans!='n') {
    F77NAME(dgtsv)(dim,1,UF->addr(),DF->addr(),LF->addr(),x.addr(),dim,
      info);
  } else {
    F77NAME(dgtsv)(dim,1,LF->addr(),DF->addr(),UF->addr(),x.addr(),dim,
      info);
  }
  CHECK_SAME(info,0);
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
}

template<> void TridiagonalMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  Vector<double,double> *UF=OPERATOR_NEW Vector<double,double>(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  int info;
  X.copy(B);
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    if (trans!='N' && trans!='n') {
      F77NAME(dgtsv)(m,n,UF->addr(),DF->addr(),LF->addr(),X.addr(),m,
        info);
    } else {
      F77NAME(dgtsv)(m,n,LF->addr(),DF->addr(),UF->addr(),X.addr(),m,
        info);
    }
  } else {
    CHECK_SAME(dim,n);
    Vector<double,double> *UF2=OPERATOR_NEW Vector<double,double>(dim-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
    F77NAME(dgttrf)(n,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
      ipiv,info);
    CHECK_TEST(info==0);
    double *t=OPERATOR_NEW_BRACKET(double,n);
    char rtrans=(trans!='N' && trans!='n' ? 'N' : 'T');
    for (int i=0;i<m;i++) {
      F77NAME(dcopy)(n,B.addr(i,0),m,t,1);
      F77NAME(dgttrs)(rtrans,n,1,LF->addr(),DF->addr(),UF->addr(),
        UF2->addr(),ipiv,t,n,info);
      F77NAME(dcopy)(n,t,1,X.addr(i,0),m);
    }
    delete [] ipiv; ipiv=0;
    delete UF2; UF2=0;
    delete [] t; t=0;
  }
  CHECK_SAME(info,0);
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
}

template class TridiagonalMatrix<double,double>;
template void testTridiagonalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double
  SymmetricTridiagonalMatrix<double,double>::safety_=double_zero_;
template<> const double
  SymmetricTridiagonalMatrix<double,double>::outofbounds_=
  double_zero_;
template<> const double
  SymmetricTridiagonalMatrix<double,double>::undefined_=
  numeric_limits<double>::infinity();

template<> SymmetricTridiagonalMatrix<double,double>::
SymmetricTridiagonalMatrix(int n,double d) : dim(n) {
  L=OPERATOR_NEW Vector<double,double>(n-1,d);
  D=OPERATOR_NEW Vector<double,double>(n,d);
}

template<> SquareMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,double> *M=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  F77NAME(dcopy)(n,D->addr(),1,M->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,M->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,M->addr(0,1),n+1);
  return M;
}

template<> void SymmetricTridiagonalMatrix<double,double>::fillWith(
double scalar) {
  *L=scalar; *D=scalar;
}

template<> double
SymmetricTridiagonalMatrix<double,double>::upperDiagonalValue(int i)
const {
  return (*L)[i];
}

template<> SymmetricTridiagonalMatrix<double,double>&
SymmetricTridiagonalMatrix<double,double>::operator=(double scalar) {
  *L=scalar; *D=scalar; return *this;
}

template<> SymmetricTridiagonalMatrix<double,double>&
SymmetricTridiagonalMatrix<double,double>::operator*=(double scalar) {
  int n=this->size(0);
  (*D)*=scalar;
  (*L)*=scalar;
  return *this;
}

template<> SymmetricTridiagonalMatrix<double,double>&
SymmetricTridiagonalMatrix<double,double>::operator/=(double scalar) {
  int n=this->size(0);
  CHECK_NONZERO(scalar)
  (*D)/=scalar;
  (*L)/=scalar;
  return *this;
}

template<> TridiagonalMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator*(double d) const {
  int n=size(0);
  TridiagonalMatrix<double,double> *S=
    OPERATOR_NEW TridiagonalMatrix<double,double>(n);
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(0,1),1);
  *S*=d;
  return S;
}

template<> TridiagonalMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator/(double d) const {
  int n=size(0);
  TridiagonalMatrix<double,double> *S=
    OPERATOR_NEW TridiagonalMatrix<double,double>(n);
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(0,1),1);
  *S/=d;
  return S;
}

template<> TridiagonalMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator+(
const TridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,double> *S=
    OPERATOR_NEW TridiagonalMatrix<double,double>(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_one_,D->addr(),1,S->addr(0,0),1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(1,0),1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(0,1),1);
  return S;
}

template<> SymmetricMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator+(
const SymmetricMatrix<double,double> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricMatrix<double,double> *T=
    OPERATOR_NEW SymmetricMatrix<double,double>(n);
  T->copy(S);
  F77NAME(daxpy)(n,double_one_,D->addr(),1,T->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,T->addr(1,0),n+1);
  return T;
}

template<> SquareMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0); 
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=double_one_;
      if (j<n-1) {
        F77NAME(dcopy)(n-j-1,M.addr(j+1,j),1,S->addr(j+1,j),1);
      }
    }
  }
  F77NAME(daxpy)(n,double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0); 
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(j+1,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(j,U.addr(0,j),1,H->addr(0,j),1);
      }
      (*H)(j,j)=double_one_;
    }
  }
  F77NAME(daxpy)(n,double_one_,D->addr(),1,H->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,H->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  S->copy(M);
  F77NAME(daxpy)(n,double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SymmetricTridiagonalMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator-(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<double,double> *S=
    OPERATOR_NEW SymmetricTridiagonalMatrix<double,double>(n);
  S->copy(*this);
  F77NAME(daxpy)(n,double_mone_,T.D->addr(),1,S->D->addr(),1);
  F77NAME(daxpy)(n-1,double_mone_,T.L->addr(),1,S->L->addr(),1);
  return S;
}

template<> TridiagonalMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator-(
const TridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,double> *S=
    OPERATOR_NEW TridiagonalMatrix<double,double>(n);
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(0,1),1);
  F77NAME(daxpy)(n,double_mone_,T.addr(0,0),1,S->addr(0,0),1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(1,0),1,S->addr(1,0),1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(0,1),1,S->addr(0,1),1);
  return S;
}

template<> TridiagonalMatrix<double,double>* operator-(
const TridiagonalMatrix<double,double> &T,
const SymmetricTridiagonalMatrix<double,double> &St) {
  int n=St.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,double> *S=
    OPERATOR_NEW TridiagonalMatrix<double,double>(n);
  F77NAME(dcopy)(n,T.addr(0,0),1,S->addr(0,0),1);
  F77NAME(dcopy)(n-1,T.addr(1,0),1,S->addr(1,0),1);
  F77NAME(dcopy)(n-1,T.addr(0,1),1,S->addr(0,1),1);
  F77NAME(daxpy)(n,double_mone_,St.diagonalAddr(0),1,S->addr(0,0),1);
  F77NAME(daxpy)(n-1,double_mone_,St.lowerDiagonalAddr(0),1,
    S->addr(1,0),1);
  F77NAME(daxpy)(n-1,double_mone_,St.lowerDiagonalAddr(0),1,
    S->addr(0,1),1);
  return S;
}

template<> SymmetricMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator-(
const SymmetricMatrix<double,double> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricMatrix<double,double> *T=
    OPERATOR_NEW SymmetricMatrix<double,double>(n,double_zero_);
  F77NAME(dcopy)(n,D->addr(),1,T->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,T->addr(1,0),n+1);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,S.addr(j,j),1,T->addr(j,j),1);
  }
  return T;
}

template<> SymmetricMatrix<double,double>* operator-(
const SymmetricMatrix<double,double> &S,
const SymmetricTridiagonalMatrix<double,double> &T) {
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,double> *M=
    OPERATOR_NEW SymmetricMatrix<double,double>(n);
  M->copy(S);
  F77NAME(daxpy)(n,double_mone_,T.diagonalAddr(0),1,M->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    M->addr(1,0),n+1);
  return M;
}

template<> SquareMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0); 
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(0,1),n+1);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(n-j,double_mone_,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)-=double_one_;
      if (j<n-1) {
        F77NAME(daxpy)(n-j-1,double_mone_,M.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const LowerTrapezoidalMatrix<double,double> &M,
const SymmetricTridiagonalMatrix<double,double> &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0); 
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=double_one_;
      if (j<n-1) {
        F77NAME(dcopy)(n-j-1,M.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  F77NAME(daxpy)(n,double_mone_,T.diagonalAddr(0),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0); 
  F77NAME(dcopy)(n,D->addr(),1,H->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,H->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,H->addr(0,1),n+1);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(j+1,double_mone_,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(j,double_mone_,U.addr(0,j),1,H->addr(0,j),1);
      }
      (*H)(j,j)-=double_one_;
    }
  }
  return H;
}

template<> UpperHessenbergMatrix<double,double>* operator-(
const UpperTrapezoidalMatrix<double,double> &U,
const SymmetricTridiagonalMatrix<double,double> &T) {
  int n=T.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0); 
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(j+1,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(dcopy)(j,U.addr(0,j),1,H->addr(0,j),1);
      }
      (*H)(j,j)=double_one_;
    }
  }
  F77NAME(daxpy)(n,double_mone_,T.diagonalAddr(0),1,H->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  F77NAME(dcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(dcopy)(n-1,L->addr(),1,S->addr(0,1),n+1);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n,double_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const Matrix<double,double> &M,
const SymmetricTridiagonalMatrix<double,double> &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  S->copy(M);
  F77NAME(daxpy)(n,double_mone_,T.diagonalAddr(0),1,S->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator*(
const SymmetricMatrix<double,double> &M) const {
// compute by bordering: note that
// [     tau    , lambda e_0^T ] [ sigma S^T ]
// [ e_0 lambda ,      T       ] [   S    S  ]
//   = [ tau sigma + lambda e_0^T s , tau s^T        + lambda e_0^T S ]
//   = [ e_0 lambda sigma +     T s , e_0 lambda s^T +            T S ]
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  (*S)(n-2,n-2)=(*D)[n-2]*M(n-2,n-2)+(*L)[n-2]*M(n-1,n-2);
  (*S)(n-1,n-2)=(*L)[n-2]*M(n-2,n-2)+(*D)[n-1]*M(n-1,n-2);
  (*S)(n-2,n-1)=(*D)[n-2]*M(n-1,n-2)+(*L)[n-2]*M(n-1,n-1);
  (*S)(n-1,n-1)=(*L)[n-2]*M(n-1,n-2)+(*D)[n-1]*M(n-1,n-1);
  for (int k=n-3;k>=0;k--) {
    F77NAME(dstmv)(n-k-1,double_one_,L->addr(k+1),D->addr(k+1),
      M.addr(k+1,k),1,double_zero_,S->addr(k+1,k),1); // T s
    F77NAME(daxpy)(n-k-1,(*L)[k],M.addr(k+1,k+1),1,S->addr(k,k+1),n);
      // lambda e_0^T S = lambda ( S e_0 )^T
    (*S)(k,k)=(*L)[k]*M(k+1,k); // lambda e_0^T s
    F77NAME(daxpy)(n-k,(*D)[k],M.addr(k,k),1,S->addr(k,k),n);
      // tau [ sigma , s^T ]
    F77NAME(daxpy)(n-k,(*L)[k],M.addr(k,k),1,S->addr(k+1,k),n);
      // e_0 lambda [ sigma , s^T ]
  }
  return S;
}

template<> SquareMatrix<double,double>* operator*(
const SymmetricMatrix<double,double> &M,
const SymmetricTridiagonalMatrix<double,double> &T) {
// compute by bordering: note that
// [ sigma s^T ] [     tau    , lambda e_0^T ]
// [   s    S  ] [ e_0 lambda ,           T  ]
//   = [ sigma tau + s^T e_0 lambda , sigma lambda e_0^T + s^T T ]
//   = [   s   tau +  S  e_o lambda ,   s   lambda e_0^T +  S  T ]
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  (*S)(n-2,n-2)=M(n-2,n-2)*T.diagonalValue(n-2)
               +M(n-1,n-2)*T.lowerDiagonalValue(n-2);
  (*S)(n-1,n-2)=M(n-1,n-2)*T.diagonalValue(n-2)
               +M(n-1,n-1)*T.lowerDiagonalValue(n-2);
  (*S)(n-2,n-1)=M(n-2,n-2)*T.upperDiagonalValue(n-2)
               +M(n-1,n-2)*T.diagonalValue(n-1);
  (*S)(n-1,n-1)=M(n-1,n-2)*T.upperDiagonalValue(n-2)
               +M(n-1,n-1)*T.diagonalValue(n-1);
  double *t=OPERATOR_NEW_BRACKET(double,n);
  for (int k=n-3;k>=0;k--) {
    F77NAME(dstmv)(n-k-1,double_one_,T.lowerDiagonalAddr(k+1),
      T.diagonalAddr(k+1),M.addr(k+1,k),1,double_zero_,t,1);
    F77NAME(dcopy)(n-k-1,t,1,S->addr(k,k+1),n); // s^T T = ( T s )^T
    F77NAME(daxpy)(n-k-1,T.lowerDiagonalValue(k),M.addr(k+1,k+1),1,
      S->addr(k+1,k),1); // S e_0 lambda
    (*S)(k,k)=M(k+1,k)*T.lowerDiagonalValue(k); // s^T e_0 lambda
    F77NAME(daxpy)(n-k,T.diagonalValue(k),M.addr(k,k),1,S->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(daxpy)(n-k,T.upperDiagonalValue(k),M.addr(k,k),1,
      S->addr(k,k+1),1);
      // [sigma ] bar(lambda)
      // [  s   ]
  }
  delete [] t; t=0;
  return S;
}

template<> Matrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &M) const {
// compute by columns: note that
// [     S_11        , e_j sigma e_0^T ] [ 0 ] = [ e_j sigma e_0^T ell ]
// [ e_0 sigma e_j^T ,    S_22         ] [ell] = [            S_22 ell ]
  int m=M.size(0),n=M.size(1);;
  CHECK_SAME(m,size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) (*S)(j-1,j)=upperDiagonalValue(j-1)*M(j,j);
        // sigma e_0^T ell
      if (j+1<m) { // S_22 ell
        F77NAME(dstmv)(m-j,double_one_,L->addr(j),D->addr(j),
          M.addr(j,j),1,double_zero_,S->addr(j,j),1);
      } else (*S)(j,j)=diagonalValue(j)*M(j,j);
    }
  } else {
//  note that
// [      sigma , lambda e_0^T ] [ 1 ] = [   sigma    + lambda e_0^T ell ]
// [ e_0 lambda ,            S ] [ell] = [ e_0 lambda +            S ell ]
    for (int j=0;j<n;j++) {
      if ( j>0 ) (*S)(j-1,j)=upperDiagonalValue(j-1);
      (*S)(j,j)=diagonalValue(j);
      if (j+1<m) {
        (*S)(j,j)+=upperDiagonalValue(j)*M(j+1,j);
        (*S)(j+1,j)=lowerDiagonalValue(j);
        if (j+2<m) {
          F77NAME(dstmv)(m-j-1,double_one_,L->addr(j+1),
            D->addr(j+1),M.addr(j+1,j),1,double_one_,
            S->addr(j+1,j),1);
        } else (*S)(j+1,j)+=diagonalValue(j+1)*M(j+1,j);
      }
    }
  }
  return S;
}

template<> Matrix<double,double>* operator*(
const LowerTrapezoidalMatrix<double,double> &M,
const SymmetricTridiagonalMatrix<double,double> &T) {
// compute by columns: note that
// L T e_j = L e_{j-1} T_{j-1,j} + L e_j T_{j,j} + L e_{j+1} T_{j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(m-j+1,T.upperDiagonalValue(j-1),M.addr(j-1,j-1),1,
          S->addr(j-1,j),1);
      }
      F77NAME(daxpy)(m-j,T.diagonalValue(j),M.addr(j,j),1,S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(m-j-1,T.lowerDiagonalValue(j),M.addr(j+1,j+1),1,
          S->addr(j+1,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        (*S)(j-1,j)=T.upperDiagonalValue(j-1);
        F77NAME(daxpy)(m-j,T.upperDiagonalValue(j-1),M.addr(j,j-1),1,
          S->addr(j,j),1);
      }
      (*S)(j,j)+=T.diagonalValue(j);
      if (j<m-1) {
        F77NAME(daxpy)(m-j-1,T.diagonalValue(j),M.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      if (j<n-1) {
        (*S)(j+1,j)+=T.lowerDiagonalValue(j);
        if (j<m-2) {
          F77NAME(daxpy)(m-j-2,T.lowerDiagonalValue(j),M.addr(j+2,j+1),1,
            S->addr(j+2,j),1);
        }
      }
    }
  }
  return S;
}

template<> Matrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator*(
const UpperTrapezoidalMatrix<double,double> &M) const {
// compute by columns: note that
// [       S_11      , e_k sigma e_0^T ] [ u ]= [ S_11 u ]
// [ e_0 sigma e_k^T ,    S_22         ] [   ]= [ e_0 sigma e_k^T u ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<m-1;j++) {
      F77NAME(dstmv)(j+1,double_one_,L->addr(),D->addr(),
        M.addr(0,j),1,double_zero_,S->addr(0,j),1);
      (*S)(j+1,j)=lowerDiagonalValue(j)*M(j,j);
    }
    for (int j=m-1;j<n;j++) {
      F77NAME(dstmv)(m,double_one_,L->addr(),D->addr(),
        M.addr(0,j),1,double_zero_,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(dstmv)(j,double_one_,L->addr(),D->addr(),
          M.addr(0,j),1,double_zero_,S->addr(0,j),1);
        (*S)(j-1,j)+=upperDiagonalValue(j-1);
        (*S)(j,j)=lowerDiagonalValue(j-1)*M(j-1,j);
      }
      (*S)(j,j)+=diagonalValue(j);
      if (j<m-1) (*S)(j+1,j)=lowerDiagonalValue(j);
    }
    for (int j=m;j<n;j++) {
      F77NAME(dstmv)(m,double_one_,L->addr(),D->addr(),
        M.addr(0,j),1,double_zero_,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<double,double>* operator*(
const UpperTrapezoidalMatrix<double,double> &M,
const SymmetricTridiagonalMatrix<double,double> &T) {
// compute by columns: note that
// U T e_j = U e_{j-1} T_{j-1,j} + U e_j T_{j,j} + U e_{j+1} T_{j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j,m),T.upperDiagonalValue(j-1),
          M.addr(0,j-1),1,S->addr(0,j),1);
      }
      F77NAME(daxpy)(min(j+1,m),T.diagonalValue(j),M.addr(0,j),1,
        S->addr(0,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(min(j+2,m),T.lowerDiagonalValue(j),
          M.addr(0,j+1),1,S->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(min(j-1,m),T.upperDiagonalValue(j-1),
          M.addr(0,j-1),1,S->addr(0,j),1);
        if (j<=m) (*S)(j-1,j)=T.upperDiagonalValue(j-1);
      }
      F77NAME(daxpy)(min(j,m),T.diagonalValue(j),M.addr(0,j),1,
        S->addr(0,j),1);
      if (j<m) (*S)(j,j)+=T.diagonalValue(j);
      if (j<n-1) {
        F77NAME(daxpy)(min(j+1,m),T.lowerDiagonalValue(j),
          M.addr(0,j+1),1,S->addr(0,j),1);
        if (j<m-1) (*S)(j+1,j)+=T.lowerDiagonalValue(j);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator*(
const SquareMatrix<double,double> &M) const {
  int n=M.size(0);
  CHECK_SAME(n,size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dstmv)(n,double_one_,L->addr(),D->addr(),
      M.addr(0,j),1,double_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator*(
const SquareMatrix<double,double> &M,
const SymmetricTridiagonalMatrix<double,double> &T) {
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(daxpy)(n,T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(daxpy)(n,T.diagonalValue(j),M.addr(0,j),1,S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(n,T.lowerDiagonalValue(j),M.addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dstmv)(m,double_one_,L->addr(),D->addr(),
      M.addr(0,j),1,double_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,
const SymmetricTridiagonalMatrix<double,double> &T) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,double> *S=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(daxpy)(m,T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(daxpy)(m,T.diagonalValue(j),M.addr(0,j),1,S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(m,T.lowerDiagonalValue(j),M.addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> Vector<double,double>*
SymmetricTridiagonalMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<double,double> *w=
    OPERATOR_NEW Vector<double,double>(n,double_zero_);
  F77NAME(dstmv)(n,double_one_,L->addr(),D->addr(),
    v.addr(),1,double_zero_,w->addr(),1);
  return w;
}

template<> double
SymmetricTridiagonalMatrix<double,double>::normFrobenius()
const {
  return F77NAME(dlanst)('F',dim,D->addr(),L->addr());
}

template<> double
SymmetricTridiagonalMatrix<double,double>::normInfinity()
const {
  return F77NAME(dlanst)('I',dim,D->addr(),L->addr());
}

template<> double
SymmetricTridiagonalMatrix<double,double>::normMaxEntry()
const {
  return F77NAME(dlanst)('M',dim,D->addr(),L->addr());
}

template<> double
SymmetricTridiagonalMatrix<double,double>::normOne() const {
  return F77NAME(dlanst)('O',dim,D->addr(),L->addr());
}

template<> double
SymmetricTridiagonalMatrix<double,double>::reciprocalConditionNumber(
char norm) const {
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  Vector<double,double> *UF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *UF2=
    OPERATOR_NEW Vector<double,double>(max(0,dim-2));
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*L);
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(dgttrf)(dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),ipiv,
    info);
  double rcond=undefined_;
  double *work=OPERATOR_NEW_BRACKET(double,2*dim);
  int *iwork=OPERATOR_NEW_BRACKET(int,dim);
  double anorm=F77NAME(dlangt)(norm,dim,L->addr(),D->addr(),L->addr());
  F77NAME(dgtcon)(norm,dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
    ipiv,anorm,rcond,work,iwork,info);
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  delete [] ipiv; ipiv=0;
  delete UF2; UF2=0;
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
  return rcond;
}

template<> void SymmetricTridiagonalMatrix<double,double>::stmv(
double alpha,const Vector<double,double> &x,double beta,
Vector<double,double> &b) const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  F77NAME(dstmv)(n,alpha,L->addr(),D->addr(),x.addr(),1,beta,b.addr(),1);
}

template<> void SymmetricTridiagonalMatrix<double,double>::stmm(
double alpha,const Matrix<double,double> &X,double beta,
Matrix<double,double> &B,char side) const {
  if (side=='L' || side=='l') {
    int m=size(0),n=X.size(1);
    CHECK_SAME(n,B.size(1));
    CHECK_SAME(m,X.size(0));
    CHECK_SAME(m,B.size(0));
    for (int j=0;j<n;j++) {
      F77NAME(dstmv)(m,alpha,L->addr(),D->addr(),X.addr(0,j),1,beta,
        B.addr(0,j),1);
    }
  } else {
    int m=X.size(0),n=X.size(1);
    CHECK_SAME(n,size(0));
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,B.size(1));
    if (abs(beta)==double_zero_) B=double_zero_;
    else if (beta!=double_one_) B*=beta;
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(daxpy)(m,upperDiagonalValue(j-1)*alpha,X.addr(0,j-1),1,
          B.addr(0,j),1);
      }
      F77NAME(daxpy)(m,diagonalValue(j)*alpha,X.addr(0,j),1,
        B.addr(0,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(m,lowerDiagonalValue(j)*alpha,X.addr(0,j+1),1,
          B.addr(0,j),1);
      }
    }
  }
}

template<> Vector<double,double>* eigenvalues(
const SymmetricTridiagonalMatrix<double,double> &T,
OrthogonalMatrix<double,double> *&Q) {
  int dim=T.size(0);
  if (Q!=0) CHECK_SAME(dim,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'I');
  Vector<double,double> *lambda=OPERATOR_NEW Vector<double,double>(dim);
  lambda->copy(*T.diagonal());
  Vector<double,double> *L_copy=OPERATOR_NEW Vector<double,double>(dim-1);
  L_copy->copy(*T.lowerDiagonal());
  double *work=OPERATOR_NEW_BRACKET(double,2*dim-2);
  int info;
  double *qa=( Q==0 ? 0 : Q->addr() );
  F77NAME(dsteqr)(jobz,dim,lambda->addr(),L_copy->addr(),qa,dim,
    work,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete L_copy; L_copy=0;
  return lambda;
}

template<> void SymmetricTridiagonalMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  Vector<double,double> *UF=OPERATOR_NEW Vector<double,double>(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*L);
  x.copy(b);
  int info;
  F77NAME(dgtsv)(dim,1,LF->addr(),DF->addr(),UF->addr(),x.addr(),dim,
    info);
  CHECK_SAME(info,0);
  delete DF; DF=0;
  delete LF; LF=0;
  delete UF; UF=0;
}

template<> void SymmetricTridiagonalMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  Vector<double,double> *UF=OPERATOR_NEW Vector<double,double>(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*L);
  int info;
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    X.copy(B);
    F77NAME(dgtsv)(m,n,LF->addr(),DF->addr(),UF->addr(),X.addr(),m,
      info);
    CHECK_SAME(info,0);
  } else {
    CHECK_SAME(dim,n);
    double *t=OPERATOR_NEW_BRACKET(double,n);
    Vector<double,double> *UF2=OPERATOR_NEW Vector<double,double>(dim-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
    F77NAME(dgttrf)(n,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
      ipiv,info);
    CHECK_TEST(info==0);
    for (int i=0;i<m;i++) {
      F77NAME(dcopy)(n,B.addr(i,0),m,t,1);
      F77NAME(dgttrs)('T',n,1,LF->addr(),DF->addr(),UF->addr(),
        UF2->addr(),ipiv,t,n,info);
      CHECK_TEST(info==0);
      F77NAME(dcopy)(n,t,1,X.addr(i,0),m);
    }
    delete [] ipiv; ipiv=0;
    delete UF2; UF2=0;
    delete [] t; t=0;
  }
  delete LF; LF=0;
  delete DF; DF=0;
  delete UF; UF=0;
}

template class SymmetricTridiagonalMatrix<double,double>;
template void testSymmetricTridiagonalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveMatrix<double,double>*
SymmetricPositiveTridiagonalMatrix<double,double>::operator+(
const SymmetricPositiveMatrix<double,double> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricPositiveMatrix<double,double> *T=
    OPERATOR_NEW SymmetricPositiveMatrix<double,double>(n);
  T->copy(S);
  F77NAME(daxpy)(n,double_one_,D->addr(),1,T->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,L->addr(),1,T->addr(1,0),n+1);
  return T;
}

template<> double SymmetricPositiveTridiagonalMatrix<double,double>
::reciprocalConditionNumber(char norm) const {
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  LF->copy(*L);
  DF->copy(*D);
  int info;
  F77NAME(dpttrf)(dim,DF->addr(),LF->addr(),info);
  CHECK_TEST(info==0);
  double rcond=undefined_;
  double *work=OPERATOR_NEW_BRACKET(double,2*dim);
  double anorm=normOne();
  F77NAME(dptcon)(dim,DF->addr(),LF->addr(),anorm,rcond,work,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete DF; DF=0;
  delete LF; LF=0;
  return rcond;
}

/*
template<> Vector<double,double>*
SymmetricPositiveTridiagonalMatrix<double,double>::eigenvalues(
OrthogonalMatrix<double,double> *&Q) const {
  if (Q!=0) CHECK_SAME(dim,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Vector<double,double> *lambda =OPERATOR_NEW Vector<double,double>(dim);
  lambda->copy(*D);
  Vector<double,double> *L_copy =OPERATOR_NEW Vector<double,double>(dim-1);
  L_copy->copy(*L);
  double *work=OPERATOR_NEW_BRACKET(double,2*dim-2);
  int info;
  double *qa=( Q==0 ? 0 : Q->addr() );
  F77NAME(dstev)(jobz,dim,lambda->addr(),L_copy->addr(),qa,dim,work,info);
  delete [] work; work=0;
  delete L_copy; L_copy=0;
  return lambda;
}
*/

template<> void SymmetricPositiveTridiagonalMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  LF->copy(*L);
  DF->copy(*D);
  x.copy(b);
  int info;
  F77NAME(dptsv)(dim,1,DF->addr(),LF->addr(),x.addr(),dim,info);
  CHECK_SAME(info,0);
  delete DF; DF=0;
  delete LF; LF=0;
}

template<> void SymmetricPositiveTridiagonalMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<double,double> *LF=OPERATOR_NEW Vector<double,double>(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  LF->copy(*L);
  DF->copy(*D);
  int info;
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    X.copy(B);
    F77NAME(dptsv)(m,n,DF->addr(),LF->addr(),X.addr(),m,info);
    CHECK_SAME(info,0);
  } else {
    CHECK_SAME(dim,n);
    double *t=OPERATOR_NEW_BRACKET(double,n);
    F77NAME(dpttrf)(n,DF->addr(),LF->addr(),info);
    CHECK_SAME(info,0);
    for (int i=0;i<m;i++) {
      F77NAME(dcopy)(n,B.addr(i,0),m,t,1);
      F77NAME(dpttrs)(n,1,DF->addr(),LF->addr(),t,n,info);
      CHECK_SAME(info,0);
      F77NAME(dcopy)(n,t,1,X.addr(i,0),m);
    }
    delete [] t; t=0;
  }
  delete DF; DF=0;
  delete LF; LF=0;
}

template class SymmetricPositiveTridiagonalMatrix<double,double>;
template void testSymmetricPositiveTridiagonalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SquareMatrix<double,double>*
DiagonalMatrix<double,double>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,double> *M=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  F77NAME(dcopy)(n,D->addr(),1,M->addr(0,0),n+1);
  return M;
}

template<> SymmetricTridiagonalMatrix<double,double>* operator+(
const DiagonalMatrix<double,double> &A,
const SymmetricTridiagonalMatrix<double,double> &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<double,double> *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<double,double>(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_one_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<double,double>*
DiagonalMatrix<double,double>::operator+(
const TridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,double> *S=OPERATOR_NEW
    TridiagonalMatrix<double,double>(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_one_,addr(),1,S->addr(0,0),1);
  return S;
}

template<> SymmetricMatrix<double,double>* operator+(
const DiagonalMatrix<double,double> &A,
const SymmetricMatrix<double,double> &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,double> *S=OPERATOR_NEW
    SymmetricMatrix<double,double>(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_one_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> UpperTriangularMatrix<double,double>*
DiagonalMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<double,double> *S=OPERATOR_NEW
    UpperTriangularMatrix<double,double>(n);
  S->copy(U);
  if (dynamic_cast<const UnitUpperTrapezoidalMatrix<double,double>*>(&U)
  ==0) {
    F77NAME(daxpy)(n,double_one_,addr(),1,S->addr(),n+1);
  } else {
    const double *di=D->addr();
    double *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di+double_one_;
  }
  return S;
}

template<> LowerTriangularMatrix<double,double>*
DiagonalMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<double,double> *S=OPERATOR_NEW
    LowerTriangularMatrix<double,double>(n);
  S->copy(L);
  if (dynamic_cast<const UnitLowerTrapezoidalMatrix<double,double>*>(&L)
  ==0) {
    F77NAME(daxpy)(n,double_one_,addr(),1,S->addr(),n+1);
  } else {
    const double *di=D->addr();
    double *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di+double_one_;
  }
  return S;
}

template<> SquareMatrix<double,double>*
DiagonalMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n);
  S->copy(M);
  F77NAME(daxpy)(n,double_one_,addr(),1,S->addr(),n+1);
  return S;
}

template<> SymmetricTridiagonalMatrix<double,double>* operator-(
const DiagonalMatrix<double,double> &A,
const SymmetricTridiagonalMatrix<double,double> &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<double,double> *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<double,double>(n);
  S->copy(T);
  *S*=double_mone_;
  F77NAME(daxpy)(n,double_one_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> SymmetricTridiagonalMatrix<double,double>* operator-(
const SymmetricTridiagonalMatrix<double,double> &T,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<double,double> *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<double,double>(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_mone_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<double,double>*
DiagonalMatrix<double,double>::operator-(
const TridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,double> *S=OPERATOR_NEW
    TridiagonalMatrix<double,double>(n);
  S->copy(T);
  *S*=double_mone_;
  F77NAME(daxpy)(n,double_one_,addr(),1,S->addr(0,0),1);
  return S;
}

template<> TridiagonalMatrix<double,double>* operator-(
const TridiagonalMatrix<double,double> &T,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,double> *S=OPERATOR_NEW
    TridiagonalMatrix<double,double>(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_mone_,A.addr(),1,S->addr(0,0),1);
  return S;
}

template<> SymmetricMatrix<double,double>* operator-(
const DiagonalMatrix<double,double> &A,
const SymmetricMatrix<double,double> &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,double> *S=OPERATOR_NEW
    SymmetricMatrix<double,double>(n);
  S->copy(T);
  *S*=double_mone_;
  F77NAME(daxpy)(n,double_one_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> SymmetricMatrix<double,double>* operator-(
const SymmetricMatrix<double,double> &T,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,double> *S=OPERATOR_NEW
    SymmetricMatrix<double,double>(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_mone_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> UpperTriangularMatrix<double,double>*
DiagonalMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<double,double> *S=OPERATOR_NEW
    UpperTriangularMatrix<double,double>(n);
  S->copy(U);
  *S*=double_mone_;
  if (dynamic_cast<const UnitUpperTrapezoidalMatrix<double,double>*>(&U)
  ==0) {
    F77NAME(daxpy)(n,double_one_,addr(),1,S->addr(),n+1);
  } else {
    const double *di=D->addr();
    double *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di-double_one_;
  }
  return S;
}

template<> UpperTriangularMatrix<double,double>* operator-(
const UpperTrapezoidalMatrix<double,double> &U,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<double,double> *S=OPERATOR_NEW
    UpperTriangularMatrix<double,double>(n);
  S->copy(U);
  if (dynamic_cast<const UnitUpperTrapezoidalMatrix<double,double>*>(&U)
  ==0) {
    F77NAME(daxpy)(n,double_mone_,A.addr(),1,S->addr(),n+1);
  } else {
    const double *di=A.addr();
    double *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=double_one_-*di;
  }
  return S;
}

template<> LowerTriangularMatrix<double,double>*
DiagonalMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<double,double> *S=OPERATOR_NEW
    LowerTriangularMatrix<double,double>(n);
  S->copy(L);
  *S*=double_mone_;
  if (dynamic_cast<const UnitLowerTrapezoidalMatrix<double,double>*>(&L)
  ==0) {
    F77NAME(daxpy)(n,double_one_,addr(),1,S->addr(),n+1);
  } else {
    const double *di=D->addr();
    double *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di-double_one_;
  }
  return S;
}

template<> LowerTriangularMatrix<double,double>* operator-(
const LowerTrapezoidalMatrix<double,double> &L,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<double,double> *S=OPERATOR_NEW
    LowerTriangularMatrix<double,double>(n);
  S->copy(L);
  if (dynamic_cast<const UnitLowerTrapezoidalMatrix<double,double>*>(&L)
  ==0) {
    F77NAME(daxpy)(n,double_mone_,A.addr(),1,S->addr(),n+1);
  } else {
    const double *di=A.addr();
    double *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=double_one_-*di;
  }
  return S;
}

template<> SquareMatrix<double,double>*
DiagonalMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n);
  S->copy(M);
  *S*=double_mone_;
  F77NAME(daxpy)(n,double_one_,addr(),1,S->addr(),n+1);
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const Matrix<double,double> &M,const DiagonalMatrix<double,double> &A)
{
  int n=A.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n);
  S->copy(M);
  F77NAME(daxpy)(n,double_mone_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> TridiagonalMatrix<double,double>* 
DiagonalMatrix<double,double>::operator*(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,double> *S=
    OPERATOR_NEW TridiagonalMatrix<double,double>(n);
  for (int i=0;i<n;i++) {
    double di=(*D)[i];
    if (i>0) (*S)(i,i-1)=di*T.lowerDiagonalValue(i-1);
    (*S)(i,i)=di*T.diagonalValue(i);
    if (i<n-1) (*S)(i,i+1)=di*T.lowerDiagonalValue(i);
  }
  return S;
}

template<> TridiagonalMatrix<double,double>* operator*(
const SymmetricTridiagonalMatrix<double,double> &T,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,double> *S=
    OPERATOR_NEW TridiagonalMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    double dj=A[j];
    if (j>0) (*S)(j-1,j)=T.lowerDiagonalValue(j-1)*dj;
    (*S)(j,j)=T.diagonalValue(j)*dj;
    if (j<n-1) (*S)(j+1,j)=T.lowerDiagonalValue(j)*dj;
  }
  return S;
}

template<> SquareMatrix<double,double>*
DiagonalMatrix<double,double>::operator*(
const SymmetricMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(dim);
  for (int i=0;i<n;i++) {
    double di=D->operator[](i);
    if (i>0) {
      const double *Tij=T.addr(i,0);
      double *Sij=S->addr(i,0);
      for (int j=0;j<i;j++,Tij+=n,Sij+=n) *Sij=di*(*Tij);
    }
    const double *Tji=T.addr(i,i);
    double *Sij=S->addr(i,i);
    for (int j=i;j<n;j++,Tji++,Sij+=n) *Sij=di*(*Tji);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator*(
const SymmetricMatrix<double,double> &T,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    double dj=A[j];
    if (j>0) {
      const double *Tji=T.addr(j,0);
      double *Sij=S->addr(0,j);
      for (int i=0;i<j;i++,Tji+=n,Sij++) *Sij=(*Tji)*dj;
    }
    const double *Tij=T.addr(j,j);
    double *Sij=S->addr(j,j);
    for (int i=j;i<n;i++,Tij++,Sij++) *Sij=(*Tij)*dj;
  }
  return S;
}

template class DiagonalMatrix<double,double>;
template void testDiagonalMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const double
  UpperHessenbergMatrix<double,double>::outofbounds_ = double_zero_;
template<> const double
  UpperHessenbergMatrix<double,double>::undefined_ =
    numeric_limits<double>::infinity();
template<> double UpperHessenbergMatrix<double,double>::safety_ =
  double_zero_;

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,double> *M=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(n,j+2),addr(0,j),1,M->addr(0,j),1);
  }
  return M;
}

template<> UpperHessenbergMatrix<double,double>&
UpperHessenbergMatrix<double,double>::operator+=(
const UpperHessenbergMatrix<double,double> &H) {
  int n=size(0);
  CHECK_SAME(n,H.size(0));
  double *colj=addr();
  const double *Hcolj=H.addr();
  for (int j=0;j<n;j++,colj+=n,Hcolj+=n) {
    F77NAME(daxpy)(min(n,j+2),double_one_,Hcolj,1,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<double,double>&
UpperHessenbergMatrix<double,double>::operator-=(
const UpperHessenbergMatrix<double,double> &H) {
  int n=size(0);
  CHECK_SAME(n,H.size(0));
  double *colj=addr();
  const double *Hcolj=H.addr();
  for (int j=0;j<n;j++,colj+=n,Hcolj+=n) {
    F77NAME(daxpy)(min(n,j+2),double_mone_,Hcolj,1,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<double,double>&
UpperHessenbergMatrix<double,double>::operator*=(double scalar) {
  int n=size(0);
  double *colj=addr();
  for (int j=0;j<n;j++,colj+=n) {
    F77NAME(dscal)(min(n,j+2),scalar,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<double,double>&
UpperHessenbergMatrix<double,double>::operator/=(double scalar) {
  int n=size(0);
  double *colj=addr();
  double s=double_one_/scalar;
  for (int j=0;j<n;j++,colj+=n) {
    F77NAME(dscal)(min(n,j+2),s,colj,1);
  }
  return *this;
}

template<> void UpperHessenbergMatrix<double,double>::copy(
const Matrix<double,double> &M) {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(j+2,n),M.addr(0,j),1,addr(0,j),1);
  }
}

/*
template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::transpose() const {
  int n=size(0);
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(j,0),n);
  }
  return S;
}
*/

template<> void UpperHessenbergMatrix<double,double>::copyFrom(
int m,const SquareMatrix<double,double> &S) {
  m=min(m,size(0));
  m=min(m,S.size(0));
  for (int j=0;j<m;j++) {
    F77NAME(dcopy)(min(j+2,m),S.addr(0,j),1,addr(0,j),1);
  }
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator+(
const DiagonalMatrix<double,double> &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(*this);
  F77NAME(daxpy)(n,double_one_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator+(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(*this);
  F77NAME(daxpy)(n,double_one_,T.diagonalAddr(0),1,H->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,T.lowerDiagonalAddr(0),1,
    H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator+(
const TridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(*this);
  F77NAME(daxpy)(n,double_one_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator+(
const SymmetricMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) F77NAME(dcopy)(n-j-1,M.addr(j+1,j),1,S->addr(j,j+1),n);
    F77NAME(daxpy)(min(j+2,n),double_one_,addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      F77NAME(daxpy)(n-j,double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)+=double_one_;
      if (j<n-1) {
        F77NAME(daxpy)(n-j-1,double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(*this);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(j+1,double_one_,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)+=double_one_;
      if (j>0) {
        F77NAME(daxpy)(j,double_one_,U.addr(0,j),1,H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(j+2,n),double_one_,addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator-(
const DiagonalMatrix<double,double> &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(*this);
  F77NAME(daxpy)(n,double_mone_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,double>* operator-(
const DiagonalMatrix<double,double> &D,
const UpperHessenbergMatrix<double,double> &U) {
  int n=U.size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(U);
  *H*=double_mone_;
  F77NAME(daxpy)(n,double_one_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator-(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(*this);
  F77NAME(daxpy)(n,double_mone_,T.diagonalAddr(0),1,H->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,double>* operator-(
const SymmetricTridiagonalMatrix<double,double> &T,
const UpperHessenbergMatrix<double,double> &U) {
  int n=U.size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(U);
  *H*=double_mone_;
  F77NAME(daxpy)(n,double_one_,T.diagonalAddr(0),1,H->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,T.lowerDiagonalAddr(0),1,
    H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator-(
const TridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(*this);
  F77NAME(daxpy)(n,double_mone_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_mone_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,double>* operator-(
const TridiagonalMatrix<double,double> &T,
const UpperHessenbergMatrix<double,double> &U) {
  int n=U.size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(U);
  *H*=double_mone_;
  F77NAME(daxpy)(n,double_one_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(daxpy)(n-1,double_one_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator-(
const SymmetricMatrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(n-j,double_mone_,M.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      F77NAME(daxpy)(n-j-1,double_mone_,M.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const SymmetricMatrix<double,double> &M,
const UpperHessenbergMatrix<double,double> &H) {
  int n=H.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      F77NAME(dcopy)(n-j-1,M.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(j+2,n),double_mone_,H.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      F77NAME(daxpy)(n-j,double_mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)-=double_one_;
      if (j<n-1) {
        F77NAME(daxpy)(n-j-1,double_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const LowerTrapezoidalMatrix<double,double> &L,
const UpperHessenbergMatrix<double,double> &H) {
  int n=H.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      F77NAME(daxpy)(min(j+2,n),double_mone_,H.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=double_one_;
      if (j<n-1) {
        F77NAME(dcopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      F77NAME(daxpy)(min(j+2,n),double_mone_,H.addr(0,j),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(*this);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(j+1,double_mone_,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)-=double_one_;
      if (j>0) {
        F77NAME(daxpy)(j,double_mone_,U.addr(0,j),1,H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> UpperHessenbergMatrix<double,double>* operator-(
const UpperTrapezoidalMatrix<double,double> &U,
const UpperHessenbergMatrix<double,double> &M) {
  int n=M.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n,double_zero_);
  H->copy(M);
  *H*=double_mone_;
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(daxpy)(j+1,double_one_,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)+=double_one_;
      if (j>0) {
        F77NAME(daxpy)(j,double_one_,U.addr(0,j),1,H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
    F77NAME(daxpy)(n,double_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const Matrix<double,double> &M,
const UpperHessenbergMatrix<double,double> &H) {
  int n=H.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n,M.addr(0,j),1,S->addr(0,j),1);
    F77NAME(daxpy)(min(j+2,n),double_mone_,H.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const UpperHessenbergMatrix<double,double> &H2) const {
  int n=size(0);
  CHECK_SAME(n,H2.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<=min(j+1,n-1);k++) {
      F77NAME(daxpy)(min(k+2,n),H2(k,j),addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const DiagonalMatrix<double,double> &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,double> *S=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n);
  S->copy(*this);
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(min(j+2,n),D[j],S->addr(0,j),1);
  }
  return S;
}

template<> UpperHessenbergMatrix<double,double>* operator*(
const DiagonalMatrix<double,double> &D,
const UpperHessenbergMatrix<double,double> &H) {
  int n=H.size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,double> *S=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(n);
  S->copy(H);
  for (int i=0;i<n;i++) {
    int j=max(0,i-1);
    F77NAME(dscal)(n-j,D[i],S->addr(i,j),n);
  }
  return S;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const SymmetricTridiagonalMatrix<double,double> &T) const {
//compute by columns: note that
// H T e_j = H e_{j-1} T_{j-1,j} + H e_j T_{j,j} + H e_{j+1} T_{j+1,j}
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(daxpy)(min(j+1,n),T.upperDiagonalValue(j-1),addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(daxpy)(min(j+2,n),T.diagonalValue(j),addr(0,j),1,
      S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(min(j+3,n),T.lowerDiagonalValue(j),addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator*(
const SymmetricTridiagonalMatrix<double,double> &T,
const UpperHessenbergMatrix<double,double> &H) {
// compute by rows: note that
// e_i^T T H
//   = T_{i,i-1} e_{i-1}^T H + T_{i,i} e_i^T H + T_{i,i+1} e_{i+1}^T H
  int n=H.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int i=0;i<n;i++) {
    if (i>0) {
      int j=max(0,i-2);
      F77NAME(daxpy)(n-j,T.lowerDiagonalValue(i-1),H.addr(i-1,j),n,
        S->addr(i,j),n);
    }
    int j=max(0,i-1);
    F77NAME(daxpy)(n-j,T.diagonalValue(i),H.addr(i,j),n,S->addr(i,j),n);
    if (i<n-1) {
      F77NAME(daxpy)(n-i,T.upperDiagonalValue(i),H.addr(i+1,i),n,
        S->addr(i,i),n);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const TridiagonalMatrix<double,double> &T) const {
//compute by columns: note that
// H T e_j = H e_{j-1} T_{j-1,j} + H e_j T_{j,j} + H e_{j+1} T_{j+1,j}
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(j-1,0);k<=min(j+1,n-1);k++) {
      F77NAME(daxpy)(min(k+2,n),T(k,j),addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator*(
const TridiagonalMatrix<double,double> &T,
const UpperHessenbergMatrix<double,double> &H) {
// compute by rows: note that
// e_i^T T H
//   = T_{i,i-1} e_{i-1}^T H + T_{i,i} e_i^T H + T_{i,i+1} e_{i+1}^T H
  int n=H.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,double> *S=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int i=0;i<n;i++) {
    for (int k=max(i-1,0);k<=min(i+1,n-1);k++) {
      int j=max(0,k-1);
      F77NAME(daxpy)(min(n-k+1,n),T(i,k),H.addr(k,j),n,S->addr(i,j),n);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const SymmetricMatrix<double,double> &S) const {
// compute by bordering: note that
// [     eta_11 h^T ] [ sigma s^T ]
// [ e_0 eta_21  H  ] [   s    S  ]
//   = [     eta_11 sigma + h^T s ,     eta_11 s^T + h^T S ]
//   = [ e_0 eta_21 sigma +  H  s , e_0 eta_21 s^T +  H  S ]
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *M=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  (*M)(n-1,n-1)=(*this)(n-1,n-1)*S(n-1,n-1);
  for (int k=n-2;k>=0;k--) {
    for (int j=k+1;j<n;j++) { // H s
      F77NAME(daxpy)(min(n-k-1,j-k+1),S(j,k),addr(k+1,j),1,
        M->addr(k+1,k),1);
    }
    F77NAME(dsymv)('L',n-k-1,double_one_,S.addr(k+1,k+1),n,
      addr(k,k+1),n,double_zero_,M->addr(k,k+1),n); // h^T S = ( S h )^T
    (*M)(k,k)=F77NAME(ddot)(n-k-1,addr(k,k+1),n,S.addr(k+1,k),1); // h^T s
    F77NAME(dger)(2,n-k,double_one_,addr(k,k),1,S.addr(k,k),1,
      M->addr(k,k),n);
  }
  return M;
}

template<> SquareMatrix<double,double>* operator*(
const SymmetricMatrix<double,double> &S,
const UpperHessenbergMatrix<double,double> &H) {
// compute by bordering: note that
// [ sigma s^T ] [     eta_11 h^T ]
// [   s    S  ] [ e_0 eta_21  H  ]
//   = [ sigma eta_11 + s^T e_0 eta_21 , sigma h^T + s^T H ]
//   = [   s   eta_11 +  S  e_0 eta_21 ,   s   h^T +  S  H ]
  int n=H.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *M=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  (*M)(n-1,n-1)=S(n-1,n-1)*H(n-1,n-1);
  for (int k=n-2;k>=0;k--) { // s^T H
    F77NAME(daxpy)(n-k-1,H(k+1,k),S.addr(k+1,k+1),1,M->addr(k+1,k),1);
      // S e_0 eta_21
    for (int j=k+1;j<n;j++) { // s^T H
      (*M)(k,j)=
        F77NAME(ddot)(min(n-k-1,j-k+1),S.addr(k+1,k),1,H.addr(k+1,j),1);
    }
    (*M)(k,k)=S(k+1,k)*H(k+1,k); // s^T e_0 eta_21
    F77NAME(dger)(n-k,n-k,double_one_,S.addr(k,k),1,H.addr(k,k),n,
      M->addr(k,k),n);
  }
  return M;
}

template<> Matrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &L) const {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<m;k++) {
        F77NAME(daxpy)(min(k+2,m),L(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(min(j+2,m),addr(0,j),1,M->addr(0,j),1);
      for (int k=j+1;k<m;k++) {
        F77NAME(daxpy)(min(k+2,m),L(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<double,double>* operator*(
const LowerTrapezoidalMatrix<double,double> &L,
const UpperHessenbergMatrix<double,double> &H) {
// compute by columns: note that
//   [ L_11      ] [ h ] = [ L_11 h ]
//   [ L_21 L_22 ] [   ] = [ L_21 h ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,H.size(0));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0 ? 'N' : 'U');
  for (int j=0;j<n;j++) {
    int k=min(j+2,n);
    F77NAME(dcopy)(k,H.addr(0,j),1,M->addr(0,j),1);
    F77NAME(dtrmv)('L','N',diag,k,L.addr(),m,M->addr(0,j),1); // L_11 h
    if (k<m) {
      F77NAME(dgemv)('N',m-k,k,double_one_,L.addr(k,0),m,H.addr(0,j),1,
        double_zero_,M->addr(k,j),1); // L_21 h
    }
  }
  return M;
}

template<> Matrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const UpperTrapezoidalMatrix<double,double> &U) const {
// compute by columns: note that
//   H [ U_1 , U_2 ] = [ H U_1 , H U_2 ]
// and that
//   [      H_11     H_12 ] [ u ] = [     H_11      u ]
//   [ e_0 eta e_k^T H_22 ] [   ] = [ e_0 eta e_k^T u ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0);
  if (U_non_unit) { // H U_1
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,m-1);k++) {
        F77NAME(daxpy)(min(k+2,m),U(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,m);k++) {
        F77NAME(daxpy)(min(k+2,m),U(k,j),addr(0,k),1,M->addr(0,j),1);
      }
      if (j<m) {
        F77NAME(daxpy)(min(j+2,m),double_one_,addr(0,j),1,M->addr(0,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<double,double>* operator*(
const UpperTrapezoidalMatrix<double,double> &U,
const UpperHessenbergMatrix<double,double> &H) {
// compute by columns: note that
//   [ U_1 , U_2 ] [     H_11          , H_12 ]
//                 [ e_0 eta e_{m-1}^T , H_22 ]
//   = [ U_1 H_11 + U_2 e_0 eta e_{m-1}^T , U_1 H_12 + U_2 H_22 ]
// and that
//   [ U_11 U_12 ] [ h ] = [ U_11 h ]
//   [      U_22 ] [   ] = [        ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,H.size(0));
  char diag=(
    dynamic_cast<const UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0
    ? 'N' : 'U');
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<m;j++) { // U_1 H_11
    int k=min(j+2,m);
    F77NAME(dcopy)(k,H.addr(0,j),1,M->addr(0,j),1);
    F77NAME(dtrmv)('U','N',diag,k,U.addr(),m,M->addr(0,j),1); // U_11 h
  }
  if (m<n) {
    F77NAME(daxpy)(m,H(m,m-1),U.addr(0,m),1,M->addr(0,m-1),1);
      // U_2 e_0 eta e_{m-1}^T
    F77NAME(dlacpy)('A',m,n-m,H.addr(0,m),m,M->addr(0,m),m);
    F77NAME(dtrmm)('L','U','N',diag,m,n-m,double_one_,U.addr(),m,
      M->addr(0,m),m); // U_1 H_12
    F77NAME(dgemm)('N','N',m,n-m,n-m,double_one_,U.addr(0,m),m,
      H.addr(m,m),n,double_one_,M->addr(0,m),m); // U_2 H_22
  }
  return M;
}

template<> SquareMatrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const SquareMatrix<double,double> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *M=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<n;k++) {
      F77NAME(daxpy)(min(k+2,n),S(k,j),addr(0,k),1,M->addr(0,j),1);
    }
  }
  return M;
}

template<> SquareMatrix<double,double>* operator*(
const SquareMatrix<double,double> &S,
const UpperHessenbergMatrix<double,double> &H) {
  int n=H.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *M=
    OPERATOR_NEW SquareMatrix<double,double>(n);
  for (int j=0;j<n;j++) {
    F77NAME(dgemv)('N',n,min(j+2,n),double_one_,S.addr(),n,H.addr(0,j),1,
      double_zero_,M->addr(0,j),1);
  }
  return M;
}

template<> Matrix<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,double> *R=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(j-1,0);k<m;k++) {
      F77NAME(daxpy)(min(k+2,m),M(k,j),addr(0,k),1,R->addr(0,j),1);
    }
  }
  return R;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,
const UpperHessenbergMatrix<double,double> &H) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,H.size(0));
  Matrix<double,double> *R=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dgemv)('N',m,min(j+2,n),double_one_,M.addr(),m,H.addr(0,j),1,
      double_zero_,R->addr(0,j),1);
  }
  return R;
}

template<> Vector<double,double>*
UpperHessenbergMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<double,double> *r=
    OPERATOR_NEW Vector<double,double>(n,double_zero_);
  for (int k=0;k<n;k++) {
    F77NAME(daxpy)(min(k+2,n),v[k],addr(0,k),1,r->addr(),1);
  }
  return r;
}

template<> void UpperHessenbergMatrix<double,double>::uhmv(double alpha,
const Vector<double,double> &x,double beta,Vector<double,double> &b,
char trans) const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  if (abs(beta)==double_zero_) b=double_zero_;
  else b*=beta;
  if (trans=='N' || trans=='n') { // b=H*x*alpha+b*beta
    double *xj=x.addr();
    for (int j=0;j<n;j++,xj++) {
      F77NAME(daxpy)(min(n,j+2),(*xj)*alpha,addr(0,j),1,b.addr(),1); 
    }
  } else { // b=H^T*x*alpha+b*beta
    double *bj=b.addr();
    for (int j=0;j<n;j++,bj++) {
      *bj+=F77NAME(ddot)(min(n,j+2),addr(0,j),1,x.addr(),1)*alpha;
    }
  }
}

template<> void UpperHessenbergMatrix<double,double>::uhmm(double alpha,
const Matrix<double,double> &X,double beta,Matrix<double,double> &B,
char side,char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (abs(beta)==double_zero_) B=double_zero_;
  else B*=beta;
  if (trans=='N' || trans=='n') {
    if (side=='L' || side=='l') { // B=H*X*alpha+B*beta
      CHECK_SAME(m,size(0));
      for (int j=0;j<n;j++) {
        for (int k=0;k<m;k++) {
          F77NAME(daxpy)(min(m,k+2),X(k,j)*alpha,addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else { // B=X*H*alpha+B*beta
      CHECK_SAME(n,size(0));
      for (int j=0;j<n;j++) {
        for (int k=0;k<=min(n-1,j+1);k++) {
          F77NAME(daxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  } else {
    if (side=='L' || side=='l') { // B=H^T*X*alpha+B*beta
      CHECK_SAME(m,size(0));
      for (int j=0;j<n;j++) {
        for (int i=0;i<m;i++) {
          B(i,j)+=F77NAME(ddot)(min(m,i+2),addr(0,i),1,X.addr(0,j),1)
                *alpha;
        }
      }
    } else { // B=X*H^T*alpha+B*beta
      CHECK_SAME(n,size(0));
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<n;k++) {
          F77NAME(daxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

template<> double UpperHessenbergMatrix<double,double>::normFrobenius()
const {
  int n=size(0);
  double *work=0;
  return F77NAME(dlanhs)('F',n,addr(),n,work);
}

template<> double UpperHessenbergMatrix<double,double>::normInfinity()
const {
  int n=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,n);
  double val=F77NAME(dlanhs)('I',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> double UpperHessenbergMatrix<double,double>::normMaxEntry()
const {
  int n=size(0);
  double *work=0;
  return F77NAME(dlanhs)('M',n,addr(),n,work);
}

template<> double UpperHessenbergMatrix<double,double>::normOne()
const {
  int n=size(0);
  double *work=0;
  return F77NAME(dlanhs)('O',n,addr(),n,work);
}

template<> Vector<double,complex<double> >*
UpperHessenbergMatrix<double,double>::eigenvalues(
SquareMatrix<double,complex<double> > *&V,
SquareMatrix<double,complex<double> > *&U) const {
  int n=size(0);
  if (V!=0) CHECK_SAME(n,V->size(0));
  if (U!=0) CHECK_SAME(n,U->size(0));
  UpperHessenbergMatrix<double,double> *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(*this);
  char job=(V==0 && U==0 ? 'E' : 'S');
  char compz=(V==0 && U==0 ? 'N' : 'I');
  double *wr=OPERATOR_NEW_BRACKET(double,n);
  double *wi=OPERATOR_NEW_BRACKET(double,n);
  OrthogonalMatrix<double,double> *Z=(V==0 && U==0 ? 0 :
    OPERATOR_NEW OrthogonalMatrix<double,double>(n,n));
  double *za=(Z==0 ? 0 : Z->addr());
  double w=double_undefined_;
  int lwork=-1;
  int info;
  F77NAME(dhseqr)(job,compz,n,1,n,H->addr(),n,wr,wi,za,n,&w,lwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dhseqr)(job,compz,n,1,n,H->addr(),n,wr,wi,za,n,work,lwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;

  Vector<double,complex<double> > *lambda =
    OPERATOR_NEW Vector<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    complex<double> &ev=(*lambda)[j];
    ev.real()=wr[j];
    ev.imag()=wi[j];
  }

  if (Z!=0) {
    char side=(V==0 ? 'R' : (U==0 ? 'L' : 'B') );
    bool *select=0;
    double *vla=0;
    SquareMatrix<double,double> *Vl=0;
    if (V!=0) {
      Vl=OPERATOR_NEW SquareMatrix<double,double>(n);
      Vl->copy(*Z);
      vla=Vl->addr();
    }
    double *vra=0;
    SquareMatrix<double,double> *Vr=0;
    if (U!=0) {
      Vr=OPERATOR_NEW SquareMatrix<double,double>(n);
      Vr->copy(*Z);
      vra=Vr->addr();
    }
    work=OPERATOR_NEW_BRACKET(double,3*n);
    int mout=-1;
    F77NAME(dtrevc)(side,'B',select,n,H->addr(),n,vla,n,vra,n,n,mout,work,
      info);
    CHECK_TEST(info==0);
    for (int j=0;j<n;) {
      int jn=j*n;
      if (wi[j]>zero_) {
        int jp1n=jn+n;
        if (V!=0) {
          for (int i=0;i<n;i++) {
            complex<double> &z=V->operator()(i,j);
            z.real()=vla[i+jn];
            z.imag()=vla[i+jp1n];
            complex<double> &zz=V->operator()(i,j+1);
            zz.real()=vla[i+jn];
            zz.imag()=-vla[i+jp1n];
          }
        }
        if (U!=0) {
          for (int i=0;i<n;i++) {
            complex<double> &z=U->operator()(i,j);
            z.real()=vra[i+jn];
            z.imag()=vra[i+jp1n];
            complex<double> &zz=U->operator()(i,j+1);
            zz.real()=vra[i+jn];
            zz.imag()=-vra[i+jp1n];
          }
        }
        j+=2;
      } else {
        if (V!=0) {
          for (int i=0;i<n;i++) {
            complex<double> &z=V->operator()(i,j);
            z.real()=vla[i+jn];
            z.imag()=zero_;
          }
        }
        if (U!=0) {
          for (int i=0;i<n;i++) {
            complex<double> &z=U->operator()(i,j);
            z.real()=vra[i+jn];
            z.imag()=zero_;
          }
        }
        j++;
      }
    }
    if (V!=0) {
      for (int j=0;j<n;j++) {
        double val=F77NAME(dznrm2)(n,V->addr(0,j),1);
        F77NAME(zdscal)(n,1./val,V->addr(0,j),1);
      }
    }
    if (U!=0) {
      for (int j=0;j<n;j++) {
        double val=F77NAME(dznrm2)(n,U->addr(0,j),1);
        F77NAME(zdscal)(n,1./val,U->addr(0,j),1);
      }
    }
    delete [] work; work=0;
    if (Vl!=0) delete Vl; Vl=0;
    if (Vr!=0) delete Vr; Vr=0;
  }
  delete [] wr; wr=0;
  delete [] wi; wi=0;
  if (Z!=0) delete Z; Z=0;
  delete H; H=0;
  return lambda;
}  

template<> int UpperHessenbergMatrix<double,double>::factor(int *ipiv) {
  int n=size(0);
  for (int k=0;k<n;k++) ipiv[k]=k;
  int info=0;
  for (int k=0;k<n-1;k++) {
    if (abs((*this)(k,k))<abs((*this)(k+1,k)) ) {
      ipiv[k]=k+1;
      F77NAME(dswap)(n-k,addr(k,k),n,addr(k+1,k),n);
    }
    if (abs((*this)(k,k))>zero_) {
      (*this)(k+1,k)/=(*this)(k,k);
      F77NAME(daxpy)(n-k-1,-(*this)(k+1,k),addr(k,k+1),n,
        addr(k+1,k+1),n);
    } else {
      info=k+1;
      break;
    }
  }
  return info;
}

template<> void UpperHessenbergMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  int n=size(0);
  CHECK_SAME(n,b.size());
  CHECK_SAME(n,x.size());
  x.copy(b);
  UpperHessenbergMatrix<double,double> *HF=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int info=HF->factor(ipiv);
  CHECK_TEST(info==0);
  if (trans!='N' && trans!='n') {
    x[0]/=(*HF)(0,0);
    for (int j=1;j<n;j++) {
      x[j]=(x[j]-F77NAME(ddot)(j,HF->addr(0,j),1,x.addr(),1))
          /(*HF)(j,j);
    }
    for (int i=n-2;i>=0;i--) {
      int ip=ipiv[i];
      double temp=x[i]-(*HF)(i+1,i)*x[i+1];
      x[i]=x[ip];
      x[ip]=temp;
    }
  } else {
    for (int i=0;i<n-1;i++) {
      int ip=ipiv[i];
      double temp=x[2*i-ip+1]-(*HF)(i+1,i)*x[ip];
      x[i]=x[ip];
      x[i+1]=temp;
    }
    for (int j=n-1;j>0;j--) {
      x[j]/=(*HF)(j,j);
      F77NAME(daxpy)(j,-x[j],HF->addr(0,j),1,x.addr(),1);
    }
    x[0]/=(*HF)(0,0);
  }
  delete ipiv; ipiv=0;
  delete HF; HF=0;
}

template<> void UpperHessenbergMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  int n=size(0);
  UpperHessenbergMatrix<double,double> *HF=
    OPERATOR_NEW UpperHessenbergMatrix<double,double>(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int info=HF->factor(ipiv);
  CHECK_TEST(info==0);
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(n,X.size(0));
    CHECK_SAME(nrhs,X.size(1));
    X.copy(B);
    if (trans!='N' && trans!='n') { // dgttsf
      for (int k=0;k<nrhs;k++) {
        X(0,k)/=(*HF)(0,0);
        for (int j=1;j<n;j++) {
          X(j,k)=(X(j,k)-F77NAME(ddot)(j,HF->addr(0,j),1,X.addr(0,k),1))
              /(*HF)(j,j);
        }
        for (int i=n-2;i>=0;i--) {
          int ip=ipiv[i];
          double temp=X(i,k)-(*HF)(i+1,i)*X(i+1,k);
          X(i,k)=X(ip,k);
          X(ip,k)=temp;
        }
      }
    } else {
      for (int k=0;k<nrhs;k++) {
        for (int i=0;i<n-1;i++) {
          int ip=ipiv[i];
          double temp=X(2*i-ip+1,k)-(*HF)(i+1,i)*X(ip,k);
          X(i,k)=X(ip,k);
          X(i+1,k)=temp;
        }
        for (int j=n-1;j>0;j--) {
          X(j,k)/=(*HF)(j,j);
          F77NAME(daxpy)(j,-X(j,k),HF->addr(0,j),1,X.addr(0,k),1);
        }
        X(0,k)/=(*HF)(0,0);
      }
    }
  } else {
    int nrhs=B.size(0);
    CHECK_SAME(n,B.size(1));
    CHECK_SAME(n,X.size(1));
    CHECK_SAME(nrhs,X.size(0));
    X.copy(B);
    if (trans!='N' && trans!='n') {
      for (int k=0;k<nrhs;k++) {
        for (int i=0;i<n-1;i++) {
          int ip=ipiv[i];
          double temp=X(k,2*i-ip+1)-(*HF)(i+1,i)*X(k,ip);
          X(k,i)=X(k,ip);
          X(k,i+1)=temp;
        }
        for (int j=n-1;j>0;j--) {
          X(k,j)/=(*HF)(j,j);
          F77NAME(daxpy)(j,-X(k,j),HF->addr(0,j),1,X.addr(k,0),nrhs);
        }
        X(k,0)/=(*HF)(0,0);
      }
    } else {
      for (int k=0;k<nrhs;k++) {
        X(k,0)/=(*HF)(0,0);
        for (int j=1;j<n;j++) {
          X(k,j)=(X(k,j)
            -F77NAME(ddot)(j,HF->addr(0,j),1,X.addr(k,0),nrhs))
            /(*HF)(j,j);
        }
        for (int i=n-2;i>=0;i--) {
          int ip=ipiv[i];
          double temp=X(k,i)-(*HF)(i+1,i)*X(k,i+1);
          X(k,i)=X(k,ip);
          X(k,ip)=temp;
        }
      }
    }
  }
  delete ipiv; ipiv=0;
  delete HF; HF=0;
}

template class UpperHessenbergMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const double BandMatrix<double,double>::outofbounds_ =
  double_zero_;
template<> double BandMatrix<double,double>::safety_ = double_zero_;

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::makeMatrix() const {
  SquareMatrix<double,double> *M=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      M->addr(ibeg,j),1);
  }
  return M;
}

template<> BandMatrix<double,double>&
BandMatrix<double,double>::operator+=(const BandMatrix<double,double> &B)
{
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  CHECK_SAME(nsup,B.nsup);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,double_one_,
      B.addr(ibeg,j),1,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<double,double>&
BandMatrix<double,double>::operator-=(const BandMatrix<double,double> &B)
{
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  CHECK_SAME(nsup,B.nsup);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<double,double>&
BandMatrix<double,double>::operator*=(double d) {
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dscal)(min(dim-1,j+nsub)-ibeg+1,d,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<double,double>&
BandMatrix<double,double>::operator/=(double d) {
  CHECK_TEST(abs(d)>double_zero_);
  double dinv=double_one_/d;
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dscal)(min(dim-1,j+nsub)-ibeg+1,dinv,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator+(const BandMatrix<double,double> &B)
const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,B.nsub),max(nsup,B.nsup),
    double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    ibeg=max(0,j-B.nsup);
    F77NAME(daxpy)(min(dim-1,j+B.nsub)-ibeg+1,double_one_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator+(
const UpperHessenbergMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    F77NAME(daxpy)(min(dim,j+2),double_one_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator+(
const DiagonalMatrix<double,double> &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<double,double> *S=
    OPERATOR_NEW BandMatrix<double,double>(*this);
  F77NAME(daxpy)(dim,double_one_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator+(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,1),max(nsup,1),double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(daxpy)(dim-1,double_one_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  F77NAME(daxpy)(dim,double_one_,T.diagonalAddr(0),1,
    S->addr(0,0),S->nt);
  F77NAME(daxpy)(dim-1,double_one_,T.lowerDiagonalAddr(0),1,
    S->addr(0,1),S->nt);
  return S;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator+(
const TridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,1),max(nsup,1),double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(daxpy)(dim,double_one_,T.diagonalAddr(),1,
    S->addr(0,0),S->nt);
  F77NAME(daxpy)(dim-1,double_one_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S->nt);
  F77NAME(daxpy)(dim-1,double_one_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S->nt);
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator+(
const SymmetricMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<dim;j++) {
    F77NAME(daxpy)(dim-j,double_one_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(daxpy)(dim-j-1,double_one_,H.addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(daxpy)(j+1,double_one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      if (j>0) {
        F77NAME(daxpy)(j,double_one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(daxpy)(dim-j,double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      (*S)(j,j)+=double_one_;
      if (j+1<dim) {
        F77NAME(daxpy)(dim-j-1,double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  S->copy(M);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,double_one_,
      addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator-(const BandMatrix<double,double> &B)
const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,B.nsub),max(nsup,B.nsup),
    double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    ibeg=max(0,j-B.nsup);
    F77NAME(daxpy)(min(B.dim-1,j+B.nsub)-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator-(
const UpperHessenbergMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    F77NAME(daxpy)(min(dim,j+2),double_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const UpperHessenbergMatrix<double,double> &H,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(n,j+2),H.addr(0,j),1,S->addr(0,j),1);
    int ibeg=max(0,j-B.supDiags());
    F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator-(
const DiagonalMatrix<double,double> &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<double,double> *S=
    OPERATOR_NEW BandMatrix<double,double>(*this);
  F77NAME(daxpy)(dim,double_mone_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> BandMatrix<double,double>* operator-(
const DiagonalMatrix<double,double> &D,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(n,B.subDiags(),B.supDiags(),double_zero_);
  F77NAME(dcopy)(n,D.addr(),1,S->addr(0,0),S->bands());
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator-(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,1),max(nsup,1),double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(daxpy)(dim-1,double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  F77NAME(daxpy)(dim,double_mone_,T.diagonalAddr(0),1,
    S->addr(0,0),S->nt);
  F77NAME(daxpy)(dim-1,double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(0,1),S->nt);
  return S;
}

template<> BandMatrix<double,double>* operator-(
const SymmetricTridiagonalMatrix<double,double> &T,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(n,max(B.subDiags(),1),max(B.supDiags(),1),
    double_zero_);
  int nb=S->bands();
  F77NAME(dcopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(1,0),nb);
  F77NAME(dcopy)(n,T.diagonalAddr(0),1,S->addr(0,0),nb);
  F77NAME(dcopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(0,1),nb);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator-(
const TridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,1),max(nsup,1),double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(daxpy)(dim,double_mone_,T.diagonalAddr(),1,
    S->addr(0,0),nt);
  F77NAME(daxpy)(dim-1,double_mone_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),nt);
  F77NAME(daxpy)(dim-1,double_mone_,T.upperDiagonalAddr(),1,
    S->addr(0,1),nt);
  return S;
}

template<> BandMatrix<double,double>* operator-(
const TridiagonalMatrix<double,double> &T,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(n,max(B.subDiags(),1),max(B.supDiags(),1),
    double_zero_);
  F77NAME(dcopy)(n,T.diagonalAddr(),1,S->addr(0,0),S->bands());
  F77NAME(dcopy)(n-1,T.lowerDiagonalAddr(),1,S->addr(1,0),S->bands());
  F77NAME(dcopy)(n-1,T.upperDiagonalAddr(),1,S->addr(0,1),S->bands());
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator-(
const SymmetricMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<dim;j++) {
    F77NAME(daxpy)(dim-j,double_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(daxpy)(dim-j-1,double_mone_,H.addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const SymmetricMatrix<double,double> &H,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,H.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      F77NAME(dcopy)(n-j-1,H.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(daxpy)(j+1,double_mone_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      if (j>0) {
        F77NAME(daxpy)(j,double_mone_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)-=double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const UpperTrapezoidalMatrix<double,double> &U,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(j+1,U.addr(0,j),1,S->addr(0,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-B.supDiags());
      if (j>0) F77NAME(dcopy)(j,U.addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)=double_one_;
      F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(daxpy)(dim-j,double_mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      (*S)(j,j)-=double_one_;
      if (j+1<dim) {
        F77NAME(daxpy)(dim-j-1,double_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const LowerTrapezoidalMatrix<double,double> &L,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=double_one_;
      if (j+1<n) F77NAME(dcopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  S->operator-=(M);
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const Matrix<double,double> &M,const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  S->copy(M);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(daxpy)(min(n-1,j+B.subDiags())-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator*(const BandMatrix<double,double> &B)
const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(dim,nsub+B.nsub,nsup+B.nsup,double_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=max(0,j-B.nsup);k<=min(dim-1,j+B.nsub);k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,B(k,j),
        addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator*(
const UpperHessenbergMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=0;k<=min(dim-1,j+1);k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,H(k,j),
        addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<double,double>* operator*(
const UpperHessenbergMatrix<double,double> &H,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(daxpy)(min(n,k+2),B(k,j),H.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator*(
const DiagonalMatrix<double,double> &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<double,double> *P=
    OPERATOR_NEW BandMatrix<double,double>(*this);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(dscal)(min(dim-1,j+nsub)-ibeg+1,D[j],
      P->addr(ibeg,j),1);
  }
  return P;
}

template<> BandMatrix<double,double>* operator*(
const DiagonalMatrix<double,double> &D,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  BandMatrix<double,double> *P=
    OPERATOR_NEW BandMatrix<double,double>(n,B.subDiags(),B.supDiags());
  P->copy(B);
  int stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    int jbeg=max(0,i-B.subDiags());
    F77NAME(dscal)(min(n-1,i+B.supDiags())-jbeg+1,D[i],
      P->addr(i,jbeg),stride);
  }
  return P;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator*(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(dim,nsub+1,nsup+1,double_zero_);
  for (int j=0;j<dim;j++) {
    if (j>0) {
      int ibeg=max(0,j-1-nsup);
      F77NAME(daxpy)(min(dim-1,nsub+j-1)-ibeg+1,
        T.upperDiagonalValue(j-1),addr(ibeg,j-1),1,P->addr(ibeg,j),1);
    }
    int ibeg=max(0,j-nsup);
    F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,T.diagonalValue(j),
      addr(ibeg,j),1,P->addr(ibeg,j),1);
    if (j<dim-1) {
      int ibeg=max(0,j+1-nsup);
      F77NAME(daxpy)(min(dim-1,j+1+nsub)-ibeg+1,
        T.lowerDiagonalValue(j),addr(ibeg,j+1),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> BandMatrix<double,double>* operator*(
const SymmetricTridiagonalMatrix<double,double> &T,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(n,B.subDiags()+1,B.supDiags()+1,
    double_zero_);
  int B_stride=B.bands()-1;
  int P_stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    if (i>0) {
      int jbeg=max(0,i-1-B.subDiags());
      F77NAME(daxpy)(min(n-1,i-1+B.supDiags())-jbeg+1,
        T.lowerDiagonalValue(i-1),B.addr(i-1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
    int jbeg=max(0,i-B.subDiags());
    F77NAME(daxpy)(min(n-1,i+B.supDiags())-jbeg+1,T.diagonalValue(i),
      B.addr(i,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    if (i<n-1) {
      int jbeg=max(0,i+1-B.subDiags());
      F77NAME(daxpy)(min(n-1,i+1+B.supDiags())-jbeg+1,
        T.upperDiagonalValue(i),B.addr(i+1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
  }
  return P;
}

template<> BandMatrix<double,double>*
BandMatrix<double,double>::operator*(
const TridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(dim,nsub+1,nsup+1,double_zero_);
  for (int j=0;j<dim;j++) {
    if (j>0) {
      int ibeg=max(0,j-1-nsup);
      F77NAME(daxpy)(min(dim-1,nsub+j-1)-ibeg+1,T(j-1,j),
        addr(ibeg,j-1),1,P->addr(ibeg,j),1);
    }
    int ibeg=max(0,j-nsup);
    F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,T(j,j),addr(ibeg,j),1,
      P->addr(ibeg,j),1);
    if (j<dim-1) {
      int ibeg=max(0,j+1-nsup);
      F77NAME(daxpy)(min(dim-1,j+1+nsub)-ibeg+1,T(j+1,j),
        addr(ibeg,j+1),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> BandMatrix<double,double>* operator*(
const TridiagonalMatrix<double,double> &T,
const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(n,B.subDiags()+1,B.supDiags()+1,
    double_zero_);
  int B_stride=B.bands()-1;
  int P_stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    if (i>0) {
      int jbeg=max(0,i-1-B.subDiags());
      F77NAME(daxpy)(min(n-1,i-1+B.supDiags())-jbeg+1,T(i,i-1),
        B.addr(i-1,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    }
    int jbeg=max(0,i-B.subDiags());
    F77NAME(daxpy)(min(n-1,i+B.supDiags())-jbeg+1,T(i,i),
      B.addr(i,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    if (i<n-1) {
      int jbeg=max(0,i+1-B.subDiags());
      F77NAME(daxpy)(min(n-1,i+1+B.supDiags())-jbeg+1,
        T(i,i+1),B.addr(i+1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
  }
  return P;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator*(
const SymmetricMatrix<double,double> &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<double,double> *P=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=0;k<dim;k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,S(k,j),addr(ibeg,k),1,
        P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<double,double>* operator*(
const SymmetricMatrix<double,double> &S,
const BandMatrix<double,double> &B) {
//compute by bordering: note that
//  [ sigma s^T ] [ beta c^T ]
//  [   s    S  ] [   b   B  ]
//  = [ sigma beta + s^T b , sigma c^T + s^T B ]
//  = [   s   beta +  S  b ,   s  c^T +   S  B ]
  int n=B.size(0),nsub=B.subDiags(),nsup=B.supDiags();
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *P=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  for (int k=n-1;k>=0;k--) {
    int iend=min(n-1,k+nsub);
    if (k<n-1) {
      // s^T B = ( B^T s )^T
      F77NAME(dgbmv)('T',n-k-1,n-k-1,nsub,nsup,
        double_one_,B.addr(k+1-nsup,k+1),B.bands(),S.addr(k+1,k),1,
        double_zero_,P->addr(k,k+1),n);
      // S b: note that
      // [ S_11 , S_21^T ] [ b ] = [ S_11 b ]
      // [ S_21 ,  S_22  ] [ 0 ] = [ S_21 b ]
      if (nsub>0) { // S_11 b
        F77NAME(dsymv)('L',min(n-k-1,nsub),double_one_,S.addr(k+1,k+1),n,
          B.addr(k+1,k),1,double_zero_,P->addr(k+1,k),1);
        if (nsub<n-k-1) { // S_21 b
          F77NAME(dgemv)('N',n-k-1-nsub,nsub,double_one_,
            S.addr(k+1+nsub,k+1),n,B.addr(k+1,k),1,double_zero_,
            P->addr(k+1+nsub,k),1);
        }
      }
    }
    if (iend>k) { // s^T b
      (*P)(k,k)=F77NAME(ddot)(iend-k,S.addr(k+1,k),1,B.addr(k+1,k),1);
    }
    int jend=min(n-1,k+nsup);
    F77NAME(dger)(n-k,jend-k+1,double_one_,S.addr(k,k),1,
      B.addr(k,k),B.bands()-1,P->addr(k,k),n);
  }
  return P;
}

template<> Matrix<double,double>*
BandMatrix<double,double>::operator*(
const UpperTrapezoidalMatrix<double,double> &U) const {
  CHECK_SAME(dim,U.size(0));
  int n=U.size(1);
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(dim,U.size(1),double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,dim-1);k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,dim);k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
      if (j<dim) {
        int ibeg=max(0,j-nsup);
        F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,double_one_,
          addr(ibeg,j),1,M->addr(ibeg,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<double,double>* operator*(
const UpperTrapezoidalMatrix<double,double> &U,
const BandMatrix<double,double> &B) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(daxpy)(min(k+1,m),B(k,j),U.addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(daxpy)(min(k,m),B(k,j),U.addr(0,k),1,M->addr(0,j),1);
        if (k<m) (*M)(k,j)+=B(k,j);
      }
    }
  }
  return M;
}

template<> Matrix<double,double>*
BandMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &L) const {
  CHECK_SAME(dim,L.size(0));
  int n=L.size(1);
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(dim,L.size(1),double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<dim;k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,double_one_,addr(ibeg,j),1,
        M->addr(ibeg,j),1);
      for (int k=j+1;k<dim;k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<double,double>* operator*(
const LowerTrapezoidalMatrix<double,double> &L,
const BandMatrix<double,double> &B) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,double> *M=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(daxpy)(m-k,B(k,j),L.addr(k,k),1,M->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        (*M)(k,j)+=B(k,j);
        if (k+1<m) {
          F77NAME(daxpy)(m-k-1,B(k,j),L.addr(k+1,k),1,M->addr(k+1,j),1);
        }
      }
    }
  }
  return M;
}

template<> SquareMatrix<double,double>*
BandMatrix<double,double>::operator*(
const SquareMatrix<double,double> &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dgbmv)('N',dim,dim,nsub,nsup,double_one_,addr(),nt,
      S.addr(0,j),1,double_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> SquareMatrix<double,double>* operator*(
const SquareMatrix<double,double> &S,const BandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(daxpy)(n,B(k,j),S.addr(0,j),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Matrix<double,double>* BandMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  CHECK_SAME(dim,M.size(0));
  int n=M.size(1);
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(dim,n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dgbmv)('N',dim,dim,nsub,nsup,double_one_,addr(),nt,
      M.addr(0,j),1,double_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,const BandMatrix<double,double> &B) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(daxpy)(m,B(k,j),M.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Vector<double,double>* BandMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  CHECK_SAME(dim,v.size());
  Vector<double,double> *w=
    OPERATOR_NEW Vector<double,double>(dim,double_zero_);
  F77NAME(dgbmv)('N',dim,dim,nsub,nsup,double_one_,addr(),nt,
    v.addr(),1,double_zero_,w->addr(),1);
  return w;
}

template<> void BandMatrix<double,double>::gbmv(double alpha,
const Vector<double,double> &x,double beta,Vector<double,double> &b,
char trans) const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  F77NAME(dgbmv)(trans,dim,dim,nsub,nsup,alpha,addr(),nt,x.addr(),1,
    beta,b.addr(),1);
}

template<> void BandMatrix<double,double>::gbmm(double alpha,
const Matrix<double,double> &X,double beta,Matrix<double,double> &B,
char side,char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      F77NAME(dgbmv)(trans,dim,dim,nsub,nsup,alpha,addr(),nt,
        X.addr(0,j),1,beta,B.addr(0,j),1);
    }
  } else {
    CHECK_SAME(n,dim);
    if (abs(beta)==double_zero_) B=double_zero_;
    else B*=beta;
    if (trans=='N' || trans=='n') {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsup);k<=min(dim-1,j+nsub);k++) {
          F77NAME(daxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsub);k<=min(dim-1,j+nsup);k++) {
          F77NAME(daxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

/*
template<> BandMatrix<double,double>*
BandMatrix<double,double>::transpose() const {
  BandMatrix<double,double> *X=
    OPERATOR_NEW BandMatrix<double,double>(dim,nsup,nsub);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsub);
    F77NAME(dcopy)(min(dim-1,j+nsup)-i+1,addr(j,i),nsub+nsup,
      X->addr(i,j),1);
  }
  return X;
}
*/

template<> double BandMatrix<double,double>::equilibrate(
Vector<double,double> &r,Vector<double,double> &c,double &rowcnd,
double &colcnd) const {
  CHECK_SAME(dim,r.size());
  CHECK_SAME(dim,c.size());
  double amax;
  int info;
  F77NAME(dgbequ)(dim,dim,nsub,nsup,addr(),nt,r.addr(),c.addr(),
    rowcnd,colcnd,amax,info);
  CHECK_TEST(info==0);
  return amax;
}

template<> double BandMatrix<double,double>::normFrobenius() const {
  double *work=0;
  return F77NAME(dlangb)('F',dim,nsub,nsup,addr(),nt,work);
}

template<> double BandMatrix<double,double>::normInfinity() const {
  double *work=OPERATOR_NEW_BRACKET(double,dim);
  double val=F77NAME(dlangb)('I',dim,nsub,nsup,addr(),nt,work);
  delete work;
  return val;
}

template<> double BandMatrix<double,double>::normMaxEntry() const {
  double *work=0;
  return F77NAME(dlangb)('M',dim,nsub,nsup,addr(),nt,work);
}

template<> double BandMatrix<double,double>::normOne() const {
  double *work=0;
  return F77NAME(dlangb)('O',dim,nsub,nsup,addr(),nt,work);
}

template<> double BandMatrix<double,double>::reciprocalConditionNumber(
char norm) const {
  BandMatrix<double,double> *BF=
    OPERATOR_NEW BandMatrix<double,double>(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(dgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);

  double anorm=(norm=='I' || norm=='i' ? normInfinity() : normOne());
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,3*dim);
  int *iwork=OPERATOR_NEW_BRACKET(int,dim);
  F77NAME(dgbcon)(norm,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,anorm,
    rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete [] iwork; iwork=0;
  delete [] ipiv; ipiv=0;
  delete BF; BF=0;
  return rcond;
}

template<> void BandMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char trans)
const {
  CHECK_SAME(dim,b.size())
  CHECK_SAME(dim,x.size())
  BandMatrix<double,double> *BF=
    OPERATOR_NEW BandMatrix<double,double>(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(dgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);

  x.copy(b);
  F77NAME(dgbtrs)(trans,dim,nsub,nsup,1,BF->addr(),BF->nt,ipiv,
    x.addr(),dim,info);
  CHECK_TEST(info==0);
  delete [] ipiv;
  delete BF; BF=0;
}

template<> void BandMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,
char trans) const {
  BandMatrix<double,double> *BF=
    OPERATOR_NEW BandMatrix<double,double>(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(dcopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(dgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(dim,B.size(0))
    CHECK_SAME(dim,X.size(0))
    CHECK_SAME(nrhs,X.size(1))
    X.copy(B);
    F77NAME(dgbtrs)(trans,dim,nsub,nsup,nrhs,BF->addr(),BF->nt,ipiv,
      X.addr(),dim,info);
    CHECK_TEST(info==0);
  } else { // dgbtrs
    int nrhs=B.size(0);
    CHECK_SAME(dim,B.size(1))
    CHECK_SAME(dim,X.size(1))
    CHECK_SAME(nrhs,X.size(0))
    X.copy(B);
    int kd=nsub+nsup+1;
    if (trans!='N' && trans!='n') {
      if (nsub>0) { // solve L Y^T = B^T: dgbtrs
        for (int j=0;j<dim-1;j++) {
          int lm=min(nsub,dim-j-1);
          int jp=ipiv[j]-1;
          if (jp!=j) F77NAME(dswap)(nrhs,X.addr(0,jp),1,X.addr(0,j),1);
          F77NAME(dger)(nrhs,lm,double_mone_,BF->addr(j+1,j),1,
            X.addr(0,j),1,X.addr(0,j+1),1);
        }
      }
      for (int i=0;i<nrhs;i++) { // solve R X^T = Y^T: dtbsv loop 20
        for (int j=dim-1;j>=0;j--) {
          if (abs(X(i,j))>double_zero_) {
            X(i,j)/=(*BF)(j,j);
            int ii=max(0,j-nsub-nsup);
            F77NAME(daxpy)(j-ii,-X(i,j),BF->addr(ii,j),1,
              X.addr(i,ii),nrhs);
          }
        }
      }
    } else {
      for (int i=0;i<nrhs;i++) { // solve R^T Y^T=B^T: dtbsv loop 100
        for (int j=0;j<dim;j++) {
          int ii=max(0,j-nsub-nsup);
          X(i,j)=(X(i,j)
            -F77NAME(ddot)(j-i,BF->addr(i,j),1,X.addr(i,ii),1))
            /(*BF)(j,j);
        }
      }
      if (nsub>0) { // solve L^T X^T = Y^T: dgbtrs loop 40
        for (int j=dim-2;j>=0;j--) {
          int lm=min(nsub,dim-1-j);
          F77NAME(dgemv)('T',nrhs,lm,double_mone_,X.addr(0,j+1),nrhs,
            BF->addr(j+1,j),1,double_one_,X.addr(0,j),nrhs);
          int jp=ipiv[j]-1;
          if (jp!=j) {
            F77NAME(dswap)(nrhs,X.addr(0,jp),1,X.addr(0,j),1);
          }
        }
      }
    }
  }

  delete [] ipiv;
  delete BF; BF=0;
}

template class BandMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const double SymmetricBandMatrix<double,double>::outofbounds_
  = double_zero_;
template<> double SymmetricBandMatrix<double,double>::safety_ =
  double_zero_;

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::makeMatrix() const {
  SquareMatrix<double,double> *M=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(dim-j,nsub+1),addr(j,j),1,M->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(dcopy)(min(dim-j-1,nsub),addr(j+1,j),1,M->addr(j,j+1),dim);
    }
  }
  return M;
}

template<> void SymmetricBandMatrix<double,double>::fillWith(double d) {
  double *colj=addr();
  for (int j=0;j<dim;j++,colj+=nt) {
    for (int i=0;i<min(dim-j,nt);i++) colj[i]=d;
  }
}

template<> double SymmetricBandMatrix<double,double>::operator()(int i,
int j) const {
  if (i-j>=0 && i-j<=nsub) return *AB->addr(i-j,j);
  if (j-i>0 && j-i<=nsub) return *AB->addr(j-i,i);
  return outofbounds_;
}

template<> SymmetricBandMatrix<double,double>::SymmetricBandMatrix(
const SymmetricTridiagonalMatrix<double,double> &T) : nsub(1), nt(2) {
  dim=T.size(0);
  AB=OPERATOR_NEW Matrix<double,double>(2,dim);
  F77NAME(dcopy)(dim,T.diagonalAddr(0),1,addr(0,0),2);
  F77NAME(dcopy)(dim-1,T.lowerDiagonalAddr(0),1,addr(1,0),2);
  (*AB)(1,dim-1)=numeric_limits<double>::infinity();
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::makeBandMatrix() const {
  BandMatrix<double,double> *B=
    OPERATOR_NEW BandMatrix<double,double>(dim,nsub,nsub);
  int stride=B->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,B->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,B->addr(j,j+1),stride);
    }
  }
  return B;
}

template<> SymmetricMatrix<double,double>*
SymmetricBandMatrix<double,double>::makeSymmetricMatrix() const {
  SymmetricMatrix<double,double> *S=
    OPERATOR_NEW SymmetricMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<double,double>&
SymmetricBandMatrix<double,double>::operator+=(
const SymmetricBandMatrix<double,double> &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  for (int j=0;j<dim;j++) {
    F77NAME(daxpy)(min(nt,dim-j),double_one_,B.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<double,double>&
SymmetricBandMatrix<double,double>::operator-=(
const SymmetricBandMatrix<double,double> &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  for (int j=0;j<dim;j++) {
    F77NAME(daxpy)(min(nt,dim-j),double_mone_,B.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<double,double>&
SymmetricBandMatrix<double,double>::operator*=(double d) {
  for (int j=0;j<dim;j++) {
    F77NAME(dscal)(min(nt,dim-j),d,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<double,double>&
SymmetricBandMatrix<double,double>::operator/=(double d) {
  CHECK_TEST(abs(d)>double_zero_);
  double dinv=double_one_/d;
  for (int j=0;j<dim;j++) {
    F77NAME(dscal)(min(nt,dim-j),dinv,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const SymmetricBandMatrix<double,double> &B) const {
  CHECK_SAME(dim,B.dim);
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(dim,max(nsub,B.nsub),double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(daxpy)(min(B.nt,dim-j),double_one_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const BandMatrix<double,double> &B) const {
  CHECK_SAME(dim,B.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,B.subDiags()),
    max(nsub,B.supDiags()),double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
        S->addr(j,j+1),stride);
    }
    int ibeg=max(0,j-B.supDiags());
    F77NAME(daxpy)(min(dim-1,j+B.subDiags())-ibeg+1,double_one_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const UpperHessenbergMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
    F77NAME(daxpy)(min(dim,j+2),double_one_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const DiagonalMatrix<double,double> &D) const {
  CHECK_SAME(dim,D.size(0));
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(*this);
  F77NAME(daxpy)(dim,double_one_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> SymmetricBandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(dim,max(nsub,1),double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  F77NAME(daxpy)(dim-1,double_one_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  F77NAME(daxpy)(dim,double_one_,T.diagonalAddr(0),1,S->addr(0,0),S->nt);
  return S;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const TridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,1),max(nsub,1),double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,S->addr(j,j+1),stride);
    }
  }
  F77NAME(daxpy)(dim,double_one_,T.diagonalAddr(),1,
    S->addr(0,0),S->bands());
  F77NAME(daxpy)(dim-1,double_one_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S->bands());
  F77NAME(daxpy)(dim-1,double_one_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S->bands());
  return S;
}

template<> SymmetricMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const SymmetricMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SymmetricMatrix<double,double> *S=makeSymmetricMatrix();
  for (int j=0;j<dim;j++) {
    F77NAME(daxpy)(dim-j,double_one_,H.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const UpperTrapezoidalMatrix<double,double> &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      F77NAME(daxpy)(j+1,double_one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      if (j>0) {
        F77NAME(daxpy)(j,double_one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)+=double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const LowerTrapezoidalMatrix<double,double> &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      F77NAME(daxpy)(dim-j,double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      (*S)(j,j)+=double_one_;
      if (j+1<dim) {
        F77NAME(daxpy)(dim-j-1,double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator+(
const Matrix<double,double> &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  S->copy(M);
  for (int j=0;j<dim;j++) {
    F77NAME(daxpy)(min(nt,dim-j),double_one_,addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(daxpy)(min(nt,dim-j)-1,double_one_,addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
  }
  return S;
}

template<> SymmetricBandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const SymmetricBandMatrix<double,double> &B) const {
  CHECK_SAME(dim,B.dim);
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(dim,max(nsub,B.nsub),double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(daxpy)(min(B.nt,dim-j),double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const BandMatrix<double,double> &B) const {
  CHECK_SAME(dim,B.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,B.subDiags()),
      max(nsub,B.supDiags()),double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
        S->addr(j,j+1),stride);
    }
    int ibeg=max(0,j-B.supDiags());
    F77NAME(daxpy)(min(dim-1,j+B.subDiags())-ibeg+1,double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>* operator-(
const BandMatrix<double,double> &B,
const SymmetricBandMatrix<double,double> &H) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(n,max(H.subDiags(),B.subDiags()),
    max(H.subDiags(),B.supDiags()),double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(dcopy)(min(n-1,j+B.subDiags())-ibeg+1,B.addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(H.bands(),n-j),double_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(min(H.bands(),n-j)-1,double_mone_,
        H.addr(j+1,j),1,S->addr(j,j+1),stride);
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const UpperHessenbergMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,S->addr(j,j+1),dim);
    }
    F77NAME(daxpy)(min(dim,j+2),double_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const UpperHessenbergMatrix<double,double> &H,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(n,j+2),H.addr(0,j),1,S->addr(0,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(min(B.bands(),n-j)-1,double_mone_,
        B.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  return S;
}

template<> SymmetricBandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const DiagonalMatrix<double,double> &D) const {
  CHECK_SAME(dim,D.size(0));
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(*this);
  F77NAME(daxpy)(dim,double_mone_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> SymmetricBandMatrix<double,double>* operator-(
const DiagonalMatrix<double,double> &D,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(n,B.subDiags(),double_zero_);
  int nt=B.bands();
  F77NAME(dcopy)(n,D.addr(),1,S->addr(0,0),nt);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(nt,n-j),double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(dim,max(nsub,1),double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  F77NAME(daxpy)(dim-1,double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  F77NAME(daxpy)(dim,double_mone_,T.diagonalAddr(0),1,S->addr(0,0),S->nt);
  return S;
}

template<> SymmetricBandMatrix<double,double>* operator-(
const SymmetricTridiagonalMatrix<double,double> &T,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(n,max(B.subDiags(),1),
    double_zero_);
  F77NAME(dcopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(1,0),S->bands());
  F77NAME(dcopy)(n,T.diagonalAddr(0),1,S->addr(0,0),S->bands());
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const TridiagonalMatrix<double,double> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(dim,max(nsub,1),max(nsub,1),double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,S->addr(j,j+1),stride);
    }
  }
  int S_nt=S->bands();
  F77NAME(daxpy)(dim,double_mone_,T.diagonalAddr(),1,
    S->addr(0,0),S_nt);
  F77NAME(daxpy)(dim-1,double_mone_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S_nt);
  F77NAME(daxpy)(dim-1,double_mone_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S_nt);
  return S;
}

template<> BandMatrix<double,double>* operator-(
const TridiagonalMatrix<double,double> &T,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,double> *S=OPERATOR_NEW
    BandMatrix<double,double>(n,max(B.subDiags(),1),max(B.subDiags(),1),
    double_zero_);
  int S_nt=S->bands();
  F77NAME(dcopy)(n,T.diagonalAddr(),1,S->addr(0,0),S_nt);
  F77NAME(dcopy)(n-1,T.lowerDiagonalAddr(),1,S->addr(1,0),S_nt);
  F77NAME(dcopy)(n-1,T.upperDiagonalAddr(),1,S->addr(0,1),S_nt);
  int stride=S->bands()-1;
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(min(B.bands(),n-j)-1,double_mone_,B.addr(j+1,j),1,
        S->addr(j,j+1),stride);
    }
  }
  return S;
}

template<> SymmetricMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const SymmetricMatrix<double,double> &H) const {
  CHECK_SAME(dim,H.size(0));
  SymmetricMatrix<double,double> *S=makeSymmetricMatrix();
  for (int j=0;j<dim;j++) {
    F77NAME(daxpy)(dim-j,double_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
//  if (j+1<dim) {
//    F77NAME(daxpy)(dim-j-1,double_mone_,H.addr(j+1,j),1,
//      S->addr(j,j+1),dim);
//  }
  }
  return S;
}

template<> SymmetricMatrix<double,double>* operator-(
const SymmetricMatrix<double,double> &H,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SymmetricMatrix<double,double> *S=OPERATOR_NEW
    SymmetricMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(n-j,H.addr(j,j),1,S->addr(j,j),1);
    F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const UpperTrapezoidalMatrix<double,double> &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      F77NAME(daxpy)(j+1,double_mone_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      if (j>0) {
        F77NAME(daxpy)(j,double_mone_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)-=double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const UpperTrapezoidalMatrix<double,double> &U,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(j+1,U.addr(0,j),1,S->addr(0,j),1);
      F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
        S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(min(B.bands(),n-j)-1,double_mone_,
          B.addr(j+1,j),1,S->addr(j,j+1),n);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) F77NAME(dcopy)(j,U.addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)=double_one_;
      F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
        S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(min(B.bands(),n-j)-1,double_mone_,
          B.addr(j+1,j),1,S->addr(j,j+1),n);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const LowerTrapezoidalMatrix<double,double> &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      F77NAME(daxpy)(dim-j,double_mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      (*S)(j,j)-=double_one_;
      if (j+1<dim) {
        F77NAME(daxpy)(dim-j-1,double_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const LowerTrapezoidalMatrix<double,double> &L,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(dcopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
        S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(min(B.bands(),n-j)-1,double_mone_,
          B.addr(j+1,j),1,S->addr(j,j+1),n);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=double_one_;
      if (j+1<n) {
        F77NAME(dcopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
        S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(daxpy)(min(B.bands(),n-j)-1,double_mone_,
          B.addr(j+1,j),1,S->addr(j,j+1),n);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator-(
const Matrix<double,double> &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(dcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(dcopy)(min(nt,dim-j)-1,addr(j+1,j),1,S->addr(j,j+1),dim);
    }
    F77NAME(daxpy)(dim,double_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,double>* operator-(
const Matrix<double,double> &M,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,double> *S=OPERATOR_NEW
    SquareMatrix<double,double>(n,double_zero_);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(daxpy)(min(B.bands(),n-j),double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(daxpy)(min(B.bands(),n-j)-1,double_mone_,
        B.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  return S;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const SymmetricBandMatrix<double,double> &B) const {
// compute by bordering: note that
// [ sigma s^T ] [ tau t^T ] = [ sigma tau + s^T t , sigma t^T + s^T T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^T +  S  T ]
  CHECK_SAME(dim,B.dim);
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(dim,nsub+B.nsub,nsub+B.nsub,double_zero_);
  int nb=min(nsub,B.nsub);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*B(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(dim-1,k+B.nsub);kk++) { // S t
        int ibeg=max(k+1,kk-B.nsub);
//      for (int i=ibeg;i<kk;i++) {
//        (*P)(i,k)+=(*this)(kk,i)*B(kk,k);
//      }
        if (kk>ibeg) {
          F77NAME(daxpy)(kk-ibeg,B(kk,k),addr(kk,ibeg),nsub,
            P->addr(ibeg,k),1);
        }
//      for (int i=max(kk,ibeg);i<=min(kk+nsub,dim-1);i++) {
//        (*P)(i,k)+=(*this)(i,kk)*B(kk,k);
//      }
        ibeg=max(kk,ibeg);
        F77NAME(daxpy)(min(kk+nsub,dim-1)-ibeg+1,B(kk,k),
          addr(ibeg,kk),1,P->addr(ibeg,k),1);
      }
      for (int j=k+1;j<=min(k+nsub+B.nsub,dim-1);j++) { // s^T T
//      for (int kk=max(k+1,j-B.nsub);kk<min(j,k+nsub);kk++) {
//        (*P)(k,j)+=(*this)(kk,k)*B(j,kk);
//      }
        int kbeg=max(k+1,j-B.nsub);
        int kend=min(j-1,k+nsub);
        if (kbeg<=kend) {
          (*P)(k,j)=F77NAME(ddot)(kend-kbeg+1,addr(kbeg,k),1,
            B.addr(j,kbeg),B.nt-1);
        }
//      for (int kk=max(kbeg,j);
//      kk<=min(dim-1,min(k+nsub,j+B.nsub));kk++) {
//        (*P)(k,j)+=(*this)(kk,k)*B(kk,j);
//      }
        kbeg=max(kbeg,j);
        kend=min(dim-1,min(k+nsub,j+B.nsub));
        if (kend>=kbeg) {
          (*P)(k,j)+=
            F77NAME(ddot)(kend-kbeg+1,addr(kbeg,k),1,B.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^T t
        F77NAME(ddot)(min(nb,dim-k-1),addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(dger)(min(dim-k,nt),min(dim-k,B.nt),double_one_,
      addr(k,k),1,B.addr(k,k),1,P->addr(k,k),P->bands()-1);
  }
  return P;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const BandMatrix<double,double> &B) const {
// compute by bordering: note that
// [ sigma s^T ] [ beta c^T ] = [ sigma beta + s^T b , sigma c^T + s^T B ]
// [   s    S  ] [  b    B  ] = [   s   beta +  S  b ,   s   c^T +  S  B ]
  CHECK_SAME(dim,B.size(0));
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(dim,nsub+B.subDiags(),nsub+B.supDiags(),
    double_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*B(dim-1,dim-1);
  int nb=min(nsub,B.subDiags());
  for (int k=dim-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(dim-1,k+B.subDiags());kk++) { // S b
        int ibeg=max(k+1,kk-B.supDiags());
        if (kk>ibeg) {
          F77NAME(daxpy)(kk-ibeg,B(kk,k),addr(kk,ibeg),nsub,
            P->addr(ibeg,k),1);
        }
        ibeg=max(kk,ibeg);
        F77NAME(daxpy)(min(kk+nsub,dim-1)-ibeg+1,B(kk,k),
          addr(ibeg,kk),1,P->addr(ibeg,k),1);
      }
      for (int j=k+1;j<=min(k+P->supDiags(),dim-1);j++) { // s^T B
        int kbeg=max(k+1,j-B.supDiags());
        int kend=min(dim-1,min(k+nsub,j+B.subDiags()));
        if (kbeg<=kend) {
          (*P)(k,j)=F77NAME(ddot)(kend-kbeg+1,addr(kbeg,k),1,
            B.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^T b
        F77NAME(ddot)(min(nb,dim-k-1),addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(dger)(min(dim-k,nt),min(dim-k,B.supDiags()+1),double_one_,
      addr(k,k),1,B.addr(k,k),1,P->addr(k,k),P->bands()-1);
  }
  return P;
}

template<> BandMatrix<double,double>* operator*(
const BandMatrix<double,double> &B,
const SymmetricBandMatrix<double,double> &S) {
  int n=S.size(0);
  CHECK_SAME(n,B.size(0));
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(n,S.subDiags()+B.subDiags(),
      S.subDiags()+B.supDiags(),double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-S.subDiags());k<=min(n-1,j+S.subDiags());k++) {
      int ibeg=max(0,k-B.supDiags());
      F77NAME(daxpy)(min(n-1,k+B.subDiags())-ibeg+1,S(k,j),
        B.addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const UpperHessenbergMatrix<double,double> &H) const {
// compute by bordering: note that
// [ sigma s^T ] [     eta_11 h^T ]
// [   s    S  ] [ e_0 eta_21  H  ]
// = [ sigma eta_11 + s^T e_0 eta_21 , sigma h^T + s^T H ]
// = [   s   eta_11 +  S  e_0 eta_21 ,   s   h^T +  S  H ]
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(dim,double_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*H(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nsub>0) {
      F77NAME(daxpy)(min(dim-k-1,nsub+1),H(k+1,k),addr(k+1,k+1),1,
        P->addr(k+1,k),1); // S  e_0 eta_21
      for (int j=k+1;j<dim;j++) { // s^T H
        (*P)(k,j)=
          F77NAME(ddot)(min(dim-k-1,min(j-k+1,nsub)),addr(k+1,k),1,
            H.addr(k+1,j),1);
      }
      (*P)(k,k)=(*this)(k+1,k)*H(k+1,k); // s^T e_0 eta_21
    }
    F77NAME(dger)(min(dim-k,nsub+1),dim-k,double_one_,addr(k,k),1,
      H.addr(k,k),dim,P->addr(k,k),dim);
  }
  return P;
}

template<> SquareMatrix<double,double>* operator*(
const UpperHessenbergMatrix<double,double> &H,
const SymmetricBandMatrix<double,double> &S) {
  int n=S.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-S.subDiags());k<=min(n-1,j+S.subDiags());k++) {
      F77NAME(daxpy)(min(n,k+2),S(k,j),H.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const DiagonalMatrix<double,double> &D) const {
  BandMatrix<double,double> *P=makeBandMatrix();
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsub);
    F77NAME(dscal)(min(dim-1,j+nsub)-ibeg+1,D[j],P->addr(ibeg,j),1);
  }
  return P;
}

template<> BandMatrix<double,double>* operator*(
const DiagonalMatrix<double,double> &D,
const SymmetricBandMatrix<double,double> &S) {
  int n=S.size(0);
  BandMatrix<double,double> *P=S.makeBandMatrix();
  int stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    int jbeg=max(0,i-S.subDiags());
    F77NAME(dscal)(min(n-1,i+S.subDiags())-jbeg+1,D[i],
      P->addr(i,jbeg),stride);
  }
  return P;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const SymmetricTridiagonalMatrix<double,double> &T) const {
// compute by bordering: note that
// [ sigma s^T ] [     tau    , lambda e_0^T ]
// [   s    S  ] [ e_0 lambda ,           T  ]
//   = [ sigma tau + s^T e_0 lambda , sigma lambda e_0^T + s^T T ]
//   = [   s   tau +  S  e_o lambda ,   s   lambda e_0^T +  S  T ]
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(dim,nsub+1,nsub+1,double_zero_);
  (*P)(dim-2,dim-2)=(*this)(dim-2,dim-2)+T.diagonalValue(dim-2);
  (*P)(dim-1,dim-2)=(*this)(dim-1,dim-1)*T.lowerDiagonalValue(dim-2);
  (*P)(dim-2,dim-1)=(*this)(dim-2,dim-2)+T.upperDiagonalValue(dim-2);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*T.diagonalValue(dim-1);
  if (nsub>0) {
    (*P)(dim-2,dim-2)+=(*this)(dim-1,dim-2)*T.lowerDiagonalValue(dim-2);
    (*P)(dim-1,dim-2)+=(*this)(dim-1,dim-2)*T.diagonalValue(dim-2);
    (*P)(dim-2,dim-1)+=(*this)(dim-1,dim-2)*T.diagonalValue(dim-1);
    (*P)(dim-1,dim-1)+=(*this)(dim-1,dim-2)*T.upperDiagonalValue(dim-2);
  }
  for (int k=dim-3;k>=0;k--) {
    if (nsub>0) {
      for (int j=k+1;j<=min(k+nsub+1,dim-1);j++) { // s^T T
        if (j-1>=k+1) (*P)(k,j)=(*this)(k,j-1)*T.upperDiagonalValue(j-1);
        if (j<=k+nsub) (*P)(k,j)+=(*this)(k,j)*T.diagonalValue(j);
        if (j+1<=k+nsub) {
          (*P)(k,j)+=(*this)(k,j+1)*T.lowerDiagonalValue(j);
        }
      }
      (*P)(k,k)=(*this)(k+1,k)*T.lowerDiagonalValue(k); // s^T e_0 lambda
    }
    F77NAME(daxpy)(min(dim-k-1,nsub+1),T.lowerDiagonalValue(k),
      addr(k+1,k+1),1,P->addr(k+1,k),1); // S e_0 lambda
    F77NAME(daxpy)(min(dim-k,nsub+1),T.diagonalValue(k),addr(k,k),1,
      P->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(daxpy)(min(dim-k,nsub+1),T.upperDiagonalValue(k),addr(k,k),1,
      P->addr(k,k+1),1);
     // [ sigma ] lambda
     // [   s   ]
  }
  return P;
}

template<> BandMatrix<double,double>* operator*(
const SymmetricTridiagonalMatrix<double,double> &T,
const SymmetricBandMatrix<double,double> &S) {
// compute by bordering: note that
// [     tau    , lambda e_0^T ] [ sigma S^T ]
// [ e_0 lambda ,      T       ] [   s    S  ]
//   = [ tau sigma + lambda e_0^T s , tau s^T        + lambda e_0^T S ]
//   = [ e_0 lambda sigma +     T s , e_0 lambda s^T +            T S ]
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=S.subDiags();
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(n,nb+1,nb+1,double_zero_);
  (*P)(n-2,n-2)=T.diagonalValue(n-2)*S(n-2,n-2);
  (*P)(n-1,n-2)=T.lowerDiagonalValue(n-2)*S(n-2,n-2);
  (*P)(n-2,n-1)=T.lowerDiagonalValue(n-2)*S(n-1,n-1);
  (*P)(n-1,n-1)=T.diagonalValue(n-1)*S(n-1,n-1);
  if (nb>0) {
    (*P)(n-2,n-2)+=T.lowerDiagonalValue(n-2)*S(n-1,n-2);
    (*P)(n-1,n-2)+=T.diagonalValue(n-1)*S(n-1,n-2);
    (*P)(n-2,n-1)+=T.diagonalValue(n-2)*S(n-1,n-2);
    (*P)(n-1,n-1)+=T.lowerDiagonalValue(n-2)*S(n-1,n-2);
  }
  int stride=P->bands()-1;
  for (int k=n-3;k>=0;k--) {
    if (nb>0) {
      for (int i=k+1;i<=min(n-1,k+nb+1);i++) { // T s
        if (i-1>=k+1) (*P)(i,k)=T.lowerDiagonalValue(i-1)*S(i-1,k); 
        if (i<=k+nb) (*P)(i,k)+=T(i,i)*S(i,k); 
        if (i+1<=min(k+nb,n-1)) {
          (*P)(i,k)+=T.lowerDiagonalValue(i)*S(i+1,k); 
        }
      }
      (*P)(k,k)=T.lowerDiagonalValue(k)*S(k+1,k);
    }
    F77NAME(daxpy)(min(n-k-1,nb+1),T.lowerDiagonalValue(k),
      S.addr(k+1,k+1),1,P->addr(k,k+1),stride);
      // lambda e_0^T S = lambda ( S e_0 )^T
    F77NAME(daxpy)(min(n-k,nb+1),T.diagonalValue(k),
      S.addr(k,k),1,P->addr(k,k),stride);
      // tau [ sigma , s^T ]
    F77NAME(daxpy)(min(n-k,nb+1),T.lowerDiagonalValue(k),S.addr(k,k),1,
      P->addr(k+1,k),stride);
      // e_0 lambda [ sigma , s^T ]
  }
  return P;
}

template<> BandMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const TridiagonalMatrix<double,double> &T) const {
// compute by bordering: note that
// [ sigma s^T ] [    tau     upsilon e_0^T ]
// [   s    S  ] [ e_0 lambda       T       ]
//   = [ sigma tau + s^T e_0 lambda , sigma upsilon e_0^T + s^T T ]
//   = [     s tau +   S e_0 lambda ,     s upsilon e_0^T +   S T ]
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(dim,nsub+1,nsub+1,double_zero_);
  (*P)(dim-2,dim-2)=(*this)(dim-2,dim-2)*T(dim-2,dim-2);
  (*P)(dim-1,dim-2)=(*this)(dim-1,dim-1)*T(dim-1,dim-2);
  (*P)(dim-2,dim-1)=(*this)(dim-2,dim-2)*T(dim-2,dim-1);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*T(dim-1,dim-1);
  if (nsub>0) {
    (*P)(dim-2,dim-2)+=(*this)(dim-1,dim-2)*T(dim-1,dim-2);
    (*P)(dim-1,dim-2)+=(*this)(dim-1,dim-2)*T(dim-2,dim-2);
    (*P)(dim-2,dim-1)+=(*this)(dim-1,dim-2)*T(dim-1,dim-1);
    (*P)(dim-1,dim-1)+=(*this)(dim-1,dim-2)*T(dim-2,dim-1);
  }
  for (int k=dim-3;k>=0;k--) {
    if (nsub>0) {
      for (int j=k+1;j<=min(k+nsub+1,dim-1);j++) { // s^T T
        if (j-1>=k+1) (*P)(k,j)=(*this)(k,j-1)*T(j-1,j);
        if (j<=k+nsub) (*P)(k,j)+=(*this)(k,j)*T(j,j);
        if (j+1<=k+nsub) (*P)(k,j)+=(*this)(k,j+1)*T(j+1,j);
      }
      (*P)(k,k)=(*this)(k+1,k)*T(k+1,k); // s^T e_0 lambda
    }
    F77NAME(daxpy)(min(dim-k-1,nsub+1),T(k+1,k),addr(k+1,k+1),1,
      P->addr(k+1,k),1); // S e_0 lambda
    F77NAME(daxpy)(min(dim-k,nsub+1),T(k,k),addr(k,k),1,P->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(daxpy)(min(dim-k,nsub+1),T(k,k+1),addr(k,k),1,
      P->addr(k,k+1),1);
      // [ sigma ] lambda
      // [   s   ]
  }
  return P;
}

template<> BandMatrix<double,double>* operator*(
const TridiagonalMatrix<double,double> &T,
const SymmetricBandMatrix<double,double> &S) {
// compute by bordering: note that
// [     tau    , upsilon e_0^T ] [ sigma S^T ]
// [ e_0 lambda ,       T       ] [   s    S  ]
//   = [ tau sigma + upsilon e_0^T s , tau s^T        + upsilon e_0^T S ]
//   = [ e_0 lambda sigma +      T s , e_0 lambda s^T +             T S ]
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=S.subDiags();
  BandMatrix<double,double> *P=OPERATOR_NEW
    BandMatrix<double,double>(n,nb+1,nb+1,double_zero_);
  (*P)(n-2,n-2)=T(n-2,n-2)*S(n-2,n-2);
  (*P)(n-1,n-2)=T(n-1,n-2)*S(n-2,n-2);
  (*P)(n-2,n-1)=T(n-2,n-1)*S(n-1,n-1);
  (*P)(n-1,n-1)=T(n-1,n-1)*S(n-1,n-1);
  if (nb>0) {
    (*P)(n-2,n-2)+=T(n-2,n-1)*S(n-1,n-2);
    (*P)(n-1,n-2)+=T(n-1,n-1)*S(n-1,n-2);
    (*P)(n-2,n-1)+=T(n-2,n-2)*S(n-1,n-2);
    (*P)(n-1,n-1)+=T(n-1,n-2)*S(n-1,n-2);
  }
  int stride=P->bands()-1;
  for (int k=n-3;k>=0;k--) {
    if (nb>0) {
      for (int i=k+1;i<=min(n-1,k+nb+1);i++) { // T s
        if (i-1>=k+1) (*P)(i,k)=T(i,i-1)*S(i-1,k); 
        if (i<=k+nb) (*P)(i,k)+=T(i,i)*S(i,k); 
        if (i+1<=min(k+nb,n-1)) (*P)(i,k)+=T(i,i+1)*S(i+1,k); 
      }
      (*P)(k,k)=T(k,k+1)*S(k+1,k);
    }
    F77NAME(daxpy)(min(n-k-1,nb+1),T(k,k+1),S.addr(k+1,k+1),1,
      P->addr(k,k+1),stride);
      // upsilon e_0^T S = upsilon ( S e_0 )^T
    F77NAME(daxpy)(min(n-k,nb+1),T(k,k),S.addr(k,k),1,
      P->addr(k,k),stride);
      // tau [ sigma , s^T ]
    F77NAME(daxpy)(min(n-k,nb+1),T(k+1,k),S.addr(k,k),1,
      P->addr(k+1,k),stride);
      // e_0 lambda [ sigma , s^T ]
  }
  return P;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const SymmetricMatrix<double,double> &S) const {
// compute by bordering: note that
// [ sigma s^T ] [ tau t^T ] = [ sigma tau + s^T t , sigma t^T + s^T T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^T +  S  T ]
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(dim,double_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*S(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nsub>0) {
      F77NAME(dsbmv)('L',dim-k-1,nsub,double_one_,addr(k+1,k+1),nt,
        S.addr(k+1,k),1,double_zero_,P->addr(k+1,k),1); // S t
      for (int j=k+1;j<dim;j++) { // s^T T
        int kend=min(j-1,k+nsub);
        if (k<kend) {
          (*P)(k,j)=F77NAME(ddot)(kend-k,addr(k+1,k),1,S.addr(j,k+1),dim);
        }
        int kbeg=max(k+1,j);
        kend=min(dim-1,k+nsub);
        if (kbeg<=kend) {
          (*P)(k,j)+=
            F77NAME(ddot)(kend-kbeg+1,addr(kbeg,k),1,S.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^T t
        F77NAME(ddot)(min(nsub,dim-k-1),addr(k+1,k),1,S.addr(k+1,k),1);
    }
    F77NAME(dger)(min(dim-k,nt),dim-k,double_one_,addr(k,k),1,
      S.addr(k,k),1,P->addr(k,k),dim);
  }
  return P;
}

template<> SquareMatrix<double,double>* operator*(
const SymmetricMatrix<double,double> &S,
const SymmetricBandMatrix<double,double> &B) {
// compute by bordering: note that
// [ sigma s^T ] [ tau t^T ] = [ sigma tau + s^T t , sigma t^T + s^T T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^T +  S  T ]
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  (*P)(n-1,n-1)=S(n-1,n-1)*B(n-1,n-1);
  int nb=B.subDiags();
  for (int k=n-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(n-1,k+nb);kk++) { // S t 
        int ibeg=k+1;
        if (kk>ibeg) {
          F77NAME(daxpy)(kk-ibeg,B(kk,k),S.addr(kk,ibeg),n,
            P->addr(ibeg,k),1);
        }
        F77NAME(daxpy)(n-kk,B(kk,k),S.addr(kk,kk),1,P->addr(kk,k),1);
      }
      for (int j=k+1;j<n;j++) { // s^T T
        int kbeg=max(k+1,j-nb);
        if (kbeg<j) {
          (*P)(k,j)=F77NAME(ddot)(j-kbeg,S.addr(kbeg,k),1,
            B.addr(j,kbeg),B.bands()-1);
        }
        kbeg=max(kbeg,j);
        (*P)(k,j)+=F77NAME(ddot)(min(n-1,j+nb)-kbeg+1,
          S.addr(kbeg,k),1,B.addr(kbeg,j),1);
      }
      (*P)(k,k)= // s^T t
        F77NAME(ddot)(min(nb,n-k-1),S.addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(dger)(n-k,min(n-k,B.bands()),double_one_,
      S.addr(k,k),1,B.addr(k,k),1,P->addr(k,k),n);
  }
  return P;
}

template<> Matrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const UpperTrapezoidalMatrix<double,double> &U) const {
  CHECK_SAME(dim,U.size(0));
  int n=U.size(1);
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(dim,n,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,dim-1);k++) {
        int ibeg=max(0,k-nsub);
        if (ibeg<k) {
          F77NAME(daxpy)(k-ibeg,U(k,j),addr(k,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,k);
        F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,dim);k++) {
        int ibeg=max(0,k-nsub);
        if (ibeg<k) {
          F77NAME(daxpy)(k-ibeg,U(k,j),addr(k,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,k);
        F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
      if (j<dim) {
        int ibeg=max(0,j-nsub);
        if (ibeg<j) {
          F77NAME(daxpy)(j-ibeg,double_one_,addr(j,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,j);
        F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,double_one_,
          addr(ibeg,j),1,P->addr(ibeg,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<double,double>* operator*(
const UpperTrapezoidalMatrix<double,double> &U,
const SymmetricBandMatrix<double,double> &B) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,double>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(daxpy)(min(k+1,m),B(k,j),U.addr(0,k),1,P->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(daxpy)(min(k,m),B(k,j),U.addr(0,k),1,P->addr(0,j),1);
        if (k<m) (*P)(k,j)+=B(k,j);
      }
    }
  }
  return P;
}

template<> Matrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const LowerTrapezoidalMatrix<double,double> &L) const {
  CHECK_SAME(dim,L.size(0));
  int n=L.size(1);
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(dim,n,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<dim;k++) {
        int ibeg=max(0,k-nsub);
        if (ibeg<k) {
          F77NAME(daxpy)(k-ibeg,L(k,j),addr(k,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,k);
        F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-nsub);
      if (ibeg<j) {
        F77NAME(daxpy)(j-ibeg,double_one_,addr(j,ibeg),nt-1,
          P->addr(ibeg,j),1);
      }
      ibeg=max(ibeg,j);
      F77NAME(daxpy)(min(dim-1,j+nsub)-ibeg+1,double_one_,addr(ibeg,j),1,
        P->addr(ibeg,j),1);
      for (int k=j+1;k<dim;k++) {
        int ibeg=max(0,k-nsub);
        if (ibeg<k) {
          F77NAME(daxpy)(k-ibeg,L(k,j),addr(k,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,k);
        F77NAME(daxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<double,double>* operator*(
const LowerTrapezoidalMatrix<double,double> &L,
const SymmetricBandMatrix<double,double> &B) {
int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,double>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(daxpy)(m-k,B(k,j),L.addr(k,k),1,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        (*P)(k,j)+=B(k,j);
        if (k+1<m) {
          F77NAME(daxpy)(m-k-1,B(k,j),L.addr(k+1,k),1,P->addr(k+1,j),1);
        }
      }
    }
  }
  return P;
}

template<> SquareMatrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const SquareMatrix<double,double> &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(dim);
  for (int j=0;j<dim;j++) {
    F77NAME(dsbmv)('L',dim,nsub,double_one_,addr(),nt,
      S.addr(0,j),1,double_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> SquareMatrix<double,double>* operator*(
const SquareMatrix<double,double> &S,
const SymmetricBandMatrix<double,double> &B) {
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,double> *P=
    OPERATOR_NEW SquareMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(daxpy)(n,B(k,j),S.addr(0,j),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Matrix<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const Matrix<double,double> &M) const {
  CHECK_SAME(dim,M.size(0));
  int n=M.size(1);
  Matrix<double,double> *P=OPERATOR_NEW Matrix<double,double>(dim,n);
  for (int j=0;j<dim;j++) {
    F77NAME(dsbmv)('L',dim,nsub,double_one_,addr(),nt,
      M.addr(0,j),1,double_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<double,double>* operator*(
const Matrix<double,double> &M,
const SymmetricBandMatrix<double,double> &B) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,double> *P=
    OPERATOR_NEW Matrix<double,double>(m,n,double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(daxpy)(m,B(k,j),M.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Vector<double,double>*
SymmetricBandMatrix<double,double>::operator*(
const Vector<double,double> &v) const {
  CHECK_SAME(dim,v.size());
  Vector<double,double> *p=OPERATOR_NEW Vector<double,double>(dim);
  F77NAME(dsbmv)('L',dim,nsub,double_one_,addr(),nt,v.addr(),1,
    double_zero_,p->addr(),1);
  return p;
}

template<> void SymmetricBandMatrix<double,double>::sbmv(double alpha,
const Vector<double,double> &x,double beta,Vector<double,double> &b)
const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  F77NAME(dsbmv)('L',dim,nsub,alpha,addr(),nt,x.addr(),1,beta,
    b.addr(),1);
}

template<> void SymmetricBandMatrix<double,double>::sbmm(double alpha,
const Matrix<double,double> &X,double beta,Matrix<double,double> &B,
char side) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      F77NAME(dsbmv)('L',dim,nsub,alpha,addr(),nt,X.addr(0,j),1,beta,
        B.addr(0,j),1);
    }
  } else {
    CHECK_SAME(n,dim);
    if (abs(beta)==double_zero_) B=double_zero_;
    else B*=beta;
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-nsub);k<=min(dim-1,j+nsub);k++) {
        F77NAME(daxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
          B.addr(0,j),1);
      }
    }
  }
}

template<> double SymmetricBandMatrix<double,double>::normFrobenius()
const {
  double *work=0;
  return F77NAME(dlansb)('F','L',dim,nsub,addr(),nt,work);
}

template<> double SymmetricBandMatrix<double,double>::normInfinity()
const {
  double *work=OPERATOR_NEW_BRACKET(double,dim);;
  double val=F77NAME(dlansb)('I','L',dim,nsub,addr(),nt,work);
  delete [] work; work=0;
  return val;
}

template<> double SymmetricBandMatrix<double,double>::normMaxEntry()
const {
  double *work=0;
  return F77NAME(dlansb)('M','L',dim,nsub,addr(),nt,work);
}

template<> double SymmetricBandMatrix<double,double>::normOne() const {
  double *work=OPERATOR_NEW_BRACKET(double,dim);;
  double val=F77NAME(dlansb)('O','L',dim,nsub,addr(),nt,work);
  delete [] work; work=0;
  return val;
}

template<> Vector<double,double>* 
SymmetricBandMatrix<double,double>::eigenvalues(
OrthogonalMatrix<double,double> *&Q) const {
  if (Q!=0) CHECK_SAME(dim,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Vector<double,double> *lambda=OPERATOR_NEW Vector<double,double>(dim);
  SymmetricBandMatrix<double,double> *copy=
    OPERATOR_NEW SymmetricBandMatrix<double,double>(*this);
  double *work=OPERATOR_NEW_BRACKET(double,3*dim-2);
  int info;
  double *qa=( Q==0 ? 0 : Q->addr() );
  F77NAME(dsbev)(jobz,'L',dim,nsub,copy->addr(),copy->bands(),
    lambda->addr(),qa,dim,work,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete copy; copy=0;
  return lambda;
}

template class SymmetricBandMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveBandMatrix<double,double>*
SymmetricPositiveBandMatrix<double,double>::operator+(
const SymmetricPositiveBandMatrix<double,double> &B) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,B.size(0));
  SymmetricPositiveBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricPositiveBandMatrix<double,double>(n,max(ns,B.subDiags()),
    double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(daxpy)(min(B.bands(),n-j),double_one_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricPositiveBandMatrix<double,double>*
SymmetricPositiveBandMatrix<double,double>::operator+(
const SymmetricPositiveTridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricPositiveBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricPositiveBandMatrix<double,double>(n,max(ns,1),
    double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=T(j,j);
    if (j+1<n) (*S)(j+1,j)+=T(j+1,j);
  }
  return S;
}

template<> SymmetricBandMatrix<double,double>*
SymmetricPositiveBandMatrix<double,double>::operator+(
const SymmetricTridiagonalMatrix<double,double> &T) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricBandMatrix<double,double> *S=OPERATOR_NEW
    SymmetricBandMatrix<double,double>(n,max(ns,1),double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=T(j,j);
    if (j+1<n) (*S)(j+1,j)+=T(j+1,j);
  }
  return S;
}

template<> SymmetricPositiveMatrix<double,double>*
SymmetricPositiveBandMatrix<double,double>::operator+(
const SymmetricPositiveMatrix<double,double> &M) const {
  int n=size(0);
  int nb=bands();
  CHECK_SAME(n,M.size(0));
  SymmetricPositiveMatrix<double,double> *S=OPERATOR_NEW
    SymmetricPositiveMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(daxpy)(n-j,double_one_,M.addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricMatrix<double,double>*
SymmetricPositiveBandMatrix<double,double>::operator+(
const SymmetricMatrix<double,double> &T) const {
  int n=size(0);
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,double> *S=OPERATOR_NEW
    SymmetricMatrix<double,double>(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(daxpy)(n-j,double_one_,T.addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> double
SymmetricPositiveBandMatrix<double,double>::equilibrate(
Vector<double,double> &s,double &scond) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,s.size());
  double amax;
  int info;
  F77NAME(dpbequ)('L',n,ns,addr(),nb,s.addr(),scond,amax,info);
  CHECK_TEST(info==0);
  return amax;
}

template<> double
SymmetricPositiveBandMatrix<double,double>::reciprocalConditionNumber()
const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  SymmetricPositiveBandMatrix<double,double> *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<double,double>(*this);
  int info;
  F77NAME(dpbtrf)('L',n,ns,BF->addr(),nb,info);
  CHECK_TEST(info==0);

  double anorm=normOne();
  double rcond=numeric_limits<double>::infinity();
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(dpbcon)('L',n,ns,BF->addr(),nb,anorm,rcond,work,iwork,info);
  delete [] work; work=0;
  delete [] iwork; iwork=0;
  delete BF; BF=0;
  return rcond;
}

template<> void SymmetricPositiveBandMatrix<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,char) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  SymmetricPositiveBandMatrix<double,double> *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<double,double>(*this);
  x.copy(b);
  int info;
  F77NAME(dpbsv)('L',n,ns,1,BF->addr(),nb,x.addr(),n,info);
  CHECK_TEST(info==0);
  delete BF; BF=0;
}

template<> void SymmetricPositiveBandMatrix<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,char side,char)
const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  SymmetricPositiveBandMatrix<double,double> *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<double,double>(*this);
  int info;
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    CHECK_SAME(nrhs,X.size(1))
    X.copy(B);
    F77NAME(dpbsv)('L',n,ns,1,BF->addr(),nb,X.addr(),n,info);
    CHECK_TEST(info==0);
  } else {
    int nrhs=B.size(0);
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    CHECK_SAME(nrhs,X.size(0))
    X.copy(B);
    F77NAME(dpbtrf)('L',n,ns,BF->addr(),nb,info);
    CHECK_TEST(info==0);
    for (int j=0;j<nrhs;j++) { // dpbtrs loop 20:
      F77NAME(dtbsv)('L','N','N',n,ns,BF->addr(),nb,X.addr(0,j),nrhs);
      F77NAME(dtbsv)('L','T','N',n,ns,BF->addr(),nb,X.addr(0,j),nrhs);
    }
  }
  delete BF; BF=0;
}

template class SymmetricPositiveBandMatrix<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SpecializedMatrix.C"
template class CompanionMatrix<double,double>; 
template class HadamardMatrix<double,double>; 
template class HankelMatrix<double,double>; 
template class HilbertMatrix<double,double>; 
template class KahanMatrix<double,double>; 
template class LauchliMatrix<double,double>; 
template class PascalMatrix<double,double>; 
template class RosserMatrix<double,double>; 
template class ToeplitzMatrix<double,double>; 
template class SymmetricToeplitzMatrix<double,double>; 
template class VandermondeMatrix<double,double>; 
template class WilkinsonMatrix<double,double>; 
template void testSpecializedMatrix(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "GaussianFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//template<> char GaussianFactorization<double,double>::norm='1';

template<>
GaussianFactorization<double,double,SquareMatrix<double,double> >::
GaussianFactorization(const SquareMatrix<double,double>& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),piv_op(po),equ_op(eo),
equed('N'),colcnd(numeric_limits<double>::infinity()),
rowcnd(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  LU=OPERATOR_NEW SquareMatrix<double,double>(n);
  LU->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    r=OPERATOR_NEW Vector<double,double>(n);
    c=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    F77NAME(dgeequ)(n,n,LU->addr(),n,r->addr(),c->addr(),rowcnd,colcnd,
      amax,info);
    if (info==0) {
      switch (equ_op) {
        case Factorization::EQUILIBRATE_ROWS:
          equed='R';
          break;
        case Factorization::EQUILIBRATE_COLUMNS:
          equed='C';
          break;
        case Factorization::EQUILIBRATE_ROWS_AND_COLUMNS:
          equed='B';
          break;
        default:
          equed='N';
      }
      F77NAME(dlaqge)(n,n,LU->addr(),n,r->addr(),c->addr(),
        rowcnd,colcnd,amax,equed);
      switch (equed) {
        case 'R':
          equ_op=Factorization::EQUILIBRATE_ROWS;
          break;
        case 'C':
          equ_op=Factorization::EQUILIBRATE_COLUMNS;
          break;
        case 'B':
          equ_op=Factorization::EQUILIBRATE_ROWS_AND_COLUMNS;
          break;
        default:
          equ_op=Factorization::NO_EQUILIBRATION;
      }
    }
  }
  anormi=LU->normInfinity();
  anormm=LU->normMaxEntry();
  anormo=LU->normOne();

  int info=0;
  switch (piv_op) {
    case Factorization::NO_PIVOTING: {
      for (int k=0;k<n;k++) {
        if (abs((*LU)(k,k))==double_zero_) {
          info=k+1;
          break;
        }
	double alpha=double_one_/(*LU)(k,k);
	int nrows=n-k-1;
	if (nrows>0) {
	  F77NAME(dscal)(nrows,alpha,LU->addr(k+1,k),1);
	  for (int j=k+1;j<n;j++) {
	    alpha=-(*LU)(k,j);
	    F77NAME(daxpy)(nrows,alpha,LU->addr(k+1,k),1,
			   LU->addr(k+1,j),1);
	  }
	}
      }
      break;
    }
    case Factorization::PIVOT_ROWS_AND_COLUMNS: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      jpiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(dgetc2)(n,LU->addr(),n,ipiv,jpiv,info);
      break;
    }
    default: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(dgetrf)(n,n,LU->addr(),n,ipiv,info);
      break;
    }
  }
  double *work=0;
  if (info>0) { // see dgesvx
    rpvgrw=F77NAME(dlantr)('M','U','N',info,info,LU->addr(),n,work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
    cerr << "zero pivot in GaussianFactorization::GaussianFactorization"
      << "\n zero pivot number,reciprocal pivot growth = " << info
      << " " << rpvgrw << endl;
    if (LU) delete LU; LU=0;
    if (r) delete r; r=0;
    if (c) delete c; c=0;
    if (ipiv) delete [] ipiv;ipiv=0;
    if (jpiv) delete [] jpiv;jpiv=0;
    A_original=0;
  } else {
    rpvgrw=F77NAME(dlantr)('M','U','N',n,n,LU->addr(),n,work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
  }
}

template<> void
GaussianFactorization<double,double,SquareMatrix<double,double> >::solve(
const Vector<double,double> &b,Vector<double,double> &x,
Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
//constructor factored   R A C = Q L U P^T
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,LU->size(0));
  if (&x!=&b) x.copy(b);
  if (to==Factorization::NO_TRANSPOSE) {
    // A x = b ==> L U P^T C^{-1} x = Q^T R b
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ri=r->addr();
      double *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
    if (ipiv!=0) F77NAME(dlaswp)(1,x.addr(),n,1,n-1,ipiv,1);
    F77NAME(dtrsm)('L','L','N','U',n,1,double_one_,LU->addr(),n,
      x.addr(),n);
    F77NAME(dtrsm)('L','U','N','N',n,1,double_one_,LU->addr(),n,
      x.addr(),n);
    if (jpiv!=0) F77NAME(dlaswp)(1,x.addr(),n,1,n-1,jpiv,-1);
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ci=c->addr();
      double *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
  } else { // A^T X = B ==> U^T L^T Q^T R^{-1} X = P^T C^T B
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ci=c->addr();
      double *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
    if (jpiv!=0) F77NAME(dlaswp)(1,x.addr(),n,1,n-1,jpiv,1);
    F77NAME(dtrsm)('L','U','T','N',n,1,double_one_,LU->addr(),n,
      x.addr(),n);
    F77NAME(dtrsm)('L','L','T','U',n,1,double_one_,LU->addr(),n,
      x.addr(),n);
    if (ipiv!=0) F77NAME(dlaswp)(1,x.addr(),n,1,n-1,ipiv,-1);
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ri=r->addr();
      double *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
  }
}

template<> void
GaussianFactorization<double,double,SquareMatrix<double,double> >::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
  CHECK_POINTER(LU);
//constructor factored   R A C = Q L U P^T
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (&X!=&B) X.copy(B);
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) {
//    A X = B ==> L U P^T C^{-1} X = Q^T R B
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ri=r->addr();
          double *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
      if (ipiv!=0) F77NAME(dlaswp)(n,X.addr(),m,1,m-1,ipiv,1);
      F77NAME(dtrsm)('L','L','N','U',m,n,double_one_,
                     LU->addr(),m,X.addr(),m);
      F77NAME(dtrsm)('L','U','N','N',m,n,double_one_,
                     LU->addr(),m,X.addr(),m);
      if (jpiv!=0) F77NAME(dlaswp)(n,X.addr(),m,1,m-1,jpiv,-1);
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ci=c->addr();
          double *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
    } else { // A^T X = B ==> U^T L^T Q^T R^{-1} X = P^T C^T B
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ci=c->addr();
          double *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
      if (jpiv!=0) F77NAME(dlaswp)(n,X.addr(),m,1,m-1,jpiv,1);
      F77NAME(dtrsm)('L','U','T','N',m,n,double_one_,
                     LU->addr(),m,X.addr(),m);
      F77NAME(dtrsm)('L','L','T','U',m,n,double_one_,
                     LU->addr(),m,X.addr(),m);
      if (ipiv!=0) F77NAME(dlaswp)(n,X.addr(),m,1,m-1,ipiv,-1);
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ri=r->addr();
          double *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
    }
  } else {
    CHECK_SAME(n,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) {
      // X A = B ==> X R^{-1} Q L U = B C P
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(dscal)(m,(*c)[j],X.addr(0,j),1);
      }
      if (jpiv!=0) {
        for (int j=n-2;j>=0;j--) {
          cout << "jpiv[" << j << "] = " << jpiv[j] << endl;
          if (j!=jpiv[j]-1) {
            F77NAME(dswap)(m,X.addr(0,j),1,X.addr(0,jpiv[j]-1),1);
          }
        }
      }
      F77NAME(dtrsm)('R','U','N','N',m,n,double_one_,
                     LU->addr(),n,X.addr(),m);
      F77NAME(dtrsm)('R','L','N','U',m,n,double_one_,
                     LU->addr(),n,X.addr(),m);
      if (ipiv!=0) {
        for (int i=0;i<n-1;i++) {
          if (i!=ipiv[i]-1) {
            F77NAME(dswap)(m,X.addr(0,i),1,X.addr(0,ipiv[i]-1),1);
          }
        }
      }
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(dscal)(m,(*r)[j],X.addr(0,j),1);
      }
    } else { // X A^T = B ==> X C^{-1} P U^T L^T = B R Q
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(dscal)(m,(*r)[j],X.addr(0,j),1);
      }
      if (ipiv!=0) {
        for (int i=n-2;i>=0;i--) {
          if (i!=ipiv[i]-1) {
            F77NAME(dswap)(m,X.addr(0,i),1,X.addr(0,ipiv[i]-1),1);
          }
        }
      }
      F77NAME(dtrsm)('R','L','T','U',m,n,double_one_,
                     LU->addr(),n,X.addr(),m);
      F77NAME(dtrsm)('R','U','T','N',m,n,double_one_,
                     LU->addr(),n,X.addr(),m);
      if (jpiv!=0) {
        for (int j=0;j<n-1;j++) {
          if (j!=jpiv[j]-1) {
            F77NAME(dswap)(m,X.addr(0,j),1,X.addr(0,jpiv[j]-1),1);
          }
        }
      }
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(dscal)(m,(*c)[j],X.addr(0,j),1);
      }
    }
  }
}

template<> double
GaussianFactorization<double,double,SquareMatrix<double,double> >::
reciprocalConditionNumber(Factorization::CONDITION_NUMBER_NORM cnn) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  double anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,4*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  F77NAME(dgecon)(norm,n,LU->addr(),n,anorm,rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void
GaussianFactorization<double,double,SquareMatrix<double,double> >::
improve(const Vector<double,double> &b,Vector<double,double> &x,
double &berr,double &ferr,Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(n+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> residual(n);
  Vector<double,double> work(n);
  Vector<double,double> v(n);
  char trans=(to==Factorization::TRANSPOSE ? 'T' : 'N');
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(dgemv)(trans,n,n,double_mone_,A_original->addr(),n,
      x.addr(),1,double_one_,residual.addr(),1);
    work.copy(b);
    F77_NAME(dla_geamv)(F77NAME(ilatrans)(trans),n,n,double_one_,
      A_original->addr(),n,x.addr(),1,double_one_,work.addr(),1);

    berr=double_zero_;
    if (to==Factorization::NO_TRANSPOSE) {
      if (r!=0) {
        for (int i=0;i<n;i++) {
          berr=(work[i]*(*r)[i]>safe2 ?
            max(berr,abs(residual[i])/work[i]) :
            max(berr,(abs(residual[i])*(*r)[i]+safe1)
                 /(work[i]*(*r)[i]+safe1)));
        }
      } else {
        for (int i=0;i<n;i++) {
          berr=(work[i]>safe2 ? max(berr,abs(residual[i])/work[i]) :
            max(berr,(abs(residual[i])+safe1)/(work[i]+safe1)));
        }
      }
    } else {
      if (c!=0) {
        for (int i=0;i<n;i++) {
          berr=(work[i]*(*c)[i]>safe2 ?
            max(berr,abs(residual[i])/work[i]) :
            max(berr,(abs(residual[i])*(*c)[i]+safe1)
                 /(work[i]*(*c)[i]+safe1)));
        }
      } else {
        for (int i=0;i<n;i++) {
          berr=(work[i]>safe2 ? max(berr,abs(residual[i])/work[i]) :
            max(berr,(abs(residual[i])+safe1)/(work[i]+safe1)));
        }
      }
    }
    if (!(berr>eps && 2.*berr<=lstres && count<=itmax)) break;
    solve(residual,residual,to);
    F77NAME(daxpy)(n,double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  double *residuali=residual.addr();
  double *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(dlacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual,to);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual,to);
    }
  }
  int i=F77NAME(idamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void
GaussianFactorization<double,double,SquareMatrix<double,double> >::
improve(const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &berr,Vector<double,double> &ferr,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
//compare to dgesvx.f
  CHECK_POINTER(LU);
  int k=LU->size(0),m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,k);
    CHECK_SAME(n,berr.size());
    CHECK_SAME(n,ferr.size());
  } else {
    CHECK_SAME(n,k);
    CHECK_SAME(m,berr.size());
    CHECK_SAME(m,ferr.size());
  }
  double nz=static_cast<double>(k+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> x(k);
  Vector<double,double> rhs(k);
  Vector<double,double> residual(k);
  Vector<double,double> work(k);
  Vector<double,double> v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  char trans='N';
  Factorization::TRANSPOSE_OPTION lto=Factorization::NO_TRANSPOSE;
  if ((so==Factorization::LEFT_SIDE)!=(to==Factorization::NO_TRANSPOSE)) {
    trans='T';
    lto=Factorization::TRANSPOSE;
  }
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(dcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(dcopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(dgemv)(trans,k,k,double_mone_,A_original->addr(),k,
        x.addr(),1,double_one_,residual.addr(),1);
      work.copy(rhs);
      F77_NAME(dla_geamv)(F77NAME(ilatrans)(trans),k,k,double_one_,
        A_original->addr(),k,x.addr(),1,double_one_,work.addr(),1);

      double &s=berr[j];
      s=double_zero_;
      if (lto==Factorization::NO_TRANSPOSE) {
        if (r!=0) {
          for (int i=0;i<k;i++) {
            s=(work[i]*(*r)[i]>safe2 ? max(s,abs(residual[i])/work[i]) :
              max(s,(abs(residual[i])*(*r)[i]+safe1)
                   /(work[i]*(*r)[i]+safe1)));
          }
        } else {
          for (int i=0;i<k;i++) {
            s=(work[i]>safe2 ? max(s,abs(residual[i])/work[i]) :
              max(s,(abs(residual[i])+safe1)/(work[i]+safe1)));
          }
        }
      } else {
        if (c!=0) {
          for (int i=0;i<k;i++) {
            s=(work[i]*(*c)[i]>safe2 ? max(s,abs(residual[i])/work[i]) :
              max(s,(abs(residual[i])*(*c)[i]+safe1)
                   /(work[i]*(*c)[i]+safe1)));
          }
        } else {
          for (int i=0;i<k;i++) {
            s=(work[i]>safe2 ? max(s,abs(residual[i])/work[i]) :
              max(s,(abs(residual[i])+safe1)/(work[i]+safe1)));
          }
        }
      }
      if (!(s>eps && 2.*s<=lstres && count<=itmax)) break;
      solve(residual,residual,lto);
      F77NAME(daxpy)(m,double_one_,residual.addr(),1,x.addr(),1);
      lstres=s;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(dcopy)(k,x.addr(),1,X.addr(0,j),1);
      } else F77NAME(dcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(dcopy)(k,x.addr(),1,X.addr(j,0),m);
      } else F77NAME(dcopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    double *residuali=residual.addr();
    double *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(dlacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual,lto);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual,lto);
      }
    }
    int i=F77NAME(idamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr[j]/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<>
GaussianFactorization<double,double,TridiagonalMatrix<double,double> >::
GaussianFactorization(const TridiagonalMatrix<double,double>& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),
piv_op(Factorization::PIVOT_ROWS),equ_op(Factorization::NO_EQUILIBRATION),
equed('N'),colcnd(numeric_limits<double>::infinity()),
rowcnd(numeric_limits<double>::infinity()) {
  int n=A.size(0);
//LU=OPERATOR_NEW SquareMatrix<double,double>(n);
//LU->copy(A);

  anormi=A.normInfinity();
  anormm=A.normMaxEntry();
  anormo=A.normOne();
  rpvgrw=double_one_; // not computed
}

template<> void
GaussianFactorization<double,double,TridiagonalMatrix<double,double> >::
solve(const Vector<double,double> &b,Vector<double,double> &x,
Factorization::TRANSPOSE_OPTION to) {
//constructor did not factor
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  LU=OPERATOR_NEW TridiagonalMatrix<double,double>(n);
  LU->copy(*A_original);
  if (&x!=&b) x.copy(b);
  int info;
  if (to==Factorization::NO_TRANSPOSE) {
    if (piv_op==Factorization::PIVOT_ROWS) {
      F77NAME(dgtsv)(n,1,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),x.addr(),n,info);
    } else {
      F77NAME(dgtsvnp)(n,1,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),x.addr(),n,info);
    }
  } else {
    if (piv_op==Factorization::PIVOT_ROWS) {
      F77NAME(dgtsv)(n,1,LU->upperDiagonalAddr(),LU->diagonalAddr(),
        LU->lowerDiagonalAddr(),x.addr(),n,info);
    } else {
      F77NAME(dgtsvnp)(n,1,LU->upperDiagonalAddr(),LU->diagonalAddr(),
        LU->lowerDiagonalAddr(),x.addr(),n,info);
    }
  }
  CHECK_TEST(info==0);
  delete LU; LU=0;
}

template<> void
GaussianFactorization<double,double,TridiagonalMatrix<double,double> >::
solve(const Matrix<double,double> &B,Matrix<double,double> &X,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
//constructor did not factor
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  LU=OPERATOR_NEW TridiagonalMatrix<double,double>(A_original->size(0));
  LU->copy(*A_original);
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) { // A X = B
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(dgtsv)(m,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(dgtsvnp)(m,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      }
    } else { // A^T X = B
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(dgtsv)(m,n,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(dgtsvnp)(m,n,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      }
    }
  } else {
    CHECK_SAME(n,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) { // X A = B ==> A^T X^T = B^T
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(dgtsvr)(n,m,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(dgtsvrnp)(n,m,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      }
    } else { // X A^T = B ==> A X^T = B^T
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(dgtsvr)(n,m,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(dgtsvrnp)(n,m,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      }
    }
  }
  CHECK_TEST(info==0);
  delete LU; LU=0;
}

template<> double
GaussianFactorization<double,double,TridiagonalMatrix<double,double> >::
reciprocalConditionNumber(Factorization::CONDITION_NUMBER_NORM cnn)
{
  LU=OPERATOR_NEW TridiagonalMatrix<double,double>(A_original->size(0));
  LU->copy(*A_original);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  double anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<double,double> *u2=OPERATOR_NEW Vector<double,double>(n-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(dgttrf)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    F77NAME(dgtcon)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,iwork,
      info);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(dgttrfnp)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    F77NAME(dgtconnp)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,iwork,info);
  }
  CHECK_SAME(info,0)
  delete LU; LU=0;
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void
GaussianFactorization<double,double,TridiagonalMatrix<double,double> >::
improve(const Vector<double,double> &b,Vector<double,double> &x,
double &berr,double &ferr,Factorization::TRANSPOSE_OPTION to) {
  LU=OPERATOR_NEW TridiagonalMatrix<double,double>(A_original->size(0));
  LU->copy(*A_original);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  char norm=(to==Factorization::NO_TRANSPOSE ? 'O' : 'I');
  double anorm=(to==Factorization::NO_TRANSPOSE ? anormo : anormi);
  char trans=(to==Factorization::NO_TRANSPOSE ? 'N' : 'T');
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  double rcond;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<double,double> *u2=OPERATOR_NEW Vector<double,double>(n-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(dgttrf)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    CHECK_TEST(info==0);
    F77NAME(dgtcon)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,iwork,
      info);
    CHECK_TEST(info==0);
    F77NAME(dgtrfs)(trans,n,1,A_original->lowerDiagonalAddr(),
      A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
      LU->lowerDiagonalAddr(),LU->diagonalAddr(),LU->upperDiagonalAddr(),
      u2->addr(),ipiv,b.addr(),n,x.addr(),n,&ferr,&berr,work,iwork,info);
    CHECK_TEST(info==0);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(dgttrfnp)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    CHECK_TEST(info==0);
    F77NAME(dgtconnp)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,iwork,info);
    CHECK_TEST(info==0);
    F77NAME(dgtrfsnp)(trans,n,1,A_original->lowerDiagonalAddr(),
      A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
      LU->lowerDiagonalAddr(),LU->diagonalAddr(),LU->upperDiagonalAddr(),
      b.addr(),n,x.addr(),n,&ferr,&berr,work,iwork,info);
    CHECK_TEST(info==0);
  }
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  delete LU; LU=0;
}

template<> void
GaussianFactorization<double,double,TridiagonalMatrix<double,double> >::
improve(const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &berr,Vector<double,double> &ferr,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
  LU=OPERATOR_NEW TridiagonalMatrix<double,double>(A_original->size(0));
  LU->copy(*A_original);
  int k=LU->size(0),m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,k);
    CHECK_SAME(n,berr.size());
    CHECK_SAME(n,ferr.size());
  } else {
    CHECK_SAME(n,k);
    CHECK_SAME(m,berr.size());
    CHECK_SAME(m,ferr.size());
  }
  char norm=(to==Factorization::NO_TRANSPOSE ? 'O' : 'I');
  double anorm=(to==Factorization::NO_TRANSPOSE ? anormo : anormi);
  char trans=(to==Factorization::NO_TRANSPOSE ? 'N' : 'T');
  double *work=OPERATOR_NEW_BRACKET(double,3*k);
  int *iwork=OPERATOR_NEW_BRACKET(int,k);
  int info;
  double rcond;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<double,double> *u2=OPERATOR_NEW Vector<double,double>(k-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,k);
    F77NAME(dgttrf)(k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    CHECK_TEST(info==0);
    F77NAME(dgtcon)(norm,k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,iwork,
      info);
    CHECK_TEST(info==0);
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(dgtrfs)(trans,k,n,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),u2->addr(),ipiv,B.addr(),m,X.addr(),m,
        ferr.addr(),berr.addr(),work,iwork,info);
    } else {
      F77NAME(dgtrfsr)(trans,k,m,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),u2->addr(),ipiv,B.addr(),m,X.addr(),m,
        ferr.addr(),berr.addr(),work,iwork,info);
    }
    CHECK_TEST(info==0);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(dgttrfnp)(k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    CHECK_TEST(info==0);
    F77NAME(dgtconnp)(norm,k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,iwork,info);
    CHECK_TEST(info==0);
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(dgtrfsnp)(trans,k,n,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),B.addr(),m,X.addr(),m,ferr.addr(),
        berr.addr(),work,iwork,info);
    } else {
      F77NAME(dgtrfsrnp)(trans,k,m,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),B.addr(),m,X.addr(),m,ferr.addr(),
        berr.addr(),work,iwork,info);
    }
    CHECK_TEST(info==0);
  }
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  delete LU; LU=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<>
GaussianFactorization<double,double,BandMatrix<double,double> >::
GaussianFactorization(const BandMatrix<double,double>& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),piv_op(po),equ_op(eo),
equed('N'),colcnd(numeric_limits<double>::infinity()),
rowcnd(numeric_limits<double>::infinity()) {
  int n=A.size(0),nsub=A.subDiags(),nsup=A.supDiags();
  LU=OPERATOR_NEW BandMatrix<double,double>(n,nsub,
    (po==Factorization::NO_PIVOTING ? nsup : nsub+nsup),double_zero_);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-nsup);
    int iend=min(n-1,j+nsub);
    F77NAME(dcopy)(iend-ibeg+1,A.addr(ibeg,j),1,LU->addr(ibeg,j),1);
  }

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    r=OPERATOR_NEW Vector<double,double>(n);
    c=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    F77NAME(dgbequ)(n,n,LU->subDiags(),LU->supDiags(),LU->addr(),
      LU->bands(),r->addr(),c->addr(),rowcnd,colcnd,amax,info);
    if (info==0) {
      switch (equ_op) {
        case Factorization::EQUILIBRATE_ROWS:
          equed='R';
          break;
        case Factorization::EQUILIBRATE_COLUMNS:
          equed='C';
          break;
        case Factorization::EQUILIBRATE_ROWS_AND_COLUMNS:
          equed='B';
          break;
        default:
          equed='N';
      }
      F77NAME(dlaqgb)(n,n,nsub,nsub+nsup,LU->addr(),LU->bands(),r->addr(),
        c->addr(),rowcnd,colcnd,amax,equed);
      switch (equed) {
        case 'R':
          equ_op=Factorization::EQUILIBRATE_ROWS;
          break;
        case 'C':
          equ_op=Factorization::EQUILIBRATE_COLUMNS;
          break;
        case 'B':
          equ_op=Factorization::EQUILIBRATE_ROWS_AND_COLUMNS;
          break;
        default:
          equ_op=Factorization::NO_EQUILIBRATION;
      }
    }
  }
  anormi=LU->normInfinity();
  anormm=LU->normMaxEntry();
  anormo=LU->normOne();

  int info=0;
  switch (piv_op) {
    case Factorization::NO_PIVOTING: {
      F77NAME(dgbtf2np)(n,n,nsub,nsup,LU->addr(),LU->bands(),info);
      break;
    }
    default: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(dgbtrf)(n,n,nsub,nsup,LU->addr(),LU->bands(),ipiv,info);
    }
  }
  double *work=0;
  if (info>0) { // see dgbsvx
    rpvgrw=F77NAME(dlantb)('M','U','N',info,min(info-1,LU->supDiags()),
      LU->addr(),LU->bands(),work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
    cerr << "zero pivot in GaussianFactorization::GaussianFactorization"
       << "\n zero pivot number,reciprocal pivot growth = " << info
       << " " << rpvgrw << endl;
    if (LU) delete LU; LU=0;
    if (r) delete r; r=0;
    if (c) delete c; c=0;
    if (ipiv) delete [] ipiv;ipiv=0;
    A_original=0;
  } else {
    rpvgrw=F77NAME(dlantb)('M','U','N',n,LU->supDiags(),LU->addr(),
      LU->bands(),work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
  }
}

template<> void
GaussianFactorization<double,double,BandMatrix<double,double> >::
solve(const Vector<double,double> &b,Vector<double,double> &x,
Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
//constructor factored R A C = Q L U
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  if (&x!=&b) x.copy(b);
  int info;
  if (to==Factorization::NO_TRANSPOSE) {
    // A x = b ==> L U C^{-1} x = Q^T R b
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ri=r->addr();
      double *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
    if (ipiv!=0) {
      F77NAME(dgbtrs)('N',n,A_original->subDiags(),A_original->supDiags(),
        1,LU->addr(),LU->bands(),ipiv,x.addr(),n,info);
    } else {
      F77NAME(dgbtrsnp)('N',n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),x.addr(),n,info);
    }
    CHECK_TEST(info==0);
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ci=c->addr();
      double *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
  } else { // A^T X = B ==> U^T L^T Q^T R^{-1} X = C^T B
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ci=c->addr();
      double *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
    if (ipiv!=0) {
      F77NAME(dgbtrs)('T',n,A_original->subDiags(),A_original->supDiags(),
        1,LU->addr(),LU->bands(),ipiv,x.addr(),n,info);
    } else {
      F77NAME(dgbtrsnp)('T',n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),x.addr(),n,info);
    }
    CHECK_TEST(info==0);
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ri=r->addr();
      double *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
  }
}

template<> void
GaussianFactorization<double,double,BandMatrix<double,double> >::
solve(const Matrix<double,double> &B,Matrix<double,double> &X,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
  CHECK_POINTER(LU);
//constructor factored   R A C = Q L U P^T
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) {
//    A X = B ==> L U C^{-1} X = Q^T R B
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ri=r->addr();
          double *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
      if (ipiv!=0) {
        F77NAME(dgbtrs)('N',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),ipiv,
          X.addr(),m,info);
      } else {
        F77NAME(dgbtrsnp)('N',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),X.addr(),m,
          info);
      }
      CHECK_TEST(info==0);
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ci=c->addr();
          double *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
    } else { // A^T X = B ==> U^T L^T Q^T R^{-1} X = C^T B
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ci=c->addr();
          double *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
      if (ipiv!=0) {
        F77NAME(dgbtrs)('T',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),ipiv,
          X.addr(),m,info);
      } else {
        F77NAME(dgbtrsnp)('T',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),X.addr(),m,
          info);
      }
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ri=r->addr();
          double *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
    }
  } else {
    CHECK_SAME(n,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) {
      // X A = B ==> X R^{-1} Q L U = B C
      //         ==> U^T L^T Q^T R^{-1} X^T = C^T B^T
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(dscal)(m,(*c)[j],X.addr(0,j),1);
      }
      Vector<double,double> *x=OPERATOR_NEW Vector<double,double>(n);
      for (int i=0;i<m;i++) {
        F77NAME(dcopy)(n,X.addr(i,0),m,x->addr(),1);
        if (ipiv!=0) {
          F77NAME(dgbtrs)('T',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
            x->addr(),n,info);
        } else {
          F77NAME(dgbtrsnp)('T',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),x->addr(),n,
            info);
        }
        F77NAME(dcopy)(n,x->addr(),1,X.addr(i,0),m);
      }
      delete x; x=0;
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(dscal)(m,(*r)[j],X.addr(0,j),1);
      }
    } else {
      // X A^T = B ==> X C^{-1} U^T L^T = B R Q
      //           ==> L U C^{-1} X^T = Q^T R B^T
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(dscal)(m,(*r)[j],X.addr(0,j),1);
      }
      Vector<double,double> *x=OPERATOR_NEW Vector<double,double>(n);
      for (int i=0;i<m;i++) {
        F77NAME(dcopy)(n,X.addr(i,0),m,x->addr(),1);
        if (ipiv!=0) {
          F77NAME(dgbtrs)('N',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
            x->addr(),n,info);
        } else {
          F77NAME(dgbtrsnp)('N',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),x->addr(),n,
            info);
        }
        F77NAME(dcopy)(n,x->addr(),1,X.addr(i,0),m);
      }
      delete x; x=0;
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(dscal)(m,(*c)[j],X.addr(0,j),1);
      }
    }
  }
}

template<> double
GaussianFactorization<double,double,BandMatrix<double,double> >::
reciprocalConditionNumber(Factorization::CONDITION_NUMBER_NORM cnn) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  double anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  if (piv_op==Factorization::PIVOT_ROWS) {
    F77NAME(dgbcon)(norm,n,A_original->subDiags(),A_original->supDiags(),
      LU->addr(),LU->bands(),ipiv,anorm,rcond,work,iwork,info);
  } else {
    F77NAME(dgbconnp)(norm,n,A_original->subDiags(),
      A_original->supDiags(),LU->addr(),LU->bands(),anorm,rcond,work,
      iwork,info);
  }
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void
GaussianFactorization<double,double,BandMatrix<double,double> >::improve(
const Vector<double,double> &b,Vector<double,double> &x,double &berr,
double &ferr,Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(A_original->bands()+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> residual(n);
  Vector<double,double> work(n);
  char trans=(to==Factorization::TRANSPOSE ? 'T' : 'N');
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(dgbmv)(trans,n,n,A_original->subDiags(),
      A_original->supDiags(),double_mone_,A_original->addr(),
      A_original->bands(),x.addr(),1,double_one_,residual.addr(),1);
    work.copy(b);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(dgbamv)(trans,n,n,A_original->subDiags(),
      A_original->supDiags(),double_one_,A_original->addr(),
      A_original->bands(),x.addr(),1,double_one_,work.addr(),1);

    berr=double_zero_;
    if (to==Factorization::NO_TRANSPOSE) {
      if (r!=0) {
        for (int i=0;i<n;i++) {
          berr=(work[i]*(*r)[i]>safe2 ?
            max(berr,abs(residual[i])/work[i]) :
            max(berr,(abs(residual[i])*(*r)[i]+safe1)
                 /(work[i]*(*r)[i]+safe1)));
        }
      } else {
        for (int i=0;i<n;i++) {
          berr=(work[i]>safe2 ? max(berr,abs(residual[i])/work[i]) :
            max(berr,(abs(residual[i])+safe1)/(work[i]+safe1)));
        }
      }
    } else {
      if (c!=0) {
        for (int i=0;i<n;i++) {
          berr=(work[i]*(*c)[i]>safe2 ?
            max(berr,abs(residual[i])/work[i]) :
            max(berr,(abs(residual[i])*(*c)[i]+safe1)
                 /(work[i]*(*c)[i]+safe1)));
        }
      } else {
        for (int i=0;i<n;i++) {
          berr=(work[i]>safe2 ? max(berr,abs(residual[i])/work[i]) :
            max(berr,(abs(residual[i])+safe1)/(work[i]+safe1)));
        }
      }
    }
    if (!(berr>eps && 2.*berr<=lstres && count<=itmax)) break;
    solve(residual,residual,to);
    F77NAME(daxpy)(n,double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }
  double *residuali=residual.addr();
  double *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  Vector<double,double> v(n);
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(dlacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual,to);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual,to);
    }
  }
  int i=F77NAME(idamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void
GaussianFactorization<double,double,BandMatrix<double,double> >::
improve(const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &berr,Vector<double,double> &ferr,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
  CHECK_POINTER(LU);
  int k=LU->size(0),m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,k);
    CHECK_SAME(n,berr.size());
    CHECK_SAME(n,ferr.size());
  } else {
    CHECK_SAME(n,k);
    CHECK_SAME(m,berr.size());
    CHECK_SAME(m,ferr.size());
  }
  double nz=static_cast<double>(A_original->bands()+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> x(k);
  Vector<double,double> rhs(k);
  Vector<double,double> residual(k);
  Vector<double,double> work(k);
  Vector<double,double> v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  char trans='N';
  Factorization::TRANSPOSE_OPTION lto=Factorization::NO_TRANSPOSE;
  if ((so==Factorization::LEFT_SIDE)!=(to==Factorization::NO_TRANSPOSE)) {
    trans='T';
    lto=Factorization::TRANSPOSE;
  }
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(dcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(dcopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(dgbmv)(trans,k,k,A_original->subDiags(),
        A_original->supDiags(),double_mone_,
        A_original->addr(),A_original->bands(),x.addr(),1,double_one_,
        residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77NAME(dgbamv)(trans,k,k,A_original->subDiags(),
        A_original->supDiags(),double_one_,A_original->addr(),
        A_original->bands(),x.addr(),1,double_one_,work.addr(),1);
      double &s=berr[j];
      s=double_zero_;
      if (lto==Factorization::NO_TRANSPOSE) {
        if (r!=0) {
          for (int i=0;i<k;i++) {
            s=(work[i]*(*r)[i]>safe2 ?
              max(s,abs(residual[i])/work[i]) :
              max(s,(abs(residual[i])*(*r)[i]+safe1)
                   /(work[i]*(*r)[i]+safe1)));
          }
        } else {
          for (int i=0;i<k;i++) {
            s=(work[i]>safe2 ? max(s,abs(residual[i])/work[i]) :
              max(s,(abs(residual[i])+safe1)/(work[i]+safe1)));
          }
        }
      } else {
        if (c!=0) {
          for (int i=0;i<k;i++) {
            s=(work[i]*(*c)[i]>safe2 ?
              max(s,abs(residual[i])/work[i]) :
              max(s,(abs(residual[i])*(*c)[i]+safe1)
                   /(work[i]*(*c)[i]+safe1)));
          }
        } else {
          for (int i=0;i<k;i++) {
            s=(work[i]>safe2 ? max(s,abs(residual[i])/work[i]) :
              max(s,(abs(residual[i])+safe1)/(work[i]+safe1)));
          }
        }
      }
      if (!(s>eps && 2.*s<=lstres && count<=itmax)) break;
      solve(residual,residual,to);
      F77NAME(daxpy)(k,double_one_,residual.addr(),1,x.addr(),1);
      lstres=s;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(dcopy)(k,x.addr(),1,X.addr(0,j),1);
      } else F77NAME(dcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(dcopy)(k,x.addr(),1,X.addr(j,0),m);
      } else F77NAME(dcopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    double *residuali=residual.addr();
    double *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(dlacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],
        kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual,to);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual,to);
      }
    }
    int i=F77NAME(idamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template class
  GaussianFactorization<double,double,SquareMatrix<double,double> >;
template void testGaussianFactorization(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "CholeskyFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> CholeskyFactorization<double,double,
SymmetricPositiveMatrix<double,double> >::CholeskyFactorization(
const SymmetricPositiveMatrix<double,double>& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equed('N'),scond(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  L=OPERATOR_NEW SymmetricPositiveMatrix<double,double>(n);
  L->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    F77NAME(dpoequb)(n,L->addr(),n,s->addr(),scond,amax,info);
    equed='N';
    F77NAME(dlaqsy)('L',n,L->addr(),n,s->addr(),scond,amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=L->normInfinity();
  anormm=L->normMaxEntry();
  anormo=L->normOne();

  int info=0;
  F77NAME(dpotrf)('L',n,L->addr(),n,info);
  double *work=OPERATOR_NEW_BRACKET(double,4*n);
  if (info!=0) { // see dposvxx
    rpvgrw=F77_NAME(dla_porpvgrw)('L',info,A_original->addr(),n,
      L->addr(),n, work);
    delete L; L=0;
    if (s!=0) delete s; s=0;
  } else {
    rpvgrw=F77_NAME(dla_porpvgrw)('L',n,A_original->addr(),n,
      L->addr(),n, work);
  }
  delete [] work; work=0;
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveMatrix<double,double> >::solve(
const Vector<double,double> &b,Vector<double,double> &x) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^T
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,L->size(0));
  if (&x!=&b) x.copy(b);
//A x = b ==> L L^T S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    double *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  int info;
  F77NAME(dpotrs)('L',n,1,L->addr(),n,x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    double *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveMatrix<double,double> >::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Factorization::SIDE_OPTION so) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^T
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (&X!=&B) X.copy(B);
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,L->size(0));
//  A X = B ==> L L^T S^{-1} X = S B
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        double *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    int info;
    F77NAME(dpotrs)('L',m,n,L->addr(),m,X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        double *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,L->size(0));
//  X A = B ==> X S^{-1} L L^T = B S
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(dscal)(m,(*s)[j],X.addr(0,j),1);
    }
    F77NAME(dtrsm)('R','L','T','N',m,n,double_one_,L->addr(),n,
      X.addr(),m);
    F77NAME(dtrsm)('R','L','N','N',m,n,double_one_,L->addr(),n,
      X.addr(),m);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(dscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> double CholeskyFactorization<double,double,
SymmetricPositiveMatrix<double,double> >::reciprocalConditionNumber() {
  CHECK_POINTER(L);
  int n=L->size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  F77NAME(dpocon)('L',n,L->addr(),n,anormo,rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveMatrix<double,double> >::improve(const Vector<double,
double> &b,Vector<double,double> &x,double &berr,double &ferr) {
  CHECK_POINTER(L);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(n+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> residual(n);
  Vector<double,double> work(n);
  Vector<double,double> v(n);
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(dsymv)('L',n,double_mone_,A_original->addr(),n,
      x.addr(),1,double_one_,residual.addr(),1);
    work.copy(b);
    F77_NAME(dla_syamv)(F77NAME(ilauplo)('L'),n,double_one_,
      A_original->addr(),n,x.addr(),1,double_one_,work.addr(),1);

    berr=double_zero_;
    if (s!=0) {
      for (int i=0;i<n;i++) {
        berr=(work[i]*(*s)[i]>safe2 ?
          max(berr,abs(residual[i])/work[i]) :
          max(berr,(abs(residual[i])*(*s)[i]+safe1)
               /(work[i]*(*s)[i]+safe1)));
      }
    } else {
      for (int i=0;i<n;i++) {
        berr=(work[i]>safe2 ? max(berr,abs(residual[i])/work[i]) :
          max(berr,(abs(residual[i])+safe1)/(work[i]+safe1)));
      }
    }
    if (!(berr>eps && 2.*berr<=lstres && count<=itmax)) break;
    solve(residual,residual);
    F77NAME(daxpy)(n,double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  double *residuali=residual.addr();
  double *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(dlacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(idamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveMatrix<double,double> >::improve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &berr,Vector<double,double> &ferr,
Factorization::SIDE_OPTION so) {
//compare to dposvx.f
  CHECK_POINTER(L);
  int k=L->size(0),m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,k);
    CHECK_SAME(n,berr.size());
    CHECK_SAME(n,ferr.size());
  } else {
    CHECK_SAME(n,k);
    CHECK_SAME(m,berr.size());
    CHECK_SAME(m,ferr.size());
  }
  double nz=static_cast<double>(k+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> x(k);
  Vector<double,double> rhs(k);
  Vector<double,double> residual(k);
  Vector<double,double> work(k);
  Vector<double,double> v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(dcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(dcopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(dsymv)('L',k,double_mone_,A_original->addr(),k,
        x.addr(),1,double_one_,residual.addr(),1);
      work.copy(rhs);
      F77_NAME(dla_syamv)(F77NAME(ilauplo)('L'),k,double_one_,
        A_original->addr(),k,x.addr(),1,double_one_,work.addr(),1);
      double &berrj=berr[j];
      berrj=double_zero_;
      if (s!=0) {
        for (int i=0;i<k;i++) {
          berrj=(work[i]*(*s)[i]>safe2 ?
            max(berrj,abs(residual[i])/work[i]) :
            max(berrj,(abs(residual[i])*(*s)[i]+safe1)
                 /(work[i]*(*s)[i]+safe1)));
        }
      } else {
        for (int i=0;i<k;i++) {
          berrj=(work[i]>safe2 ? max(berrj,abs(residual[i])/work[i]) :
            max(berrj,(abs(residual[i])+safe1)/(work[i]+safe1)));
        }
      }
      if (!(berrj>eps && 2.*berrj<=lstres && count<=itmax)) break;
      solve(residual,residual);
      F77NAME(daxpy)(k,double_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(dcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      F77NAME(dcopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    double *residuali=residual.addr();
    double *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(dlacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(idamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr[j]/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> CholeskyFactorization<double,double,
SymmetricPositiveTridiagonalMatrix<double,double> >::
CholeskyFactorization(
const SymmetricPositiveTridiagonalMatrix<double,double>& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equ_op(Factorization::NO_EQUILIBRATION),equed('N'),
scond(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  anormi=A.normInfinity();
  anormm=A.normMaxEntry();
  anormo=A.normOne();
  rpvgrw=double_one_; // not computed
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveTridiagonalMatrix<double,double> >::solve(
const Vector<double,double> &b,Vector<double,double> &x) {
//constructor did not factor
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<double,double>(n);
  L->copy(*A_original);
  if (&x!=&b) x.copy(b);
  int info;
  F77NAME(dptsv)(n,1,L->diagonalAddr(0),L->lowerDiagonalAddr(0),
    x.addr(),n,info);
  CHECK_TEST(info==0);
  delete L; L=0;
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveTridiagonalMatrix<double,double> >::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Factorization::SIDE_OPTION so) {
//constructor did not factor
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<double,double>(
    A_original->size(0));
  L->copy(*A_original);
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,L->size(0));
    F77NAME(dptsv)(m,n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),
      X.addr(),m,info);
  } else {
    CHECK_SAME(n,L->size(0));
    F77NAME(dpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
    for (int j=0;j<m;j++) {
      for (int i=1;i<n;i++) {
        X(j,i)-=X(j,i-1)*L->lowerDiagonalValue(i-1);
      }
      X(j,n-1)/=L->diagonalValue(n-1);
      for (int i=n-2;i>=0;i--) {
        X(j,i)=X(j,i)/L->diagonalValue(i)
              -X(j,i+1)*L->lowerDiagonalValue(i);
      }
    }
  }
  CHECK_TEST(info==0);
  delete L; L=0;
}

template<> double CholeskyFactorization<double,double,
SymmetricPositiveTridiagonalMatrix<double,double> >::
reciprocalConditionNumber() {
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<double,double>(
    A_original->size(0));
  L->copy(*A_original);
  int n=L->size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,n);
  int info=0;
  F77NAME(dpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  F77NAME(dptcon)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,work,info);
  CHECK_SAME(info,0)
  delete L; L=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveTridiagonalMatrix<double,double> >::improve(
const Vector<double,double> &b,Vector<double,double> &x,
double &berr,double &ferr) {
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<double,double>(
    A_original->size(0));
  L->copy(*A_original);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double *work=OPERATOR_NEW_BRACKET(double,2*n);
  int info;
  double rcond;
  F77NAME(dpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  CHECK_TEST(info==0);
  F77NAME(dptcon)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,work,info);
  CHECK_TEST(info==0);
  F77NAME(dptrfs)(n,1,A_original->diagonalAddr(0),
    A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
    L->lowerDiagonalAddr(0),b.addr(),n,x.addr(),n,&ferr,&berr,work,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete L; L=0;
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveTridiagonalMatrix<double,double> >::improve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &berr,Vector<double,double> &ferr,
Factorization::SIDE_OPTION so) {
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<double,double>(
    A_original->size(0));
  L->copy(*A_original);
  int k=L->size(0),m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,k);
    CHECK_SAME(n,berr.size());
    CHECK_SAME(n,ferr.size());
  } else {
    CHECK_SAME(n,k);
    CHECK_SAME(m,berr.size());
    CHECK_SAME(m,ferr.size());
  }
  double *work=OPERATOR_NEW_BRACKET(double,2*k);
  int info;
  double rcond;
  F77NAME(dpttrf)(k,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  CHECK_TEST(info==0);
  F77NAME(dptcon)(k,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,work,info);
  CHECK_TEST(info==0);
  if (so==Factorization::LEFT_SIDE) {
    F77NAME(dptrfs)(k,n,A_original->diagonalAddr(0),
      A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
      L->lowerDiagonalAddr(0),B.addr(),m,X.addr(),m,ferr.addr(),
      berr.addr(),work,info);
  } else {
    F77NAME(dptrfsr)(k,m,A_original->diagonalAddr(0),
      A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
      L->lowerDiagonalAddr(0),B.addr(),m,X.addr(),m,ferr.addr(),
      berr.addr(),work,info);
  }
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete L; L=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> CholeskyFactorization<double,double,
SymmetricPositiveBandMatrix<double,double> >::CholeskyFactorization(
const SymmetricPositiveBandMatrix<double,double>& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equ_op(eo),equed('N'),scond(numeric_limits<double>::infinity()) {
  int n=A.size(0),nsub=A.subDiags();
  L=OPERATOR_NEW
    SymmetricPositiveBandMatrix<double,double>(n,nsub,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(n-j,nsub+1),A.addr(j,j),1,L->addr(j,j),1);
  }

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    F77NAME(dpbequ)('L',n,L->subDiags(),L->addr(),L->bands(),s->addr(),
      scond,amax,info);
    equed='N';
    F77NAME(dlaqsb)('L',n,nsub,L->addr(),L->bands(),s->addr(),scond,
      amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=L->normInfinity();
  anormm=L->normMaxEntry();
  anormo=L->normOne();

  int info=0;
  F77NAME(dpbtrf)('L',n,nsub,L->addr(),L->bands(),info);
  double *work=0;
  if (info>0) { // see dgbsvx
    rpvgrw=F77NAME(dlansb)('M','L',info,min(info-1,L->subDiags()),
      L->addr(),L->bands(),work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
    cerr << "zero pivot in CholeskyFactorization::CholeskyFactorization"
       << "\n zero pivot number,reciprocal pivot growth = " << info
       << " " << rpvgrw << endl;
    if (L) delete L; L=0;
    if (s) delete s; s=0;
    A_original=0;
  } else {
    rpvgrw=F77NAME(dlansb)('M','L',n,L->subDiags(),L->addr(),
      L->bands(),work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
  }
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveBandMatrix<double,double> >::solve(
const Vector<double,double> &b,Vector<double,double> &x) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^T
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  if (&x!=&b) x.copy(b);
  int info;
//A x = b ==> L l^T S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    double *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  F77NAME(dpbtrs)('L',n,L->subDiags(),1,L->addr(),L->bands(),
    x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    double *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveBandMatrix<double,double> >::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Factorization::SIDE_OPTION so) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^T
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,L->size(0));
//  A X = B ==> L L^T S^{-1} X = S B
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        double *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    F77NAME(dpbtrs)('L',m,L->subDiags(),n,L->addr(),L->bands(),
      X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        double *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,L->size(0));
//  X A = B ==> X R^{-1} Q L U = B C
//          ==> U^T L^T Q^T R^{-1} X^T = C^T B^T
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(dscal)(m,(*s)[j],X.addr(0,j),1);
    }
    Vector<double,double> *x=OPERATOR_NEW Vector<double,double>(n);
    for (int i=0;i<m;i++) {
      F77NAME(dcopy)(n,X.addr(i,0),m,x->addr(),1);
      F77NAME(dpbtrs)('L',n,L->subDiags(),1,L->addr(),L->bands(),
        x->addr(),n,info);
      F77NAME(dcopy)(n,x->addr(),1,X.addr(i,0),m);
    }
    delete x; x=0;
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(dscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> double CholeskyFactorization<double,double,
SymmetricPositiveBandMatrix<double,double> >::reciprocalConditionNumber()
{
  CHECK_POINTER(L);
  int n=L->size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  F77NAME(dpbcon)('L',n,L->subDiags(),L->addr(),L->bands(),anormo,
    rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveBandMatrix<double,double> >::improve(
const Vector<double,double> &b,Vector<double,double> &x,double &berr,
double &ferr) {
  CHECK_POINTER(L);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(A_original->bands()+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> residual(n);
  Vector<double,double> work(n);
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(dsbmv)('L',n,A_original->subDiags(),double_mone_,
      A_original->addr(),A_original->bands(),x.addr(),1,
      double_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(dsbamv)('L',n,A_original->subDiags(),double_one_,
      A_original->addr(),A_original->bands(),x.addr(),1,
      double_one_,work.addr(),1);

    berr=double_zero_;
    if (s!=0) {
      for (int i=0;i<n;i++) {
        berr=(work[i]*(*s)[i]>safe2 ?
          max(berr,abs(residual[i])/work[i]) :
          max(berr,(abs(residual[i])*(*s)[i]+safe1)
               /(work[i]*(*s)[i]+safe1)));
      }
    } else {
      for (int i=0;i<n;i++) {
        berr=(work[i]>safe2 ? max(berr,abs(residual[i])/work[i]) :
          max(berr,(abs(residual[i])+safe1)/(work[i]+safe1)));
      }
    }
    if (!(berr>eps && 2.*berr<=lstres && count<=itmax)) break;
    solve(residual,residual);
    F77NAME(daxpy)(n,double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }
  double *residuali=residual.addr();
  double *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  Vector<double,double> v(n);
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(dlacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(idamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void CholeskyFactorization<double,double,
SymmetricPositiveBandMatrix<double,double> >::improve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &berr,Vector<double,double> &ferr,
Factorization::SIDE_OPTION so) {
  CHECK_POINTER(L);
  int k=L->size(0),m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,k);
    CHECK_SAME(n,berr.size());
    CHECK_SAME(n,ferr.size());
  } else {
    CHECK_SAME(n,k);
    CHECK_SAME(m,berr.size());
    CHECK_SAME(m,ferr.size());
  }
  double nz=static_cast<double>(A_original->bands()+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> x(k);
  Vector<double,double> rhs(k);
  Vector<double,double> residual(k);
  Vector<double,double> work(k);
  Vector<double,double> v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(dcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(dcopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(dsbmv)('L',k,A_original->subDiags(),double_mone_,
        A_original->addr(),A_original->bands(),x.addr(),1,double_one_,
        residual.addr(),1);
      for (int i=0;i<n;i++) work[i]=abs(rhs[i]);
      F77NAME(dsbamv)('L',k,A_original->subDiags(),double_one_,
        A_original->addr(),A_original->bands(),x.addr(),1,double_one_,
        work.addr(),1);
      double &berrj=berr[j];
      berrj=double_zero_;
      if (s!=0) {
        for (int i=0;i<k;i++) {
          berrj=(work[i]*(*s)[i]>safe2 ?
            max(berrj,abs(residual[i])/work[i]) :
            max(berrj,(abs(residual[i])*(*s)[i]+safe1)
                 /(work[i]*(*s)[i]+safe1)));
        }
      } else {
        for (int i=0;i<k;i++) {
          berrj=(work[i]>safe2 ? max(berrj,abs(residual[i])/work[i]) :
            max(berrj,(abs(residual[i])+safe1)/(work[i]+safe1)));
        }
      }
      if (!(berrj>eps && 2.*berrj<=lstres && count<=itmax)) break;
      solve(residual,residual);
      F77NAME(daxpy)(k,double_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(dcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      F77NAME(dcopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    double *residuali=residual.addr();
    double *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(dlacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],
        kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(idamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template class CholeskyFactorization<double,double,
  SymmetricPositiveMatrix<double,double> >;
template void testCholeskyFactorization(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "MDMtFactorization.C"
template<> MDMtFactorization<double,double>::MDMtFactorization(
const SymmetricMatrix<double,double>& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),MD(0),ipiv(0),s(0),
equed('N'),scond(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  MD=OPERATOR_NEW SymmetricMatrix<double,double>(n);
  MD->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    double *work=OPERATOR_NEW_BRACKET(double,3*n);
    F77NAME(dsyequb)('L',n,MD->addr(),n,s->addr(),scond,amax,work,info);
    equed='N';
    F77NAME(dlaqsy)('L',n,MD->addr(),n,s->addr(),scond,amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
    delete [] work; work=0;
  }
  anormi=MD->normInfinity();
  anormm=MD->normMaxEntry();
  anormo=MD->normOne();

  int info=0;
  ipiv=OPERATOR_NEW_BRACKET(int,n);
  double workq=double_zero_;
  int lwork=-1;
  F77NAME(dsytrf)('L',n,MD->addr(),n,ipiv,&workq,lwork,info);
  lwork=max(2*n,static_cast<int>(workq));
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dsytrf)('L',n,MD->addr(),n,ipiv,work,lwork,info);
  if (info!=0) { // see dsysvxx
    rpvgrw=F77_NAME(dla_syrpvgrw)('L',n,info,A_original->addr(),n,
      MD->addr(),n,ipiv,work);
    delete MD; MD=0;
    if (s!=0) delete s; s=0;
    if (ipiv!=0) delete ipiv; ipiv=0;
  } else {
    rpvgrw=F77_NAME(dla_syrpvgrw)('L',n,n,A_original->addr(),n,
      MD->addr(),n,ipiv,work);
  }
  delete [] work; work=0;
}

template<> void MDMtFactorization<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x) {
  CHECK_POINTER(MD);
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,MD->size(0));
  if (&x!=&b) x.copy(b);
//A x = b ==> M D M^T S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    double *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  int info;
  F77NAME(dsytrs)('L',n,1,MD->addr(),n,ipiv,x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    double *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void MDMtFactorization<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Factorization::SIDE_OPTION so) {
  CHECK_POINTER(MD);
//constructor factored P^T S A S P = M D M^T
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (&X!=&B) X.copy(B);
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,MD->size(0));
//  A X = B ==> M D M^T S^{-1} X = S B
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        double *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    int info;
    F77NAME(dsytrs)('L',m,n,MD->addr(),m,ipiv,X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        double *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,MD->size(0));
//  X A = B ==> X S^{-1} P M D M^T = B S P
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(dscal)(m,(*s)[j],X.addr(0,j),1);
    }
    Vector<double,double> x(n);
    for (int i=0;i<m;i++) {
      double *Xij=X.addr(i,0); 
      double *xj=x.addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *xj=*Xij;
      int info;
      F77NAME(dsytrs)('L',n,1,MD->addr(),n,ipiv,x.addr(),n,info);
      CHECK_TEST(info==0);
      Xij=X.addr(i,0); 
      xj=x.addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *Xij=*xj;
    }
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(dscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> double
MDMtFactorization<double,double>::reciprocalConditionNumber() {
  CHECK_POINTER(MD);
  int n=MD->size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  F77NAME(dsycon)('L',n,MD->addr(),n,ipiv,anormo,rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void MDMtFactorization<double,double>::improve(
const Vector<double,double> &b,Vector<double,double> &x,double &berr,
double &ferr) {
  CHECK_POINTER(MD);
  int n=MD->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(n+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> residual(n);
  Vector<double,double> work(n);
  Vector<double,double> v(n);
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(dsymv)('L',n,double_mone_,A_original->addr(),n,
      x.addr(),1,double_one_,residual.addr(),1);
    work.copy(b);
    F77_NAME(dla_syamv)(F77NAME(ilauplo)('L'),n,double_one_,
      A_original->addr(),n,x.addr(),1,double_one_,work.addr(),1);

    berr=double_zero_;
    if (s!=0) {
      for (int i=0;i<n;i++) {
        berr=(work[i]*(*s)[i]>safe2 ?
          max(berr,abs(residual[i])/work[i]) :
          max(berr,(abs(residual[i])*(*s)[i]+safe1)
               /(work[i]*(*s)[i]+safe1)));
      }
    } else {
      for (int i=0;i<n;i++) {
        berr=(work[i]>safe2 ? max(berr,abs(residual[i])/work[i]) :
          max(berr,(abs(residual[i])+safe1)/(work[i]+safe1)));
      }
    }
    if (!(berr>eps && 2.*berr<=lstres && count<=itmax)) break;
    solve(residual,residual);
    F77NAME(daxpy)(n,double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  double *residuali=residual.addr();
  double *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(dlacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(idamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void MDMtFactorization<double,double>::improve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &berr,Vector<double,double> &ferr,
Factorization::SIDE_OPTION so) {
//compare to dsysvx.f
  CHECK_POINTER(MD);
  int k=MD->size(0),m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,k);
    CHECK_SAME(n,berr.size());
    CHECK_SAME(n,ferr.size());
  } else {
    CHECK_SAME(n,k);
    CHECK_SAME(m,berr.size());
    CHECK_SAME(m,ferr.size());
  }
  double nz=static_cast<double>(k+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,double> x(k);
  Vector<double,double> rhs(k);
  Vector<double,double> residual(k);
  Vector<double,double> work(k);
  Vector<double,double> v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(dcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(dcopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(dcopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(dsymv)('L',k,double_mone_,A_original->addr(),k,
        x.addr(),1,double_one_,residual.addr(),1);
      work.copy(rhs);
      F77_NAME(dla_syamv)(F77NAME(ilauplo)('L'),k,double_one_,
        A_original->addr(),k,x.addr(),1,double_one_,work.addr(),1);
      double &berrj=berr[j];
      berrj=double_zero_;
      if (s!=0) {
        for (int i=0;i<k;i++) {
          berrj=(work[i]*(*s)[i]>safe2 ?
            max(berrj,abs(residual[i])/work[i]) :
            max(berrj,(abs(residual[i])*(*s)[i]+safe1)
                 /(work[i]*(*s)[i]+safe1)));
        }
      } else {
        for (int i=0;i<k;i++) {
          berrj=(work[i]>safe2 ? max(berrj,abs(residual[i])/work[i]) :
            max(berrj,(abs(residual[i])+safe1)/(work[i]+safe1)));
        }
      }
      if (!(berrj>eps && 2.*berrj<=lstres && count<=itmax)) break;
      solve(residual,residual);
      F77NAME(daxpy)(k,double_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(dcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      F77NAME(dcopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    double *residuali=residual.addr();
    double *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(dlacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(idamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr[j]/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template class MDMtFactorization<double,double>;
template void testMDMtFactorization(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "HouseholderQRFactorization.C"
template<>
HouseholderQRFactorization<double,double>::HouseholderQRFactorization(
const Matrix<double,double> &A,Factorization::PIVOT_OPTION po) :
QR(0),tau(0),jpvt(0),piv_op(po),A_original(&A),iascl(0),
ascl(double_one_),anrm(numeric_limits<double>::infinity()) {
//TRACER_CALL(t,"HouseholderQRFactorization::HouseholderQRFactorization");
  int m=A.size(0),n=A.size(1);
  QR=OPERATOR_NEW Matrix<double,double>(m,n);
  QR->copy(A);
  double w=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);
  anrm=F77NAME(dlange)('M',m,n,A.addr(),m,&w);
  if (anrm>double_zero_ && anrm<smlnum) {
    ascl=smlnum;
    iascl=1;
  } else if (anrm>bignum) {
    ascl=bignum;
    iascl=2;
  } else if (anrm==double_zero_) {
    delete QR; QR=0;
    return;
  }
  int info=0;
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,anrm,ascl,m,n,QR->addr(),m,info);
    CHECK_SAME(info,0);
  }

  if (piv_op==Factorization::NO_PIVOTING) {
    int lwork=-1;
    if (m>=n) {
      tau=OPERATOR_NEW Vector<double,double>(n);
      F77NAME(dgeqrf)(m,n,QR->addr(),m,tau->addr(),&w,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dgeqrf)(m,n,QR->addr(),m,tau->addr(),work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work;
    } else {
      tau=OPERATOR_NEW Vector<double,double>(m);
      F77NAME(dgelqf)(m,n,QR->addr(),m,tau->addr(),&w,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dgelqf)(m,n,QR->addr(),m,tau->addr(),work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work;
    }
  } else {
    int mn=min(m,n);
    int nb=F77NAME(ilaenv)(1,"DGEQRF"," ",m,n,-1,1);
    int lwkmin=mn+max(2*mn,n+1);
    int lwork=max(lwkmin,mn+2*n+nb*(n+1));
    double *work=OPERATOR_NEW_BRACKET(double,lwork);
    tau=OPERATOR_NEW Vector<double,double>(mn);
    jpvt=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(dgeqp3)(m,n,QR->addr(),m,jpvt,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
  }
}

template<> OrthogonalMatrix<double,double>* 
HouseholderQRFactorization<double,double>::orthogonalPart() const {
//TRACER_CALL(t,"HouseholderQRFactorization::orthogonalPart");
  int m=QR->size(0),n=QR->size(1);
  if (m>=n) {
    OrthogonalMatrix<double,double> *Q=
      OPERATOR_NEW OrthogonalMatrix<double,double>(m,m);
    Q->copyFrom('A',m,n,*QR);

    double w=numeric_limits<double>::infinity();
    int lwork=-1;
    int info;
    F77NAME(dorgqr)(m,m,n,Q->addr(),m,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    double *work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dorgqr)(m,m,n,Q->addr(),m,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    return Q;
  } else {
    OrthogonalMatrix<double,double> *Q=
      OPERATOR_NEW OrthogonalMatrix<double,double>(n,n);
    Q->copyFrom('A',m,n,*QR);

    double w=numeric_limits<double>::infinity();
    int lwork=-1;
    int info;
    F77NAME(dorglq)(n,n,m,Q->addr(),n,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    double *work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dorglq)(n,n,m,Q->addr(),n,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    return Q;
  }
}

template<> double HouseholderQRFactorization<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::solve");
  int m=QR->size(0),n=QR->size(1);
  double w=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(dlange)('M',brow,1,b.addr(),b.size(),&w);
  int ibscl=0;
  double bscl=double_one_;
  if (bnrm>double_zero_ && bnrm<smlnum) {
    ibscl=1;
    bscl=smlnum;
  } else if (bnrm > bignum) {
    ibscl=2;
    bscl=bignum;
  }

  int info;
  int lwork=-1;
  int scllen=0;
  double residual_norm=numeric_limits<double>::infinity();
  if (piv_op==Factorization::NO_PIVOTING) {
    if (m>=n) {
      if (tr==Factorization::NO_TRANSPOSE) { // min || b - A x ||
//      if A=Q[R], r=b-Ax and 0=A'r then
//            [0]
//        [v]=Q'r=Q'b-[R]x=[y]-[Rx] and 0=[R',0]Q'r=R'v ==> v=0, w=z
//        [w]         [0]  [z] [0]
//      so compute [y]=Q'b, solve Rx=y
//                 [z]
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        Vector<double,double> xtmp(m);
        xtmp.copy(b);
        if (ibscl!=0) {
          F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      Q' b = [ y \\ z ]:
        F77NAME(dormqr)('L','T',m,1,n,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,&w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dormqr)('L','T',m,1,n,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,work,lwork,info);
        CHECK_SAME(info,0)

//      solve R x = y:
        F77NAME(dtrsv)('U','N','N',n,QR->addr(),m,xtmp.addr(),1);
        delete [] work; work=0;
        x.copyFrom(n,xtmp);
        residual_norm=F77NAME(dnrm2)(m-n,xtmp.addr(n),1);
        scllen=n;
      } else { // min || x || s.t. A' x = b
//      if A=Q[R], b=A'x and x=Ax then
//            [0]
//        b=[R',0]Q'x=[R',0][v]=R'v and [v]=Q'x=[R]s ==> w=0
//                          [w]         [w]     [0]
//      so solve R'v=b, compute x=Q[v]
//                                 [0]
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())

        x=double_zero_;
        x.copyFrom(m,b);
        if (ibscl!=0) {
          F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,x.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      solve R' v = b
        F77NAME(dtrsv)('U','T','N',n,QR->addr(),m,x.addr(),1);
//      Q [v\\0]
        F77NAME(dormqr)('L','N',m,1,n,QR->addr(),m,tau->addr(),
                        x.addr(),m,&w,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dormqr)('L','N',m,1,n,QR->addr(),m,tau->addr(),
                        x.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=double_zero_;
        scllen=m;
      }
    } else {
      if (tr==Factorization::NO_TRANSPOSE) {// min || x || s.t. A x = b
//      if A=[L,0] Q, b=Ax and x=A'x then
//        b=[L,0]Qx=[L,0][v]=Lv and [v]=Qx=[L']s ==> w=0
//                       [w]        [w]    [0 ]
//      so solve Lv=b, compute x=Q'[v]
//                                 [0]
        CHECK_SAME(m,b.size());
        CHECK_SAME(n,x.size());
        x=double_zero_;
        x.copyFrom(m,b);
        if (ibscl!=0) {
          F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,x.addr(),m,info);
          CHECK_SAME(info,0)
        }
        F77NAME(dtrsv)('L','N','N',m,QR->addr(),m,x.addr(),1);
        F77NAME(dormlq)('L','T',n,1,m,QR->addr(),m,tau->addr(),
                        x.addr(),n,&w,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dormlq)('L','T',n,1,m,QR->addr(),m,tau->addr(),
                        x.addr(),n,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=double_zero_;
        scllen=n;
      } else { // min || b - A' x ||
//      if A=[L,0]Q, r=b-A'x and 0=Ar then
//        [v]=Qr=Qb-[L']x=[y]-[L'x] and 0=[L,0]Qr=Lv ==> v=0, w=z
//        [w]       [0 ]  [z] [ 0 ]
//      so compute [y]=Qb, solve L'x=y
//                 [z]
        CHECK_SAME(n,b.size());
        CHECK_SAME(m,x.size());
        Vector<double,double> xtmp(n);
        xtmp.copy(b);
        if (ibscl!=0) {
          F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
          CHECK_SAME(info,0)
        }
        F77NAME(dormlq)('L','N',m,1,m,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,&w,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dormlq)('L','N',m,1,m,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        F77NAME(dtrsv)('L','T','N',m,QR->addr(),m,xtmp.addr(),1);
        x.copyFrom(m,xtmp);
        residual_norm=F77NAME(dnrm2)(n-m,xtmp.addr(m),1);
        scllen=m;
      }
    }
  } else {
    int mn=min(m,n);
    if (tr==Factorization::NO_TRANSPOSE) {// min || x || s.t. A x = b
      CHECK_SAME(m,b.size())
      CHECK_SAME(n,x.size())
      Vector<double,double> xtmp(m);
      xtmp.copy(b);
      if (ibscl!=0) {
        F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    Q' b = [ y \\ z ]:
      F77NAME(dormqr)('L','T',m,1,mn,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,&w,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dormqr)('L','T',m,1,mn,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
//    solve R x' = y:
      F77NAME(dtrsv)('U','N','N',mn,QR->addr(),m,xtmp.addr(),1);
//    permute
      for (int i=0;i<n;i++) x[jpvt[i]-1]=xtmp[i];
      residual_norm=F77NAME(dnrm2)(m-mn,xtmp.addr(mn),1);
    } else {
      CHECK_SAME(n,b.size())
      CHECK_SAME(m,x.size())
      Vector<double,double> xtmp(n);
//    P'b:
      for (int i=0;i<n;i++) xtmp[i]=b[jpvt[i]-1];
      if (ibscl!=0) {
        F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    solve R' v = P' b
      F77NAME(dtrsv)('U','T','N',mn,QR->addr(),m,xtmp.addr(),1);
      x=double_zero_;
      x.copyFrom(mn,xtmp);
//    Q [v\\0]
      double w=numeric_limits<double>::infinity();
      lwork=-1;
      F77NAME(dormqr)('L','N',m,1,mn,QR->addr(),m,tau->addr(),
                      x.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)
      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dormqr)('L','N',m,1,mn,QR->addr(),m,tau->addr(),
                      x.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      if (mn<n) {
        residual_norm=F77NAME(dnrm2)(n-mn,xtmp.addr(mn),1);
      } else residual_norm=double_zero_;
    }
  }
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(dlascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
  return residual_norm;
}

template<> void HouseholderQRFactorization<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &residual_norm,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::solve");
  int m=QR->size(0),n=QR->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  CHECK_SAME(k,residual_norm.size())
  double w=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(dlange)('M',brow,k,B.addr(),B.size(0),&w);
  int ibscl=0;
  double bscl=double_one_;
  if (bnrm>double_zero_ && bnrm<smlnum) {
    ibscl=1;
    bscl=smlnum;
  } else if (bnrm > bignum) {
    ibscl=2;
    bscl=bignum;
  }

  int info;
  int lwork=-1;
  int scllen=0;
  if (piv_op==Factorization::NO_PIVOTING) {
    if (m>=n) {
      if (tr==Factorization::NO_TRANSPOSE) { // min || B - A x ||
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        Matrix<double,double> Xtmp(m,k);
        Xtmp.copy(B);
        if (ibscl!=0) {
          F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      Q' B = [ Y \\ Z ]:
        F77NAME(dormqr)('L','T',m,k,n,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,&w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dormqr)('L','T',m,k,n,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,work,lwork,info);
        CHECK_SAME(info,0)

//      solve R X = Y:
        F77NAME(dtrsm)('L','U','N','N',n,k,double_one_,QR->addr(),m,
          Xtmp.addr(),m);
        delete [] work; work=0;
        X.copyFrom('A',n,k,Xtmp);
        for (int j=0;j<k;j++) {
          residual_norm[j]=F77NAME(dnrm2)(m-n,Xtmp.addr(n,j),1);
        }
        scllen=n;
      } else { // min || x || s.t. A^t x = B
        CHECK_SAME(n,B.size(0))
        CHECK_SAME(m,X.size(0))

        X=double_zero_;
        X.copyFrom('A',m,k,B);
        if (ibscl!=0) {
          F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,k,X.addr(),m,info);
          CHECK_SAME(info,0)
        }
        F77NAME(dtrsm)('L','U','T','N',n,k,double_one_,QR->addr(),m,
          X.addr(),m);
        F77NAME(dormqr)('L','N',m,k,n,QR->addr(),m,tau->addr(),
                        X.addr(),m,&w,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dormqr)('L','N',m,k,n,QR->addr(),m,tau->addr(),
                        X.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=double_zero_;
        scllen=m;
      }
    } else {
      if (tr==Factorization::NO_TRANSPOSE) { // min || x || s.t. A x = B
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        X=double_zero_;
        X.copyFrom('A',m,k,B);
        if (ibscl!=0) {
          F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,k,X.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      solve L X = B:
        F77NAME(dtrsm)('L','L','N','N',m,k,double_one_,QR->addr(),m,
          X.addr(),n);
//      Q' V
        F77NAME(dormlq)('L','T',n,k,m,QR->addr(),m,tau->addr(),
                        X.addr(),n,&w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dormlq)('L','T',n,k,m,QR->addr(),m,tau->addr(),
                        X.addr(),n,work,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=double_zero_;
        scllen=n;
      } else { // // min || B - A' X ||
        CHECK_SAME(n,B.size(0))
        CHECK_SAME(m,X.size(0))

        Matrix<double,double> Xtmp(n,k);
        Xtmp.copy(B);
        if (ibscl!=0) {
          F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
          CHECK_SAME(info,0)
        }
        F77NAME(dormlq)('L','N',m,k,m,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,&w,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        double *work=OPERATOR_NEW_BRACKET(double,lwork);
        F77NAME(dormlq)('L','N',m,k,m,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;

        F77NAME(dtrsm)('L','L','T','N',m,k,double_one_,QR->addr(),m,
          Xtmp.addr(),n);
        X=double_zero_;
        X.copyFrom('A',n,k,Xtmp);
        residual_norm=double_zero_;
        for (int j=0;j<k;j++) {
          residual_norm[j]=F77NAME(dnrm2)(n-m,Xtmp.addr(m,j),1);
        }
        scllen=m;
      }
    }
  } else {
//  TRACER_CALL(t,"HouseholderQRFactorization::solve");
    int mn=min(m,n);
    if (tr==Factorization::NO_TRANSPOSE) {
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<double,double> Xtmp(m,k);
      Xtmp.copy(B);
      if (ibscl!=0) {
        F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    Q' B = [ Y \\ Z ]:
      F77NAME(dormqr)('L','T',m,k,mn,QR->addr(),m,tau->addr(),
                      Xtmp.addr(),m,&w,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dormqr)('L','T',m,k,mn,QR->addr(),m,tau->addr(),
                      Xtmp.addr(),m,work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
//    solve R V = Y:
      F77NAME(dtrsv)('U','N','N',mn,QR->addr(),m,Xtmp.addr(),1);
//    P V:
      for (int j=0;j<k;j++) {
        for (int i=0;i<n;i++) X(jpvt[i]-1,j)=Xtmp(i,j);
        residual_norm[j]=F77NAME(dnrm2)(m-mn,Xtmp.addr(mn,j),1);
      }
    } else {
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<double,double> Xtmp(n,k);
//    P' B:
      for (int j=0;j<k;j++) {
        for (int i=0;i<n;i++) Xtmp(i,j)=B(jpvt[i]-1,j);
      }
      if (ibscl!=0) {
        F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    solve R' V = P' B
      F77NAME(dtrsv)('U','T','N',mn,QR->addr(),m,Xtmp.addr(),1);
      X=double_zero_;
      X.copyFrom('A',mn,k,Xtmp);
//    Q [V\\0]
      w=numeric_limits<double>::infinity();
      lwork=-1;
      F77NAME(dormqr)('L','N',m,k,mn,QR->addr(),m,tau->addr(),
                      X.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)
      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dormqr)('L','N',m,k,mn,QR->addr(),m,tau->addr(),
                      X.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      if (mn<n) {
        for (int j=0;j<k;j++) {
          residual_norm[j]=F77NAME(dnrm2)(n-mn,Xtmp.addr(mn,j),1);
        }
      } else residual_norm=double_zero_;
    }
  }
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(dlascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
}

template<> void HouseholderQRFactorization<double,double>::improve(
const Vector<double,double> &b,Vector<double,double> &x,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::improve");
  CHECK_TEST(piv_op==Factorization::NO_PIVOTING)
  int m=QR->size(0),n=QR->size(1);
  double w=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(dlange)('M',brow,1,b.addr(),b.size(),&w);
  int ibscl=0;
  double bscl=double_one_;
  if (bnrm>double_zero_ && bnrm<smlnum) {
    ibscl=1;
    bscl=smlnum;
  } else if (bnrm > bignum) {
    ibscl=2;
    bscl=bignum;
  }

  int info;
  int lwork=-1;
  int scllen=0;
  if (ibscl!=0) {
    F77NAME(dlascl)('G',0,0,bnrm,bscl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,ascl,anrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (m>=n) {
    if (tr==Factorization::NO_TRANSPOSE) { // min || b - A x ||
//    if A=Q[R], db=b-r-Ax=dr+Adx and dc=-A'r=A'dr then
//          [0]
//      [dy]=Q'db=Q'dr+[R]dx=[dv]+[Rdx] ==> dw=dz and Rdx=dy-dv
//      [dz]           [0]   [dw] [0  ]
//      dc=[R',0]Q'dr=[R',0][dv]=R'dv
//                          [dz]
//    so compute [dy\\dz]=Q'(b-r-Ax), solve R'(-dv)=A'r, Rdx=dy+(-dv),
//      dr=Q[dv\\dz]
      CHECK_SAME(m,b.size())
      CHECK_SAME(n,x.size())
      Vector<double,double> r(m);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    r=b-Ax
      F77NAME(dgemv)('N',m,n,double_mone_,A_original->addr(),m,x.addr(),1,
        double_one_,r.addr(),1);
      Vector<double,double> c(n,double_zero_);
//    -dc=A'r
      F77NAME(dgemv)('T',m,n,double_one_,A_original->addr(),m,r.addr(),1,
        double_zero_,c.addr(),1);
//    solve R'(-dv) =-dc:
      F77NAME(dtrsv)('U','T','N',n,QR->addr(),m,c.addr(),1);
//    solve Rdx=-v:
      F77NAME(dtrsv)('U','N','N',n,QR->addr(),m,c.addr(),1);
      x+=c;
      scllen=n;
    } else { // min || x || s.t. A^t x = b
//    if A=Q[R], da=-x+As=dx-Ads and db=b-A'x=A'dx then
//          [0]
//      db=[R',0]Q'dx=[R',0][dv]=R'dv
//                          [dw]
//      [df]=Q'da=Q'dx-[R']ds=[dv]-Rds ==> dw=dg, Rds=dv-df
//      [dg]           [0 ]   [dw]
//    so compute [df\\dg]=Q'(-x+As), solve R'dv=b-A'x,
//      compute dx=Q[dv\\dg], solve Rds=dv-df
      CHECK_SAME(n,b.size())
      CHECK_SAME(m,x.size())
      Vector<double,double> r(n);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    r=b-A'x
      F77NAME(dgemv)('T',m,n,double_mone_,A_original->addr(),m,x.addr(),1,
        double_one_,r.addr(),1);
//    solve R'dv=r:
      Vector<double,double> c(m,double_zero_);
      c.copyFrom(n,r);
      F77NAME(dtrsv)('U','T','N',n,QR->addr(),m,c.addr(),1);
//    dx=Q[dv\\0]
      F77NAME(dormqr)('L','T',m,1,n,QR->addr(),m,tau->addr(),
                      c.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dormqr)('L','T',m,1,n,QR->addr(),m,tau->addr(),
                      c.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      x+=c;
      scllen=m;
    }
  } else {
    if (tr==Factorization::NO_TRANSPOSE) { // min || x || s.t. A x = b
//    if A=[L,0]Q, da=-x+A's=dx-A'ds and db=b-Ax=Adx then
//      db=[L,0]Qdx=[L,0][dv]=Ldv
//                       [dw]
//      [df]=Qda=Qdx-[L]ds=[dv]-L'ds ==> dw=dg, L'ds=dv-df
//      [dg]         [0]   [dw]
//    so compute [df\\dg]=Q(-x+A's), solve Ldv=b-Ax,
//      compute dx=Q'[dv\\dg], solve L'ds=dv-df
      CHECK_SAME(m,b.size())
      CHECK_SAME(n,x.size())
      Vector<double,double> v(n,double_zero_);
      v.copyFrom(m,b);
      if (ibscl!=0) {
        F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,v.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    r=b-Ax
      F77NAME(dgemv)('N',m,n,double_mone_,A_original->addr(),m,x.addr(),1,
        double_one_,v.addr(),1);
//    solve Ldv=r:
      F77NAME(dtrsv)('L','N','N',m,QR->addr(),m,v.addr(),1);
//    dx=A'[dv\\dg]:
      F77NAME(dormlq)('L','T',n,1,m,QR->addr(),m,tau->addr(),
                      v.addr(),n,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dormlq)('L','T',n,1,m,QR->addr(),m,tau->addr(),
                      v.addr(),n,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      x+=v;
      scllen=n;
    } else { // min || b - A' x ||
//    if A=[L,0]Q, db=b-r-A'x=dr+A'dx and dc=-Ar=Adr then
//      [dy]=Qdb=Qdr+[L']dx=[dv]+[L'dx] ==> dw=dz and L'dx=dy-dv
//      [dz]         [0 ]   [dw] [0   ]
//      dc=[L,0]Qdr=[L,0][dv]=Ldv
//                       [dz]
//    so compute [dy\\dz]=Q(b-r-A'x), solve L(-dv)=Ar, L'dx=dy+(-dv),
//      dr=Q'[dv\\dz]
      CHECK_SAME(n,b.size())
      CHECK_SAME(m,x.size())
      Vector<double,double> r(n);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    r=b-A'x
      F77NAME(dgemv)('T',m,n,double_mone_,A_original->addr(),m,x.addr(),1,
        double_one_,r.addr(),1);
      Vector<double,double> c(m,double_zero_);
//    -dc=Ar
      F77NAME(dgemv)('N',m,n,double_one_,A_original->addr(),m,r.addr(),1,
        double_zero_,c.addr(),1);
//    solve L(-dv)=-dc:
      F77NAME(dtrsv)('L','N','N',m,QR->addr(),m,c.addr(),1);
//    solve L'dx=-dv:
      F77NAME(dtrsv)('L','T','N',m,QR->addr(),m,c.addr(),1);
      x+=c;
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(dlascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
}

template<> void HouseholderQRFactorization<double,double>::improve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::improve");
  CHECK_TEST(piv_op==Factorization::NO_PIVOTING)
  int m=QR->size(0),n=QR->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  double w=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);
  
  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(dlange)('M',brow,k,B.addr(),B.size(0),&w);
  int ibscl=0;
  double bscl=double_one_;
  if (bnrm>double_zero_ && bnrm<smlnum) {
    ibscl=1;
    bscl=smlnum;
  } else if (bnrm > bignum) {
    ibscl=2;
    bscl=bignum;
  }

  int info;
  int lwork=-1;
  int scllen=0;
  if (ibscl!=0) {
    F77NAME(dlascl)('G',0,0,bnrm,bscl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,ascl,anrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (m>=n) {
    if (tr==Factorization::NO_TRANSPOSE) { // min || B - A x ||
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<double,double> R(m,k);
      R.copy(B);
//    R=B-AX
      F77NAME(dgemm)('N','N',m,k,n,double_mone_,A_original->addr(),m,
        X.addr(),n,double_one_,R.addr(),m);
      Matrix<double,double> C(n,k,double_zero_);
//    -dC=A'R
      F77NAME(dgemm)('T','N',n,k,m,double_one_,A_original->addr(),m,
        R.addr(),m,double_zero_,C.addr(),n);
//    solve R'(-dV)=-dC:
      F77NAME(dtrsm)('L','U','T','N',n,k,double_one_,QR->addr(),m,
        C.addr(),n);
//    solve RdX=-V:
      F77NAME(dtrsm)('L','U','N','N',n,k,double_one_,QR->addr(),m,
        C.addr(),n);
      X+=C;
      scllen=n;
    } else { // min || X || s.t. A^t X = B
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<double,double> R(n,k);
      R.copy(B);
//    R=B-A'X
      F77NAME(dgemm)('T','N',n,k,m,double_mone_,A_original->addr(),m,
        X.addr(),m,double_one_,R.addr(),n);
//    solve R'dV=R:
      Matrix<double,double> C(m,k,double_zero_);
      C.copyFrom('A',n,k,R);
      F77NAME(dtrsm)('L','U','T','N',n,k,double_one_,QR->addr(),m,
        C.addr(),n);
//    dX=Q[dV\\0]
      F77NAME(dormqr)('L','T',m,k,n,QR->addr(),m,tau->addr(),
                      C.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dormqr)('L','T',m,k,n,QR->addr(),m,tau->addr(),
                      C.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      X+=C;
      scllen=m;
    }
  } else {
    if (tr==Factorization::NO_TRANSPOSE) { // min || X || s.t. A X = B
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<double,double> V(n,k,double_zero_);
      V.copyFrom('A',m,k,B);
//    R=B-AX
      F77NAME(dgemm)('N','N',m,k,n,double_mone_,A_original->addr(),m,
        X.addr(),n,double_one_,V.addr(),m);
//    solve LdV=R:
      F77NAME(dtrsm)('L','L','N','N',m,k,double_one_,QR->addr(),m,
        V.addr(),n);
//    dX=Q'[dV\\0]:
      F77NAME(dormlq)('L','T',n,k,m,QR->addr(),m,tau->addr(),
                      V.addr(),n,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      double *work=OPERATOR_NEW_BRACKET(double,lwork);
      F77NAME(dormlq)('L','T',n,k,m,QR->addr(),m,tau->addr(),
                      V.addr(),n,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      X+=V;
      scllen=n;
    } else { // min || B - A' X ||
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<double,double> R(n,k);
      R.copy(B);
//    R=B-A'X
      F77NAME(dgemm)('T','N',n,k,m,double_mone_,A_original->addr(),m,
        X.addr(),m,double_one_,R.addr(),n);
      Matrix<double,double> C(m,k,double_zero_);
//    (-dC)=A R
      F77NAME(dgemm)('N','N',m,k,n,double_one_,A_original->addr(),m,
        R.addr(),n,double_zero_,C.addr(),m);
//    L(-dV)=-dC
      F77NAME(dtrsm)('L','L','N','N',m,k,double_one_,QR->addr(),m,
        C.addr(),m);
//    L'dX=-dV
      F77NAME(dtrsm)('L','L','T','N',m,k,double_one_,QR->addr(),m,
        C.addr(),m);
      X+=C;
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(dlascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
}

template class HouseholderQRFactorization<double,double>;
//template void testHouseholderQRFactorization(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "CompleteOrthogonalDecomposition.C"
template<> CompleteOrthogonalDecomposition<double,double>::
CompleteOrthogonalDecomposition(const Matrix<double,double> &A,
double rcond) : the_rank(0),URV(0),utau(0),vtau(0),jpvt(0),iascl(0),
ascl(double_one_),anrm(numeric_limits<double>::infinity()),A_original(&A)
{
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::CompleteOrthogonalDecomposition");
  int m=A.size(0),n=A.size(1);
  URV=OPERATOR_NEW Matrix<double,double>(m,n);
  URV->copy(A);
  double w=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);
  anrm=F77NAME(dlange)('M',m,n,A.addr(),m,&w);
  if (anrm>double_zero_ && anrm<smlnum) {
    ascl=smlnum;
    iascl=1;
  } else if (anrm>bignum) {
    ascl=bignum;
    iascl=2;
  } else if (anrm==double_zero_) {
    delete URV; URV=0;
    return;
  }
  int info=0;
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,anrm,ascl,m,n,URV->addr(),m,info);
    CHECK_SAME(info,0);
  }

  int mn=min(m,n);
  int nb1=F77NAME(ilaenv)(1,"DGEQRF"," ",m,n,-1,1);
  int nb2=F77NAME(ilaenv)(1,"DGERQF"," ",m,n,-1,1);
  int nb=max(nb1,nb2);
  int lwkmin=mn+max(2*mn,n+1);
  int lwork=max(lwkmin,mn+2*n+nb*(n+1));
  double *work=OPERATOR_NEW_BRACKET(double,lwork);

  utau=OPERATOR_NEW Vector<double,double>(mn);
  jpvt=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(dgeqp3)(m,n,URV->addr(),m,jpvt,utau->addr(),work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  double smax=abs((*URV)(0,0));
  if (smax<=double_zero_) return;
  double smin=smax;

  the_rank=1;
  Vector<double,double> xmin(mn);
  xmin[0]=double_one_;
  Vector<double,double> xmax(mn);
  xmax[0]=double_one_;
  while (the_rank<mn) {
    int i=the_rank+1;
    double sminpr=numeric_limits<double>::infinity(),
      s1=numeric_limits<double>::infinity(),
      c1=numeric_limits<double>::infinity();
//                 (job,j,x,sest,w,gamma,sestpr,s,c)
    F77NAME(dlaic1)(2,the_rank,xmin.addr(),smin,URV->addr(0,i-1),
      (*URV)(i-1,i-1),sminpr,s1,c1);
    double smaxpr=numeric_limits<double>::infinity(),
      s2=numeric_limits<double>::infinity(),
      c2=numeric_limits<double>::infinity();
    F77NAME(dlaic1)(1,the_rank,xmax.addr(),smax,URV->addr(0,i-1),
      (*URV)(i-1,i-1),smaxpr,s2,c2);
    if (smaxpr*rcond>sminpr) break;
    for (int i=0;i<the_rank;i++) {
      xmin[i]=s1*xmin[i];
      xmax[i]=s2*xmax[i];
    }
    xmin[the_rank-1]=c1;
    xmax[the_rank-1]=c2;
    smin=sminpr;
    smax=smaxpr;
    the_rank++;
  }
  if (the_rank<n) {
    vtau=OPERATOR_NEW Vector<double,double>(the_rank);
    w=numeric_limits<double>::infinity();
    int lwork=-1;
    F77NAME(dtzrzf)(the_rank,n,URV->addr(),m,vtau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dtzrzf)(the_rank,n,URV->addr(),m,vtau->addr(),work,lwork,
      info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
  }
}

template<> OrthogonalMatrix<double,double>* 
CompleteOrthogonalDecomposition<double,double>::leftOrthogonalPart()
const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::leftOrthogonalPart");
  int m=URV->size(0),n=URV->size(1);
  OrthogonalMatrix<double,double> *Q=
    OPERATOR_NEW OrthogonalMatrix<double,double>(m,m);
  Q->copyFrom('A',m,the_rank,*URV);

  double w=numeric_limits<double>::infinity();
  int lwork=-1;
  int info;
  F77NAME(dorgqr)(m,m,the_rank,Q->addr(),m,utau->addr(),&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dorgqr)(m,m,the_rank,Q->addr(),m,utau->addr(),work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  return Q;
}

template<> OrthogonalMatrix<double,double>* 
CompleteOrthogonalDecomposition<double,double>
::rightOrthogonalPartTransposed() const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::rightOrthogonalPartTransposed");
  int m=URV->size(0),n=URV->size(1);
  OrthogonalMatrix<double,double> *Q=
    OPERATOR_NEW OrthogonalMatrix<double,double>(n,n);

  double w=numeric_limits<double>::infinity();
  int lwork=-1;
  int info;
  F77NAME(dormrz)('L','T',n,n,the_rank,n-the_rank,URV->addr(),m,
    vtau->addr(), Q->addr(),n,&w,lwork,info);
  CHECK_SAME(info,0)
  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dormrz)('L','T',n,n,the_rank,n-the_rank,URV->addr(),m,
    vtau->addr(), Q->addr(),n,work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  return Q;
}

template<> double CompleteOrthogonalDecomposition<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve");
  int m=URV->size(0),n=URV->size(1);
  double w=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(dlange)('M',brow,1,b.addr(),b.size(),&w);
  int ibscl=0;
  double bscl=double_one_;
  if (bnrm>double_zero_ && bnrm<smlnum) {
    ibscl=1;
    bscl=smlnum;
  } else if (bnrm > bignum) {
    ibscl=2;
    bscl=bignum;
  }

  int info;
  int mn=min(m,n);
  int lwork=-1;
  int scllen=0;
  double residual_norm=numeric_limits<double>::infinity();
  if (tr==Factorization::NO_TRANSPOSE) { // min || b - A x || min || x ||
//  if A=U[R 0]V'P', r=b-Ax, 0=A'r and X=A's then
//        [0 0]
//    [f]=U'r=U'b-[R 0]V'P'x=[y]-[R 0][v]=[y-Rv] ==> g=z
//    [g]         [0 0]      [z] [0 0][w] [z]
//  and
//    0=A'r=PV[R' 0]U'r=PV[R' 0][f]=PV[R'f] ==> f=0
//            [0  0]      [0  0][g]   [0  ]
//  and
//    [v]=V'P'x=[R' 0]U's=[R' 0][p]=[R'p] ==> w=0
//    [w]       [0  0]    [0  0][q] [0  ]
//  so compute [y\\z]=U'b, solve Rv=y, compute x=PV[v\\0]
    CHECK_SAME(m,b.size())
    CHECK_SAME(n,x.size())
    Vector<double,double> xtmp(m);
    xtmp.copy(b);
    if (ibscl!=0) {
      F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
      CHECK_SAME(info,0)
    }
//  U' b = [ y \\ z ]:
    F77NAME(dormqr)('L','T',m,1,mn,URV->addr(),m,utau->addr(),
                    xtmp.addr(),m,&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    double *work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dormqr)('L','T',m,1,mn,URV->addr(),m,utau->addr(),
                    xtmp.addr(),m,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R v = y:
    F77NAME(dtrsv)('U','N','N',the_rank,URV->addr(),m,xtmp.addr(),1);
    Vector<double,double> px(n,double_zero_);
    px.copyFrom(the_rank,xtmp);
    double w=numeric_limits<double>::infinity();
    lwork=-1;
//  P'x=V[v\\0]:
    F77NAME(dormrz)('L','T',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(), px.addr(),n,&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dormrz)('L','T',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),px.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  permute
    for (int i=0;i<n;i++) x[jpvt[i]-1]=px[i];
    residual_norm=F77NAME(dnrm2)(m-the_rank,xtmp.addr(the_rank),1);
  } else { // min || x || min || b - A' x ||
//  if A=U[R 0]V'P', r=b-A'x, 0=Ar and x=Ax then
//          [0 0]
//    [f]=V'P'r=V'P'b-[R' 0]U'x=[y]-[R' 0][v]=[y-Rv] ==> g=z, R'v=y
//    [g]             [0  0]    [z] [0  0][w] [z   ]
//  and
//    0=Ar=U[R 0]V'P'r=U[R 0][f]=U[Rf] ==> f=0
//          [0 0]       [0 0][g]  [0 ]
//  and
//    [v]=Ux=[R 0]V'P's=[R 0][p]=[Rp] ==> w=0
//    [w]    [0 0]      [0 0][q] [0 ]
//  so compute [y\\z]=V'P'b, solve R'v=y, compute x=U[v\\0]
    CHECK_SAME(n,b.size())
    CHECK_SAME(m,x.size())
    Vector<double,double> xtmp(n);
//  P'b:
    for (int i=0;i<n;i++) xtmp[i]=b[jpvt[i]-1];
    if (ibscl!=0) {
      F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
      CHECK_SAME(info,0)
    }
    double w=numeric_limits<double>::infinity();
    lwork=-1;
//  [y\\z]=V'(P'b):
    F77NAME(dormrz)('L','N',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),xtmp.addr(),n,&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    double *work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dormrz)('L','N',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),xtmp.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R' v = y
    F77NAME(dtrsv)('U','T','N',the_rank,URV->addr(),m,xtmp.addr(),1);

    x=double_zero_;
    x.copyFrom(the_rank,xtmp);
//  U [v\\0]
    w=numeric_limits<double>::infinity();
    lwork=-1;
    F77NAME(dormqr)('L','N',m,1,mn,URV->addr(),m,utau->addr(),
                    x.addr(),m,&w,lwork,info );
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dormqr)('L','N',m,1,mn,URV->addr(),m,utau->addr(),
                    x.addr(),m,work,lwork,info );
    CHECK_SAME(info,0)
    delete [] work; work=0;
    residual_norm=F77NAME(dnrm2)(n-the_rank,xtmp.addr(the_rank),1);
  }
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(dlascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
  return residual_norm;
}

template<> void CompleteOrthogonalDecomposition<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,
Vector<double,double> &residual_norm,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve");
#ifdef DEBUG
//cout << "\tB = " << endl;
//B.printOn(cout);
#endif
  int m=URV->size(0),n=URV->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  CHECK_SAME(k,residual_norm.size())
  double w=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(dlange)('M',brow,k,B.addr(),B.size(0),&w);
  int ibscl=0;
  double bscl=double_one_;
  if (bnrm>double_zero_ && bnrm<smlnum) {
    ibscl=1;
    bscl=smlnum;
  } else if (bnrm > bignum) {
    ibscl=2;
    bscl=bignum;
  }

  int info;
  int mn=min(m,n);
  int lwork=-1;
  int scllen=0;
  residual_norm=numeric_limits<double>::infinity();
  if (tr==Factorization::NO_TRANSPOSE) { // min || B - A X || min || X ||
//  TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve NO_TRANSPOSE");
//  if A=U[R 0]V'P', R=B-Ax, 0=A'R and X=A'S then
//        [0 0]
//    [F]=U'R=U'B-[R 0]V'P'X=[Y]-[R 0][V]=[Y-Rv] ==> G=Z
//    [G]         [0 0]      [Z] [0 0][W] [Z]
//  and
//    0=A'R=PV[R' 0]U'R=PV[R' 0][F]=PV[R'F] ==> F=0
//            [0  0]      [0  0][G]   [0  ]
//  and
//    [V]=V'P'X=[R' 0]U'S=[R' 0][p]=[R'p] ==> W=0
//    [W]       [0  0]    [0  0][q] [0  ]
//  so compute [Y\\Z]=U'B, solve Rv=Y, compute X=PV[V\\0]
    CHECK_SAME(m,B.size(0))
    CHECK_SAME(n,X.size(0))
    Matrix<double,double> Xtmp(m,k);
    Xtmp.copy(B);
    if (ibscl!=0) {
      F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
      CHECK_SAME(info,0)
    }
#ifdef DEBUG
//  cout << "\tXtmp after dlascl = " << endl;
//  Xtmp.printOn(cout);
#endif
//  U' B = [ Y \\ Z ]:
    F77NAME(dormqr)('L','T',m,k,mn,URV->addr(),m,utau->addr(),
                    Xtmp.addr(),m,&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    double *work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dormqr)('L','T',m,k,mn,URV->addr(),m,utau->addr(),
                    Xtmp.addr(),m,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
#ifdef DEBUG
//  cout << "\tXtmp after dormqr = " << endl;
//  Xtmp.printOn(cout);
#endif
//  solve R V = Y:
    F77NAME(dtrsm)('L','U','N','N',the_rank,k,double_one_,URV->addr(),m,
      Xtmp.addr(),m);
#ifdef DEBUG
//  cout << "\tXtmp after dtrsm = " << endl;
//  Xtmp.printOn(cout);
#endif
    Matrix<double,double> PX(n,k,double_zero_);
    PX.copyFrom('A',the_rank,k,Xtmp);
    double w=numeric_limits<double>::infinity();
    lwork=-1;
//  P'X=V[V\\0]:
    F77NAME(dormrz)('L','T',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),PX.addr(),n,&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dormrz)('L','T',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),PX.addr(),n,work,lwork,info);
#ifdef DEBUG
//  cout << "\tPX = " << endl;
//  PX.printOn(cout);
#endif
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  permute
    for (int j=0;j<k;j++) {
      for (int i=0;i<n;i++) X(jpvt[i]-1,j)=PX(i,j);
      residual_norm[j]=F77NAME(dnrm2)(m-the_rank,Xtmp.addr(the_rank,j),1);
    }
#ifdef DEBUG
//  cout << "\tX = " << endl;
//  X.printOn(cout);
#endif
  } else { // min || X || min || B - A' X ||
//  if A=U[R 0]V'P', R=B-A'X, 0=Ar and X=AS then
//          [0 0]
//    [F]=V'P'R=V'P'B-[R' 0]U'X=[Y]-[R' 0][V]=[Y-Rv] ==> G=Z, R'V=Y
//    [G]             [0  0]    [Z] [0  0][W] [Z   ]
//  and
//    0=Ar=U[R 0]V'P'R=U[R 0][F]=U[Rf] ==> F=0
//          [0 0]       [0 0][G]  [0 ]
//  and
//    [V]=Ux=[R 0]V'P'S=[R 0][p]=[Rp] ==> W=0
//    [W]    [0 0]      [0 0][q] [0 ]
//  so compute [Y\\Z]=V'P'B, solve R'V=Y, compute X=U[V\\0]
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(m,X.size(0))
    Matrix<double,double> Xtmp(n,k);
//  P'B:
    for (int j=0;j<k;j++) {
      for (int i=0;i<n;i++) Xtmp(i,j)=B(jpvt[i]-1,j);
    }
    if (ibscl!=0) {
      F77NAME(dlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
      CHECK_SAME(info,0)
    }
    double w=numeric_limits<double>::infinity();
    lwork=-1;
//  [Y\\Z]=V'(P'B):
    F77NAME(dormrz)('L','N',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),Xtmp.addr(),n,&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    double *work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dormrz)('L','N',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),Xtmp.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R' V = B
    F77NAME(dtrsm)('L','U','T','N',the_rank,k,double_one_,URV->addr(),m,
      Xtmp.addr(),n);

    X=double_zero_;
    X.copyFrom('A',the_rank,k,Xtmp);
//  U [V\\0]
    w=numeric_limits<double>::infinity();
    lwork=-1;
    F77NAME(dormqr)('L','N',m,k,mn,URV->addr(),m,utau->addr(),
                    X.addr(),m,&w,lwork,info );
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dormqr)('L','N',m,k,mn,URV->addr(),m,utau->addr(),
                    X.addr(),m,work,lwork,info );
    CHECK_SAME(info,0)
    delete [] work; work=0;
    for (int j=0;j<k;j++) {
      residual_norm[j]=F77NAME(dnrm2)(n-the_rank,Xtmp.addr(the_rank,j),1);
    }
  }
  if (iascl!=0) {
    F77NAME(dlascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(dlascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
}

template class CompleteOrthogonalDecomposition<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "GramSchmidtQRFactorization.C"
template<> double GramSchmidtQRFactorization<double,double>::tol=
  double_one_/sqrt(2.);

template<> void 
GramSchmidtQRFactorization<double,double>::reorthogonalize(int j) {
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,R->size(0))
  CHECK_BOUNDS(j,0,n)
  for (int k=0;k<j;k++) {
    double dot=F77NAME(ddot)(m,Q->addr(0,k),1,Q->addr(0,j),1)/(*R)(k,k);
    (*R)(k,j)+=dot;
    F77NAME(daxpy)(m,-dot,Q->addr(0,k),1,Q->addr(0,j),1);
  }
  (*R)(j,j)=F77NAME(ddot)(m,Q->addr(0,j),1,Q->addr(0,j),1);
}

template<> 
GramSchmidtQRFactorization<double,double>::GramSchmidtQRFactorization(
const Matrix<double,double> &A) : Q(0),R(0),A_original(&A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::GramSchmidtQRFactorization");
  int m=A.size(0),n=A.size(1);
  assert(m>=n);
  Q=OPERATOR_NEW OrthogonalMatrix<double,double>(m,n);
  R=OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(n,n);
  Q->copy(A);

  double *norm=OPERATOR_NEW_BRACKET(double,n);
  for (int j=0;j<n;j++) {
    norm[j]=F77NAME(ddot)(m,Q->addr(0,j),1,Q->addr(0,j),1);
  }
  double tol2=tol*tol;
  for (int k=0;k<n;k++) {
    (*R)(k,k)=F77NAME(ddot)(m,Q->addr(0,k),1,Q->addr(0,k),1);
    if ((*R)(k,k)<tol2*norm[k]) reorthogonalize(k);
    if ((*R)(k,k)>double_zero_ && k<n-1) {
      F77NAME(dgemv)('T',m,n-k-1,double_one_/(*R)(k,k),Q->addr(0,k+1),m,
        Q->addr(0,k),1,double_zero_,R->addr(k,k+1),n);
      F77NAME(dger)(m,n-k-1,double_mone_,Q->addr(0,k),1,R->addr(k,k+1),n,
        Q->addr(0,k+1),m);
    }
  }
  delete [] norm;
}

template<> void
GramSchmidtQRFactorization<double,double>::solveOverdetermined(
const Vector<double,double> &b,Vector<double,double> &x,
Vector<double,double> &residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,b.size());
  CHECK_SAME(n,x.size());
  CHECK_SAME(m,residual.size());
  residual.copy(b);
  for (int j=0;j<n;j++) {
    double dot=-F77NAME(ddot)(m,Q->addr(0,j),1,residual.addr(),1)
              /(*R)(j,j);
    x[j]=-dot;
    F77NAME(daxpy)(m,dot,Q->addr(0,j),1,residual.addr(),1);
  }
  F77NAME(dtrsv)('U','N','U',n,R->addr(),n,x.addr(),1);
}

template<> void
GramSchmidtQRFactorization<double,double>::solveOverdetermined(
const Matrix<double,double> &B,Matrix<double,double> &X,
Matrix<double,double> &Residual) const { 
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,X.size(0));
  CHECK_SAME(k,X.size(1));
  CHECK_SAME(m,Residual.size(0));
  CHECK_SAME(k,Residual.size(1));
  Residual.copy(B);
  for (int j=0;j<n;j++) {
    F77NAME(dgemv)('T',m,k,double_one_/(*R)(j,j),Residual.addr(),m,
      Q->addr(0,j),1,double_zero_,X.addr(j,0),n);
    F77NAME(dger)(m,k,double_mone_,Q->addr(0,j),1,X.addr(j,0),n,
      Residual.addr(),m);
  }
  F77NAME(dtrsm)('L','U','N','U',n,k,double_one_,R->addr(),n,X.addr(),n);
}

template<> void
GramSchmidtQRFactorization<double,double>::solveUnderdetermined(
const Vector<double,double> &b,Vector<double,double> &x) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveUnderdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,b.size());
  CHECK_SAME(m,x.size());

  Vector<double,double> *v=OPERATOR_NEW Vector<double,double>(n);
  v->copy(b);

  F77NAME(dtrsv)('U','T','U',n,R->addr(),n,v->addr(),1);
  x=double_zero_;
  for (int j=0;j<n;j++) {
    double omega=((*v)[j]-F77NAME(ddot)(m,Q->addr(0,j),1,x.addr(),1))
      / (*R)(j,j);
    F77NAME(daxpy)(m,omega,Q->addr(0,j),1,x.addr(),1);
  }
  delete v; v=0;
}

template<>
void GramSchmidtQRFactorization<double,double>::solveUnderdetermined(
const Matrix<double,double> &B,Matrix<double,double> &X) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveUnderdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(n,B.size(0));
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(k,X.size(1));

  Matrix<double,double> *V=OPERATOR_NEW Matrix<double,double>(n,k);
  V->copy(B);
  F77NAME(dtrsm)('L','U','T','U',n,k,double_one_,R->addr(),n,
    V->addr(),n);
  X=double_zero_;
  double *work=OPERATOR_NEW_BRACKET(double,k);
  for (int j=0;j<n;j++) {
    F77NAME(dgemv)('T',m,k,double_one_,X.addr(),m,Q->addr(0,j),1,
      double_zero_,work,1);
    for (int kk=0;kk<k;kk++) work[kk]=((*V)(j,kk)-work[kk])/(*R)(j,j);
    F77NAME(dger)(m,k,double_one_,Q->addr(0,j),1,work,1,X.addr(),m);
  }
  delete [] work; work=0;
  delete V; V=0;
}

template<> void
GramSchmidtQRFactorization<double,double>::improveOverdetermined(
const Vector<double,double> &b,Vector<double,double> &x,
Vector<double,double> &residual) const { 
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveOverdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,b.size());
  CHECK_SAME(n,x.size());
  CHECK_SAME(m,residual.size());
  Vector<double,double> *delta_b=OPERATOR_NEW Vector<double,double>(m);
  Vector<double,double> *delta_z=OPERATOR_NEW Vector<double,double>(n,0.);
  Vector<double,double> *delta_y=OPERATOR_NEW Vector<double,double>(n,0.);
  F77NAME(dgemv)('N',m,n,double_mone_,A_original->addr(),m,x.addr(),1,
    double_zero_,delta_b->addr(),1);
  F77NAME(dgemv)('T',m,n,double_mone_,A_original->addr(),m,
    residual.addr(),1,double_zero_,delta_z->addr(),1);
  for (int i=0;i<m;i++) (*delta_b)[i]+=b[i]-residual[i];
  F77NAME(dtrsv)('U','T','U',n,R->addr(),n,delta_z->addr(),1);
  for (int j=0;j<n;j++) {
    (*delta_y)[j]=(F77NAME(ddot)(m,Q->addr(0,j),1,delta_b->addr(),1)
                 -(*delta_z)[j])/(*R)(j,j);
    F77NAME(daxpy)(m,-(*delta_y)[j],Q->addr(0,j),1,delta_b->addr(),1);
  }
  F77NAME(dtrsv)('U','N','U',n,R->addr(),n,delta_y->addr(),1);
  x+=(*delta_y);
  residual+=(*delta_b);
  delete delta_y; delta_y=0;
  delete delta_z; delta_z=0;
  delete delta_b; delta_b=0;
}

template<> void
GramSchmidtQRFactorization<double,double>::improveOverdetermined(
const Matrix<double,double> &B,Matrix<double,double> &X,
Matrix<double,double> &Residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,X.size(0));
  CHECK_SAME(k,X.size(1));
  CHECK_SAME(m,Residual.size(0));
  CHECK_SAME(k,Residual.size(1));
  Matrix<double,double> *delta_B=OPERATOR_NEW Matrix<double,double>(m,k);
  Matrix<double,double> *delta_Z=
    OPERATOR_NEW Matrix<double,double>(n,k,0.);
  Matrix<double,double> *delta_Y=
    OPERATOR_NEW Matrix<double,double>(n,k,0.);
  F77NAME(dgemm)('N','N',m,k,n,double_mone_,A_original->addr(),m,
    X.addr(),n,double_zero_,delta_B->addr(),m);
  F77NAME(dgemm)('T','N',n,k,m,double_mone_,A_original->addr(),m,
    Residual.addr(),m,double_zero_,delta_Z->addr(),n);
  for (int j=0;j<k;j++) {
    for (int i=0;i<m;i++) (*delta_B)(i,j)+=B(i,j)-Residual(i,j);
  }
  F77NAME(dtrsm)('L','U','T','U',n,k,double_one_,R->addr(),n,
    delta_Z->addr(),n);
  double *work=OPERATOR_NEW_BRACKET(double,k);
  for (int j=0;j<n;j++) {
    F77NAME(dgemv)('T',m,k,double_one_,delta_B->addr(),m,Q->addr(0,j),1,
      double_zero_,work,1);
    for (int kk=0;kk<k;kk++) {
      (*delta_Y)(j,kk)=(work[kk]-(*delta_Z)(j,kk))/(*R)(j,j);
    }
    F77NAME(dger)(m,k,double_mone_,Q->addr(0,j),1,delta_Y->addr(j,0),n,
      delta_B->addr(),m);
  }
  delete [] work; work=0;
  F77NAME(dtrsm)('L','U','N','U',n,k,double_one_,R->addr(),n,
    delta_Y->addr(),n);
  X+=(*delta_Y);
  Residual+=(*delta_B);
  delete delta_Y; delta_Y=0;
  delete delta_Z; delta_Z=0;
  delete delta_B; delta_B=0;
}

template<>
void GramSchmidtQRFactorization<double,double>::improveUnderdetermined(
const Vector<double,double> &b,Vector<double,double> &x) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveUnderdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,b.size());
  CHECK_SAME(m,x.size());

  Vector<double,double> *y=OPERATOR_NEW Vector<double,double>(n);
  Vector<double,double> *delta_x=OPERATOR_NEW Vector<double,double>(m);
  Vector<double,double> *delta_by=OPERATOR_NEW Vector<double,double>(n);
  y->copy(b);
  delta_x->copy(x);
  delta_by->copy(b);
  F77NAME(dtrsv)('U','T','U',n,R->addr(),n,y->addr(),1); // R^T y=b
  for (int j=0;j<n;j++) (*y)[j]/=(*R)(j,j);
  F77NAME(dgemv)('N',m,n,double_one_,Q->addr(),m,y->addr(),1,
    double_mone_,delta_x->addr(),1); // delta_c = Q Sigma^{-2} y - x
  F77NAME(dgemv)('T',m,n,double_mone_,A_original->addr(),m,x.addr(),1,
    double_one_,delta_by->addr(),1); // delta_b = b - A x
  F77NAME(dtrsv)('U','T','U',n,R->addr(),n,delta_by->addr(),1);
    // R^T delta_y = delta_b
  for (int j=0;j<n;j++) {
    double omega=(F77NAME(ddot)(m,Q->addr(0,j),1,delta_x->addr(),1)
      -(*delta_by)[j])/(*R)(j,j);
    F77NAME(daxpy)(m,-omega,Q->addr(0,j),1,delta_x->addr(),1);
  }
  x+=(*delta_x);
  delete delta_by; delta_by=0;
  delete delta_x; delta_x=0;
  delete y; y=0;
}

template<>
void GramSchmidtQRFactorization<double,double>::improveUnderdetermined(
const Matrix<double,double> &B,Matrix<double,double> &X) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveUnderdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(n,B.size(0));
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(k,X.size(1));

  Matrix<double,double> *Y=OPERATOR_NEW Matrix<double,double>(n,k);
  Matrix<double,double> *delta_X=OPERATOR_NEW Matrix<double,double>(m,k);
  Matrix<double,double> *delta_BY=OPERATOR_NEW Matrix<double,double>(n,k);
  Y->copy(B);
  delta_X->copy(X);
  delta_BY->copy(B);

  F77NAME(dtrsm)('L','U','T','U',n,k,double_one_,R->addr(),n,Y->addr(),n);
    // R^T Y = B
  for (int j=0;j<n;j++) {
    F77NAME(dscal)(k,double_one_/(*R)(j,j),Y->addr(j,0),n);
  }
  F77NAME(dgemm)('N','N',m,k,n,double_one_,Q->addr(),m,Y->addr(),n,
    double_mone_,delta_X->addr(),m); // delta_C = Q Sigma^{-2} Y - X
  F77NAME(dgemm)('T','N',n,k,m,double_mone_,A_original->addr(),m,
    X.addr(),m,double_one_,delta_BY->addr(),n); // delta_B = B - A X
  F77NAME(dtrsm)('L','U','T','U',n,k,double_one_,R->addr(),n,
    delta_BY->addr(),n);
    // R^T delta_Y = delta_B
  double *work=OPERATOR_NEW_BRACKET(double,k);
  for (int j=0;j<n;j++) {
    F77NAME(dgemv)('T',m,k,double_one_,delta_X->addr(),m,Q->addr(0,j),1,
      double_zero_,work,1); // w^T = delta_X^T q_j
    for (int kk=0;kk<k;kk++) {
      work[kk]=(work[kk]-(*delta_BY)(j,kk))/(*R)(j,j);
    }
    F77NAME(dger)(m,k,double_mone_,Q->addr(0,j),1,work,1,
      delta_X->addr(),m);
  }
  X+=(*delta_X);
  delete [] work; work=0;
  delete delta_BY; delta_BY=0;
  delete delta_X; delta_X=0;
  delete Y; Y=0;
}

template<> void GramSchmidtQRFactorization<double,double>::exchangeColumn(
int j,int jAin,int jAout,const Matrix<double,double> &A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::exchangeColumn");
#ifdef DEBUG
//cout << "\tj,jAin,jAout = " << j << " " << jAin << " " << jAout << endl;
//cout << "\tA = " << endl;
//A.printOn(cout);
//printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(j,0,n)
  CHECK_BOUNDS(jAin,0,A.size(1))
  CHECK_BOUNDS(jAout,0,A.size(1))
  Vector<double,double> dif(m);
  for (int i=0;i<m;i++) dif[i]=A(i,jAin)-A(i,jAout);
#ifdef DEBUG
//cout << "\tdif = " << endl;
//dif.printOn(cout);
#endif
  Vector<double,double> QTdif(n,double_zero_);
  Q->gemv(double_one_,dif,double_zero_,QTdif,'T');
#ifdef DEBUG
//cout << "\tQTdif = " << endl;
//QTdif.printOn(cout);
#endif
  Vector<double,double> subdiag(n-1,double_zero_);
  for (int k=n-2;k>j;k--) {
    double c,s;
    F77NAME(drotg)(QTdif[k],QTdif[k+1],c,s);
#ifdef DEBUG
//  cout << "\tk,c,s = " << k << " " << c << " " << s << endl;
#endif
    int ncols=n-k-1;
    F77NAME(drot)(1,R->addr(k,k),n,subdiag.addr(k),1,c,s);
    F77NAME(drot)(n-k-1,R->addr(k,k+1),n,R->addr(k+1,k+1),n,c,s);
    F77NAME(drot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
#ifdef DEBUG
//  printOn(cout);
#endif
  }
  F77NAME(daxpy)(j+1,double_one_,QTdif.addr(),1,R->addr(0,j),1);
  if (j+1<n) {
    subdiag[j]=QTdif[j+1];
    for (int k=j;k<n-1;k++) {
      double c,s;
      F77NAME(drotg)((*R)(k,k),subdiag[k],c,s);
#ifdef DEBUG
//    cout << "\tk,c,s = " << k << " " << c << " " << s << endl;
#endif
      F77NAME(drot)(n-k-1,R->addr(k,k+1),n,R->addr(k+1,k+1),n,c,s);
      F77NAME(drot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
#ifdef DEBUG
//    printOn(cout);
#endif
    }
  }
#ifdef DEBUG
//printOn(cout);
#endif
}

template<> void GramSchmidtQRFactorization<double,double>::exchangeRow(
int i,int iAin,int iAout,const Matrix<double,double> &A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::exchangeRow");
#ifdef DEBUG
//cout << "\ti,iAin,iAout = " << i << " " << iAin << " " << iAout << endl;
//cout << "\tA = " << endl;
//A.printOn(cout);
//printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(iAin,0,A.size(0));
  CHECK_BOUNDS(iAout,0,A.size(0));
  Vector<double,double> dif(n);
  for (int j=0;j<n;j++) dif[j]=A(iAin,j)-A(iAout,j);
#ifdef DEBUG
//cout << "\tdif = " << endl;
//dif.printOn(cout);
#endif
  Vector<double,double> axis(m,double_zero_);
  axis[i]=double_one_;
  Vector<double,double> QTaxis(n,double_zero_);
  Q->gemv(double_one_,axis,double_zero_,QTaxis,'T');
#ifdef DEBUG
//cout << "\tQTaxis = " << endl;
//QTaxis.printOn(cout);
#endif
  Vector<double,double> subdiag(n-1,double_zero_);
  for (int k=n-2;k>=0;k--) {
    double c,s;
    F77NAME(drotg)(QTaxis[k],QTaxis[k+1],c,s);
#ifdef DEBUG
//  cout << "\tk,c,s = " << k << " " << c << " " << s << endl;
#endif
    int ncols=n-k-1;
    F77NAME(drot)(1,R->addr(k,k),n,subdiag.addr(k),1,c,s);
    F77NAME(drot)(n-k-1,R->addr(k,k+1),n,R->addr(k+1,k+1),n,c,s);
    F77NAME(drot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
#ifdef DEBUG
//  printOn(cout);
#endif
  }
  F77NAME(daxpy)(n,QTaxis[0],dif.addr(),1,R->addr(),n);
  for (int k=0;k<n-1;k++) {
    double c,s;
    F77NAME(drotg)((*R)(k,k),subdiag[k],c,s);
#ifdef DEBUG
//    cout << "\tk,c,s = " << k << " " << c << " " << s << endl;
#endif
      F77NAME(drot)(n-k-1,R->addr(k,k+1),n,R->addr(k+1,k+1),n,c,s);
      F77NAME(drot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
#ifdef DEBUG
//    printOn(cout);
#endif
    }
#ifdef DEBUG
//printOn(cout);
#endif
}

template<> void GramSchmidtQRFactorization<double,double>::dropColumn(
int j) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropColumn");
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(j,0,n)
  for (int k=j+1;k<n;k++) {
    double c,s;
    F77NAME(drotg)(R->operator()(k,k),R->operator()(j,k),c,s);
    int ncols=n-k-1;
    if (ncols>0) {
      F77NAME(drot)(ncols,R->addr(k,k+1),n,R->addr(j,k+1),n,c,s);
    }
    F77NAME(drot)(m,Q->addr(0,k),1,Q->addr(0,j),1,c,s);
  }

  for (int k=j+1;k<n;k++) {
    for (int i=0;i<j;i++) (*R)(i,k-1)=(*R)(i,k);
    for (int i=j;i<k;i++) (*R)(i,k-1)=(*R)(i+1,k);
    F77NAME(dcopy)(m,Q->addr(0,k),1,Q->addr(0,k-1),1);
  }

  OrthogonalMatrix<double,double> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,double>(m,n-1);
  new_Q->copyFrom('A',m,n-1,*Q);
  delete Q;
  Q=new_Q;
  UpperTrapezoidalMatrix<double,double> *new_R
    =OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(n-1,n-1);
  new_R->copyFrom(n-1,n-1,*R);
  delete R;
  R=new_R;
}

template<> void GramSchmidtQRFactorization<double,double>::addColumn(
int j,const Matrix<double,double> &A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::addColumn");
#ifdef DEBUG
//printOn(cout);
//cout << "\tj = " << j << endl;
//cout << "\tA = " << endl;
//A.printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,A.size(0))

  OrthogonalMatrix<double,double> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,double>(m,n+1);
  new_Q->copyFrom('A',m,n,*Q);
  delete Q;
  Q=new_Q;
  F77NAME(dcopy)(m,A.addr(0,j),1,Q->addr(0,n),1);

  UpperTrapezoidalMatrix<double,double> *new_R
    =OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(n+1,n+1);
  new_R->copyFrom(n,n,*R);
  delete R;
  R=new_R;

  double norm=F77NAME(dnrm2)(m,Q->addr(0,n),1);
  for (j=0;j<n;j++) {
    double dot=-F77NAME(ddot)(m,Q->addr(0,j),1,Q->addr(0,n),1);
    (*R)(j,n)=-dot;
    F77NAME(daxpy)(m,dot,Q->addr(0,j),1,Q->addr(0,n),1);
  }
  (*R)(n,n)=F77NAME(dnrm2)(m,Q->addr(0,n),1);
  if ((*R)(n,n)<tol*norm) reorthogonalize(n);
  double a=double_one_/(*R)(n,n);
  F77NAME(dscal)(m,a,Q->addr(0,n),1);
#ifdef DEBUG
//printOn(cout);
#endif
}

template<> void GramSchmidtQRFactorization<double,double>::dropRow(int i)
{
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropRow");
#ifdef DEBUG
//printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(i,0,m)
//permute row i to last
  for (int k=i+1;k<m;k++) {
    F77NAME(dswap)(n,Q->addr(k-1,0),m,Q->addr(k,0),m);
  }
#ifdef DEBUG
//Q->printOn(cout);
#endif
//the following from Daniel, Gragg, Kaufman and Stewart fails if
//  the last axis vector was already a column of Q
  Vector<double,double> work_column(m);
  int Mmin1=m-1;
  F77NAME(dgemv)('N',Mmin1,n,-double_one_,Q->addr(),m,Q->addr(m-1,0),m,
		 double_zero_,work_column.addr(),1);
  work_column[m-1]=sqrt(double_one_
	  -F77NAME(ddot)(n,Q->addr(m-1,0),n,Q->addr(m-1,0),m));
  double alpha=double_one_/work_column[m-1];
  F77NAME(dscal)(Mmin1,alpha,work_column.addr(),1);

  Vector<double,double> work_row(n,double_zero_);
  for (int j=n-1;j>=0;j--) {
    double c,s;
    F77NAME(drotg)(work_column[Mmin1],Q->operator()(Mmin1,j),c,s);
    F77NAME(drot)(Mmin1,work_column.addr(),1,Q->addr(0,j),1,c,s);
    int ncols=n-j;
    F77NAME(drot)(ncols,work_row.addr(j),1,R->addr(j,j),n,c,s);
  }

  OrthogonalMatrix<double,double> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,double>(Mmin1,n);
  new_Q->copyFrom('A',Mmin1,n,*Q);
  delete Q;
  Q=new_Q;
#ifdef DEBUG
//printOn(cout);
#endif
}

template<> void
GramSchmidtQRFactorization<double,double>::dropRowAndColumn(int i,int j)
{
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropRowAndColumn");
#ifdef DEBUG
//cout << "\ti,j = " << i << " " << j << endl;
//printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,n)
  Vector<double,double> rwork(n-1);
  UpperTrapezoidalMatrix<double,double> *new_R
    =OPERATOR_NEW UpperTrapezoidalMatrix<double,double>(n-1,n-1);
  for (int k=0;k<j;k++) {
    memcpy(new_R->addr(0,k),R->addr(0,k),(k+1)*sizeof(double));
  }
  for (int k=j+1;k<n;k++) {
    memcpy(new_R->addr(0,k-1),R->addr(0,k),k*sizeof(double));
    rwork[k-1]=(*R)(k,k);
  }
  for (int k=j;k<n-1;k++) {
    double c,s;
    F77NAME(drotg)((*new_R)(k,k),rwork[k],c,s);
    int ncols=n-k-2;
    if (ncols>0) {
      F77NAME(drot)(ncols,new_R->addr(k,k+1),n-1,
        new_R->addr(k+1,k+1),n-1,c,s);
    }
    F77NAME(drot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
  }
  OrthogonalMatrix<double,double> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,double>(m-1,n-1);
  for (int k=0;k<i;k++) {
    F77NAME(dcopy)(n-1,Q->addr(k,0),m,new_Q->addr(k,0),m-1);
  }
  for (int k=i+1;k<m;k++) {
    F77NAME(dcopy)(n-1,Q->addr(k-1,0),m,new_Q->addr(k-1,0),m-1);
  }
  Vector<double,double> qwork(m-1); // last column of Q
  memcpy(qwork.addr(),Q->addr(0,n-1),(m-1)*sizeof(double));
  rwork=double_zero_;
  for (int k=n-2;k>=0;k--) {
    double c,s;
    F77NAME(drotg)((*Q)(n-1,n-1),(*Q)(n-1,k),c,s);
    F77NAME(drot)(m-1,qwork.addr(),1,new_Q->addr(0,k),1,c,s);
    F77NAME(drot)(n-1-k,rwork.addr(k),1,new_R->addr(k,k),n-1,c,s);
  }
  delete Q;
  Q=new_Q;
  delete R;
  R=new_R;
#ifdef DEBUG
//cout << "\n\tfactorization at end" << endl;
//printOn(cout);
#endif
}

template<> void GramSchmidtQRFactorization<double,double>::addRow(int i,
const Matrix<double,double> &A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::addRow");
#ifdef DEBUG
//printOn(cout);
//cout << "\ti = " << i << endl;
//cout << "\tA = " << endl;
//A.printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,A.size(1))

  OrthogonalMatrix<double,double> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,double>(m+1,n);
  new_Q->copyFrom('A',m,n,*Q);
  delete Q;
  Q=new_Q;
  for (int j=0;j<n;j++) (*Q)(m,j)=double_zero_;

  Vector<double,double> work_row(n);
  F77NAME(dcopy)(n,A.addr(i,0),A.size(0),work_row.addr(),1);

  Vector<double,double> work_column(m+1,double_zero_);
  work_column[m]=double_one_;

  int Mp1=m+1;
  for (int j=0;j<n;j++) {
    double c,s;
    F77NAME(drotg)((*R)(j,j),work_row[j],c,s);
    int ncols=n-j-1;
    if (ncols>0) {
      F77NAME(drot)(ncols,R->addr(j,j+1),n,work_row.addr(j+1),1,c,s);
    }
    F77NAME(drot)(Mp1,Q->addr(0,j),1,work_column.addr(),1,c,s);
  }
#ifdef DEBUG
//printOn(cout);
#endif
}

template void testGramSchmidtQRFactorization(double,double);
template class GramSchmidtQRFactorization<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SingularValueDecomposition.C"
template<> SingularValueDecomposition<double,double>::
SingularValueDecomposition(const Matrix<double,double> &A) : U(0),
Vtranspose(0),s(0),A_original(&A) {
//TRACER_CALL(t,"SingularValueDecomposition::SingularValueDecomposition");
  int m=A.size(0),n=A.size(1);
  int mn=min(m,n);
  U=OPERATOR_NEW OrthogonalMatrix<double,double>(m,mn);
  Vtranspose=OPERATOR_NEW OrthogonalMatrix<double,double>(mn,n);
  s=OPERATOR_NEW Vector<double,double>(mn);
  Matrix<double,double> Acopy(m,n);
  Acopy.copy(A);

  int info=0;
  double w=numeric_limits<double>::infinity();
  int lwork=-1;
  F77NAME(dgesvd)('S','S',m,n,Acopy.addr(),m,s->addr(),U->addr(),m,
    Vtranspose->addr(),mn,&w,lwork,info);
  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dgesvd)('S','S',m,n,Acopy.addr(),m,s->addr(),U->addr(),m,
    Vtranspose->addr(),mn,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;
}

template<> void SingularValueDecomposition<double,double>::solve(
const Vector<double,double> &b,Vector<double,double> &x,
double rcond,Factorization::TRANSPOSE_OPTION to) const {
//TRACER_CALL(t,"SingularValueDecomposition::solve");
  CHECK_TEST(rcond<double_one_);
  int m=A_original->size(0),n=A_original->size(1);
  int mn=min(m,n);
  if (rcond<double_zero_) rcond=F77NAME(dlamch)('P');
  double thr=max(rcond*(*s)[0],F77NAME(dlamch)('S'));
  int rank=0;
  for (int i=0;i<mn;i++) {
    if ((*s)[i]>thr) rank++;
    else break;
  }
  double *y=OPERATOR_NEW_BRACKET(double,rank);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(dgemv)('C',m,rank,double_one_,U->addr(),m,b.addr(),1,
      double_zero_,y,1);
    for (int i=0;i<rank;i++) y[i]/=(*s)[i];
    F77NAME(dgemv)('C',rank,n,double_one_,Vtranspose->addr(),mn,y,1,
      double_zero_,x.addr(),1);
  } else {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(dgemv)('N',rank,n,double_one_,Vtranspose->addr(),mn,
      b.addr(),1,double_zero_,y,1);
    for (int i=0;i<rank;i++) y[i]/=(*s)[i];
    F77NAME(dgemv)('N',m,rank,double_one_,U->addr(),m,y,1,double_zero_,
      x.addr(),1);
  }
  delete [] y; y=0;
}

template<> void SingularValueDecomposition<double,double>::regularize(
const Vector<double,double> &b,Vector<double,double> &x,
double ridge,Factorization::TRANSPOSE_OPTION to) const {
  CHECK_TEST(ridge>=double_zero_);
  int m=A_original->size(0),n=A_original->size(1);
  int mn=min(m,n);
  double sr=sqrt(ridge);
  double *y=OPERATOR_NEW_BRACKET(double,mn);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(dgemv)('C',m,mn,double_one_,U->addr(),m,b.addr(),1,
      double_zero_,y,1);
    for (int i=0;i<mn;i++) {
      double si=(*s)[i];
      double scale=max(si,sr);
      si/=scale;
      double ri=sr/scale;
      y[i]*=si/(scale*(si*si+ri*ri));
    }
    F77NAME(dgemv)('C',mn,n,double_one_,Vtranspose->addr(),mn,y,1,
      double_zero_,x.addr(),1);
  } else {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(dgemv)('N',mn,n,double_one_,Vtranspose->addr(),mn,b.addr(),1,
      double_zero_,y,1);
    for (int i=0;i<mn;i++) {
      double si=(*s)[i];
      double scale=max(si,sr);
      si/=scale;
      double ri=sr/scale;
      y[i]*=si/(scale*(si*si+ri*ri));
    }
    F77NAME(dgemv)('N',m,mn,double_one_,U->addr(),m,y,1,double_zero_,
      x.addr(),1);
  }
  delete [] y; y=0;
}

template<> void SingularValueDecomposition<double,double>::solve(
const Matrix<double,double> &B,Matrix<double,double> &X,
double rcond,Factorization::TRANSPOSE_OPTION to) const {
//TRACER_CALL(t,"SingularValueDecomposition::solve");
  CHECK_TEST(rcond<double_one_);
  int m=A_original->size(0),n=A_original->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1));
  int mn=min(m,n);
  if (rcond<double_zero_) rcond=F77NAME(dlamch)('P');
  double thr=max(rcond*(*s)[0],F77NAME(dlamch)('S'));
  int rank=0;
  for (int i=0;i<mn;i++) {
    if ((*s)[i]>thr) rank++;
    else break;
  }
  double *Y=OPERATOR_NEW_BRACKET(double,rank*k);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,X.size(0));
    F77NAME(dgemm)('C','N',rank,k,m,double_one_,U->addr(),m,B.addr(),m,
      double_zero_,Y,rank);
    for (int i=0;i<rank;i++) F77NAME(drscl)(k,(*s)[i],Y+i,rank);
    F77NAME(dgemm)('C','N',n,k,rank,double_one_,Vtranspose->addr(),mn,
      Y,rank,double_zero_,X.addr(),n);
  } else {
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(m,X.size(0));
    F77NAME(dgemm)('N','N',rank,k,n,double_one_,Vtranspose->addr(),mn,
      B.addr(),n,double_zero_,Y,rank);
    for (int i=0;i<rank;i++) F77NAME(drscl)(k,(*s)[i],Y+i,rank);
    F77NAME(dgemm)('N','N',m,k,rank,double_one_,U->addr(),m,Y,rank,
      double_zero_,X.addr(),m);
  }
  delete [] Y; Y=0;
}

template<> void SingularValueDecomposition<double,double>::regularize(
const Matrix<double,double> &B,Matrix<double,double> &X,
double ridge,Factorization::TRANSPOSE_OPTION to) const {
  CHECK_TEST(ridge>=double_zero_);
  int m=A_original->size(0),n=A_original->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1));
  int mn=min(m,n);
  double sr=sqrt(ridge);
  double *Y=OPERATOR_NEW_BRACKET(double,mn*k);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,X.size(0));
    F77NAME(dgemm)('C','N',mn,k,m,double_one_,U->addr(),m,B.addr(),m,
      double_zero_,Y,mn);
    for (int i=0;i<mn;i++) {
      double si=(*s)[i];
      double scale=max(abs(si),sr);
      si/=scale;
      double ri=sr/scale;
      F77NAME(dscal)(k,si/(scale*(si*si+ri*ri)),Y+i,mn);
    }
    F77NAME(dgemm)('C','N',n,k,mn,double_one_,Vtranspose->addr(),mn,Y,mn,
      double_zero_,X.addr(),n);
  } else {
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(m,X.size(0));
    F77NAME(dgemm)('N','N',mn,k,n,double_one_,Vtranspose->addr(),mn,
      B.addr(),n,double_zero_,Y,mn);
    for (int i=0;i<mn;i++) {
      double si=(*s)[i];
      double scale=max(abs(si),sr);
      si/=scale;
      double ri=sr/scale;
      F77NAME(dscal)(k,si/(scale*(si*si+ri*ri)),Y+i,mn);
    }
    F77NAME(dgemm)('N','N',m,k,mn,double_one_,U->addr(),m,Y,mn,
      double_zero_,X.addr(),m);
  }
  delete [] Y; Y=0;
}

template class SingularValueDecomposition<double,double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> CholeskyFactorization<double,double>::CholeskyFactorization(
const SymmetricPositiveMatrix<double,double> &A) : anorm(double_zero_) {
//TRACER_CALL(t,"CholeskyFactorization::CF");
  int n=A.size(0);
  int inc=1;
  for (int j=0;j<n;j++) {
    int k=n-j;
    double col_sum=F77NAME(dasum)(k,A.addr(j,j),inc);
    col_sum+=F77NAME(dasum)(j,A.addr(j,0),n);
    if (col_sum>anorm) anorm=col_sum;
  }
  SymmetricPositiveMatrix<double,double> *Acopy=
    OPERATOR_NEW SymmetricPositiveMatrix<double,double>(n);
  Acopy->copy(A);
  int info;
  char uplo='L';
  F77NAME(dpotrf)(uplo,n,Acopy->addr(),n,info);
  CHECK_SAME(info,0)
  L=OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(n,n);
  for (int j=0;j<n;j++) {
    for (int i=j;i<n;i++) (*L)(i,j)=(*Acopy)(i,j);
  }
  delete Acopy;
}

template<> Matrix<double,double>*
CholeskyFactorization<double,double>::solveAXeqB(
const Matrix<double,double> &B) const {
//TRACER_CALL(t,"CholeskyFactorization::solveAXeqB");
  int n=L->size(0); 
  CHECK_SAME(n,B.size(0))
  int k=B.size(1);

  Matrix<double,double> *X=OPERATOR_NEW Matrix<double,double>(n,k);
  X->copy(B);
  char uplo='L';
  int info;
  F77NAME(dpotrs)(uplo,n,k,L->addr(),n,X->addr(),n,info);
  CHECK_SAME(info,0)
  return X;
}

template<> double
CholeskyFactorization<double,double>::conditionNumber() const {
  int n=L->size(0);
  double rcond=double_zero_;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  char uplo='L';
  F77NAME(dpocon)(uplo,n,L->addr(),n,anorm,rcond,work,iwork,
                  info);
  CHECK_SAME(info,0)
  delete [] iwork;
  delete [] work;
  if (rcond>double_zero_) return double_one_/rcond;
  return numeric_limits<double>::infinity();
}

template class CholeskyFactorization<double,double>;
//template ostream& operator<<(ostream&,
//  const CholeskyFactorization<double,double>&);
template void testCholeskyFactorization(double,double);
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> MDMtFactorization<double,double>::MDMtFactorization(
const SymmetricMatrix<double,double> &A) : anorm(double_zero_) {
//TRACER_CALL(t,"MDMtFactorization::MDMtF");
  int n=A.size(0);
  ipiv=OPERATOR_NEW_BRACKET(int,n);
  int inc=1;
  for (int j=0;j<n;j++) {
    int k=n-j;
    double col_sum=F77NAME(dasum)(k,A.addr(j,j),inc);
    col_sum+=F77NAME(dasum)(j,A.addr(j,0),n);
    if (col_sum>anorm) anorm=col_sum;
  }
  SymmetricMatrix<double,double> *Acopy=
    OPERATOR_NEW SymmetricMatrix<double,double>(n);
  Acopy->copy(A);
  int info;
  char uplo='L';
  double w=numeric_limits<double>::infinity();
  int lwork=-1;
  F77NAME(dsytrf)(uplo,n,Acopy->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dsytrf)(uplo,n,Acopy->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<double,double>(n,n);
  for (int j=0;j<n;j++) {
    for (int i=j;i<n;i++) (*L)(i,j)=(*Acopy)(i,j);
  }
  delete Acopy;
}

template<> Matrix<double,double>*
MDMtFactorization<double,double>::solveAXeqB(
const Matrix<double,double> &B) const {
//TRACER_CALL(t,"MDMtFactorization::solveAXeqB");
  int n=L->size(0);
  CHECK_SAME(n,B.size(0))
  int k=B.size(1);

  Matrix<double,double> *X=OPERATOR_NEW Matrix<double,double>(n,k);
  X->copy(B);
  char uplo='L';
  int info;
  F77NAME(dsytrs)(uplo,n,k,L->addr(),n,ipiv,X->addr(),n,info);
  CHECK_SAME(info,0)
  return X;
}

template<> double MDMtFactorization<double,double>::conditionNumber()
const {
  int n=L->size(0);
  double rcond=double_zero_;
  double *work=OPERATOR_NEW_BRACKET(double,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  char uplo='L';
  F77NAME(dsycon)(uplo,n,L->addr(),n,ipiv,anorm,rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork;
  delete [] work;
  if (rcond>double_zero_) return double_one_/rcond;
  return numeric_limits<double>::infinity();
}

template class MDMtFactorization<double,double>;
//template ostream& operator<<(ostream&,
//  const MDMtFactorization<double,double>&);
template void testMDMtFactorization(double,double);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<>
SingularValueDecomposition<double,double>::SingularValueDecomposition(
const Matrix<double,double> &A) {
//TRACER_CALL(t,"SingularValueDecomposition::SVD");
  int m=A.size(0),n=A.size(1),minmn=min(m,n);
  Sigma=OPERATOR_NEW Vector<double,double>(minmn);
  int info=0;
  double *work=0;
  int lwork=5*max(m,n);
  if (m>=n) {
    U=OPERATOR_NEW Matrix<double,double>(m,n);
    U->copy(A);
    V_transpose=OPERATOR_NEW Matrix<double,double>(n,n);
    char jobu='O';
    char jobvt='S';
//  int lwork=max(1,max( 3*n + m, 5*n - 4));
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dgesvd)(jobu,jobvt,m,n,U->addr(),m,Sigma->addr(),
      0,m,V_transpose->addr(),n,work,lwork,info);
  } else {
    U=OPERATOR_NEW Matrix<double,double>(m,m);
    V_transpose=OPERATOR_NEW Matrix<double,double>(m,n);
    V_transpose->copy(A);
    char jobu='S';
    char jobvt='O';
//  int lwork=max(1,max( 3*m + n, 5*m - 4));
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dgesvd)(jobu,jobvt,m,n,V_transpose->addr(),m,
      Sigma->addr(),U->addr(),m,0,m,work,lwork,info);
  }
  CHECK_SAME(info,0)
  if (work) delete work;
}

template<> double SingularValueDecomposition<double,double>::conditionNumber(
double cutoff) const {
  CHECK_POSITIVE(cutoff)
  int r=rank(cutoff);
  if (r<1) return numeric_limits<double>::infinity();
  if (r==1) return double_one_;
  return (*Sigma)[0]/(*Sigma)[r-1];
}

template<> Matrix<double,double>* SingularValueDecomposition<double,double>::solveAXeqB(
const Matrix<double,double> &B,double cutoff,transpose_option option) const {
//TRACER_CALL(t,"SingularValueDecomposition::solveAXeqB");
  int m=U->size(0),n=V_transpose->size(1),minmn=Sigma->size();
  int K=B.size(1);
  int r=rank(cutoff);
  double alpha=double_one_;
  double beta=double_zero_;
  char trans='N';
  switch (option) {
    case no_matrix_transpose: {
      CHECK_SAME(m,B.size(0))
      Matrix<double,double> C(r,K);
      char transa='T';
      F77NAME(dgemm)(transa,trans,r,K,m,alpha,U->addr(),m,
        B.addr(),m,beta,C.addr(),r);
      for (int i=0;i<r;i++) {
        double temp=1./(*Sigma)[i];
        F77NAME(dscal)(K,temp,C.addr(i,0),r);
      }
      Matrix<double,double> *X=OPERATOR_NEW Matrix<double,double>(n,K);
      F77NAME(dgemm)(transa,trans,n,K,r,alpha,
        V_transpose->addr(),minmn,C.addr(),r,beta,X->addr(),n);
      return X;
      break;
    }
    case matrix_transpose: {
      CHECK_SAME(n,B.size(0))
      Matrix<double,double> C(r,K);
      F77NAME(dgemm)(trans,trans,r,K,n,alpha,V_transpose->addr(),
        n,B.addr(),n,beta,C.addr(),r);
      for (int i=0;i<r;i++) {
        double temp=1./(*Sigma)[i];
        F77NAME(dscal)(K,temp,C.addr(i,0),r);
      }
      Matrix<double,double> *X=OPERATOR_NEW Matrix<double,double>(m,K);
      F77NAME(dgemm)(trans,trans,m,K,r,alpha,U->addr(),m,
        C.addr(),r,beta,X->addr(),m);
      return X;
      break;
    }
    default:
      return 0;
  }
}

template class SingularValueDecomposition<double,double>;
//template ostream& operator<<(ostream&,
//  const SingularValueDecomposition<double,double>&);
template void testSingularValueDecomposition(double,double);
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//#include "LinearProgram.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//template<> double VirtualLinearProgram<double>::zero_ = double_zero_;
//template<> double VirtualLinearProgram<double>::one_ = double_one_;
//template<> double VirtualLinearProgram<double>::huge_ = 
//  numeric_limits<double>::infinity();

//template class VirtualLinearProgram<double>;
//template ostream& operator<<(ostream&,
//  const VirtualLinearProgram<double>&);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> int SFLinearProgram<double>::smallestDualSlack() {
//TRACER_CALL(t,"SFLinearProgram::smallestDualSlack");
  int M=A->size(0),N=A->size(1);
  int jmin=0;
  double rmin=double_zero_;
  for (int j=0;j<N-M;j++) {
    int colj=column[j+M];
    double rj=(*c)[colj]-F77NAME(ddot)(M,y->addr(),1,A->addr(0,colj),1);
    if (rj < rmin) { jmin=j; rmin=rj; }
    (*r)[j]=rj;
  }
#ifdef DEBUG
//cout << "\tr = " << endl;
//r->printOn(cout);
#endif
  return F77NAME(idmin)(N-M,r->addr(),1)-1;
}
*/

/*
//if perturbation in single component c_i between lower[i] and upper[i],
//then optimal x unchanged
template<> void SFLinearProgram<double>::costBounds(
Vector<double,double> &lower,Vector<double,double> &upper) const {
  int M=A->size(0),N=A->size(1);
  CHECK_SAME(N,lower.size())
  CHECK_SAME(N,upper.size())

//perturbation in basic variable
  for (int i=0;i<M;i++) {
    double low=-huge_;
    double high=huge_;
    Vector<double,double> axisi(M,this->zero_);
    axisi[i]=this->one_;
    Vector<double,double> hi(M);
    QR->solveUnderdetermined(axisi,hi);
    for (int j=M;j<N;j++) {
      int colj=column[j];
      double denom=F77NAME(ddot)(M,hi.addr(),1,A->addr(0,colj),1);
      double rj=(*r)[j-M];
      if (denom<this->zero_ && denom*low>rj) low=rj/denom;
      else if (denom>this->zero_ && denom*high>rj) high=rj/denom;
    }
    int coli=column[i];
    lower[coli]=(*c)[coli]+low;
    upper[coli]=(*c)[coli]+high;
  }
//perturbation in non-basic variable
  for (int j=M;j<N;j++) {
    int colj=column[j];
    lower[colj]=(*c)[colj]-(*r)[j-M];
    upper[colj]=huge_;
  }
}
*/

/*
template<> void SFLinearProgram<double>::arrayBounds(int i,int j,
double &lower,double &upper) const {
//TRACER_CALL(t,"SFLinearProgram::arrayBounds");
  int m=A->size(0),n=A->size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,n)
  char safe='S';
  double sfmin=F77NAME(dlamch)(&safe);
  upper=double_one_/sfmin;
  lower=-upper;
  int incA=1;

//perturbation in basic column of A
  double yi=(*y_trans)(0,i);
  double xj=(*x)(j,0);
  for (int jj=0;jj<m;jj++) {
    if (column[jj]==j) {
      double ainv_ji=(*basic_inverse)(jj,i);
      if (fabs(ainv_ji)>sfmin) {
	double bound=-double_one_/ainv_ji;
	if (ainv_ji>double_zero_) lower=bound;
	else upper=bound;
      }
      for (int k=0;k<m;k++) {
        if (k!=jj) {
	  double xk=(*x)(column[k],0);
	  double denom=(*basic_inverse)(k,i)*xj-ainv_ji*xk;
	  if (fabs(denom)>xk*sfmin) {
	    double bound=xk/denom;
	    if (denom>double_zero_ && bound<upper) upper=bound;
	    if (denom<double_zero_ && bound>lower) lower=bound;
	  }
	}
      }
      for (int kk=m;kk<n;kk++) {
	double rk=(*r_trans)(0,kk-m);
	double term=F77NAME(ddot)(m,basic_inverse->addr(jj,0),m,
			   A->addr(0,column[kk]),incA);
	double denom=yi*term-ainv_ji*rk;
	if (fabs(denom)>rk*sfmin) {
	  double bound=rk/denom;
	  if (denom>double_zero_ && bound<upper) upper=bound;
	  if (denom<double_zero_ && bound>lower) lower=bound;
	}
      }
      lower += (*A)(i,j);
      upper += (*A)(i,j);
      return;
    }
  }
//perturbation in non-basic column of A
  { for (int jj=m;jj<n;jj++) {
    if (column[jj]==j) {
      double rj=(*r_trans)(0,jj-m);
      if (fabs(yi)>rj*sfmin) {
	double bound=(*A)(i,j)+rj/yi;
	if (yi>double_zero_ && bound<upper) upper=bound;
	if (yi<double_zero_ && bound>lower) bound=lower;
      }
    }
  } }
}
*/
//template class SFLinearProgram<double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> LinearProgram<double>::LinearProgram(
const Matrix<double,double> &Ain,const Vector<double,double> &bin,
const Vector<double,double> &cin) : 
VirtualLinearProgram<double>(Ain,bin,cin),column(0),row(0),QR(0),
bbasic(0),cbasic(0),xbasic(0),ybasic(0),r(0),s(0),basic_number(0) {
//TRACER_CALL(t,"LinearProgram::LinearProgram");
  int M=A->size(0),N=A->size(1);
#ifdef DEBUG
//cout << "\tA = " << endl;
//Ain.printOn(cout);
//cout << "\tb = " << endl;
//bin.printOn(cout);
//cout << "\tidsumn(b) = " << F77NAME(idsumn)(M,bin.addr(),1) << endl;
//cout << "\tc = " << endl;
//cin.printOn(cout);
//cout << "\tidsump(c) = " << F77NAME(idsump)(N,cin.addr(),1) << endl;
#endif

  row=OPERATOR_NEW_BRACKET(int,M);
  for (int i=0;i<M;i++) row[i]=i;
  column=OPERATOR_NEW_BRACKET(int,N);
  for (int j=0;j<N;j++) column[j]=j;
  if (F77NAME(idsumn)(M,bin.addr(),1)>=M) {
    current_status=PRIMAL_FEASIBLE; // x=0 is feasible, all nonbasic
  }
  if (F77NAME(idsump)(N,cin.addr(),1)>=N) {
    current_status=DUAL_FEASIBLE; // y=0 is feasible, all nonbasic
  }
#ifdef DEBUG
//cout << "\tcurrent_status = " << current_status << endl;
#endif
  if (current_status!=UNKNOWN) computeSolution();
}
*/

/*
template<> void LinearProgram<double>::specifyBasicVariables(
int basic_no,int *basic_row_index,int *basic_column_index) {
//TRACER_CALL(t,"LinearProgram::specifyBasicVariables");
#ifdef DEBUG
//cout << "\tbasic_no = " << basic_no << endl;
//cout << "\tbasic_row_index = " << endl;
//for (int i=0;i<basic_no;i++) cout << basic_row_index[i] << " ";
//cout << endl;
//cout << "\tbasic_column_index = " << endl;
//for (int j=0;j<basic_no;j++) cout << basic_column_index[j] << " ";
//cout << endl;
#endif
  CHECK_TEST(this->current_status==UNKNOWN);
  CHECK_TEST(basic_no>0);
  this->basic_number=basic_no;
  int M=this->A->size(0),N=this->A->size(1);
  CHECK_TEST(basic_no<=min(M,N));
#ifdef DEBUG
  for (int i=1;i<basic_no;i++) {
    CHECK_TEST(basic_row_index[i-1]<basic_row_index[i]);
    CHECK_TEST(basic_column_index[i-1]<basic_column_index[i]);
  }
#endif
  SquareMatrix<double,double> Abasic(basic_no);
  bbasic=OPERATOR_NEW Vector<double,double>(basic_no);
  cbasic=OPERATOR_NEW Vector<double,double>(basic_no);
  for (int j=0;j<basic_no;j++) {
    (*bbasic)[j]=(*this->b)[basic_row_index[j]];
    (*cbasic)[j]=(*this->c)[basic_column_index[j]];
    for (int i=0;i<basic_no;i++) {
      Abasic(i,j)=(*this->A)(basic_row_index[i],basic_column_index[j]);
    }
  }
  QR=OPERATOR_NEW GramSchmidtQRFactorization<double,double>(Abasic);

  for (int j=0;j<basic_no;j++) {
    int t=column[j];
    column[j]=basic_column_index[j];
    column[basic_column_index[j]]=t;
    t=row[j];
    row[j]=basic_row_index[j];
    row[basic_row_index[j]]=t;
  }
  computeSolution();

//find cost reduction vector
  if (F77NAME(idsump)(basic_no,xbasic->addr(),1)>=basic_no) {
    if (basic_no>=M) {
      this->current_status=PRIMAL_FEASIBLE;
      return;
    }
    if (F77NAME(idsump)(M-basic_no,s->addr(),1)>=M-basic_no) {
      this->current_status=PRIMAL_FEASIBLE;
      return;
    }
  }
  if (F77NAME(idsump)(basic_no,ybasic->addr(),1)>=basic_no) {
    if (basic_no>=N) {
      this->current_status=DUAL_FEASIBLE;
      return;
    }
    if (F77NAME(idsump)(N-basic_no,r->addr(),1)>=N-basic_no) {
      this->current_status=DUAL_FEASIBLE;
      return;
    }
  }
  CHECK_TEST(0); // basic indices not basic
#ifdef DEBUG
//printOn(cout);
#endif
}
*/

/*
template<> STATUS_OPTION LinearProgram<double>::primalSimplexStep() {
//TRACER_CALL(t,"LinearProgram::primalSimplexStep");
  if (this->current_status==OPTIMAL) return OPTIMAL;
  CHECK_TEST(this->current_status==PRIMAL_FEASIBLE);
#ifdef DEBUG
//printOn(cout);
#endif

//find largest feasible simplex step
  int M=this->A->size(0),N=this->A->size(1);
  int i_basic=-1,j_non_basic=-1;
  double rmin=this->huge_,ymin=this->huge_;
  if (basic_number>0) {
    i_basic=F77NAME(idmin)(basic_number,ybasic->addr(),1)-1;
    ymin=(*ybasic)[i_basic];
  }
  if (basic_number<N) {
    j_non_basic=F77NAME(idmin)(N-basic_number,r->addr(),1)-1;
    rmin=(*r)[j_non_basic];
    j_non_basic+=basic_number;
  }
#ifdef DEBUG
//cout << "\ti_basic,ymin = " << i_basic << " " << ymin << endl;
//cout << "\tj_non_basic,rmin = " << j_non_basic << " " << rmin << endl;
#endif

  if (min(ymin,rmin)>=this->zero_) return this->current_status=OPTIMAL;
  if (ymin<rmin) {
//  TRACER_CALL(t,"LinearProgram::primalSimplexStep ymin<rmin");
    Vector<double,double> axis(basic_number,double_zero_);
    axis[i_basic]=-this->one_;
    Vector<double,double> gb(basic_number,this->zero_);
    Vector<double,double> resid(basic_number);
    QR->solveOverdetermined(axis,gb,resid);
    double gb_norm=gb.nrm2()*static_cast<double>(basic_number);
    double basic_ratio=double_undefined_;
    int j_basic=basicPrimalPivot(gb,basic_ratio);

    double non_basic_ratio=double_undefined_;
    int i_non_basic=-1;
    if (basic_number<M) {
      Vector<double,double> gn(M-basic_number);
      for (int i=basic_number;i<M;i++) {
        double gni=this->zero_;
        int rowi=row[i];
        double Ai_norm=double_zero_;
        for (int j=0;j<basic_number;j++) {
          double Aij=(*A)(rowi,column[j]);
          gni+=Aij*gb[j];
          Ai_norm+=Aij*Aij;
        }
        if (abs(gni)<gb_norm*sqrt(Ai_norm)*DBL_EPSILON) gni=double_zero_;
        gn[i-basic_number]=gni;
      }
      i_non_basic=nonBasicPrimalPivot(gn,non_basic_ratio);
    }
#ifdef DEBUG
//  cout << "i_basic = " << i_basic << endl;
//  cout << "basic_ratio,non_basic_ratio = " << basic_ratio << " "
//       << non_basic_ratio << endl;
//  cout << "j_basic = " << j_basic << endl;
//  cout << "i_non_basic = " << i_non_basic << endl;
#endif

    if (basic_ratio<non_basic_ratio) drop(i_basic,j_basic);
    else if (i_non_basic<M) switchRows(i_basic,i_non_basic);
    else this->current_status=UNBOUNDED;
  } else {
//  TRACER_CALL(t,"LinearProgram::primalSimplexStep ymin>=rmin");
    int coljn=column[j_non_basic];
    Vector<double,double> *hb=0;
    int j_basic=-1,i_non_basic=-1;
    double basic_ratio=double_undefined_;
    double non_basic_ratio=double_undefined_;
    double hb_norm=double_zero_;
    if (basic_number>0) {
      Vector<double,double> rhs(basic_number);
      for (int i=0;i<basic_number;i++) {
        rhs[i]=(*A)(row[i],coljn);
      }
      hb=OPERATOR_NEW Vector<double,double>(basic_number);
      Vector<double,double> resid(basic_number);
      QR->solveOverdetermined(rhs,*hb,resid);
      j_basic=basicPrimalPivot(*hb,basic_ratio);
      hb_norm=hb->nrm2()*static_cast<double>(basic_number);
    }
    if (basic_number<M) {
      Vector<double,double> hn(M-basic_number);
      for (int i=basic_number;i<M;i++) {
        int rowi=row[i];
        double hni=-(*this->A)(rowi,coljn);
        double Ai_norm=double_zero_;
        for (int j=0;j<basic_number;j++) {
          double Aij=(*A)(rowi,column[j]);
          hni+=Aij*(*hb)[j];
          Ai_norm+=Aij*Aij;
        }
        if (abs(hni)<(abs((*A)(rowi,coljn))+hb_norm*Ai_norm)*DBL_EPSILON){
          hni=double_zero_;
        }
        hn[i-basic_number]=hni;
      }
      i_non_basic=nonBasicPrimalPivot(hn,non_basic_ratio);
    }
    if (hb) delete hb; hb=0;
#ifdef DEBUG
//  cout << "basic_ratio,non_basic_ratio = " << basic_ratio << " "
//       << non_basic_ratio << endl;
//  cout << "j_non_basic = " << j_non_basic << endl;
//  cout << "j_basic = " << j_basic << endl;
//  cout << "i_non_basic = " << i_non_basic << endl;
#endif

    if (basic_ratio<non_basic_ratio) switchColumns(j_basic,j_non_basic);
    else if (i_non_basic<M) add(i_non_basic,j_non_basic);
    else this->current_status=UNBOUNDED;
  }
  computeSolution();
  return this->current_status;
}
*/

/*
template<> STATUS_OPTION LinearProgram<double>::dualSimplexStep() {
//TRACER_CALL(t,"LinearProgram::dualSimplexStep");
  if (current_status==OPTIMAL) return OPTIMAL;
  CHECK_TEST(current_status==DUAL_FEASIBLE);
#ifdef DEBUG
//printOn(cout);
#endif
  int M=A->size(0),N=A->size(1);
//find largest feasible simplex step
  double smin=double_undefined_,xmin=double_undefined_;
  int i_non_basic=-1,j_basic=-1;
  if (basic_number<M) {
    i_non_basic=F77NAME(idmin)(M-basic_number,s->addr(),1)-1;
    smin=(*s)[i_non_basic];
    i_non_basic+=basic_number;
  }
  if (basic_number>0) {
    j_basic=F77NAME(idmin)(basic_number,xbasic->addr(),1)-1;
    xmin=(*xbasic)[j_basic];
  }
#ifdef DEBUG
//cout << "\tj_basic,xmin = " << j_basic << " " << xmin << endl;
//cout << "\ti_non_basic,smin = " << i_non_basic << " " << smin << endl;
#endif
  if (min(xmin,smin)>=double_zero_) return current_status=OPTIMAL;

  if (xmin<smin) {
//  TRACER_CALL(t,"LinearProgram::dualSimplexStep xmin<smin");
    Vector<double,double> axis(basic_number,double_zero_);
    axis[j_basic]=double_one_;
    Vector<double,double> gb(basic_number);
    QR->solveUnderdetermined(axis,gb);
    double gb_norm=gb.nrm2()*static_cast<double>(basic_number);
    double basic_ratio=double_undefined_;
    int i_basic=basicDualPivot(gb,basic_ratio);

    double non_basic_ratio=double_undefined_;
    int j_non_basic=-1;
    if (basic_number<N) {
      Vector<double,double> gn(N-basic_number);
      for (int j=basic_number;j<N;j++) {
        double gnj=double_zero_;
        int colj=column[j];
        double Aj_norm=double_zero_;
        for (int i=0;i<basic_number;i++) {
          double Aij=(*A)(row[i],colj);
          gnj+=(gb)[i]*Aij;
          Aj_norm+=Aij*Aij;
        }
        if (abs(gnj)<gb_norm*sqrt(Aj_norm)*DBL_EPSILON) gnj=double_zero_;
        gn[j-basic_number]=-gnj;
      }
      j_non_basic=nonBasicDualPivot(gn,non_basic_ratio);
    }
#ifdef DEBUG
//  cout << "j_basic = " << j_basic << endl;
//  cout << "basic_ratio,non_basic_ratio = " << basic_ratio << " "
//       << non_basic_ratio << endl;
//  cout << "i_basic = " << i_basic << endl;
//  cout << "j_non_basic = " << j_non_basic << endl;
#endif

    if (basic_ratio<non_basic_ratio) drop(i_basic,j_basic);
    else if (j_non_basic<N) switchColumns(j_basic,j_non_basic);
    else current_status=UNBOUNDED;
  } else {
//  TRACER_CALL(t,"LinearProgram::dualSimplexStep xmin>=smin");
    double basic_ratio=double_undefined_;
    double non_basic_ratio=double_undefined_;
    int i_basic=-1,j_non_basic=-1;
    int rowin=row[i_non_basic];
    Vector<double,double> *hb=0;
    double hb_norm=double_zero_;
    if (basic_number>0) {
      Vector<double,double> rhs(basic_number);
      for (int j=0;j<basic_number;j++) {
        rhs[j]=(*A)(rowin,column[j]);
      }
      hb=OPERATOR_NEW Vector<double,double>(basic_number);
      QR->solveUnderdetermined(rhs,*hb);
      i_basic=basicDualPivot(*hb,basic_ratio);
      hb_norm=hb->nrm2()*static_cast<double>(basic_number);
    }
    if (basic_number<N) {
      Vector<double,double> hn(N-basic_number);
      for (int j=basic_number;j<N;j++) {
        int colj=column[j];
        double hnj=(*A)(rowin,colj);
        double Aj_norm=double_zero_;
        for (int i=0;i<basic_number;i++) {
          double Aij=(*A)(row[i],colj);
          hnj-=(*hb)[i]*Aij;
          Aj_norm+=Aij*Aij;
        }
        if (abs(hnj)<(abs((*A)(rowin,colj))+hb_norm*Aj_norm)*DBL_EPSILON){
          hnj=double_zero_;
        }
        hn[j-basic_number]=hnj;
      }
      j_non_basic=nonBasicDualPivot(hn,non_basic_ratio);
    }
    if (hb) delete hb; hb=0;
#ifdef DEBUG
//  cout << "basic_ratio,non_basic_ratio = " << basic_ratio << " "
//       << non_basic_ratio << endl;
//  cout << "i_non_basic in permuted order = " << i_non_basic << endl;
//  cout << "i_basic in permuted order = " << i_basic << endl;
//  cout << "j_non_basic in permuted order = " << j_non_basic << endl;
#endif
    if (basic_ratio<non_basic_ratio) switchRows(i_basic,i_non_basic);
    else if (j_non_basic<N) add(i_non_basic,j_non_basic);
    else current_status=UNBOUNDED;
  }
  computeSolution();
  return current_status;
}
*/

/*
template<> STATUS_OPTION 
LinearProgram<double>::findPrimalBasicFeasibleGuess() {
//TRACER_CALL(t,"LinearProgram::findPrimalBasicFeasibleGuess");
//if we have already dinked with this program, get out now
  if (current_status==OPTIMAL || current_status==PRIMAL_FEASIBLE) {
    return current_status;
  }
  CHECK_TEST(current_status!=INFEASIBLE && current_status!=UNBOUNDED);

//introduce artificial variables for positive b's
#ifdef DEBUG
//cout << "introducing artificial variables" << endl;
#endif
  int m=A->size(0),n=A->size(1);
  int ib=F77NAME(idsumn)(m,b->addr(),1);
  Matrix<double,double> *temp_A=
    OPERATOR_NEW Matrix<double,double>(m,n+m-ib,double_zero_);
  temp_A->copyFrom('A',m,n,*A);
  Vector<double,double> *temp_c=
    OPERATOR_NEW Vector<double,double>(n+m-ib,double_zero_);
  int pos_b=n;
  for (int i=0;i<m;i++) {
    if ((*b)[i]>double_zero_) {
      (*temp_A)(i,pos_b)=double_one_;
      (*temp_c)[pos_b]=double_one_;
      pos_b++;
    }
  }
  LinearProgram LP(*temp_A,*b,*temp_c);
//int *basic_row_index=OPERATOR_NEW_BRACKET(int,m-ib);
//int ip=0;
//for (int i=0;i<m;i++) {
//  if ((*b)[i]>double_zero_) basic_row_index[ip++]=i;
//}
//int *basic_column_index=OPERATOR_NEW_BRACKET(int,m-ib);
//for (int j=n;j<n+m-ib;j++) basic_column_index[j-n]=j;
//LP.specifyBasicVariables(m-ib,basic_row_index,basic_column_index);

  while (LP.current_status==DUAL_FEASIBLE) LP.simplexStep();
#ifdef DEBUG
//cout << "after solving LP to find initial feasible guess" << endl;
//cout << "LP:" << endl;
//LP.printOn(cout);
#endif
  CHECK_TEST(LP.current_status==OPTIMAL);

  basic_number=LP.basic_number;
  int jj=0;
  for (int j=0;j<n+m-ib;j++) {
    int colj=LP.column[j];
    if (colj<n) column[jj++]=colj;
  }
  memcpy(row,LP.row,m*sizeof(int));
  delete QR; QR=0;
  QR=OPERATOR_NEW GramSchmidtQRFactorization<double,double>(*LP.QR);

  computeSolution();
  current_status=PRIMAL_FEASIBLE;
#ifdef DEBUG
//cout << "after setting initial feasible guess" << endl;
//printOn(cout);
#endif
  delete temp_A;
  delete temp_c;
  return current_status;
}
*/

/*
template<> STATUS_OPTION 
LinearProgram<double>::findDualBasicFeasibleGuess() {
//TRACER_CALL(t,"LinearProgram::findDualBasicFeasibleGuess");
#ifdef DEBUG
//cout << "current_status = " << current_status << endl;
#endif
  if (current_status==OPTIMAL || current_status==DUAL_FEASIBLE) {
    return current_status;
  }
  CHECK_TEST(current_status==UNKNOWN);


#ifdef DEBUG
//printOn(cout);
//cout << "\n\n\tintroducing artificial variables" << endl;
#endif
  int m=A->size(0),n=A->size(1);
  int ic=F77NAME(idsump)(n,c->addr(),1);
  Matrix<double,double> *temp_A=
    OPERATOR_NEW Matrix<double,double>(m+n-ic,n,double_zero_);
  temp_A->copyFrom('A',m,n,*A);
  Vector<double,double> *temp_b=
    OPERATOR_NEW Vector<double,double>(m+n-ic,double_zero_);
  int neg_c=m;
  for (int j=0;j<n;j++) {
    if ((*c)[j]<double_zero_) {
      (*temp_A)(neg_c,j)=-double_one_;
      (*temp_b)[neg_c]=-double_one_;
      neg_c++;
    }
  }
  LinearProgram LP(*temp_A,*temp_b,*c);
#ifdef DEBUG
//LP.printOn(cout);
#endif
//int *basic_column_index=OPERATOR_NEW_BRACKET(int,n-ic);
//int jn=0;
//for (int j=0;j<n;j++) {
//  if ((*c)[j]<double_zero_) basic_column_index[jn++]=j;
//}
//int *basic_row_index=OPERATOR_NEW_BRACKET(int,n-ic);
//for (int i=m;i<m+n-ic;i++) basic_row_index[i-m]=i;
//{
//TRACER_CALL(t,"LinearProgram::findDualBasicFeasibleGuess");
//LP.specifyBasicVariables(n-ic,basic_row_index,basic_column_index);
#ifdef DEBUG
//LP.printOn(cout);
#endif
//}
//delete basic_column_index; basic_column_index=0;
//delete basic_row_index; basic_row_index=0;

  while (LP.current_status==PRIMAL_FEASIBLE) LP.simplexStep();
#ifdef DEBUG
//cout << "after solving LP to find initial feasible guess" << endl;
//cout << "LP:" << endl;
//LP.printOn(cout);
#endif
  CHECK_TEST(LP.current_status==OPTIMAL);

  basic_number=LP.basic_number;
  int ii=0;
  for (int i=0;i<m+n-ic;i++) {
    int rowi=LP.row[i];
    if (rowi<m) row[ii++]=rowi;
  }
  memcpy(column,LP.column,n*sizeof(int));
  delete QR; QR=0;
  QR=OPERATOR_NEW GramSchmidtQRFactorization<double,double>(*LP.QR);

  computeSolution();
  current_status=DUAL_FEASIBLE;
#ifdef DEBUG
//cout << "after setting initial feasible guess" << endl;
//printOn(cout);
#endif
  delete temp_A;
  delete temp_b;
  return current_status;
}
*/

/*
template class LinearProgram<double>;
//template LinearProgram<double>::LinearProgram(const Matrix<double,double>&,const Matrix<double,double>&,const Matrix<double,double>&);
//template double LinearProgram<double>::currentValue() const;
//template status_option LinearProgram<double>::simplexStep();
//template void LinearProgram<double>::costBounds(Matrix<double,double>&,
//   Matrix<double,double>&) const;
//template void LinearProgram<double>::constraintBounds(Matrix<double,double>&,
//   Matrix<double,double>&) const;
//template void LinearProgram<double>::printOn(ostream&) const;
*/

/*
template<> STATUS_OPTION LinearProgram<double>::findBasicFeasibleGuess() {
//TRACER_CALL(t,"LinearProgram::findBasicFeasibleGuess");
//if we have already dinked with this program, get out now
  CHECK_TEST(current_status!=INFEASIBLE && current_status!=UNBOUNDED);
  if (current_status!=UNKNOWN) return current_status;
  int m=A->size(0),n=A->size(1);
  int ib=F77NAME(idsumn)(m,b->addr(),1);
  int ic=F77NAME(idsump)(n,c->addr(),1);
#ifdef DEBUG
//cout << "ib,ic = " << ib << " " << ic << endl;
#endif
  if (m-ib<=n-ic) return findPrimalBasicFeasibleGuess();
  else return findDualBasicFeasibleGuess();
}
*/
