#include <complex>
#include <float.h>
#include <limits>
#include <math.h>
#include "Vector.H"
double double_mone_=-1.;
double double_one_=1.;
double double_zero_=0.;
double double_undefined_=numeric_limits<double>::infinity();
complex<double> complex_double_mone_=complex<double>(-1.,0.);
complex<double> complex_double_one_=complex<double>(1.,0.);
complex<double> complex_double_zero_=complex<double>(0.,0.);
complex<double> complex_double_undefined_=
  complex<double>(double_undefined_,double_undefined_);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  void F77NAME(daxpy)(const int &n,const double &a,const double *x,
    const int &incx,double *y,const int &incy);
  void F77NAME(dlabad)(double &small,double &large);
  double F77NAME(dlamch)(const char&);
  double F77NAME(dzasum)(const int&,const complex<double>*,const int&);
  double F77NAME(dznrm2)(const int &n,complex<double> *x,
    const int &incx);
  int F77NAME(ilaenv)(const int &ispec,const char *name,const char *opts,
    const int &n1,const int &n2,const int &n3,const int &n4);
  int F77NAME(ilatrans)(const char &trans);
  int F77NAME(ilauplo)(const char &trans);
  void F77NAME(zaxpy)(const int &n,const complex<double> &a,
    const complex<double> *x,const int &incx,complex<double> *y,
    const int &incy);
  void F77NAME(zcopy)(const int &n,const complex<double> *sx,
    const int &incx,complex<double> *sy,const int &incy);
  complex<double> F77NAME(zdotc)(const int &n,const complex<double> *sx,
    const int &incx,const complex<double> *sy,const int &incy);
  complex<double> F77NAME(zdotu)(const int &n,const complex<double> *sx,
    const int &incx,const complex<double> *sy,const int &incy);
  void F77NAME(zdrscl)(const int &n,const double &sa,complex<double> *sx,
    const int &incx);
  void F77NAME(zgbcon)(const char &norm,const int &n,const int &kl,
    const int &ku,const complex<double> *AB,const int &ldab,int *ipiv,
    const double &anorm,double &rcond,complex<double> *work,
    double *rwork,int &info);
  void F77NAME(zgbconnp)(const char &norm,const int &n,const int &kl,
    const int &ku,const complex<double> *AB,const int &ldab,
    const double &anorm,double &rcond,complex<double> *work,
    double *rwork,int &info);
  void F77NAME(zgbequ)(const int &m,const int &n,const int &kl,
    const int &ku,const complex<double> *AB,const int &ldab,double *r,
    double *c,double &rowcnd,double &colcnd,double &amax,int &info);
  void F77NAME(zgbamv)(const char &trans,const int &m,const int &n,
    const int &kl,const int &ku,const double &alpha,
    const complex<double> *A,const int &lda,const complex<double> *x,
    const int &incx,const double &beta,double *y,
    const int &incy);
  void F77NAME(zgbmv)(const char &trans,const int &m,const int &n,
    const int &kl,const int &ku,const complex<double> &alpha,
    const complex<double> *A,const int &lda,const complex<double> *x,
    const int &incx,const complex<double> &beta,complex<double> *y,
    const int &incy);
  void F77NAME(zgbtf2np)(const int &m,const int &n,const int &kl,
    const int &ku,complex<double> *AB,const int &ldab,int &info);
  void F77NAME(zgbtrf)(const int &m,const int &n,const int &kl,
    const int &ku,complex<double> *AB,const int &ldab,int *ipiv,
    int &info);
  void F77NAME(zgbtrs)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const complex<double> *AB,
    const int &ldab,const int *ipiv,complex<double> *b,const int &ldb,
    int &info);
  void F77NAME(zgbtrsnp)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const complex<double> *AB,
    const int &ldab,complex<double> *b,const int &ldb,int &info);
  void F77NAME(zdscal)(const int &n,const double &a,
    complex<double> *x,const int &incx);
  void F77NAME(zgeamv)(const int &trans,const int &m,const int &n,
    double &alpha,const complex<double> *A,const int &lda,
    const complex<double> *x,const int &incx,const double &beta,
    double *y,const int &incy);
  void F77NAME(zgecon)(const char &norm,const int &n,
    const complex<double> *A,const int &lda,const double &anorm,
    double &rcond,complex<double> *work,double *rwork,int &info);
  void F77NAME(zgeequ)(const int &m,const int &n,
    const complex<double> *A,const int &lda,double *r,double *c,
    double &rowcnd,double &rolcnd,double &amax,int &info);
  void F77NAME(zgeev)(const char &jobvl,const char &jobvr,const int &n,
    complex<double> *A,const int &lda,complex<double> *w,
    complex<double> *vl,const int &ldvl,complex<double> *vr,
    const int &ldvr,complex<double> *work,const int &lwork,
    double *rwork,int &info);
  void F77NAME(zgels)(const char &trans,const int &m,const int &n,
    const int &nrhs,complex<double> *a,const int &lda,complex<double> *b,
    const int &ldb,complex<double> *work,const int &lwork,
    const int &info);
  void F77NAME(zgemm)(const char &transa,const char &transb,
    const int &m,const int &n,const int &k, 
    const complex<double> &alpha,const complex<double> *A,
    const int &lda,const complex<double> *B,const int &ldb,
    const complex<double> &beta,complex<double> *C,const int &ldc);
  void F77NAME(zgemv)(const char &trans,const int &m,const int &n,
    const complex<double> &alpha,const complex<double> *a,
    const int &lda,const complex<double> *x,const int &incx,
    const complex<double> &beta,complex<double> *y,const int &incy);
  void F77NAME(zgelqf)(const int &m,const int &n,complex<double> *A,
    const int &lda,complex<double> *tau,complex<double> *work,
    int &lwork,int &info);
  void F77NAME(zgeqp3)(const int &m,const int &n,complex<double> *A,
    const int &lda,int *jpvt,complex<double> *tau,complex<double> *work,
    const int &lwork,double *rwork,int &info);
  void F77NAME(zgeqrf)(const int &m,const int &n,complex<double> *A,
    const int &lda,complex<double> *tau,complex<double> *work,
    int &lwork,int &info);
  void F77NAME(zgerc)(const int&,const int&,const complex<double>&,
    const complex<double>*,const int&,const complex<double>*,const int&,
    complex<double>*,const int&);
  void F77NAME(zgeru)(const int&,const int&,const complex<double>&,
    const complex<double>*,const int&,const complex<double>*,const int&,
    complex<double>*,const int&);
  void F77NAME(zgesv)(int &n,int &nrhs,complex<double> *a,int &lda,
    int *ipiv,complex<double> *b, int &ldb, int &info);
  void F77NAME(zgesvd)(const char &jobu,const char &jobvt,const int &m, 
    const int &n,const complex<double> *A,const int &lda,
    double *S,complex<double> *U,const int &ldu,complex<double> *Vt,
    const int &ldvt,complex<double> *work,const int &lwork,
    double *rwork,int &info);
  void F77NAME(zgesc2)(const int &n,const complex<double> *A,
    const int &lda,complex<double> *rhs,const int *ipiv,const int *jpiv,
    complex<double> &scale);
  void F77NAME(zgetc2)(const int &n,complex<double> *A,const int &lda,
    int *ipiv,int *jpiv,int &info);
  void F77NAME(zgetrf)(int &m,int &n,complex<double> *a,int &lda,
    int *ipiv,int &info);
  void F77NAME(zgetrs)(const char &trans,const int &n,const int &nrhs,
    complex<double> *a,const int &lda,int *ipiv,complex<double> *b,
    const int &ldb,int &info);
  void F77NAME(zgetri)(int &n,complex<double> *a,int &lda,int *ipiv,
    complex<double> *work,int &lwork,int &info);
  void F77NAME(zgtcon)(const char &norm,const int &n,
    const complex<double> *L,const complex<double> *D,
    const complex<double> *U,const complex<double> *U2,
    const int *ipiv,const double &anorm,double &rcond,
    complex<double> *work,int &info);
  void F77NAME(zgtconnp)(const char &norm,const int &n,
    const complex<double> *L,const complex<double> *D,
    const complex<double> *U,const double &anorm,double &rcond,
    complex<double> *work,int &info);
  void F77NAME(zgtmv)(const int &n,const complex<double> &alpha,
    const complex<double> *L,const complex<double> *D,
    const complex<double> *U,const complex<double> *x,const int &incx,
    const complex<double> &beta,complex<double> *y,const int &incy);
  void F77NAME(zgtrfs)(const char &trans,const int &n,const int &nrhs,
    const complex<double> *dl,const complex<double> *d,
    const complex<double> *du,const complex<double> *dlf,
    const complex<double> *df,const complex<double> *duf,
    const complex<double> *du2,const int *ipiv,const complex<double> *B,
    const int &ldb,complex<double> *X,const int &ldx,double *ferr,
    double *berr,complex<double> *work,double *rwork,int &info);
  void F77NAME(zgtrfsnp)(const char &trans,const int &n,const int &nrhs,
    const complex<double> *dl,const complex<double> *d,
    const complex<double> *du,const complex<double> *dlf,
    const complex<double> *df,const complex<double> *duf,
    const complex<double> *B,const int &ldb,complex<double> *X,
    const int &ldx,double *ferr,double *berr,complex<double> *work,
    double *rwork,int &info);
  void F77NAME(zgtrfsr)(const char &trans,const int &n,const int &nrhs,
    const complex<double> *dl,const complex<double> *d,
    const complex<double> *du,const complex<double> *dlf,
    const complex<double> *df,const complex<double> *duf,
    const complex<double> *du2,const int *ipiv,const complex<double> *B,
    const int &ldb,complex<double> *X,const int &ldx,double *ferr,
    double *berr,complex<double> *work,double *rwork,int &info);
  void F77NAME(zgtrfsrnp)(const char &trans,const int &n,const int &nrhs,
    const complex<double> *dl,const complex<double> *d,
    const complex<double> *du,const complex<double> *dlf,
    const complex<double> *df,const complex<double> *duf,
    const complex<double> *B,const int &ldb,complex<double> *X,
    const int &ldx,double *ferr,double *berr,complex<double> *work,
    double *rwork,int &info);
  void F77NAME(zhtmv)(const int &n,const complex<double> &alpha,
    const complex<double> *L,const double *D,const complex<double> *x,
    const int &incx,const complex<double> &beta,complex<double> *y,
    const int &incy);
  void F77NAME(zsteqr)(const char &compz,const int &n,double *D,
    complex<double> *E,complex<double> *Z,const int &ldz,double *work,
    int &info);
  void F77NAME(zgtsv)(const int &n,const int &nrhs,complex<double> *dl,
    complex<double> *d,complex<double> *du,complex<double> *b,
    const int &ldb,int &info);
  void F77NAME(zgtsvnp)(const int &n,const int &nrhs,complex<double> *dl,
    complex<double> *d,complex<double> *du,complex<double> *b,
    const int &ldb,int &info);
  void F77NAME(zgtsvr)(const int &n,const int &nrhs,complex<double> *dl,
    complex<double> *d,complex<double> *du,complex<double> *b,
    const int &ldb,int &info);
  void F77NAME(zgtsvrnp)(const int &n,const int &nrhs,complex<double> *dl,
    complex<double> *d,complex<double> *du,complex<double> *b,
    const int &ldb,int &info);
  void F77NAME(zgttrf)(const int &n,complex<double> *L,
    complex<double> *D,complex<double> *U,complex<double> *U2,int *ipiv,
    int &info);
  void F77NAME(zgttrfnp)(const int &n,complex<double> *L,
    complex<double> *D,complex<double> *U,int &info);
  void F77NAME(zgttrs)(const char &trans,const int &n,const int &nrhs,
    const complex<double> *L,const complex<double> *D,
    const complex<double> *U,const complex<double> *U2,
    const int *ipiv,complex<double> *B,const int &ldb,int &info);
  void F77NAME(zhbev)(const char &jobz,const char &uplo,const int &n,
    const int &kd,complex<double> *AB,const int &ldab,double *w,
    complex<double> *Z,const int &ldz,complex<double> *work,double *rwork,
    int &info);
  void F77NAME(zhbamv)(const char &uplo,const int &n,const int &k,
    const double &alpha,const complex<double> *A,const int &lda,
    const complex<double> *x,const int &incx,const double &beta,
    double *y,const int &incy);
  void F77NAME(zhbmv)(const char &uplo,const int &n,const int &k,
    const complex<double> &alpha,const complex<double> *A,const int &lda,
    const complex<double> *x,const int &incx,const complex<double> &beta,
    complex<double> *y,const int &incy);
  void F77NAME(zhecon)(const char &uplo,const int &n,
    const complex<double> *A,const int &lda,const int *ipiv,
    const double &anorm,double &rcond,complex<double> *work,int &info);
  void F77NAME(zheequb)(const char &uplo,const int &n,
    const complex<double> *A,const int &lda,double *s,double &scond,
    double &amax,complex<double> *work,int &info);
  void F77NAME(zheev)(const char &jobz,const char &uplo,const int &n,
    complex<double> *A,const int &lda,double *w,complex<double> *work,
    const int &lwork,double *rwork,int &info);
  void F77NAME(zhemm)(const char &side,const char &uplo,const int &m,
    const int &n,const complex<double> &alpha,const complex<double> *A,
    const int &lda,const complex<double> *B,const int &ldb,
    const complex<double> &beta,complex<double> *C,const int &ldc);
  void F77NAME(zhemv)(const char &uplo,const int &n,
    const complex<double> &alpha,const complex<double> *A,
    const int &lda,const complex<double> *x,const int &incx,
    const complex<double> &beta,complex<double> *y,const int &incy);
  void F77NAME(zher)(const char &uplo,const int &n,
    const complex<double> &alpha,const complex<double> *X,
    const int &incx,complex<double> *A,const int &lda);
  void F77NAME(zher2)(const char &uplo,const int &n,
    const complex<double> &alpha,const complex<double> *X,
    const int &incx,const complex<double> *y,const int &ldy,
    complex<double> *A,const int &lda);
  void F77NAME(zher2k)(const char &uplo,const char &trans,const int &n,
    const int &k,const complex<double> &alpha,const complex<double> *A,
    const int &lda,const complex<double> *B,const int &ldb,
    const complex<double> &beta,complex<double> *C,const int &ldc);
  void F77NAME(zherk)(const char &uplo,const char &trans,const int &n,
    const int &k,const complex<double> &alpha,const complex<double> *A,
    const int &lda,const complex<double> &beta,complex<double> *C,
    const int &ldc);
  void F77NAME(zhetrf)(const char &uplo,const int &n,complex<double> *a,
    const int &lda,int *ipiv,complex<double> *work,const int &lwork,
    int &info);
  void F77NAME(zhetri)(const char &uplo,const int &n,complex<double> *a,
    const int &lda,int *ipiv,complex<double> *work,int &info);
  void F77NAME(zhetrs)(const char &uplo,const int &n,const int &nrhs,
    complex<double> *a,const int &lda,int *ipiv,complex<double> *b,
    const int &ldb,int &info);
  void F77NAME(zhseqr)(const char &job,const char &compz,const int &n,
    const int &ilo,const int &ihi,complex<double> *H,const int &ldh,
    complex<double> *w,complex<double> *Z,const int &ldz,
    complex<double> *work,const int &lwork,int &info);
  void F77NAME(zlacn2)(const int &n,complex<double> *v,
    complex<double> *x,double &est,int &kase,int *isave);
  void F77NAME(zlacpy)(const char &uplo,const int &m,const int &n,
    const complex<double> *A,const int &lda,complex<double> *B,
    const int &ldb);
  void F77NAME(zladiv)(const complex<double>&,const complex<double>&,
    const complex<double>&,const complex<double>&,complex<double>&,
    complex<double>&);
  void F77NAME(zlaic1)(const int &job,const int &j,
    const complex<double> *x,const double &sest,const complex<double> *w,
    const complex<double> &gamma,double &sestpr,complex<double> &s,
    complex<double> &c);
  double F77NAME(zlangb)(const char &norm,const int &n,const int &kl,
    const int &ku,const complex<double> *AB,const int &ldab,
    double *work);
  double F77NAME(zlange)(const char &norm,const int &m,const int &n,
    const complex<double> *A,const int &lda,double *work);
  double F77NAME(zlanhb)(const char &norm,const char &uplo,const int &n,
    const int &k,const complex<double> *AB,const int &ldab,double *work);
  double F77NAME(zlanhs)(const char &norm,const int &n,
    const complex<double> *A,const int &lda,double *work);
  double F77NAME(zlangt)(const char &norm,const int &n,
    const complex<double> *L,const complex<double> *D,
    const complex<double> *U);
  double F77NAME(zlanhe)(const char &norm,const char &uplo,const int &n,
    const complex<double> *A,const int &lda,double *work);
  double F77NAME(zlanht)(const char &norm,const int &n,
    const double *D,const complex<double> *E);
  double F77NAME(zlansb)(const char &norm,const char &uplo,
    const int &n,const int &k,const complex<double> *AB,const int &ldab,
    double *work);
  double F77NAME(zlantb)(const char &norm,const char &uplo,
    const char &diag,const int &n,const int &k,const complex<double> *AB,
    const int &ldab,double *work);
  double F77NAME(zlantr)(const char &norm,const char &uplo,
    const char &diag,const int &m,const int &n,const complex<double> *A,
    const int &lda,double *work);
  void F77NAME(zlaqgb)(const int &m,const int &n,const int &kl,
    const int &ku,complex<double> *AB,const int &ldab,double *r,
    double *c,const double &rowcnd,const double &colcnd,
    const double &amax,char &equed);
  void F77NAME(zlaqge)(const int &m,const int &n,complex<double> *A,
    const int &lda,const double *r,const double *c,const double &rowcnd,
    const double &colcnd,const double &amax,char &equed);
  void F77NAME(zlaqhe)(const char &uplo,const int &n,complex<double> *A,
    const int &lda,const double *s,const double &scond,const double &amax,
    char &equed);
  void F77NAME(zlaqsb)(const char &uplo,const int &n,const int &kd,
    complex<double> *AB,const int &ldab,const double *s,
    const double &scond,const double &amax,char &equed);
  void F77NAME(zlaqsy)(const char &uplo,const int &n,complex<double> *A,
    const int &lda,double *s,const double &scond,const double &amax,
    char &equed);
  void F77NAME(zlargv)(const int &n,complex<double> *x,const int &incx,
    complex<double> *y,const int &incy,double *c,
    const int &incc);
  void F77NAME(zlartv)(const int &n,complex<double> *x,const int &incx,
    complex<double> *y,const int &incy,const double *c,
    const complex<double> *s,const int &incc);
  void F77NAME(zlascl)(const char &type,const int &kl,const int &ku,
    const double &cfrom,const double &cto,const int &m,const int &n,
    complex<double> *A,const int &lda,int &info);
  void F77NAME(zlaset)(const char &uplo,const int &m,const int &n,
    const complex<double> &alpha,const complex<double> &beta,
    complex<double> *A,const int &lda);
  void F77NAME(zlaswp)(const int &n,complex<double> *A,const int &lda,
    const int &k1,const int &k2,const int *ipiv,const int &incx);
  void F77_NAME(zla_heamv)(const int &uplo,const int &n,
    const double &alpha,const complex<double> *A,const int &lda,
    const complex<double> *x,const int &incx,const double &beta,
    double *y,const int &incy);
  double F77_NAME(zla_porpvgrw)(const char &uplo,const int &ncols,
    const complex<double> *A,const int &lda,const complex<double> *AF,
    const int &ldaf,complex<double> *work);
  double F77_NAME(zla_syrpvgrw)(const char &uplo,const int &n,
    const int &info,const complex<double> *A,const int &lda,
    const complex<double> *AF,const int &ldaf,const int *ipiv,
    complex<double> *work);
  void F77NAME(zorgqr)(int &m,int &n,int &k,complex<double> *a,int &lda,
    const complex<double> *tau,complex<double> *work,int &lwork,
    int &info);
  void F77NAME(zpbcon)(const char &uplo,const int &n,const int &kd,
    const complex<double> *AB,const int &ldab,const double &anorm,
    double &rcond,complex<double> *work,double *rwork,int &info);
  void F77NAME(zpbequ)(const char &uplo,const int &n,const int &kd,
    const complex<double> *AB,const int &ldab,double *s,double &scond,
    double &amax,int &info);
  void F77NAME(zpbsv)(const char &uplo,const int &n,const int &kd,
    const int &nrhs,complex<double> *AB,const int &ldab,
    complex<double> *B,const int &ldb,int &info);
  void F77NAME(zpbtrf)(const char &uplo,const int &n,const int &kd,
    complex<double> *AB,const int &ldab,int &info);
  void F77NAME(zpbtrs)(const char &uplo,int &n,const int &kd,
    const int &nrhs,const complex<double> *AB,const int &ldab,
    complex<double> *b,const int &ldb,int &info);
  void F77NAME(zpocon)(const char &uplo,const int &n,
    const complex<double> *a,const int &lda,const complex<double> &anorm,
    double &rcond,complex<double> *work,double *rwork,int &info);
  void F77NAME(zpoequb)(const int &n,const complex<double> *A,
    const int &lda,double *s,double &scond,double &amax,int &info);
  void F77NAME(zpotrf)(const char &uplo,const int &n,complex<double> *A,
    const int &lda,int &info);
  void F77NAME(zpotri)(const char &uplo,const int &n,complex<double> *A,
    const int &lda,int &info);
  void F77NAME(zpotrs)(const char &uplo,const int &n,const int &nrhs,
    complex<double> *A,const int &lda,complex<double> *B,const int &ldb,
    int &info);
  void F77NAME(zptcon)(const int &n,const double *d,
    const complex<double> *e,const double &anorm,double &rcond,
    double *rwork,int &info);
  void F77NAME(zptrfs)(const char &uplo,const int &n,const int &nrhs,
    const double *d,const complex<double> *E,const double *df,
    const complex<double> *ef,const complex<double> *B,const int &ldb,
    complex<double> *X,const int &ldx,double *ferr,double *berr,
    complex<double> *work,double *rwork,int &info);
  void F77NAME(zptrfsr)(const char &uplo,const int &n,const int &nrhs,
    const double *d,const complex<double> *E,const double *df,
    const complex<double> *ef,const complex<double> *B,const int &ldb,
    complex<double> *X,const int &ldx,double *ferr,double *berr,
    complex<double> *work,double *rwork,int &info);
  void F77NAME(zptsv)(const int &n,const int &nrhs,double *d,
    complex<double> *e,complex<double> *b,const int &ldb,int &info);
  void F77NAME(zpttrf)(const int &n,double *d,complex<double> *e,
    int &info);
  void F77NAME(zpttrs)(const char &uplo,const int &n,const int &nrhs,
    double *d,complex<double> *e,complex<double> *B,const int &ldb,
    int &info);
  void F77NAME(zrot)(const int &n,complex<double> *sx,const int &incx,
    complex<double> *sy,const int &incy,const double &c,
    const complex<double> &s);
  void F77NAME(zrotg)(complex<double> &sa,complex<double> &sb,
    double &c,complex<double> &s);
  void F77NAME(zscal)(const int &n,const complex<double> &a,
    complex<double> *x,const int &incx);
  void F77NAME(zswap)(const int &n,complex<double> *sx,const int &incx,
    complex<double> *sy,const int &incy);
  void F77NAME(ztbsv)(const char &uplo,const char &trans,const char &diag,
    const int &n,const int &k,const complex<double> *A,const int &lda,
    complex<double> *x,const int &incx);
  void F77NAME(ztrcon)(const char &norm,const char &uplo,
    const char &diag,const int &n,const complex<double> *A,const int &lda,
    double &rcond,complex<double> *work,double *rwork,int &info);
  void F77NAME(ztrevc)(const char &side,const char &howmny,
    const bool *select,const int &n,complex<double> *T,const int &ldt,
    complex<double> *VL,const int &ldvl,complex<double> *VR,
    const int &ldvr,const int &mm,int &m,complex<double> *work,
    double *rwork,int &info);
  void F77NAME(ztrmv)(const char &uplo,const char &trans,
    const char &diag,const int &n,const complex<double> *A,const int &lda,
    complex<double> *x,const int &incx);
  void F77NAME(ztrsv)(const char&,const char&,const char&,const int &n,
    const complex<double> *A,const int &lda,complex<double> *x,
    const int &incx);
  void F77NAME(ztrmm)(const char &side,const char &uplo,
    const char &transa,const char &diag,const int &m,const int &n,
    const complex<double> &alpha,const complex<double> *A,const int &lda,
    complex<double> *B,const int &ldb);
  void F77NAME(ztrsm)(const char &side,const char &uplo,
    const char &transa,const char &diag,const int &m,const int &n,
    const complex<double> &alpha,const complex<double> *a,const int &lda,
    complex<double> *b,const int &ldb);
  void F77NAME(ztrtri)(const char &uplo,const char &diag,const int &n,
    complex<double> *A,const int &lda,int &info);
  void F77NAME(ztzrzf)(const int &m,const int &n,complex<double> *A,
    const int &lda,complex<double> *tau,complex<double> *work,
    const int &lwork,int &info);
  void F77NAME(zunglq)(const int &m,const int &n,const int &k,
    complex<double> *A,const int &lda,const complex<double> *tau,
    complex<double> *work,const int &lwork,int &info);
  void F77NAME(zungqr)(const int &m,const int &n,const int &k,
    complex<double> *A,const int &lda,const complex<double> *tau,
    complex<double> *work,const int &lwork,int &info);
  void F77NAME(zunmlq)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,complex<double> *A,const int &lda,
    complex<double> *tau,complex<double> *C,const int &ldc,
    complex<double> *work,int &lwork,int &info);
  void F77NAME(zunmqr)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,complex<double> *A,const int &lda,
    complex<double> *tau,complex<double> *C,const int &ldc,
    complex<double> *work,int &lwork,int &info);
  void F77NAME(zunmrz)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,const int &l,const complex<double> *A,
    const int &lda,const complex<double> *tau,complex<double> *C,
    const int &ldc,complex<double> *work,const int &lwork,int &info);
  int F77NAME(ilazlc)(const int &m,const int &n,const complex<double> *A,
    const int &lda);
  int F77NAME(ilazlr)(const int &m,const int &n,const complex<double> *A,
    const int &lda);
  int F77NAME(izamax)(const int &n,const complex<double> *x,
    const int &incx);
  int F77NAME(izmin)(const int &n,const complex<double> *a,
    const int &inca);
  int F77NAME(izsumn)(const int &n,const complex<double> *a,
    const int &inca);
  int F77NAME(izsump)(const int &n,const complex<double> *a,
    const int &inca);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Vector.C"
#ifdef DEBUG
template<> const complex<double>
  Vector<double,complex<double> >::undefined_(
  numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
#endif
template<> const complex<double>
  Vector<double,complex<double> >::zero_(double_zero_,double_zero_);
template<> const complex<double>
  Vector<double,complex<double> >::one_(double_one_,double_zero_);
template<> const complex<double>
  Vector<double,complex<double> >::mone_(double_mone_,double_zero_);

template<> int Vector<double,complex<double> >::amax() const { 
  return F77NAME(izamax)(sz,data,1)-1;
}
template<> double Vector<double,complex<double> >::asum() const {
  return F77NAME(dzasum)(sz,data,1);
}
template<> void Vector<double,complex<double> >::axpy(complex<double> a,
const Vector<double,complex<double> > &x) {
  F77NAME(zaxpy)(min(sz,x.sz),a,x.data,1,data,1);
}
template<> complex<double> Vector<double,complex<double> >::dot(
const Vector<double,complex<double> > &x) const {
  OBSOLETE(0);
}
template<> complex<double> Vector<double,complex<double> >::dotc(
const Vector<double,complex<double> > &x) const {
  return F77NAME(zdotc)(min(sz,x.sz),x.data,1,data,1);
}
template<> complex<double> Vector<double,complex<double> >::dotu(
const Vector<double,complex<double> > &x) const {
  return F77NAME(zdotu)(min(sz,x.sz),x.data,1,data,1);
}
template<> double Vector<double,complex<double> >::nrm2() const {
  return F77NAME(dznrm2)(sz,data,1);
}
template<> void Vector<double,complex<double> >::rot(
Vector<double,complex<double> > &x,double c,double s) {
  F77NAME(zrot)(min(sz,x.sz),x.data,1,data,1,c,s);
}
template<> void Vector<double,complex<double> >::scal(complex<double> a)
{
  F77NAME(zscal)(sz,a,data,1);
}
template<> void Vector<double,complex<double> >::swap(
Vector<double,complex<double> > &x) {
  F77NAME(zswap)(min(sz,x.sz),x.data,1,data,1);
}
template<> void Vector<double,complex<double> >::largv(
Vector<double,complex<double> > &y,Vector<double,double> &c) {
  int n=size();
  CHECK_SAME(n,y.size());
  CHECK_SAME(n,c.size());
  F77NAME(zlargv)(n,data,1,y.data,1,c.addr(),1);
}
template<> void Vector<double,complex<double> >::lartv(
Vector<double,complex<double> > &y,
const Vector<double,double> &c,
const Vector<double,complex<double> > &s) {
  int n=size();
  CHECK_SAME(n,y.size());
  CHECK_SAME(n,c.size());
  CHECK_SAME(n,s.size());
  F77NAME(zlartv)(n,addr(),1,y.addr(),1,c.addr(),s.addr(),1);
}

template class Vector<double,complex<double> >; 
template void testVector(const double&,const complex<double>&);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Matrix.C"
#ifdef DEBUG
template<> const complex<double>
  Matrix<double,complex<double> >::undefined_(
  numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
#endif
template<> const complex<double>
  Matrix<double,complex<double> >::zero_(double_zero_,double_zero_);
template<> const complex<double>
  Matrix<double,complex<double> >::one_(double_one_,double_zero_);
template<> const complex<double>
  Matrix<double,complex<double> >::mone_(double_mone_,double_zero_);

/*
template<> Matrix<double,complex<double> >*
Matrix<double,complex<double> >::transpose() const {
  int m=size(0),n=size(1);
  Matrix<double,complex<double> > *T=
    OPERATOR_NEW Matrix<double,complex<double> >(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(m,addr(0,j),1,T->addr(j,0),n);
  }
  return T;
}

template<> Matrix<double,complex<double> >*
Matrix<double,complex<double> >::conjugateTranspose() const {
  int m=size(0),n=size(1);
  Matrix<double,complex<double> > *T=
    OPERATOR_NEW Matrix<double,complex<double> >(n,m);
  for (int j=0;j<n;j++) {
    const complex<double> *col_j=addr(0,j);
    complex<double> *T_row_j=T->addr(j,0);
    for (int i=0;i<m;i++,col_j++,T_row_j+=n) {
      *T_row_j=conj(*col_j); 
    }
  }
  return T;
}
*/

template<> void Matrix<double,complex<double> >::interchangeColumns(
int i,int j) {
  int m=size(0),n=size(1);
  CHECK_BOUNDS(i,0,n)
  CHECK_BOUNDS(j,0,n)
  F77NAME(zswap)(m,addr(0,i),1,addr(0,j),1);
}

template<> void Matrix<double,complex<double> >::interchangeRows(
int i,int j) {
  int m=size(0),n=size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,m)
  F77NAME(zswap)(n,addr(i,0),m,addr(j,0),m);
}

template<> void Matrix<double,complex<double> >::gemv(
complex<double> alpha,const Vector<double,complex<double> > &x,
complex<double> beta,Vector<double,complex<double> > &y,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(x.size(),m);
    CHECK_SAME(y.size(),n);
    F77NAME(zgemv)(trans,m,n,alpha,addr(),m,x.addr(),1,beta,y.addr(),1);
  } else {
    CHECK_SAME(x.size(),n);
    CHECK_SAME(y.size(),m);
    F77NAME(zgemv)(trans,m,n,alpha,addr(),m,x.addr(),1,beta,y.addr(),1);
  }
}

template<> void Matrix<double,complex<double> >::ger(
complex<double> alpha,const Vector<double,complex<double> > &x,
const Vector<double,complex<double> > &y) {
  OBSOLETE("inappropriate for this class");
}

template<> void Matrix<double,complex<double> >::gerc(
complex<double> alpha,const Vector<double,complex<double> > &x,
const Vector<double,complex<double> > &y) {
  int m=size(0),n=size(1);
  CHECK_SAME(x.size(),m);
  CHECK_SAME(y.size(),n);
  F77NAME(zgerc)(m,n,alpha,x.addr(),1,y.addr(),1,addr(),m);
}

template<> void Matrix<double,complex<double> >::geru(
complex<double> alpha,const Vector<double,complex<double> > &x,
const Vector<double,complex<double> > &y) {
  int m=size(0),n=size(1);
  CHECK_SAME(x.size(),m);
  CHECK_SAME(y.size(),n);
  F77NAME(zgeru)(m,n,alpha,x.addr(),1,y.addr(),1,addr(),m);
}

template<> void Matrix<double,complex<double> >::gemm(
complex<double> alpha,const Matrix<double,complex<double> > &A,
const Matrix<double,complex<double> > &B,complex<double> beta,
char transa,char transb) {
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
  F77NAME(zgemm)(transa,transb,m,n,k,alpha,A.addr(),A.size(0),
    B.addr(),B.size(0),beta,addr(),m);
}

template<> double Matrix<double,complex<double> >::equilibrate(
Vector<double,double> &r,Vector<double,double> &c,
double &rowcnd,double &colcnd) const {
  int m=size(0),n=size(1);
  double amax=numeric_limits<double>::infinity();
  int info;
  F77NAME(zgeequ)(m,n,addr(),m,r.addr(),c.addr(),rowcnd,colcnd,
    amax,info);
  CHECK_SAME(info,0);
  return amax;
}

template<> void Matrix<double,complex<double> >::copyFrom(char uplo,
int m,int n,const Matrix<double,complex<double> > &A) {
  int s0=size(0),as0=A.size(0);
  m=min(m,min(s0,as0));
  n=min(n,min(size(1),A.size(1)));
  if (uplo=='A' || uplo=='a') {
    F77NAME(zlacpy)(uplo,m,n,A.addr(),as0,addr(),s0);
//  otherwise, dlacpy only copies a triangular part:
  } else if (uplo=='L' || uplo=='l') {
    for (int j=0;j<n;j++) {
      if (j<m) {
        F77NAME(zcopy)(m-j,A.addr(j,j),1,addr(j,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+1,m),A.addr(0,j),1,addr(0,j),1);
    }
  }
}

template<> void Matrix<double,complex<double> >::scale(char type,
int kl,int ku,double denominator,double numerator) {
  int m=size(0),info;
  F77NAME(zlascl)(type,kl,ku,denominator,numerator,m,size(1),addr(),
    m,info);
  CHECK_SAME(info,0);
}

template<> void Matrix<double,complex<double> >::set(char uplo,
complex<double> offdiag,complex<double> diag) {
  int m=size(0);
  F77NAME(zlaset)(uplo,m,size(1),offdiag,diag,addr(),m);
}

template<> int Matrix<double,complex<double> >::lastNonzeroColumn()
const {
  int m=size(0);
  return F77NAME(ilazlc)(m,size(1),addr(),m);
}

template<> int Matrix<double,complex<double> >::lastNonzeroRow() const {
  int m=size(0);
  return F77NAME(ilazlr)(m,size(1),addr(),m);
}

template<> double Matrix<double,complex<double> >::normFrobenius()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlange)('F',m,size(1),addr(),m,work);
}

template<> double Matrix<double,complex<double> >::normInfinity() const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double val=F77NAME(zlange)('I',m,size(1),addr(),m,work);
  delete [] work;
  return val;
}

template<> double Matrix<double,complex<double> >::normMaxEntry() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlange)('M',m,size(1),addr(),m,work);
}

template<> double Matrix<double,complex<double> >::normOne() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlange)('1',m,size(1),addr(),m,work);
}

template<> Matrix<double,complex<double> >*
Matrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(M.size(0),k);
  Matrix<double,complex<double> > *P=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zgemm)('N','N',m,n,k,complex_double_one_,addr(),m,M.addr(),k,
    complex_double_zero_,P->addr(),m);
  return P;
}

template<> Vector<double,complex<double> >*
Matrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(v.size(),n);
  Vector<double,complex<double> > *P=
    OPERATOR_NEW Vector<double,complex<double> >(m);
  F77NAME(zgemv)('N',m,n,complex_double_one_,addr(),m,v.addr(),1,
    complex_double_zero_,P->addr(),1);
  return P;
}

template<> void Matrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  int m=size(0), n=size(1);
  Matrix<double,complex<double> > AF(*this);
  int info;
  if (m==n) {
    CHECK_SAME(m,b.size())
    CHECK_SAME(n,x.size())
    x.copy(b);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(zgetrf)(m,m,AF.addr(),m,ipiv,info);
    if (info==0) {
      F77NAME(zgetrs)(trans,m,1,AF.addr(),m,ipiv,x.addr(),m,info);
    }
    CHECK_SAME(info,0)
    delete [] ipiv;
  } else {
    complex<double> w(numeric_limits<double>::infinity(),
      numeric_limits<double>::infinity());
    int lwork=-1;
    bool transposed=(trans!='N' && trans!='n');
    if (m>n) {
      if (transposed) {
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())
        x.copyFrom(n,b);
        F77NAME(zgels)(trans,m,n,1,AF.addr(),m,x.addr(),m,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zgels)(trans,m,n,1,AF.addr(),m,x.addr(),m,work,
          lwork,info);
        CHECK_SAME(info,0)
        delete [] work;
      } else {
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        Vector<double,complex<double> > xtmp(m);
        xtmp.copy(b);
        F77NAME(zgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),m,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),m,work,
          lwork,info);
        CHECK_SAME(info,0)

        x.copyFrom(n,xtmp);
        delete [] work;
      }
    } else {
      if (transposed) {
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())
        Vector<double,complex<double> > xtmp(n);
        xtmp.copy(b);
        F77NAME(zgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),n,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),n,work,
          lwork,info);
        CHECK_SAME(info,0)

        x.copyFrom(m,xtmp);
        delete [] work;
      } else {
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        x.copyFrom(m,b);
        F77NAME(zgels)(trans,m,n,1,AF.addr(),m,x.addr(),n,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zgels)(trans,m,n,1,AF.addr(),m,x.addr(),n,work,
          lwork,info);
        CHECK_SAME(info,0)
        delete [] work;
      }
    }
  }
}

template<> void Matrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  int m=size(0), n=size(1);
  Matrix<double,complex<double> > AF(*this);
  int info;
  bool left_side=(side=='L' || side=='l');
  bool transposed=(trans!='N' && trans!='n');
  if (m==n) {
    X.copy(B);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(zgetrf)(m,m,AF.addr(),m,ipiv,info);
    if (info==0) {
      if (left_side) {
        int k=B.size(1);
        CHECK_SAME(k,X.size(1))
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        F77NAME(zgetrs)(trans,m,k,AF.addr(),m,ipiv,X.addr(),m,info);
      } else {
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        CHECK_SAME(m,B.size(1))
        CHECK_SAME(n,X.size(1))
        if (transposed) {
          for (int j=0;j<m-1;j++) {
            if (ipiv[j]-1!=j) {
              F77NAME(zswap)(k,X.addr(0,j),1,X.addr(0,ipiv[j]-1),1);
            }
          }
          F77NAME(ztrsm)('R','L','T','U',k,m,complex_double_one_,
            AF.addr(),m,X.addr(),k);
          F77NAME(ztrsm)('R','U','T','N',k,m,complex_double_one_,
            AF.addr(),m,X.addr(),k);
        } else {
          F77NAME(ztrsm)('R','U','N','N',k,m,complex_double_one_,
            AF.addr(),m,X.addr(),k);
          F77NAME(ztrsm)('R','L','N','U',k,m,complex_double_one_,
            AF.addr(),m,X.addr(),k);
          for (int j=m-2;j>=0;j--) {
            if (ipiv[j]-1!=j) {
              F77NAME(zswap)(k,X.addr(0,j),1,X.addr(0,ipiv[j]-1),1);
            }
          }
        }
      }
    }
    CHECK_SAME(info,0)
    delete [] ipiv; ipiv=0;
  } else {
    complex<double> w(numeric_limits<double>::infinity(),
      numeric_limits<double>::infinity());
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
          F77NAME(zgels)(trans,m,n,k,AF.addr(),m,X.addr(),m,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w.real());
          complex<double> *work=
            OPERATOR_NEW_BRACKET(complex<double>,lwork);
          F77NAME(zgels)(trans,m,n,k,AF.addr(),m,X.addr(),m,work,
            lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        } else {
          CHECK_SAME(m,B.size(0))
          CHECK_SAME(n,X.size(0))
          Matrix<double,complex<double> > Xtmp(B);
          F77NAME(zgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),m,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w.real());
          complex<double> *work=
            OPERATOR_NEW_BRACKET(complex<double>,lwork);
          F77NAME(zgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),m,work,
            lwork,info);
          CHECK_SAME(info,0)

          X.copyFrom('A',n,k,Xtmp);
          delete [] work; work=0;
        }
      } else { // m>n, right side
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        complex<double> *tau=OPERATOR_NEW_BRACKET(complex<double>,n);
        F77NAME(zgeqrf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zgeqrf)(m,n,AF.addr(),m,tau,work,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        if (transposed) {
          CHECK_SAME(m,B.size(1))
          CHECK_SAME(n,X.size(1))
          Matrix<double,complex<double> > Xtmp(B);

          lwork=-1;
          F77NAME(zunmqr)('R','N',k,m,n,AF.addr(),m,tau,Xtmp.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w.real());
          work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
          F77NAME(zunmqr)('R','N',k,m,n,AF.addr(),m,tau,Xtmp.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
          X.copyFrom('A',k,n,Xtmp);

          F77NAME(ztrsm)('R','U','T','N',k,n,complex_double_one_,
            AF.addr(),m,X.addr(),k);
        } else {
          CHECK_SAME(n,B.size(1))
          CHECK_SAME(m,X.size(1))
          X=zero_;
          X.copyFrom('A',k,n,B);
          F77NAME(ztrsm)('R','U','N','N',k,n,complex_double_one_,
            AF.addr(),m,X.addr(),k);
          lwork=-1;
          F77NAME(zunmqr)('R','C',k,m,n,AF.addr(),m,tau,X.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w.real());
          work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
          F77NAME(zunmqr)('R','C',k,m,n,AF.addr(),m,tau,X.addr(),k,
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
          Matrix<double,complex<double> > Xtmp(B);
          F77NAME(zgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),n,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w.real());
          complex<double> *work=
            OPERATOR_NEW_BRACKET(complex<double>,lwork);
          F77NAME(zgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),n,work,
            lwork,info);
          CHECK_SAME(info,0)

          X.copyFrom('A',m,k,Xtmp);
          delete [] work; work=0;
        } else {
          CHECK_SAME(m,B.size(0))
          CHECK_SAME(n,X.size(0))
          X.copyFrom('A',m,k,B);
          F77NAME(zgels)(trans,m,n,k,AF.addr(),m,X.addr(),n,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w.real());
          complex<double> *work=
            OPERATOR_NEW_BRACKET(complex<double>,lwork);
          F77NAME(zgels)(trans,m,n,k,AF.addr(),m,X.addr(),n,work,
            lwork,info);
          CHECK_SAME(info,0)
          delete [] work;
        }
      } else { // m<n, right side
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        complex<double> *tau=OPERATOR_NEW_BRACKET(complex<double>,n);
        F77NAME(zgelqf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)
        lwork=static_cast<int>(w.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zgelqf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        if (transposed) {
          CHECK_SAME(m,B.size(1))
          CHECK_SAME(n,X.size(1))
          X=zero_;
          X.copyFrom('A',k,m,B);
          F77NAME(ztrsm)('R','L','T','N',k,m,complex_double_one_,
            AF.addr(),m,X.addr(),k);
          lwork=-1;
          F77NAME(zunmlq)('R','N',k,n,m,AF.addr(),m,tau,X.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w.real());
          work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
          F77NAME(zunmlq)('R','N',k,n,m,AF.addr(),m,tau,X.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        } else {
          CHECK_SAME(n,B.size(1))
          CHECK_SAME(m,X.size(1))
          Matrix<double,complex<double> > Xtmp(B);
          lwork=-1;
          F77NAME(zunmlq)('R','C',k,n,m,AF.addr(),m,tau,Xtmp.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w.real());
          work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
          F77NAME(zunmlq)('R','C',k,n,m,AF.addr(),m,tau,Xtmp.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
          X.copyFrom('A',k,m,Xtmp);

          F77NAME(ztrsm)('R','L','N','N',k,m,complex_double_one_,
            AF.addr(),m,X.addr(),k);
        }
        delete [] tau; tau=0;
      }
    }
  }
}

template class Matrix<double,complex<double> >; 
template void testMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SquareMatrix.C"

template<> SquareMatrix<double,complex<double> >*
SquareMatrix<double,complex<double> >::operator*(
const SquareMatrix<double,complex<double> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *P=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  F77NAME(zgemm)('N','N',n,n,n,complex_double_one_,addr(),n,S.addr(),n,
    complex_double_zero_,P->addr(),n);
  return P;
}

template<> Matrix<double,complex<double> >*
SquareMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *P=(m==n ?
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n) :
    OPERATOR_NEW Matrix<double,complex<double> >(m,n));
  F77NAME(zgemm)('N','N',m,n,m,complex_double_one_,addr(),m,M.addr(),m,
    complex_double_zero_,P->addr(),m);
  return P;
}

/*
template<> SquareMatrix<double,complex<double> >*
SquareMatrix<double,complex<double> >::transpose() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *T=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n,addr(0,j),1,T->addr(j,0),n);
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >*
SquareMatrix<double,complex<double> >::conjugateTranspose() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *T=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    const complex<double> *col_j=addr(0,j);
    complex<double> *T_row_j=T->addr(j,0);
    for (int i=0;i<n;i++,col_j++,T_row_j+=n) {
      *T_row_j=conj(*col_j);
    }
  }
  return T;
}
*/

template<> double 
SquareMatrix<double,complex<double> >::reciprocalConditionNumber(
char norm) const {
  int n=size(0),info;
  double rcond;
  double *rwork=OPERATOR_NEW_BRACKET(double,2*n);
  double anorm=F77NAME(zlange)(norm,n,n,addr(),n,rwork);

  SquareMatrix<double,complex<double> > *AF=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(zgetrf)(n,n,AF->addr(0,0),n,ipiv,info);
  CHECK_SAME(info,0)
  delete [] ipiv; ipiv=0;

  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(zgecon)(norm,n,AF->addr(),n,anorm,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] iwork; iwork=0;
  delete [] rwork; rwork=0;
  delete[] work; work=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SquareMatrix<double,complex<double> >*
SquareMatrix<double,complex<double> >::inverse() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *Ainv=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  *Ainv = *this;
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(zgetrf)(n,n,Ainv->addr(0,0),n,ipiv,info);
  CHECK_SAME(info,0)

  complex<double> w(numeric_limits<double>::infinity(),
    numeric_limits<double>::infinity());
  int lwork=-1;
  F77NAME(zgetri)(n,Ainv->addr(0,0),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zgetri)(n,Ainv->addr(0,0),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)

  delete [] ipiv;
  delete [] work;
  return Ainv;
}
*/

template<> Vector<double,complex<double> >*
SquareMatrix<double,complex<double> >::eigenvalues(
SquareMatrix<double,complex<double> > *&V,
SquareMatrix<double,complex<double> > *&U) const {
  int n=size(0);
  if (V!=0) CHECK_SAME(n,V->size(0));
  if (U!=0) CHECK_SAME(n,U->size(0));
  char jobvl=(V==0 ? 'N' : 'V');
  char jobvr=(U==0 ? 'N' : 'V');
  complex<double> *vl=(V==0 ? 0 : V->addr());
  complex<double> *vr=(U==0 ? 0 : U->addr());
  SquareMatrix<double,complex<double> > AF(n);
  AF.copy(*this);
  Vector<double,complex<double> > *lambda =
    OPERATOR_NEW Vector<double,complex<double> >(n);
  complex<double> w;
  double *rwork=OPERATOR_NEW_BRACKET(double,2*n);
  int lwork=-1,info;
  F77NAME(zgeev)(jobvl,jobvr,n,AF.addr(),n,lambda->addr(),vl,n,vr,n,&w,
    lwork,rwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zgeev)(jobvl,jobvr,n,AF.addr(),n,lambda->addr(),vl,n,vr,n,
    work,lwork,rwork,info);
  CHECK_TEST(info==0);
  delete [] work;
  delete [] rwork;
  return lambda;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const SquareMatrix<double,complex<double> > &S) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<double,complex<double> > *P=(m==n ?
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n) :
    OPERATOR_NEW Matrix<double,complex<double> >(m,n));
  F77NAME(zgemm)('N','N',m,n,n,complex_double_one_,M.addr(),m,S.addr(),n,
    complex_double_zero_,P->addr(),m);
  return P;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "TrapezoidalMatrix.C"
template<> complex<double>
  TrapezoidalMatrix<double,complex<double> >::outofbounds_=
  complex<double>(double_zero_,double_zero_);
template<> complex<double>
  TrapezoidalMatrix<double,complex<double> >::safety_=
  complex<double>(double_zero_,double_zero_);

template class TrapezoidalMatrix<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> Matrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::makeMatrix() const {
  int m=size(0),n=size(1);
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(this)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(m-j,addr(j,j),1,M->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*M)(j,j)=complex_double_one_;
      if (j+1<m) F77NAME(zcopy)(m-j-1,addr(j+1,j),1,M->addr(j+1,j),1);
    }
  }
  return M;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >& 
LowerTrapezoidalMatrix<double,complex<double> >::operator+=(
const LowerTrapezoidalMatrix<double,complex<double> > &L) {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0))
  CHECK_SAME(n,L.size(1))
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m-j,one_,L.addr(j,j),1,addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(zaxpy)(m-j-1,one_,L.addr(j+1,j),1,addr(j+1,j),1);
      }
      (*this)(j,j)+=complex_double_one_;
    }
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >&
LowerTrapezoidalMatrix<double,complex<double> >::operator-=(
const LowerTrapezoidalMatrix<double,complex<double> > &L) {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0))
  CHECK_SAME(n,L.size(1))
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m-j,mone_,L.addr(j,j),1,addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(zaxpy)(m-j-1,mone_,L.addr(j+1,j),1,addr(j+1,j),1);
      }
      (*this)(j,j)-=complex_double_one_;
    }
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >& 
LowerTrapezoidalMatrix<double,complex<double> >::operator*=(
complex<double> d) {
  int m=size(0),n=size(1);
  for (int j=0;j<n;j++) {
    F77NAME(zscal)(m-j,d,addr(j,j),1);
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >&
LowerTrapezoidalMatrix<double,complex<double> >::operator/=(
complex<double> d) {
  int m=size(0),n=size(1);
  complex<double> dinv=one_/d;
  for (int j=0;j<n;j++) {
    F77NAME(zscal)(m-j,dinv,addr(j,j),1);
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<double,complex<double> > *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(m,n);
  S->copy(*this);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m-j,complex_double_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m-j-1,complex_double_one_,L.addr(j+1,j),1,
        S->addr(j+1,j),1);
      (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *S=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(m-j,complex_double_one_,addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<double,complex<double> > *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(m,n);
  S->copy(*this);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m-j,complex_double_mone_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m-j-1,complex_double_mone_,L.addr(j+1,j),1,
        S->addr(j+1,j),1);
      (*S)(j,j)-=complex_double_one_;
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *D=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zlaset)('U',m,n,complex_double_zero_,complex_double_zero_,
    D->addr(),m);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(m-j,addr(j,j),1,D->addr(j,j),1);
    F77NAME(zaxpy)(m,complex_double_mone_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >* 
LowerTrapezoidalMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
// to compute the jth column of the product
//   [ L_11   0  ] [ 0 ] = [    0   ]
//   [ L_21 L_22 ] [ m ] = [ L_22 m ]
//   [ L_31 L_32 ]       = [ L_32 m ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  LowerTrapezoidalMatrix<double,complex<double> > *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(m,n);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(k-j,L.addr(j,j),1,P->addr(j,j),1);
      F77NAME(ztrmv)('L','N','N',k-j,addr(j,j),m,P->addr(j,j),1);
      if (m>k) {
        F77NAME(zgemv)('N',m-k,k-j,complex_double_one_,addr(k,j),m,
          L.addr(j,j),1,complex_double_zero_,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*P)(j,j)=(*this)(j,j);
      if (j<k-1) {
        F77NAME(zcopy)(k-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
        F77NAME(ztrmv)('L','N','N',k-j-1,addr(j+1,j+1),m,
          P->addr(j+1,j),1);
        F77NAME(zaxpy)(k-j-1,complex_double_one_,addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
      if (m>k) {
        if (j<k-1) {
          F77NAME(zgemv)('N',m-k,k-j-1,complex_double_one_,addr(k,j+1),m,
            L.addr(j+1,j),1,complex_double_zero_,P->addr(k,j),1);
        }
        F77NAME(zaxpy)(m-k,complex_double_one_,addr(k,j),1,
          P->addr(k,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<double,complex<double> > *P=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zlacpy)('A',k,n,M.addr(),k,P->addr(),m);
  F77NAME(ztrmm)('L','L','N','N',k,n,complex_double_one_,addr(),m,
    P->addr(),m);
  if (m>k) {
    F77NAME(zgemm)('N','N',m-k,n,k,complex_double_one_,addr(k,0),m,
      M.addr(),k,complex_double_zero_,P->addr(k,0),m);
  }
  return P;
}

template<> Vector<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<double,complex<double> > *p=
    OPERATOR_NEW Vector<double,complex<double> >(m);
  F77NAME(zcopy)(n,v.addr(),1,p->addr(),1);
  F77NAME(ztrmv)('L','N','N',n,addr(),m,p->addr(),1);
  if (m>n) {
    F77NAME(zgemv)('N',m-n,n,complex_double_one_,addr(n,0),m,v.addr(),1,
      complex_double_zero_,p->addr(n),1);
  }
  return p;
}

/*
template<> UpperTrapezoidalMatrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::transpose() const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<double,complex<double> > *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(m-j,addr(j,j),1,U->addr(j,j),n);
  }
  return U;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::conjugateTranspose()
const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<double,complex<double> > *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(n,m);
  for (int j=0;j<n;j++) {
    const complex<double> *col_j=addr(j,j);
    complex<double> *U_row_j=U->addr(j,j);
    for (int i=j;i<m;i++,col_j++,U_row_j+=n) {
      *U_row_j=conj(*col_j);
    }
  }
  return U;
}
*/

template<> Vector<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::trmv(
const Vector<double,complex<double> > &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<double,complex<double> > *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<double,complex<double> >(n);
    F77NAME(zcopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(ztrmv)('L',trans,'N',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(zgemv)(trans,m-n,n,one_,addr(n,0),m,x.addr(n),1,one_,
        p->addr(),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<double,complex<double> >(m);
    F77NAME(zcopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(ztrmv)('L','N','N',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(zgemv)('N',m-n,n,one_,addr(n,0),m,x.addr(),1,zero_,
        p->addr(n),1);
    }
  }
  return p;
}

template<> Matrix<double,complex<double> >*
LowerTrapezoidalMatrix<double,complex<double> >::trmm(
const Matrix<double,complex<double> > &X,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<double,complex<double> > *P=0;
  if (side=='L' || side=='l') {
    int k=X.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,X.size(0));
      P=OPERATOR_NEW Matrix<double,complex<double> >(n,k);
      F77NAME(zlacpy)('A',n,k,X.addr(),m,P->addr(),n);
      F77NAME(ztrmm)('L','L',trans,'N',n,k,one_,addr(),m,P->addr(),n);
      if (m>n) {
        F77NAME(zgemm)(trans,'N',n,k,m-n,one_,addr(n,0),m,X.addr(n,0),m,
          one_,P->addr(),n);
      }
    } else {
      CHECK_SAME(n,X.size(0));
      P=OPERATOR_NEW Matrix<double,complex<double> >(m,k);
      F77NAME(zlacpy)('A',n,k,X.addr(),n,P->addr(),m);
      F77NAME(ztrmm)('L','L','N','N',n,k,one_,addr(),m,P->addr(),m);
      if (m>n) {
        F77NAME(zgemm)('N','N',m-n,k,n,one_,addr(n,0),m,X.addr(),n,
          zero_,P->addr(n,0),m);
      }
    }
  } else {
    int k=X.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,X.size(1));
      P=OPERATOR_NEW Matrix<double,complex<double> >(k,m);
      F77NAME(zlacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(ztrmm)('R','L',trans,'N',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(zgemm)('N',trans,k,m-n,n,one_,X.addr(),k,addr(n,0),m,
          zero_,P->addr(0,n),k);
      }
    } else {
      CHECK_SAME(m,X.size(1));
      P=OPERATOR_NEW Matrix<double,complex<double> >(k,n);
      F77NAME(zlacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(ztrmm)('R','L','N','N',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(zgemm)('N','N',k,n,m-n,one_,X.addr(0,n),k,addr(n,0),m,
          one_,P->addr(),k);
      }
    }
  }
  return P;
}

template<> void LowerTrapezoidalMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(zcopy)(n,b.addr(),1,x.addr(),1);
    if (m>n) { // use trailing entries of x as free variables
      F77NAME(zgemv)(trans,m-n,n,mone_,addr(n,0),m,x.addr(n),1,one_,
        x.addr(),1);
    }
    F77NAME(ztrsv)('L',trans,'N',n,addr(),m,x.addr(),1);
  } else { // calling routine will have to check consistency conditions
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(zcopy)(n,b.addr(),1,x.addr(),1);
    F77NAME(ztrsv)('L','N','N',n,addr(),m,x.addr(),1);
  }
}

template<> void LowerTrapezoidalMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(zlacpy)('A',n,k,B.addr(),n,X.addr(),n);
      if (m>n) { // use these entries of X as free variables
        F77NAME(zgemm)(trans,'N',n,k,m-n,mone_,addr(n,0),m,X.addr(n,0),m,
          one_,X.addr(),m);
      }
      F77NAME(ztrsm)('L','L',trans,'N',n,k,one_,addr(),m,X.addr(),m);
    } else { // calling routine will have to check consistency conditions
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(zlacpy)('A',n,k,B.addr(),m,X.addr(),n);
      F77NAME(ztrsm)('L','L','N','N',n,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(zlacpy)('A',k,n,B.addr(),k,X.addr(),k);
      F77NAME(ztrsm)('R','L',trans,'N',k,n,one_,addr(),m,X.addr(),k);
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(zlacpy)('A',k,n,B.addr(),k,X.addr(),k);
      if (m>n) { // use these entries of X as free variables
        F77NAME(zgemm)('N','N',k,n,m-n,mone_,X.addr(0,n),n,addr(n,0),m,
          one_,X.addr(),k);
      }
      F77NAME(ztrsm)('R','L','N','N',k,n,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> double
LowerTrapezoidalMatrix<double,complex<double> >::normFrobenius() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('F','L','N',m,size(1),addr(),m,work);
}

template<> double
LowerTrapezoidalMatrix<double,complex<double> >::normInfinity() const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double result=F77NAME(zlantr)('I','L','N',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> double
LowerTrapezoidalMatrix<double,complex<double> >::normMaxEntry() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('M','L','N',m,size(1),addr(),m,work);
}

template<> double
LowerTrapezoidalMatrix<double,complex<double> >::normOne() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('O','L','N',m,size(1),addr(),m,work);
}

template<> Matrix<double,complex<double> >* operator+(
const Matrix<double,complex<double> > &M,
const LowerTrapezoidalMatrix<double,complex<double> > &L) {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *S=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  S->copy(M);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m-j,complex_double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(zaxpy)(m-j-1,complex_double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const LowerTrapezoidalMatrix<double,complex<double> > &L) {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *D=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  D->copy(M);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m-j,complex_double_mone_,L.addr(j,j),1,
        D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(zaxpy)(m-j-1,complex_double_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
      (*D)(j,j)-=complex_double_one_;
    }
  }
  return D;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const LowerTrapezoidalMatrix<double,complex<double> > &L) {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  char diagL=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0 ?
    'N' : 'U');
  int m=M.size(0),k=L.size(0),n=L.size(1);
  CHECK_SAME(k,M.size(1));
  Matrix<double,complex<double> > *P=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zlacpy)('A',m,n,M.addr(),m,P->addr(),m);
  F77NAME(ztrmm)('R','L','N',diagL,m,n,complex_double_one_,L.addr(),k,
    P->addr(),m);
  if (k>n) {
    F77NAME(zgemm)('N','N',m,n,k-n,complex_double_one_,M.addr(0,n),m,
      L.addr(n,0),k,complex_double_one_,P->addr(),m);
  }
  return P;
}
template class LowerTrapezoidalMatrix<double,complex<double> >;
template void testLowerTrapezoidalMatrix(double,complex<double> );

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double
LowerTriangularMatrix<double,complex<double> >::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info;
  F77NAME(ztrcon)(norm,'L','N',n,addr(),n,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork;
  delete [] work;
  return rcond;
}

/*
template<> LowerTriangularMatrix<double,complex<double> >*
LowerTriangularMatrix<double,complex<double> >::inverse() const {
  LowerTriangularMatrix<double,complex<double> > *L=
    OPERATOR_NEW LowerTriangularMatrix<double,complex<double> >(*this);
  int n=size(0);
  int info;
  F77NAME(ztrtri)('L','N',n,L->addr(),n,info);
  CHECK_TEST(info==0);
  return L;
}
*/

template class LowerTriangularMatrix<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void
UnitLowerTrapezoidalMatrix<double,complex<double> >::copyFrom(int m,int n,
const Matrix<double,complex<double> > &L) {
  m=min(m,min(size(0),L.size(0)));
  n=min(n,min(size(1),L.size(1)));
  for (int j=0;j<n;j++) {
    if (j+1<m) {
      F77NAME(zcopy)(m-j-1,L.addr(j+1,j),1,addr(j+1,j),1);
    }
  }
}

/*
template<> UnitUpperTrapezoidalMatrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::transpose() const {
  int m=size(0),n=size(1);
  UnitUpperTrapezoidalMatrix<double,complex<double> > *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<double,complex<double> >(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(m-j-1,addr(j+1,j),1,U->addr(j,j+1),n);
  }
  return U;
}

template<> UnitUpperTrapezoidalMatrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::conjugateTranspose()
const {
  int m=size(0),n=size(1);
  UnitUpperTrapezoidalMatrix<double,complex<double> > *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<double,complex<double> >(n,m);
  for (int j=0;j<n;j++) {
    const complex<double> *col_j=addr(j+1,j);
    complex<double> *U_row_j=U->addr(j,j+1);
    for (int i=j+1;i<m;i++,col_j++,U_row_j+=n) {
      *U_row_j=conj(*col_j);
    }
  }
  return U;
}
*/

template<> LowerTrapezoidalMatrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<double,complex<double> > *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(m,n);
  if (L_non_unit) {
    S->copy(L);
    for (int j=0;j<n;j++) {
      (*S)(j,j)+=complex_double_one_;
      F77NAME(zaxpy)(m-j-1,complex_double_one_,addr(j+1,j),1,
        S->addr(j+1,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=2.;
      if (j<m-1) {
        F77NAME(zcopy)(m-j-1,addr(j+1,j),1,S->addr(j+1,j),1);
        F77NAME(zaxpy)(m-j-1,complex_double_one_,addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *S=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    (*S)(j,j)+=complex_double_one_;
    F77NAME(zaxpy)(m-j-1,complex_double_one_,addr(j+1,j),1,
      S->addr(j+1,j),1);
  }
  return S;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<double,complex<double> > *dif=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(m,n);
  for (int j=0;j<n;j++) {
    (*dif)(j,j)=
      (L_non_unit ? complex_double_one_-L(j,j) : complex_double_zero_);
    if (j<m-1) {
      F77NAME(zcopy)(m-j-1,addr(j+1,j),1,dif->addr(j+1,j),1);
      F77NAME(zaxpy)(m-j-1,complex_double_mone_,L.addr(j+1,j),1,
        dif->addr(j+1,j),1);
    }
  }
  return dif;
}

template<> Matrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *dif=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zlaset)('U',m,n,complex_double_zero_,complex_double_zero_,
    dif->addr(),m);
  for (int j=0;j<n;j++) {
    (*dif)(j,j)=complex_double_one_-M(j,j);
    F77NAME(zcopy)(m-j-1,addr(j+1,j),1,dif->addr(j+1,j),1);
    F77NAME(zaxpy)(m-j-1,complex_double_mone_,M.addr(j+1,j),1,
      dif->addr(j+1,j),1);
  }
  return dif;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator*(
complex<double> d) const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<double,complex<double> > *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(m,n);
  for (int j=0;j<n;j++) {
    (*P)(j,j)=d;
    if (j<m-1) {
      F77NAME(zcopy)(m-j-1,addr(j+1,j),1,P->addr(j+1,j),1);
      F77NAME(zscal)(m-j-1,d,P->addr(j+1,j),1);
    }
  }
  return P;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator/(
complex<double> d) const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<double,complex<double> > *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(m,n);
  complex<double> dinv=complex_double_one_/d;
  for (int j=0;j<n;j++) {
    (*P)(j,j)=dinv;
    if (j<m-1) {
      F77NAME(zcopy)(m-j-1,addr(j+1,j),1,P->addr(j+1,j),1);
      F77NAME(zscal)(m-j-1,dinv,P->addr(j+1,j),1);
    }
  }
  return P;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
// to compute the jth column of the product
//   [ L_11   0  ] [ 0 ] = [    0   ]
//   [ L_21 L_22 ] [ m ] = [ L_22 m ]
//   [ L_31 L_32 ]       = [ L_32 m ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  LowerTrapezoidalMatrix<double,complex<double> > *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(m,n);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(k-j,L.addr(j,j),1,P->addr(j,j),1);
      F77NAME(ztrmv)('L','N','U',k-j,addr(j,j),m,P->addr(j,j),1);
      if (m>k) {
        F77NAME(zgemv)('N',m-k,k-j,complex_double_one_,addr(k,j),m,
          L.addr(j,j),1,complex_double_zero_,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*P)(j,j)=complex_double_one_;
      if (j<k-1) {
        F77NAME(zcopy)(k-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
        F77NAME(ztrmv)('L','N','U',k-j-1,addr(j+1,j+1),m,
          P->addr(j+1,j),1);
        F77NAME(zaxpy)(k-j-1,complex_double_one_,addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
      if (m>k) {
        if (j<k-1) {
          F77NAME(zgemv)('N',m-k,k-j-1,complex_double_one_,addr(k,j+1),m,
            L.addr(j+1,j),1,complex_double_zero_,P->addr(k,j),1);
        }
        F77NAME(zaxpy)(m-k,complex_double_one_,addr(k,j),1,
          P->addr(k,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<double,complex<double> > *P=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n,zero_);
  F77NAME(zlacpy)('A',k,n,M.addr(),k,P->addr(),m);
  F77NAME(ztrmm)('L','L','N','U',k,n,complex_double_one_,addr(),m,
    P->addr(),m);
  if (m>k) {
    F77NAME(zgemm)('N','N',m-k,n,k,complex_double_one_,addr(k,0),m,
      M.addr(),k,complex_double_one_,P->addr(k,0),m);
  }
  return P;
}

template<> Vector<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<double,complex<double> > *p=
    OPERATOR_NEW Vector<double,complex<double> >(m,zero_);
  F77NAME(zcopy)(n,v.addr(),1,p->addr(),1);
  F77NAME(ztrmv)('L','N','U',n,addr(),m,p->addr(),1);
  if (m>n) {
    F77NAME(zgemv)('N',m-n,n,complex_double_one_,addr(n,0),m,
      v.addr(),1,complex_double_one_,p->addr(n),m);
  }
  return p;
}

template<> Vector<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::trmv(
const Vector<double,complex<double> > &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<double,complex<double> > *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<double,complex<double> >(n);
    F77NAME(zcopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(ztrmv)('L',trans,'U',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(zgemv)(trans,m-n,n,one_,addr(n,0),m,x.addr(n),1,one_,
        p->addr(),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<double,complex<double> >(m);
    F77NAME(zcopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(ztrmv)('L','N','U',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(zgemv)('N',m-n,n,one_,addr(n,0),m,x.addr(),1,zero_,
        p->addr(n),1);
    }
  }
  return p;
}

template<> Matrix<double,complex<double> >*
UnitLowerTrapezoidalMatrix<double,complex<double> >::trmm(
const Matrix<double,complex<double> > &X,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<double,complex<double> > *P=0;
  if (side=='L' || side=='l') {
    int k=X.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,X.size(0));
      P=OPERATOR_NEW Matrix<double,complex<double> >(n,k);
      F77NAME(zlacpy)('A',n,k,X.addr(),m,P->addr(),n);
      F77NAME(ztrmm)('L','L',trans,'U',n,k,one_,addr(),m,P->addr(),n);
      if (m>n) {
        F77NAME(zgemm)(trans,'N',n,k,m-n,one_,addr(n,0),m,X.addr(n,0),m,
          one_,P->addr(),n);
      }
    } else {
      CHECK_SAME(n,X.size(0));
      P=OPERATOR_NEW Matrix<double,complex<double> >(m,k);
      F77NAME(zlacpy)('A',n,k,X.addr(),n,P->addr(),m);
      F77NAME(ztrmm)('L','L','N','U',n,k,one_,addr(),m,P->addr(),m);
      if (m>n) {
        F77NAME(zgemm)('N','N',m-n,k,n,one_,addr(n,0),m,X.addr(),n,
          zero_,P->addr(n,0),m);
      }
    }
  } else {
    int k=X.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,X.size(1));
      P=OPERATOR_NEW Matrix<double,complex<double> >(k,m);
      F77NAME(zlacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(ztrmm)('R','L',trans,'U',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(zgemm)('N',trans,k,m-n,n,one_,X.addr(),k,addr(n,0),m,
          zero_,P->addr(0,n),k);
      }
    } else {
      CHECK_SAME(m,X.size(1));
      P=OPERATOR_NEW Matrix<double,complex<double> >(k,n);
      F77NAME(zlacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(ztrmm)('R','L','N','U',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(zgemm)('N','N',k,n,m-n,one_,X.addr(0,n),k,addr(n,0),m,
          one_,P->addr(),k);
      }
    }
  }
  return P;
}

template<> void
UnitLowerTrapezoidalMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(zcopy)(n,b.addr(),1,x.addr(),1);
    if (m>n) { // use trailing entries of x as free variables
      F77NAME(zgemv)(trans,m-n,n,mone_,addr(n,0),m,x.addr(n),1,one_,
        x.addr(),1);
    }
    F77NAME(ztrsv)('L',trans,'U',n,addr(),m,x.addr(),1);
  } else { // calling routine will have to check consistency conditions
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(zcopy)(n,b.addr(),1,x.addr(),1);
    F77NAME(ztrsv)('L','N','U',n,addr(),m,x.addr(),1);
  }
}

template<> void
UnitLowerTrapezoidalMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(zlacpy)('A',n,k,B.addr(),n,X.addr(),n);
      if (m>n) { // use these entries of X as free variables
        F77NAME(zgemm)(trans,'N',n,k,m-n,mone_,addr(n,0),m,X.addr(n,0),m,
          one_,X.addr(),m);
      }
      F77NAME(ztrsm)('L','L',trans,'U',n,k,one_,addr(),m,X.addr(),m);
    } else { // calling routine will have to check consistency conditions
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(zlacpy)('A',n,k,B.addr(),m,X.addr(),n);
      F77NAME(ztrsm)('L','L','N','U',n,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(zlacpy)('A',k,n,B.addr(),k,X.addr(),k);
      F77NAME(ztrsm)('R','L',trans,'U',k,n,one_,addr(),m,X.addr(),k);
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(zlacpy)('A',k,n,B.addr(),k,X.addr(),k);
      if (m>n) { // use these entries of X as free variables
        F77NAME(zgemm)('N','N',k,n,m-n,mone_,X.addr(0,n),n,addr(n,0),m,
          one_,X.addr(),k);
      }
      F77NAME(ztrsm)('R','L','N','U',k,n,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> double
UnitLowerTrapezoidalMatrix<double,complex<double> >::normFrobenius()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('F','L','U',m,size(1),addr(),m,work);
}

template<> double
UnitLowerTrapezoidalMatrix<double,complex<double> >::normInfinity()
const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double result=F77NAME(zlantr)('I','L','U',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> double
UnitLowerTrapezoidalMatrix<double,complex<double> >::normMaxEntry()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('M','L','U',m,size(1),addr(),m,work);
}

template<> double
UnitLowerTrapezoidalMatrix<double,complex<double> >::normOne() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('O','L','U',m,size(1),addr(),m,work);
}

template class UnitLowerTrapezoidalMatrix<double,complex<double> >;
template void testUnitLowerTrapezoidalMatrix(double,complex<double> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double UnitLowerTriangularMatrix<double,complex<double> >
::reciprocalConditionNumber(char norm) const {
  int n=size(0);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info;
  F77NAME(ztrcon)(norm,'L','U',n,addr(),n,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork;
  delete [] work;
  return rcond;
}

/*
template<> UnitLowerTriangularMatrix<double,complex<double> >*
UnitLowerTriangularMatrix<double,complex<double> >::inverse() const {
  UnitLowerTriangularMatrix<double,complex<double> > *result=
    OPERATOR_NEW UnitLowerTriangularMatrix<double,complex<double> >(*this);
  int n=size(0);
  int info;
  F77NAME(ztrtri)('L','U',n,result->addr(),n,info);
  CHECK_TEST(info==0);
  return result;
}
*/

template class UnitLowerTriangularMatrix<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> Matrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::makeMatrix() const {
  int m=size(0),n=size(1);
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(this)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+1,m),addr(0,j),1,M->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j,m),addr(0,j),1,M->addr(0,j),1);
      if (j<m) (*M)(j,j)=complex_double_one_;
    }
  }
  return M;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >& 
UpperTrapezoidalMatrix<double,complex<double> >::operator+=(
const UpperTrapezoidalMatrix<double,complex<double> > &M) {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0))
  CHECK_SAME(n,M.size(1))
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(j+1,m),one_,M.addr(0,j),1,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >&
UpperTrapezoidalMatrix<double,complex<double> >::operator-=(
const UpperTrapezoidalMatrix<double,complex<double> > &M) {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0))
  CHECK_SAME(n,M.size(1))
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(j+1,m),mone_,M.addr(0,j),1,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >& 
UpperTrapezoidalMatrix<double,complex<double> >::operator*=(
complex<double> d) {
  int m=size(0),n=size(1);
  for (int j=0;j<n;j++) {
    F77NAME(zscal)(min(j+1,m),d,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >&
UpperTrapezoidalMatrix<double,complex<double> >::operator/=(
complex<double> d) {
  int m=size(0),n=size(1);
  complex<double> dinv=complex_double_one_/d;
  for (int j=0;j<n;j++) {
    F77NAME(zscal)(min(j+1,m),dinv,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<double,complex<double> > *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(m,n);
  if (U_non_unit) {
    S->copy(U);
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(min(j+1,m),complex_double_one_,addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      if (j>0) {
        F77NAME(zaxpy)(min(j,m),complex_double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(m,zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(j+1,addr(0,j),1,S->addr(0,j),1);
      F77NAME(zaxpy)(m-j,complex_double_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(j+1,addr(0,j),1,S->addr(0,j),1);
      if (j<m-1) {
        F77NAME(zaxpy)(m-j-1,complex_double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *S=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(m,j+1),complex_double_one_,addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<double,complex<double> > *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(m,n);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      F77NAME(zaxpy)(min(j+1,m),complex_double_mone_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      if (j>0) {
        F77NAME(zaxpy)(min(j,m),complex_double_mone_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)-=complex_double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *D=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(j+1,addr(0,j),1,D->addr(0,j),1);
      F77NAME(zaxpy)(m-j,complex_double_mone_,L.addr(j,j),1,
        D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(j+1,addr(0,j),1,D->addr(0,j),1);
      if (j<m-1) {
        F77NAME(zaxpy)(m-j-1,complex_double_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
      (*D)(j,j)-=complex_double_one_;
    }
  }
  return D;
}

template<> Matrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *D=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(j+1,addr(0,j),1,D->addr(0,j),1);
    F77NAME(zaxpy)(m,complex_double_mone_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >* 
UpperTrapezoidalMatrix<double,complex<double> >::operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
// [ U_1 , U_2 ] [ V_11 , V_12 , V_13 ]
//               [   0  , V_22 , V_23 ]
//   = [ U_1 V_11 , U_1 V_12 + U_2 V_22 , U_1 V_13 + U_2 V_23 ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  int m=size(0),k=size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  UpperTrapezoidalMatrix<double,complex<double> > *P=OPERATOR_NEW
    UpperTrapezoidalMatrix<double,complex<double> >(m,n,
      complex_double_zero_);
  if (U_non_unit) {
    for (int j=0;j<m;j++) {
      F77NAME(zcopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ztrmv)('U','N','N',j+1,addr(),m,P->addr(0,j),1);
    }
    for (int j=m;j<k;j++) {
      F77NAME(zcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ztrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
      F77NAME(zgemv)('N',m,j-m+1,complex_double_one_,addr(0,m),m,
        U.addr(m,j),1,complex_double_one_,P->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(zcopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(ztrmv)('U','N','N',j,addr(),m,P->addr(0,j),1);
        F77NAME(zaxpy)(j,complex_double_one_,addr(0,j),1,P->addr(0,j),1);
      }
      (*P)(j,j)=(*this)(j,j);
    }
    for (int j=m;j<k;j++) {
      F77NAME(zcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ztrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
      if (j>m) {
        F77NAME(zgemv)('N',m,j-m,complex_double_one_,addr(0,m),m,
          U.addr(m,j),1,complex_double_one_,P->addr(0,j),1);
      }
      F77NAME(zaxpy)(m,complex_double_one_,addr(0,j),1,P->addr(0,j),1);
    }
  }
  for (int j=k;j<n;j++) {
    F77NAME(zcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
    F77NAME(ztrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
    F77NAME(zgemv)('N',m,k-m,complex_double_one_,addr(0,m),m,
      U.addr(m,j),1,complex_double_one_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  if (m>=n) {
//  note that
//  [ U_11 U_12 U_13 ] [ L_1 ] = [ U_11 L_1 + U_12 L_2 + U_13 L_3 ]
//  [      U_22 U_23 ] [ L_2 ] = [            U_22 L_2 + U_23 L_3 ]
//                     [ L_3 ]
    if (L_non_unit) { // U_11 L_1:
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(zgemv)('N',j,n-j,complex_double_one_,addr(0,j),m,
            L.addr(j,j),1,complex_double_zero_,P->addr(0,j),1);
        }
        F77NAME(zcopy)(n-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(ztrmv)('U','N','N',n-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(zaxpy)(j,complex_double_one_,addr(0,j),1,
            P->addr(0,j),1);
          if (j<n-1) {
            F77NAME(zgemv)('N',j,n-j-1,complex_double_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,complex_double_zero_,P->addr(0,j),1);
          }
        }
        if (j<n-1) {
          F77NAME(zcopy)(n-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(ztrmv)('U','N','N',n-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(zdotu)(n-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (n<k) { // U_12 L_2 + U_13 L_3:
      F77NAME(zgemm)('N','N',n,n,k-n,double_one_,addr(0,n),m,
        L.addr(n,0),k,double_one_,P->addr(0,0),m);
    }
    if (n<m) {
      F77NAME(zlacpy)('A',m-n,n,L.addr(n,0),k,P->addr(n,0),m);
      F77NAME(ztrmm)('L','U','N','N',m-n,n,complex_double_one_,
        addr(n,n),m,P->addr(n,0),m); // U_22 L_2
      if (m<k) {
        F77NAME(zgemm)('N','N',m-n,n,k-m,complex_double_one_,addr(n,m),m,
          L.addr(m,0),k,complex_double_one_,P->addr(n,0),m); // U_23 L_3
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
          F77NAME(zgemv)('N',j,m-j,complex_double_one_,addr(0,j),m,
            L.addr(j,j),1,complex_double_zero_,P->addr(0,j),1);
        }
        F77NAME(zcopy)(m-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(ztrmv)('U','N','N',m-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(zaxpy)(j,complex_double_one_,addr(0,j),1,
            P->addr(0,j),1);
          if (j<m-1) {
            F77NAME(zgemv)('N',j,m-j-1,complex_double_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,complex_double_zero_,P->addr(0,j),1);
          }
        }
        if (j<m-1) {
          F77NAME(zcopy)(m-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(ztrmv)('U','N','N',m-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(zdotu)(m-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (m<k) { // U_2 L_21 + U_3 L_31 
      F77NAME(zgemm)('N','N',m,m,k-m,complex_double_one_,addr(0,m),m,
        L.addr(m,0),k,complex_double_one_,P->addr(0,0),m);
    }
    char diagL=(L_non_unit ? 'N' : 'U');
    F77NAME(zlacpy)('A',m,n-m,addr(0,m),m,P->addr(0,m),m);
    F77NAME(ztrmm)('R','L','N',diagL,m,n-m,complex_double_one_,
      L.addr(m,m),k,P->addr(0,m),m); // U_2 L_22
    if (n<k) {
      F77NAME(zgemm)('N','N',m,n-m,k-n,complex_double_one_,addr(0,n),m,
        L.addr(n,m),k,complex_double_one_,P->addr(0,m),m); // U_3 L_32
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
// [ U_1 U_2 ] [ M_1 ] = U_1 M_1 + U_2 M_2
//             [ M_2 ]
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<double,complex<double> > *P=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zlacpy)('A',m,n,M.addr(),k,P->addr(),m);
  F77NAME(ztrmm)('L','U','N','N',m,n,complex_double_one_,addr(),m,
    P->addr(),m); // U_1 M_1
  if (k>m) {
    F77NAME(zgemm)('N','N',m,n,k-m,complex_double_one_,addr(0,m),m,
      M.addr(m,0),k,complex_double_one_,P->addr(),m); // U_2 M_2
  }
  return P;
}

template<> Vector<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<double,complex<double> > *p=
    OPERATOR_NEW Vector<double,complex<double> >(m);
  F77NAME(zcopy)(m,v.addr(),1,p->addr(),1);
  F77NAME(ztrmv)('U','N','N',m,addr(),m,p->addr(),1);
  F77NAME(zgemv)('N',m,n-m,complex_double_one_,addr(0,m),m,v.addr(m),1,
    complex_double_one_,p->addr(),1);
  return p;
}

/*
template<> LowerTrapezoidalMatrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::transpose() const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<double,complex<double> > *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(n,m);
  for (int i=0;i<m;i++) {
    F77NAME(zcopy)(n-i,addr(i,i),m,L->addr(i,i),1);
  }
  return L;
}

template<> LowerTrapezoidalMatrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::conjugateTranspose()
const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<double,complex<double> > *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<double,complex<double> >(n,m);
  for (int i=0;i<m;i++) {
    const complex<double> *row_i=addr(i,i);
    complex<double> *L_col_i=L->addr(i,i);
    for (int j=i;j<n;j++,row_i+=m,L_col_i++) {
      *L_col_i=conj(*row_i);
    }
  }
  return L;
}
*/

template<> Vector<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::trmv(
const Vector<double,complex<double> > &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<double,complex<double> > *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<double,complex<double> >(n);
    F77NAME(zcopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(ztrmv)('U',trans,'N',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(zgemv)(trans,m,n-m,one_,addr(0,m),m,x.addr(),1,zero_,
        p->addr(m),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<double,complex<double> >(m);
    F77NAME(zcopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(ztrmv)('U','N','N',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(zgemv)('N',m,n-m,one_,addr(0,m),m,x.addr(m),1,one_,
        p->addr(),1);
    }
  }
  return p;
}

template<> Matrix<double,complex<double> >*
UpperTrapezoidalMatrix<double,complex<double> >::trmm(
const Matrix<double,complex<double> > &M,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<double,complex<double> > *P=0;
  if (side=='L' || side=='l') {
    int k=M.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,M.size(0));
      P=OPERATOR_NEW Matrix<double,complex<double> >(n,k);
      F77NAME(zlacpy)('A',m,k,M.addr(),m,P->addr(),n);
      F77NAME(ztrmm)('L','U',trans,'N',m,k,one_,addr(),m,P->addr(),n);
      if (n>m) {
        F77NAME(zgemm)(trans,'N',n-m,k,m,one_,addr(0,m),m,M.addr(),m,
          zero_,P->addr(m,0),n);
      }
    } else {
      CHECK_SAME(n,M.size(0));
      P=OPERATOR_NEW Matrix<double,complex<double> >(m,k);
      F77NAME(zlacpy)('A',m,k,M.addr(),n,P->addr(),m);
      F77NAME(ztrmm)('L','U','N','N',m,k,one_,addr(),m,P->addr(),m);
      if (n>m) {
        F77NAME(zgemm)('N','N',m,k,n-m,one_,addr(0,m),m,M.addr(m,0),n,
          one_,P->addr(),m);
      }
    }
  } else {
    int k=M.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,M.size(1));
      P=OPERATOR_NEW Matrix<double,complex<double> >(k,m);
      F77NAME(zlacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(ztrmm)('R','U',trans,'N',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(zgemm)('N',trans,k,m,n-m,one_,M.addr(0,m),k,addr(0,m),m,
          one_,P->addr(),k);
      }
    } else {
      CHECK_SAME(m,M.size(1));
      P=OPERATOR_NEW Matrix<double,complex<double> >(k,n);
      F77NAME(zlacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(ztrmm)('R','U','N','N',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(zgemm)('N','N',k,n-m,m,one_,M.addr(),k,addr(0,m),m,
          zero_,P->addr(0,m),k);
      }
    }
  }
  return P;
}

template<> double
UpperTrapezoidalMatrix<double,complex<double> >::normFrobenius() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('F','U','N',m,size(1),addr(),m,work);
}

template<> double
UpperTrapezoidalMatrix<double,complex<double> >::normInfinity() const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double result=F77NAME(zlantr)('I','U','N',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> double
UpperTrapezoidalMatrix<double,complex<double> >::normMaxEntry() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('M','U','N',m,size(1),addr(),m,work);
}

template<> double
UpperTrapezoidalMatrix<double,complex<double> >::normOne() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('O','U','N',m,size(1),addr(),m,work);
}

template<> void UpperTrapezoidalMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  char diag=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(this) ==0 ?
    'N' : 'U');
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    // calling routine will have to check consistency conditions
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(zcopy)(m,b.addr(),1,x.addr(),1);
    F77NAME(ztrsv)('U',trans,'N',m,addr(),m,x.addr(),1);
  } else {
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(zcopy)(m,b.addr(),1,x.addr(),1);
    if (n>m) { // use trailing entries of x as free variables
      F77NAME(zgemv)('N',m,n-m,mone_,addr(0,m),m,x.addr(m),1,one_,
        x.addr(),1);
    }
    F77NAME(ztrsv)('U','N','N',m,addr(),m,x.addr(),1);
  }
}

template<> void UpperTrapezoidalMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  char diag=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(this) ==0 ?
    'N' : 'U');
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(zlacpy)('A',m,k,B.addr(),n,X.addr(),m);
      F77NAME(ztrsm)('L','U',trans,'N',m,k,one_,addr(),m,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(zlacpy)('A',m,k,B.addr(),m,X.addr(),n);
      if (n>m) { // use these entries of X as free variables
        F77NAME(zgemm)('N','N',m,k,n-m,mone_,addr(0,m),m,X.addr(m,0),n,
          one_,X.addr(),n);
      }
      F77NAME(ztrsm)('L','U','N','N',m,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(zlacpy)('A',k,m,B.addr(),k,X.addr(),k);
      if (n>m) { // use these entries of X as free variables
        F77NAME(zgemm)('N',trans,k,m,n-m,mone_,X.addr(0,m),k,addr(0,m),m,
          one_,X.addr(),k);
      }
      F77NAME(ztrsm)('R','U',trans,'N',k,m,one_,addr(),m,X.addr(),k);
    } else { // calling routine will have to check consistency
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(zlacpy)('A',k,m,B.addr(),k,X.addr(),k);
      F77NAME(ztrsm)('R','U','N','N',k,m,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> SquareMatrix<double,complex<double> >* operator+(
const UnitLowerTrapezoidalMatrix<double,complex<double> > &L,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    (*S)(j,j)=complex_double_one_;
    if (j<m-1) {
      F77NAME(zcopy)(m-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
    }
    if (U_non_unit) {
      F77NAME(zaxpy)(j+1,complex_double_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    if (L_non_unit) {
      F77NAME(zcopy)(m-j,L.addr(j,j),1,S->addr(j,j),1);
    } else {
      (*S)(j,j)=complex_double_one_;
      if (j<m-1) {
        F77NAME(zcopy)(m-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
    }
    if (U_non_unit) {
      F77NAME(zaxpy)(j+1,complex_double_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >* operator+(
const Matrix<double,complex<double> > &M,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(not_trapezoidal);
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  Matrix<double,complex<double> > *S=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  S->copy(M);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(min(j+1,m),complex_double_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j,m),complex_double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const UnitLowerTrapezoidalMatrix<double,complex<double> > &L,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  SquareMatrix<double,complex<double> > *D=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    (*D)(j,j)=complex_double_one_;
    if (j<m-1) {
      F77NAME(zcopy)(m-j-1,L.addr(j+1,j),1,D->addr(j+1,j),1);
    }
    if (U_non_unit) {
      F77NAME(zaxpy)(j+1,complex_double_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_mone_,U.addr(0,j),1,
          D->addr(0,j),1);
      }
      (*D)(j,j)-=complex_double_one_;
    }
  }
  return D;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  SquareMatrix<double,complex<double> > *D=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    if (L_non_unit) {
      F77NAME(zcopy)(m-j,L.addr(j,j),1,D->addr(j,j),1);
    } else {
      (*D)(j,j)=complex_double_one_;
      if (j<m-1) {
        F77NAME(zcopy)(m-j-1,L.addr(j+1,j),1,D->addr(j+1,j),1);
      }
    }
    if (U_non_unit) {
      F77NAME(zaxpy)(j+1,complex_double_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_mone_,U.addr(0,j),1,
          D->addr(0,j),1);
      }
      (*D)(j,j)-=complex_double_one_;
    }
  }
  return D;
}

template<> Matrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(not_trapezoidal);
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  Matrix<double,complex<double> > *D=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  D->copy(M);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(min(j+1,m),complex_double_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j,m),complex_double_mone_,U.addr(0,j),1,
          D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)-=complex_double_one_;
    }
  }
  return D;
}

template<> Matrix<double,complex<double> >* operator*(
const UnitLowerTrapezoidalMatrix<double,complex<double> > &L,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
//  Note that
//  [ L_1 ] [ U_1 U_2 ] = [ L_1 U_1 , L_1 U_2 ]
//  [ L_2 ]             = [ L_2 U_1 , L_2 U_2 ]
//  and that
//  [ L_11      ] [ u ] = [ L_11 u ]
//  [ L_21 L_22 ] [   ] = [ L_21 u ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U) ==0);
  int m=L.size(0),k=L.size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<k;j++) { // L_1 U_1
    if (U_non_unit) {
      F77NAME(zcopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ztrmv)('L','N','U',j+1,L.addr(),m,P->addr(0,j),1);
      if (j+1<k) {
        F77NAME(zgemv)('N',k-j-1,j+1,complex_double_one_,L.addr(j+1,0),m,
          U.addr(0,j),1,complex_double_zero_,P->addr(j+1,j),1);
      }
    } else {
      if (j>0) {
        F77NAME(zcopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(ztrmv)('L','N','U',j,L.addr(),m,P->addr(0,j),1);
        (*P)(j,j)=F77NAME(zdotu)(j,L.addr(j,0),m,U.addr(0,j),1);
        if (j+1<k) {
          F77NAME(zgemv)('N',k-j-1,j,complex_double_one_,L.addr(j+1,0),m,
            U.addr(0,j),1,complex_double_zero_,P->addr(j+1,j),1);
        }
      }
      (*P)(j,j)+=L(j,j);
      if (j+1<k) {
        F77NAME(zaxpy)(k-j-1,complex_double_one_,L.addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
    }
  }
  if (n>k) { // L_1 U_2
    F77NAME(zlacpy)('A',k,n-k,U.addr(0,k),k,P->addr(0,k),m);
    F77NAME(ztrmm)('L','L','N','U',k,n-k,complex_double_one_,L.addr(),m,
      P->addr(0,k),m);
  }
  if (m>k) {
    char diagU=(U_non_unit ? 'N' : 'U');
    F77NAME(zlacpy)('A',m-k,k,L.addr(k,0),m,P->addr(k,0),m);
    F77NAME(ztrmm)('R','U','N',diagU,m-k,k,complex_double_one_,U.addr(),k,
      P->addr(k,0),m); // L_2 U_1
    if (n>k) { // L_2 U_2
      F77NAME(zgemm)('N','N',m-k,n-k,k,complex_double_one_,L.addr(k,0),m,
        U.addr(0,k),k,complex_double_zero_,P->addr(k,k),m);
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >* operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
//  Note that
//  [ L_1 ] [ U_1 U_2 ] = [ L_1 U_1 , L_1 U_2 ]
//  [ L_2 ]             = [ L_2 U_1 , L_2 U_2 ]
//  and that
//  [ L_11      ] [ u ] = [ L_11 u ]
//  [ L_21 L_22 ] [   ] = [ L_21 u ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  char diagL=(L_non_unit ? 'N' : 'U');
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  int m=L.size(0),k=L.size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<k;j++) { // L_1 U_1
    if (U_non_unit) {
      F77NAME(zcopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ztrmv)('L','N',diagL,j+1,L.addr(),m,P->addr(0,j),1);
      if (j+1<k) {
        F77NAME(zgemv)('N',k-j-1,j+1,complex_double_one_,L.addr(j+1,0),m,
          U.addr(0,j),1,complex_double_zero_,P->addr(j+1,0),1);
      }
    } else {
      if (j>0) {
        F77NAME(zcopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(ztrmv)('L','N',diagL,j,L.addr(),m,P->addr(0,j),1);
        (*P)(j,j)=F77NAME(zdotu)(j,L.addr(j,0),m,U.addr(0,j),1);
        if (j+1<k) {
          F77NAME(zgemv)('N',k-j-1,j,complex_double_one_,L.addr(j+1,0),m,
            U.addr(0,j),1,complex_double_zero_,P->addr(j+1,j),1);
        }
      }
      (*P)(j,j)+=(L_non_unit ? L(j,j) : complex_double_one_);
      if (j+1<k) {
        F77NAME(zaxpy)(k-j-1,complex_double_one_,L.addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
    }
  }
  if (n>k) { // L_1 U_2
    F77NAME(zlacpy)('A',k,n-k,U.addr(0,k),k,P->addr(0,k),m);
    F77NAME(ztrmm)('L','L','N',diagL,k,n-k,complex_double_one_,L.addr(),m,
      P->addr(0,k),m);
  }
  if (m>k) {
    char diagU=(U_non_unit ? 'N' : 'U');
    F77NAME(zlacpy)('A',m-k,k,L.addr(k,0),m,P->addr(k,0),m);
    F77NAME(ztrmm)('R','U','N',diagU,m-k,k,complex_double_one_,U.addr(),k,
      P->addr(k,0),m); // L_2 U_1
    if (n>k) { // L_2 U_2
      F77NAME(zgemm)('N','N',m-k,n-k,k,complex_double_one_,L.addr(k,0),m,
        U.addr(0,k),k,complex_double_zero_,P->addr(k,k),m); 
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const UpperTrapezoidalMatrix<double,complex<double> > &U) {
// Note that
// M [ U_1 , U_2 ] = [ M U_1 , M U_2 ]
  char diagU=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U) ==0 ?
    'N' : 'U');
  int m=M.size(0),k=U.size(0),n=U.size(1);
  CHECK_SAME(k,M.size(1));
  Matrix<double,complex<double> > *P=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zlacpy)('A',m,k,M.addr(),m,P->addr(),m);
  F77NAME(ztrmm)('R','U','N',diagU,m,k,complex_double_one_,U.addr(),k,
    P->addr(),m); // M U_1
  if (n>k) { // M U_2
    F77NAME(zgemm)('N','N',m,n-k,k,complex_double_one_,M.addr(),m,
      U.addr(0,k),k,complex_double_zero_,P->addr(0,k),m);
  }
  return P;
}

template class UpperTrapezoidalMatrix<double,complex<double> >;
template void testUpperTrapezoidalMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double
UpperTriangularMatrix<double,complex<double> >::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info;
  F77NAME(ztrcon)(norm,'U','N',n,addr(),n,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork;
  delete [] work;
  return rcond;
}

/*
template<> UpperTriangularMatrix<double,complex<double> >*
UpperTriangularMatrix<double,complex<double> >::inverse() const {
  UpperTriangularMatrix<double,complex<double> > *I=
    OPERATOR_NEW UpperTriangularMatrix<double,complex<double> >(*this);
  int n=size(0);
  int info;
  F77NAME(ztrtri)('U','N',n,I->addr(),n,info);
  CHECK_TEST(info==0);
  return I;
}
*/

template class UpperTriangularMatrix<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void
UnitUpperTrapezoidalMatrix<double,complex<double> >::copyFrom(int m,int n,
const Matrix<double,complex<double> > &U) {
  m=min(m,min(size(0),U.size(0)));
  n=min(n,min(size(1),U.size(1)));
  for (int j=1;j<n;j++) {
    F77NAME(zcopy)(min(j,m),U.addr(0,j),1,addr(0,j),1);
  }
}

template<> UpperTrapezoidalMatrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<double,complex<double> > *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(m,n);
  if (U_non_unit) {
    S->copy(U);
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j,m),complex_double_one_,addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=complex_double_one_;
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zcopy)(min(j,m),addr(0,j),1,S->addr(0,j),1);
        F77NAME(zaxpy)(min(j,m),complex_double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)=2.;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(m,complex_double_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zcopy)(j,addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)=complex_double_one_;
      F77NAME(zaxpy)(m-j,complex_double_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zcopy)(j,addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)=2.;
      if (m<j-1) {
        F77NAME(zaxpy)(m-j-1,complex_double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *S=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(zaxpy)(min(j,m),complex_double_one_,addr(0,j),1,
        S->addr(0,j),1);
    }
    if (j<m) (*S)(j,j)+=complex_double_one_;
  }
  return S;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<double,complex<double> > *D=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(m,n);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zcopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)=complex_double_one_;
      F77NAME(zaxpy)(min(j+1,m),complex_double_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zcopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
        F77NAME(zaxpy)(min(j,m),complex_double_mone_,addr(0,j),1,
          D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)=complex_double_zero_;
    }
  }
  return D;
}

template<> SquareMatrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *D=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(m,complex_double_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zcopy)(j,addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)=complex_double_one_;
      F77NAME(zaxpy)(m-j,complex_double_mone_,L.addr(j,j),1,
        D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zcopy)(j,addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)=complex_double_zero_;
      if (m<j-1) {
        F77NAME(zaxpy)(m-j-1,complex_double_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
    }
  }
  return D;
}

template<> Matrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *D=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(zcopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
    }
    if (j<m) (*D)(j,j)=complex_double_one_;
    F77NAME(zaxpy)(m,complex_double_one_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator*(
complex<double> d) const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<double,complex<double> > *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(m,n);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(zcopy)(min(j,m),addr(0,j),1,P->addr(0,j),1);
      F77NAME(zscal)(min(j,m),d,P->addr(0,j),1);
    }
    if (j<m) (*P)(j,j)=d;
  }
  return P;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator/(
complex<double> d) const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<double,complex<double> > *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(m,n);
  complex<double> dinv=complex_double_one_/d;
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(zcopy)(min(j,m),addr(0,j),1,P->addr(0,j),1);
      F77NAME(zscal)(min(j,m),dinv,P->addr(0,j),1);
    }
    if (j<m) (*P)(j,j)=dinv;
  }
  return P;
}

template<> UpperTrapezoidalMatrix<double,complex<double> >* 
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
// [ U_1 , U_2 ] [ V_11 , V_12 , V_13 ]
//               [   0  , V_22 , V_23 ]
//   = [ U_1 V_11 , U_1 V_12 + U_2 V_22 , U_1 V_13 + U_2 V_23 ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  int m=size(0),k=size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  UpperTrapezoidalMatrix<double,complex<double> > *P=OPERATOR_NEW
    UpperTrapezoidalMatrix<double,complex<double> >(m,n,
      complex_double_zero_);
  if (U_non_unit) {
    for (int j=0;j<m;j++) {
      F77NAME(zcopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ztrmv)('U','N','U',j+1,addr(),m,P->addr(0,j),1);
    }
    for (int j=m;j<k;j++) {
      F77NAME(zcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ztrmv)('U','N','U',m,addr(),m,P->addr(0,j),1);
      F77NAME(zgemv)('N',m,j-m+1,complex_double_one_,addr(0,m),m,
        U.addr(m,j),1,complex_double_one_,P->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(zcopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(ztrmv)('U','N','U',j,addr(),m,P->addr(0,j),1);
        F77NAME(zaxpy)(j,complex_double_one_,addr(0,j),1,P->addr(0,j),1);
      }
      (*P)(j,j)=complex_double_one_;
    }
    for (int j=m;j<k;j++) {
      F77NAME(zcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ztrmv)('U','N','U',m,addr(),m,P->addr(0,j),1);
      if (j>m) {
        F77NAME(zgemv)('N',m,j-m,complex_double_one_,addr(0,m),m,
          U.addr(m,j),1,complex_double_one_,P->addr(0,j),1);
      }
      F77NAME(zaxpy)(m,complex_double_one_,addr(0,j),1,P->addr(0,j),1);
    }
  }
  for (int j=k;j<n;j++) {
    F77NAME(zcopy)(m,U.addr(0,j),1,P->addr(0,j),1);
    F77NAME(ztrmv)('U','N','U',m,addr(),m,P->addr(0,j),1);
    F77NAME(zgemv)('N',m,k-m,complex_double_one_,addr(0,m),m,
      U.addr(m,j),1,complex_double_one_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  if (m>=n) {
//  note that
//  [ U_11 U_12 U_13 ] [ L_1 ] = [ U_11 L_1 + U_12 L_2 + U_13 L_3 ]
//  [      U_22 U_23 ] [ L_2 ] = [            U_22 L_2 + U_23 L_3 ]
//                     [ L_3 ]
    if (L_non_unit) { // U_11 L_1:
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(zgemv)('N',j,n-j,complex_double_one_,addr(0,j),m,
            L.addr(j,j),1,complex_double_zero_,P->addr(0,j),1);
        }
        F77NAME(zcopy)(n-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(ztrmv)('U','N','U',n-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(zaxpy)(j,complex_double_one_,addr(0,j),1,
            P->addr(0,j),1);
          if (j<n-1) {
            F77NAME(zgemv)('N',j,n-j-1,complex_double_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,complex_double_zero_,P->addr(0,j),1);
          }
        }
        if (j<n-1) {
          F77NAME(zcopy)(n-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(ztrmv)('U','N','U',n-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(zdotu)(n-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=complex_double_one_;
      }
    }
    if (n<k) { // U_12 L_2 + U_13 L_3:
      F77NAME(zgemm)('N','N',n,n,k-n,double_one_,addr(0,n),m,
        L.addr(n,0),k,double_one_,P->addr(0,0),m);
    }
    if (n<m) {
      F77NAME(zlacpy)('A',m-n,n,L.addr(n,0),k,P->addr(n,0),m);
      F77NAME(ztrmm)('L','U','N','U',m-n,n,complex_double_one_,
        addr(n,n),m,P->addr(n,0),m); // U_22 L_2
      if (m<k) {
        F77NAME(zgemm)('N','N',m-n,n,k-m,complex_double_one_,addr(n,m),m,
          L.addr(m,0),k,complex_double_one_,P->addr(n,0),m); // U_23 L_3
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
          F77NAME(zgemv)('N',j,m-j,complex_double_one_,addr(0,j),m,
            L.addr(j,j),1,complex_double_zero_,P->addr(0,j),1);
        }
        F77NAME(zcopy)(m-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(ztrmv)('U','N','U',m-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(zaxpy)(j,complex_double_one_,addr(0,j),1,
            P->addr(0,j),1);
          if (j<m-1) {
            F77NAME(zgemv)('N',j,m-j-1,complex_double_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,complex_double_zero_,P->addr(0,j),1);
          }
        }
        if (j<m-1) {
          F77NAME(zcopy)(m-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(ztrmv)('U','N','U',m-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(zdotu)(m-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (m<k) { // U_2 L_21 + U_3 L_31 
      F77NAME(zgemm)('N','N',m,m,k-m,complex_double_one_,addr(0,m),m,
        L.addr(m,0),k,complex_double_one_,P->addr(0,0),m);
    }
    char diagL=(L_non_unit ? 'N' : 'U');
    F77NAME(zlacpy)('A',m,n-m,addr(0,m),m,P->addr(0,m),m);
    F77NAME(ztrmm)('R','L','U',diagL,m,n-m,complex_double_one_,
      L.addr(m,m),k,P->addr(0,m),m); // U_2 L_22
    if (n<k) {
      F77NAME(zgemm)('N','N',m,n-m,k-n,complex_double_one_,addr(0,n),m,
        L.addr(n,m),k,complex_double_one_,P->addr(0,m),m); // U_3 L_32
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
// [ U_1 U_2 ] [ M_1 ] = U_1 M_1 + U_2 M_2
//             [ M_2 ]
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<double,complex<double> > *P=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zlacpy)('A',m,n,M.addr(),k,P->addr(),m);
  F77NAME(ztrmm)('L','U','N','U',m,n,complex_double_one_,addr(),m,
    P->addr(),m); // U_1 M_1
  if (k>m) {
    F77NAME(zgemm)('N','N',m,n,k-m,complex_double_one_,addr(0,m),m,
      M.addr(m,0),k,complex_double_one_,P->addr(),m); // U_2 M_2
  }
  return P;
}

template<> Vector<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<double,complex<double> > *p=
    OPERATOR_NEW Vector<double,complex<double> >(m);
  F77NAME(zcopy)(m,v.addr(),1,p->addr(),1);
  F77NAME(ztrmv)('U','N','U',m,addr(),m,p->addr(),1);
  F77NAME(zgemv)('N',m,n-m,complex_double_one_,addr(0,m),m,v.addr(m),1,
    complex_double_one_,p->addr(),1);
  return p;
}

/*
template<> UnitLowerTrapezoidalMatrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::transpose() const {
  int m=size(0),n=size(1);
  UnitLowerTrapezoidalMatrix<double,complex<double> > *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<double,complex<double> >(n,m);
  for (int i=0;i<m;i++) {
    F77NAME(zcopy)(n-i-1,addr(i,i+1),m,L->addr(i+1,i),1);
  }
  return L;
}

template<> UnitLowerTrapezoidalMatrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::conjugateTranspose()
const {
  int m=size(0),n=size(1);
  UnitLowerTrapezoidalMatrix<double,complex<double> > *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<double,complex<double> >(n,m);
  for (int i=0;i<m;i++) {
    const complex<double> *row_i=addr(i,i+1);
    complex<double> *L_col_i=L->addr(i+1,i);
    for (int j=i+1;j<n;j++,row_i+=m,L_col_i++) {
      *L_col_i=conj(*row_i);
    }
  }
  return L;
}
*/

template<> Vector<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::trmv(
const Vector<double,complex<double> > &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<double,complex<double> > *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<double,complex<double> >(n);
    F77NAME(zcopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(ztrmv)('U',trans,'U',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(zgemv)(trans,m,n-m,one_,addr(0,m),m,x.addr(),1,zero_,
        p->addr(m),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<double,complex<double> >(m);
    F77NAME(zcopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(ztrmv)('U','N','U',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(zgemv)('N',m,n-m,one_,addr(0,m),m,x.addr(m),1,one_,
        p->addr(),1);
    }
  }
  return p;
}

template<> Matrix<double,complex<double> >*
UnitUpperTrapezoidalMatrix<double,complex<double> >::trmm(
const Matrix<double,complex<double> > &M,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<double,complex<double> > *P=0;
  if (side=='L' || side=='l') {
    int k=M.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,M.size(0));
      P=OPERATOR_NEW Matrix<double,complex<double> >(n,k);
      F77NAME(zlacpy)('A',m,k,M.addr(),m,P->addr(),n);
      F77NAME(ztrmm)('L','U',trans,'U',m,k,one_,addr(),m,P->addr(),n);
      if (n>m) {
        F77NAME(zgemm)(trans,'N',n-m,k,m,one_,addr(0,m),m,M.addr(),m,
          zero_,P->addr(m,0),n);
      }
    } else {
      CHECK_SAME(n,M.size(0));
      P=OPERATOR_NEW Matrix<double,complex<double> >(m,k);
      F77NAME(zlacpy)('A',m,k,M.addr(),n,P->addr(),m);
      F77NAME(ztrmm)('L','U','N','U',m,k,one_,addr(),m,P->addr(),m);
      if (n>m) {
        F77NAME(zgemm)('N','N',m,k,n-m,one_,addr(0,m),m,M.addr(m,0),n,
          one_,P->addr(),m);
      }
    }
  } else {
    int k=M.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,M.size(1));
      P=OPERATOR_NEW Matrix<double,complex<double> >(k,m);
      F77NAME(zlacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(ztrmm)('R','U',trans,'U',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(zgemm)('N',trans,k,m,n-m,one_,M.addr(0,m),k,addr(0,m),m,
          one_,P->addr(),k);
      }
    } else {
      CHECK_SAME(m,M.size(1));
      P=OPERATOR_NEW Matrix<double,complex<double> >(k,n);
      F77NAME(zlacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(ztrmm)('R','U','N','U',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(zgemm)('N','N',k,n-m,m,one_,M.addr(),k,addr(0,m),m,
          zero_,P->addr(0,m),k);
      }
    }
  }
  return P;
}

template<> double
UnitUpperTrapezoidalMatrix<double,complex<double> >::normFrobenius()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('F','U','U',m,size(1),addr(),m,work);
}

template<> double
UnitUpperTrapezoidalMatrix<double,complex<double> >::normInfinity()
const {
  int m=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,m);
  double result=F77NAME(zlantr)('I','U','U',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> double
UnitUpperTrapezoidalMatrix<double,complex<double> >::normMaxEntry()
const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('M','U','U',m,size(1),addr(),m,work);
}

template<> double
UnitUpperTrapezoidalMatrix<double,complex<double> >::normOne() const {
  int m=size(0);
  double *work=0;
  return F77NAME(zlantr)('O','U','U',m,size(1),addr(),m,work);
}

template<> void
UnitUpperTrapezoidalMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    // calling routine will have to check consistency conditions
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(zcopy)(m,b.addr(),1,x.addr(),1);
    F77NAME(ztrsv)('U',trans,'U',m,addr(),m,x.addr(),1);
  } else {
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(zcopy)(m,b.addr(),1,x.addr(),1);
    if (n>m) { // use trailing entries of x as free variables
      F77NAME(zgemv)('N',m,n-m,mone_,addr(0,m),m,x.addr(m),1,one_,
        x.addr(),1);
    }
    F77NAME(ztrsv)('U','N','U',m,addr(),m,x.addr(),1);
  }
}

template<> void
UnitUpperTrapezoidalMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<double,complex<double> >*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(zlacpy)('A',m,k,B.addr(),n,X.addr(),m);
      F77NAME(ztrsm)('L','U',trans,'U',m,k,one_,addr(),m,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(zlacpy)('A',m,k,B.addr(),m,X.addr(),n);
      if (n>m) { // use these entries of X as free variables
        F77NAME(zgemm)('N','N',m,k,n-m,mone_,addr(0,m),m,X.addr(m,0),n,
          one_,X.addr(),n);
      }
      F77NAME(ztrsm)('L','U','N','U',m,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(zlacpy)('A',k,m,B.addr(),k,X.addr(),k);
      if (n>m) { // use these entries of X as free variables
        F77NAME(zgemm)('N',trans,k,m,n-m,mone_,X.addr(0,m),k,addr(0,m),m,
          one_,X.addr(),k);
      }
      F77NAME(ztrsm)('R','U',trans,'U',k,m,one_,addr(),m,X.addr(),k);
    } else { // calling routine will have to check consistency
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(zlacpy)('A',k,m,B.addr(),k,X.addr(),k);
      F77NAME(ztrsm)('R','U','N','U',k,m,one_,addr(),m,X.addr(),k);
    }
  }
}

template class UnitUpperTrapezoidalMatrix<double,complex<double> >;
template void testUnitUpperTrapezoidalMatrix(double,complex<double> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> double
UnitUpperTriangularMatrix<double,complex<double> >::
reciprocalConditionNumber(char norm) const {
  int n=size(0);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info;
  F77NAME(ztrcon)(norm,'U','U',n,addr(),n,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork;
  delete [] work;
  return rcond;
}

/*
template<> UnitUpperTriangularMatrix<double,complex<double> >*
UnitUpperTriangularMatrix<double,complex<double> >::inverse() const {
  UnitUpperTriangularMatrix<double,complex<double> > *I=
    OPERATOR_NEW UnitUpperTriangularMatrix<double,complex<double> >(*this);
  int n=size(0);
  int info;
  F77NAME(ztrtri)('U','U',n,I->addr(),n,info);
  CHECK_TEST(info==0);
  return I;
}
*/

template class UnitUpperTriangularMatrix<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "OrthogonalMatrix.C"
/*
template<> OrthogonalMatrix<double,complex<double> >*
OrthogonalMatrix<double,complex<double> >::transpose() const {
  int m=size(0),n=size(1);
  OrthogonalMatrix<double,complex<double> > *Q=
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(m,addr(0,j),1,Q->addr(j,0),n);
  }
  return Q;
}

template<> OrthogonalMatrix<double,complex<double> >*
OrthogonalMatrix<double,complex<double> >::conjugateTranspose() const {
  int m=size(0),n=size(1);
  OrthogonalMatrix<double,complex<double> > *Q=
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(n,m);
  for (int j=0;j<n;j++) {
    const complex<double> *col_j=addr(0,j);
    complex<double> *Q_row_j=Q->addr(j,0);
    for (int i=0;i<m;i++,col_j++,Q_row_j+=n) {
      *Q_row_j=conj(*col_j);
    }
  }
  return Q;
}
*/

/*
template<> void OrthogonalMatrix<double,complex<double> >::solve(
const UpperTrapezoidalMatrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  int m=size(0), n=size(1);
  bool B_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&B)==0);
  X=complex_double_zero_;
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<=min(j,n-1);i++) {
            F77NAME(zaxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<min(j,n);i++) {
            F77NAME(zaxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
          if (j<n) {
            F77NAME(zaxpy)(m,complex_double_one_,addr(0,j),1,
              X.addr(0,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(zdotc)(min(j+1,m),addr(0,i),1,B.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(zdotc)(min(j,m),addr(0,i),1,B.addr(0,j),1);
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
            X(i,j)=F77NAME(zdotu)(m-i,B.addr(i,i),k,addr(i,j),1);
          }
        }
      } else {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=(*this)(i,j)
              +F77NAME(zdotu)(m-i-1,B.addr(i,i+1),k,addr(i+1,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      if (B_non_unit) {
        for (int i=0;i<k;i++) {
          for (int j=i;j<n;j++) {
            F77NAME(zaxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      } else {
        for (int i=0;i<k;i++) {
          F77NAME(zaxpy)(m,double_one_,addr(0,i),1,X.addr(i,0),k);
          for (int j=i+1;j<n;j++) {
            F77NAME(zaxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      }
    }
  }
}
*/

/*
template<> void OrthogonalMatrix<double,complex<double> >::solve(
const LowerTrapezoidalMatrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  int m=size(0), n=size(1);
  bool B_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&B)==0);
  X=complex_double_zero_;
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=j;i<n;i++) {
            F77NAME(zaxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          F77NAME(zaxpy)(m,complex_double_one_,addr(0,j),1,X.addr(0,j),1);
          for (int i=j+1;i<n;i++) {
            F77NAME(zaxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(zdotc)(m-j,addr(j,i),1,B.addr(j,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=(*this)(j,i)
              +F77NAME(zdotc)(m-j-1,addr(j+1,i),1,B.addr(j+1,j),1);
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
            X(i,j)=F77NAME(zdotu)(min(i+1,m),B.addr(i,0),k,addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=F77NAME(zdotu)(min(i,m),B.addr(i,0),k,addr(0,j),1);
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
            F77NAME(zaxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      } else {
        for (int i=0;i<k;i++) {
          for (int j=0;j<min(i,n);j++) {
            F77NAME(zaxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
          if (i<n) {
            F77NAME(zaxpy)(m,complex_double_one_,addr(0,i),1,
              X.addr(i,0),k);
          }
        }
      }
    }
  }
}
*/

template<> void OrthogonalMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  int m=size(0), n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans=='C' || trans=='c') { // X=Q*B is smallest solution
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(zgemm)('N','N',m,k,n,double_one_,addr(),m,B.addr(),n,
        complex_double_zero_,X.addr(),m);
    } else if (trans=='T' || trans=='t') {
    // X=conj(Q)*B=conj(Q*conj(B)) is smallest
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      X=complex_double_zero_;
      for (int kk=0;kk<k;kk++) {
        for (int j=0;j<n;j++) {
          F77NAME(zaxpy)(m,conj(B(j,kk)),addr(0,j),1,X.addr(0,kk),1);
        }
        complex<double> *Xik=X.addr(0,kk);
        for (int i=0;i<m;i++,Xik++) *Xik=conj(*Xik);
      }
    } else { // X=Q^H*B minimizes residual
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
//    F77NAME(zgemm)('C','N',n,k,m,double_one_,addr(),m,B.addr(),m,
//      complex_double_zero_,X.addr(),n);
      Vector<double,complex<double> > *r=
        OPERATOR_NEW Vector<double,complex<double> >(m);
      for (int l=0;l<k;l++) {
        F77NAME(zcopy)(m,B.addr(0,l),1,r->addr(),1);
        for (int j=0;j<n;j++) {
          complex<double> &Xjl=X(j,l);
          Xjl=F77NAME(zdotc)(m,addr(0,j),1,r->addr(),1);
          F77NAME(zaxpy)(m,-Xjl,addr(0,j),1,r->addr(),1);
        }
      }
      delete r; r=0;
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans=='C' || trans=='c') { // X=B*Q minimizes residual
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
//    F77NAME(zgemm)('N','N',k,n,m,double_one_,B.addr(),k,addr(),m,
//      complex_double_zero_,X.addr(),k);
      Vector<double,complex<double> > *r=
        OPERATOR_NEW Vector<double,complex<double> >(m);
      for (int l=0;l<k;l++) {
        for (int i=0;i<m;i++) (*r)[i]=conj(B(l,i));
        for (int j=0;j<n;j++) {
          complex<double> Xlj=F77NAME(zdotc)(m,addr(0,j),1,r->addr(),1);
          F77NAME(zaxpy)(m,-Xlj,addr(0,j),1,r->addr(),1);
          X(l,j)=conj(Xlj);
        }
      }
      delete r; r=0;
    } else if (trans=='T' || trans=='t') {
    // X=B*conj(Q) minimizes resid
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
//    X=complex_double_zero_;
//    for (int j=0;j<n;j++) {
//      for (int i=0;i<m;i++) {
//        F77NAME(zaxpy)(k,conj((*this)(i,j)),B.addr(0,i),1,
//          X.addr(0,j),1);
//      }
//    }
      Vector<double,complex<double> > *r=
        OPERATOR_NEW Vector<double,complex<double> >(m);
      for (int l=0;l<k;l++) {
        F77NAME(zcopy)(m,B.addr(l,0),k,r->addr(),1);
        for (int j=0;j<n;j++) {
          complex<double> &Xlj=X(l,j);
          Xlj=F77NAME(zdotc)(m,addr(0,j),1,r->addr(),1);
          F77NAME(zaxpy)(m,-Xlj,addr(0,j),1,r->addr(),1);
        }
      }
      delete r; r=0;
    } else { // X=B*Q^H is smallest solution
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(zgemm)('N','C',k,m,n,double_one_,B.addr(),k,addr(),m,
        complex_double_zero_,X.addr(),k);
    }
  }
}

template<> void OrthogonalMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  int m=size(0), n=size(1);
  if (trans=='C' || trans=='c') { // x=Q*b is smallest solution
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(zgemv)('N',m,n,double_one_,addr(),m,b.addr(),1,
      complex_double_zero_,x.addr(),1);
  } else if (trans=='T' || trans=='t') {
  // x=conj(Q)*b=conj(Q*conj(b)) is smallest
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    x=complex_double_zero_;
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(m,conj(b[j]),addr(0,j),1,x.addr(),1);
    }
    complex<double> *xi=x.addr();
    for (int i=0;i<m;i++,xi++) *xi=conj(*xi);
  } else { // x=Q^H*b minimizes residual
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(zgemv)('C',m,n,double_one_,addr(),m,b.addr(),1,
      complex_double_zero_,x.addr(),1);
  }
}

template class OrthogonalMatrix<double,complex<double> >;
template void testOrthogonalMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SymmetricMatrix.C"
template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,addr(j,j),1,M->addr(j,j),1);
    if (j+1<n) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Mji=M->addr(j,j+1);
      for (int i=j+1;i<n;i++,Aij++,Mji+=n) *Mji=conj(*Aij);
    }
  }
  return M;
}

template<> void SymmetricMatrix<double,complex<double> >::fillWith(
complex<double> d) {
  Matrix<double,complex<double> >::set('L',d,d.real());
}

template<> complex<double>
SymmetricMatrix<double,complex<double> >::operator()(int i,int j) const {
  return (j<=i ? *(this->addr(i,j)) : conj(*(this->addr(j,i))) );
}

template<> SymmetricMatrix<double,complex<double> >&
SymmetricMatrix<double,complex<double> >::operator+=(
const SymmetricMatrix<double,complex<double> > &S) {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0))
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_one_,S.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<double,complex<double> >&
SymmetricMatrix<double,complex<double> >::operator-=(
const SymmetricMatrix<double,complex<double> > &S) {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0))
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_mone_,S.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<double,complex<double> >&
SymmetricMatrix<double,complex<double> >::operator*=(
complex<double> scalar) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    F77NAME(zdscal)(n-j,scalar.real(),addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<double,complex<double> >&
SymmetricMatrix<double,complex<double> >::operator/=(
complex<double> scalar) {
  int n=this->size(0);
  CHECK_NONZERO(abs(scalar))
  for (int j=0;j<n;j++) {
    F77NAME(zdscal)(n-j,double_one_/scalar.real(),addr(j,j),1);
  }
  return *this;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  int n=this->size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *col_j=addr(j+1,j);
      complex<double> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(min(j+1,n),one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j,n),one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      if (j<n) (*S)(j,j)+=one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int n=this->size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *col_j=addr(j+1,j);
      complex<double> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(n-j,one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(zaxpy)(n-j-1,one_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      (*S)(j,j)+=one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator+(
const SquareMatrix<double,complex<double> > &S) const {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *T=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  T->copy(S);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,one_,addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *col_j=addr(j+1,j);
      complex<double> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,T_row_j+=n) *T_row_j+=conj(*col_j);
    }
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  int n=this->size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,one_,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *col_j=addr(j+1,j);
      complex<double> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j+=conj(*col_j);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  int n=this->size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *col_j=addr(j+1,j);
      complex<double> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(min(j+1,n),mone_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j,n),mone_,U.addr(0,j),1,S->addr(0,j),1);
      }
      if (j<n) (*S)(j,j)+=mone_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const UnitUpperTrapezoidalMatrix<double,complex<double> > &U,
const SymmetricMatrix<double,complex<double> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,complex<double> > *T=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(zcopy)(min(j,n),U.addr(0,j),1,T->addr(0,j),1);
    }
    if (j<n) (*T)(j,j)=complex_double_one_;
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *S_col_j=S.addr(j+1,j);
      complex<double> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const SymmetricMatrix<double,complex<double> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,complex<double> > *T=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+1,n),U.addr(0,j),1,T->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zcopy)(min(j,n),U.addr(0,j),1,T->addr(0,j),1);
      }
      if (j<n) (*T)(j,j)=complex_double_one_;
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *S_col_j=S.addr(j+1,j);
      complex<double> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int n=this->size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *col_j=addr(j+1,j);
      complex<double> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(n-j,mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(zaxpy)(n-j-1,mone_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      (*S)(j,j)+=mone_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const UnitLowerTrapezoidalMatrix<double,complex<double> > &L,
const SymmetricMatrix<double,complex<double> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *T=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    if (j<n-1) {
      F77NAME(zcopy)(n-j-1,L.addr(j+1,j),1,T->addr(j+1,j),1);
    }
    (*T)(j,j)=complex_double_one_;
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *S_col_j=S.addr(j+1,j);
      complex<double> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const SymmetricMatrix<double,complex<double> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *T=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(n-j,L.addr(j,j),1,T->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(zcopy)(n-j-1,L.addr(j+1,j),1,T->addr(j+1,j),1);
      }
      (*T)(j,j)=complex_double_one_;
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *S_col_j=S.addr(j+1,j);
      complex<double> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator-(
const SquareMatrix<double,complex<double> > &S) const {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *T=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *col_j=addr(j+1,j);
      complex<double> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,T_row_j+=n) *T_row_j=conj(*col_j);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n,mone_,S.addr(0,j),1,T->addr(0,j),1);
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const SquareMatrix<double,complex<double> > &S,
const SymmetricMatrix<double,complex<double> > &SS) {
  int n=SS.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *T=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  T->copy(S);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_mone_,SS.addr(j,j),1,
      T->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *S_col_j=S.addr(j+1,j);
      complex<double> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  int n=this->size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *col_j=addr(j+1,j);
      complex<double> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n,mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const SymmetricMatrix<double,complex<double> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *T=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  T->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *S_col_j=S.addr(j+1,j);
      complex<double> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator*(complex<double> d)
const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    const complex<double> *col_j=addr(j,j);
    complex<double> *S_col_j=S->addr(j,j);
    *S_col_j=(*col_j)*d;
    if (j<n-1) {
      S_col_j++;
      col_j++;
      complex<double> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_col_j++,S_row_j+=n) {
        *S_col_j=(*col_j)*d;
        *S_row_j=conj(*col_j)*d;
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator/(complex<double> d)
const {
  CHECK_NONZERO(abs(d));
  return operator*(complex_double_one_/d);
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator*(
const SymmetricMatrix<double,complex<double> > &S) const {
// compute by bordering: note that
// [ sigma s^H ] [ tau t^H ] = [ sigma tau + s^H t , sigma t^H + s^H T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^H +  S  T ]
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *T=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-1;k>=0;k--) {
    if (k<n-1) {
      F77NAME(zhemv)('L',n-k-1,complex_double_one_,addr(k+1,k+1),n,
        S.addr(k+1,k),1,complex_double_zero_,T->addr(k+1,k),1); // S t
      complex<double> *tt=t+(k+1);
      F77NAME(zhemv)('L',n-k-1,complex_double_one_,S.addr(k+1,k+1),n,
        addr(k+1,k),1,complex_double_zero_,tt,1); // s^H T = ( T s )^H
      complex<double> *tk=T->addr(k,k+1);
      for (int i=k+1;i<n;i++,tk+=n,tt++) *tk=conj(*tt);
      (*T)(k,k)=F77NAME(zdotc)(n-k-1,addr(k+1,k),1,S.addr(k+1,k),1);//s^Ht
    }
    F77NAME(zgerc)(n-k,n-k,complex_double_one_,addr(k,k),1,S.addr(k,k),1,
      T->addr(k,k),n);
  }
  delete [] t; t=0;
  return T;
}

template<> Matrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
// compute by bordering: note that
// S [ U_1 , U_2 ] = [ S U_1 , S U_2 ]
// and that
// [ sigma s^H ] [ upsilon u^T ] = [ sigma upsilon , sigma u^T + s^H U ]
// [   s    S  ] [          U  ] [ [   s   upsilon ,   s   u^T +  S  U ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  char diag=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0
    ? 'N' : 'U');
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,m);
  for (int k=m-1;k>=0;k--) { // S U_1
    if (k<m-1) {
      F77NAME(zcopy)(m-k-1,addr(k+1,k),1,t+(k+1),1);
      complex<double> *tt=t+(k+1);
      F77NAME(ztrmv)('U','C',diag,m-k-1,U.addr(k+1,k+1),m,tt,1);
      complex<double> *mk=M->addr(k,k+1);
      for (int j=k+1;j<m;j++,mk+=m,tt++) *mk=conj(*tt); // s^H U
    }
    if (diag=='N') {
      F77NAME(zgeru)(m-k,m-k,complex_double_one_,addr(k,k),1,
        U.addr(k,k),m,M->addr(k,k),m);
    } else {
      F77NAME(zaxpy)(m-k,complex_double_one_,addr(k,k),1,M->addr(k,k),1);
      F77NAME(zgeru)(m-k,m-k-1,complex_double_one_,addr(k,k),1,
        U.addr(k,k+1),m,M->addr(k,k+1),m);
    }
  }
  if (n>m) { // S U_2
    F77NAME(zhemm)('L','L',m,n-m,complex_double_one_,addr(),m,
      U.addr(0,m),m,complex_double_zero_,M->addr(0,m),m);
  }
  delete [] t; t=0;
  return M;
}

template<> Matrix<double,complex<double> >* operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const SymmetricMatrix<double,complex<double> > &S) {
//compute by bordering: note that
// [ U_1 U_2 ] [ S_11 S_21^H ]
//             [ S_21  S_22  ]
//   = [ U_1 S_11 + U_2 S_21 , U_1 S_21^H + U_2 S_22 ]
// and that
// [ U    u    ] [  S  bar(s)] = [ U S +   u    s^T , U bar(s) + u sigma ]
// [   upsilon ] [ s^T sigma ] = [      upsilon s^T ,    upsilon   sigma ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  char diag=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0
    ? 'N' : 'U');
  for (int k=0;k<m;k++) { // U_1 S_11
    if (k>0) {
      complex<double> *mk=M->addr(0,k);
      const complex<double> *sk=S.addr(k,0);
      for (int i=0;i<k;i++,mk++,sk+=n) *mk=conj(*sk);
      F77NAME(ztrmv)('U','N',diag,k,U.addr(),m,M->addr(0,k),1); //U bar(s)
    }
    if (diag=='N') {
      F77NAME(zgeru)(k+1,k+1,complex_double_one_,U.addr(0,k),1,
        S.addr(k,0),n,M->addr(),m);
    } else {
      F77NAME(zgeru)(k,k+1,complex_double_one_,U.addr(0,k),1,
        S.addr(k,0),n,M->addr(),m);
      F77NAME(zaxpy)(k+1,complex_double_one_,S.addr(k,0),n,
        M->addr(k,0),m);
    }
  }
  if (n>m) {
    F77NAME(zgemm)('N','N',m,m,n-m,complex_double_one_,U.addr(0,m),m,
      S.addr(m,0),n,complex_double_one_,M->addr(),m); // U_2 S_21
    for (int j=m;j<n;j++) {
      complex<double> *mj=M->addr(0,j);
      const complex<double> *sj=S.addr(j,0);
      for (int i=0;i<m;i++,mj++,sj+=n) *mj=conj(*sj);
    }
    F77NAME(ztrmm)('L','U','N',diag,m,n-m,complex_double_one_,U.addr(),m,
      M->addr(0,m),m); // U_1 S_21^T
    F77NAME(zhemm)('R','L',m,n-m,complex_double_one_,S.addr(m,m),n,
      U.addr(0,m),m,complex_double_one_,M->addr(0,m),m); // U_2 S_22
  }
  return M;
}

template<> Matrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
// compute by bordering: note that
// [ S_11 S_21^H ] [ L_1 ] = [ S_11 L_1 + S_21^H L_2 ]
// [ S_21  S_22  ] [ L_2 ] = [ S_21 L_1 +  S_22  L_2 ]
// and that
// [  S  bar(s)] [   L          ] = [  S  L +   s   ell^T ,bar(s) lambda ]
// [ s^T sigma ] [ ell^T lambda ] = [ s^T L + sigma ell^T , sigma lambda 
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(1));
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0
    ? 'N' : 'U');
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,m);
  for (int k=0;k<n;k++) { // S_11 L_1
    if (k>0) {
      F77NAME(zcopy)(k,addr(k,0),m,M->addr(k,0),m);
      F77NAME(ztrmv)('L','T',diag,k,L.addr(),m,M->addr(k,0),m); // s^T L
    }
    const complex<double> *sk=addr(k,0);
    complex<double> *tt=t;
    for (int i=0;i<=k;i++,sk+=m,tt++) *tt=conj(*sk);
    if (diag=='N') {
      F77NAME(zgeru)(k+1,k+1,complex_double_one_,t,1,L.addr(k,0),m,
        M->addr(),m);
    } else {
      F77NAME(zgeru)(k+1,k,complex_double_one_,t,1,L.addr(k,0),m,
        M->addr(),m);
      F77NAME(zaxpy)(k+1,complex_double_one_,t,1,M->addr(0,k),1);
    }
  }
  if (m>n) {
    F77NAME(zgemm)('C','N',n,n,m-n,complex_double_one_,addr(n,0),m,
      L.addr(n,0),m,complex_double_one_,M->addr(),m); // S_21^H L_2
    F77NAME(zlacpy)('A',m-n,n,addr(n,0),m,M->addr(n,0),m);
    F77NAME(ztrmm)('R','L','N',diag,m-n,n,complex_double_one_,L.addr(),m,
      M->addr(n,0),m); // S_21 L_1
    F77NAME(zhemm)('L','L',m-n,n,complex_double_one_,addr(n,n),m,
      L.addr(n,0),m,complex_double_one_,M->addr(n,0),m); // S_22 L_2
  }
  delete [] t; t=0;
  return M;
}

template<> Matrix<double,complex<double> >* operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const SymmetricMatrix<double,complex<double> > &S) {
// compute by bordering: note that
// [ L_1 ] S = [ L_1 S ]
// [ L_2 ]   = [ L_2 S ]
// and that
// [ lambda   ] [ sigma s^H ] = [ lambda sigma       , lambda s^H       ]
// [  ell   L ] [   s    S  ] = [  ell   sigma + L s ,   ell  s^H + L S ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0
    ? 'N' : 'U');
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-1;k>=0;k--) { // L_1 S
    if (k<n-1) {
      F77NAME(zcopy)(n-k-1,S.addr(k+1,k),1,M->addr(k+1,k),1);
      F77NAME(ztrmv)('L','N',diag,n-k-1,L.addr(k+1,k+1),m,
        M->addr(k+1,k),1); // L s
    }
    complex<double> *tk=t+k;
    F77NAME(zcopy)(n-k,S.addr(k,k),1,tk,1);
    if (diag=='N') {
      F77NAME(zgerc)(n-k,n-k,complex_double_one_,L.addr(k,k),1,
        tk,1,M->addr(k,k),m);
    } else {
      F77NAME(zgerc)(n-k-1,n-k,complex_double_one_,L.addr(k+1,k),1,
        tk,1,M->addr(k+1,k),m);
      for (int i=k;i<n;i++) t[i]=conj(t[i]);
      F77NAME(zaxpy)(n-k,complex_double_one_,tk,1,M->addr(k,k),m);
    }
  }
  if (m>n) { // L_2 S
    F77NAME(zhemm)('R','L',m-n,n,complex_double_one_,S.addr(),n,
      L.addr(n,0),m,complex_double_zero_,M->addr(n,0),m);
  } 
  delete [] t; t=0;
  return M;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator*(
const SquareMatrix<double,complex<double> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *T=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  F77NAME(zhemm)('L','L',n,n,complex_double_one_,addr(),n,S.addr(),n,
    complex_double_zero_,T->addr(),n);
  return T;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SquareMatrix<double,complex<double> > &A,
const SymmetricMatrix<double,complex<double> > &B) {
// compute by bordering: note that
// [ mu v^T ] [ sigma s^H ] = [ mu sigma + v^T s , mu s^H + v^T S ]
// [  m  M  ] [   s    S  ] = [  m sigma +  M  s ,  m s^H +  M  S ]
  int n=A.size(0);
  CHECK_SAME(n,B.size(0));
  SquareMatrix<double,complex<double> > *T=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  complex<double> *v=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-1;k>=0;k--) {
    if (k+1<n) {
      F77NAME(zgemv)('N',n-k-1,n-k-1,double_one_,A.addr(k+1,k+1),n,
        B.addr(k+1,k),1,double_zero_,T->addr(k+1,k),1); // M s
      complex<double> *ti=t;
      const complex<double> *Aik=A.addr(k,k+1);
      for (int i=k+1;i<n;i++,ti++,Aik+=n) *ti=conj(*Aik);
      F77NAME(zhemv)('L',n-k-1,double_one_,B.addr(k+1,k+1),n,
        t,1,double_zero_,v,1);
      complex<double> *vi=v;
      complex<double> *Tki=T->addr(k,k+1);
      for (int i=k+1;i<n;i++,vi++,Tki+=n) *Tki=conj(*vi);
      //v^T S = (S bar(v))^H
      (*T)(k,k)=
        F77NAME(zdotu)(n-k-1,A.addr(k,k+1),n,B.addr(k+1,k),1); // v^T s
    }
    F77NAME(zgerc)(n-k,n-k,double_one_,A.addr(k,k),1,B.addr(k,k),1,
      T->addr(k,k),n);
  }
  delete [] t; t=0;
  delete [] v; v=0;
  return T;
}

template<> Matrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  int m=size(0),n=M.size(1);
  CHECK_SAME(m,M.size(0));
  Matrix<double,complex<double> > *T=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zhemm)('L','L',m,n,complex_double_one_,addr(),m,M.addr(),m,
    complex_double_zero_,T->addr(),m);
  return T;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const SymmetricMatrix<double,complex<double> > &S) {
  int m=M.size(0),n=S.size(0);
  CHECK_SAME(n,M.size(1));
  Matrix<double,complex<double> > *T=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  F77NAME(zhemm)('R','L',m,n,complex_double_one_,S.addr(),n,M.addr(),m,
    complex_double_zero_,T->addr(),m);
  return T;
}

template<> Vector<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<double,complex<double> > *w=
    OPERATOR_NEW Vector<double,complex<double> >(n);
  F77NAME(zhemv)('L',n,complex_double_one_,addr(),n,v.addr(),1,
    complex_double_zero_,w->addr(),1);
  return w;
}

// y := A * x * alpha + y * beta
template<> void SymmetricMatrix<double,complex<double> >::symv(
complex<double> alpha,const Vector<double,complex<double> > &x,
complex<double> beta,Vector<double,complex<double> > &y) const {
  int n=this->size(0);
  F77NAME(zhemv)('L',n,alpha,addr(),n,x.addr(),1,beta,y.addr(),1);
}

// A += x * alpha * x^T
template<> void SymmetricMatrix<double,complex<double> >::syr(
complex<double> alpha,const Vector<double,complex<double> > &x) {
  int n=this->size(0);
  F77NAME(zher)('L',n,alpha,x.addr(),1,addr(),n);
}

// A += x * alpha * y^T + y * alpha * x^T
template<> void SymmetricMatrix<double,complex<double> >::syr2(
complex<double> alpha,const Vector<double,complex<double> > &x,
const Vector<double,complex<double> > &y) {
  int n=this->size(0);
  F77NAME(zher2)('L',n,alpha,x.addr(),1,y.addr(),1,addr(),n);
}

// C := A * alpha * B + C * beta
template<> void SymmetricMatrix<double,complex<double> >::symm(
complex<double> alpha,const Matrix<double,complex<double> > &B,
complex<double> beta,Matrix<double,complex<double> > &C,char side) const {
  int m=C.size(0),n=C.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,size(0));
    F77NAME(zhemm)('L','L',m,n,alpha,addr(),m,B.addr(),m,beta,
      C.addr(),m);
  } else {
    CHECK_SAME(n,size(0));
    F77NAME(zhemm)('R','L',m,n,alpha,addr(),n,B.addr(),m,beta,
      C.addr(),m);
  }
}

// C := A * alpha * A^T + C * beta
template<> void SymmetricMatrix<double,complex<double> >::syrk(
complex<double> alpha,const Matrix<double,complex<double> > &A,
complex<double> beta,char transa) {
  int m=A.size(0),n=A.size(1);
  if (transa!='N' && transa!='n') {
    CHECK_SAME(n,size(0));
    F77NAME(zherk)('L',transa,n,m,alpha,A.addr(),m,beta,addr(),n);
  } else {
    CHECK_SAME(m,size(0));
    F77NAME(zherk)('L',transa,m,n,alpha,A.addr(),m,beta,addr(),m);
  }
}

// C := A * alpha * B^T + B * alpha + A^T + C * beta
template<> void SymmetricMatrix<double,complex<double> >::syr2k(
complex<double> alpha,const Matrix<double,complex<double> > &A,
const Matrix<double,complex<double> > &B,complex<double> beta,
char transab) {
  int m=A.size(0),n=A.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (transab!='N' && transab!='n') {
    CHECK_SAME(n,size(0));
    F77NAME(zher2k)('L',transab,n,m,alpha,A.addr(),m,B.addr(),m,beta,
      addr(),n);
  } else {
    CHECK_SAME(m,size(0));
    F77NAME(zher2k)('L',transab,m,n,alpha,A.addr(),m,B.addr(),m,beta,
      addr(),m);
  }
}

/*
// y := abs(A) * abs(x) * alpha + abs(y) * beta
template<> void SymmetricMatrix<double,complex<double> >::syamv(double alpha,
const Vector<double,complex<double> > &x,double beta,Vector<double,complex<double> > &y)
const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,y.size());
  F77_NAME(dla_syamv)('L',n,alpha,addr(),n,x.addr(),1,beta,y.addr(),1);
}
*/

template<> double SymmetricMatrix<double,complex<double> >::equilibrate(
Vector<double,double> &s,double &scond) const {
  int n=size(0);
  CHECK_SAME(n,s.size());
  double amax=numeric_limits<double>::infinity();
  int info;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,3*n);
  F77NAME(zheequb)('L',n,addr(),n,s.addr(),scond,amax,work,info);
  CHECK_SAME(info,0);
  delete [] work;
  return amax;
}

template<> double
SymmetricMatrix<double,complex<double> >::normFrobenius() const {
  int n=size(0);
  double *work=0;
  return F77NAME(zlanhe)('F','L',n,addr(),n,work);
}

template<> double
SymmetricMatrix<double,complex<double> >::normInfinity() const {
  int n=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,n);
  double val=F77NAME(zlanhe)('I','L',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> double
SymmetricMatrix<double,complex<double> >::normMaxEntry() const {
  int n=size(0);
  double *work=0;
  return F77NAME(zlanhe)('M','L',n,addr(),n,work);
}

template<> double
SymmetricMatrix<double,complex<double> >::normOne() const {
  int n=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,n);
  double val=F77NAME(zlanhe)('O','L',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> double
SymmetricMatrix<double,complex<double> >::reciprocalConditionNumber()
const {
  int n=size(0);

  SymmetricMatrix<double,complex<double> > *AF=
    OPERATOR_NEW SymmetricMatrix<double,complex<double> >(*this);
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  complex<double> w=complex_double_undefined_;
  F77NAME(zhetrf)('L',n,AF->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0);

  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zhetrf)('L',n,AF->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;

  double *rwork=OPERATOR_NEW_BRACKET(double,2*n);
  double anorm=F77NAME(zlanhe)('O','L',n,addr(),n,rwork);
  delete rwork; rwork=0;

  double rcond;
  work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  F77NAME(zhecon)('L',n,AF->addr(),n,ipiv,anorm,rcond,work,info);
  delete [] ipiv; ipiv=0;
  delete [] work; work=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SymmetricMatrix<double,complex<double> >*
SymmetricMatrix<double,complex<double> >::inverse() const {
  int n=size(0);
  SymmetricMatrix<double,complex<double> > *Ainv=
    OPERATOR_NEW SymmetricMatrix<double,complex<double> >(*this);
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  complex<double> w=undefined_;
  F77NAME(zhetrf)('L',n,Ainv->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0);

  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zhetrf)('L',n,Ainv->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;

  work=OPERATOR_NEW_BRACKET(complex<double>,n);
  F77NAME(zhetri)('L',n,Ainv->addr(),n,ipiv,work,info);
  CHECK_SAME(info,0)

  delete [] ipiv;
  delete [] work;
  return Ainv;
}
*/

template<> Vector<double,double>*
SymmetricMatrix<double,complex<double> >::eigenvalues(
OrthogonalMatrix<double,complex<double> > *&Q) const {
  int n=size(0);
  if (Q!=0) CHECK_SAME(n,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Q->Matrix<double,complex<double> >::copyFrom('L',n,n,*this);
  Vector<double,double> *lambda=OPERATOR_NEW Vector<double,double>(n);
  complex<double> w;
  int lwork=-1,info;
  double *rwork=OPERATOR_NEW_BRACKET(double,max(1,3*n-2));
  F77NAME(zheev)(jobz,'L',n,Q->addr(),n,lambda->addr(),&w,lwork,rwork,
    info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zheev)(jobz,'L',n,Q->addr(),n,lambda->addr(),work,lwork,rwork,
    info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete [] rwork; rwork=0;

  return lambda;
}

template<> void SymmetricMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char) const {
  int n=size(0);
  SymmetricMatrix<double,complex<double> > AF(*this);
  int info;
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  x.copy(b);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  complex<double> w=complex_double_undefined_;
  F77NAME(zhetrf)('L',n,AF.addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zhetrf)('L',n,AF.addr(),n,ipiv,work,lwork,info);
  delete [] work; work=0;

  if (info==0) {
    F77NAME(zhetrs)('L',n,1,AF.addr(),n,ipiv,x.addr(),n,info);
  }
  CHECK_SAME(info,0)
  delete [] ipiv; ipiv=0;
}

template<> void SymmetricMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char) const {
  int n=size(0);
  SymmetricMatrix<double,complex<double> > AF(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  complex<double> w=complex_double_undefined_;
  int info;
  F77NAME(zhetrf)('L',n,AF.addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zhetrf)('L',n,AF.addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;

  bool left_side=(side=='L' || side=='l');
  if (left_side) {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1))
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    X.copy(B);
    F77NAME(zhetrs)('L',n,1,AF.addr(),n,ipiv,X.addr(),n,info);
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0))
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
    for (int i=0;i<k;i++) {
      complex<double> *tt=t;
      const complex<double> *bi=B.addr(i,0);
      for (int j=0;j<n;j++,tt++,bi+=k) *tt=conj(*bi);
      F77NAME(zhetrs)('L',n,1,AF.addr(),n,ipiv,t,n,info);
      CHECK_SAME(info,0)
      tt=t;
      complex<double> *xi=X.addr(i,0);
      for (int j=0;j<n;j++,tt++,xi+=k) *xi=conj(*tt);
    }
    delete [] t; t=0;
  }
  delete [] ipiv; ipiv=0;
}

template class SymmetricMatrix<double,complex<double> >;
template void testSymmetricMatrix(double,complex<double> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveMatrix<double,complex<double> >&
SymmetricPositiveMatrix<double,complex<double> >::operator*=(
complex<double> scalar) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    F77NAME(zdscal)(n-j,abs(scalar),addr(j,j),1);
  }
  return *this;
}

template<> SymmetricPositiveMatrix<double,complex<double> >&
SymmetricPositiveMatrix<double,complex<double> >::operator/=(
complex<double> scalar) {
  int n=this->size(0);
  CHECK_NONZERO(abs(scalar))
  for (int j=0;j<n;j++) {
    F77NAME(zdscal)(n-j,double_one_/abs(scalar),addr(j,j),1);
  }
  return *this;
}

template<> double
SymmetricPositiveMatrix<double,complex<double> >::equilibrate(
Vector<double,double> &s,double &scond) const {
  int n=size(0);
  CHECK_SAME(n,s.size());
  double amax=numeric_limits<double>::infinity();
  int info;
  F77NAME(zpoequb)(n,addr(),n,s.addr(),scond,amax,info);
  CHECK_SAME(info,0);
  return amax;
}

template<> double
SymmetricPositiveMatrix<double,complex<double> >::reciprocalConditionNumber(
) const {
  int n=size(0);
  SymmetricPositiveMatrix<double,complex<double> > *AF=
    OPERATOR_NEW SymmetricPositiveMatrix<double,complex<double> >(*this);
  int info;
  F77NAME(zpotrf)('L',n,AF->addr(),n,info);
  CHECK_SAME(info,0)

  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  double anorm=F77NAME(zlanhe)('O','L',n,addr(),n,rwork);
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double rcond;
  F77NAME(zpocon)('L',n,addr(),n,anorm,rcond,work,rwork,info);
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SymmetricPositiveMatrix<double,complex<double> >*
SymmetricPositiveMatrix<double,complex<double> >::inverse() const {
  int n=size(0);
  SymmetricPositiveMatrix<double,complex<double> > *Ainv=
    OPERATOR_NEW SymmetricPositiveMatrix<double,complex<double> >(*this);
  int info;
  F77NAME(zpotrf)('L',n,Ainv->addr(),n,info);
  CHECK_SAME(info,0);

  F77NAME(zpotri)('L',n,Ainv->addr(),n,info);
  CHECK_SAME(info,0)
  return Ainv;
}
*/

template<> void SymmetricPositiveMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char) const {
  int n=size(0);
  SymmetricPositiveMatrix<double,complex<double> > AF(*this);
  int info;
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  x.copy(b);
  F77NAME(zpotrf)('L',n,AF.addr(),n,info);
  CHECK_SAME(info,0)

  F77NAME(zpotrs)('L',n,1,AF.addr(),n,x.addr(),n,info);
  CHECK_SAME(info,0)
}

template<> void SymmetricPositiveMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char) const {
  int n=size(0);
  SymmetricPositiveMatrix<double,complex<double> > AF(*this);
  int info;
  F77NAME(zpotrf)('L',n,AF.addr(),n,info);
  CHECK_SAME(info,0)

  bool left_side=(side=='L' || side=='l');
  if (left_side) {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1))
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    X.copy(B);
    F77NAME(zpotrs)('L',n,1,AF.addr(),n,X.addr(),n,info);
    CHECK_SAME(info,0)
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0))
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
    for (int i=0;i<k;i++) {
      complex<double> *tt=t;
      const complex<double> *bi=B.addr(i,0);
      for (int j=0;j<n;j++,tt++,bi+=k) *tt=conj(*bi);
      F77NAME(zpotrs)('L',n,1,AF.addr(),n,t,n,info);
      CHECK_SAME(info,0)
      tt=t;
      complex<double> *xi=X.addr(i,0);
      for (int j=0;j<n;j++,tt++,xi+=k) *xi=conj(*tt);
    }
    delete [] t; t=0;
  }
}

template class SymmetricPositiveMatrix<double,complex<double> >;
template void testSymmetricPositiveMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "BandMatrix.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> complex<double>
  TridiagonalMatrix<double,complex<double> >::safety_ =
  complex_double_zero_;
template<> const complex<double>
  TridiagonalMatrix<double,complex<double> >::outofbounds_ =
  complex_double_zero_;
template<> const complex<double>
  TridiagonalMatrix<double,complex<double> >::undefined_ =
  complex<double>(numeric_limits<double>::infinity(),
    numeric_limits<double>::infinity());

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  F77NAME(zcopy)(n,D->addr(),1,M->addr(0,0),n+1);
  F77NAME(zcopy)(n-1,L->addr(),1,M->addr(1,0),n+1);
  F77NAME(zcopy)(n-1,U->addr(),1,M->addr(0,1),n+1);
  return M;
}

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator+(
const SymmetricMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int k=0;k<n;k++) {
    F77NAME(zcopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) F77NAME(zcopy)(n-k-1,M.addr(k+1,k),1,S->addr(k,k+1),n);
  }
  F77NAME(zaxpy)(n,complex_double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(zcopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)=complex_double_one_;
      if (k<n-1) F77NAME(zcopy)(n-k-1,M.addr(k+1,k),1,S->addr(k+1,k),1);
    }
  }
  F77NAME(zaxpy)(n,complex_double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<double,complex<double> > *S=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(zcopy)(k+1,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(zcopy)(k,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)=complex_double_one_;
    }
  }
  F77NAME(zaxpy)(n,complex_double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  S->copy(M);
  F77NAME(zaxpy)(n,complex_double_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator-(
const SymmetricMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  F77NAME(zcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(zcopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    F77NAME(zaxpy)(n-k,double_mone_,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) {
      F77NAME(zaxpy)(n-k-1,double_mone_,M.addr(k+1,k),1,S->addr(k,k+1),n);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const SymmetricMatrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) { 
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  for (int k=0;k<n;k++) {
    F77NAME(zcopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) F77NAME(zcopy)(n-k-1,M.addr(k+1,k),1,S->addr(k,k+1),n);
  }
  F77NAME(zaxpy)(n,complex_double_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  F77NAME(zcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(zcopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(zaxpy)(n-k,double_mone_,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)-=double_one_;
      if (k<n-1) {
        F77NAME(zaxpy)(n-k-1,double_mone_,M.addr(k+1,k),1,
          S->addr(k+1,k),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) { 
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(zcopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)=double_one_;
      if (k<n-1) {
        F77NAME(zcopy)(n-k-1,M.addr(k+1,k),1,S->addr(k+1,k),1);
      }
    }
  }
  F77NAME(zaxpy)(n,complex_double_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<double,complex<double> > *S=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  F77NAME(zcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(zcopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(zaxpy)(k+1,double_mone_,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(zaxpy)(k,double_mone_,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)-=double_one_;
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >* operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<double,complex<double> > *S=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(zcopy)(k+1,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(zcopy)(k,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)=double_one_;
    }
  }
  F77NAME(zaxpy)(n,complex_double_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  F77NAME(zcopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(zcopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  *S-=M;
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  S->copy(M);
  F77NAME(zaxpy)(n,complex_double_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator*(
const SymmetricMatrix<double,complex<double> > &M) const {
// compute by bordering: note that
// [    tau     upsilon e_0^T ] [ sigma s^H ]
// [ e_0 lambda       T       ] [   s    S  ]
//   = [ tau sigma + upsilon e_0^T s , tau s^H + upsilon e_0^T S ]
//   = [ e_0 lambda sigma +      T s , e_0 lambda s^H + T S      ]
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  (*S)(n-2,n-2)=(*D)[n-2]*M(n-2,n-2)+(*U)[n-2]*M(n-1,n-2);
  (*S)(n-1,n-2)=(*L)[n-2]*M(n-2,n-2)+(*D)[n-1]*M(n-1,n-2);
  (*S)(n-2,n-1)=(*D)[n-2]*conj(M(n-1,n-2))+(*U)[n-2]*M(n-1,n-1);
  (*S)(n-1,n-1)=(*L)[n-2]*conj(M(n-1,n-2))+(*D)[n-1]*M(n-1,n-1);
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-3;k>=0;k--) {
    F77NAME(zgtmv)(n-k-1,complex_double_one_,L->addr(k+1),D->addr(k+1),
      U->addr(k+1),M.addr(k+1,k),1,complex_double_zero_,S->addr(k+1,k),1);
      // T s
    complex<double> *ti=t;
    const complex<double> *Mik=M.addr(k+1,k+1);
    for (int i=k+1;i<n;i++,ti++,Mik++) *ti=conj(*Mik);
    F77NAME(zaxpy)(n-k-1,(*U)[k],t,1,S->addr(k,k+1),n); // upsilon e_0^T S
    (*S)(k,k)=(*U)[k]*M(k+1,k); // upsilon e_0^T s
    ti=t;
    Mik=M.addr(k,k);
    for (int i=k;i<n;i++,ti++,Mik++) *ti=conj(*Mik);
    F77NAME(zaxpy)(n-k,(*D)[k],t,1,S->addr(k,k),n); // tau [ sigma , s^H ]
    F77NAME(zaxpy)(n-k,(*L)[k],t,1,S->addr(k+1,k),n);
      // e_0 lambda [ sigma, s^H ]
  }
  delete [] t; t=0;
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SymmetricMatrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) {
// compute by bordering: note that
// [ sigma s^H ] [    tau     upsilon e_0^T ]
// [   s    S  ] [ e_0 lambda       T       ]
//   = [ sigma tau + s^H e_0 lambda , sigma upsilon e_0^T + s^H T ]
//   = [     s tau +   S e_0 lambda ,     s upsilon e_0^T +   S T ]
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  (*S)(n-2,n-2)=M(n-2,n-2)*T(n-2,n-2)+conj(M(n-1,n-2))*T(n-1,n-2);
  (*S)(n-1,n-2)=M(n-1,n-2)*T(n-2,n-2)+     M(n-1,n-1) *T(n-1,n-2);
  (*S)(n-2,n-1)=M(n-2,n-2)*T(n-2,n-1)+conj(M(n-1,n-2))*T(n-1,n-1);
  (*S)(n-1,n-1)=M(n-1,n-2)*T(n-2,n-1)+     M(n-1,n-1) *T(n-1,n-1);
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-3;k>=0;k--) {
    complex<double> *tj=t;
    const complex<double> *Mjk=M.addr(k+1,k);
    for (int j=k+1;j<n;j++,tj++,Mjk++) *tj=conj(*Mjk);
    F77NAME(zgtmv)(n-k-1,double_one_,T.addr(k+1,k+2),
      T.addr(k+1,k+1),T.addr(k+2,k+1),t,1,double_zero_,S->addr(k,k+1),n);
      // s^H T = ( T^T bar(s) )^T
    F77NAME(zaxpy)(n-k-1,T(k+1,k),M.addr(k+1,k+1),1,S->addr(k+1,k),1);
      // S e_0 lambda
    (*S)(k,k)=conj(M(k+1,k))*T(k+1,k); // s^H e_0 lambda
    F77NAME(zaxpy)(n-k,T(k,k),M.addr(k,k),1,S->addr(k,k),1);
      // [ sigma ] tau
      // [   s   ]
    F77NAME(zaxpy)(n-k,T(k,k+1),M.addr(k,k),1,S->addr(k,k+1),1);
      // [ sigma ] upsilon
      // [   s   ]
  }
  delete [] t; t=0;
  return S;
}

template<> Matrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &M) const {
// to compute jth column of product:
//   [      T_11        e_j tau_12 e_0^T ] [  0  ]
//   [ e_0 tau_21 e_j^t      T_22        ] [ ell ]
//     = [   e_j tau_12 e_0^T ell ]
//     = [        T_22 ell        ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *S=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n,double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  if (M_non_unit) {
    F77NAME(zgtmv)(m,complex_double_one_,L->addr(),D->addr(),U->addr(),
      M.addr(),1,complex_double_zero_,S->addr(),1);
    for (int j=1;j<n;j++) {
      (*S)(j-1,j)=(*U)[j-1]*M(j,j); // e_j upsilon e_0^T ell
      if (j<m-1) { // T_22 ell
        F77NAME(zgtmv)(m-j,complex_double_one_,L->addr(j),D->addr(j),
          U->addr(j),M.addr(j,j),1,complex_double_zero_,S->addr(j,j),1);
      } else (*S)(j,j)=(*D)[j]*M(j,j);
    }
  } else {
// note that
// [     tau    , upsilon e_0^T ] [ 1 ] = [ tau + upsilon e_0^T ell ]
// [ e_0 lambda ,        T      ] [ell] = [ e_0 lambda      + T ell ]
    (*S)(0,0)=(*D)[0]+(*U)[0]*M(1,0);
    (*S)(1,0)=(*L)[0];
    F77NAME(zgtmv)(m-1,complex_double_one_,L->addr(1),D->addr(1),
      U->addr(1),M.addr(1,0),1,complex_double_one_,S->addr(1,0),1);
    for (int j=1;j<n;j++) {
      (*S)(j-1,j)=(*U)[j-1]; // e_j upsilon e_0^T ell
      (*S)(j,j)=(*D)[j]; // tau
      if (j<m-1) {
        (*S)(j,j)+=(*U)[j]*M(j+1,j); // upsilon e_0^T ell
        (*S)(j+1,j)=(*L)[j]; // e_0 lambda
      }
      if (j<m-2) { // T ell
        F77NAME(zgtmv)(m-j-1,complex_double_one_,L->addr(j+1),
          D->addr(j+1),U->addr(j+1),M.addr(j+1,j),1,
          complex_double_one_,S->addr(j+1,j),1);
      } else if (j<m-1) (*S)(j+1,j)+=(*D)[j+1]*M(j+1,j);
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >* operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) {
// compute by columns: note that
// L T e_j = L e_{j-1} T_{j-1,j} + L e_j T_{j,j} + L e_{j+1{ T){j+1,j}
//         = L e_{j-1} U[j-1] + L e_j D[j] + L e_{j+1} L[j]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(m-j+1,T(j-1,j),M.addr(j-1,j-1),1,
          S->addr(j-1,j),1);
      }
      F77NAME(zaxpy)(m-j,T(j,j),M.addr(j,j),1,S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(zaxpy)(m-j-1,T(j+1,j),M.addr(j+1,j+1),1,S->addr(j+1,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        (*S)(j-1,j)=T(j-1,j);
        F77NAME(zaxpy)(m-j,T(j-1,j),M.addr(j,j-1),1,S->addr(j,j),1);
      }
      (*S)(j,j)+=T(j,j);
      if (j<m-1) {
        F77NAME(zaxpy)(m-j-1,T(j,j),M.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      if (j<n-1) {
        (*S)(j+1,j)+=T(j+1,j);
        if (j<m-2) {
          F77NAME(zaxpy)(m-j-2,T(j+1,j),M.addr(j+2,j+1),1,
            S->addr(j+2,j),1);
        }
      }
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &M) const {
// compute by columns: note that
// T [ U_1 , U_2 ] = [ T U_1 , T U_2 ]
// and that
// [      T_11       ,e_j upsilon e_0^T ][ u ] = [             T_11 u ]
// [ e_0 lambda e_j^T,      T_22        ][ 0 ] = [ e_0 lambda e_j^T u ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  if (M_non_unit) {
    (*S)(0,0)=(*D)[0]*M(0,0);
    (*S)(1,0)=(*L)[0]*M(0,0);
    for (int j=1;j<m-1;j++) { // T U_1
      F77NAME(zgtmv)(j+1,complex_double_one_,L->addr(),D->addr(),
        U->addr(),M.addr(0,j),1,complex_double_zero_,S->addr(0,j),1);
      (*S)(j+1,j)=(*L)[j]*M(j,j);
    }
    for (int j=m-1;j<n;j++) { // T U_2
      F77NAME(zgtmv)(m,complex_double_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,complex_double_zero_,S->addr(0,j),1);
    }
  } else {
    (*S)(0,0)=(*D)[0];
    (*S)(1,0)=(*L)[0];
    for (int j=1;j<m;j++) {
      F77NAME(zgtmv)(j,complex_double_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,complex_double_zero_,S->addr(0,j),1);
      (*S)(j-1,j)+=(*U)[j-1];
      (*S)(j,j)=(*L)[j-1]*M(j-1,j)+(*D)[j];
      if (j<m-1) (*S)(j+1,j)=(*L)[j];
    }
    for (int j=m;j<n;j++) {
      F77NAME(zgtmv)(m,complex_double_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,complex_double_zero_,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >* operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) {
// compute by columns: note that
// U T e_j = U e_{j-1} T_{j-1,j} + U e_j T_{j,j} + U e_{j+1{ T){j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j,m),T(j-1,j),M.addr(0,j-1),1,S->addr(0,j),1);
      }
      F77NAME(zaxpy)(min(j+1,m),T(j,j),M.addr(0,j),1,S->addr(0,j),1);
      if (j<n-1) {
        F77NAME(zaxpy)(min(j+2,m),T(j+1,j),M.addr(0,j+1),1,
          S->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j-1,m),T(j-1,j),M.addr(0,j-1),1,
          S->addr(0,j),1);
        if (j<=m) (*S)(j-1,j)=T(j-1,j);
      }
      F77NAME(zaxpy)(min(j,m),T(j,j),M.addr(0,j),1,S->addr(0,j),1);
      if (j<m) (*S)(j,j)+=T(j,j);
      if (j<n-1) {
        F77NAME(zaxpy)(min(j+1,m),T(j+1,j),M.addr(0,j+1),1,
          S->addr(0,j),1);
        if (j<m-1) (*S)(j+1,j)+=T(j+1,j);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator*(
const SquareMatrix<double,complex<double> > &M) const {
  int m=M.size(0);
  CHECK_SAME(m,size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(m,complex_double_zero_);
  for (int j=0;j<m;j++) {
    F77NAME(zgtmv)(m,double_one_,L->addr(),D->addr(),U->addr(),
      M.addr(0,j),1,double_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SquareMatrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) {
  int m=M.size(0);
  CHECK_SAME(m,T.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(m,complex_double_zero_);
  for (int j=0;j<m;j++) {
    for (int k=max(0,j-1);k<=min(m-1,j+1);k++) {
      F77NAME(zaxpy)(m,T(k,j),M.addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zgtmv)(m,double_one_,L->addr(),D->addr(),U->addr(),
      M.addr(0,j),1,double_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const TridiagonalMatrix<double,complex<double> > &T) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
      F77NAME(zaxpy)(m,T(k,j),M.addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Vector<double,complex<double> >*
TridiagonalMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<double,complex<double> > *w=OPERATOR_NEW
    Vector<double,complex<double> >(n,complex_double_zero_);
  F77NAME(zgtmv)(n,complex_double_one_,L->addr(),D->addr(),U->addr(),
    v.addr(),1,complex_double_zero_,w->addr(),1);
  return w;
}

template<> double
TridiagonalMatrix<double,complex<double> >::normFrobenius() const {
  return F77NAME(zlangt)('F',dim,L->addr(),D->addr(),U->addr());
}

template<>
double TridiagonalMatrix<double,complex<double> >::normInfinity() const {
  return F77NAME(zlangt)('I',dim,L->addr(),D->addr(),U->addr());
}

template<>
double TridiagonalMatrix<double,complex<double> >::normMaxEntry() const {
  return F77NAME(zlangt)('M',dim,L->addr(),D->addr(),U->addr());
}

template<>
double TridiagonalMatrix<double,complex<double> >::normOne() const {
  return F77NAME(zlangt)('O',dim,L->addr(),D->addr(),U->addr());
}

template<> double
TridiagonalMatrix<double,complex<double> >::reciprocalConditionNumber(
char norm) const {
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *DF=
    OPERATOR_NEW Vector<double,complex<double> >(dim);
  Vector<double,complex<double> > *UF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *UF2=
    OPERATOR_NEW Vector<double,complex<double> >(max(0,dim-2));
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(zgttrf)(dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),ipiv,
    info);
  double rcond=numeric_limits<double>::infinity();
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*dim);
  double anorm=F77NAME(zlangt)(norm,dim,L->addr(),D->addr(),U->addr());
  F77NAME(zgtcon)(norm,dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
    ipiv,anorm,rcond,work,info);
  delete [] work; work=0;
  delete [] ipiv; ipiv=0;
  delete UF2; UF2=0;
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
  return rcond;
}

template<> void TridiagonalMatrix<double,complex<double> >::gtmv(
complex<double> alpha,const Vector<double,complex<double> > &x,
complex<double> beta,Vector<double,complex<double> > &b,char trans)
const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  if (trans=='N' || trans=='n') {
    F77NAME(zgtmv)(n,alpha,L->addr(),D->addr(),U->addr(),
      x.addr(),1,beta,b.addr(),1);
  } else if (trans=='T' || trans=='t') {
    F77NAME(zgtmv)(n,alpha,U->addr(),D->addr(),L->addr(),
      x.addr(),1,beta,b.addr(),1);
  } else {
    if (abs(beta)==double_zero_) b=complex_double_zero_;
    else b*=beta;
    if (abs(alpha)==double_zero_) return;
    if (abs(x[0])>double_zero_) {
      complex<double> temp=x[0]*alpha;
      b[0]+=conj((*D)[0])*temp;
      b[1]+=conj((*U)[0])*temp;
    }
    for (int i=1;i<n-1;i++) {
      if (abs(x[i])>double_zero_) {
        complex<double> temp=x[i]*alpha;
        b[i-1]+=conj((*L)[i-1])*temp;
        b[i]+=conj((*D)[i])*temp;
        b[i+1]+=conj((*U)[i])*temp;
      }
    }
    if (abs(x[n-1])>double_zero_) {
      complex<double> temp=x[n-1]*alpha;
      b[n-2]+=conj((*L)[n-2])*temp;
      b[n-1]+=conj((*D)[n-1])*temp;
    }
  }
}

template<> void TridiagonalMatrix<double,complex<double> >::gtmm(
complex<double> alpha,const Matrix<double,complex<double> > &X,
complex<double> beta,Matrix<double,complex<double> > &B,char side,
char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,size(0));
    if (trans=='N' || trans=='n') { // A X alpha + B beta
      for (int j=0;j<n;j++) {
        F77NAME(zgtmv)(m,alpha,L->addr(),D->addr(),U->addr(),
          X.addr(0,j),1,beta,B.addr(0,j),1);
      }
    } else if (trans=='T' || trans=='t') { // A^T X alpha + B beta
      for (int j=0;j<n;j++) {
        F77NAME(zgtmv)(m,alpha,U->addr(),D->addr(),L->addr(),
          X.addr(0,j),1,beta,B.addr(0,j),1);
      }
    } else { // A^H X alpha + B beta
      if (abs(beta)==double_zero_) B=complex_double_zero_;
      else B*=beta;
      if (abs(alpha)==double_zero_) return;
      for (int j=0;j<n;j++) {
        if (abs(X(0,j))>double_zero_) {
          complex<double> temp=X(0,j)*alpha;
          B(0,j)+=conj((*D)[0])*temp;
          B(1,j)+=conj((*U)[0])*temp;
        }
        for (int i=1;i<m-1;i++) {
          if (abs(X(i,j))>double_zero_) {
            complex<double> temp=X(i,j)*alpha;
            B(i-1,j)+=conj((*L)[i-1])*temp;
            B(i,j)+=conj((*D)[i])*temp;
            B(i+1,j)+=conj((*U)[i])*temp;
          }
        }
        if (abs(X(m-1,j))>double_zero_) {
          complex<double> temp=X(m-1,j)*alpha;
          B(m-2,j)+=conj((*L)[m-2])*temp;
          B(m-1,j)+=conj((*D)[m-1])*temp;
        }
      }
    }
  } else {
    CHECK_SAME(n,size(0));
    if (abs(beta)==double_zero_) B=double_zero_;
    else if (beta!=double_one_) B*=beta;
    if (trans=='N' || trans=='n') { // X A alpha + B beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(zaxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else if (trans=='T' || trans=='t') { // X A^T alpha + B beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(zaxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else { // X A^H alpha + B beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(zaxpy)(m,conj((*this)(j,k))*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

template<> void TridiagonalMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *DF=
    OPERATOR_NEW Vector<double,complex<double> >(dim);
  Vector<double,complex<double> > *UF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  x.copy(b);
  int info;
  if (trans!='N' && trans!='n') {
    F77NAME(zgtsv)(dim,1,UF->addr(),DF->addr(),LF->addr(),x.addr(),dim,
      info);
  } else {
    F77NAME(zgtsv)(dim,1,LF->addr(),DF->addr(),UF->addr(),x.addr(),dim,
      info);
  }
  CHECK_SAME(info,0);
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
}

template<> void TridiagonalMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *DF=
    OPERATOR_NEW Vector<double,complex<double> >(dim);
  Vector<double,complex<double> > *UF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  int info;
  X.copy(B);
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    if (trans!='N' && trans!='n') {
      F77NAME(zgtsv)(m,n,UF->addr(),DF->addr(),LF->addr(),X.addr(),m,
        info);
    } else {
      F77NAME(zgtsv)(m,n,LF->addr(),DF->addr(),UF->addr(),X.addr(),m,
        info);
    }
  } else {
    CHECK_SAME(dim,n);
    Vector<double,complex<double> > *UF2=
      OPERATOR_NEW Vector<double,complex<double> >(dim-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
    F77NAME(zgttrf)(n,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
      ipiv,info);
    CHECK_TEST(info==0);
    complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
    char rtrans=(trans!='N' && trans!='n' ? 'N' : 'T');
    for (int i=0;i<m;i++) {
      F77NAME(zcopy)(n,B.addr(i,0),m,t,1);
      F77NAME(zgttrs)(rtrans,n,1,LF->addr(),DF->addr(),UF->addr(),
        UF2->addr(),ipiv,t,n,info);
      F77NAME(zcopy)(n,t,1,X.addr(i,0),m);
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
template class TridiagonalMatrix<double,complex<double> >;
template void testTridiagonalMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> complex<double>
  SymmetricTridiagonalMatrix<double,complex<double> >::safety_=
  complex_double_zero_;
template<> const complex<double>
  SymmetricTridiagonalMatrix<double,complex<double> >::outofbounds_=
  complex_double_zero_;
template<> const complex<double>
  SymmetricTridiagonalMatrix<double,complex<double> >::undefined_=
  complex<double>(numeric_limits<double>::infinity(),
    numeric_limits<double>::infinity());

template<> SymmetricTridiagonalMatrix<double,complex<double> >::
SymmetricTridiagonalMatrix(int n,complex<double> d) : dim(n) {
  L=OPERATOR_NEW Vector<double,complex<double> >(n-1,d);
  D=OPERATOR_NEW Vector<double,double>(n,d.real());
}

template<> SquareMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  F77NAME(zcopy)(n-1,L->addr(),1,M->addr(1,0),n+1);
  const double *Aii=D->addr();
  complex<double> *Mii=M->addr(0,0);
  for (int i=0;i<dim;i++,Aii++,Mii+=n+1) *Mii=*Aii;
  const complex<double> *Aip1i=L->addr();
  complex<double> *Miip1=M->addr(0,1);
  for (int i=0;i<dim-1;i++,Aip1i++,Miip1+=n+1) *Miip1=conj(*Aip1i);
  return M;
}

template<> void
SymmetricTridiagonalMatrix<double,complex<double> >::fillWith(
complex<double> scalar) {
  *L=scalar; *D=scalar.real();
}

template<> complex<double>
SymmetricTridiagonalMatrix<double,complex<double> >::
upperDiagonalValue(int i) const {
  return conj((*L)[i]);
}

template<> SymmetricTridiagonalMatrix<double,complex<double> >&
SymmetricTridiagonalMatrix<double,complex<double> >::operator=(
complex<double> scalar) {
  *L=scalar; *D=scalar.real(); return *this;
}

template<> SymmetricTridiagonalMatrix<double,complex<double> >&
SymmetricTridiagonalMatrix<double,complex<double> >::operator*=(
complex<double> scalar) {
  int n=this->size(0);
  (*D)*=scalar.real();
  (*L)*=scalar;
  return *this;
}

template<> SymmetricTridiagonalMatrix<double,complex<double> >&
SymmetricTridiagonalMatrix<double,complex<double> >::operator/=(
complex<double> scalar) {
  int n=this->size(0);
  CHECK_NONZERO(scalar.real())
  (*D)/=scalar.real();
  (*L)/=scalar;
  return *this;
}

template<> TridiagonalMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator*(
complex<double> d) const {
  int n=size(0);
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  double *di=D->addr();
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) *Sii=*di;
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),1);
  complex<double> *li=L->addr();
  complex<double> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1=conj(*li);
  *S*=d;
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator/(
complex<double> d) const {
  int n=size(0);
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  double *di=D->addr();
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) *Sii=*di;
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),1);
  complex<double> *li=L->addr();
  complex<double> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1=conj(*li);
  *S/=d;
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator+(
const TridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  S->copy(T);
  double *di=D->addr();
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) *Sii+=*di;
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,S->addr(1,0),1);
  complex<double> *li=L->addr();
  complex<double> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1+=conj(*li);
  return S;
}

template<> SymmetricMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator+(
const SymmetricMatrix<double,complex<double> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricMatrix<double,complex<double> > *T=
    OPERATOR_NEW SymmetricMatrix<double,complex<double> >(n);
  T->copy(S);
  double *di=D->addr();
  complex<double> *Tii=T->addr(0,0);
  for (int i=0;i<n;i++,di++,Tii+=n+1) {
    *Tii+=complex<double>(*di,double_zero_);
  }
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,T->addr(1,0),n+1);
  return T;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0); 
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    complex<double> *Sjj=S->addr(0,0);
    for (int j=0;j<n;j++,Sjj+=n+1) *Sjj=double_one_;
    for (int j=0;j<n-1;j++) {
      F77NAME(zcopy)(n-j-1,M.addr(j+1,j),1,S->addr(j+1,j),1);
    }
  }
  double *di=D->addr();
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii+=n+1) {
    *Sii+=complex<double>(*di,double_zero_);
  }
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,S->addr(1,0),n+1);
  complex<double> *lj=L->addr();
  complex<double> *Sjjp1=S->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Sjjp1+=n+1) *Sjjp1+=conj(*lj);
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
   UpperHessenbergMatrix<double,complex<double> >(n,complex_double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0); 
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(j+1,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=1;j<n;j++) F77NAME(zcopy)(j,U.addr(0,j),1,H->addr(0,j),1);
    complex<double> *Hjj=H->addr(0,0);
    for (int j=0;j<n;j++,Hjj+=n+1) *Hjj=complex_double_one_;
  }
  double *di=D->addr();
  complex<double> *Hii=H->addr(0,0);
  for (int i=0;i<n;i++,di++,Hii+=n+1) {
    *Hii+=complex<double>(*di,double_zero_);
  }
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,H->addr(1,0),n+1);
  complex<double> *lj=L->addr();
  complex<double> *Hjjp1=H->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Hjjp1+=n+1) *Hjjp1+=conj(*lj);
  return H;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  S->copy(M);
  double *di=D->addr();
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii+=*di;
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,S->addr(1,0),n+1);
  complex<double> *li=L->addr();
  complex<double> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1+=n+1) *Siip1+=conj(*li);
  return S;
}

template<> SymmetricTridiagonalMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW SymmetricTridiagonalMatrix<double,complex<double> >(n);
  S->copy(*this);
  F77NAME(daxpy)(n,double_mone_,T.D->addr(),1,S->D->addr(),1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.L->addr(),1,S->L->addr(),1);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator-(
const TridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  double *di=D->addr();
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) *Sii=*di;
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),1);
  complex<double> *li=L->addr();
  complex<double> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1=conj(*li);
  F77NAME(zaxpy)(n,complex_double_mone_,T.addr(0,0),1,S->addr(0,0),1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(1,0),1,S->addr(1,0),1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(0,1),1,S->addr(0,1),1);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >* operator-(
const TridiagonalMatrix<double,complex<double> > &T,
const SymmetricTridiagonalMatrix<double,complex<double> > &St) {
  int n=St.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  S->copy(T);
  double *di=St.diagonalAddr(0);
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) {
    *Sii-=complex<double>(*di,double_zero_);
  }
  F77NAME(zaxpy)(n-1,complex_double_mone_,St.lowerDiagonalAddr(0),1,
    S->addr(1,0),1);
  complex<double> *li=St.lowerDiagonalAddr(0);
  complex<double> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1-=conj(*li);
  return S;
}

template<> SymmetricMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator-(
const SymmetricMatrix<double,complex<double> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricMatrix<double,complex<double> > *T=OPERATOR_NEW
    SymmetricMatrix<double,complex<double> >(n,complex_double_zero_);
  double *di=D->addr();
  complex<double> *Tii=T->addr(0,0);
  for (int i=0;i<n;i++,di++,Tii+=n+1) {
    *Tii=complex<double>(*di,double_zero_);
  }
  F77NAME(zcopy)(n-1,L->addr(),1,T->addr(1,0),n+1);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,complex_double_mone_,S.addr(j,j),1,T->addr(j,j),1);
  }
  return T;
}

template<> SymmetricMatrix<double,complex<double> >* operator-(
const SymmetricMatrix<double,complex<double> > &S,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,complex<double> > *M=
    OPERATOR_NEW SymmetricMatrix<double,complex<double> >(n);
  M->copy(S);
  double *di=T.diagonalAddr(0);
  complex<double> *Mii=M->addr(0,0);
  for (int i=0;i<n;i++,di++,Mii+=n+1) {
    *Mii-=complex<double>(*di,double_zero_);
  }
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.lowerDiagonalAddr(0),1,
    M->addr(1,0),n+1);
  return M;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0); 
  double *di=D->addr();
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di;
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  complex<double> *li=L->addr();
  complex<double> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1+=n+1) *Siip1=conj(*li);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(n-j,complex_double_mone_,M.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    Sii=S->addr(0,0);
    for (int i=0;i<n;i++,Sii+=n+1) *Sii-=double_one_;
    for (int j=0;j<n-1;j++) {
      F77NAME(zaxpy)(n-j-1,complex_double_mone_,M.addr(j+1,j),1,
        S->addr(j+1,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &M,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0); 
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    complex<double> *Sjj=S->addr(0,0);
    for (int j=0;j<n;j++,Sjj+=n+1) *Sjj=complex_double_one_;
    for (int j=0;j<n-1;j++) {
      F77NAME(zcopy)(n-j-1,M.addr(j+1,j),1,
        S->addr(j+1,j),1);
    }
  }
  double *dj=T.diagonalAddr(0);
  complex<double> *Sjj=S->addr(0,0);
  for (int j=0;j<n;j++,dj++,Sjj+=n+1) *Sjj-=*dj;
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),n+1);
  complex<double> *lj=T.lowerDiagonalAddr(0);
  complex<double> *Sjjp1=S->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Sjjp1+=n+1) *Sjjp1-=conj(*lj);
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
   UpperHessenbergMatrix<double,complex<double> >(n,complex_double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0); 
  double *dj=D->addr();
  complex<double> *Hjj=H->addr(0,0);
  for (int j=0;j<n;j++,dj++,Hjj+=n+1) *Hjj=*dj;
  F77NAME(zcopy)(n-1,L->addr(),1,H->addr(1,0),n+1);
  complex<double> *lj=L->addr(0);
  complex<double> *Hjjp1=H->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Hjjp1+=n+1) *Hjjp1=conj(*lj);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(j+1,complex_double_mone_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  } else {
    Hjj=H->addr(0,0);
    for (int j=0;j<n;j++,Hjj+=n+1) *Hjj-=double_one_;
    for (int j=1;j<n;j++) {
      F77NAME(zaxpy)(j,complex_double_mone_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  }
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >* operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
   UpperHessenbergMatrix<double,complex<double> >(n,complex_double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0); 
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(j+1,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    complex<double> *Hjj=H->addr(0,0);
    for (int j=0;j<n;j++,Hjj+=n+1) *Hjj=double_one_;
    for (int j=1;j<n;j++) F77NAME(zcopy)(j,U.addr(0,j),1,H->addr(0,j),1);
  }
  double *dj=T.diagonalAddr(0);
  complex<double> *Hjj=H->addr(0,0);
  for (int j=0;j<n;j++,dj++,Hjj+=n+1) {
    *Hjj-=complex<double>(*dj,double_zero_);
  }
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  complex<double> *lj=T.lowerDiagonalAddr(0);
  complex<double> *Hjjp1=H->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Hjjp1+=n+1) *Hjjp1-=conj(*lj);
  return H;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  double *dj=D->addr();
  complex<double> *Sjj=S->addr(0,0);
  for (int j=0;j<n;j++,dj++,Sjj+=n+1) *Sjj=*dj;
  F77NAME(zcopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  complex<double> *lj=L->addr();
  complex<double> *Sjjp1=S->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Sjjp1+=n+1) *Sjjp1=conj(*lj);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n,complex_double_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  S->copy(M);
  double *di=T.diagonalAddr(0);
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii-=*di;
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),n+1);
  complex<double> *li=T.lowerDiagonalAddr(0);
  complex<double> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1+=n+1) *Siip1-=conj(*li);
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator*(
const SymmetricMatrix<double,complex<double> > &M) const {
// compute by bordering: note that
// [     tau    , bar(lambda) e_0^T ] [ sigma S^H ]
// [ e_0 lambda ,           T       ] [   S    S  ]
//   = [ tau sigma + bar(lambda) e_0^T s , tau s^H + bar(lambda) e_0^T S ]
//   = [ e_0 lambda sigma +          T s , e_0 lambda s^H +          T S ]
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  (*S)(n-2,n-2)=(*D)[n-2]*     M(n-2,n-2) +conj((*L)[n-2])*M(n-1,n-2);
  (*S)(n-1,n-2)=(*L)[n-2]*     M(n-2,n-2) +     (*D)[n-1] *M(n-1,n-2);
  (*S)(n-2,n-1)=(*D)[n-2]*conj(M(n-1,n-2))+conj((*L)[n-2])*M(n-1,n-1);
  (*S)(n-1,n-1)=(*L)[n-2]*conj(M(n-1,n-2))+     (*D)[n-1] *M(n-1,n-1);
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-3;k>=0;k--) {
    F77NAME(zhtmv)(n-k-1,double_one_,L->addr(k+1),D->addr(k+1),
      M.addr(k+1,k),1,double_zero_,S->addr(k+1,k),1); // T s
    complex<double> *ti=t;
    const complex<double> *Mik=M.addr(k+1,k+1);
    for (int i=k+1;i<n;i++,ti++,Mik++) *ti=conj(*Mik);
    F77NAME(zaxpy)(n-k-1,conj((*L)[k]),t,1,S->addr(k,k+1),n);
      // bar(lambda) e_0^T S = bar(lambda) ( S e_0 )^H
    (*S)(k,k)=conj((*L)[k])*M(k+1,k); // bar(lambda) e_0^T s
    ti=t;
    Mik=M.addr(k,k);
    for (int i=k;i<n;i++,ti++,Mik++) *ti=conj(*Mik);
    F77NAME(zaxpy)(n-k,(*D)[k],t,1,S->addr(k,k),n);
      // tau [ sigma , s^H ]
    F77NAME(zaxpy)(n-k,(*L)[k],t,1,S->addr(k+1,k),n);
      // e_0 lamba [ sigma , s^H ]
  }
  delete [] t; t=0;
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SymmetricMatrix<double,complex<double> > &M,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
// compute by bordering: note that
// [ sigma s^H ] [     tau    , bar(lambda) e_0^T ]
// [   s    S  ] [ e_0 lambda ,           T       ]
//   = [ sigma tau + s^H e_0 lambda , sigma bar(lambda) e_0^T + s^H T ]
//   = [   s   tau +  S  e_o lambda ,   s   bar(lambda) e_0^T +  S  T ]
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  (*S)(n-2,n-2)=M(n-2,n-2)*T.diagonalValue(n-2)
               +conj(M(n-1,n-2))*T.lowerDiagonalValue(n-2);
  (*S)(n-1,n-2)=M(n-1,n-2)*T.diagonalValue(n-2)
               +     M(n-1,n-1) *T.lowerDiagonalValue(n-2);
  (*S)(n-2,n-1)=M(n-2,n-2)*T.upperDiagonalValue(n-2)
               +conj(M(n-1,n-2))*T.diagonalValue(n-1);
  (*S)(n-1,n-1)=M(n-1,n-2)*T.upperDiagonalValue(n-2)
               +     M(n-1,n-1) *T.diagonalValue(n-1);
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-3;k>=0;k--) {
    F77NAME(zhtmv)(n-k-1,complex_double_one_,T.lowerDiagonalAddr(k+1),
      T.diagonalAddr(k+1),M.addr(k+1,k),1,complex_double_zero_,t,1);
    complex<double> *tj=t;
    complex<double> *Skj=S->addr(k,k+1);
    for (int j=k+1;j<n;j++,tj++,Skj+=n) *Skj=conj(*tj);
      // s^H T = ( T s )^H
    F77NAME(zaxpy)(n-k-1,T.lowerDiagonalValue(k),M.addr(k+1,k+1),1,
      S->addr(k+1,k),1); // S e_0 lambda
    (*S)(k,k)=conj(M(k+1,k))*T.lowerDiagonalValue(k); // s^H e_0 lambda
    F77NAME(zaxpy)(n-k,T.diagonalValue(k),M.addr(k,k),1,S->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(zaxpy)(n-k,T.upperDiagonalValue(k),M.addr(k,k),1,
      S->addr(k,k+1),1);
      // [sigma ] bar(lambda)
      // [  s   ]
  }
  delete [] t; t=0;
  return S;
}

template<> Matrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &M) const {
// compute by columns: note that
// [     S_11        , e_j bar(sigma) e_0^T ] [ 0 ]
// [ e_0 sigma e_j^T ,         S_22         ] [ell]
//   = [ e_j bar(sigma) e_0^T ell ]
//   = [ S_22 ell ]
  int m=M.size(0),n=M.size(1);;
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) (*S)(j-1,j)=upperDiagonalValue(j-1)*M(j,j);
        // sigma e_0^T ell
      if (j+1<m) { // S_22 ell
        F77NAME(zhtmv)(m-j,complex_double_one_,L->addr(j),D->addr(j),
          M.addr(j,j),1,complex_double_zero_,S->addr(j,j),1);
      } else (*S)(j,j)=diagonalValue(j)*M(j,j);
    }
  } else {
//  note that
// [      sigma , bar(lambda) e_0^T ] [ 1 ]
// [ e_0 lambda ,                 S ] [ell]
//   = [ sigma + bar(lambda) e_0^T ell ]
//   = [ e_0 lambda + S ell ]
    for (int j=0;j<n;j++) {
      if (j>0) (*S)(j-1,j)=upperDiagonalValue(j-1);
      (*S)(j,j)=diagonalValue(j);
      if (j+1<m) {
        (*S)(j,j)+=upperDiagonalValue(j)*M(j+1,j);
        (*S)(j+1,j)=lowerDiagonalValue(j);
        if (j+2<m) {
          F77NAME(zhtmv)(m-j-1,complex_double_one_,L->addr(j+1),
            D->addr(j+1),M.addr(j+1,j),1,complex_double_zero_,
            S->addr(j+1,j),1);
        } else (*S)(j+1,j)+=diagonalValue(j+1)*M(j+1,j);
      }
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >* operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &M,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
// compute by columns: note that
// L T e_j = L e_{j-1} T_{j-1,j} + L e_j T_{j,j} + L e_{j+1} T_{j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(m-j+1,T.upperDiagonalValue(j-1),M.addr(j-1,j-1),1,
          S->addr(j-1,j),1);
      }
      F77NAME(zaxpy)(m-j,T.diagonalValue(j),M.addr(j,j),1,S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(zaxpy)(m-j-1,T.lowerDiagonalValue(j),M.addr(j+1,j+1),1,
          S->addr(j+1,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        (*S)(j-1,j)=T.upperDiagonalValue(j-1);
        F77NAME(zaxpy)(m-j,T.upperDiagonalValue(j-1),M.addr(j,j-1),1,
          S->addr(j,j),1);
      }
      (*S)(j,j)+=T.diagonalValue(j);
      if (j<m-1) {
        F77NAME(zaxpy)(m-j-1,T.diagonalValue(j),M.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      if (j<n-1) {
        (*S)(j+1,j)+=T.lowerDiagonalValue(j);
        if (j<m-2) {
          F77NAME(zaxpy)(m-j-2,T.lowerDiagonalValue(j),M.addr(j+2,j+1),1,
            S->addr(j+2,j),1);
        }
      }
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &M) const {
// compute by columns: note that
// [       S_11      , e_k bar(sigma) e_0^T ] [ u ]= [ S_11 u ]
// [ e_0 sigma e_k^T ,         S_22         ] [   ]= [ e_0 sigma e_k^T u ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<m-1;j++) {
      F77NAME(zhtmv)(j+1,complex_double_one_,L->addr(),D->addr(),
        M.addr(0,j),1,complex_double_zero_,S->addr(0,j),1);
      (*S)(j+1,j)=lowerDiagonalValue(j)*M(j,j);
    }
    for (int j=m-1;j<n;j++) {
      F77NAME(zhtmv)(m,complex_double_one_,L->addr(),D->addr(),
        M.addr(0,j),1,complex_double_zero_,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(zhtmv)(j,complex_double_one_,L->addr(),D->addr(),
          M.addr(0,j),1,complex_double_zero_,S->addr(0,j),1);
        (*S)(j-1,j)+=upperDiagonalValue(j-1);
        (*S)(j,j)=lowerDiagonalValue(j-1)*M(j-1,j);
      }
      (*S)(j,j)+=diagonalValue(j);
      if (j+1<m) (*S)(j+1,j)=lowerDiagonalValue(j);
    }
    for (int j=m;j<n;j++) {
      F77NAME(zhtmv)(m,complex_double_one_,L->addr(),D->addr(),
        M.addr(0,j),1,complex_double_zero_,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >* operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &M,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
// compute by columns: note that
// U T e_j = U e_{j-1} T_{j-1,j} + U e_j T_{j,j} + U e_{j+1} T_{j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j,m),T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
          S->addr(0,j),1);
      }
      F77NAME(zaxpy)(min(j+1,m),T.diagonalValue(j),M.addr(0,j),1,
        S->addr(0,j),1);
      if (j<n-1) {
        F77NAME(zaxpy)(min(j+2,m),T.lowerDiagonalValue(j),M.addr(0,j+1),1,
          S->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(min(j-1,m),T.upperDiagonalValue(j-1),
          M.addr(0,j-1),1,S->addr(0,j),1);
        if (j<=m) (*S)(j-1,j)=T.upperDiagonalValue(j-1);
      }
      F77NAME(zaxpy)(min(j,m),T.diagonalValue(j),M.addr(0,j),1,
        S->addr(0,j),1);
      if (j<m) (*S)(j,j)+=T.diagonalValue(j);
      if (j<n-1) {
        F77NAME(zaxpy)(min(j+1,m),T.lowerDiagonalValue(j),M.addr(0,j+1),1,
          S->addr(0,j),1);
        if (j<m-1) (*S)(j+1,j)+=T.lowerDiagonalValue(j);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator*(
const SquareMatrix<double,complex<double> > &M) const {
  int n=M.size(0);
  CHECK_SAME(n,size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zhtmv)(n,double_one_,L->addr(),D->addr(),
      M.addr(0,j),1,double_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SquareMatrix<double,complex<double> > &M,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(zaxpy)(n,T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(zaxpy)(n,T.diagonalValue(j),M.addr(0,j),1,S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(zaxpy)(n,T.lowerDiagonalValue(j),M.addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zhtmv)(m,double_one_,L->addr(),D->addr(),
      M.addr(0,j),1,double_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<double,complex<double> > *S=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(zaxpy)(m,T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(zaxpy)(m,T.diagonalValue(j),M.addr(0,j),1,S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(zaxpy)(m,T.lowerDiagonalValue(j),M.addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> Vector<double,complex<double> >*
SymmetricTridiagonalMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<double,complex<double> > *w=OPERATOR_NEW
    Vector<double,complex<double> >(n,complex_double_zero_);
  F77NAME(zhtmv)(n,complex_double_one_,L->addr(),D->addr(),
    v.addr(),1,double_zero_,w->addr(),1);
  return w;
}

template<> double
SymmetricTridiagonalMatrix<double,complex<double> >::normFrobenius()
const {
  return F77NAME(zlanht)('F',dim,D->addr(),L->addr());
}

template<> double
SymmetricTridiagonalMatrix<double,complex<double> >::normInfinity()
const {
  return F77NAME(zlanht)('I',dim,D->addr(),L->addr());
}

template<> double
SymmetricTridiagonalMatrix<double,complex<double> >::normMaxEntry()
const {
  return F77NAME(zlanht)('M',dim,D->addr(),L->addr());
}

template<> double
SymmetricTridiagonalMatrix<double,complex<double> >::normOne() const {
  return F77NAME(zlanht)('O',dim,D->addr(),L->addr());
}

template<> double SymmetricTridiagonalMatrix<double,complex<double> >
::reciprocalConditionNumber(char norm) const {
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *UF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *UF2=
    OPERATOR_NEW Vector<double,complex<double> >(max(0,dim-2));

  complex<double> *zd=OPERATOR_NEW_BRACKET(complex<double>,dim);
  const double *di=diagonalAddr(0);
  complex<double> *zdi=zd;
  for (int i=0;i<dim;i++,di++,zdi++) {
    *zdi=complex<double>(*di,double_zero_);
  }
  LF->copy(*L);
  UF->copy(*L);
  double anorm=F77NAME(zlangt)(norm,dim,L->addr(),zd,L->addr());

  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(zgttrf)(dim,LF->addr(),zd,UF->addr(),UF2->addr(),ipiv,info);
  double rcond=numeric_limits<double>::infinity();
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*dim);
  F77NAME(zgtcon)(norm,dim,LF->addr(),zd,UF->addr(),UF2->addr(),
    ipiv,anorm,rcond,work,info);
  delete [] work; work=0;
  delete [] ipiv; ipiv=0;
  delete UF2; UF2=0;
  delete UF; UF=0;
  delete [] zd; zd=0;
  delete LF; LF=0;
  return rcond;
}

template<> void SymmetricTridiagonalMatrix<double,complex<double> >::stmv(
complex<double> alpha,const Vector<double,complex<double> > &x,
complex<double> beta,Vector<double,complex<double> > &b) const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  F77NAME(zhtmv)(n,alpha,L->addr(),D->addr(),x.addr(),1,beta,b.addr(),1);
}

template<> void SymmetricTridiagonalMatrix<double,complex<double> >::stmm(
complex<double> alpha,const Matrix<double,complex<double> > &X,
complex<double> beta,Matrix<double,complex<double> > &B,char side) const {
  if (side=='L' || side=='l') {
    int m=size(0),n=X.size(1);
    CHECK_SAME(n,B.size(1));
    CHECK_SAME(m,X.size(0));
    CHECK_SAME(m,B.size(0));
    for (int j=0;j<n;j++) {
      F77NAME(zhtmv)(m,alpha,L->addr(),D->addr(),X.addr(0,j),1,beta,
        B.addr(0,j),1);
    }
  } else {
    int m=X.size(0),n=X.size(1);
    CHECK_SAME(n,size(0));
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,B.size(1));
    if (abs(beta)==double_zero_) B=complex_double_zero_;
    else if (beta!=double_one_) B*=beta;
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(zaxpy)(m,upperDiagonalValue(j-1)*alpha,X.addr(0,j-1),1,
          B.addr(0,j),1);
      }
      F77NAME(zaxpy)(m,diagonalValue(j)*alpha,X.addr(0,j),1,
        B.addr(0,j),1);
      if (j<n-1) {
        F77NAME(zaxpy)(m,lowerDiagonalValue(j)*alpha,X.addr(0,j+1),1,
          B.addr(0,j),1);
      }
    }
  }
}

template<> void
SymmetricTridiagonalMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,Vector<double,
complex<double> > &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *UF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *DF=
    OPERATOR_NEW Vector<double,complex<double> >(dim);
  double *di=D->addr();
  complex<double> *dfi=DF->addr();
  for (int i=0;i<dim;i++,di++,dfi++) {
    *dfi=complex<double>(*di,double_zero_);
  }
  const complex<double> *li=L->addr();
  complex<double> *ufi=UF->addr();
  for (int i=0;i<dim-1;i++,li++,ufi++) *ufi=conj(*li);
  LF->copy(*L);
  x.copy(b);
  int info;
  F77NAME(zgtsv)(dim,1,LF->addr(),DF->addr(),UF->addr(),x.addr(),dim,
    info);
  CHECK_SAME(info,0);
  delete DF; DF=0;
  delete LF; LF=0;
  delete UF; UF=0;
}

template<> void
SymmetricTridiagonalMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *UF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,complex<double> > *DF=
    OPERATOR_NEW Vector<double,complex<double> >(dim);
  double *di=D->addr();
  complex<double> *dfi=DF->addr();
  for (int i=0;i<dim;i++,di++,dfi++) {
    *dfi=complex<double>(*di,double_zero_);
  }
  const complex<double> *li=L->addr();
  complex<double> *ufi=UF->addr();
  for (int i=0;i<dim-1;i++,li++,ufi++) *ufi=conj(*li);
  LF->copy(*L);
  int info;
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    X.copy(B);
    F77NAME(zgtsv)(m,n,LF->addr(),DF->addr(),UF->addr(),X.addr(),m,info);
    CHECK_SAME(info,0);
  } else {
    CHECK_SAME(dim,n);
    complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
    Vector<double,complex<double> > *UF2=
      OPERATOR_NEW Vector<double,complex<double> >(dim-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
    F77NAME(zgttrf)(n,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),ipiv,
      info);
    CHECK_TEST(info==0);
    for (int i=0;i<m;i++) {
      F77NAME(zcopy)(n,B.addr(i,0),m,t,1);
      F77NAME(zgttrs)('T',n,1,LF->addr(),DF->addr(),UF->addr(),
        UF2->addr(),ipiv,t,n,info);
      CHECK_TEST(info==0);
      F77NAME(zcopy)(n,t,1,X.addr(i,0),m);
    }
    delete [] ipiv; ipiv=0;
    delete UF2; UF2=0;
    delete [] t; t=0;
  }
  delete DF; DF=0;
  delete LF; LF=0;
  delete UF; UF=0;
}

template class SymmetricTridiagonalMatrix<double,complex<double> >;
template void testSymmetricTridiagonalMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveMatrix<double,complex<double> >*
SymmetricPositiveTridiagonalMatrix<double,complex<double> >::operator+(
const SymmetricPositiveMatrix<double,complex<double> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricPositiveMatrix<double,complex<double> > *T=
    OPERATOR_NEW SymmetricPositiveMatrix<double,complex<double> >(n);
  T->copy(S);
  double *di=D->addr();
  complex<double> *Tii=T->addr(0,0);
  for (int i=0;i<n;i++,di++,Tii+=n+1) {
    *Tii+=complex<double>(*di,double_zero_);
  }
  F77NAME(zaxpy)(n-1,complex_double_one_,L->addr(),1,T->addr(1,0),n+1);
  return T;
}

template<> double
SymmetricPositiveTridiagonalMatrix<double,complex<double> >
::reciprocalConditionNumber(char norm) const {
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  LF->copy(*L);
  DF->copy(*D);
  int info;
  F77NAME(zpttrf)(dim,DF->addr(),LF->addr(),info);
  CHECK_TEST(info==0);
  double rcond=numeric_limits<double>::infinity();
  double *work=OPERATOR_NEW_BRACKET(double,2*dim);
  double anorm=normOne();
  F77NAME(zptcon)(dim,DF->addr(),LF->addr(),anorm,rcond,work,info);
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

template<> void
SymmetricPositiveTridiagonalMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  LF->copy(*L);
  DF->copy(*D);
  x.copy(b);
  int info;
  F77NAME(zptsv)(dim,1,DF->addr(),LF->addr(),x.addr(),dim,info);
  CHECK_SAME(info,0);
  delete DF; DF=0;
  delete LF; LF=0;
}

template<> void
SymmetricPositiveTridiagonalMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<double,complex<double> > *LF=
    OPERATOR_NEW Vector<double,complex<double> >(dim-1);
  Vector<double,double> *DF=OPERATOR_NEW Vector<double,double>(dim);
  LF->copy(*L);
  DF->copy(*D);
  int info;
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    X.copy(B);
    F77NAME(zptsv)(m,n,DF->addr(),LF->addr(),X.addr(),m,info);
    CHECK_SAME(info,0);
  } else {
    CHECK_SAME(dim,n);
    complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
    F77NAME(zpttrf)(n,DF->addr(),LF->addr(),info);
    CHECK_SAME(info,0);
    for (int i=0;i<m;i++) {
      const complex<double> *Bij=B.addr(i,0);
      complex<double> *tj=t;
      for (int j=0;j<n;j++,Bij++,tj++) *tj=conj(*Bij);
      F77NAME(zpttrs)('L',n,1,DF->addr(),LF->addr(),t,n,info);
      CHECK_SAME(info,0);
      complex<double> *Xij=X.addr(i,0);
      tj=t;
      for (int j=0;j<n;j++,Xij++,tj++) *Xij=conj(*tj);
    }
    delete [] t; t=0;
  }
  delete DF; DF=0;
  delete LF; LF=0;
}

template class
  SymmetricPositiveTridiagonalMatrix<double,complex<double> >;
template void
  testSymmetricPositiveTridiagonalMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SquareMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  F77NAME(zcopy)(n,D->addr(),1,M->addr(0,0),n+1);
  return M;
}

template<> SymmetricTridiagonalMatrix<double,complex<double> >*
operator+(const DiagonalMatrix<double,double> &A,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<double,complex<double> >(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_one_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >*
operator+(const DiagonalMatrix<double,complex<double> > &A,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  F77NAME(zcopy)(n-1,T.lowerDiagonalAddr(0),1,S->lowerDiagonalAddr(),1);
  const complex<double> *di=A.addr();
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Sii=S->diagonalAddr();
  for (int i=0;i<n;i++,di++,Tii++,Sii++) *Sii=*di+*Tii;
  const complex<double> *Tip1i=T.lowerDiagonalAddr(0);
  complex<double> *Siip1=S->upperDiagonalAddr();
  for (int i=0;i<n-1;i++,Tip1i++,Siip1++) *Siip1=conj(*Tip1i);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator+(
const TridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=OPERATOR_NEW
    TridiagonalMatrix<double,complex<double> >(n);
  S->copy(T);
  F77NAME(zaxpy)(n,complex_double_one_,addr(),1,S->addr(0,0),1);
  return S;
}

template<> SymmetricMatrix<double,complex<double> >* operator+(
const DiagonalMatrix<double,double> &A,
const SymmetricMatrix<double,complex<double> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,complex<double> > *S=
    OPERATOR_NEW SymmetricMatrix<double,complex<double> >(n);
  S->copy(T);
  const double *di=A.addr();
  complex<double> *Sii=S->addr();
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii+=*di;
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator+(
const DiagonalMatrix<double,complex<double> > &A,
const SymmetricMatrix<double,complex<double> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,T.addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=A[j];
    if (j>0) {
      const complex<double> *Tji=T.addr(j,0);
      complex<double> *Sij=S->addr(0,j);
      for (int i=0;i<j;i++,Tji+=n,Sij++) *Sij=conj(*Tji);
    }
  }
  return S;
}

template<> UpperTriangularMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<double,complex<double> > *S=OPERATOR_NEW
    UpperTriangularMatrix<double,complex<double> >(n);
  S->copy(U);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    F77NAME(zaxpy)(n,complex_double_one_,addr(),1,S->addr(),n+1);
  } else {
    const complex<double> *di=D->addr();
    complex<double> *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di+double_one_;
  }
  return S;
}

template<> LowerTriangularMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<double,complex<double> > *S=OPERATOR_NEW
    LowerTriangularMatrix<double,complex<double> >(n);
  S->copy(L);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    F77NAME(zaxpy)(n,double_one_,addr(),1,S->addr(),n+1);
  } else {
    const complex<double> *di=D->addr();
    complex<double> *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di+double_one_;
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n);
  S->copy(M);
  F77NAME(zaxpy)(n,complex_double_one_,addr(),1,S->addr(),n+1);
  return S;
}

template<> SymmetricTridiagonalMatrix<double,complex<double> >* operator-(
const DiagonalMatrix<double,double> &A,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<double,complex<double> >(n);
  S->copy(T);
  *S*=double_mone_;
  F77NAME(daxpy)(n,double_one_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >* operator-(
const DiagonalMatrix<double,complex<double> > &A,
const SymmetricTridiagonalMatrix<double,complex<double> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  const complex<double> *di=A.addr();
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Sii=S->diagonalAddr();
  for (int i=0;i<n;i++,di++,Tii++,Sii++) *Sii=*di-*Tii;
  const complex<double> *Tip1i=T.lowerDiagonalAddr(0);
  complex<double> *Siip1=S->upperDiagonalAddr();
  complex<double> *Sip1i=S->lowerDiagonalAddr();
  for (int i=0;i<n-1;i++,Tip1i++,Siip1++,Sip1i++) {
    *Sip1i=-*Tip1i;
    *Siip1=-conj(*Tip1i);
  }
  return S;
}

template<> SymmetricTridiagonalMatrix<double,complex<double> >* operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<double,complex<double> >(n);
  S->copy(T);
  F77NAME(daxpy)(n,double_mone_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >* operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const DiagonalMatrix<double,complex<double> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  F77NAME(zcopy)(n-1,T.lowerDiagonalAddr(0),1,S->lowerDiagonalAddr(),1);
  const complex<double> *di=A.addr();
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Sii=S->diagonalAddr();
  for (int i=0;i<n;i++,di++,Tii++,Sii++) *Sii=*Tii-*di;
  const complex<double> *Tip1i=T.lowerDiagonalAddr(0);
  complex<double> *Siip1=S->upperDiagonalAddr();
  for (int i=0;i<n-1;i++,Tip1i++,Siip1++) *Siip1=conj(*Tip1i);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator-(
const TridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=OPERATOR_NEW
    TridiagonalMatrix<double,complex<double> >(n);
  S->copy(T);
  *S*=double_mone_;
  F77NAME(zaxpy)(n,complex_double_one_,addr(),1,S->addr(0,0),1);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >* operator-(
const TridiagonalMatrix<double,complex<double> > &T,
const DiagonalMatrix<double,complex<double> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=OPERATOR_NEW
    TridiagonalMatrix<double,complex<double> >(n);
  S->copy(T);
  F77NAME(zaxpy)(n,complex_double_mone_,A.addr(),1,S->addr(0,0),1);
  return S;
}

template<> SymmetricMatrix<double,complex<double> >* operator-(
const DiagonalMatrix<double,double> &A,
const SymmetricMatrix<double,complex<double> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricMatrix<double,complex<double> >(n);
  S->copy(T);
  *S*=complex_double_mone_;
  const double *di=A.addr();
  complex<double> *Sii=S->addr();
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii+=*di;
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const DiagonalMatrix<double,complex<double> > &A,
const SymmetricMatrix<double,complex<double> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    (*S)(j,j)=A[j]-T(j,j);
    if (j+1<n) {
      const complex<double> *Tij=T.addr(j+1,j);
      complex<double> *Sij=S->addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Tij++,Sij++,Sji+=n) {
        *Sij=-*Tij;
        *Sji=-conj(*Tij);
      }
    }
  }
  return S;
}

template<> SymmetricMatrix<double,complex<double> >* operator-(
const SymmetricMatrix<double,complex<double> > &T,
const DiagonalMatrix<double,double> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricMatrix<double,complex<double> >(n);
  S->copy(T);
  const double *di=A.addr();
  complex<double> *Sii=S->addr();
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii-=*di;
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const SymmetricMatrix<double,complex<double> > &T,
const DiagonalMatrix<double,complex<double> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    (*S)(j,j)=T(j,j)-A[j];
    if (j+1<n) {
      const complex<double> *Tij=T.addr(j+1,j);
      complex<double> *Sij=S->addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Tij++,Sij++,Sji+=n) {
        *Sij=*Tij;
        *Sji=conj(*Tij);
      }
    }
  }
  return S;
}

template<> UpperTriangularMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<double,complex<double> > *S=OPERATOR_NEW
    UpperTriangularMatrix<double,complex<double> >(n);
  S->copy(U);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U) ==0) {
    *S*=double_mone_;
    F77NAME(zaxpy)(n,complex_double_one_,addr(),1,S->addr(),n+1);
  } else {
    const complex<double> *dj=D->addr();
    for (int j=0;j<n;j++,dj++) {
      complex<double> *Sij=S->addr(0,j);
      const complex<double> *Uij=U.addr(0,j);
      for (int i=0;i<j;i++,Sij++,Uij++) *Sij=-*Uij;
      *Sij=*dj-double_one_;
    }
  }
  return S;
}

template<> UpperTriangularMatrix<double,complex<double> >* operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const DiagonalMatrix<double,complex<double> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<double,complex<double> > *S=OPERATOR_NEW
    UpperTriangularMatrix<double,complex<double> >(n);
  S->copy(U);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    F77NAME(zaxpy)(n,complex_double_mone_,A.addr(),1,S->addr(),n+1);
  } else {
    const complex<double> *di=A.addr();
    complex<double> *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=double_one_-*di;
  }
  return S;
}

template<> LowerTriangularMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<double,complex<double> > *S=OPERATOR_NEW
    LowerTriangularMatrix<double,complex<double> >(n);
  S->copy(L);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    *S*=complex_double_mone_;
    F77NAME(zaxpy)(n,complex_double_one_,addr(),1,S->addr(),n+1);
  } else {
    const complex<double> *dj=D->addr();
    for (int j=0;j<n;j++,dj++) {
      complex<double> *Sij=S->addr(j,j);
      *Sij=*dj-double_one_;
      if (j+1<n) {
        Sij++;
        const complex<double> *Lij=L.addr(j+1,j);
        for (int i=j+1;i<n;i++,Sij++,Lij++) *Sij=-*Lij;
      }
    }
  }
  return S;
}

template<> LowerTriangularMatrix<double,complex<double> >* operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const DiagonalMatrix<double,complex<double> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<double,complex<double> > *S=
    OPERATOR_NEW LowerTriangularMatrix<double,complex<double> >(n);
  S->copy(L);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    F77NAME(zaxpy)(n,complex_double_mone_,A.addr(),1,S->addr(),n+1);
  } else {
    const complex<double> *di=A.addr();
    complex<double> *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=double_one_-*di;
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  S->copy(M);
  *S*=complex_double_mone_;
  F77NAME(zaxpy)(n,complex_double_one_,addr(),1,S->addr(),n+1);
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const DiagonalMatrix<double,complex<double> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  S->copy(M);
  F77NAME(zaxpy)(n,complex_double_mone_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >* 
DiagonalMatrix<double,complex<double> >::operator*(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  for (int i=0;i<n;i++) {
    complex<double> di=(*D)[i];
    if (i>0) (*S)(i,i-1)=di*T.lowerDiagonalValue(i-1);
    (*S)(i,i)=di*T.diagonalValue(i);
    if (i<n-1) (*S)(i,i+1)=di*conj(T.lowerDiagonalValue(i));
  }
  return S;
}

template<> TridiagonalMatrix<double,complex<double> >* operator*(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const DiagonalMatrix<double,complex<double> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<double,complex<double> > *S=
    OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    complex<double> dj=A[j];
    if (j>0) (*S)(j-1,j)=conj(T.lowerDiagonalValue(j-1))*dj;
    (*S)(j,j)=T.diagonalValue(j)*dj;
    if (j<n-1) (*S)(j+1,j)=T.lowerDiagonalValue(j)*dj;
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
DiagonalMatrix<double,complex<double> >::operator*(
const SymmetricMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(dim);
  for (int i=0;i<n;i++) {
    complex<double> di=D->operator[](i);
    if (i>0) {
      const complex<double> *Tij=T.addr(i,0);
      complex<double> *Sij=S->addr(i,0);
      for (int j=0;j<i;j++,Tij+=n,Sij+=n) *Sij=di*(*Tij);
    }
    const complex<double> *Tji=T.addr(i,i);
    complex<double> *Sij=S->addr(i,i);
    for (int j=i;j<n;j++,Tji++,Sij+=n) *Sij=di*conj(*Tji);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SymmetricMatrix<double,complex<double> > &T,
const DiagonalMatrix<double,complex<double> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    complex<double> dj=A[j];
    if (j>0) {
      const complex<double> *Tji=T.addr(j,0);
      complex<double> *Sij=S->addr(0,j);
      for (int i=0;i<j;i++,Tji+=n,Sij++) *Sij=conj(*Tji)*dj;
    }
    const complex<double> *Tij=T.addr(j,j);
    complex<double> *Sij=S->addr(j,j);
    for (int i=j;i<n;i++,Tij++,Sij++) *Sij=(*Tij)*dj;
  }
  return S;
}

template class DiagonalMatrix<double,complex<double> >;
template void testDiagonalMatrix(double,complex<double> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const complex<double>
  UpperHessenbergMatrix<double,complex<double> >::outofbounds_ =
  complex_double_zero_;
template<> const complex<double> UpperHessenbergMatrix<double,
  complex<double> >::undefined_(numeric_limits<double>::infinity(),
    numeric_limits<double>::infinity());
template<> complex<double>
  UpperHessenbergMatrix<double,complex<double> >::safety_ =
  complex_double_zero_;

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(n,j+2),addr(0,j),1,M->addr(0,j),1);
  }
  return M;
}

template<> UpperHessenbergMatrix<double,complex<double> >&
UpperHessenbergMatrix<double,complex<double> >::operator+=(
const UpperHessenbergMatrix<double,complex<double> > &H) {
  int n=size(0);
  CHECK_SAME(n,H.size(0));
  complex<double> *colj=addr();
  const complex<double> *Hcolj=H.addr();
  for (int j=0;j<n;j++,colj+=n,Hcolj+=n) {
    F77NAME(zaxpy)(min(n,j+2),complex_double_one_,Hcolj,1,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<double,complex<double> >&
UpperHessenbergMatrix<double,complex<double> >::operator-=(
const UpperHessenbergMatrix<double,complex<double> > &H) {
  int n=size(0);
  CHECK_SAME(n,H.size(0));
  complex<double> *colj=addr();
  const complex<double> *Hcolj=H.addr();
  for (int j=0;j<n;j++,colj+=n,Hcolj+=n) {
    F77NAME(zaxpy)(min(n,j+2),complex_double_mone_,Hcolj,1,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<double,complex<double> >&
UpperHessenbergMatrix<double,complex<double> >::operator*=(
complex<double> scalar) {
  int n=size(0);
  complex<double> *colj=addr();
  for (int j=0;j<n;j++,colj+=n) {
    F77NAME(zscal)(min(n,j+2),scalar,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<double,complex<double> >&
UpperHessenbergMatrix<double,complex<double> >::operator/=(
complex<double> scalar) {
  int n=size(0);
  complex<double> *colj=addr();
  complex<double> s=complex_double_one_/scalar;
  for (int j=0;j<n;j++,colj+=n) {
    F77NAME(zscal)(min(n,j+2),s,colj,1);
  }
  return *this;
}

template<> void UpperHessenbergMatrix<double,complex<double> >::copy(
const Matrix<double,complex<double> > &M) {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(j+2,n),M.addr(0,j),1,addr(0,j),1);
  }
}

/*
template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::transpose() const {
  int n=size(0);
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(j,0),n);
  }
  return S;
}
*/

template<> void UpperHessenbergMatrix<double,complex<double> >::copyFrom(
int m,const SquareMatrix<double,complex<double> > &S) {
  m=min(m,size(0));
  m=min(m,S.size(0));
  for (int j=0;j<m;j++) {
    F77NAME(zcopy)(min(j+2,m),S.addr(0,j),1,addr(0,j),1);
  }
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator+(
const DiagonalMatrix<double,complex<double> > &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
      complex_double_zero_);
  H->copy(*this);
  F77NAME(zaxpy)(n,complex_double_one_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator+(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(*this);
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Hii=H->addr(0,0);
  for (int i=0;i<n;i++,Tii++,Hii+=n+1) *Hii+=*Tii;
  F77NAME(zaxpy)(n-1,complex_double_one_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  const complex<double> *Tip1i=T.lowerDiagonalAddr(0);
  complex<double> *Hiip1=H->addr(0,1);
  for (int i=0;i<n-1;i++,Tip1i++,Hiip1+=n+1) *Hiip1+=conj(*Tip1i);
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator+(
const TridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(*this);
  F77NAME(zaxpy)(n,complex_double_one_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator+(
const SymmetricMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      const complex<double> *Mij=M.addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Mij++,Sji+=n) *Sji=conj(*Mij);
    }
    F77NAME(zaxpy)(min(j+2,n),double_one_,addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      F77NAME(zaxpy)(n-j,complex_double_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)+=complex_double_one_;
      if (j<n-1) {
        F77NAME(zaxpy)(n-j-1,complex_double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(*this);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(j+1,complex_double_one_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)+=complex_double_one_;
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_one_,U.addr(0,j),1,
          H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(j+2,n),complex_double_one_,addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator-(
const DiagonalMatrix<double,complex<double> > &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(*this);
  F77NAME(zaxpy)(n,complex_double_mone_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >* operator-(
const DiagonalMatrix<double,complex<double> > &D,
const UpperHessenbergMatrix<double,complex<double> > &U) {
  int n=U.size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(U);
  *H*=complex_double_mone_;
  F77NAME(zaxpy)(n,complex_double_one_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(*this);
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Hii=H->addr(0,0);
  for (int i=0;i<n;i++,Tii++,Hii+=n+1) *Hii-=*Tii;
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  const complex<double> *Tip1i=T.lowerDiagonalAddr(0);
  complex<double> *Hiip1=H->addr(0,1);
  for (int i=0;i<n-1;i++,Tip1i++,Hiip1+=n+1) *Hiip1-=conj(*Tip1i);
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >* operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const UpperHessenbergMatrix<double,complex<double> > &U) {
  int n=U.size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(U);
  *H*=complex_double_mone_;
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Hii=H->addr(0,0);
  for (int i=0;i<n;i++,Tii++,Hii+=n+1) *Hii+=*Tii;
  F77NAME(zaxpy)(n-1,complex_double_one_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  const complex<double> *Tip1i=T.lowerDiagonalAddr(0);
  complex<double> *Hiip1=H->addr(0,1);
  for (int i=0;i<n-1;i++,Tip1i++,Hiip1+=n+1) *Hiip1+=conj(*Tip1i);
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator-(
const TridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(*this);
  F77NAME(zaxpy)(n,complex_double_mone_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_mone_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >* operator-(
const TridiagonalMatrix<double,complex<double> > &T,
const UpperHessenbergMatrix<double,complex<double> > &U) {
  int n=U.size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(U);
  *H*=complex_double_mone_;
  F77NAME(zaxpy)(n,complex_double_one_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(zaxpy)(n-1,complex_double_one_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator-(
const SymmetricMatrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(n-j,double_mone_,M.addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<double> *Mij=M.addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Mij++,Sji+=n) *Sji-=conj(*Mij);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const SymmetricMatrix<double,complex<double> > &M,
const UpperHessenbergMatrix<double,complex<double> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      const complex<double> *Mij=M.addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Mij++,Sji+=n) *Sji=conj(*Mij);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(j+2,n),complex_double_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      F77NAME(zaxpy)(n-j,complex_double_mone_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)-=complex_double_one_;
      if (j<n-1) {
        F77NAME(zaxpy)(n-j-1,complex_double_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const UpperHessenbergMatrix<double,complex<double> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      F77NAME(zaxpy)(min(j+2,n),complex_double_mone_,H.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=complex_double_one_;
      if (j<n-1) {
        F77NAME(zcopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      F77NAME(zaxpy)(min(j+2,n),complex_double_mone_,H.addr(0,j),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,
    complex_double_zero_);
  H->copy(*this);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(j+1,complex_double_mone_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)-=complex_double_one_;
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_mone_,U.addr(0,j),1,
          H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> UpperHessenbergMatrix<double,complex<double> >* operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const UpperHessenbergMatrix<double,complex<double> > &M) {
  int n=M.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<double,complex<double> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<double,complex<double> >(n,double_zero_);
  H->copy(M);
  *H*=complex_double_mone_;
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(j+1,complex_double_one_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)+=complex_double_one_;
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_one_,U.addr(0,j),1,
          H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
    F77NAME(zaxpy)(n,complex_double_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const UpperHessenbergMatrix<double,complex<double> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n,M.addr(0,j),1,S->addr(0,j),1);
    F77NAME(zaxpy)(min(j+2,n),complex_double_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const UpperHessenbergMatrix<double,complex<double> > &H2) const {
  int n=size(0);
  CHECK_SAME(n,H2.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<=min(j+1,n-1);k++) {
      F77NAME(zaxpy)(min(k+2,n),H2(k,j),addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const DiagonalMatrix<double,complex<double> > &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,complex<double> > *S=
    OPERATOR_NEW UpperHessenbergMatrix<double,complex<double> >(n);
  S->copy(*this);
  for (int j=0;j<n;j++) {
    F77NAME(zscal)(min(j+2,n),D[j],S->addr(0,j),1);
  }
  return S;
}

template<> UpperHessenbergMatrix<double,complex<double> >* operator*(
const DiagonalMatrix<double,complex<double> > &D,
const UpperHessenbergMatrix<double,complex<double> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<double,complex<double> > *S=
    OPERATOR_NEW UpperHessenbergMatrix<double,complex<double> >(n);
  S->copy(H);
  for (int i=0;i<n;i++) {
    int j=max(0,i-1);
    F77NAME(zscal)(n-j,D[i],S->addr(i,j),n);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
//compute by columns: note that
// H T e_j = H e_{j-1} T_{j-1,j} + H e_j T_{j,j} + H e_{j+1} T_{j+1,j}
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(zaxpy)(min(j+1,n),T.upperDiagonalValue(j-1),addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(zaxpy)(min(j+2,n),T.diagonalValue(j),addr(0,j),1,
      S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(zaxpy)(min(j+3,n),T.lowerDiagonalValue(j),addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const UpperHessenbergMatrix<double,complex<double> > &H) {
// compute by rows: note that
// e_i^T T H
//   = T_{i,i-1} e_{i-1}^T H + T_{i,i} e_i^T H + T_{i,i+1} e_{i+1}^T H
  int n=H.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int i=0;i<n;i++) {
    if (i>0) {
      int j=max(0,i-2);
      F77NAME(zaxpy)(n-j,T.lowerDiagonalValue(i-1),H.addr(i-1,j),n,
        S->addr(i,j),n);
    }
    int j=max(0,i-1);
    F77NAME(zaxpy)(n-j,T.diagonalValue(i),H.addr(i,j),n,S->addr(i,j),n);
    if (i<n-1) {
      F77NAME(zaxpy)(n-i,T.upperDiagonalValue(i),H.addr(i+1,i),n,
        S->addr(i,i),n);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const TridiagonalMatrix<double,complex<double> > &T) const {
//compute by columns: note that
// H T e_j = H e_{j-1} T_{j-1,j} + H e_j T_{j,j} + H e_{j+1} T_{j+1,j}
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(j-1,0);k<=min(j+1,n-1);k++) {
      F77NAME(zaxpy)(min(k+2,n),T(k,j),addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const TridiagonalMatrix<double,complex<double> > &T,
const UpperHessenbergMatrix<double,complex<double> > &H) {
// compute by rows: note that
// e_i^T T H
//   = T_{i,i-1} e_{i-1}^T H + T_{i,i} e_i^T H + T_{i,i+1} e_{i+1}^T H
  int n=H.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int i=0;i<n;i++) {
    for (int k=max(i-1,0);k<=min(i+1,n-1);k++) {
      int j=max(0,k-1);
      F77NAME(zaxpy)(min(n-k+1,n),T(i,k),H.addr(k,j),n,S->addr(i,j),n);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const SymmetricMatrix<double,complex<double> > &S) const {
// compute by bordering: note that
// [     eta_11 h^T ] [ sigma s^H ]
// [ e_0 eta_21  H  ] [   s    S  ]
//   = [     eta_11 sigma + h^T s ,     eta_11 s^H + h^T S ]
//   = [ e_0 eta_21 sigma +  H  s , e_0 eta_21 s^H +  H  S ]
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  (*M)(n-1,n-1)=(*this)(n-1,n-1)*S(n-1,n-1);
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  complex<double> *v=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-2;k>=0;k--) {
    for (int j=k+1;j<n;j++) { // H s
      F77NAME(zaxpy)(min(n-k-1,j-k+1),S(j,k),addr(k+1,j),1,
        M->addr(k+1,k),1);
    }
    const complex<double> *hkj=addr(k,k+1);
    complex<double> *tj=t;
    for (int j=k+1;j<n;j++,hkj+=n,tj++) *tj=conj(*hkj);
    F77NAME(zhemv)('L',n-k-1,double_one_,S.addr(k+1,k+1),n,
      t,1,double_zero_,v,1);
    complex<double> *vj=v;
    complex<double> *Mkj=M->addr(k,k+1);
    for (int j=k+1;j<n;j++,vj++,Mkj+=n) *Mkj=conj(*vj);
      //h^T S=(S bar(h))^H
    (*M)(k,k)=F77NAME(zdotu)(n-k-1,addr(k,k+1),n,S.addr(k+1,k),1);//h^T s
    F77NAME(zgerc)(2,n-k,complex_double_one_,addr(k,k),1,S.addr(k,k),1,
      M->addr(k,k),n);
  }
  delete t; t=0;
  delete v; v=0;
  return M;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SymmetricMatrix<double,complex<double> > &S,
const UpperHessenbergMatrix<double,complex<double> > &H) {
// compute by bordering: note that
// [ sigma s^H ] [     eta_11 h^T ]
// [   s    S  ] [ e_0 eta_21  H  ]
//   = [ sigma eta_11 + s^H e_0 eta_21 , sigma h^T + s^H H ]
//   = [   s   eta_11 +  S  e_0 eta_21 ,   s   h^T +  S  H ]
  int n=H.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  (*M)(n-1,n-1)=S(n-1,n-1)*H(n-1,n-1);
  for (int k=n-2;k>=0;k--) { // s^T H
    F77NAME(zaxpy)(n-k-1,H(k+1,k),S.addr(k+1,k+1),1,M->addr(k+1,k),1);
      // S e_0 eta_21
    for (int j=k+1;j<n;j++) { // s^H H
      (*M)(k,j)=
        F77NAME(zdotc)(min(n-k-1,j-k+1),S.addr(k+1,k),1,H.addr(k+1,j),1);
    }
    (*M)(k,k)=S(k,k+1)*H(k+1,k); // s^H e_0 eta_21
    F77NAME(zgeru)(n-k,n-k,complex_double_one_,S.addr(k,k),1,
      H.addr(k,k),n,M->addr(k,k),n);
  }
  return M;
}

template<> Matrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *M=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n,double_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<m;k++) {
        F77NAME(zaxpy)(min(k+2,m),L(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(min(j+2,m),addr(0,j),1,M->addr(0,j),1);
      for (int k=j+1;k<m;k++) {
        F77NAME(zaxpy)(min(k+2,m),L(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<double,complex<double> >* operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const UpperHessenbergMatrix<double,complex<double> > &H) {
// compute by columns: note that
//   [ L_11      ] [ h ] = [ L_11 h ]
//   [ L_21 L_22 ] [   ] = [ L_21 h ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,H.size(0));
  Matrix<double,complex<double> > *M=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n,double_zero_);
  char diag=(dynamic_cast<const UnitLowerTrapezoidalMatrix<double,
    complex<double> >*>(&L)==0 ? 'N' : 'U');
  for (int j=0;j<n;j++) {
    int k=min(j+2,n);
    F77NAME(zcopy)(k,H.addr(0,j),1,M->addr(0,j),1);
    F77NAME(ztrmv)('L','N',diag,k,L.addr(),m,M->addr(0,j),1); // L_11 h
    if (k<m) {
      F77NAME(zgemv)('N',m-k,k,complex_double_one_,L.addr(k,0),m,
        H.addr(0,j),1,complex_double_zero_,M->addr(k,j),1); // L_21 h
    }
  }
  return M;
}

template<> Matrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
// compute by columns: note that
//   H [ U_1 , U_2 ] = [ H U_1 , H U_2 ]
// and that
//   [      H_11     H_12 ] [ u ] = [     H_11      u ]
//   [ e_0 eta e_k^T H_22 ] [   ] = [ e_0 eta e_k^T u ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *M=
    OPERATOR_NEW Matrix<double,complex<double> >(m,n,double_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0);
  if (U_non_unit) { // H U_1
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,m-1);k++) {
        F77NAME(zaxpy)(min(k+2,m),U(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,m);k++) {
        F77NAME(zaxpy)(min(k+2,m),U(k,j),addr(0,k),1,M->addr(0,j),1);
      }
      if (j<m) {
        F77NAME(zaxpy)(min(j+2,m),double_one_,addr(0,j),1,M->addr(0,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<double,complex<double> >* operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const UpperHessenbergMatrix<double,complex<double> > &H) {
// compute by columns: note that
//   [ U_1 , U_2 ] [     H_11          , H_12 ]
//                 [ e_0 eta e_{m-1}^T , H_22 ]
//   = [ U_1 H_11 + U_2 e_0 eta e_{m-1}^T , U_1 H_12 + U_2 H_22 ]
// and that
//   [ U_11 U_12 ] [ h ] = [ U_11 h ]
//   [      U_22 ] [   ] = [        ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,H.size(0));
  char diag=(dynamic_cast<const UnitUpperTrapezoidalMatrix<double,
    complex<double> >*>(&U)==0 ? 'N' : 'U');
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<m;j++) { // U_1 H_11
    int k=min(j+2,m);
    F77NAME(zcopy)(k,H.addr(0,j),1,M->addr(0,j),1);
    F77NAME(ztrmv)('U','N',diag,k,U.addr(),m,M->addr(0,j),1); // U_11 h
  }
  if (m<n) {
    F77NAME(zaxpy)(m,H(m,m-1),U.addr(0,m),1,M->addr(0,m-1),1);
      // U_2 e_0 eta e_{m-1}^T
    F77NAME(zlacpy)('A',m,n-m,H.addr(0,m),m,M->addr(0,m),m);
    F77NAME(ztrmm)('L','U','N',diag,m,n-m,complex_double_one_,U.addr(),m,
      M->addr(0,m),m); // U_1 H_12
    F77NAME(zgemm)('N','N',m,n-m,n-m,complex_double_one_,U.addr(0,m),m,
      H.addr(m,m),n,complex_double_one_,M->addr(0,m),m); // U_2 H_22
  }
  return M;
}

template<> SquareMatrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const SquareMatrix<double,complex<double> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<n;k++) {
      F77NAME(zaxpy)(min(k+2,n),S(k,j),addr(0,k),1,M->addr(0,j),1);
    }
  }
  return M;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SquareMatrix<double,complex<double> > &S,
const UpperHessenbergMatrix<double,complex<double> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *M=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(zgemv)('N',n,min(j+2,n),complex_double_one_,S.addr(),n,
      H.addr(0,j),1,complex_double_zero_,M->addr(0,j),1);
  }
  return M;
}

template<> Matrix<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<double,complex<double> > *R=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<m;k++) {
      F77NAME(zaxpy)(min(m,k+2),M(k,j),addr(0,k),1,R->addr(0,j),1);
    }
  }
  return R;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const UpperHessenbergMatrix<double,complex<double> > &H) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,H.size(0));
  Matrix<double,complex<double> > *R=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<=min(n-1,j+1);k++) {
      F77NAME(zaxpy)(m,H(k,j),M.addr(0,k),1,R->addr(0,j),1);
    }
  }
  return R;
}

template<> Vector<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<double,complex<double> > *r=
    OPERATOR_NEW Vector<double,complex<double> >(n,complex_double_zero_);
  for (int k=0;k<n;k++) {
    F77NAME(zaxpy)(min(n,k+2),v[k],addr(0,k),1,r->addr(),1);
  }
  return r;
}

template<> void UpperHessenbergMatrix<double,complex<double> >::uhmv(
complex<double> alpha,const Vector<double,complex<double> > &x,
complex<double> beta,Vector<double,complex<double> > &b,char trans)
const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  if (abs(beta)==double_zero_) b=complex_double_zero_;
  else b*=beta;
  if (trans=='N' || trans=='n') { // b=H*x*alpha+b*beta
    complex<double> *xj=x.addr();
    for (int j=0;j<n;j++,xj++) {
      F77NAME(zaxpy)(min(n,j+2),(*xj)*alpha,addr(0,j),1,b.addr(),1); 
    }
  } else if (trans=='T' || trans=='t') { // b=H^T*x*alpha+b*beta
    complex<double> *bj=b.addr();
    for (int j=0;j<n;j++,bj++) {
      *bj+=F77NAME(zdotu)(min(n,j+2),addr(0,j),1,x.addr(),1)*alpha;
    }
  } else { // b=H^H*x*alpha+b*beta
    complex<double> *bj=b.addr();
    for (int j=0;j<n;j++,bj++) {
      *bj+=F77NAME(zdotc)(min(n,j+2),addr(0,j),1,x.addr(),1)*alpha;
    }
  }
}

template<> void UpperHessenbergMatrix<double,complex<double> >::uhmm(
complex<double> alpha,const Matrix<double,complex<double> > &X,
complex<double> beta,Matrix<double,complex<double> > &B,char side,
char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (abs(beta)==double_zero_) B=complex_double_zero_;
  else B*=beta;
  if (side=='L' || side=='l') {
    CHECK_SAME(m,size(0));
    if (trans=='N' || trans=='n') { // B=H*X*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int k=0;k<m;k++) {
          F77NAME(zaxpy)(min(m,k+2),X(k,j)*alpha,addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else if (trans=='T' || trans=='t') { // B=H^T*X*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int i=0;i<m;i++) {
          B(i,j)+=F77NAME(zdotu)(min(m,i+2),addr(0,i),1,X.addr(0,j),1)
                *alpha;
        }
      }
    } else { // B=H^H*X*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int i=0;i<m;i++) {
          B(i,j)+=F77NAME(zdotc)(min(m,i+2),addr(0,i),1,X.addr(0,j),1)
                *alpha;
        }
      }
    }
  } else {
    CHECK_SAME(n,size(0));
    if (trans=='N' || trans=='n') { // B=X*H*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int k=0;k<=min(n-1,j+1);k++) {
          F77NAME(zaxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else if (trans=='T' || trans=='t') { // // B=X*H^T*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<n;k++) {
          F77NAME(zaxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else { // B=X*H^H*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<n;k++) {
          F77NAME(zaxpy)(m,conj((*this)(j,k))*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

template<> double
UpperHessenbergMatrix<double,complex<double> >::normFrobenius() const {
  int n=size(0);
  double *work=0;
  return F77NAME(zlanhs)('F',n,addr(),n,work);
}

template<> double
UpperHessenbergMatrix<double,complex<double> >::normInfinity() const {
  int n=size(0);
  double *work=OPERATOR_NEW_BRACKET(double,n);
  double val=F77NAME(zlanhs)('I',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> double
UpperHessenbergMatrix<double,complex<double> >::normMaxEntry() const {
  int n=size(0);
  double *work=0;
  return F77NAME(zlanhs)('M',n,addr(),n,work);
}

template<> double
UpperHessenbergMatrix<double,complex<double> >::normOne() const {
  int n=size(0);
  double *work=0;
  return F77NAME(zlanhs)('O',n,addr(),n,work);
}

template<> Vector<double,complex<double> >*
UpperHessenbergMatrix<double,complex<double> >::eigenvalues(
SquareMatrix<double,complex<double> > *&V,
SquareMatrix<double,complex<double> > *&U) const {
  int n=size(0);
  if (V!=0) CHECK_SAME(n,V->size(0));
  if (U!=0) CHECK_SAME(n,U->size(0));
  UpperHessenbergMatrix<double,complex<double> > *H=
    OPERATOR_NEW UpperHessenbergMatrix<double,complex<double> >(*this);
  char job=(V==0 && U==0 ? 'E' : 'S');
  char compz=(V==0 && U==0 ? 'N' : 'I');
  Vector<double,complex<double> > *lambda =
    OPERATOR_NEW Vector<double,complex<double> >(n);
  OrthogonalMatrix<double,complex<double> > *Z=(V==0 && U==0 ? 0 :
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(n,n));
  complex<double> *za=(Z==0 ? 0 : Z->addr());
  complex<double> ww(numeric_limits<double>::infinity(),
    numeric_limits<double>::infinity());
  int lwork=-1;
  int info;
  F77NAME(zhseqr)(job,compz,n,1,n,H->addr(),n,lambda->addr(),za,n,&ww,
    lwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(ww.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zhseqr)(job,compz,n,1,n,H->addr(),n,lambda->addr(),za,n,work,
    lwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;

  if (Z!=0) {
    char side=(V==0 ? 'R' : (U==0 ? 'L' : 'B') );
    bool *select=0;
    complex<double> *vla=0;
//  SquareMatrix<double,complex<double> > *Vl=0;
    if (V!=0) {
//    Vl=OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
      V->copy(*Z);
      vla=V->addr();
    }
    complex<double> *vra=0;
//  SquareMatrix<double,complex<double> > *Vr=0;
    if (U!=0) {
//    Vr=OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
      U->copy(*Z);
      vra=U->addr();
    }
    work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
    int mout=-1;
    double *rwork=OPERATOR_NEW_BRACKET(double,n);
    F77NAME(ztrevc)(side,'B',select,n,H->addr(),n,vla,n,vra,n,n,mout,work,
      rwork,info);
    CHECK_TEST(info==0);
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
    delete [] rwork; rwork=0;
    delete [] work; work=0;
//  if (Vl!=0) delete Vl; Vl=0;
//  if (Vr!=0) delete Vr; Vr=0;
  }
  if (Z!=0) delete Z; Z=0;
  delete H; H=0;
  return lambda;
}  

template<> int UpperHessenbergMatrix<double,complex<double> >::factor(
int *ipiv) {
  int n=size(0);
  for (int k=0;k<n;k++) ipiv[k]=k;
  int info=0;
  for (int k=0;k<n-1;k++) {
    if (abs((*this)(k,k))<abs((*this)(k+1,k)) ) {
      ipiv[k]=k+1;
      F77NAME(zswap)(n-k,addr(k,k),n,addr(k+1,k),n);
    }
    if (abs((*this)(k,k))>double_zero_) {
      (*this)(k+1,k)/=(*this)(k,k);
      F77NAME(zaxpy)(n-k-1,-(*this)(k+1,k),addr(k,k+1),n,
        addr(k+1,k+1),n);
    } else {
      info=k+1;
      break;
    }
  }
  return info;
}

template<> void UpperHessenbergMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  int n=size(0);
  CHECK_SAME(n,b.size());
  CHECK_SAME(n,x.size());
  x.copy(b);
  UpperHessenbergMatrix<double,complex<double> > *HF=
    OPERATOR_NEW UpperHessenbergMatrix<double,complex<double> >(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int info=HF->factor(ipiv);
  CHECK_TEST(info==0);
  if (trans!='N' && trans!='n') {
  // H = Q^T L U and H^H x = b ==> U^T L^T Q conj(x) = conj(b)
  //              or H^T x = b ==> U^T * L^T Q x = b
    if (trans=='C' || trans=='c') {
      for (int j=0;j<n;j++) x[j]=conj(x[j]);
    }
    x[0]/=(*HF)(0,0); // forward-solve U^T y = conj(b);
    for (int j=1;j<n;j++) {
      x[j]=(x[j]-F77NAME(zdotu)(j,HF->addr(0,j),1,x.addr(),1))
          /(*HF)(j,j);
    }
    for (int i=n-2;i>=0;i--) { // back-solve L^T Q z = y
      int ip=ipiv[i];
      complex<double> temp=x[i]-(*HF)(i+1,i)*x[i+1];
      x[i]=x[ip];
      x[ip]=temp;
    }
    if (trans=='C' || trans=='c') {
      for (int j=0;j<n;j++) x[j]=conj(x[j]); // x = conj(y)
    }
  } else { // H = Q^T L U and H x = b ==> Q^T L U x = b
    for (int i=0;i<n-1;i++) { // forward-solve Q^T L  = b
      int ip=ipiv[i];
      complex<double> temp=x[2*i-ip+1]-(*HF)(i+1,i)*x[ip];
      x[i]=x[ip];
      x[i+1]=temp;
    }
    for (int j=n-1;j>0;j--) { // back-solve U x = y
      x[j]/=(*HF)(j,j);
      F77NAME(zaxpy)(j,-x[j],HF->addr(0,j),1,x.addr(),1);
    }
    x[0]/=(*HF)(0,0);
  }
  delete ipiv; ipiv=0;
  delete HF; HF=0;
}

template<> void UpperHessenbergMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  int n=size(0);
  UpperHessenbergMatrix<double,complex<double> > *HF=
    OPERATOR_NEW UpperHessenbergMatrix<double,complex<double> >(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int info=HF->factor(ipiv);
  CHECK_TEST(info==0);
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(n,X.size(0));
    CHECK_SAME(nrhs,X.size(1));
    X.copy(B);
    if (trans!='N' && trans!='n') {
    // H = Q^T L U and H^H X = B ==> U^T L^T Q conj(X) = conj(B)
    //              or H^T X = B ==> U^T L^T Q X = B
      for (int k=0;k<nrhs;k++) {
        complex<double> *X0k=X.addr(0,k);
        complex<double> *Xjk=X0k;
        if (trans=='C' || trans=='c') {
          for (int j=0;j<n;j++,Xjk+=nrhs) *Xjk=conj(*Xjk);
        }
        *X0k/=(*HF)(0,0); // forward-solve U^T Y = conj(B)
        Xjk=X0k+nrhs;
        for (int j=1;j<n;j++,Xjk+=nrhs) {
          *Xjk=(*Xjk-F77NAME(zdotu)(j,HF->addr(0,j),1,X0k,1))/(*HF)(j,j);
        }
        for (int i=n-2;i>=0;i--) { // back-solve L^T Q Z = Y
          int ip=ipiv[i];
          complex<double> temp=X(i,k)-(*HF)(i+1,i)*X(i+1,k);
          X(i,k)=X(ip,k);
          X(ip,k)=temp;
        }
        if (trans=='C' || trans=='c') {
          Xjk=X0k;
          for (int j=0;j<n;j++,Xjk+=nrhs) *Xjk=conj(*Xjk); // X = conj(Z)
        }
      }
    } else {
    // H = Q^T L U and H X = B ==> Q^T L U X = B
      for (int k=0;k<nrhs;k++) {
        for (int i=0;i<n-1;i++) { // forward-solve Q^T L Y = B
          int ip=ipiv[i];
          complex<double> temp=X(2*i-ip+1,k)-(*HF)(i+1,i)*X(ip,k);
          X(i,k)=X(ip,k);
          X(i+1,k)=temp;
        }
        for (int j=n-1;j>0;j--) { // back-solve U X = Y
          X(j,k)/=(*HF)(j,j);
          F77NAME(zaxpy)(j,-X(j,k),HF->addr(0,j),1,X.addr(0,k),1);
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
    // H = Q^T L U and X H^H = B ==> Q^T L U X^H = B^H
    //              or X H^T = B ==> Q^T L U X^T = B^T
      for (int k=0;k<nrhs;k++) {
        complex<double> *Xk0=X.addr(k,0);
        complex<double> *Xki=Xk0;
        if (trans=='C' || trans=='c') {
          for (int i=0;i<n;i++,Xki+=nrhs) *Xki=conj(*Xki);
        }
        for (int i=0;i<n-1;i++) { // forward-solve Q^T L Y = B^H
          int ip=ipiv[i];
          complex<double> temp=X(k,2*i-ip+1)-(*HF)(i+1,i)*X(k,ip);
          X(k,i)=X(k,ip);
          X(k,i+1)=temp;
        }
        complex<double> *Xkj=X.addr(k,n-1);
        for (int j=n-1;j>0;j--,Xkj-=nrhs) { // back-solve U Z = Y
          *Xkj/=(*HF)(j,j);
          F77NAME(zaxpy)(j,-(*Xkj),HF->addr(0,j),1,Xk0,nrhs);
        }
        (*Xk0)/=(*HF)(0,0);
        if (trans=='C' || trans=='c') {
          Xki=Xk0;
          for (int i=0;i<n;i++,Xki+=nrhs) *Xki=conj(*Xki); // X^H = Z;
        }
      }
    } else {
    // H = Q^T L U and X H = B ==> U^T L^T Q X^T = B^T
      for (int k=0;k<nrhs;k++) {
        complex<double> *Xk0=X.addr(k,0);
        complex<double> *Xkj=Xk0;
        *Xkj/=(*HF)(0,0); // forward-solve U^T Y = B^T
        for (int j=1;j<n;j++,Xkj+=nrhs) {
          *Xkj=(*Xkj-F77NAME(zdotu)(j,HF->addr(0,j),1,Xk0,nrhs))
              /(*HF)(j,j);
        }
        for (int i=n-2;i>=0;i--) { // back-solve L^T Q X^T = Y
          int ip=ipiv[i];
          complex<double> temp=X(k,i)-(*HF)(i+1,i)*X(k,i+1);
          X(k,i)=X(k,ip);
          X(k,ip)=temp;
        }
      }
    }
  }
  delete ipiv; ipiv=0;
  delete HF; HF=0;
}

template class UpperHessenbergMatrix<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const complex<double>
  BandMatrix<double,complex<double> >::outofbounds_=complex_double_zero_;
template<> complex<double> BandMatrix<double,complex<double> >::safety_ =
  complex_double_zero_;

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::makeMatrix() const {
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      M->addr(ibeg,j),1);
  }
  return M;
}

template<> BandMatrix<double,complex<double> >&
BandMatrix<double,complex<double> >::operator+=(
const BandMatrix<double,complex<double> > &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  CHECK_SAME(nsup,B.nsup);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,complex_double_one_,
      B.addr(ibeg,j),1,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<double,complex<double> >&
BandMatrix<double,complex<double> >::operator-=(
const BandMatrix<double,complex<double> > &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  CHECK_SAME(nsup,B.nsup);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<double,complex<double> >&
BandMatrix<double,complex<double> >::operator*=(complex<double> d) {
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zscal)(min(dim-1,j+nsub)-ibeg+1,d,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<double,complex<double> >&
BandMatrix<double,complex<double> >::operator/=(complex<double> d) {
  CHECK_TEST(abs(d)>double_zero_);
  complex<double> dinv=complex_double_one_/d;
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zscal)(min(dim-1,j+nsub)-ibeg+1,dinv,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator+(
const BandMatrix<double,complex<double> > &B) const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,B.nsub),
    max(nsup,B.nsup),complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    ibeg=max(0,j-B.nsup);
    F77NAME(zaxpy)(min(dim-1,j+B.nsub)-ibeg+1,complex_double_one_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* 
BandMatrix<double,complex<double> >::operator+(
const UpperHessenbergMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    F77NAME(zaxpy)(min(dim,j+2),complex_double_one_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator+(
const DiagonalMatrix<double,complex<double> > &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<double,complex<double> > *S=
    OPERATOR_NEW BandMatrix<double,complex<double> >(*this);
  F77NAME(zaxpy)(dim,complex_double_one_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator+(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,1),max(nsup,1),
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(zaxpy)(dim-1,complex_double_one_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),nt);
  complex<double> *Sjj=S->addr(0,0);
  const double *Tjj=T.diagonalAddr(0);
  for (int j=0;j<dim;j++,Sjj+=S->nt,Tjj++) *Sjj+=*Tjj;
  complex<double> *Sjjp1=S->addr(0,1);
  const complex<double> *Tjm1j=T.lowerDiagonalAddr(0);
  for (int j=0;j<dim-1;j++,Sjjp1+=S->nt,Tjm1j++) *Sjjp1+=conj(*Tjm1j);
  return S;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator+(
const TridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,1),max(nsup,1),
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(zaxpy)(dim,complex_double_one_,T.diagonalAddr(),1,
    S->addr(0,0),S->nt);
  F77NAME(zaxpy)(dim-1,complex_double_one_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S->nt);
  F77NAME(zaxpy)(dim-1,complex_double_one_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S->nt);
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator+(
const SymmetricMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<dim;j++) {
    F77NAME(zaxpy)(dim-j,complex_double_one_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(zaxpy)(dim-j-1,complex_double_one_,H.addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(zaxpy)(j+1,complex_double_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(zaxpy)(dim-j,complex_double_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      (*S)(j,j)+=complex_double_one_;
      if (j+1<dim) {
        F77NAME(zaxpy)(dim-j-1,complex_double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  S->copy(M);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,complex_double_one_,
      addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const BandMatrix<double,complex<double> > &B) const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,B.nsub),
    max(nsup,B.nsup),complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    ibeg=max(0,j-B.nsup);
    F77NAME(zaxpy)(min(B.dim-1,j+B.nsub)-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const UpperHessenbergMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    F77NAME(zaxpy)(min(dim,j+2),complex_double_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const UpperHessenbergMatrix<double,complex<double> > &H,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(n,j+2),H.addr(0,j),1,S->addr(0,j),1);
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const DiagonalMatrix<double,complex<double> > &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<double,complex<double> > *S=
    OPERATOR_NEW BandMatrix<double,complex<double> >(*this);
  F77NAME(zaxpy)(dim,complex_double_mone_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> BandMatrix<double,complex<double> >* operator-(
const DiagonalMatrix<double,complex<double> > &D,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,B.subDiags(),B.supDiags(),
    complex_double_zero_);
  F77NAME(zcopy)(n,D.addr(),1,S->addr(0,0),S->bands());
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,1),max(nsup,1),
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(zaxpy)(dim-1,complex_double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  complex<double> *Sjj=S->addr(0,0);
  const double *Tjj=T.diagonalAddr(0);
  for (int j=0;j<dim;j++,Sjj+=S->nt,Tjj++) *Sjj-=*Tjj;
  complex<double> *Sjjp1=S->addr(0,1);
  const complex<double> *Tjp1j=T.lowerDiagonalAddr(0);
  for (int j=0;j<dim-1;j++,Sjjp1+=S->nt,Tjp1j++) *Sjjp1+=conj(*Tjp1j);
  return S;
}

template<> BandMatrix<double,complex<double> >* operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,max(B.subDiags(),1),
    max(B.supDiags(),1),complex_double_zero_);
  int nb=S->bands();
  F77NAME(zcopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(1,0),nb);
  complex<double> *Sjj=S->addr(0,0);
  const double *Tjj=T.diagonalAddr(0);
  for (int j=0;j<n;j++,Sjj+=nb,Tjj++) *Sjj=*Tjj;
  complex<double> *Sjjp1=S->addr(0,1);
  const complex<double> *Tjp1j=T.lowerDiagonalAddr(0);
  for (int j=0;j<n-1;j++,Sjjp1+=nb,Tjp1j++) *Sjjp1=conj(*Tjp1j);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const TridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,1),max(nsup,1),
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(zaxpy)(dim,complex_double_mone_,T.diagonalAddr(),1,
    S->addr(0,0),nt);
  F77NAME(zaxpy)(dim-1,complex_double_mone_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),nt);
  F77NAME(zaxpy)(dim-1,complex_double_mone_,T.upperDiagonalAddr(),1,
    S->addr(0,1),nt);
  return S;
}

template<> BandMatrix<double,complex<double> >* operator-(
const TridiagonalMatrix<double,complex<double> > &T,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,max(B.subDiags(),1),
    max(B.supDiags(),1),complex_double_zero_);
  F77NAME(zcopy)(n,T.diagonalAddr(),1,S->addr(0,0),S->bands());
  F77NAME(zcopy)(n-1,T.lowerDiagonalAddr(),1,S->addr(1,0),S->bands());
  F77NAME(zcopy)(n-1,T.upperDiagonalAddr(),1,S->addr(0,1),S->bands());
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const SymmetricMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<dim;j++) {
    F77NAME(zaxpy)(dim-j,complex_double_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim) {
      complex<double> *Sji=S->addr(j,j+1);
      const complex<double> *Hij=H.addr(j+1,j);
      for (int i=j+1;i<dim;i++,Sji+=dim,Hij++) *Sji-=conj(*Hij);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const SymmetricMatrix<double,complex<double> > &H,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,H.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      complex<double> *Sji=S->addr(j,j+1);
      const complex<double> *Hij=H.addr(j+1,j);
      for (int i=j+1;i<n;i++,Sji+=n,Hij++) *Sji=conj(*Hij);
    }
  }
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(zaxpy)(j+1,complex_double_mone_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_mone_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)-=complex_double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(j+1,U.addr(0,j),1,S->addr(0,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-B.supDiags());
      if (j>0) F77NAME(zcopy)(j,U.addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)=complex_double_one_;
      F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(zaxpy)(dim-j,complex_double_mone_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      (*S)(j,j)-=complex_double_one_;
      if (j+1<dim) {
        F77NAME(zaxpy)(dim-j-1,complex_double_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=complex_double_one_;
      if (j+1<n) F77NAME(zcopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  S->operator-=(M);
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  S->copy(M);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zaxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const BandMatrix<double,complex<double> > &B) const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,nsub+B.nsub,nsup+B.nsup,
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=max(0,j-B.nsup);k<=min(dim-1,j+B.nsub);k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,B(k,j),
        addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const UpperHessenbergMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=0;k<=min(dim-1,j+1);k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,H(k,j),
        addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const UpperHessenbergMatrix<double,complex<double> > &H,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(zaxpy)(min(n,k+2),B(k,j),H.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const DiagonalMatrix<double,complex<double> > &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<double,complex<double> > *P=
    OPERATOR_NEW BandMatrix<double,complex<double> >(*this);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(zscal)(min(dim-1,j+nsub)-ibeg+1,D[j],
      P->addr(ibeg,j),1);
  }
  return P;
}

template<> BandMatrix<double,complex<double> >* operator*(
const DiagonalMatrix<double,complex<double> > &D,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,B.subDiags(),B.supDiags());
  P->copy(B);
  int stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    int jbeg=max(0,i-B.subDiags());
    int jend=min(n-1,i+B.supDiags());
    F77NAME(zscal)(min(n-1,i+B.supDiags())-jbeg+1,D[i],
      P->addr(i,jbeg),stride);
  }
  return P;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,nsub+1,nsup+1,
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    if (j>0) {
      int ibeg=max(0,j-1-nsup);
      F77NAME(zaxpy)(min(dim-1,nsub+j-1)-ibeg+1,
        T.upperDiagonalValue(j-1),addr(ibeg,j-1),1,P->addr(ibeg,j),1);
    }
    int ibeg=max(0,j-nsup);
    F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,T.diagonalValue(j),
      addr(ibeg,j),1,P->addr(ibeg,j),1);
    if (j<dim-1) {
      int ibeg=max(0,j+1-nsup);
      F77NAME(zaxpy)(min(dim-1,j+1+nsub)-ibeg+1,
        T.lowerDiagonalValue(j),addr(ibeg,j+1),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> BandMatrix<double,complex<double> >* operator*(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,B.subDiags()+1,B.supDiags()+1,
    complex_double_zero_);
  int B_stride=B.bands()-1;
  int P_stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    if (i>0) {
      int jbeg=max(0,i-1-B.subDiags());
      F77NAME(zaxpy)(min(n-1,i-1+B.supDiags())-jbeg+1,
        T.lowerDiagonalValue(i-1),B.addr(i-1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
    int jbeg=max(0,i-B.subDiags());
    F77NAME(zaxpy)(min(n-1,i+B.supDiags())-jbeg+1,T.diagonalValue(i),
      B.addr(i,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    if (i<n-1) {
      int jbeg=max(0,i+1-B.subDiags());
      F77NAME(zaxpy)(min(n-1,i+1+B.supDiags())-jbeg+1,
        T.upperDiagonalValue(i),B.addr(i+1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
  }
  return P;
}

template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const TridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,nsub+1,nsup+1,
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    if (j>0) {
      int ibeg=max(0,j-1-nsup);
      F77NAME(zaxpy)(min(dim-1,nsub+j-1)-ibeg+1,T(j-1,j),
        addr(ibeg,j-1),1,P->addr(ibeg,j),1);
    }
    int ibeg=max(0,j-nsup);
    F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,T(j,j),addr(ibeg,j),1,
      P->addr(ibeg,j),1);
    if (j<dim-1) {
      int ibeg=max(0,j+1-nsup);
      F77NAME(zaxpy)(min(dim-1,j+1+nsub)-ibeg+1,T(j+1,j),
        addr(ibeg,j+1),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> BandMatrix<double,complex<double> >* operator*(
const TridiagonalMatrix<double,complex<double> > &T,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,B.subDiags()+1,B.supDiags()+1,
    complex_double_zero_);
  int B_stride=B.bands()-1;
  int P_stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    if (i>0) {
      int jbeg=max(0,i-1-B.subDiags());
      F77NAME(zaxpy)(min(n-1,i-1+B.supDiags())-jbeg+1,T(i,i-1),
        B.addr(i-1,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    }
    int jbeg=max(0,i-B.subDiags());
    F77NAME(zaxpy)(min(n-1,i+B.supDiags())-jbeg+1,T(i,i),
      B.addr(i,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    if (i<n-1) {
      int jbeg=max(0,i+1-B.subDiags());
      F77NAME(zaxpy)(min(n-1,i+1+B.supDiags())-jbeg+1,
        T(i,i+1),B.addr(i+1,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    }
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const SymmetricMatrix<double,complex<double> > &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=0;k<dim;k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,S(k,j),addr(ibeg,k),1,
        P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SymmetricMatrix<double,complex<double> > &S,
const BandMatrix<double,complex<double> > &B) {
//compute by bordering: note that
//  [ sigma s^H ] [ beta c^T ]
//  [   s    S  ] [   b   B  ]
//  = [ sigma beta + s^H b , sigma c^T + s^H B ]
//  = [   s   beta +  S  b ,   s  c^T +   S  B ]
  int n=B.size(0),nsub=B.subDiags(),nsup=B.supDiags();
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  complex<double> *t=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int k=n-1;k>=0;k--) {
    int iend=min(n-1,k+nsub);
    if (k<n-1) {
      // s^H B = ( B^H s )^H
      F77NAME(zgbmv)('C',n-k-1,n-k-1,nsub,nsup,
        complex_double_one_,B.addr(k+1-nsup,k+1),B.bands(),
        S.addr(k+1,k),1,complex_double_zero_,t+k+1,1);
      for (int j=k+1;j<n;j++) (*P)(k,j)=conj(t[j]);
      // S b: note that
      // [ S_11 , S_21^H ] [ b ] = [ S_11 b ]
      // [ S_21 ,  S_22  ] [ 0 ] = [ S_21 b ]
      if (nsub>0) { // S_11 b
        F77NAME(zhemv)('L',min(n-k-1,nsub),complex_double_one_,
          S.addr(k+1,k+1),n,B.addr(k+1,k),1,complex_double_zero_,
          P->addr(k+1,k),1);
        if (nsub<n-k-1) { // S_21 b
          F77NAME(zgemv)('N',n-k-1-nsub,nsub,complex_double_one_,
            S.addr(k+1+nsub,k+1),n,B.addr(k+1,k),1,complex_double_zero_,
            P->addr(k+1+nsub,k),1);
        }
      }
    }
    if (iend>k) { // s^H b
      (*P)(k,k)=F77NAME(zdotc)(iend-k,S.addr(k+1,k),1,B.addr(k+1,k),1);
    }
    int jend=min(n-1,k+nsup);
    F77NAME(zgeru)(n-k,jend-k+1,double_one_,S.addr(k,k),1,
      B.addr(k,k),B.bands()-1,P->addr(k,k),n);
  }
  delete [] t; t=0;
  return P;
}

template<> Matrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  CHECK_SAME(dim,U.size(0));
  int n=U.size(1);
  Matrix<double,complex<double> > *M=
    OPERATOR_NEW Matrix<double,complex<double> >(dim,U.size(1),
    complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,dim-1);k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,dim);k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
      if (j<dim) {
        int ibeg=max(0,j-nsup);
        F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,double_one_,
          addr(ibeg,j),1,M->addr(ibeg,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<double,complex<double> >* operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const BandMatrix<double,complex<double> > &B) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(zaxpy)(min(k+1,m),B(k,j),U.addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(zaxpy)(min(k,m),B(k,j),U.addr(0,k),1,M->addr(0,j),1);
        if (k<m) (*M)(k,j)+=B(k,j);
      }
    }
  }
  return M;
}

template<> Matrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  CHECK_SAME(dim,L.size(0));
  int n=L.size(1);
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(dim,L.size(1),complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<dim;k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,complex_double_one_,
        addr(ibeg,j),1,M->addr(ibeg,j),1);
      for (int k=j+1;k<dim;k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<double,complex<double> >* operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const BandMatrix<double,complex<double> > &B) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,complex<double> > *M=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(zaxpy)(m-k,B(k,j),L.addr(k,k),1,M->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        (*M)(k,j)+=B(k,j);
        if (k+1<m) {
          F77NAME(zaxpy)(m-k-1,B(k,j),L.addr(k+1,k),1,M->addr(k+1,j),1);
        }
      }
    }
  }
  return M;
}

template<> SquareMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const SquareMatrix<double,complex<double> > &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zgbmv)('N',dim,dim,nsub,nsup,complex_double_one_,addr(),nt,
      S.addr(0,j),1,complex_double_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SquareMatrix<double,complex<double> > &S,
const BandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(zaxpy)(n,B(k,j),S.addr(0,j),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  CHECK_SAME(dim,M.size(0));
  int n=M.size(1);
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(dim,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zgbmv)('N',dim,dim,nsub,nsup,complex_double_one_,addr(),nt,
      M.addr(0,j),1,complex_double_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const BandMatrix<double,complex<double> > &B) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(zaxpy)(m,B(k,j),M.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Vector<double,complex<double> >*
BandMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  CHECK_SAME(dim,v.size());
  Vector<double,complex<double> > *w=
    OPERATOR_NEW Vector<double,complex<double> >(dim,double_zero_);
  F77NAME(zgbmv)('N',dim,dim,nsub,nsup,complex_double_one_,addr(),nt,
    v.addr(),1,complex_double_zero_,w->addr(),1);
  return w;
}

template<> void BandMatrix<double,complex<double> >::gbmv(
complex<double> alpha,const Vector<double,complex<double> > &x,
complex<double> beta,Vector<double,complex<double> > &b,char trans)
const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  F77NAME(zgbmv)(trans,dim,dim,nsub,nsup,alpha,addr(),nt,x.addr(),1,
    beta,b.addr(),1);
}

template<> void BandMatrix<double,complex<double> >::gbmm(
complex<double> alpha,const Matrix<double,complex<double> > &X,
complex<double> beta,Matrix<double,complex<double> > &B,char side,
char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') { // B=op(A)*X*alpha+B*beta
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      F77NAME(zgbmv)(trans,dim,dim,nsub,nsup,alpha,addr(),nt,
        X.addr(0,j),1,beta,B.addr(0,j),1);
    }
  } else { // B=X*op(A)*alpha+B*beta
    CHECK_SAME(n,dim);
    if (abs(beta)==double_zero_) B=double_zero_;
    else B*=beta;
    if (trans=='N' || trans=='n') {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsup);k<=min(dim-1,j+nsub);k++) {
          F77NAME(zaxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else if (trans=='C' || trans=='c') {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsub);k<=min(dim-1,j+nsup);k++) {
          F77NAME(zaxpy)(m,conj((*this)(j,k))*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsub);k<=min(dim-1,j+nsup);k++) {
          F77NAME(zaxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

/*
template<> BandMatrix<double,complex<double> >*
BandMatrix<double,complex<double> >::transpose() const {
  BandMatrix<double,complex<double> > *X=
    OPERATOR_NEW BandMatrix<double,complex<double> >(dim,nsup,nsub);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsub);
    F77NAME(dcopy)(min(dim-1,j+nsup)-i+1,addr(j,i),nsub+nsup,
      X->addr(i,j),1);
  }
  return X;
}
*/

template<> double BandMatrix<double,complex<double> >::equilibrate(
Vector<double,double> &r,Vector<double,double> &c,double &rowcnd,
double &colcnd) const {
  CHECK_SAME(dim,r.size());
  CHECK_SAME(dim,c.size());
  double amax;
  int info;
  F77NAME(zgbequ)(dim,dim,nsub,nsup,addr(),nt,r.addr(),c.addr(),
    rowcnd,colcnd,amax,info);
  CHECK_TEST(info==0);
  return amax;
}

template<> double BandMatrix<double,complex<double> >::normFrobenius()
const {
  double *work=0;
  return F77NAME(zlangb)('F',dim,nsub,nsup,addr(),nt,work);
}

template<> double BandMatrix<double,complex<double> >::normInfinity()
const {
  double *work=OPERATOR_NEW_BRACKET(double,dim);
  double val=F77NAME(zlangb)('I',dim,nsub,nsup,addr(),nt,work);
  delete work;
  return val;
}

template<> double BandMatrix<double,complex<double> >::normMaxEntry()
const {
  double *work=0;
  return F77NAME(zlangb)('M',dim,nsub,nsup,addr(),nt,work);
}

template<> double BandMatrix<double,complex<double> >::normOne() const {
  double *work=0;
  return F77NAME(zlangb)('O',dim,nsub,nsup,addr(),nt,work);
}

template<> double
BandMatrix<double,complex<double> >::reciprocalConditionNumber(char norm)
const {
  BandMatrix<double,complex<double> > *BF=
    OPERATOR_NEW BandMatrix<double,complex<double> >(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(zgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);

  double anorm=(norm=='I' || norm=='i' ? normInfinity() : normOne());
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*dim);
  double *rwork=OPERATOR_NEW_BRACKET(double,dim);
  F77NAME(zgbcon)(norm,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,anorm,
    rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete [] rwork; rwork=0;
  delete [] ipiv; ipiv=0;
  delete BF; BF=0;
  return rcond;
}

template<> void BandMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char trans) const {
  CHECK_SAME(dim,b.size())
  CHECK_SAME(dim,x.size())
  BandMatrix<double,complex<double> > *BF=
    OPERATOR_NEW BandMatrix<double,complex<double> >(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(zgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);

  x.copy(b);
  F77NAME(zgbtrs)(trans,dim,nsub,nsup,1,BF->addr(),BF->nt,ipiv,
    x.addr(),dim,info);
  CHECK_TEST(info==0);
  delete [] ipiv;
  delete BF; BF=0;
}

template<> void BandMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char trans) const {
  BandMatrix<double,complex<double> > *BF=
    OPERATOR_NEW BandMatrix<double,complex<double> >(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(zcopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(zgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(dim,B.size(0))
    CHECK_SAME(dim,X.size(0))
    CHECK_SAME(nrhs,X.size(1))
    X.copy(B);
    F77NAME(zgbtrs)(trans,dim,nsub,nsup,nrhs,BF->addr(),BF->nt,ipiv,
      X.addr(),dim,info);
    CHECK_TEST(info==0);
  } else { // zgbtrs
    int nrhs=B.size(0);
    CHECK_SAME(dim,B.size(1))
    CHECK_SAME(dim,X.size(1))
    CHECK_SAME(nrhs,X.size(0))
    X.copy(B);
//  int kd=nsub+nsup+1;
    if (trans!='N' && trans!='n') {
    // A = Q^T L U and B = X A^H = X U^H L^H Q ==> Q^T L U X^H = B^H
    //              or B = X A^T = X U^T L^T Q ==> Q^T L U X^T = B^T
      if (trans=='C' || trans=='c') {
        for (int j=0;j<dim;j++) {
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<nrhs;i++,Xij++) *Xij=conj(*Xij);
        }
      }
      if (nsub>0) { // solve Q^T L Y^T = B^T: zgbtrs
        for (int j=0;j<dim-1;j++) {
          int lm=min(nsub,dim-j-1);
          int jp=ipiv[j]-1;
          if (jp!=j) F77NAME(zswap)(nrhs,X.addr(0,jp),1,X.addr(0,j),1);
          F77NAME(zgeru)(nrhs,lm,double_mone_,BF->addr(j+1,j),1,
            X.addr(0,j),1,X.addr(0,j+1),1);
        }
      }
      for (int i=0;i<nrhs;i++) { // solve U X^T = Y^T: ztbsv loop 20
        for (int j=dim-1;j>=0;j--) {
          if (abs(X(i,j))>double_zero_) {
            X(i,j)/=(*BF)(j,j);
            int ii=max(0,j-nsub-nsup);
            F77NAME(zaxpy)(j-ii,-X(i,j),BF->addr(ii,j),1,
              X.addr(i,ii),nrhs);
          }
        }
      }
      if (trans=='C' || trans=='c') {
        for (int j=0;j<dim;j++) {
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<nrhs;i++,Xij++) *Xij=conj(*Xij);
        }
      }
    } else {
    // A = Q^T L U and B = X A = X Q^T L U ==> U^T L^T Q X^T = B^T
      for (int i=0;i<nrhs;i++) { // solve U^T Y^T=B^T: ztbsv loop 100
        for (int j=0;j<dim;j++) {
          int ii=max(0,j-nsub-nsup);
          X(i,j)=(X(i,j)
            -F77NAME(zdotu)(j-i,BF->addr(i,j),1,X.addr(i,ii),1))
            /(*BF)(j,j);
        }
      }
      if (nsub>0) { // solve L^T Q X^T = Y^T: zgbtrs loop 40
        for (int j=dim-2;j>=0;j--) {
          int lm=min(nsub,dim-1-j);
          F77NAME(zgemv)('T',nrhs,lm,double_mone_,X.addr(0,j+1),nrhs,
            BF->addr(j+1,j),1,double_one_,X.addr(0,j),nrhs);
          int jp=ipiv[j]-1;
          if (jp!=j) {
            F77NAME(zswap)(nrhs,X.addr(0,jp),1,X.addr(0,j),1);
          }
        }
      }
    }
  }

  delete [] ipiv;
  delete BF; BF=0;
}

template class BandMatrix<double,complex<double> >;
//template void testBandMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const complex<double>
  SymmetricBandMatrix<double,complex<double> >::outofbounds_
  = complex_double_zero_;
template<> complex<double>
  SymmetricBandMatrix<double,complex<double> >::safety_ =
  complex_double_zero_;

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::makeMatrix() const {
  SquareMatrix<double,complex<double> > *M=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(dim-j,nsub+1),addr(j,j),1,M->addr(j,j),1);
    if (j+1<dim) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Mji=M->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Mji+=dim) {
        *Mji=conj(*Aij);
      }
    }
  }
  return M;
}

template<> void SymmetricBandMatrix<double,complex<double> >::fillWith(
complex<double> d) {
  complex<double> *colj=addr();
  for (int j=0;j<dim;j++,colj+=nt) {
    colj[0]=d.real();
    for (int i=1;i<min(dim-j,nt);i++) colj[i]=d;
  }
}

template<> complex<double>
SymmetricBandMatrix<double,complex<double> >::operator()(int i,int j)
const {
  if (i-j>=0 && i-j<=nsub) return *AB->addr(i-j,j);
  if (j-i>0 && j-i<=nsub) return conj(*AB->addr(j-i,i));
  return outofbounds_;
}

template<>
SymmetricBandMatrix<double,complex<double> >::SymmetricBandMatrix(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) : nsub(1),
nt(2) {
  dim=T.size(0);
  AB=OPERATOR_NEW Matrix<double,complex<double> >(2,dim);
  complex<double> *dii=addr(0,0);
  double *Tii=T.diagonalAddr(0);
  for (int i=0;i<dim;i++,dii+=2,Tii++) *dii=*Tii;
  F77NAME(zcopy)(dim-1,T.lowerDiagonalAddr(0),1,addr(1,0),2);
  (*AB)(1,dim-1)=complex_double_undefined_;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::makeBandMatrix() const {
  BandMatrix<double,complex<double> > *B=
    OPERATOR_NEW BandMatrix<double,complex<double> >(dim,nsub,nsub);
  int stride=B->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,B->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Bji=B->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Bji+=stride) {
        *Bji=conj(*Aij);
      }
    }
  }
  return B;
}

template<> SymmetricMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::makeSymmetricMatrix()
const {
  SymmetricMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >&
SymmetricBandMatrix<double,complex<double> >::operator+=(
const SymmetricBandMatrix<double,complex<double> > &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  for (int j=0;j<dim;j++) {
    F77NAME(zaxpy)(min(nt,dim-j),complex_double_one_,B.addr(j,j),1,
      addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<double,complex<double> >&
SymmetricBandMatrix<double,complex<double> >::operator-=(
const SymmetricBandMatrix<double,complex<double> > &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  for (int j=0;j<dim;j++) {
    F77NAME(zaxpy)(min(nt,dim-j),complex_double_mone_,B.addr(j,j),1,
      addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<double,complex<double> >&
SymmetricBandMatrix<double,complex<double> >::operator*=(
complex<double> d) {
  for (int j=0;j<dim;j++) {
    F77NAME(zscal)(min(nt,dim-j),d,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<double,complex<double> >&
SymmetricBandMatrix<double,complex<double> >::operator/=(
complex<double> d) {
  CHECK_TEST(abs(d)>double_zero_);
  complex<double> dinv=complex_double_one_/d;
  for (int j=0;j<dim;j++) {
    F77NAME(zscal)(min(nt,dim-j),dinv,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const SymmetricBandMatrix<double,complex<double> > &B) const {
  CHECK_SAME(dim,B.dim);
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(dim,max(nsub,B.nsub),
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(zaxpy)(min(B.nt,dim-j),complex_double_one_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const BandMatrix<double,complex<double> > &B) const {
  CHECK_SAME(dim,B.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,B.subDiags()),
    max(nsub,B.supDiags()),complex_double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=stride) {
        *Sji=conj(*Aij);
      }
    }
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zaxpy)(min(dim-1,j+B.subDiags())-ibeg+1,complex_double_one_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const UpperHessenbergMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
        *Sji=conj(*Aij);
      }
    }
    F77NAME(zaxpy)(min(dim,j+2),complex_double_one_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const DiagonalMatrix<double,complex<double> > &D) const {
  CHECK_SAME(dim,D.size(0));
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(*this);
  F77NAME(zaxpy)(dim,complex_double_one_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(dim,max(nsub,1),
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  F77NAME(zaxpy)(dim-1,complex_double_one_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<dim;i++,Tii++,Sii+=S->nt) *Sii+=*Tii;
  return S;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const TridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,1),max(nsub,1),
    complex_double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=stride) {
        *Sji=conj(*Aij);
      }
    }
  }
  F77NAME(zaxpy)(dim,complex_double_one_,T.diagonalAddr(),1,
    S->addr(0,0),S->bands());
  F77NAME(zaxpy)(dim-1,complex_double_one_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S->bands());
  F77NAME(zaxpy)(dim-1,complex_double_one_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S->bands());
  return S;
}

template<> SymmetricMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const SymmetricMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SymmetricMatrix<double,complex<double> > *S=makeSymmetricMatrix();
  for (int j=0;j<dim;j++) {
    F77NAME(zaxpy)(dim-j,complex_double_one_,H.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<double> *Aij=addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      F77NAME(zaxpy)(j+1,complex_double_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<double> *Aij=addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)+=complex_double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<double> *Aij=addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      F77NAME(zaxpy)(dim-j,complex_double_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<double> *Aij=addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      (*S)(j,j)+=complex_double_one_;
      if (j+1<dim) {
        F77NAME(zaxpy)(dim-j-1,complex_double_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator+(
const Matrix<double,complex<double> > &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  S->copy(M);
  for (int j=0;j<dim;j++) {
    F77NAME(zaxpy)(min(nt,dim-j),complex_double_one_,addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
        *Sji+=conj(*Aij);
      }
    }
  }
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const SymmetricBandMatrix<double,complex<double> > &B) const {
  CHECK_SAME(dim,B.dim);
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(dim,max(nsub,B.nsub),
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(zaxpy)(min(B.nt,dim-j),complex_double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const BandMatrix<double,complex<double> > &B) const {
  CHECK_SAME(dim,B.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,B.subDiags()),
    max(nsub,B.supDiags()),complex_double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=stride) {
        *Sji=conj(*Aij);
      }
    }
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zaxpy)(min(dim-1,j+B.subDiags())-ibeg+1,complex_double_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >* operator-(
const BandMatrix<double,complex<double> > &B,
const SymmetricBandMatrix<double,complex<double> > &H) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,max(H.subDiags(),B.subDiags()),
    max(H.subDiags(),B.supDiags()),complex_double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(zcopy)(min(n-1,j+B.subDiags())-ibeg+1,B.addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(H.bands(),n-j),complex_double_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<n && H.subDiags()>0) {
      const complex<double> *Hij=H.addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(n-1,j+H.subDiags());i++,Hij++,Sji+=stride) {
        *Sji-=conj(*Hij);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const UpperHessenbergMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
        *Sji=conj(*Aij);
      }
    }
    F77NAME(zaxpy)(min(dim,j+2),complex_double_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const UpperHessenbergMatrix<double,complex<double> > &H,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(n,j+2),H.addr(0,j),1,S->addr(0,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<n && B.subDiags()>0) {
      const complex<double> *Bij=B.addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
        *Sji-=conj(*Bij);
      }
    }
  }
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const DiagonalMatrix<double,complex<double> > &D) const {
  CHECK_SAME(dim,D.size(0));
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(*this);
  F77NAME(zaxpy)(dim,complex_double_mone_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >* operator-(
const DiagonalMatrix<double,complex<double> > &D,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(n,B.subDiags(),
    complex_double_zero_);
  int nt=B.bands();
  F77NAME(zcopy)(n,D.addr(),1,S->addr(0,0),nt);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(nt,n-j),complex_double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(dim,max(nsub,1),
    complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  F77NAME(zaxpy)(dim-1,complex_double_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<dim;i++,Tii++,Sii+=S->nt) *Sii-=*Tii;
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >* operator-(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(n,max(B.subDiags(),1),
    complex_double_zero_);
  F77NAME(zcopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(1,0),S->bands());
  const double *Tii=T.diagonalAddr(0);
  complex<double> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,Tii++,Sii+=S->bands()) *Sii=*Tii;
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const TridiagonalMatrix<double,complex<double> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,max(nsub,1),max(nsub,1),
    complex_double_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=stride) {
        *Sji=conj(*Aij);
      }
    }
  }
  int S_nt=S->bands();
  F77NAME(zaxpy)(dim,complex_double_mone_,T.diagonalAddr(),1,
    S->addr(0,0),S_nt);
  F77NAME(zaxpy)(dim-1,complex_double_mone_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S_nt);
  F77NAME(zaxpy)(dim-1,complex_double_mone_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S_nt);
  return S;
}

template<> BandMatrix<double,complex<double> >* operator-(
const TridiagonalMatrix<double,complex<double> > &T,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<double,complex<double> > *S=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,max(B.subDiags(),1),
    max(B.subDiags(),1),complex_double_zero_);
  int S_nt=S->bands();
  F77NAME(zcopy)(n,T.diagonalAddr(),1,S->addr(0,0),S_nt);
  F77NAME(zcopy)(n-1,T.lowerDiagonalAddr(),1,S->addr(1,0),S_nt);
  F77NAME(zcopy)(n-1,T.upperDiagonalAddr(),1,S->addr(0,1),S_nt);
  int stride=S->bands()-1;
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<n && B.subDiags()>0) {
      const complex<double> *Bij=B.addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=stride) {
        *Sji-=conj(*Bij);
      }
    }
  }
  return S;
}

template<> SymmetricMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const SymmetricMatrix<double,complex<double> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SymmetricMatrix<double,complex<double> > *S=makeSymmetricMatrix();
  for (int j=0;j<dim;j++) {
    F77NAME(zaxpy)(dim-j,complex_double_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricMatrix<double,complex<double> >* operator-(
const SymmetricMatrix<double,complex<double> > &H,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SymmetricMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(n-j,H.addr(j,j),1,S->addr(j,j),1);
    F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<double> *Aij=addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      F77NAME(zaxpy)(j+1,complex_double_mone_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<double> *Aij=addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      if (j>0) {
        F77NAME(zaxpy)(j,complex_double_mone_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)-=complex_double_one_;
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(j+1,U.addr(0,j),1,S->addr(0,j),1);
    }
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,
        B.addr(j,j),1,S->addr(j,j),1);
      if (j+1<n && B.subDiags()>0) {
        const complex<double> *Bij=B.addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
          *Sji-=conj(*Bij);
        }
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) F77NAME(zcopy)(j,U.addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)=complex_double_one_;
    }
    for (int j=0;j<n;j++) {
      F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,
        B.addr(j,j),1,S->addr(j,j),1);
      if (j+1<n && B.subDiags()>0) {
        const complex<double> *Bij=B.addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
          *Sji-=conj(*Bij);
        }
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<double> *Aij=addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      F77NAME(zaxpy)(dim-j,complex_double_mone_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<double> *Aij=addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      (*S)(j,j)-=complex_double_one_;
      if (j+1<dim) {
        F77NAME(zaxpy)(dim-j-1,complex_double_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(zcopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,
        B.addr(j,j),1,S->addr(j,j),1);
      if (j+1<n && B.subDiags()>0) {
        const complex<double> *Bij=B.addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
          *Sji-=conj(*Bij);
        }
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=complex_double_one_;
      if (j+1<n) {
        F77NAME(zcopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,
        B.addr(j,j),1,S->addr(j,j),1);
      if (j+1<n && B.subDiags()>0) {
        const complex<double> *Bij=B.addr(j+1,j);
        complex<double> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
          *Sji-=conj(*Bij);
        }
      }
    }
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator-(
const Matrix<double,complex<double> > &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(zcopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<double> *Aij=addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
        *Sji=conj(*Aij);
      }
    }
    F77NAME(zaxpy)(dim,complex_double_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<double,complex<double> >* operator-(
const Matrix<double,complex<double> > &M,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<double,complex<double> > *S=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<n && B.subDiags()>0) {
      const complex<double> *Bij=B.addr(j+1,j);
      complex<double> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
        *Sji-=conj(*Bij);
      }
    }
  }
  return S;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const SymmetricBandMatrix<double,complex<double> > &B) const {
// compute by bordering: note that
// [ sigma s^H ] [ tau t^H ] = [ sigma tau + s^H t , sigma t^H + s^H T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^H +  S  T ]
  CHECK_SAME(dim,B.dim);
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,nsub+B.nsub,nsub+B.nsub,
    complex_double_zero_);
  int nb=min(nsub,B.nsub);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*B(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(dim-1,k+B.nsub);kk++) { // S t
        int ibeg=max(k+1,kk-B.nsub);
        if (kk>ibeg) {
          complex<double> Bkkk=B(kk,k);
          const complex<double> *Akki=addr(kk,ibeg);
          complex<double> *Pik=P->addr(ibeg,k);
          for (int i=ibeg;i<kk;i++,Akki+=nsub,Pik++) {
            *Pik+=conj(*Akki)*Bkkk;
          }
        }
        ibeg=max(kk,ibeg);
        F77NAME(zaxpy)(min(kk+nsub,dim-1)-ibeg+1,B(kk,k),
          addr(ibeg,kk),1,P->addr(ibeg,k),1);
      }
      for (int j=k+1;j<=min(k+nsub+B.nsub,dim-1);j++) { // s^H T
        int kbeg=max(k+1,j-B.nsub);
        int kend=min(j-1,k+nsub);
        if (kbeg<=kend) {
          (*P)(k,j)=conj(F77NAME(zdotu)(kend-kbeg+1,addr(kbeg,k),1,
            B.addr(j,kbeg),B.nt-1));
        }
        kbeg=max(kbeg,j);
        kend=min(dim-1,min(k+nsub,j+B.nsub));
        if (kend>=kbeg) {
          (*P)(k,j)+=
            F77NAME(zdotc)(kend-kbeg+1,addr(kbeg,k),1,B.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^H t
        F77NAME(zdotc)(min(nb,dim-k-1),addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(zgerc)(min(dim-k,nt),min(dim-k,B.nt),complex_double_one_,
      addr(k,k),1,B.addr(k,k),1,P->addr(k,k),P->bands()-1);
  }
  return P;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const BandMatrix<double,complex<double> > &B) const {
// compute by bordering: note that
// [ sigma s^H ] [ beta c^T ] = [ sigma beta + s^H b , sigma c^T + s^H B ]
// [   s    S  ] [  b    B  ] = [   s   beta +  S  b ,   s   c^T +  S  B ]
  CHECK_SAME(dim,B.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,nsub+B.subDiags(),
    nsub+B.supDiags(),complex_double_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*B(dim-1,dim-1);
  int nb=min(nsub,B.subDiags());
  for (int k=dim-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(dim-1,k+B.subDiags());kk++) { // S b
        int ibeg=max(k+1,kk-B.supDiags());
        if (kk>ibeg) {
          complex<double> Bkkk=B(kk,k);
          const complex<double> *Akki=addr(kk,ibeg);
          complex<double> *Pik=P->addr(ibeg,k);
          for (int i=ibeg;i<kk;i++,Akki+=nsub,Pik++) {
            *Pik+=conj(*Akki)*Bkkk;
          }
        }
        ibeg=max(kk,ibeg);
        F77NAME(zaxpy)(min(kk+nsub,dim-1)-ibeg+1,B(kk,k),
          addr(ibeg,kk),1,P->addr(ibeg,k),1);
      }
      for (int j=k+1;j<=min(k+P->supDiags(),dim-1);j++) { // s^H B
        int kbeg=max(k+1,j-B.supDiags());
        int kend=min(dim-1,min(k+nsub,j+B.subDiags()));
        if (kbeg<=kend) {
          (*P)(k,j)=F77NAME(zdotc)(kend-kbeg+1,addr(kbeg,k),1,
            B.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^H b
        F77NAME(zdotc)(min(nb,dim-k-1),addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(zgeru)(min(dim-k,nt),min(dim-k,B.supDiags()+1),
      complex_double_one_,addr(k,k),1,B.addr(k,k),1,
      P->addr(k,k),P->bands()-1);
  }
  return P;
}

template<> BandMatrix<double,complex<double> >* operator*(
const BandMatrix<double,complex<double> > &B,
const SymmetricBandMatrix<double,complex<double> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,B.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,S.subDiags()+B.subDiags(),
      S.subDiags()+B.supDiags(),complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-S.subDiags());k<=min(n-1,j+S.subDiags());k++) {
      int ibeg=max(0,k-B.supDiags());
      F77NAME(zaxpy)(min(n-1,k+B.subDiags())-ibeg+1,S(k,j),
        B.addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const UpperHessenbergMatrix<double,complex<double> > &H) const {
// compute by bordering: note that
// [ sigma s^H ] [     eta_11 h^T ]
// [   s    S  ] [ e_0 eta_21  H  ]
// = [ sigma eta_11 + s^H e_0 eta_21 , sigma h^T + s^H H ]
// = [   s   eta_11 +  S  e_0 eta_21 ,   s   h^T +  S  H ]
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*H(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nsub>0) {
      F77NAME(zaxpy)(min(dim-k-1,nsub+1),H(k+1,k),addr(k+1,k+1),1,
        P->addr(k+1,k),1); // S  e_0 eta_21
      for (int j=k+1;j<dim;j++) { // s^H H
        (*P)(k,j)=
          F77NAME(zdotc)(min(dim-k-1,min(j-k+1,nsub)),addr(k+1,k),1,
            H.addr(k+1,j),1);
      }
      (*P)(k,k)=conj((*this)(k+1,k))*H(k+1,k); // s^H e_0 eta_21
    }
    F77NAME(zgeru)(min(dim-k,nsub+1),dim-k,complex_double_one_,
      addr(k,k),1,H.addr(k,k),dim,P->addr(k,k),dim);
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const UpperHessenbergMatrix<double,complex<double> > &H,
const SymmetricBandMatrix<double,complex<double> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-S.subDiags());k<=min(n-1,j+S.subDiags());k++) {
      F77NAME(zaxpy)(min(n,k+2),S(k,j),H.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const DiagonalMatrix<double,complex<double> > &D) const {
  BandMatrix<double,complex<double> > *P=makeBandMatrix();
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsub);
    F77NAME(zscal)(min(dim-1,j+nsub)-ibeg+1,D[j],P->addr(ibeg,j),1);
  }
  return P;
}

template<> BandMatrix<double,complex<double> >* operator*(
const DiagonalMatrix<double,complex<double> > &D,
const SymmetricBandMatrix<double,complex<double> > &S) {
  int n=S.size(0);
  BandMatrix<double,complex<double> > *P=S.makeBandMatrix();
  int stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    int jbeg=max(0,i-S.subDiags());
    F77NAME(zscal)(min(n-1,i+S.subDiags())-jbeg+1,D[i],
      P->addr(i,jbeg),stride);
  }
  return P;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
// compute by bordering: note that
// [ sigma s^H ] [     tau    , conj(lambda) e_0^T ]
// [   s    S  ] [ e_0 lambda ,                 T  ]
//   = [ sigma tau + s^H e_0 lambda , sigma conj(lambda) e_0^T + s^H T ]
//   = [   s   tau +  S  e_o lambda ,   s   conj(lambda) e_0^T +  S  T ]
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,nsub+1,nsub+1,
    complex_double_zero_);
  (*P)(dim-2,dim-2)=(*this)(dim-2,dim-2)+T.diagonalValue(dim-2);
  (*P)(dim-1,dim-2)=(*this)(dim-1,dim-1)*T.lowerDiagonalValue(dim-2);
  (*P)(dim-2,dim-1)=(*this)(dim-2,dim-2)*T.upperDiagonalValue(dim-2);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*T.diagonalValue(dim-1);
  if (nsub>0) {
    (*P)(dim-2,dim-2)+=(*this)(dim-2,dim-1)*T.lowerDiagonalValue(dim-2);
    (*P)(dim-1,dim-2)+=(*this)(dim-1,dim-2)*T.diagonalValue(dim-2);
    (*P)(dim-2,dim-1)+=(*this)(dim-2,dim-1)*T.diagonalValue(dim-1);
    (*P)(dim-1,dim-1)+=(*this)(dim-1,dim-2)*T.upperDiagonalValue(dim-2);
  }
  for (int k=dim-3;k>=0;k--) {
    if (nsub>0) {
      for (int j=k+1;j<=min(k+nsub+1,dim-1);j++) { // s^H T
        if (j-1>=k+1) (*P)(k,j)=(*this)(k,j-1)*T.upperDiagonalValue(j-1);
        if (j<=k+nsub) (*P)(k,j)+=(*this)(k,j)*T.diagonalValue(j);
        if (j+1<=k+nsub) {
          (*P)(k,j)+=(*this)(k,j+1)*T.lowerDiagonalValue(j);
        }
      }
      (*P)(k,k)=(*this)(k,k+1)*T.lowerDiagonalValue(k); // s^H e_0 lambda
    }
    F77NAME(zaxpy)(min(dim-k-1,nsub+1),T.lowerDiagonalValue(k),
      addr(k+1,k+1),1,P->addr(k+1,k),1); // S e_0 lambda
    F77NAME(zaxpy)(min(dim-k,nsub+1),T.diagonalValue(k),addr(k,k),1,
      P->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(zaxpy)(min(dim-k,nsub+1),T.upperDiagonalValue(k),addr(k,k),1,
      P->addr(k,k+1),1);
     // [ sigma ] conj(lambda)
     // [   s   ]
  }
  return P;
}

template<> BandMatrix<double,complex<double> >* operator*(
const SymmetricTridiagonalMatrix<double,complex<double> > &T,
const SymmetricBandMatrix<double,complex<double> > &S) {
// compute by bordering: note that
// [     tau    , conj(lambda) e_0^T ] [ sigma s^H ]
// [ e_0 lambda ,            T       ] [   s    S  ]
// = [ tau sigma + conj(lambda) e_0^T s, tau s^H  + conj(lambda) e_0^T S ]
// = [ e_0 lambda sigma +          T s , e_0 lambda s^H +            T S ]
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=S.subDiags();
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,nb+1,nb+1,complex_double_zero_);
  (*P)(n-2,n-2)=T.diagonalValue(n-2)*S(n-2,n-2);
  (*P)(n-1,n-2)=T.lowerDiagonalValue(n-2)*S(n-2,n-2);
  (*P)(n-2,n-1)=T.upperDiagonalValue(n-2)*S(n-1,n-1);
  (*P)(n-1,n-1)=T.diagonalValue(n-1)*S(n-1,n-1);
  if (nb>0) {
    (*P)(n-2,n-2)+=T.upperDiagonalValue(n-2)*S(n-1,n-2);
    (*P)(n-1,n-2)+=T.diagonalValue(n-1)*S(n-1,n-2);
    (*P)(n-2,n-1)+=T.diagonalValue(n-2)*S(n-2,n-1);
    (*P)(n-1,n-1)+=T.lowerDiagonalValue(n-2)*S(n-2,n-1);
  }
  int stride=P->bands()-1;
  for (int k=n-3;k>=0;k--) {
    if (nb>0) {
      for (int i=k+1;i<=min(n-1,k+nb+1);i++) { // T s
        if (i-1>=k+1) (*P)(i,k)=T.lowerDiagonalValue(i-1)*S(i-1,k); 
        if (i<=k+nb) (*P)(i,k)+=T.diagonalValue(i)*S(i,k); 
        if (i+1<=min(k+nb,n-1)) {
          (*P)(i,k)+=T.upperDiagonalValue(i)*S(i+1,k); 
        }
      }
      (*P)(k,k)=T.upperDiagonalValue(k)*S(k+1,k);
    }
    complex<double> Tkkp1=T.upperDiagonalValue(k);
    const complex<double> *Sjkp1=S.addr(k+1,k+1);
    complex<double> *Pkj=P->addr(k,k+1);
    for (int j=k+1;j<=min(k+1+nb,n-1);j++,Sjkp1++,Pkj+=stride) {
      *Pkj+=Tkkp1*conj(*Sjkp1);
    }
      // conj(lambda) e_0^T S = conj(lambda) ( S e_0 )^H
    double Tkk=T.diagonalValue(k);
    const complex<double> *Sjk=S.addr(k,k);
    Pkj=P->addr(k,k);
    for (int j=k;j<=min(n-1,k+nb);j++,Sjk++,Pkj+=stride) {
      *Pkj+=Tkk*conj(*Sjk);
    }
      // tau [ sigma , s^H ]
    complex<double> Tkp1k=T.lowerDiagonalValue(k);
    Sjk=S.addr(k,k);
    complex<double> *Pkp1j=P->addr(k+1,k);
    for (int j=k;j<=min(n-1,k+nb);j++,Sjk++,Pkp1j+=stride) {
      *Pkp1j+=Tkp1k*conj(*Sjk);
    }
      // e_0 lambda [ sigma , s^H ]
  }
  return P;
}

template<> BandMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const TridiagonalMatrix<double,complex<double> > &T) const {
// compute by bordering: note that
// [ sigma s^H ] [    tau     upsilon e_0^T ]
// [   s    S  ] [ e_0 lambda       T       ]
//   = [ sigma tau + s^H e_0 lambda , sigma upsilon e_0^T + s^H T ]
//   = [     s tau +   S e_0 lambda ,     s upsilon e_0^T +   S T ]
  CHECK_SAME(dim,T.size(0));
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(dim,nsub+1,nsub+1,
    complex_double_zero_);
  (*P)(dim-2,dim-2)=(*this)(dim-2,dim-2)*T(dim-2,dim-2);
  (*P)(dim-1,dim-2)=(*this)(dim-1,dim-1)*T(dim-1,dim-2);
  (*P)(dim-2,dim-1)=(*this)(dim-2,dim-2)*T(dim-2,dim-1);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*T(dim-1,dim-1);
  if (nsub>0) {
    (*P)(dim-2,dim-2)+=(*this)(dim-2,dim-1)*T(dim-1,dim-2);
    (*P)(dim-1,dim-2)+=(*this)(dim-1,dim-2)*T(dim-2,dim-2);
    (*P)(dim-2,dim-1)+=(*this)(dim-2,dim-1)*T(dim-1,dim-1);
    (*P)(dim-1,dim-1)+=(*this)(dim-1,dim-2)*T(dim-2,dim-1);
  }
  for (int k=dim-3;k>=0;k--) {
    if (nsub>0) {
      for (int j=k+1;j<=min(k+nsub+1,dim-1);j++) { // s^H T
        if (j-1>=k+1) (*P)(k,j)=(*this)(k,j-1)*T(j-1,j);
        if (j<=k+nsub) (*P)(k,j)+=(*this)(k,j)*T(j,j);
        if (j+1<=k+nsub) (*P)(k,j)+=(*this)(k,j+1)*T(j+1,j);
      }
      (*P)(k,k)=(*this)(k,k+1)*T(k+1,k); // s^H e_0 lambda
    }
    F77NAME(zaxpy)(min(dim-k-1,nsub+1),T(k+1,k),addr(k+1,k+1),1,
      P->addr(k+1,k),1); // S e_0 lambda
    F77NAME(zaxpy)(min(dim-k,nsub+1),T(k,k),addr(k,k),1,P->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(zaxpy)(min(dim-k,nsub+1),T(k,k+1),addr(k,k),1,
      P->addr(k,k+1),1);
      // [ sigma ] lambda
      // [   s   ]
  }
  return P;
}

template<> BandMatrix<double,complex<double> >* operator*(
const TridiagonalMatrix<double,complex<double> > &T,
const SymmetricBandMatrix<double,complex<double> > &S) {
// compute by bordering: note that
// [     tau    , upsilon e_0^T ] [ sigma s^H ]
// [ e_0 lambda ,       T       ] [   s    S  ]
//   = [ tau sigma + upsilon e_0^T s , tau s^H        + upsilon e_0^T S ]
//   = [ e_0 lambda sigma +      T s , e_0 lambda s^H +             T S ]
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=S.subDiags();
  BandMatrix<double,complex<double> > *P=OPERATOR_NEW
    BandMatrix<double,complex<double> >(n,nb+1,nb+1,complex_double_zero_);
  (*P)(n-2,n-2)=T(n-2,n-2)*S(n-2,n-2);
  (*P)(n-1,n-2)=T(n-1,n-2)*S(n-2,n-2);
  (*P)(n-2,n-1)=T(n-2,n-1)*S(n-1,n-1);
  (*P)(n-1,n-1)=T(n-1,n-1)*S(n-1,n-1);
  if (nb>0) {
    (*P)(n-2,n-2)+=T(n-2,n-1)*S(n-1,n-2);
    (*P)(n-1,n-2)+=T(n-1,n-1)*S(n-1,n-2);
    (*P)(n-2,n-1)+=T(n-2,n-2)*S(n-2,n-1);
    (*P)(n-1,n-1)+=T(n-1,n-2)*S(n-2,n-1);
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
    complex<double> Tkkp1=T(k,k+1);
    const complex<double> *Sjkp1=S.addr(k+1,k+1);
    complex<double> *Pkj=P->addr(k,k+1);
    for (int j=k+1;j<=min(n-1,k+nb+1);j++,Sjkp1++,Pkj+=stride) {
      *Pkj+=Tkkp1*conj(*Sjkp1);
    }
      // upsilon e_0^T S = upsilon ( S e_0 )^T
    complex<double> Tkk=T(k,k);
    const complex<double> *Sjk=S.addr(k,k);
    Pkj=P->addr(k,k);
    for (int j=k;j<=min(n-1,k+nb);j++,Sjk++,Pkj+=stride) {
      *Pkj+=Tkk*conj(*Sjk);
    }
      // tau [ sigma , s^H ]
    complex<double> Tkp1k=T(k+1,k);
    Sjk=S.addr(k,k);
    complex<double> *Pkp1j=P->addr(k+1,k);
    for (int j=k;j<=min(n-1,k+nb);j++,Sjk++,Pkp1j+=stride) {
      *Pkp1j+=Tkp1k*conj(*Sjk);
    }
      // e_0 lambda [ sigma , s^H ]
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const SymmetricMatrix<double,complex<double> > &S) const {
// compute by bordering: note that
// [ sigma s^H ] [ tau t^H ] = [ sigma tau + s^H t , sigma t^H + s^H T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^H +  S  T ]
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(dim,complex_double_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*S(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nsub>0) {
      F77NAME(zhbmv)('L',dim-k-1,nsub,complex_double_one_,
        addr(k+1,k+1),nt,S.addr(k+1,k),1,complex_double_zero_,
        P->addr(k+1,k),1); // S t
      for (int j=k+1;j<dim;j++) { // s^H T
        int kend=min(j-1,k+nsub);
        if (k<kend) {
          (*P)(k,j)=
            conj(F77NAME(zdotu)(kend-k,addr(k+1,k),1,S.addr(j,k+1),dim));
        }
        int kbeg=max(k+1,j);
        kend=min(dim-1,k+nsub);
        if (kbeg<=kend) {
          (*P)(k,j)+=
            F77NAME(zdotc)(kend-kbeg+1,addr(kbeg,k),1,S.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^H t
        F77NAME(zdotc)(min(nsub,dim-k-1),addr(k+1,k),1,S.addr(k+1,k),1);
    }
    F77NAME(zgerc)(min(dim-k,nt),dim-k,complex_double_one_,addr(k,k),1,
      S.addr(k,k),1,P->addr(k,k),dim);
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SymmetricMatrix<double,complex<double> > &S,
const SymmetricBandMatrix<double,complex<double> > &B) {
// compute by bordering: note that
// [ sigma s^H ] [ tau t^H ] = [ sigma tau + s^H t , sigma t^H + s^H T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^H +  S  T ]
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  (*P)(n-1,n-1)=S(n-1,n-1)*B(n-1,n-1);
  int nb=B.subDiags();
  for (int k=n-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(n-1,k+nb);kk++) { // S t 
        int ibeg=k+1;
        complex<double> Bkkk=B(kk,k);
        if (kk>ibeg) {
          const complex<double> *Skki=S.addr(kk,ibeg);
          complex<double> *Pik=P->addr(ibeg,k);
          for (int i=ibeg;i<kk;i++,Skki+=n,Pik++) {
            *Pik+=conj(*Skki)*Bkkk;
          }
        }
        F77NAME(zaxpy)(n-kk,Bkkk,S.addr(kk,kk),1,P->addr(kk,k),1);
      }
      for (int j=k+1;j<n;j++) { // s^H T
        int kbeg=max(k+1,j-nb);
        if (kbeg<j) {
          (*P)(k,j)=conj(F77NAME(zdotu)(j-kbeg,S.addr(kbeg,k),1,
            B.addr(j,kbeg),B.bands()-1));
        }
        kbeg=max(kbeg,j);
        (*P)(k,j)+=F77NAME(zdotc)(min(n-1,j+nb)-kbeg+1,
          S.addr(kbeg,k),1,B.addr(kbeg,j),1);
      }
      (*P)(k,k)= // s^H t
        F77NAME(zdotc)(min(nb,n-k-1),S.addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(zgerc)(n-k,min(n-k,B.bands()),complex_double_one_,
      S.addr(k,k),1,B.addr(k,k),1,P->addr(k,k),n);
  }
  return P;
}

template<> Matrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U) const {
  CHECK_SAME(dim,U.size(0));
  int n=U.size(1);
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(dim,n,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,dim-1);k++) {
        int ibeg=max(0,k-nsub);
        complex<double> Ukj=U(k,j);
        if (ibeg<k) {
          const complex<double> *Aki=addr(k,ibeg);
          complex<double> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<k;i++,Aki+=nsub,Pij++) {
            *Pij+=conj(*Aki)*Ukj;
          }
        }
        ibeg=max(ibeg,k);
        F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,Ukj,addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,dim);k++) {
        int ibeg=max(0,k-nsub);
        complex<double> Ukj=U(k,j);
        if (ibeg<k) {
          const complex<double> *Aki=addr(k,ibeg);
          complex<double> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<k;i++,Aki+=nsub,Pij++) {
            *Pij+=conj(*Aki)*Ukj;
          }
        }
        ibeg=max(ibeg,k);
        F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,Ukj,addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
      if (j<dim) {
        int ibeg=max(0,j-nsub);
        if (ibeg<j) {
          const complex<double> *Aji=addr(j,ibeg);
          complex<double> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<j;i++,Aji+=nsub,Pij++) {
            *Pij+=conj(*Aji);
          }
        }
        ibeg=max(ibeg,j);
        F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,complex_double_one_,
          addr(ibeg,j),1,P->addr(ibeg,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >* operator*(
const UpperTrapezoidalMatrix<double,complex<double> > &U,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<double,complex<double> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(zaxpy)(min(k+1,m),B(k,j),U.addr(0,k),1,P->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(zaxpy)(min(k,m),B(k,j),U.addr(0,k),1,P->addr(0,j),1);
        if (k<m) (*P)(k,j)+=B(k,j);
      }
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L) const {
  CHECK_SAME(dim,L.size(0));
  int n=L.size(1);
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(dim,n,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<dim;k++) {
        int ibeg=max(0,k-nsub);
        complex<double> Lkj=L(k,j);
        if (ibeg<k) {
          const complex<double> *Aki=addr(k,ibeg);
          complex<double> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<k;i++,Aki+=nsub,Pij++) {
            *Pij+=conj(*Aki)*Lkj;
          }
        }
        ibeg=max(ibeg,k);
        F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,Lkj,addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-nsub);
      if (ibeg<j) {
        const complex<double> *Aji=addr(j,ibeg);
        complex<double> *Pij=P->addr(ibeg,j);
        for (int i=ibeg;i<j;i++,Aji+=nsub,Pij++) {
          *Pij+=conj(*Aji);
        }
      }
      ibeg=max(ibeg,j);
      F77NAME(zaxpy)(min(dim-1,j+nsub)-ibeg+1,complex_double_one_,
        addr(ibeg,j),1,P->addr(ibeg,j),1);
      for (int k=j+1;k<dim;k++) {
        int ibeg=max(0,k-nsub);
        complex<double> Lkj=L(k,j);
        if (ibeg<k) {
          const complex<double> *Aki=addr(k,ibeg);
          complex<double> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<k;i++,Aki+=nsub,Pij++) {
            *Pij+=conj(*Aki)*Lkj;
          }
        }
        ibeg=max(ibeg,k);
        F77NAME(zaxpy)(min(dim-1,k+nsub)-ibeg+1,Lkj,addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >* operator*(
const LowerTrapezoidalMatrix<double,complex<double> > &L,
const SymmetricBandMatrix<double,complex<double> > &B) {
int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<double,complex<double> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(zaxpy)(m-k,B(k,j),L.addr(k,k),1,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        (*P)(k,j)+=B(k,j);
        if (k+1<m) {
          F77NAME(zaxpy)(m-k-1,B(k,j),L.addr(k+1,k),1,P->addr(k+1,j),1);
        }
      }
    }
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const SquareMatrix<double,complex<double> > &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<double,complex<double> > *P=
    OPERATOR_NEW SquareMatrix<double,complex<double> >(dim);
  for (int j=0;j<dim;j++) {
    F77NAME(zhbmv)('L',dim,nsub,complex_double_one_,addr(),nt,
      S.addr(0,j),1,complex_double_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> SquareMatrix<double,complex<double> >* operator*(
const SquareMatrix<double,complex<double> > &S,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<double,complex<double> > *P=OPERATOR_NEW
    SquareMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(zaxpy)(n,B(k,j),S.addr(0,j),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Matrix<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const Matrix<double,complex<double> > &M) const {
  CHECK_SAME(dim,M.size(0));
  int n=M.size(1);
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(dim,n);
  for (int j=0;j<dim;j++) {
    F77NAME(zhbmv)('L',dim,nsub,complex_double_one_,addr(),nt,
      M.addr(0,j),1,complex_double_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<double,complex<double> >* operator*(
const Matrix<double,complex<double> > &M,
const SymmetricBandMatrix<double,complex<double> > &B) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<double,complex<double> > *P=OPERATOR_NEW
    Matrix<double,complex<double> >(m,n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(zaxpy)(m,B(k,j),M.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Vector<double,complex<double> >*
SymmetricBandMatrix<double,complex<double> >::operator*(
const Vector<double,complex<double> > &v) const {
  CHECK_SAME(dim,v.size());
  Vector<double,complex<double> > *p=
    OPERATOR_NEW Vector<double,complex<double> >(dim);
  F77NAME(zhbmv)('L',dim,nsub,complex_double_one_,addr(),nt,v.addr(),1,
    complex_double_zero_,p->addr(),1);
  return p;
}

template<> void SymmetricBandMatrix<double,complex<double> >::sbmv(
complex<double> alpha,const Vector<double,complex<double> > &x,
complex<double> beta,Vector<double,complex<double> > &b) const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  F77NAME(zhbmv)('L',dim,nsub,alpha,addr(),nt,x.addr(),1,beta,
    b.addr(),1);
}

template<> void SymmetricBandMatrix<double,complex<double> >::sbmm(
complex<double> alpha,const Matrix<double,complex<double> > &X,
complex<double> beta,Matrix<double,complex<double> > &B,char side)
const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      F77NAME(zhbmv)('L',dim,nsub,alpha,addr(),nt,X.addr(0,j),1,beta,
        B.addr(0,j),1);
    }
  } else {
    CHECK_SAME(n,dim);
    if (abs(beta)==complex_double_zero_) B=complex_double_zero_;
    else B*=beta;
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-nsub);k<=min(dim-1,j+nsub);k++) {
        F77NAME(zaxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
          B.addr(0,j),1);
      }
    }
  }
}

template<> double
SymmetricBandMatrix<double,complex<double> >::normFrobenius() const {
  double *work=0;
  return F77NAME(zlanhb)('F','L',dim,nsub,addr(),nt,work);
}

template<> double
SymmetricBandMatrix<double,complex<double> >::normInfinity() const {
  double *work=OPERATOR_NEW_BRACKET(double,dim);;
  double val=F77NAME(zlanhb)('I','L',dim,nsub,addr(),nt,work);
  delete [] work; work=0;
  return val;
}

template<> double
SymmetricBandMatrix<double,complex<double> >::normMaxEntry() const {
  double *work=0;
  return F77NAME(zlanhb)('M','L',dim,nsub,addr(),nt,work);
}

template<> double
SymmetricBandMatrix<double,complex<double> >::normOne() const {
  double *work=OPERATOR_NEW_BRACKET(double,dim);;
  double val=F77NAME(zlanhb)('O','L',dim,nsub,addr(),nt,work);
  delete [] work; work=0;
  return val;
}

template<> Vector<double,double>* 
SymmetricBandMatrix<double,complex<double> >::eigenvalues(
OrthogonalMatrix<double,complex<double> > *&Q) const {
  if (Q!=0) CHECK_SAME(dim,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Vector<double,double> *lambda=OPERATOR_NEW Vector<double,double>(dim);
  SymmetricBandMatrix<double,complex<double> > *copy=
    OPERATOR_NEW SymmetricBandMatrix<double,complex<double> >(*this);
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,dim);
  double *rwork=OPERATOR_NEW_BRACKET(double,max(1,3*dim-2));
  int info;
  complex<double> *qa=( Q==0 ? 0 : Q->addr() );
  F77NAME(zhbev)(jobz,'L',dim,nsub,copy->addr(),copy->bands(),
    lambda->addr(),qa,dim,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete [] rwork; rwork=0;
  delete copy; copy=0;
  return lambda;
}
template class SymmetricBandMatrix<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveBandMatrix<double,complex<double> >*
SymmetricPositiveBandMatrix<double,complex<double> >::operator+(
const SymmetricPositiveBandMatrix<double,complex<double> > &B) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,B.size(0));
  SymmetricPositiveBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricPositiveBandMatrix<double,complex<double> >(n,
    max(ns,B.subDiags()),complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(zaxpy)(min(B.bands(),n-j),complex_double_one_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricPositiveBandMatrix<double,complex<double> >*
SymmetricPositiveBandMatrix<double,complex<double> >::operator+(
const SymmetricPositiveTridiagonalMatrix<double,complex<double> > &T)
const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricPositiveBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricPositiveBandMatrix<double,complex<double> >(n,max(ns,1),
    complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=T(j,j);
    if (j+1<n) (*S)(j+1,j)+=T(j+1,j);
  }
  return S;
}

template<> SymmetricBandMatrix<double,complex<double> >*
SymmetricPositiveBandMatrix<double,complex<double> >::operator+(
const SymmetricTridiagonalMatrix<double,complex<double> > &T) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricBandMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricBandMatrix<double,complex<double> >(n,max(ns,1),
    complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=T(j,j);
    if (j+1<n) (*S)(j+1,j)+=T(j+1,j);
  }
  return S;
}

template<> SymmetricPositiveMatrix<double,complex<double> >*
SymmetricPositiveBandMatrix<double,complex<double> >::operator+(
const SymmetricPositiveMatrix<double,complex<double> > &M) const {
  int n=size(0);
  int nb=bands();
  CHECK_SAME(n,M.size(0));
  SymmetricPositiveMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricPositiveMatrix<double,complex<double> >(n,
    complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(zaxpy)(n-j,complex_double_one_,M.addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricMatrix<double,complex<double> >*
SymmetricPositiveBandMatrix<double,complex<double> >::operator+(
const SymmetricMatrix<double,complex<double> > &T) const {
  int n=size(0);
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<double,complex<double> > *S=OPERATOR_NEW
    SymmetricMatrix<double,complex<double> >(n,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(zaxpy)(n-j,complex_double_one_,T.addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> double
SymmetricPositiveBandMatrix<double,complex<double> >::equilibrate(
Vector<double,double> &s,double &scond) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,s.size());
  double amax;
  int info;
  F77NAME(zpbequ)('L',n,ns,addr(),nb,s.addr(),scond,amax,info);
  CHECK_TEST(info==0);
  return amax;
}

template<> double
SymmetricPositiveBandMatrix<double,complex<double> >
::reciprocalConditionNumber() const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  SymmetricPositiveBandMatrix<double,complex<double> > *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<double,complex<double> >(
    *this);
  int info;
  F77NAME(zpbtrf)('L',n,ns,BF->addr(),nb,info);
  CHECK_TEST(info==0);

  double anorm=normOne();
  double rcond=numeric_limits<double>::infinity();
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  F77NAME(zpbcon)('L',n,ns,BF->addr(),nb,anorm,rcond,work,rwork,info);
  delete [] work; work=0;
  delete [] rwork; rwork=0;
  delete BF; BF=0;
  return rcond;
}

template<> void
SymmetricPositiveBandMatrix<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,char) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  SymmetricPositiveBandMatrix<double,complex<double> > *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<double,complex<double> >(
    *this);
  x.copy(b);
  int info;
  F77NAME(zpbsv)('L',n,ns,1,BF->addr(),nb,x.addr(),n,info);
  CHECK_TEST(info==0);
  delete BF; BF=0;
}

template<> void
SymmetricPositiveBandMatrix<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,char side,char) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  SymmetricPositiveBandMatrix<double,complex<double> > *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<double,complex<double> >(
    *this);
  int info;
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    CHECK_SAME(nrhs,X.size(1))
    X.copy(B);
    F77NAME(zpbsv)('L',n,ns,1,BF->addr(),nb,X.addr(),n,info);
    CHECK_TEST(info==0);
  } else {
    int nrhs=B.size(0);
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    CHECK_SAME(nrhs,X.size(0))
    X.copy(B);
    F77NAME(zpbtrf)('L',n,ns,BF->addr(),nb,info);
    CHECK_TEST(info==0);
    for (int j=0;j<nrhs;j++) { // dpbtrs loop 20:
      F77NAME(ztbsv)('L','N','N',n,ns,BF->addr(),nb,X.addr(0,j),nrhs);
      F77NAME(ztbsv)('L','T','N',n,ns,BF->addr(),nb,X.addr(0,j),nrhs);
    }
  }
  delete BF; BF=0;
}

template class SymmetricPositiveBandMatrix<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SpecializedMatrix.C"
template class CompanionMatrix<double,complex<double> >;
template class HadamardMatrix<double,complex<double> >;
template class HankelMatrix<double,complex<double> >;
template class HilbertMatrix<double,complex<double> >;
template class KahanMatrix<double,complex<double> >;
template class LauchliMatrix<double,complex<double> >;
template class PascalMatrix<double,complex<double> >;
template class RosserMatrix<double,complex<double> >;
template class ToeplitzMatrix<double,complex<double> >;
template class SymmetricToeplitzMatrix<double,complex<double> >;
template class VandermondeMatrix<double,complex<double> >;
template class WilkinsonMatrix<double,complex<double> >;
template void testSpecializedMatrix(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "GaussianFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> GaussianFactorization<double,complex<double>,
SquareMatrix<double,complex<double> > >::GaussianFactorization(
const SquareMatrix<double,complex<double> >& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),piv_op(po),equ_op(eo),
equed('N'),colcnd(numeric_limits<double>::infinity()),
rowcnd(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  LU=OPERATOR_NEW SquareMatrix<double,complex<double> >(n);
  LU->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    r=OPERATOR_NEW Vector<double,double>(n);
    c=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    F77NAME(zgeequ)(n,n,LU->addr(),n,r->addr(),c->addr(),rowcnd,colcnd,
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
      F77NAME(zlaqge)(n,n,LU->addr(),n,r->addr(),c->addr(),
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
	complex<double> alpha=complex_double_one_/(*LU)(k,k);
	int nrows=n-k-1;
	if (nrows>0) {
	  F77NAME(zscal)(nrows,alpha,LU->addr(k+1,k),1);
	  for (int j=k+1;j<n;j++) {
	    alpha=-(*LU)(k,j);
	    F77NAME(zaxpy)(nrows,alpha,LU->addr(k+1,k),1,
			   LU->addr(k+1,j),1);
	  }
	}
      }
      break;
    }
    case Factorization::PIVOT_ROWS_AND_COLUMNS: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      jpiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(zgetc2)(n,LU->addr(),n,ipiv,jpiv,info);
      break;
    }
    default: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(zgetrf)(n,n,LU->addr(),n,ipiv,info);
      break;
    }
  }
  double *work=0;
  if (info>0) { // see dgesvx
    rpvgrw=F77NAME(zlantr)('M','U','N',info,info,LU->addr(),n,work);
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
    rpvgrw=F77NAME(zlantr)('M','U','N',n,n,LU->addr(),n,work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
  }
}

template<> void GaussianFactorization<double,complex<double>,
SquareMatrix<double,complex<double> > >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,Factorization::TRANSPOSE_OPTION to) {
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
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
    if (ipiv!=0) F77NAME(zlaswp)(1,x.addr(),n,1,n-1,ipiv,1);
    F77NAME(ztrsm)('L','L','N','U',n,1,complex_double_one_,LU->addr(),n,
      x.addr(),n);
    F77NAME(ztrsm)('L','U','N','N',n,1,complex_double_one_,LU->addr(),n,
      x.addr(),n);
    if (jpiv!=0) F77NAME(zlaswp)(1,x.addr(),n,1,n-1,jpiv,-1);
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ci=c->addr();
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
  } else {
    // A^T X = B ==> U^T L^T Q^T R^{-1} X = P^T C B
    // A^H X = B ==> U^H L^H Q^T R^{-1} X = P^T C B
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ci=c->addr();
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
    if (jpiv!=0) F77NAME(zlaswp)(1,x.addr(),n,1,n-1,jpiv,1);
    if (to==Factorization::TRANSPOSE) {
      F77NAME(ztrsm)('L','U','T','N',n,1,complex_double_one_,LU->addr(),n,
        x.addr(),n);
      F77NAME(ztrsm)('L','L','T','U',n,1,complex_double_one_,LU->addr(),n,
        x.addr(),n);
    } else {
      F77NAME(ztrsm)('L','U','C','N',n,1,complex_double_one_,LU->addr(),n,
        x.addr(),n);
      F77NAME(ztrsm)('L','L','C','U',n,1,complex_double_one_,LU->addr(),n,
        x.addr(),n);
    }
    if (ipiv!=0) F77NAME(zlaswp)(1,x.addr(),n,1,n-1,ipiv,-1);
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ri=r->addr();
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
  }
}

template<> void GaussianFactorization<double,complex<double>,
SquareMatrix<double,complex<double> > >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Factorization::TRANSPOSE_OPTION to,
Factorization::SIDE_OPTION so) {
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
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
      if (ipiv!=0) F77NAME(zlaswp)(n,X.addr(),m,1,m-1,ipiv,1);
      F77NAME(ztrsm)('L','L','N','U',m,n,complex_double_one_,
                     LU->addr(),m,X.addr(),m);
      F77NAME(ztrsm)('L','U','N','N',m,n,complex_double_one_,
                     LU->addr(),m,X.addr(),m);
      if (jpiv!=0) F77NAME(zlaswp)(n,X.addr(),m,1,m-1,jpiv,-1);
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ci=c->addr();
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
    } else {
      // A^T X = B ==> U^T L^T Q^T R^{-1} X = P^T C B
      // A^H X = B ==> U^H L^H Q^T R^{-1} X = P^T C B
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ci=c->addr();
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
      if (jpiv!=0) F77NAME(zlaswp)(n,X.addr(),m,1,m-1,jpiv,1);
      if (to==Factorization::TRANSPOSE) {
        F77NAME(ztrsm)('L','U','T','N',m,n,complex_double_one_,
                       LU->addr(),m,X.addr(),m);
        F77NAME(ztrsm)('L','L','T','U',m,n,complex_double_one_,
                       LU->addr(),m,X.addr(),m);
      } else {
        F77NAME(ztrsm)('L','U','C','N',m,n,complex_double_one_,
                       LU->addr(),m,X.addr(),m);
        F77NAME(ztrsm)('L','L','C','U',m,n,complex_double_one_,
                       LU->addr(),m,X.addr(),m);
      }
      if (ipiv!=0) F77NAME(zlaswp)(n,X.addr(),m,1,m-1,ipiv,-1);
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ri=r->addr();
          complex<double> *Xij=X.addr(0,j);
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
        for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*c)[j],X.addr(0,j),1);
      }
      if (jpiv!=0) {
        for (int j=n-2;j>=0;j--) {
          cout << "jpiv[" << j << "] = " << jpiv[j] << endl;
          if (j!=jpiv[j]-1) {
            F77NAME(zswap)(m,X.addr(0,j),1,X.addr(0,jpiv[j]-1),1);
          }
        }
      }
      F77NAME(ztrsm)('R','U','N','N',m,n,complex_double_one_,
                     LU->addr(),n,X.addr(),m);
      F77NAME(ztrsm)('R','L','N','U',m,n,complex_double_one_,
                     LU->addr(),n,X.addr(),m);
      if (ipiv!=0) {
        for (int i=0;i<n-1;i++) {
          if (i!=ipiv[i]-1) {
            F77NAME(zswap)(m,X.addr(0,i),1,X.addr(0,ipiv[i]-1),1);
          }
        }
      }
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*r)[j],X.addr(0,j),1);
      }
    } else {
      // X A^T = B ==> X C^{-1} P U^T L^T = B R Q
      // X A^H = B ==> X C^{-1} P U^H L^H = B R Q
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*r)[j],X.addr(0,j),1);
      }
      if (ipiv!=0) {
        for (int i=n-2;i>=0;i--) {
          if (i!=ipiv[i]-1) {
            F77NAME(zswap)(m,X.addr(0,i),1,X.addr(0,ipiv[i]-1),1);
          }
        }
      }
      if (to==Factorization::TRANSPOSE) {
        F77NAME(ztrsm)('R','L','T','U',m,n,complex_double_one_,
                       LU->addr(),n,X.addr(),m);
        F77NAME(ztrsm)('R','U','T','N',m,n,complex_double_one_,
                       LU->addr(),n,X.addr(),m);
      } else {
        F77NAME(ztrsm)('R','L','C','U',m,n,complex_double_one_,
                       LU->addr(),n,X.addr(),m);
        F77NAME(ztrsm)('R','U','C','N',m,n,complex_double_one_,
                       LU->addr(),n,X.addr(),m);
      }
      if (jpiv!=0) {
        for (int j=0;j<n-1;j++) {
          if (j!=jpiv[j]-1) {
            F77NAME(zswap)(m,X.addr(0,j),1,X.addr(0,jpiv[j]-1),1);
          }
        }
      }
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*c)[j],X.addr(0,j),1);
      }
    }
  }
}

template<> double GaussianFactorization<double,complex<double>,
SquareMatrix<double,complex<double> > >::reciprocalConditionNumber(
Factorization::CONDITION_NUMBER_NORM cnn) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  double anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,2*n);
  int info=0;
  F77NAME(zgecon)(norm,n,LU->addr(),n,anorm,rcond,work,rwork,info);
  CHECK_SAME(info,0)
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void GaussianFactorization<double,complex<double>,
SquareMatrix<double,complex<double> > >::improve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double &berr,double &ferr,
Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(n+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,complex<double> > residual(n);
  Vector<double,double> work(n);
  Vector<double,complex<double> > v(n);
  char trans='N';
  if (to==Factorization::TRANSPOSE) trans='T';
  else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(zgemv)(trans,n,n,complex_double_mone_,A_original->addr(),n,
      x.addr(),1,complex_double_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(zgeamv)(trans,n,n,double_one_,A_original->addr(),n,
      x.addr(),1,double_one_,work.addr(),1);

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
    F77NAME(zaxpy)(n,complex_double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  complex<double> *residuali=residual.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(zlacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual,to);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual,to);
    }
  }
  int i=F77NAME(izamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void GaussianFactorization<double,complex<double>,
SquareMatrix<double,complex<double> > >::improve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Vector<double,double> &berr,
Vector<double,double> &ferr,Factorization::TRANSPOSE_OPTION to,
Factorization::SIDE_OPTION so) {
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
  Vector<double,complex<double> > x(k);
  Vector<double,complex<double> > rhs(k);
  Vector<double,complex<double> > residual(k);
  Vector<double,double> work(k);
  Vector<double,complex<double> > v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  char trans='N';
  Factorization::TRANSPOSE_OPTION lto=Factorization::NO_TRANSPOSE;
  if (so==Factorization::LEFT_SIDE) {
    if (to==Factorization::TRANSPOSE) trans='T';
    else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
    lto=to;
  } else {
    if (to==Factorization::NO_TRANSPOSE) {
      trans='T';
      lto=Factorization::TRANSPOSE;
    } else {
      trans='N';
      lto=Factorization::NO_TRANSPOSE;
    }
  }
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(zcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(zcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        const complex<double> *Xji=X.addr(j,0);
        complex<double> *xi=x.addr();
        const complex<double> *Bji=B.addr(j,0);
        complex<double> *rhsi=rhs.addr();
        for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
          *xi=conj(*Xji);
          *rhsi=conj(*Bji);
        }
      } else {
        F77NAME(zcopy)(k,X.addr(j,0),m,x.addr(),1);
        F77NAME(zcopy)(k,B.addr(j,0),m,rhs.addr(),1);
      }
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(zgemv)(trans,k,k,complex_double_mone_,A_original->addr(),k,
        x.addr(),1,complex_double_one_,residual.addr(),double_one_);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77NAME(zgeamv)(trans,k,k,double_one_,A_original->addr(),k,
        x.addr(),1,double_one_,work.addr(),1);

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
      F77NAME(zaxpy)(k,complex_double_one_,residual.addr(),1,x.addr(),1);
      lstres=s;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(zcopy)(k,x.addr(),1,X.addr(0,j),1);
      } else F77NAME(zcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<double> *Xji=X.addr(j,0);
        const complex<double> *xi=x.addr();
        for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
      } else {
        F77NAME(zcopy)(k,x.addr(),1,X.addr(j,0),m);
      }
    }

    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(zlacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual,lto);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual,lto);
      }
    }
    int i=F77NAME(izamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr[j]/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> GaussianFactorization<double,complex<double>,
TridiagonalMatrix<double,complex<double> > >::GaussianFactorization(
const TridiagonalMatrix<double,complex<double> >& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),
piv_op(Factorization::PIVOT_ROWS),equ_op(Factorization::NO_EQUILIBRATION),
equed('N'),colcnd(numeric_limits<double>::infinity()),
rowcnd(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  anormi=A.normInfinity();
  anormm=A.normMaxEntry();
  anormo=A.normOne();
  rpvgrw=double_one_; // not computed
}

template<> void GaussianFactorization<double,complex<double>,
TridiagonalMatrix<double,complex<double> > >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,Factorization::TRANSPOSE_OPTION to) {
//constructor did not factor
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  LU=OPERATOR_NEW TridiagonalMatrix<double,complex<double> >(n);
  LU->copy(*A_original);
  if (&x!=&b) x.copy(b);
  int info;
  if (to==Factorization::NO_TRANSPOSE) {
    if (piv_op==Factorization::PIVOT_ROWS) {
      F77NAME(zgtsv)(n,1,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),x.addr(),n,info);
    } else {
      F77NAME(zgtsvnp)(n,1,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),x.addr(),n,info);
    }
  } else {
    if (to==Factorization::CONJUGATE_TRANSPOSE) {
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,xi++) *xi=conj(*xi);
    }
    if (piv_op==Factorization::PIVOT_ROWS) {
      F77NAME(zgtsv)(n,1,LU->upperDiagonalAddr(),LU->diagonalAddr(),
        LU->lowerDiagonalAddr(),x.addr(),n,info);
    } else {
      F77NAME(zgtsvnp)(n,1,LU->upperDiagonalAddr(),LU->diagonalAddr(),
        LU->lowerDiagonalAddr(),x.addr(),n,info);
    }
    if (to==Factorization::CONJUGATE_TRANSPOSE) {
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,xi++) *xi=conj(*xi);
    }
  }
  CHECK_TEST(info==0);
  delete LU; LU=0;
}

template<> void GaussianFactorization<double,complex<double>,
TridiagonalMatrix<double,complex<double> > >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Factorization::TRANSPOSE_OPTION to,
Factorization::SIDE_OPTION so) {
//constructor did not factor
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  LU=OPERATOR_NEW
    TridiagonalMatrix<double,complex<double> >(A_original->size(0));
  LU->copy(*A_original);
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) { // A X = B
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(zgtsv)(m,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(zgtsvnp)(m,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      }
    } else { // A^T X = B, or A^H X = B <==> A^T conj(X) = conj(B)
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<double> *Xij=X.addr();
        for (int ij=0;ij<m*n;ij++,Xij++) *Xij=conj(*Xij);
      }
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(zgtsv)(m,n,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(zgtsvnp)(m,n,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      }
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<double> *Xij=X.addr();
        for (int ij=0;ij<m*n;ij++,Xij++) *Xij=conj(*Xij);
      }
    }
  } else {
    CHECK_SAME(n,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) { // X A = B ==> A^T X^T = B^T
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(zgtsvr)(n,m,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(zgtsvrnp)(n,m,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      }
    } else {
      // X A^T = B ==> A X^T = B^T
      // X A^H = B ==> A X^H = B^H
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<double> *Xij=X.addr();
        for (int ij=0;ij<m*n;ij++,Xij++) *Xij=conj(*Xij);
      }
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(zgtsvr)(n,m,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(zgtsvrnp)(n,m,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      }
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<double> *Xij=X.addr();
        for (int ij=0;ij<m*n;ij++,Xij++) *Xij=conj(*Xij);
      }
    }
  }
  CHECK_TEST(info==0);
  delete LU; LU=0;
}

template<> double GaussianFactorization<double,complex<double>,
TridiagonalMatrix<double,complex<double> > >::reciprocalConditionNumber(
Factorization::CONDITION_NUMBER_NORM cnn) {
  LU=OPERATOR_NEW
    TridiagonalMatrix<double,complex<double> >(A_original->size(0));
  LU->copy(*A_original);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  double anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  int info=0;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<double,complex<double> > *u2=
      OPERATOR_NEW Vector<double,complex<double> >(n-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(zgttrf)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    F77NAME(zgtcon)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,info);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(zgttrfnp)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    F77NAME(zgtconnp)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,info);
  }
  CHECK_SAME(info,0)
  delete LU; LU=0;
  delete [] work; work=0;
  return rcond;
}

template<> void GaussianFactorization<double,complex<double>,
TridiagonalMatrix<double,complex<double> > >::improve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double &berr,double &ferr,
Factorization::TRANSPOSE_OPTION to) {
  LU=OPERATOR_NEW
    TridiagonalMatrix<double,complex<double> >(A_original->size(0));
  LU->copy(*A_original);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  char norm=(to==Factorization::NO_TRANSPOSE ? 'O' : 'I');
  double anorm=(to==Factorization::NO_TRANSPOSE ? anormo : anormi);
  char trans='N';
  if (to==Factorization::TRANSPOSE) trans='T';
  else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info;
  double rcond;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<double,complex<double> > *u2=
      OPERATOR_NEW Vector<double,complex<double> >(n-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(zgttrf)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    CHECK_TEST(info==0);
    F77NAME(zgtcon)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,info);
    CHECK_TEST(info==0);
    F77NAME(zgtrfs)(trans,n,1,A_original->lowerDiagonalAddr(),
      A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
      LU->lowerDiagonalAddr(),LU->diagonalAddr(),LU->upperDiagonalAddr(),
      u2->addr(),ipiv,b.addr(),n,x.addr(),n,&ferr,&berr,work,rwork,info);
    CHECK_TEST(info==0);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(zgttrfnp)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    CHECK_TEST(info==0);
    F77NAME(zgtconnp)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,info);
    CHECK_TEST(info==0);
    F77NAME(zgtrfsnp)(trans,n,1,A_original->lowerDiagonalAddr(),
      A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
      LU->lowerDiagonalAddr(),LU->diagonalAddr(),LU->upperDiagonalAddr(),
      b.addr(),n,x.addr(),n,&ferr,&berr,work,rwork,info);
    CHECK_TEST(info==0);
  }
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  delete LU; LU=0;
}

template<> void GaussianFactorization<double,complex<double>,
TridiagonalMatrix<double,complex<double> > >::improve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Vector<double,double> &berr,
Vector<double,double> &ferr,Factorization::TRANSPOSE_OPTION to,
Factorization::SIDE_OPTION so) {
  LU=OPERATOR_NEW
    TridiagonalMatrix<double,complex<double> >(A_original->size(0));
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
  char trans='N';
  if (to==Factorization::TRANSPOSE) trans='T';
  else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*k);
  double *rwork=OPERATOR_NEW_BRACKET(double,k);
  int info;
  double rcond;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<double,complex<double> > *u2=
      OPERATOR_NEW Vector<double,complex<double> >(k-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,k);
    F77NAME(zgttrf)(k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    CHECK_TEST(info==0);
    F77NAME(zgtcon)(norm,k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,info);
    CHECK_TEST(info==0);
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(zgtrfs)(trans,k,n,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),u2->addr(),ipiv,B.addr(),m,X.addr(),m,
        ferr.addr(),berr.addr(),work,rwork,info);
    } else {
      F77NAME(zgtrfsr)(trans,k,m,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),u2->addr(),ipiv,B.addr(),m,X.addr(),m,
        ferr.addr(),berr.addr(),work,rwork,info);
    }
    CHECK_TEST(info==0);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(zgttrfnp)(k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    CHECK_TEST(info==0);
    F77NAME(zgtconnp)(norm,k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,info);
    CHECK_TEST(info==0);
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(zgtrfsnp)(trans,k,n,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),B.addr(),m,X.addr(),m,ferr.addr(),
        berr.addr(),work,rwork,info);
    } else {
      F77NAME(zgtrfsrnp)(trans,k,m,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),B.addr(),m,X.addr(),m,ferr.addr(),
        berr.addr(),work,rwork,info);
    }
    CHECK_TEST(info==0);
  }
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  delete LU; LU=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> GaussianFactorization<double,complex<double>,
BandMatrix<double,complex<double> > >::GaussianFactorization(
const BandMatrix<double,complex<double> >& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),piv_op(po),equ_op(eo),
equed('N'),colcnd(numeric_limits<double>::infinity()),
rowcnd(numeric_limits<double>::infinity()) {
  int n=A.size(0),nsub=A.subDiags(),nsup=A.supDiags();
  LU=OPERATOR_NEW BandMatrix<double,complex<double> >(n,nsub,
    (po==Factorization::NO_PIVOTING ? nsup : nsub+nsup),double_zero_);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-nsup);
    int iend=min(n-1,j+nsub);
    F77NAME(zcopy)(iend-ibeg+1,A.addr(ibeg,j),1,LU->addr(ibeg,j),1);
  }

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    r=OPERATOR_NEW Vector<double,double>(n);
    c=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    F77NAME(zgbequ)(n,n,LU->subDiags(),LU->supDiags(),LU->addr(),
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
      F77NAME(zlaqgb)(n,n,nsub,nsub+nsup,LU->addr(),LU->bands(),
        r->addr(),c->addr(),rowcnd,colcnd,amax,equed);
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
      F77NAME(zgbtf2np)(n,n,nsub,nsup,LU->addr(),LU->bands(),info);
      break;
    }
    default: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(zgbtrf)(n,n,nsub,nsup,LU->addr(),LU->bands(),ipiv,info);
    }
  }
  double *work=0;
  if (info>0) { // see dgbsvx
    rpvgrw=F77NAME(zlantb)('M','U','N',info,min(info-1,LU->supDiags()),
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
    rpvgrw=F77NAME(zlantb)('M','U','N',n,LU->supDiags(),LU->addr(),
      LU->bands(),work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
  }
}

template<> void GaussianFactorization<double,complex<double>,
BandMatrix<double,complex<double> > >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,Factorization::TRANSPOSE_OPTION to) {
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
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
    if (ipiv!=0) {
      F77NAME(zgbtrs)('N',n,A_original->subDiags(),A_original->supDiags(),
        1,LU->addr(),LU->bands(),ipiv,x.addr(),n,info);
    } else {
      F77NAME(zgbtrsnp)('N',n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),x.addr(),n,info);
    }
    CHECK_TEST(info==0);
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ci=c->addr();
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
  } else {
    // A^T X = B ==> U^T L^T Q^T R^{-1} X = C B
    // A^H X = B ==> U^H L^H Q^T R^{-1} X = C B
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ci=c->addr();
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
    char trans=(to==Factorization::TRANSPOSE ? 'T' : 'C');
    if (ipiv!=0) {
      F77NAME(zgbtrs)(trans,n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
        x.addr(),n,info);
    } else {
      F77NAME(zgbtrsnp)(trans,n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),x.addr(),n,info);
    }
    CHECK_TEST(info==0);
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const double *ri=r->addr();
      complex<double> *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
  }
}

template<> void GaussianFactorization<double,complex<double>,
BandMatrix<double,complex<double> > >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Factorization::TRANSPOSE_OPTION to,
Factorization::SIDE_OPTION so) {
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
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
      if (ipiv!=0) {
        F77NAME(zgbtrs)('N',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),ipiv,
          X.addr(),m,info);
      } else {
        F77NAME(zgbtrsnp)('N',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),X.addr(),m,
          info);
      }
      CHECK_TEST(info==0);
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ci=c->addr();
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
    } else {
      // A^T X = B ==> U^T L^T Q^T R^{-1} X = C B
      // A^H X = B ==> U^H L^H Q^T R^{-1} X = C B
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ci=c->addr();
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
      char trans=(to==Factorization::TRANSPOSE ? 'T' : 'C');
      if (ipiv!=0) {
        F77NAME(zgbtrs)(trans,m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),ipiv,
          X.addr(),m,info);
      } else {
        F77NAME(zgbtrsnp)(trans,m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),X.addr(),m,
          info);
      }
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const double *ri=r->addr();
          complex<double> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
    }
  } else {
    CHECK_SAME(n,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) {
      // X A = B ==> A^T X^T = B^T
      //         ==> U^T L^T Q^T R^{-1} X^T = ( B C )^T
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*c)[j],X.addr(0,j),1);
      }
      Vector<double,complex<double> > *x=
        OPERATOR_NEW Vector<double,complex<double> >(n);
      for (int i=0;i<m;i++) {
        F77NAME(zcopy)(n,X.addr(i,0),m,x->addr(),1);
        if (ipiv!=0) {
          F77NAME(zgbtrs)('T',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
            x->addr(),n,info);
        } else {
          F77NAME(zgbtrsnp)('T',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),x->addr(),n,
            info);
        }
        F77NAME(zcopy)(n,x->addr(),1,X.addr(i,0),m);
      }
      delete x; x=0;
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*r)[j],X.addr(0,j),1);
      }
    } else {
      // X A^T = B ==> A X^T = B^T
      //           ==> L U C^{-1} X^T = Q^T ( B R )^T
      // X A^H = B ==> A X^H = B^H
      //           ==> L U C^{-1} X^H = Q^T ( B R )^H
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*r)[j],X.addr(0,j),1);
      }
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<double> *Xij=X.addr();
        for (int i=0;i<m*n;i++,Xij++) *Xij=conj(*Xij);
      }
      Vector<double,complex<double> > *x=
        OPERATOR_NEW Vector<double,complex<double> >(n);
      for (int i=0;i<m;i++) {
        F77NAME(zcopy)(n,X.addr(i,0),m,x->addr(),1);
        if (ipiv!=0) {
          F77NAME(zgbtrs)('N',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
            x->addr(),n,info);
        } else {
          F77NAME(zgbtrsnp)('N',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),x->addr(),n,
            info);
        }
        F77NAME(zcopy)(n,x->addr(),1,X.addr(i,0),m);
      }
      delete x; x=0;
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*c)[j],X.addr(0,j),1);
      }
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<double> *Xij=X.addr();
        for (int i=0;i<m*n;i++,Xij++) *Xij=conj(*Xij);
      }
    }
  }
}

template<> double GaussianFactorization<double,complex<double>,
BandMatrix<double,complex<double> > >::reciprocalConditionNumber(
Factorization::CONDITION_NUMBER_NORM cnn) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  double anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info=0;
  if (piv_op==Factorization::PIVOT_ROWS) {
    F77NAME(zgbcon)(norm,n,A_original->subDiags(),A_original->supDiags(),
      LU->addr(),LU->bands(),ipiv,anorm,rcond,work,rwork,info);
  } else {
    F77NAME(zgbconnp)(norm,n,A_original->subDiags(),
      A_original->supDiags(),LU->addr(),LU->bands(),anorm,rcond,work,
      rwork,info);
  }
  CHECK_SAME(info,0)
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void GaussianFactorization<double,complex<double>,
BandMatrix<double,complex<double> > >::improve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double &berr,double &ferr,
Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(A_original->bands()+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,complex<double> > residual(n);
  Vector<double,double> work(n);
  char trans='N';
  if (to==Factorization::TRANSPOSE) trans='T';
  else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(zgbmv)(trans,n,n,A_original->subDiags(),
      A_original->supDiags(),double_mone_,A_original->addr(),
      A_original->bands(),x.addr(),1,double_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(zgbamv)(trans,n,n,A_original->subDiags(),
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
    F77NAME(zaxpy)(n,double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  Vector<double,complex<double> > v(n);
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(zlacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual,to);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual,to);
    }
  }
  int i=F77NAME(izamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void GaussianFactorization<double,complex<double>,
BandMatrix<double,complex<double> > >::improve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Vector<double,double> &berr,
Vector<double,double> &ferr,Factorization::TRANSPOSE_OPTION to,
Factorization::SIDE_OPTION so) {
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
  Vector<double,complex<double> > x(k);
  Vector<double,complex<double> > rhs(k);
  Vector<double,complex<double> > residual(k);
  Vector<double,double> work(k);
  Vector<double,complex<double> > v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  char trans='N';
  Factorization::TRANSPOSE_OPTION lto=Factorization::NO_TRANSPOSE;
  if (so==Factorization::LEFT_SIDE) {
    if (to==Factorization::TRANSPOSE) trans='T';
    else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
    lto=to;
  } else {
    if (to==Factorization::NO_TRANSPOSE) {
      trans='T';
      lto=Factorization::TRANSPOSE;
    } else {
      trans='N';
      lto=Factorization::NO_TRANSPOSE;
    }
  }
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(zcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(zcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        const complex<double> *Xji=X.addr(j,0);
        complex<double> *xi=x.addr();
        const complex<double> *Bji=B.addr(j,0);
        complex<double> *rhsi=rhs.addr();
        for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
          *xi=conj(*Xji);
          *rhsi=conj(*Bji);
        }
      } else {
        F77NAME(zcopy)(k,X.addr(j,0),m,x.addr(),1);
        F77NAME(zcopy)(k,B.addr(j,0),m,rhs.addr(),1);
      }
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(zgbmv)(trans,k,k,A_original->subDiags(),
        A_original->supDiags(),double_mone_,
        A_original->addr(),A_original->bands(),x.addr(),1,double_one_,
        residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77NAME(zgbamv)(trans,k,k,A_original->subDiags(),
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
      solve(residual,residual,lto);
      F77NAME(zaxpy)(k,complex_double_one_,residual.addr(),1,x.addr(),1);
      lstres=s;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(zcopy)(k,x.addr(),1,X.addr(0,j),1);
      } else F77NAME(zcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<double> *Xji=X.addr(j,0);
        const complex<double> *xi=x.addr();
        for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
      } else {
        F77NAME(zcopy)(k,x.addr(),1,X.addr(j,0),m);
      }
    }

    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(zlacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual,lto);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual,lto);
      }
    }
    int i=F77NAME(izamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template class GaussianFactorization<double,complex<double>,
  SquareMatrix<double,complex<double> > >;
template void testGaussianFactorization(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "CholeskyFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> CholeskyFactorization<double,complex<double>,
SymmetricPositiveMatrix<double,complex<double> > >::CholeskyFactorization(
const SymmetricPositiveMatrix<double,complex<double> >& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equed('N'),scond(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  L=OPERATOR_NEW SymmetricPositiveMatrix<double,complex<double> >(n);
  L->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,3*n);
    int info;
    F77NAME(zheequb)('L',n,L->addr(),n,s->addr(),scond,amax,work,info);
    delete [] work; work=0;
    equed='N';
    F77NAME(zlaqhe)('L',n,L->addr(),n,s->addr(),scond,amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=L->normInfinity();
  anormm=L->normMaxEntry();
  anormo=L->normOne();

  int info=0;
  F77NAME(zpotrf)('L',n,L->addr(),n,info);
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,4*n);
  if (info!=0) { // see dposvxx
    rpvgrw=F77_NAME(zla_porpvgrw)('L',info,A_original->addr(),n,
      L->addr(),n, work);
    delete L; L=0;
    if (s!=0) delete s; s=0;
  } else {
    rpvgrw=F77_NAME(zla_porpvgrw)('L',n,A_original->addr(),n,
      L->addr(),n, work);
  }
  delete [] work; work=0;
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveMatrix<double,complex<double> > >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x) {
  CHECK_POINTER(L);
//constructor factored S A S = M D M^T
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,L->size(0));
  if (&x!=&b) x.copy(b);
//A x = b ==> L L^T S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    complex<double> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  int info;
  F77NAME(zpotrs)('L',n,1,L->addr(),n,x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    complex<double> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveMatrix<double,complex<double> > >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Factorization::SIDE_OPTION so) {
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
        complex<double> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    int info;
    F77NAME(zpotrs)('L',m,n,L->addr(),m,X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        complex<double> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,L->size(0));
//  X A = B ==> X S^{-1} L L^T = B S
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*s)[j],X.addr(0,j),1);
    }
    F77NAME(ztrsm)('R','L','C','N',m,n,complex_double_one_,L->addr(),n,
      X.addr(),m);
    F77NAME(ztrsm)('R','L','N','N',m,n,complex_double_one_,L->addr(),n,
      X.addr(),m);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> double CholeskyFactorization<double,complex<double>,
SymmetricPositiveMatrix<double,complex<double> > >::
reciprocalConditionNumber() {
  CHECK_POINTER(L);
  int n=L->size(0);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info=0;
  F77NAME(zpocon)('L',n,L->addr(),n,anormo,rcond,work,rwork,info);
  CHECK_SAME(info,0)
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveMatrix<double,complex<double> > >::improve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double &berr,double &ferr) {
  CHECK_POINTER(L);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(n+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,complex<double> > residual(n);
  Vector<double,double> work(n);
  Vector<double,complex<double> > v(n);
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(zhemv)('L',n,complex_double_mone_,A_original->addr(),n,
      x.addr(),1,complex_double_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77_NAME(zla_heamv)(F77NAME(ilauplo)('L'),n,double_one_,
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
    F77NAME(zaxpy)(n,complex_double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  complex<double> *residuali=residual.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(zlacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(izamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isave; isave=0;
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveMatrix<double,complex<double> > >::improve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Vector<double,double> &berr,
Vector<double,double> &ferr,Factorization::SIDE_OPTION so) {
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
  Vector<double,complex<double> > x(k);
  Vector<double,complex<double> > rhs(k);
  Vector<double,complex<double> > residual(k);
  Vector<double,double> work(k);
  Vector<double,complex<double> > v(k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(zcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(zcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      const complex<double> *Xji=X.addr(j,0);
      complex<double> *xi=x.addr();
      const complex<double> *Bji=B.addr(j,0);
      complex<double> *rhsi=rhs.addr();
      for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
        *xi=conj(*Xji);
        *rhsi=conj(*Bji);
      }
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(zhemv)('L',k,complex_double_mone_,A_original->addr(),k,
        x.addr(),1,complex_double_one_,residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77_NAME(zla_heamv)(F77NAME(ilauplo)('L'),k,double_one_,
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
      F77NAME(zaxpy)(k,complex_double_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(zcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      complex<double> *Xji=X.addr(j,0);
      const complex<double> *xi=x.addr();
      for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
    }

    complex<double> *residuali=residual.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(zlacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(izamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr[j]/=lstres;
  }
  delete [] isave; isave=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> CholeskyFactorization<double,complex<double>,
SymmetricPositiveTridiagonalMatrix<double,complex<double> > >::
CholeskyFactorization(
const SymmetricPositiveTridiagonalMatrix<double,complex<double> >& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equ_op(Factorization::NO_EQUILIBRATION),equed('N'),
scond(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  anormi=A.normInfinity();
  anormm=A.normMaxEntry();
  anormo=A.normOne();
  rpvgrw=double_one_; // not computed
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveTridiagonalMatrix<double,complex<double> > >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x) {
//constructor did not factor
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<double,complex<double> >(n);
  L->copy(*A_original);
  if (&x!=&b) x.copy(b);
  int info;
  F77NAME(zptsv)(n,1,L->diagonalAddr(0),L->lowerDiagonalAddr(0),
    x.addr(),n,info);
  CHECK_TEST(info==0);
  delete L; L=0;
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveTridiagonalMatrix<double,complex<double> > >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Factorization::SIDE_OPTION so) {
//constructor did not factor
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<double,complex<double> >(
    A_original->size(0));
  L->copy(*A_original);
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,L->size(0));
    F77NAME(zptsv)(m,n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),
      X.addr(),m,info);
  } else {
    CHECK_SAME(n,L->size(0));
    F77NAME(zpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
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

template<> double CholeskyFactorization<double,complex<double>,
SymmetricPositiveTridiagonalMatrix<double,complex<double> > >::
reciprocalConditionNumber() {
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<double,complex<double> >(
    A_original->size(0));
  L->copy(*A_original);
  int n=L->size(0);
  double rcond;
  double *work=OPERATOR_NEW_BRACKET(double,n);
  int info=0;
  F77NAME(zpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  F77NAME(zptcon)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,work,info);
  CHECK_SAME(info,0)
  delete L; L=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveTridiagonalMatrix<double,complex<double> > >::improve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double &berr,double &ferr) {
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<double,complex<double> >(
    A_original->size(0));
  L->copy(*A_original);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info;
  double rcond;
  F77NAME(zpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  CHECK_TEST(info==0);
  F77NAME(zptcon)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,rwork,info);
  CHECK_TEST(info==0);
  F77NAME(zptrfs)('L',n,1,A_original->diagonalAddr(0),
    A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
    L->lowerDiagonalAddr(0),b.addr(),n,x.addr(),n,&ferr,&berr,
    work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  delete L; L=0;
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveTridiagonalMatrix<double,complex<double> > >::improve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Vector<double,double> &berr,
Vector<double,double> &ferr,Factorization::SIDE_OPTION so) {
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<double,complex<double> >(
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
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,k);
  double *rwork=OPERATOR_NEW_BRACKET(double,k);
  int info;
  double rcond;
  F77NAME(zpttrf)(k,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  CHECK_TEST(info==0);
  F77NAME(zptcon)(k,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,rwork,info);
  CHECK_TEST(info==0);
  if (so==Factorization::LEFT_SIDE) {
    F77NAME(zptrfs)('L',k,n,A_original->diagonalAddr(0),
      A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
      L->lowerDiagonalAddr(0),B.addr(),m,X.addr(),m,ferr.addr(),
      berr.addr(),work,rwork,info);
  } else {
    F77NAME(zptrfsr)('L',k,m,A_original->diagonalAddr(0),
      A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
      L->lowerDiagonalAddr(0),B.addr(),m,X.addr(),m,ferr.addr(),
      berr.addr(),work,rwork,info);
  }
  CHECK_TEST(info==0);
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  delete L; L=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> CholeskyFactorization<double,complex<double>,
SymmetricPositiveBandMatrix<double,complex<double> > >::
CholeskyFactorization(
const SymmetricPositiveBandMatrix<double,complex<double> >& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equ_op(eo),equed('N'),scond(numeric_limits<double>::infinity()) {
  int n=A.size(0),nsub=A.subDiags();
  L=OPERATOR_NEW SymmetricPositiveBandMatrix<double,complex<double> >(n,
    nsub,complex_double_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(zcopy)(min(n-j,nsub+1),A.addr(j,j),1,L->addr(j,j),1);
  }

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    F77NAME(zpbequ)('L',n,L->subDiags(),L->addr(),L->bands(),s->addr(),
      scond,amax,info);
    equed='N';
    F77NAME(zlaqsb)('L',n,nsub,L->addr(),L->bands(),s->addr(),scond,
      amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=L->normInfinity();
  anormm=L->normMaxEntry();
  anormo=L->normOne();

  int info=0;
  F77NAME(zpbtrf)('L',n,nsub,L->addr(),L->bands(),info);
  double *work=0;
  if (info>0) { // see dgbsvx
    rpvgrw=F77NAME(zlansb)('M','L',info,min(info-1,L->subDiags()),
      L->addr(),L->bands(),work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
    cerr << "zero pivot in CholeskyFactorization::CholeskyFactorization"
       << "\n zero pivot number,reciprocal pivot growth = " << info
       << " " << rpvgrw << endl;
    if (L) delete L; L=0;
    if (s) delete s; s=0;
    A_original=0;
  } else {
    rpvgrw=F77NAME(zlansb)('M','L',n,L->subDiags(),L->addr(),
      L->bands(),work);
    rpvgrw=(rpvgrw==double_zero_ ? double_one_ : anormm / rpvgrw);
  }
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveBandMatrix<double,complex<double> > >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^H
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  if (&x!=&b) x.copy(b);
  int info;
//A x = b ==> L L^H S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    complex<double> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  {
  F77NAME(zpbtrs)('L',n,L->subDiags(),1,L->addr(),L->bands(),
    x.addr(),n,info);
  }
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    complex<double> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveBandMatrix<double,complex<double> > >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Factorization::SIDE_OPTION so) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^H
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
        complex<double> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    F77NAME(zpbtrs)('L',m,L->subDiags(),n,L->addr(),L->bands(),
      X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        complex<double> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,L->size(0));
//  X A = B ==> X S^{-1} L L^H = B S
//          ==> L L^H S^{-1} X^H = S B^H
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*s)[j],X.addr(0,j),1);
    }
    Vector<double,complex<double> > *x=
      OPERATOR_NEW Vector<double,complex<double> >(n);
    for (int i=0;i<m;i++) {
      complex<double> *Xij=X.addr(i,0);
      complex<double> *xj=x->addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *xj=conj(*Xij);
      F77NAME(zpbtrs)('L',n,L->subDiags(),1,L->addr(),L->bands(),
        x->addr(),n,info);
      xj=x->addr();
      Xij=X.addr(i,0);
      for (int j=0;j<n;j++,Xij+=m,xj++) *Xij=conj(*xj);
    }
    delete x; x=0;
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> double CholeskyFactorization<double,complex<double>,
SymmetricPositiveBandMatrix<double,complex<double> > >::
reciprocalConditionNumber() {
  CHECK_POINTER(L);
  int n=L->size(0);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  double *rwork=OPERATOR_NEW_BRACKET(double,n);
  int info=0;
  F77NAME(zpbcon)('L',n,L->subDiags(),L->addr(),L->bands(),anormo,
    rcond,work,rwork,info);
  CHECK_SAME(info,0)
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveBandMatrix<double,complex<double> > >::improve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double &berr,double &ferr) {
  CHECK_POINTER(L);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(A_original->bands()+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,complex<double> > residual(n);
  Vector<double,double> work(n);
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(zhbmv)('L',n,A_original->subDiags(),complex_double_mone_,
      A_original->addr(),A_original->bands(),x.addr(),1,
      complex_double_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(zhbamv)('L',n,A_original->subDiags(),double_one_,
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
    F77NAME(zaxpy)(n,complex_double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  Vector<double,complex<double> > v(n);
  int kase=0;
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(zlacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(izamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isave; isave=0;
}

template<> void CholeskyFactorization<double,complex<double>,
SymmetricPositiveBandMatrix<double,complex<double> > >::improve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,
Vector<double,double> &berr,
Vector<double,double> &ferr,Factorization::SIDE_OPTION so) {
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
  Vector<double,complex<double> > x(k);
  Vector<double,complex<double> > rhs(k);
  Vector<double,complex<double> > residual(k);
  Vector<double,double> work(k);
  Vector<double,complex<double> > v(k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(zcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(zcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      const complex<double> *Xji=X.addr(j,0);
      complex<double> *xi=x.addr();
      const complex<double> *Bji=B.addr(j,0);
      complex<double> *rhsi=rhs.addr();
      for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
        *xi=conj(*Xji);
        *rhsi=conj(*Bji);
      }
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(zhbmv)('L',k,A_original->subDiags(),complex_double_mone_,
        A_original->addr(),A_original->bands(),x.addr(),1,
        complex_double_one_,residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77NAME(zhbamv)('L',k,A_original->subDiags(),double_one_,
        A_original->addr(),A_original->bands(),x.addr(),1,
        double_one_,work.addr(),1);
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
      F77NAME(zaxpy)(k,complex_double_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(zcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      complex<double> *Xji=X.addr(j,0);
      const complex<double> *xi=x.addr();
      for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
    }

    complex<double> *residuali=residual.addr();
    double *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(zlacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(izamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr/=lstres;
  }
  delete [] isave; isave=0;
}

template class CholeskyFactorization<double,complex<double>,
  SymmetricPositiveMatrix<double,complex<double> > >;
template void testCholeskyFactorization(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "MDMtFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> MDMtFactorization<double,complex<double> >::MDMtFactorization(
const SymmetricMatrix<double,complex<double> >& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),MD(0),ipiv(0),s(0),
equed('N'),scond(numeric_limits<double>::infinity()) {
  int n=A.size(0);
  MD=OPERATOR_NEW SymmetricMatrix<double,complex<double> >(n);
  MD->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<double,double>(n);
    double amax;
    int info;
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,3*n);
    F77NAME(zheequb)('L',n,MD->addr(),n,s->addr(),scond,amax,work,info);
    delete [] work; work=0;
    equed='N';
    F77NAME(zlaqsy)('L',n,MD->addr(),n,s->addr(),scond,amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=MD->normInfinity();
  anormm=MD->normMaxEntry();
  anormo=MD->normOne();

  int info=0;
  ipiv=OPERATOR_NEW_BRACKET(int,n);
  complex<double> workq=complex_double_zero_;
  int lwork=-1;
  F77NAME(zhetrf)('L',n,MD->addr(),n,ipiv,&workq,lwork,info);
  lwork=max(2*n,static_cast<int>(workq.real()));
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zhetrf)('L',n,MD->addr(),n,ipiv,work,lwork,info);
  if (info!=0) { // see dposvxx
    rpvgrw=F77_NAME(zla_syrpvgrw)('L',n,info,A_original->addr(),n,
      MD->addr(),n,ipiv,work);
    delete MD; MD=0;
    if (s!=0) delete s; s=0;
    if (ipiv!=0) delete ipiv; ipiv=0;
  } else {
    rpvgrw=F77_NAME(zla_syrpvgrw)('L',n,n,A_original->addr(),n,
      MD->addr(),n,ipiv,work);
  }
  delete [] work; work=0;
}

template<> void MDMtFactorization<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x) {
  CHECK_POINTER(MD);
//constructor factored S A S = M D M^H
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,MD->size(0));
  if (&x!=&b) x.copy(b);
//A x = b ==> M D M^H S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    complex<double> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  int info;
  F77NAME(zhetrs)('L',n,1,MD->addr(),n,ipiv,x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const double *si=s->addr();
    complex<double> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void MDMtFactorization<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Factorization::SIDE_OPTION so) {
  CHECK_POINTER(MD);
//constructor factored S A S = M D M^H
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  if (&X!=&B) X.copy(B);
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,MD->size(0));
//  A X = B ==> M D M^H S^{-1} X = S B
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        complex<double> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    int info;
    F77NAME(zhetrs)('L',m,n,MD->addr(),m,ipiv,X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const double *si=s->addr();
        complex<double> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,MD->size(0));
//  X A = B ==> X S^{-1} M D M^T = B S
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*s)[j],X.addr(0,j),1);
    }
    Vector<double,complex<double> > x(n);
    for (int i=0;i<m;i++) {
      complex<double> *Xij=X.addr(i,0);
      complex<double> *xj=x.addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *xj=conj(*Xij);
      int info;
      F77NAME(zhetrs)('L',n,1,MD->addr(),n,ipiv,x.addr(),n,info);
      CHECK_TEST(info==0);
      Xij=X.addr(i,0);
      xj=x.addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *Xij=conj(*xj);
    }
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(zdscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> double
MDMtFactorization<double,complex<double> >::reciprocalConditionNumber() {
  CHECK_POINTER(MD);
  int n=MD->size(0);
  double rcond;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,2*n);
  int info=0;
  F77NAME(zhecon)('L',n,MD->addr(),n,ipiv,anormo,rcond,work,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  return rcond;
}

template<> void MDMtFactorization<double,complex<double> >::improve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double &berr,double &ferr) {
  CHECK_POINTER(MD);
  int n=MD->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  double nz=static_cast<double>(n+1);
  double eps=F77NAME(dlamch)('E');
  double safe1=nz*F77NAME(dlamch)('S');
  double safe2=safe1/eps;
  int itmax=5;
  Vector<double,complex<double> > residual(n);
  Vector<double,double> work(n);
  Vector<double,complex<double> > v(n);
  double lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(zhemv)('L',n,complex_double_mone_,A_original->addr(),n,
      x.addr(),1,complex_double_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77_NAME(zla_heamv)(F77NAME(ilauplo)('L'),n,double_one_,
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
    F77NAME(zaxpy)(n,complex_double_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  complex<double> *residuali=residual.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(zlacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(izamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=double_zero_) ferr/=lstres;
  delete [] isave; isave=0;
}

template<> void MDMtFactorization<double,complex<double> >::improve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Vector<double,double> &berr,
Vector<double,double> &ferr,Factorization::SIDE_OPTION so) {
//compare to dposvx.f
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
  Vector<double,complex<double> > x(k);
  Vector<double,complex<double> > rhs(k);
  Vector<double,complex<double> > residual(k);
  Vector<double,double> work(k);
  Vector<double,complex<double> > v(k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(zcopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(zcopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      const complex<double> *Xji=X.addr(j,0);
      complex<double> *xi=x.addr();
      const complex<double> *Bji=B.addr(j,0);
      complex<double> *rhsi=rhs.addr();
      for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
        *xi=conj(*Xji);
        *rhsi=conj(*Bji);
      }
    }
    double lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(zhemv)('L',k,complex_double_mone_,A_original->addr(),k,
        x.addr(),1,complex_double_one_,residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77_NAME(zla_heamv)(F77NAME(ilauplo)('L'),k,double_one_,
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
      F77NAME(zaxpy)(k,complex_double_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(zcopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      complex<double> *Xji=X.addr(j,0);
      const complex<double> *xi=x.addr();
      for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
    }

    complex<double> *residuali=residual.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(zlacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(izamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=double_zero_) ferr[j]/=lstres;
  }
  delete [] isave; isave=0;
}

template class MDMtFactorization<double,complex<double> >;
template void testMDMtFactorization(double,complex<double>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "HouseholderQRFactorization.C"
template<> HouseholderQRFactorization<double,complex<double> >
::HouseholderQRFactorization(const Matrix<double,complex<double> > &A,
Factorization::PIVOT_OPTION po) : QR(0),tau(0),jpvt(0),piv_op(po),
A_original(&A),iascl(0),ascl(double_one_),
anrm(numeric_limits<double>::infinity()) {
//TRACER_CALL(t,"HouseholderQRFactorization::HouseholderQRFactorization");
  int m=A.size(0),n=A.size(1);
  QR=OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  QR->copy(A);
  double dw=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);
  anrm=F77NAME(zlange)('M',m,n,A.addr(),m,&dw);
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
    F77NAME(zlascl)('G',0,0,anrm,ascl,m,n,QR->addr(),m,info);
    CHECK_SAME(info,0);
  }

  if (piv_op==Factorization::NO_PIVOTING) {
    complex<double> zw=complex_double_undefined_;
    int lwork=-1;
    if (m>=n) {
      tau=OPERATOR_NEW Vector<double,complex<double> >(n);
      F77NAME(zgeqrf)(m,n,QR->addr(),m,tau->addr(),&zw,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zgeqrf)(m,n,QR->addr(),m,tau->addr(),work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
    } else {
      tau=OPERATOR_NEW Vector<double,complex<double> >(m);
      F77NAME(zgelqf)(m,n,QR->addr(),m,tau->addr(),&zw,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zgelqf)(m,n,QR->addr(),m,tau->addr(),work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
    }
  } else {
    int mn=min(m,n);
    int nb=F77NAME(ilaenv)(1,"ZGEQRF"," ",m,n,-1,1);
    int lwkmin=mn+max(2*mn,n+1);
    int lwork=max(lwkmin,mn+2*n+nb*(n+1));
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    tau=OPERATOR_NEW Vector<double,complex<double> >(mn);
    jpvt=OPERATOR_NEW_BRACKET(int,n);
    double *rwork=OPERATOR_NEW_BRACKET(double,2*n);
    F77NAME(zgeqp3)(m,n,QR->addr(),m,jpvt,tau->addr(),work,lwork,
      rwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    delete [] rwork; rwork=0;
  }
}

template<> OrthogonalMatrix<double,complex<double> >* 
HouseholderQRFactorization<double,complex<double> >::orthogonalPart()
const {
//TRACER_CALL(t,"HouseholderQRFactorization::orthogonalPart");
  int m=QR->size(0),n=QR->size(1);
  complex<double> w=complex_double_undefined_;
  if (m>=n) {
    OrthogonalMatrix<double,complex<double> > *Q=
      OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(m,m);
    Q->copyFrom('A',m,n,*QR);

    int lwork=-1;
    int info;
    F77NAME(zungqr)(m,m,n,Q->addr(),m,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w.real());
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zungqr)(m,m,n,Q->addr(),m,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    return Q;
  } else {
    OrthogonalMatrix<double,complex<double> > *Q=
      OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(n,n);
    Q->copyFrom('A',m,n,*QR);

    int lwork=-1;
    int info;
    F77NAME(zunglq)(n,n,m,Q->addr(),n,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w.real());
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunglq)(n,n,m,Q->addr(),n,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    return Q;
  }
}

template<> double
HouseholderQRFactorization<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,Factorization::TRANSPOSE_OPTION tr)
const {
//TRACER_CALL(t,"HouseholderQRFactorization::solve");
  int m=QR->size(0),n=QR->size(1);
  double dw=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(zlange)('M',brow,1,b.addr(),b.size(),&dw);
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
  complex<double> zw=complex_double_undefined_;
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
        Vector<double,complex<double> > xtmp(m);
        xtmp.copy(b);
        if (ibscl!=0) {
          F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      Q' b = [ y \\ z ]:
        F77NAME(zunmqr)('L','C',m,1,n,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,&zw,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zunmqr)('L','C',m,1,n,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,work,lwork,info);
        CHECK_SAME(info,0)

//      solve R x = y:
        F77NAME(ztrsv)('U','N','N',n,QR->addr(),m,xtmp.addr(),1);
        delete [] work; work=0;
        x.copyFrom(n,xtmp);
        residual_norm=F77NAME(dznrm2)(m-n,xtmp.addr(n),1);
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

        x=complex_double_zero_;
        x.copyFrom(m,b);
        if (ibscl!=0) {
          F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,x.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      solve R' v = b
        F77NAME(ztrsv)('U','C','N',n,QR->addr(),m,x.addr(),1);
//      Q [v\\0]
        F77NAME(zunmqr)('L','N',m,1,n,QR->addr(),m,tau->addr(),
                        x.addr(),m,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zunmqr)('L','N',m,1,n,QR->addr(),m,tau->addr(),
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
        x=complex_double_zero_;
        x.copyFrom(m,b);
        if (ibscl!=0) {
          F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,x.addr(),m,info);
          CHECK_SAME(info,0)
        }
        F77NAME(ztrsv)('L','N','N',m,QR->addr(),m,x.addr(),1);
        F77NAME(zunmlq)('L','C',n,1,m,QR->addr(),m,tau->addr(),
                        x.addr(),n,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zunmlq)('L','C',n,1,m,QR->addr(),m,tau->addr(),
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
        Vector<double,complex<double> > xtmp(n);
        xtmp.copy(b);
        if (ibscl!=0) {
          F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
          CHECK_SAME(info,0)
        }
        F77NAME(zunmlq)('L','N',m,1,m,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zunmlq)('L','N',m,1,m,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        F77NAME(ztrsv)('L','C','N',m,QR->addr(),m,xtmp.addr(),1);
        x.copyFrom(m,xtmp);
        residual_norm=F77NAME(dznrm2)(n-m,xtmp.addr(m),1);
        scllen=m;
      }
    }
  } else {
    int mn=min(m,n);
    if (tr==Factorization::NO_TRANSPOSE) {
      CHECK_SAME(m,b.size())
      CHECK_SAME(n,x.size())
      Vector<double,complex<double> > xtmp(m);
      xtmp.copy(b);
      if (ibscl!=0) {
        F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    Q' b = [ y \\ z ]:
      F77NAME(zunmqr)('L','C',m,1,mn,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,&zw,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zunmqr)('L','C',m,1,mn,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
//    solve R x' = y:
      F77NAME(ztrsv)('U','N','N',mn,QR->addr(),m,xtmp.addr(),1);
//    permute
      for (int i=0;i<n;i++) x[jpvt[i]-1]=xtmp[i];
      residual_norm=F77NAME(dznrm2)(m-mn,xtmp.addr(mn),1);
    } else {
      CHECK_SAME(n,b.size())
      CHECK_SAME(m,x.size())
      Vector<double,complex<double> > xtmp(n);
//    P'b:
      for (int i=0;i<n;i++) xtmp[i]=b[jpvt[i]-1];
      if (ibscl!=0) {
        F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    solve R' v = b
      F77NAME(ztrsv)('U','C','N',mn,QR->addr(),m,xtmp.addr(),1);
      x=complex_double_zero_;
      x.copyFrom(mn,xtmp);
//    Q [v\\0]
      lwork=-1;
      F77NAME(zunmqr)('L','N',m,1,mn,QR->addr(),m,tau->addr(),
                      x.addr(),m,&zw,lwork,info );
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zunmqr)('L','N',m,1,mn,QR->addr(),m,tau->addr(),
                      x.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      if (mn<n) {
        residual_norm=F77NAME(dznrm2)(n-mn,xtmp.addr(mn),1);
      } else residual_norm=double_zero_;
    }
  }
  if (iascl!=0) {
    F77NAME(zlascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(zlascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  return residual_norm;
}

template<> void
HouseholderQRFactorization<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,
Vector<double,complex<double> > &residual_norm,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::solve");
  int m=QR->size(0),n=QR->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  CHECK_SAME(k,residual_norm.size())
  double dw=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(zlange)('M',brow,k,B.addr(),B.size(0),&dw);
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
  complex<double> zw=complex_double_undefined_;
  int lwork=-1;
  int scllen=0;
  if (piv_op==Factorization::NO_PIVOTING) {
    if (m>=n) {
      if (tr==Factorization::NO_TRANSPOSE) { // min || B - A x ||
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        Matrix<double,complex<double> > Xtmp(m,k);
        Xtmp.copy(B);
        if (ibscl!=0) {
          F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      Q' B = [ Y \\ Z ]:
        F77NAME(zunmqr)('L','C',m,k,n,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,&zw,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zunmqr)('L','C',m,k,n,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,work,lwork,info);
        CHECK_SAME(info,0)

//      solve R X = Y:
        F77NAME(ztrsm)('L','U','N','N',n,k,complex_double_one_,
          QR->addr(),m,Xtmp.addr(),m);
        delete [] work; work=0;
        X.copyFrom('A',n,k,Xtmp);
        for (int j=0;j<k;j++) {
          residual_norm[j]=F77NAME(dznrm2)(m-n,Xtmp.addr(n,j),1);
        }
        scllen=n;
      } else { // min || x || s.t. A^t x = B
        CHECK_SAME(n,B.size(0))
        CHECK_SAME(m,X.size(0))

        X=complex_double_zero_;
        X.copyFrom('A',m,k,B);
        if (ibscl!=0) {
          F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,k,X.addr(),m,info);
          CHECK_SAME(info,0)
        }
        F77NAME(ztrsm)('L','U','C','N',n,k,complex_double_one_,
          QR->addr(),m,X.addr(),m);
        F77NAME(zunmqr)('L','N',m,k,n,QR->addr(),m,tau->addr(),
                        X.addr(),m,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zunmqr)('L','N',m,k,n,QR->addr(),m,tau->addr(),
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
        X=complex_double_zero_;
        X.copyFrom('A',m,k,B);
        if (ibscl!=0) {
          F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,k,X.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      solve L X = B:
        F77NAME(ztrsm)('L','L','N','N',m,k,complex_double_one_,
          QR->addr(),m,X.addr(),n);
//      Q' V
        F77NAME(zunmlq)('L','C',n,k,m,QR->addr(),m,tau->addr(),
                        X.addr(),n,&zw,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zunmlq)('L','C',n,k,m,QR->addr(),m,tau->addr(),
                        X.addr(),n,work,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=double_zero_;
        scllen=n;
      } else { // // min || B - A' X ||
        CHECK_SAME(n,B.size(0))
        CHECK_SAME(m,X.size(0))

        Matrix<double,complex<double> > Xtmp(n,k);
        Xtmp.copy(B);
        if (ibscl!=0) {
          F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
          CHECK_SAME(info,0)
        }
        F77NAME(zunmlq)('L','N',m,k,m,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
        F77NAME(zunmlq)('L','N',m,k,m,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;

        F77NAME(ztrsm)('L','L','C','N',m,k,complex_double_one_,
          QR->addr(),m,Xtmp.addr(),n);
        X.copyFrom('A',m,k,Xtmp);
        residual_norm=double_zero_;
        for (int j=0;j<k;j++) {
          residual_norm[j]=F77NAME(dznrm2)(n-m,Xtmp.addr(m,j),1);
        }
        scllen=m;
      }
    }
  } else {
    int mn=min(m,n);
    if (tr==Factorization::NO_TRANSPOSE) {
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<double,complex<double> > Xtmp(m,k);
      Xtmp.copy(B);
      if (ibscl!=0) {
        F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    Q' B = [ Y \\ Z ]:
      F77NAME(zunmqr)('L','C',m,k,mn,QR->addr(),m,tau->addr(),
      Xtmp.addr(),m,&zw,lwork,info);
      CHECK_SAME(info,0)
      
      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zunmqr)('L','C',m,k,mn,QR->addr(),m,tau->addr(),
      Xtmp.addr(),m,work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
//    solve R V = Y:
      F77NAME(ztrsv)('U','N','N',mn,QR->addr(),m,Xtmp.addr(),1);
//    P V:
      for (int j=0;j<k;j++) {
        for (int i=0;i<n;i++) X(jpvt[i]-1,j)=Xtmp(i,j);
        residual_norm[j]=
          F77NAME(dznrm2)(m-mn,Xtmp.addr(mn,j),1);
      }
    } else {
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<double,complex<double> > Xtmp(n,k);
//    P'B:
      for (int j=0;j<k;j++) {
        for (int i=0;i<n;i++) Xtmp(i,j)=B(jpvt[i]-1,j);
      }
      if (ibscl!=0) {
        F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    solve R' V = P' B
      F77NAME(ztrsv)('U','C','N',mn,QR->addr(),m,Xtmp.addr(),1);
      X=complex_double_zero_;
      X.copyFrom('A',mn,k,Xtmp);
//    Q [V\\0]
      lwork=-1;
      F77NAME(zunmqr)('L','N',m,k,mn,QR->addr(),m,tau->addr(),
        X.addr(),m,&zw,lwork,info );
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zunmqr)('L','N',m,k,mn,QR->addr(),m,tau->addr(),
        X.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      if (mn<n) {
        for (int j=0;j<k;j++) {
          residual_norm[j]=
            F77NAME(dznrm2)(n-mn,Xtmp.addr(mn,j),1);
        }
      } else residual_norm=double_zero_;
    }
  }
  if (iascl!=0) {
    F77NAME(zlascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(zlascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
}

template<> void
HouseholderQRFactorization<double,complex<double> >::improve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::improve");
  CHECK_TEST(piv_op==Factorization::NO_PIVOTING)
  int m=QR->size(0),n=QR->size(1);
  double dw=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(zlange)('M',brow,1,b.addr(),b.size(),&dw);
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
    F77NAME(zlascl)('G',0,0,bnrm,bscl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (iascl!=0) {
    F77NAME(zlascl)('G',0,0,ascl,anrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  complex<double> zw=complex_double_undefined_;
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
      Vector<double,complex<double> > r(m);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    r=b-Ax
      F77NAME(zgemv)('N',m,n,complex_double_mone_,A_original->addr(),m,
        x.addr(),1,complex_double_one_,r.addr(),1);
      Vector<double,complex<double> > c(n,complex_double_zero_);
//    -dc=A'r
      F77NAME(zgemv)('C',m,n,complex_double_one_,A_original->addr(),m,
        r.addr(),1,complex_double_zero_,c.addr(),1);
//    solve R'(-dv) =-dc:
      F77NAME(ztrsv)('U','C','N',n,QR->addr(),m,c.addr(),1);
//    solve Rdx=-v:
      F77NAME(ztrsv)('U','N','N',n,QR->addr(),m,c.addr(),1);
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
      Vector<double,complex<double> > r(n);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    r=b-A'x
      F77NAME(zgemv)('C',m,n,complex_double_mone_,A_original->addr(),m,
        x.addr(),1,complex_double_one_,r.addr(),1);
//    solve R'dv=r:
      Vector<double,complex<double> > c(m,complex_double_zero_);
      c.copyFrom(n,r);
      F77NAME(ztrsv)('U','C','N',n,QR->addr(),m,c.addr(),1);
//    dx=Q[dv\\0]
      F77NAME(zunmqr)('L','C',m,1,n,QR->addr(),m,tau->addr(),
                      c.addr(),m,&zw,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zunmqr)('L','C',m,1,n,QR->addr(),m,tau->addr(),
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
      Vector<double,complex<double> > v(n,complex_double_zero_);
      v.copyFrom(m,b);
      if (ibscl!=0) {
        F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,v.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    r=b-Ax
      F77NAME(zgemv)('N',m,n,complex_double_mone_,A_original->addr(),m,
        x.addr(),1,complex_double_one_,v.addr(),1);
//    solve Ldv=r:
      F77NAME(ztrsv)('L','N','N',m,QR->addr(),m,v.addr(),1);
//    dx=A'[dv\\dg]:
      F77NAME(zunmlq)('L','C',n,1,m,QR->addr(),m,tau->addr(),
                      v.addr(),n,&zw,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zunmlq)('L','C',n,1,m,QR->addr(),m,tau->addr(),
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
      Vector<double,complex<double> > r(n);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    r=b-A'x
      F77NAME(zgemv)('C',m,n,complex_double_mone_,A_original->addr(),m,
        x.addr(),1,complex_double_one_,r.addr(),1);
      Vector<double,complex<double> > c(m,complex_double_zero_);
//    -dc=Ar
      F77NAME(zgemv)('N',m,n,complex_double_one_,A_original->addr(),m,
        r.addr(),1,complex_double_zero_,c.addr(),1);
//    solve L(-dv)=-dc:
      F77NAME(ztrsv)('L','N','N',m,QR->addr(),m,c.addr(),1);
//    solve L'dx=-dv:
      F77NAME(ztrsv)('L','C','N',m,QR->addr(),m,c.addr(),1);
      x+=c;
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(zlascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(zlascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
}

template<> void
HouseholderQRFactorization<double,complex<double> >
::improve(const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,Factorization::TRANSPOSE_OPTION tr)
const {
//TRACER_CALL(t,"HouseholderQRFactorization::improve");
  CHECK_TEST(piv_op==Factorization::NO_PIVOTING)
  int m=QR->size(0),n=QR->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  double dw=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);
  
  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(zlange)('M',brow,k,B.addr(),B.size(0),&dw);
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
    F77NAME(zlascl)('G',0,0,bnrm,bscl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (iascl!=0) {
    F77NAME(zlascl)('G',0,0,ascl,anrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  complex<double> zw=complex_double_undefined_;
  if (m>=n) {
    if (tr==Factorization::NO_TRANSPOSE) { // min || B - A x ||
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<double,complex<double> > R(m,k);
      R.copy(B);
//    R=B-AX
      F77NAME(zgemm)('N','N',m,k,n,complex_double_mone_,
        A_original->addr(),m,X.addr(),n,complex_double_one_,R.addr(),m);
      Matrix<double,complex<double> > C(n,k,complex_double_zero_);
//    -dC=A'R
      F77NAME(zgemm)('C','N',n,k,m,complex_double_one_,
        A_original->addr(),m,R.addr(),m,complex_double_zero_,C.addr(),n);
//    solve R'(-dV)=-dC:
      F77NAME(ztrsm)('L','U','C','N',n,k,complex_double_one_,QR->addr(),m,
        C.addr(),n);
//    solve RdX=-V:
      F77NAME(ztrsm)('L','U','N','N',n,k,complex_double_one_,QR->addr(),m,
        C.addr(),n);
      X+=C;
      scllen=n;
    } else { // min || X || s.t. A^t X = B
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<double,complex<double> > R(n,k);
      R.copy(B);
//    R=B-A'X
      F77NAME(zgemm)('C','N',n,k,m,complex_double_mone_,
        A_original->addr(),m,X.addr(),m,complex_double_one_,R.addr(),n);
//    solve R'dV=R:
      Matrix<double,complex<double> > C(m,k,complex_double_zero_);
      C.copyFrom('A',n,k,R);
      F77NAME(ztrsm)('L','U','C','N',n,k,complex_double_one_,QR->addr(),m,
        C.addr(),n);
//    dX=Q[dV\\0]
      F77NAME(zunmqr)('L','C',m,k,n,QR->addr(),m,tau->addr(),
                      C.addr(),m,&zw,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zunmqr)('L','C',m,k,n,QR->addr(),m,tau->addr(),
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
      Matrix<double,complex<double> > V(n,k,complex_double_zero_);
      V.copyFrom('A',m,k,B);
//    R=B-AX
      F77NAME(zgemm)('N','N',m,k,n,complex_double_mone_,
        A_original->addr(),m,X.addr(),n,complex_double_one_,V.addr(),m);
//    solve LdV=R:
      F77NAME(ztrsm)('L','L','N','N',m,k,complex_double_one_,QR->addr(),m,
        V.addr(),n);
//    dX=Q'[dV\\0]:
      F77NAME(zunmlq)('L','C',n,k,m,QR->addr(),m,tau->addr(),
                      V.addr(),n,&zw,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(zw.real());
      complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
      F77NAME(zunmlq)('L','C',n,k,m,QR->addr(),m,tau->addr(),
                      V.addr(),n,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      X+=V;
      scllen=n;
    } else { // min || B - A' X ||
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<double,complex<double> > R(n,k);
      R.copy(B);
//    R=B-A'X
      F77NAME(zgemm)('C','N',n,k,m,complex_double_mone_,
        A_original->addr(),m,X.addr(),m,complex_double_one_,R.addr(),n);
      Matrix<double,complex<double> > C(m,k,complex_double_zero_);
//    (-dC)=A R
      F77NAME(zgemm)('N','N',m,k,n,complex_double_one_,
        A_original->addr(),m,R.addr(),n,complex_double_zero_,C.addr(),m);
//    L(-dV)=-dC
      F77NAME(ztrsm)('L','L','N','N',m,k,complex_double_one_,QR->addr(),m,
        C.addr(),m);
//    L'dX=-dV
      F77NAME(ztrsm)('L','L','C','N',m,k,complex_double_one_,QR->addr(),m,
        C.addr(),m);
      X+=C;
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(zlascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(zlascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
}

template class HouseholderQRFactorization<double,complex<double> >;
//template void testHouseholderQRFactorization(double,complex<double> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "CompleteOrthogonalDecomposition.C"
template<> CompleteOrthogonalDecomposition<double,complex<double> >::
CompleteOrthogonalDecomposition(const Matrix<double,complex<double> > &A,
double rcond) : the_rank(0),URV(0),utau(0),vtau(0),jpvt(0),iascl(0),
ascl(double_one_),anrm(numeric_limits<double>::infinity()),
A_original(&A) {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::CompleteOrthogonalDecomposition");
  int m=A.size(0),n=A.size(1);
  URV=OPERATOR_NEW Matrix<double,complex<double> >(m,n);
  URV->copy(A);
  double dw=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);
  anrm=F77NAME(zlange)('M',m,n,A.addr(),m,&dw);
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
    F77NAME(zlascl)('G',0,0,anrm,ascl,m,n,URV->addr(),m,info);
    CHECK_SAME(info,0);
  }

  int mn=min(m,n);
  int nb1=F77NAME(ilaenv)(1,"DGEQRF"," ",m,n,-1,1);
  int nb2=F77NAME(ilaenv)(1,"DGERQF"," ",m,n,-1,1);
  int nb=max(nb1,nb2);
  int lwkmin=mn+max(2*mn,n+1);
  int lwork=max(lwkmin,mn+2*n+nb*(n+1));
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);

  utau=OPERATOR_NEW Vector<double,complex<double> >(mn);
  jpvt=OPERATOR_NEW_BRACKET(int,n);
  double *rwork=OPERATOR_NEW_BRACKET(double,2*n);
  F77NAME(zgeqp3)(m,n,URV->addr(),m,jpvt,utau->addr(),work,lwork,
    rwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  delete [] rwork; rwork=0;
  double smax=abs((*URV)(0,0));
  if (smax<=double_zero_) return;
  double smin=smax;

  the_rank=1;
  Vector<double,complex<double> > xmin(mn);
  xmin[0]=complex_double_one_;
  Vector<double,complex<double> > xmax(mn);
  xmax[0]=complex_double_one_;
  while (the_rank<mn) {
    int i=the_rank+1;
    double sminpr=numeric_limits<double>::infinity();
    complex<double> s1=complex_double_undefined_,
      c1=complex_double_undefined_;
//                 (job,j,x,sest,w,gamma,sestpr,s,c)
    F77NAME(zlaic1)(2,the_rank,xmin.addr(),smin,URV->addr(0,i-1),
      (*URV)(i-1,i-1),sminpr,s1,c1);
    double smaxpr=numeric_limits<double>::infinity();
    complex<double> s2=complex_double_undefined_,
      c2=complex_double_undefined_;
    F77NAME(zlaic1)(1,the_rank,xmax.addr(),smax,URV->addr(0,i-1),
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
    vtau=OPERATOR_NEW Vector<double,complex<double> >(the_rank);
    complex<double> zw=complex_double_undefined_;
    int lwork=-1;
    F77NAME(ztzrzf)(the_rank,n,URV->addr(),m,vtau->addr(),&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(ztzrzf)(the_rank,n,URV->addr(),m,vtau->addr(),work,lwork,
      info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
  }
}

template<> OrthogonalMatrix<double,complex<double> >* 
CompleteOrthogonalDecomposition<double,complex<double> >
::leftOrthogonalPart() const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::leftOrthogonalPart");
  int m=URV->size(0),n=URV->size(1);
  OrthogonalMatrix<double,complex<double> > *Q=
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(m,m);
  Q->copyFrom('A',m,the_rank,*URV);

  complex<double> w=complex_double_undefined_;
  int lwork=-1;
  int info;
  F77NAME(zungqr)(m,m,the_rank,Q->addr(),m,utau->addr(),&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zungqr)(m,m,the_rank,Q->addr(),m,utau->addr(),work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  return Q;
}

template<> OrthogonalMatrix<double,complex<double> >* 
CompleteOrthogonalDecomposition<double,complex<double> >
::rightOrthogonalPartTransposed() const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::rightOrthogonalPartTransposed");
  int m=URV->size(0),n=URV->size(1);
  OrthogonalMatrix<double,complex<double> > *Q=
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(n,n);

  int lwork=-1;
  int info;
  complex<double> zw=complex_double_undefined_;
  F77NAME(zunmrz)('L','C',n,n,the_rank,n-the_rank,URV->addr(),m,
    vtau->addr(),Q->addr(),n,&zw,lwork,info);
  CHECK_SAME(info,0)
  lwork=static_cast<int>(zw.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  CHECK_SAME(info,0)
  F77NAME(zunmrz)('L','C',n,n,the_rank,n-the_rank,URV->addr(),m,
    vtau->addr(),Q->addr(),n,work,lwork,info);
  delete [] work; work=0;
  return Q;
}

template<> double
CompleteOrthogonalDecomposition<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,Factorization::TRANSPOSE_OPTION tr)
const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve");
  int m=URV->size(0),n=URV->size(1);
  double dw=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(zlange)('M',brow,1,b.addr(),b.size(),&dw);
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
  complex<double> zw=complex_double_undefined_;
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
    Vector<double,complex<double> > xtmp(m);
    xtmp.copy(b);
    if (ibscl!=0) {
      F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
      CHECK_SAME(info,0)
    }
//  U' b = [ y \\ z ]:
    F77NAME(zunmqr)('L','C',m,1,mn,URV->addr(),m,utau->addr(),
                    xtmp.addr(),m,&zw,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(zw.real());
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunmqr)('L','C',m,1,mn,URV->addr(),m,utau->addr(),
                    xtmp.addr(),m,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R v = y:
    F77NAME(ztrsv)('U','N','N',the_rank,URV->addr(),m,xtmp.addr(),1);
    Vector<double,complex<double> > px(n,complex_double_zero_);
    px.copyFrom(the_rank,xtmp);
    lwork=-1;
//  P'x=V[v\\0]:
    F77NAME(zunmrz)('L','C',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),px.addr(),n,&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunmrz)('L','C',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),px.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  permute
    for (int i=0;i<n;i++) x[jpvt[i]-1]=px[i];
    residual_norm=F77NAME(dznrm2)(m-the_rank,xtmp.addr(the_rank),1);
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
    Vector<double,complex<double> > xtmp(n);
//  P'b:
    for (int i=0;i<n;i++) xtmp[i]=b[jpvt[i]-1];
    if (ibscl!=0) {
      F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
      CHECK_SAME(info,0)
    }
    lwork=-1;
//  [y\\z]=V'(P'b):
    F77NAME(zunmrz)('L','N',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),xtmp.addr(),n,&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunmrz)('L','N',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),xtmp.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R' v = b
    F77NAME(ztrsv)('U','C','N',the_rank,URV->addr(),m,xtmp.addr(),1);

    x=complex_double_zero_;
    x.copyFrom(the_rank,xtmp);
//  U [v\\0]
    lwork=-1;
    F77NAME(zunmqr)('L','N',m,1,mn,URV->addr(),m,utau->addr(),
                    x.addr(),m,&zw,lwork,info );
    CHECK_SAME(info,0)

    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunmqr)('L','N',m,1,mn,URV->addr(),m,utau->addr(),
                    x.addr(),m,work,lwork,info );
    CHECK_SAME(info,0)
    delete [] work; work=0;
    residual_norm=F77NAME(dznrm2)(n-the_rank,xtmp.addr(the_rank),1);
  }
  if (iascl!=0) {
    F77NAME(zlascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(zlascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
  return residual_norm;
}

template<> void
CompleteOrthogonalDecomposition<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,
Vector<double,double> &residual_norm,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve");
  int m=URV->size(0),n=URV->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  CHECK_SAME(k,residual_norm.size())
  double dw=numeric_limits<double>::infinity();
  double smlnum=F77NAME(dlamch)('S')/F77NAME(dlamch)('P');
  double bignum=double_one_/smlnum;
  F77NAME(dlabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  double bnrm=F77NAME(zlange)('M',brow,k,B.addr(),B.size(0),&dw);
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
  complex<double> zw=complex_double_undefined_;
  if (tr==Factorization::NO_TRANSPOSE) { // min || B - A X || min || X ||
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
    Matrix<double,complex<double> > Xtmp(m,k);
    Xtmp.copy(B);
    if (ibscl!=0) {
      F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
      CHECK_SAME(info,0)
    }
//  U' B = [ Y \\ Z ]:
    F77NAME(zunmqr)('L','C',m,k,mn,URV->addr(),m,utau->addr(),
                    Xtmp.addr(),m,&zw,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(zw.real());
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunmqr)('L','C',m,k,mn,URV->addr(),m,utau->addr(),
                    Xtmp.addr(),m,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R V = Y:
    F77NAME(ztrsm)('L','U','N','N',the_rank,k,complex_double_one_,
      URV->addr(),m,Xtmp.addr(),m);
    Matrix<double,complex<double> > PX(n,k,complex_double_zero_);
    PX.copyFrom('A',the_rank,k,Xtmp);
    lwork=-1;
//  P'X=V[V\\0]:
    F77NAME(zunmrz)('L','C',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),PX.addr(),n,&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunmrz)('L','C',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),PX.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  permute
    for (int j=0;j<k;j++) {
      for (int i=0;i<n;i++) X(jpvt[i]-1,j)=PX(i,j);
      residual_norm[j]=
        F77NAME(dznrm2)(m-the_rank,Xtmp.addr(the_rank,j),1);
    }
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
    Matrix<double,complex<double> > Xtmp(n,k);
//  P'B:
    for (int j=0;j<k;j++) {
      for (int i=0;i<n;i++) Xtmp(i,j)=B(jpvt[i]-1,j);
    }
    if (ibscl!=0) {
      F77NAME(zlascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
      CHECK_SAME(info,0)
    }
    lwork=-1;
//  [Y\\Z]=V'(P'B):
    F77NAME(zunmrz)('L','N',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),Xtmp.addr(),n,&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunmrz)('L','N',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),Xtmp.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R' V = B
    F77NAME(ztrsm)('L','U','C','N',the_rank,k,complex_double_one_,
      URV->addr(),m,Xtmp.addr(),n);

    X=complex_double_zero_;
    X.copyFrom('A',the_rank,k,Xtmp);
//  U [V\\0]
    lwork=-1;
    F77NAME(zunmqr)('L','N',m,k,mn,URV->addr(),m,utau->addr(),
                    X.addr(),m,&zw,lwork,info );
    CHECK_SAME(info,0)

    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
    F77NAME(zunmqr)('L','N',m,k,mn,URV->addr(),m,utau->addr(),
                    X.addr(),m,work,lwork,info );
    CHECK_SAME(info,0)
    delete [] work; work=0;
    for (int j=0;j<k;j++) {
      residual_norm[j]=
        F77NAME(dznrm2)(n-the_rank,Xtmp.addr(the_rank,j),1);
    }
  }
  if (iascl!=0) {
    F77NAME(zlascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(zlascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
}

template class CompleteOrthogonalDecomposition<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "GramSchmidtQRFactorization.C"
template<> double
  GramSchmidtQRFactorization<double,complex<double> >::tol=
  double_one_/sqrt(2.);

template<> void
GramSchmidtQRFactorization<double,complex<double> >::reorthogonalize(
int j) {
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,R->size(0))
  CHECK_BOUNDS(j,0,n)
  for (int k=0;k<j;k++) {
    complex<double> dot=F77NAME(zdotc)(m,Q->addr(0,k),1,Q->addr(0,j),1)
      /(*R)(k,k);
    (*R)(k,j)+=dot;
    F77NAME(zaxpy)(m,-dot,Q->addr(0,k),1,Q->addr(0,j),1);
  }
  (*R)(j,j)=F77NAME(zdotc)(m,Q->addr(0,j),1,Q->addr(0,j),1);
}

template<> GramSchmidtQRFactorization<double,complex<double> >
::GramSchmidtQRFactorization(const Matrix<double,complex<double> > &A) :
A_original(&A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::GramSchmidtQRFactorization");
  int m=A.size(0),n=A.size(1);
  assert(m>=n);
  Q=OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(m,n);
  R=OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(n,n);
  Q->copy(A);

  double *norm=OPERATOR_NEW_BRACKET(double,n);
  for (int j=0;j<n;j++) {
    norm[j]=real(F77NAME(zdotc)(m,Q->addr(0,j),1,Q->addr(0,j),1));
  }
  double tol2=tol*tol;
  for (int k=0;k<n;k++) {
    complex<double> length=
      F77NAME(zdotc)(m,Q->addr(0,k),1,Q->addr(0,k),1);
    double rlength=real(length);
    (*R)(k,k)=rlength;
    if (rlength<tol2*norm[k]) reorthogonalize(k);
    if (rlength>double_zero_ &&k<n-1) {
      F77NAME(zgemv)('C',m,n-k-1,complex_double_one_/(*R)(k,k),
        Q->addr(0,k+1),m,Q->addr(0,k),1,complex_double_zero_,
        R->addr(k,k+1),n);
      for (int j=k+1;j<n;j++) (*R)(k,j)=conj((*R)(k,j));
      F77NAME(zgeru)(m,n-k-1,complex_double_mone_,Q->addr(0,k),1,
        R->addr(k,k+1),n,Q->addr(0,k+1),m);
    }
  }
  delete [] norm;
}

template<> void
GramSchmidtQRFactorization<double,complex<double> >::solveOverdetermined(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,
Vector<double,complex<double> > &residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,b.size());
  CHECK_SAME(n,x.size());
  CHECK_SAME(m,residual.size());
  residual.copy(b);
  for (int j=0;j<n;j++) {
    complex<double> dot=
      -F77NAME(zdotc)(m,Q->addr(0,j),1,residual.addr(),1)/(*R)(j,j);
    x[j]=-dot;
    F77NAME(zaxpy)(m,dot,Q->addr(0,j),1,residual.addr(),1);
  }
  F77NAME(ztrsv)('U','N','U',n,R->addr(),n,x.addr(),1);
}

template<> void
GramSchmidtQRFactorization<double,complex<double> >::solveOverdetermined(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,
Matrix<double,complex<double> > &Residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,X.size(0));
  CHECK_SAME(k,X.size(1));
  CHECK_SAME(m,Residual.size(0));
  CHECK_SAME(k,Residual.size(1));
  Residual.copy(B);
  for (int j=0;j<n;j++) {
    F77NAME(zgemv)('C',m,k,complex_double_one_/(*R)(j,j),
      Residual.addr(),m,Q->addr(0,j),1,
      complex_double_zero_,X.addr(j,0),n);
    for (int kk=0;kk<k;kk++) X(j,kk)=conj(X(j,kk));
    F77NAME(zgeru)(m,k,complex_double_mone_,Q->addr(0,j),1,X.addr(j,0),n,
      Residual.addr(),m);
  }
  F77NAME(ztrsm)('L','U','N','U',n,k,complex_double_one_,R->addr(),n,
    X.addr(),n);
}

template<> void
GramSchmidtQRFactorization<double,complex<double> >::solveUnderdetermined(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveUnderdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,b.size());
  CHECK_SAME(m,x.size());

  Vector<double,complex<double> > *v=
    OPERATOR_NEW Vector<double,complex<double> >(n);
  v->copy(b);

  F77NAME(ztrsv)('U','C','U',n,R->addr(),n,v->addr(),1);
  x=complex_double_zero_;
  for (int j=0;j<n;j++) {
    complex<double> omega=((*v)[j]
      -F77NAME(zdotc)(m,Q->addr(0,j),1,x.addr(),1))/(*R)(j,j);
    F77NAME(zaxpy)(m,omega,Q->addr(0,j),1,x.addr(),1);
  }
  delete v; v=0;
}

template<> void
GramSchmidtQRFactorization<double,complex<double> >::solveUnderdetermined(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveUnderdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(n,B.size(0));
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(k,X.size(1));

  Matrix<double,complex<double> > *V=
    OPERATOR_NEW Matrix<double,complex<double> >(n,k);
  V->copy(B);
  F77NAME(ztrsm)('L','U','C','U',n,k,complex_double_one_,R->addr(),n,
    V->addr(),n);
  X=complex_double_zero_;
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,k);
  for (int j=0;j<n;j++) {
    F77NAME(zgemv)('C',m,k,complex_double_one_,X.addr(),m,Q->addr(0,j),1,
      complex_double_zero_,work,1);
    for (int kk=0;kk<k;kk++) {
      work[kk]=((*V)(j,kk)-conj(work[kk]))/(*R)(j,j);
    }
    F77NAME(zgeru)(m,k,complex_double_one_,Q->addr(0,j),1,work,1,
      X.addr(),m);
  }
  delete [] work; work=0;
  delete V; V=0;
}

template<> void GramSchmidtQRFactorization<double,complex<double> >
::improveOverdetermined(const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,
Vector<double,complex<double> > &residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveOverdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,b.size());
  CHECK_SAME(n,x.size());
  CHECK_SAME(m,residual.size());
  Vector<double,complex<double> > *delta_b=
    OPERATOR_NEW Vector<double,complex<double> >(m);
  Vector<double,complex<double> > *delta_z=
    OPERATOR_NEW Vector<double,complex<double> >(n,0.);
  Vector<double,complex<double> > *delta_y=
    OPERATOR_NEW Vector<double,complex<double> >(n,0.);
  F77NAME(zgemv)('N',m,n,complex_double_mone_,A_original->addr(),m,
    x.addr(),1,complex_double_zero_,delta_b->addr(),1);
  F77NAME(zgemv)('C',m,n,complex_double_mone_,A_original->addr(),m,
    residual.addr(),1,complex_double_zero_,delta_z->addr(),1);
  for (int i=0;i<m;i++) (*delta_b)[i]+=b[i]-residual[i];
  F77NAME(ztrsv)('U','C','U',n,R->addr(),n,delta_z->addr(),1);
  for (int j=0;j<n;j++) {
    (*delta_y)[j]=(F77NAME(zdotc)(m,Q->addr(0,j),1,delta_b->addr(),1)
                  -(*delta_z)[j])/(*R)(j,j);
    F77NAME(zaxpy)(m,-(*delta_y)[j],Q->addr(0,j),1,delta_b->addr(),1);
  }
  F77NAME(ztrsv)('U','N','U',n,R->addr(),n,delta_y->addr(),1);
  x+=(*delta_y);
  residual+=(*delta_b);
  delete delta_y; delta_y=0;
  delete delta_z; delta_z=0;
  delete delta_b; delta_b=0;
}

template<> void GramSchmidtQRFactorization<double,complex<double> >
::improveOverdetermined(const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,
Matrix<double,complex<double> > &Residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,X.size(0));
  CHECK_SAME(k,X.size(1));
  CHECK_SAME(m,Residual.size(0));
  CHECK_SAME(k,Residual.size(1));
  Matrix<double,complex<double> > *delta_B=
    OPERATOR_NEW Matrix<double,complex<double> >(m,k);
  Matrix<double,complex<double> > *delta_Z=
    OPERATOR_NEW Matrix<double,complex<double> >(n,k,0.);
  Matrix<double,complex<double> > *delta_Y=
    OPERATOR_NEW Matrix<double,complex<double> >(n,k,0.);
  F77NAME(zgemm)('N','N',m,k,n,complex_double_mone_,A_original->addr(),m,
    X.addr(),n,complex_double_zero_,delta_B->addr(),m);
  F77NAME(zgemm)('C','N',n,k,m,complex_double_mone_,A_original->addr(),m,
    Residual.addr(),m,complex_double_zero_,delta_Z->addr(),n);
  for (int j=0;j<k;j++) {
    for (int i=0;i<m;i++) (*delta_B)(i,j)+=B(i,j)-Residual(i,j);
  }
  F77NAME(ztrsm)('L','U','C','U',n,k,complex_double_one_,R->addr(),n,
    delta_Z->addr(),n);
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,k);
  for (int j=0;j<n;j++) {
    F77NAME(zgemv)('C',m,k,complex_double_one_,delta_B->addr(),m,
      Q->addr(0,j),1,complex_double_zero_,work,1);
    for (int kk=0;kk<k;kk++) {
      (*delta_Y)(j,kk)=(conj(work[kk])-(*delta_Z)(j,kk))/(*R)(j,j);
    }
    F77NAME(zgeru)(m,k,complex_double_mone_,Q->addr(0,j),1,
      delta_Y->addr(j,0),n,delta_B->addr(),m);
  }
  delete [] work; work=0;
  F77NAME(ztrsm)('L','U','N','U',n,k,complex_double_one_,R->addr(),n,
    delta_Y->addr(),n);
  X+=(*delta_Y);
  Residual+=(*delta_B);
  delete delta_Y; delta_Y=0;
  delete delta_Z; delta_Z=0;
  delete delta_B; delta_B=0;
}

template<> void GramSchmidtQRFactorization<double,complex<double> >
::improveUnderdetermined(const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveUnderdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,b.size());
  CHECK_SAME(m,x.size());

  Vector<double,complex<double> > *y=
    OPERATOR_NEW Vector<double,complex<double> >(n);
  Vector<double,complex<double> > *delta_x=
    OPERATOR_NEW Vector<double,complex<double> >(m);
  Vector<double,complex<double> > *delta_by=
    OPERATOR_NEW Vector<double,complex<double> >(n);
  y->copy(b);
  delta_x->copy(x);
  delta_by->copy(b);
  F77NAME(ztrsv)('U','C','U',n,R->addr(),n,y->addr(),1);
  for (int j=0;j<n;j++) (*y)[j]/=(*R)(j,j);
  F77NAME(zgemv)('N',m,n,complex_double_one_,Q->addr(),m,y->addr(),1,
    complex_double_mone_,delta_x->addr(),1);
  F77NAME(zgemv)('C',m,n,complex_double_mone_,A_original->addr(),m,
    x.addr(),1,complex_double_one_,delta_by->addr(),1);
  F77NAME(ztrsv)('U','C','U',n,R->addr(),n,delta_by->addr(),1);
  for (int j=0;j<n;j++) {
    complex<double> omega=
      (F77NAME(zdotc)(m,Q->addr(0,j),1,delta_x->addr(),1)-(*delta_by)[j])
      /(*R)(j,j);
    F77NAME(zaxpy)(m,-omega,Q->addr(0,j),1,delta_x->addr(),1);
  }
  x+=(*delta_x);
  delete delta_by; delta_by=0;
  delete delta_x; delta_x=0;
  delete y; y=0;
}

template<> void GramSchmidtQRFactorization<double,complex<double> >
::improveUnderdetermined(const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveUnderdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(n,B.size(0));
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(k,X.size(1));

  Matrix<double,complex<double> > *Y=
    OPERATOR_NEW Matrix<double,complex<double> >(n,k);
  Matrix<double,complex<double> > *delta_X=
    OPERATOR_NEW Matrix<double,complex<double> >(m,k);
  Matrix<double,complex<double> > *delta_BY=
    OPERATOR_NEW Matrix<double,complex<double> >(n,k);
  Y->copy(B);
  delta_X->copy(X);
  delta_BY->copy(B);

  F77NAME(ztrsm)('L','U','C','U',n,k,complex_double_one_,R->addr(),n,
    Y->addr(),n);
  for (int j=0;j<n;j++) {
    F77NAME(zdscal)(k,double_one_/real((*R)(j,j)),Y->addr(j,0),n);
  }
  F77NAME(zgemm)('N','N',m,k,n,complex_double_one_,Q->addr(),m,
    Y->addr(),n,complex_double_mone_,delta_X->addr(),m);
  F77NAME(zgemm)('C','N',n,k,m,complex_double_mone_,A_original->addr(),m,
    X.addr(),m,complex_double_one_,delta_BY->addr(),n);
  F77NAME(ztrsm)('L','U','C','U',n,k,complex_double_one_,R->addr(),n,
    delta_BY->addr(),n);
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,k);
  for (int j=0;j<n;j++) {
    F77NAME(zgemv)('C',m,k,complex_double_one_,delta_X->addr(),m,
      Q->addr(0,j),1,complex_double_zero_,work,1);
    for (int kk=0;kk<k;kk++) {
      work[kk]=(conj(work[kk])-(*delta_BY)(j,kk))/(*R)(j,j);
    }
    F77NAME(zgeru)(m,k,complex_double_mone_,Q->addr(0,j),1,work,1,
      delta_X->addr(),m);
  }
  X+=(*delta_X);
  delete [] work; work=0;
  delete delta_BY; delta_BY=0;
  delete delta_X; delta_X=0;
  delete Y; Y=0;
}

template<> void
GramSchmidtQRFactorization<double,complex<double> >::dropColumn(int j) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropColumn");
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(j,0,n)
  for (int k=j+1;k<n;k++) {
    double c;
    complex<double> s;
    F77NAME(zrotg)(R->operator()(k,k),R->operator()(j,k),c,s);
    int ncols=n-k-1;
    if (ncols>0) {
      F77NAME(zrot)(ncols,R->addr(k,k+1),n,R->addr(j,k+1),n,c,s);
    }
    F77NAME(zrot)(m,Q->addr(0,k),1,Q->addr(0,j),1,c,s);
  }

  for (int k=j+1;k<n;k++) {
    for (int i=0;i<j;i++) (*R)(i,k-1)=(*R)(i,k);
    for (int i=j;i<k;i++) (*R)(i,k-1)=(*R)(i+1,k);
    F77NAME(zcopy)(m,Q->addr(0,k),1,Q->addr(0,k-1),1);
  }

  OrthogonalMatrix<double,complex<double> > *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(m,n-1);
  new_Q->copyFrom('A',m,n-1,*Q);
  delete Q;
  Q=new_Q;
  UpperTrapezoidalMatrix<double,complex<double> > *new_R=OPERATOR_NEW
    UpperTrapezoidalMatrix<double,complex<double> >(n-1,n-1);
  new_R->copyFrom(n-1,n-1,*R);
  delete R;
  R=new_R;
}

template<> void
GramSchmidtQRFactorization<double,complex<double> >::addColumn(int j,
const Matrix<double,complex<double> > &A) {
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,A.size(0))

  OrthogonalMatrix<double,complex<double> > *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(m,n+1);
  new_Q->copyFrom('A',m,n,*Q);
  delete Q;
  Q=new_Q;
  F77NAME(zcopy)(m,A.addr(0,j),1,Q->addr(0,n),1);

  UpperTrapezoidalMatrix<double,complex<double> > *new_R=
    OPERATOR_NEW UpperTrapezoidalMatrix<double,complex<double> >(n+1,n+1);
  new_R->copyFrom(n,n,*R);
  delete R;
  R=new_R;

  double norm=F77NAME(dznrm2)(m,Q->addr(0,n),1);
  for (j=0;j<n;j++) {
    complex<double>  dot=
      -F77NAME(zdotc)(m,Q->addr(0,j),1,Q->addr(0,n),1);
    (*R)(j,n)=-dot;
    F77NAME(zaxpy)(m,dot,Q->addr(0,j),1,Q->addr(0,n),1);
  }
  double nrm=F77NAME(dznrm2)(m,Q->addr(0,n),1);
  (*R)(n,n)=nrm;
  if (nrm<tol*norm) reorthogonalize(n);
  complex<double> a=complex_double_one_/(*R)(n,n);
  F77NAME(zscal)(m,a,Q->addr(0,n),1);
}

template<> void
GramSchmidtQRFactorization<double,complex<double> >::dropRow(int i) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropRow");
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(i,0,m)
  for (int k=i+1;k<m;k++) {
    F77NAME(zswap)(n,Q->addr(k-1,0),m,Q->addr(k,0),m);
  }

  complex<double> *work_column=OPERATOR_NEW_BRACKET(complex<double>,m);
  for (int i=0;i<m;i++) work_column[i]=complex_double_zero_;
  int Mmin1=m-1;
  for (int j=0;j<n;j++) {
    F77NAME(zaxpy)(Mmin1,conj((*Q)(m-1,j)),Q->addr(0,j),1,work_column,1);
  }
  F77NAME(zdscal)(Mmin1,double_mone_,work_column,1);

  double nrm=sqrt(double_one_-pow(F77NAME(dznrm2)(n,Q->addr(m-1,0),m),2));
  work_column[m-1]=nrm;
  double alpha=double_one_/nrm;
  F77NAME(zdscal)(Mmin1,alpha,work_column,1);

  complex<double> *work_row=OPERATOR_NEW_BRACKET(complex<double>,n);
  for (int j=0;j<n;j++) work_row[j]=complex_double_zero_;
  for (int j=n-1;j>=0;j--) {
    double c;
    complex<double> s;
    F77NAME(zrotg)(work_column[Mmin1],Q->operator()(Mmin1,j),c,s);
    F77NAME(zrot)(Mmin1,work_column,1,Q->addr(0,j),1,c,s);
    int ncols=n-j;
    F77NAME(zrot)(ncols,work_row+j,1,R->addr(j,j),n,c,s);
  }

  OrthogonalMatrix<double,complex<double> > *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(Mmin1,n);
  new_Q->copyFrom('A',Mmin1,n,*Q);
  delete Q;
  Q=new_Q;

  delete [] work_row;
  delete [] work_column;
}

template<> void
GramSchmidtQRFactorization<double,complex<double> >::addRow(int i,
const Matrix<double,complex<double> > &A) {
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,A.size(1))

  OrthogonalMatrix<double,complex<double> > *new_Q=
    OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(m+1,n);
  new_Q->copyFrom('A',m,n,*Q);
  delete Q;
  Q=new_Q;
  for (int j=0;j<n;j++) (*Q)(m,j)=complex_double_zero_;

  complex<double> *work_row=OPERATOR_NEW_BRACKET(complex<double>,n);
  int incA=A.size(0);
  F77NAME(zcopy)(n,A.addr(i,0),incA,work_row,1);

  complex<double> *work_column=OPERATOR_NEW_BRACKET(complex<double>,m+1);
  for (int k=0;k<m;k++) work_column[k]=complex_double_zero_;
  work_column[m]=complex_double_one_;

  int Mp1=m+1;
  for (int j=0;j<n;j++) {
    double c;
    complex<double> s;
    F77NAME(zrotg)(R->operator()(j,j),work_row[j],c,s);
    int ncols=n-j-1;
    if (ncols>0) {
      F77NAME(zrot)(ncols,R->addr(j,j+1),n,work_row+j+1,1,c,s);
    }
    F77NAME(zrot)(Mp1,Q->addr(0,j),1,work_column,1,c,s);
  }
  delete [] work_column;
  delete [] work_row;
}

template class GramSchmidtQRFactorization<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SingularValueDecomposition.C"
template<> SingularValueDecomposition<double,complex<double> >::
SingularValueDecomposition(const Matrix<double,complex<double> > &A) :
U(0),Vtranspose(0),s(0),A_original(&A) {
//TRACER_CALL(t,"SingularValueDecomposition::SingularValueDecomposition");
  int m=A.size(0),n=A.size(1);
  int mn=min(m,n);
  U=OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(m,mn);
  Vtranspose=OPERATOR_NEW OrthogonalMatrix<double,complex<double> >(mn,n);
  s=OPERATOR_NEW Vector<double,double>(mn);
  Matrix<double,complex<double> > Acopy(m,n);
  Acopy.copy(A);

  int info=0;
  complex<double> w=complex_double_undefined_;
  int lwork=-1;
  double *rwork=OPERATOR_NEW_BRACKET(double,5*mn);
  F77NAME(zgesvd)('S','S',m,n,Acopy.addr(),m,s->addr(),U->addr(),m,
    Vtranspose->addr(),mn,&w,lwork,rwork,info);
  lwork=static_cast<int>(w.real());
  complex<double> *work=OPERATOR_NEW_BRACKET(complex<double>,lwork);
  F77NAME(zgesvd)('S','S',m,n,Acopy.addr(),m,s->addr(),U->addr(),m,
    Vtranspose->addr(),mn,work,lwork,rwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;
  delete [] rwork; rwork=0;
}

template<> void
SingularValueDecomposition<double,complex<double> >::solve(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double rcond,
Factorization::TRANSPOSE_OPTION to) const {
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
  complex<double> *y=OPERATOR_NEW_BRACKET(complex<double>,rank);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(zgemv)('C',m,rank,double_one_,U->addr(),m,b.addr(),1,
      double_zero_,y,1);
    for (int i=0;i<rank;i++) y[i]/=(*s)[i];
    F77NAME(zgemv)('C',rank,n,double_one_,Vtranspose->addr(),mn,y,1,
      double_zero_,x.addr(),1);
  } else {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(zgemv)('N',rank,n,double_one_,Vtranspose->addr(),mn,
      b.addr(),1,double_zero_,y,1);
    for (int i=0;i<rank;i++) y[i]/=(*s)[i];
    F77NAME(zgemv)('N',m,rank,double_one_,U->addr(),m,y,1,double_zero_,
      x.addr(),1);
  }
  delete [] y; y=0;
}

template<> void
SingularValueDecomposition<double,complex<double> >::regularize(
const Vector<double,complex<double> > &b,
Vector<double,complex<double> > &x,double ridge,
Factorization::TRANSPOSE_OPTION to) const {
  CHECK_TEST(ridge>=double_zero_);
  int m=A_original->size(0),n=A_original->size(1);
  int mn=min(m,n);
  double sr=sqrt(ridge);
  complex<double> *y=OPERATOR_NEW_BRACKET(complex<double>,mn);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(zgemv)('C',m,mn,double_one_,U->addr(),m,b.addr(),1,
      double_zero_,y,1);
    for (int i=0;i<mn;i++) {
      double si=(*s)[i];
      double scale=max(si,sr);
      si/=scale;
      double ri=sr/scale;
      y[i]*=si/(scale*(si*si+ri*ri));
    }
    F77NAME(zgemv)('C',mn,n,double_one_,Vtranspose->addr(),mn,y,1,
      double_zero_,x.addr(),1);
  } else {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(zgemv)('N',mn,n,double_one_,Vtranspose->addr(),mn,b.addr(),1,
      double_zero_,y,1);
    for (int i=0;i<mn;i++) {
      double si=(*s)[i];
      double scale=max(abs(si),sr);
      si/=scale;
      double ri=sr/scale;
      y[i]*=si/(scale*(si*si+ri*ri));
    }
    F77NAME(zgemv)('N',m,mn,double_one_,U->addr(),m,y,1,double_zero_,
      x.addr(),1);
  }
  delete [] y; y=0;
}

template<> void
SingularValueDecomposition<double,complex<double> >::solve(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,double rcond,
Factorization::TRANSPOSE_OPTION to) const {
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
  complex<double> *Y=OPERATOR_NEW_BRACKET(complex<double>,rank*k);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,X.size(0));
    F77NAME(zgemm)('C','N',rank,k,m,double_one_,U->addr(),m,B.addr(),m,
      double_zero_,Y,rank);
    for (int i=0;i<rank;i++) F77NAME(zdrscl)(k,(*s)[i],Y+i,rank);
    F77NAME(zgemm)('C','N',n,k,rank,double_one_,Vtranspose->addr(),mn,
      Y,rank,double_zero_,X.addr(),n);
  } else {
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(m,X.size(0));
    F77NAME(zgemm)('N','N',rank,k,n,double_one_,Vtranspose->addr(),mn,
      B.addr(),n,double_zero_,Y,rank);
    for (int i=0;i<rank;i++) F77NAME(zdrscl)(k,(*s)[i],Y+i,rank);
    F77NAME(zgemm)('N','N',m,k,rank,double_one_,U->addr(),m,Y,rank,
      double_zero_,X.addr(),m);
  }
  delete [] Y; Y=0;
}

template<> void
SingularValueDecomposition<double,complex<double> >::regularize(
const Matrix<double,complex<double> > &B,
Matrix<double,complex<double> > &X,double ridge,
Factorization::TRANSPOSE_OPTION to) const {
  CHECK_TEST(ridge>=double_zero_);
  int m=A_original->size(0),n=A_original->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1));
  int mn=min(m,n);
  double sr=sqrt(ridge);
  complex<double> *Y=OPERATOR_NEW_BRACKET(complex<double>,mn*k);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,X.size(0));
    F77NAME(zgemm)('C','N',mn,k,m,double_one_,U->addr(),m,B.addr(),m,
      double_zero_,Y,mn);
    for (int i=0;i<mn;i++) {
      double si=(*s)[i];
      double scale=max(abs(si),sr);
      si/=scale;
      double ri=sr/scale;
      F77NAME(zdscal)(k,si/(scale*(si*si+ri*ri)),Y+i,mn);
    }
    F77NAME(zgemm)('C','N',n,k,mn,double_one_,Vtranspose->addr(),mn,Y,mn,
      double_zero_,X.addr(),n);
  } else {
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(m,X.size(0));
    F77NAME(zgemm)('N','N',mn,k,n,double_one_,Vtranspose->addr(),mn,
      B.addr(),n,double_zero_,Y,mn);
    for (int i=0;i<mn;i++) {
      double si=(*s)[i];
      double scale=max(abs(si),sr);
      si/=scale;
      double ri=sr/scale;
      F77NAME(zdscal)(k,si/(scale*(si*si+ri*ri)),Y+i,mn);
    }
    F77NAME(zgemm)('N','N',m,k,mn,double_one_,U->addr(),m,Y,mn,
      double_zero_,X.addr(),m);
  }
  delete [] Y; Y=0;
}

template class SingularValueDecomposition<double,complex<double> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
#ifdef __GNUC__
  template class CholeskyFactorization<double>;
  template ostream& operator<<(ostream&,
                               const CholeskyFactorization<double>&);
  template void testCholeskyFactorization(double);
#endif

template<> CholeskyFactorization<double>::CholeskyFactorization(
const SymmetricPositiveMatrix<double> &A) : anorm(double_zero_) {
//TRACER_CALL(t,"CholeskyFactorization::CF");
  int n=A.size(0);
  int inc=1;
  for (int j=0;j<n;j++) {
    int k=n-j;
    double col_sum=F77NAME(dasum)(&k,A.addr(j,j),&inc);
    col_sum+=F77NAME(dasum)(&j,A.addr(j,0),&n);
    if (col_sum>anorm) anorm=col_sum;
  }
  SymmetricPositiveMatrix<double> *Acopy=
    OPERATOR_NEW SymmetricPositiveMatrix<double>(A);
  int info;
  char uplo='L';
  F77NAME(dpotrf)(&uplo,&n,Acopy->addr(),&n,&info);
  CHECK_SAME(int,info,0)
  L=OPERATOR_NEW LowerTrapezoidalMatrix<double>(n,n);
  for (int j=0;j<n;j++) {
    for (int i=j;i<n;i++) (*L)(i,j)=(*Acopy)(i,j);
  }
  delete Acopy;
}

template<> Matrix<double>* CholeskyFactorization<double>::solveAXeqB(
const Matrix<double> &B) const {
//TRACER_CALL(t,"CholeskyFactorization::solveAXeqB");
  int n=L->size(0); 
  CHECK_SAME(int,n,B.size(0))
  int k=B.size(1);

  Matrix<double> *X=OPERATOR_NEW Matrix<double>(B);
  char uplo='L';
  int info;
  F77NAME(dpotrs)(&uplo,&n,&k,L->addr(),&n,X->addr(),&n,&info);
  CHECK_SAME(int,info,0)
  return X;
}

template<> double CholeskyFactorization<double>::conditionNumber() const {
  int n=L->size(0);
  double rcond=double_zero_;
  double *work=OPERATOR_NEW_BRACKET(double,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  char uplo='L';
  F77NAME(dpocon)(&uplo,&n,L->addr(),&n,&anorm,&rcond,work,iwork,
                  &info);
  CHECK_SAME(int,info,0)
  delete [] iwork;
  delete [] work;
  if (rcond>double_zero_) return double_one_/rcond;
  return numeric_limits<double>::infinity();
}
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
#ifdef __GNUC__
  template class MDMtFactorization<double>;
  template ostream& operator<<(ostream&,
                               const MDMtFactorization<double>&);
  template void testMDMtFactorization(double);
#endif

template<> MDMtFactorization<double>::MDMtFactorization(
const SymmetricMatrix<double> &A) : anorm(double_zero_) {
//TRACER_CALL(t,"MDMtFactorization::MDMtF");
  int n=A.size(0);
  ipiv=OPERATOR_NEW_BRACKET(int,n);
  int inc=1;
  for (int j=0;j<n;j++) {
    int k=n-j;
    double col_sum=F77NAME(dasum)(&k,A.addr(j,j),&inc);
    col_sum+=F77NAME(dasum)(&j,A.addr(j,0),&n);
    if (col_sum>anorm) anorm=col_sum;
  }
  SymmetricMatrix<double> *Acopy=OPERATOR_NEW SymmetricMatrix<double>(A);
  int info;
  char uplo='L';
  int nb=LaEnvBlockSize("DSYTRF",A);
  int lwork=n*nb;
  double *work=OPERATOR_NEW_BRACKET(double,lwork);
  F77NAME(dsytrf)(&uplo,&n,Acopy->addr(),&n,ipiv->addr(),work,&lwork,
                  &info);
  CHECK_SAME(int,info,0)
  delete [] work;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<double>(n,n);
  for (int j=0;j<n;j++) {
    for (int i=j;i<n;i++) (*L)(i,j)=(*Acopy)(i,j);
  }
  delete Acopy;
}

template<> Matrix<double>* MDMtFactorization<double>::solveAXeqB(
const Matrix<double> &B) const {
//TRACER_CALL(t,"MDMtFactorization::solveAXeqB");
  int n=L->size(0);
  CHECK_SAME(int,n,B.size(0))
  int k=B.size(1);

  Matrix<double> *X=OPERATOR_NEW Matrix<double>(B);
  char uplo='L';
  int info;
  F77NAME(dsytrs)(&uplo,&n,&k,L->addr(),&n,ipiv->addr(),X->addr(),&n,
                  &info);
  CHECK_SAME(int,info,0)
  return X;
}

template<> double MDMtFactorization<double>::conditionNumber() const {
  int n=L->size(0);
  double rcond=double_zero_;
  double *work=OPERATOR_NEW_BRACKET(double,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  char uplo='L';
  F77NAME(dsycon)(&uplo,&n,L->addr(),&n,ipiv->addr(),&anorm,&rcond,work,
                  iwork,&info);
  CHECK_SAME(int,info,0)
  delete [] iwork;
  delete [] work;
  if (rcond>double_zero_) return double_one_/rcond;
  return numeric_limits<double>::infinity();
}
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
#ifdef __GNUC__
  template class SingularValueDecomposition<double>;
  template ostream& operator<<(ostream&,const SingularValueDecomposition<double>&);
  template void testSingularValueDecomposition(double);
#endif

template<> SingularValueDecomposition<double>::SingularValueDecomposition(
const Matrix<double> &A) {
//TRACER_CALL(t,"SingularValueDecomposition::SVD");
  int m=A.size(0),n=A.size(1),minmn=min(m,n);
  Sigma=OPERATOR_NEW Vector<double>(minmn);
  int info=0;
  double *work=0;
  int lwork=5*max(m,n);
  if (m>=n) {
    U=OPERATOR_NEW Matrix<double>(A);
    V_transpose=OPERATOR_NEW Matrix<double>(n,n);
    char jobu='O';
    char jobvt='S';
//  int lwork=max(1,max( 3*n + m, 5*n - 4));
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dgesvd)(&jobu,&jobvt,&m,&n,U->addr(),&m,Sigma->addr(),
      0,&m,V_transpose->addr(),&n,work,&lwork,&info);
  } else {
    U=OPERATOR_NEW Matrix<double>(m,m);
    V_transpose=OPERATOR_NEW Matrix<double>(A);
    char jobu='S';
    char jobvt='O';
//  int lwork=max(1,max( 3*m + n, 5*m - 4));
    work=OPERATOR_NEW_BRACKET(double,lwork);
    F77NAME(dgesvd)(&jobu,&jobvt,&m,&n,V_transpose->addr(),&m,
      Sigma->addr(),U->addr(),&m,0,&m,work,&lwork,&info);
  }
  CHECK_SAME(int,info,0)
  if (work) delete work;
}

template<> double SingularValueDecomposition<double>::conditionNumber(
double cutoff) const {
  CHECK_POSITIVE(cutoff)
  int r=rank(cutoff);
  if (r<1) return numeric_limits<double>::infinity();
  if (r==1) return double_one_;
  return (*Sigma)[0]/(*Sigma)[r-1];
}

template<> Matrix<double>* SingularValueDecomposition<double>::solveAXeqB(
const Matrix<double> &B,double cutoff,transpose_option option) const {
//TRACER_CALL(t,"SingularValueDecomposition::solveAXeqB");
  int m=U->size(0),n=V_transpose->size(1),minmn=Sigma->size();
  int K=B.size(1);
  int r=rank(cutoff);
  double alpha=double_one_;
  double beta=double_zero_;
  char trans='N';
  switch (option) {
    case no_matrix_transpose: {
      CHECK_SAME(int,m,B.size(0))
      Matrix<double> C(r,K);
      char transa='T';
      F77NAME(dgemm)(&transa,&trans,&r,&K,&m,&alpha,U->addr(),&m,
        B.addr(),&m,&beta,C.addr(),&r);
      for (int i=0;i<r;i++) {
        double temp=1./(*Sigma)[i];
        F77NAME(dscal)(&K,&temp,C.addr(i,0),&r);
      }
      Matrix<double> *X=OPERATOR_NEW Matrix<double>(n,K);
      F77NAME(dgemm)(&transa,&trans,&n,&K,&r,&alpha,
        V_transpose->addr(),&minmn,C.addr(),&r,&beta,X->addr(),&n);
      return X;
      break;
    }
    case matrix_transpose: {
      CHECK_SAME(int,n,B.size(0))
      Matrix<double> C(r,K);
      F77NAME(dgemm)(&trans,&trans,&r,&K,&n,&alpha,V_transpose->addr(),
        &n,B.addr(),&n,&beta,C.addr(),&r);
      for (int i=0;i<r;i++) {
        double temp=1./(*Sigma)[i];
        F77NAME(dscal)(&K,&temp,C.addr(i,0),&r);
      }
      Matrix<double> *X=OPERATOR_NEW Matrix<double>(m,K);
      F77NAME(dgemm)(&trans,&trans,&m,&K,&r,&alpha,U->addr(),&m,
        C.addr(),&r,&beta,X->addr(),&m);
      return X;
      break;
    }
    default:
      return 0;
  }
}
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//#include "LinearProgram.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
#ifdef __GNUC__
  template class VirtualLinearProgram<double>;
  template ostream& operator<<(ostream&,
                               const VirtualLinearProgram<double>&);
#endif
template<> double VirtualLinearProgram<double>::zero_ = double_zero_;
template<> double VirtualLinearProgram<double>::one_ = double_one_;
template<> double VirtualLinearProgram<double>::huge_ = numeric_limits<double>::infinity();

template<> double VirtualLinearProgram<double>::currentValue() const {
  switch (current_status) {
    case unknown:
    case infeasible:
    case unbounded:
      return huge_;
    default:
      break;
  }
  int m=A.size(0),n=A.size(1);
  if (current_status==primal_feasible || n<m) {
    int incc=c_trans.size(0);
    int incx=1;
    return F77NAME(ddot)(&n,c_trans.addr(),&incc,x->addr(),&incx);
  } else {
    int incy=y_trans->size(0);
    int incb=1;
    return F77NAME(ddot)(&m,y_trans->addr(),&incy,b.addr(),&incb);
  }
}
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
#ifdef __GNUC__
  template class SFLinearProgram<double>;
  template ostream& operator<<(ostream&,const SFLinearProgram<double>&);
#endif

template<> SFLinearProgram<double>::SFLinearProgram(
const Matrix<double> &Ain,const Matrix<double> &bin,
const Matrix<double> &cin) : VirtualLinearProgram<double>(Ain,bin,cin) {
  int m=A.size(0),n=A.size(1);
  assert(m<=n);

  int incb=1;
  int imin=F77NAME(idmin)(&m,b.addr(),&incb)-1;
//if (b(imin,0)<zero_) throw bad_constraint_vector;
  assert(b(imin,0)>=zero_);

  r_trans = OPERATOR_NEW Matrix<double>(1,n-m);
  h=OPERATOR_NEW Matrix<double>(m,1);
  column=OPERATOR_NEW_BRACKET(int,n);
  for (int j=0;j<n;j++) (*column)[j]=j;
  new_basic=n;
}

template<> status_option 
SFLinearProgram<double> ::findBasicFeasibleGuess() {
//TRACER_CALL(t,"SFLinearProgram::findBasicFeasibleGuess");
//if we have already dinked with this program, get out now
  switch (current_status) {
    case primal_feasible:
    case optimal:
      return current_status;
    case dual_feasible:
    case infeasible:
    case unbounded:
      assert(0);
    default:
      break;
  }

//see if we can luck out with the first m unknowns
  int m=A.size(0),n=A.size(1),incA=1;
  SquareMatrix<double> *active=OPERATOR_NEW SquareMatrix<double>(m,m);
  A.copyInto(*active);
  basic_inverse=inverse(*active);
  delete active;

  Matrix<double> *x_basic = productOf(*basic_inverse,b);
  int incx=1;
  int imin=F77NAME(idmin)(&m,x_basic->addr(),&incx)-1;
  double xmin=(*x_basic)(imin,0);
  x->fillWith(double_zero_);
  x_basic->copyInto(*x);
  delete x_basic;

  if (xmin<double_zero_) {
//  augment the program and find a basic feasible solution
    Matrix<double> *augmented_A=
      OPERATOR_NEW Matrix<double>(m,m+n,double(double_zero_));
    Matrix<double> *augmented_c_trans=
      OPERATOR_NEW Matrix<double>(1,m+n,double(double_zero_));
    int incaA=1;
    for (int j=0;j<m;j++) (*augmented_A)(j,j)=double_one_;
    { for (int j=0;j<n;j++) {
      F77NAME(dcopy)(&m,A.addr(0,j),&incA,augmented_A->addr(0,j+m),
		     &incaA);
      (*augmented_c_trans)(0,j+m)=c_trans(0,j);
    } }

//  initialize and solve the augmented program
    SFLinearProgram<double> LP(*augmented_A,b,*augmented_c_trans);
    LP.basic_inverse=OPERATOR_NEW SquareMatrix<double>(m,m);
    augmented_A->copyInto(*LP.basic_inverse);
    LP.x->fillWith(double_zero_);
    b.copyInto(*(LP.x));
    LP.y_trans->fillWith(double_zero_);
    c_trans.copyInto(*(LP.r_trans));
    int incr=LP.r_trans->size(0);
    LP.new_basic=F77NAME(idmin)(&n,r_trans->addr(),&incr)-1;
    LP.current_status=primal_feasible;
    while (LP.simplexStep() == primal_feasible) {;}

//  if no augmented variables are basic, store basic variables here
//  we could have trouble if an augmented variable were basic and zero
    x->fillWith(double_zero_);
    { for (int j=0;j<m;j++) {
      int lp_colj=(*(LP.column))[j];
      if (lp_colj<m) return current_status=infeasible;
      int colj=lp_colj-m;
      (*column)[j]=colj;
      (*x)(colj,0)=(*(LP.x))(lp_colj,0);
    } }
//  store order of non-basic variables
    int i=m;
    { for (int j=m;j<m+n;j++) {
      int lp_colj=(*(LP.column))[j];
      if (lp_colj>=m) {
	(*column)[i]=lp_colj-m;
	i=i+1;
      }
    } }
    basic_inverse->copy(*LP.basic_inverse);
    delete augmented_c_trans;
    delete augmented_A;
  }

//find dual solution
  for (int j=0;j<m;j++) {
    double yj=double_zero_;
    for (int k=0;k<m;k++) {
      yj += c_trans(0,(*column)[k]) * (*basic_inverse)(k,j);
    }
    (*y_trans)(0,j)=yj;
  }

//find cost reduction vector
  int jmin=largestCostReduction();
  current_status=
    ((*r_trans)(0,jmin)<double_zero_ ? primal_feasible : optimal);
  new_basic=jmin+m;

  return current_status;
}

template<> int SFLinearProgram<double>::largestCostReduction() {
//TRACER_CALL(t,"SFLinearProgram::largestCostReduction");
  int m=A.size(0),n=A.size(1);
  int jmin=0; double rmin=double_zero_;
  int incy=y_trans->size(0);
  int incA=1;
  for (int j=0;j<n-m;j++) {
    int colj=(*column)[j+m];
    double rj=c_trans(0,colj)
          -F77NAME(ddot)(&m,y_trans->addr(),&incy,A.addr(0,colj),&incA);
    if (rj < rmin) { jmin=j; rmin=rj; }
    (*r_trans)(0,j)=rj;
  }
  int NminM=n-m;
  int incr=r_trans->size(0);
  return F77NAME(idmin)(&NminM,r_trans->addr(),&incr)-1;
}

template<> status_option SFLinearProgram<double>::simplexStep() {
//TRACER_CALL(t,"SFLinearProgram::simplexStep");
//if we have already dinked with this program, get out now
  switch (current_status) {
    case optimal:
      return current_status;
    case infeasible:
    case unbounded:
      assert(0);
    default:
      break;
  }
  int m=A.size(0),n=A.size(1);

//find largest feasible simplex step
  int imin=n; 
  int colmin=n;
  double epsmin=double_zero_;
  char trans='n';
  double alpha=double_one_;
  double beta=double_zero_;
  int inch=1;
  int incA=1;
  int col_new_basic=(*column)[new_basic];
  F77NAME(dgemv)(&trans,&m,&m,&alpha,basic_inverse->addr(),&m,
		 A.addr(0,col_new_basic),&incA,&beta,h->addr(),&inch);
  for (int i=0;i<m;i++) {
    double hi=(*h)(i,0);
    int coli=(*column)[i];
    if (hi>double_zero_) {
      double xi=(*x)(coli,0);
      if (imin>=n || xi < epsmin*hi) { 
	imin=i; colmin=coli; epsmin=xi/hi;
      } else if (xi<=double_zero_ && coli<colmin) { 
	imin=i; colmin=coli; epsmin=double_zero_; 
      }
    }
  }
#ifdef DEBUG
//cout << "\th = A_1^inv * A_2 * e_jmin" << endl;
//cout << *h << endl;
//cout << "\timin = " << imin << endl;
//cout << "\tepsmin = " << epsmin << endl;
#endif
  if (imin>=n) return current_status=unbounded;

//update solution
  { for (int i=0;i<m;i++) {
    (*x)((*column)[i],0) -= (*h)(i,0) * epsmin;
  } }
  (*x)((*column)[imin],0)=double_zero_;
  (*x)(col_new_basic,0)=epsmin;
#ifdef DEBUG
//cout << "\tupdated x" << endl;
//cout << *x << endl;
#endif

//update dual solution
  double hmin=(*h)(imin,0);
  alpha=(*r_trans)(0,new_basic-m)/hmin;
  int incy=1;
  F77NAME(daxpy)(&m,&alpha,basic_inverse->addr(imin,0),&m,
		 y_trans->addr(),&incy);
#ifdef DEBUG
//cout << "\tupdated y_trans" << endl;
//cout << *y_trans << endl;
#endif

//update inverse of active constraint matrix
  { for (int i=0;i<m;i++) {
    if (i!=imin) {
      alpha=-(*h)(i,0)/hmin;
      F77NAME(daxpy)(&m,&alpha,basic_inverse->addr(imin,0),&m,
		     basic_inverse->addr(i,0),&m);
    }
  } }
  alpha=double_one_/hmin;
  F77NAME(dscal)(&m,&alpha,basic_inverse->addr(imin,0),&m);
#ifdef DEBUG
//cout << "\tupdated basic_inverse = " << endl;
//cout << *basic_inverse << endl;
#endif

//switch columns
  int col=(*column)[imin];
  (*column)[imin]=(*column)[new_basic];
  (*column)[new_basic]=col;
#ifdef DEBUG
//cout << "\tupdated column = " << endl;
//cout << *column << endl;
#endif

//find cost reduction vector
  int jmin=largestCostReduction();
  current_status=
    ((*r_trans)(0,jmin)<double_zero_ ? primal_feasible : optimal);
  new_basic=jmin+m;
#ifdef DEBUG
//cout << "\tupdated r_trans = " << endl;
//cout << *r_trans << endl;
//cout << "new_basic = " << (*column)[new_basic] << endl;
#endif

  return current_status;
}

template<> void SFLinearProgram<double>::costBounds(
Matrix<double> &lower_trans,Matrix<double> &upper_trans) const {
  int m=A.size(0),n=A.size(1);
  CHECK_SAME(int,n,lower_trans.size(1))
  CHECK_SAME(int,n,upper_trans.size(1))

  int inc0=1;
  int inc1=A.size(0);
//perturbation in basic variable
  for (int i=0;i<m;i++) {
    double lo=-huge_;
    double hi=huge_;
    for (int j=m;j<n;j++) {
      int colj=(*column)[j];
      double denom=F77NAME(ddot)(&m,basic_inverse->addr(i,0),&inc1,
				A.addr(0,colj),&inc0);
      double rj=(*r_trans)(0,j-m);
      if (denom<double_zero_ && denom*lo>rj) lo=rj/denom;
      else if (denom>double_zero_ && denom*hi>rj) hi=rj/denom;
    }
    int coli=(*column)[i];
    lower_trans(0,coli)=c_trans(0,coli)+lo;
    upper_trans(0,coli)=c_trans(0,coli)+hi;
  }
//perturbation in non-basic variable
  for (int j=m;j<n;j++) {
    int colj=(*column)[j];
    lower_trans(0,colj)=c_trans(0,colj)-(*r_trans)(0,j-m);
    upper_trans(0,colj)=huge_;
  }
}

template<> double SFLinearProgram<double>::costSensitivity(int j,
Matrix<double> &dual_derivative) const {
  int m=A.size(0);
  CHECK_SAME(int,m,dual_derivative.size(1))
  for (int i=0;i<m;i++) {
    if ((*column)[i]==j) {
      for (int k=0;k<m;k++) dual_derivative(0,k)=(*basic_inverse)(i,k);
      return (*x)(j,0);
    }
  }
  dual_derivative.fillWith(double_zero_);
  return double_zero_;
}

template<> void SFLinearProgram<double>::arrayBounds(int i,int j,
double &lower,double &upper) const {
//TRACER_CALL(t,"SFLinearProgram::arrayBounds");
  int m=A.size(0),n=A.size(1);
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
    if ((*column)[jj]==j) {
      double ainv_ji=(*basic_inverse)(jj,i);
      if (fabs(ainv_ji)>sfmin) {
	double bound=-double_one_/ainv_ji;
	if (ainv_ji>double_zero_) lower=bound;
	else upper=bound;
      }
      for (int k=0;k<m;k++) {
        if (k!=jj) {
	  double xk=(*x)((*column)[k],0);
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
	double term=F77NAME(ddot)(&m,basic_inverse->addr(jj,0),&m,
			   A.addr(0,(*column)[kk]),&incA);
	double denom=yi*term-ainv_ji*rk;
	if (fabs(denom)>rk*sfmin) {
	  double bound=rk/denom;
	  if (denom>double_zero_ && bound<upper) upper=bound;
	  if (denom<double_zero_ && bound>lower) lower=bound;
	}
      }
      lower += A(i,j);
      upper += A(i,j);
      return;
    }
  }
//perturbation in non-basic column of A
  { for (int jj=m;jj<n;jj++) {
    if ((*column)[jj]==j) {
      double rj=(*r_trans)(0,jj-m);
      if (fabs(yi)>rj*sfmin) {
	double bound=A(i,j)+rj/yi;
	if (yi>double_zero_ && bound<upper) upper=bound;
	if (yi<double_zero_ && bound>lower) bound=lower;
      }
    }
  } }
}

template<> void testSFLinearProgram(double scalar) {
  scalar=double_one_;
//Giapetto problem, p. 45 of Winston 
  Matrix<double> A(3,5); 
  Matrix<double> b(3,1); 
  Matrix<double> c(1,5); 
  c(0,0)=0.; c(0,1)=0.; c(0,2)=0.; c(0,3)=-3.; c(0,4)=-2.;
  A(0,0)=1.; A(0,1)=0.; A(0,2)=0.; A(0,3)= 2.; A(0,4)= 1.; b(0,0)=100.;
  A(1,0)=0.; A(1,1)=1.; A(1,2)=0.; A(1,3)= 1.; A(1,4)= 1.; b(1,0)= 80.;
  A(2,0)=0.; A(2,1)=0.; A(2,2)=1.; A(2,3)= 1.; A(2,4)= 0.; b(2,0)= 40.;

//
//Dorian problem, p. 130 of Winston 
//Matrix<double> A(4,7); 
//Matrix<double> b(4,1); 
//Matrix<double> c(1,7); 
//c(0,0)=0.; c(0,1)=0.; c(0,2)=0.; c(0,3)=0.; c(0,4)=-60.; c(0,5)=-30. ; c(0,6)=-20. ;
//A(0,0)=1.; A(0,1)=0.; A(0,2)=0.; A(0,3)=0.; A(0,4)=  8.; A(0,5)=  6. ; A(0,6)= 1. ; b(0,0)=48.;
//A(1,0)=0.; A(1,1)=1.; A(1,2)=0.; A(1,3)=0.; A(1,4)=  4.; A(1,5)=  2. ; A(1,6)= 1.5; b(1,0)=20.;
//A(2,0)=0.; A(2,1)=0.; A(2,2)=1.; A(2,3)=0.; A(2,4)=  2.; A(2,5)=  1.5; A(2,6)= 0.5; b(2,0)= 8.;
//A(3,0)=0.; A(3,1)=0.; A(3,2)=0.; A(3,3)=1.; A(3,4)=  0.; A(3,5)=  1. ; A(3,6)= 0. ; b(3,0)= 5.;

//first two unknowns not basic, feasible
//Matrix<double> A(2,5); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,5); 
//c(0,0)= 5.; c(0,1)=0.; c(0,2)=0.; c(0,3)= 3.; c(0,4)=-2.;
//A(0,0)=-6.; A(0,1)=0.; A(0,2)=1.; A(0,3)=-2.; A(0,4)= 2.; b(0,0)= 6.;
//A(1,0)=-3.; A(1,1)=1.; A(1,2)=0.; A(1,3)= 6.; A(1,4)= 3.; b(1,0)=15.;

//problem unbounded, p. 146
//Matrix<double> A(2,6); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,6); 
//c(0,0)=0.; c(0,1)=0.; c(0,2)=-36.; c(0,3)=-30.; c(0,4)= 3.; c(0,5)= 4.;
//A(0,0)=1.; A(0,1)=0.; A(0,2)=  1.; A(0,3)=  1.; A(0,4)=-1.; A(0,5)= 0.; b(0,0)= 5.;
//A(1,0)=0.; A(1,1)=1.; A(1,2)=  6.; A(1,3)=  5.; A(1,4)= 0.; A(1,5)=-1.; b(1,0)=10.;

//basic feasible solution is degenerate, p. 159
//Matrix<double> A(2,4);
//Matrix<double> b(2,1); 
//Matrix<double> c(1,4); 
//c(0,0)=0.; c(0,1)=0.; c(0,2)=-5.; c(0,3)=-2.;
//A(0,0)=1.; A(0,1)=0.; A(0,2)= 1.; A(0,3)= 1.; b(0,0)=6.;
//A(1,0)=0.; A(1,1)=1.; A(1,2)= 1.; A(1,3)=-1.; b(1,0)=0.;

//problem cycles without lexico-graphic ordering, p. 161
//Matrix<double> A(2,6);
//Matrix<double> b(2,1); 
//Matrix<double> c(1,6); 
//c(0,0)=0.; c(0,1)=0.; c(0,2)=-2.  ; c(0,3)=-3.; c(0,4)= 1.   ; c(0,5)=12.;
//A(0,0)=1.; A(0,1)=0.; A(0,2)=-2.  ; A(0,3)=-9.; A(0,4)= 1.   ; A(0,5)= 9.; b(0,0)=0.;
//A(1,0)=0.; A(1,1)=1.; A(1,2)=1./3.; A(1,3)= 1.; A(1,4)=-1./3.; A(1,5)=-2.; b(1,0)=0.;
//

  int m=A.size(0),n=A.size(1);
  SFLinearProgram<double> LP(A,b,c);
  LP.findBasicFeasibleGuess();

  while (LP.currentStatus()==primal_feasible) {
    LP.simplexStep();
    cout << "\nstatus = " << LP.currentStatus()
         << endl;
    cout << "current value = " << LP.currentValue() << endl;
    cout << "current solution = " << LP.currentPrimalSolution() << endl;
    char dummy[80];
    cout << "hit RETURN to continue" << endl;
    cin.getline(dummy,80);
  }

  cout << "\n\nafter iteration, status = " << LP.currentStatus()
       << endl;
  cout << "final value = " << LP.currentValue() << endl;
  cout << "final solution = " << LP.currentPrimalSolution() << endl;

  Matrix<double> *c_lower=OPERATOR_NEW Matrix<double>(1,n);
  Matrix<double> *c_upper=OPERATOR_NEW Matrix<double>(1,n);
  LP.costBounds(*c_lower,*c_upper);
  cout << "lower bounds on cost vector:" << endl;
  cout << *c_lower << endl;
  cout << "upper bounds on cost vector:" << endl;
  cout << *c_upper << endl;

  Matrix<double> *dual_derivative=OPERATOR_NEW Matrix<double>(1,m);
  for (int j=0;j<n;j++) {
    double dzdc=LP.costSensitivity(j,*dual_derivative);
    cout << "\nj,dzdc = " << j << " " << dzdc << endl;
    cout << "dual_derivative = " << *dual_derivative << endl;
  }

  Matrix<double> *b_lower=OPERATOR_NEW Matrix<double>(m,1);
  Matrix<double> *b_upper=OPERATOR_NEW Matrix<double>(m,1);
  LP.constraintBounds(*b_lower,*b_upper);
  cout << "lower bounds on constraint vector:" << endl;
  cout << *b_lower << endl;
  cout << "upper bounds on constraint vector:" << endl;
  cout << *b_upper << endl;

  Matrix<double> *primal_derivative=OPERATOR_NEW Matrix<double>(n,1);
  for (int i=0;i<m;i++) {
    double dzdb=LP.constraintSensitivity(i,*primal_derivative);
    cout << "\ni,dzdb = " << i << " " << dzdb << endl;
    cout << "primal_derivative = " << *primal_derivative << endl;
  }

  double A_lower,A_upper;
  cout << "lower bounds on array:" << endl;
  { for (int i=0;i<m;i++) {
    for (int j=0;j<n;j++) {
      LP.arrayBounds(i,j,A_lower,A_upper);
      cout << A_lower << " ";
    }
    cout << endl;
  } }
  cout << "upper bounds on array:" << endl;
  { for (int i=0;i<m;i++) {
    for (int j=0;j<n;j++) {
      LP.arrayBounds(i,j,A_lower,A_upper);
      cout << A_upper << " ";
    }
    cout << endl;
  } }

  { for (int i=0;i<m;i++) {
    for (int j=0;j<n;j++) {
      double dzda=
      LP.arraySensitivity(i,j,*primal_derivative,*dual_derivative);
      cout << "\ni,j,dzda = " << i << " " << j << " " << dzda <<endl;
      cout << "primal_derivative = " << endl;
      cout << *primal_derivative << endl;
      cout << "dual_derivative = " << endl;
      cout << *dual_derivative << endl;
    }
  } }
}
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
#ifdef __GNUC__
  template class LinearProgram<double>;
//template LinearProgram<double>::LinearProgram(const Matrix<double>&,const Matrix<double>&,const Matrix<double>&);
//template double LinearProgram<double>::currentValue() const;
//template status_option LinearProgram<double>::simplexStep();
//template void LinearProgram<double>::costBounds(Matrix<double>&,
//   Matrix<double>&) const;
//template void LinearProgram<double>::constraintBounds(Matrix<double>&,
//   Matrix<double>&) const;
//template void LinearProgram<double>::printOn(ostream&) const;
  template ostream& operator<<(ostream&,const LinearProgram<double>&);
#endif

template<> status_option LinearProgram<double>::findBasicFeasibleGuess() {
//TRACER_CALL(t,"LinearProgram::findBasicFeasibleGuess");
//if we have already dinked with this program, get out now
  switch (current_status) {
    case dual_feasible:
    case primal_feasible:
    case optimal:
      return current_status;
    case infeasible:
    case unbounded:
      assert(0);
    default:
      break;
  }
  int m=A.size(0),n=A.size(1);
  int incb=1,incc=c_trans.size(0);
  int ib=F77NAME(idsumn)(&m,b.addr(),&incb);
  int ic=F77NAME(idsump)(&n,c_trans.addr(),&incc);
#ifdef DEBUG
  cout << "ib,ic = " << ib << " " << ic << endl;
#endif
  if (m-ib<=n-ic) return findPrimalBasicFeasibleGuess();
  else return findDualBasicFeasibleGuess();
}

template<> status_option 
LinearProgram<double>::findPrimalBasicFeasibleGuess() {
//TRACER_CALL(t,"LinearProgram::findPrimalBasicFeasibleGuess");
//if we have already dinked with this program, get out now
  switch (current_status) {
    case primal_feasible:
    case optimal:
      return current_status;
    case infeasible:
    case unbounded:
      assert(0);
    default:
      break;
  }

  basic_number=0;
//see if x=0 is feasible for primal
  int m=A.size(0),n=A.size(1);
  int incb=1;
  int ib=F77NAME(idsumn)(&m,b.addr(),&incb);
  if (ib>=m) return current_status=primal_feasible;

//introduce artificial variables for positive b's
#ifdef DEBUG
  cout << "introducing artificial variables" << endl;
#endif
  Matrix<double> *temp_A=OPERATOR_NEW Matrix<double>(m,n+m-ib,double_zero_);
  A.copyInto(*temp_A);
  Matrix<double> *temp_c_trans=OPERATOR_NEW Matrix<double>(1,n+m-ib,double_zero_);
  int pos_b=n;
  for (int i=0;i<m;i++) {
    if (b(i,0)>double_zero_) {
      (*temp_A)(i,pos_b)=double_one_;
      (*temp_c_trans)(0,pos_b)=double_one_;
      pos_b++;
    }
  }
  LinearProgram temp_LP(*temp_A,b,*temp_c_trans);
  temp_LP.basic_number=pos_b-n;
  delete temp_LP.basic_inverse;
  temp_LP.basic_inverse=OPERATOR_NEW SquareMatrix<double>(temp_LP.basic_number,0);
  basic_inverse->fillWith(double_zero_);
  for (int k=0;k<temp_LP.basic_number;k++) {
    (*temp_LP.basic_inverse)(k,k)=double_one_;
  }
  temp_LP.x->fillWith(double_zero_);
  pos_b=n;
  { for (int i=0;i<m;i++) {
    if (b(i,0)>double_zero_) (*temp_LP.x)(pos_b++,0)=b(i,0);
  } }
  temp_LP.current_status=primal_feasible;
#ifdef DEBUG
  cout << "before solving LP to find initial feasible guess" << endl;
  cout << "temp_LP:" << endl;
  cout << temp_LP << endl;
#endif
  while (temp_LP.current_status==primal_feasible) {
    temp_LP.simplexStep();
  }
#ifdef DEBUG
  cout << "after solving LP to find initial feasible guess" << endl;
  cout << "temp_LP:" << endl;
  cout << temp_LP << endl;
#endif

  current_status=primal_feasible;
  basic_number=temp_LP.basic_number;
  basic_inverse->copy(*temp_LP.basic_inverse);
  x->fillWith(double_zero_);
  int pos=0;
  for (int j=0;j<basic_number;j++) {
    int temp_col=(*temp_LP.column)[j];
    if (temp_col>=n) {
      current_status=infeasible;
      break;
    } else {
      int col=(*column)[pos]=temp_col;
      (*x)(col,0)=(*temp_LP.x)(temp_col,0);
      pos++;
    }
  }
#ifdef DEBUG
  cout << "after setting initial feasible guess" << endl;
  printOn(cout);
#endif
  delete temp_A;
  delete temp_c_trans;
  return current_status;
}

template<> status_option 
LinearProgram<double>::findDualBasicFeasibleGuess() {
//TRACER_CALL(t,"LinearProgram::findDualBasicFeasibleGuess");
//if we have already dinked with this program, get out now
  switch (current_status) {
    case dual_feasible:
    case optimal:
      return current_status;
    case infeasible:
    case unbounded:
      assert(0);
    default:
      break;
  }

  basic_number=0;
//see if y=0 is feasible for dual
  int m=A.size(0),n=A.size(1);
  int incc=c_trans.size(0);
  int ic=F77NAME(idsump)(&n,c_trans.addr(),&incc);
  if (ic>=n) return current_status=dual_feasible;

//introduce artificial variables for negative c's
#ifdef DEBUG
  cout << "introducing artificial variables" << endl;
#endif
  Matrix<double> *temp_A=OPERATOR_NEW Matrix<double>(m+n-ic,n,double_zero_);
  A.copyInto(*temp_A);
  Matrix<double> *temp_b=OPERATOR_NEW Matrix<double>(m+n-ic,1,double_zero_);
  int neg_c=m;
  for (int j=0;j<n;j++) {
    if (c_trans(0,j)<double_zero_) {
      (*temp_A)(neg_c,j)=-double_one_;
      (*temp_b)(neg_c,0)=double_one_;
      neg_c++;
    }
  }
  LinearProgram temp_LP(*temp_A,*temp_b,c_trans);
  temp_LP.basic_number=neg_c;
  delete temp_LP.basic_inverse;
  temp_LP.basic_inverse=OPERATOR_NEW SquareMatrix<double>(temp_LP.basic_number,0);
  basic_inverse->fillWith(double_zero_);
  for (int k=0;k<temp_LP.basic_number;k++) {
    (*temp_LP.basic_inverse)(k,k)=double_one_;
  }
  temp_LP.y_trans->fillWith(double_zero_);
  neg_c=m;
  { for (int j=0;j<n;j++) {
    if (c_trans(0,j)<double_zero_) (*temp_LP.y_trans)(0,neg_c++)=-c_trans(0,j);
  } }
  temp_LP.current_status=dual_feasible;
#ifdef DEBUG
  cout << "before solving LP to find initial feasible guess" << endl;
  cout << "temp_LP:" << endl;
  cout << temp_LP << endl;
#endif
  while (temp_LP.current_status==dual_feasible) {
    temp_LP.simplexStep();
  }
#ifdef DEBUG
  cout << "after solving LP to find initial feasible guess" << endl;
  cout << "temp_LP:" << endl;
  cout << temp_LP << endl;
#endif

  current_status=dual_feasible;
  basic_number=temp_LP.basic_number;
  basic_inverse->copy(*temp_LP.basic_inverse);
  y_trans->fillWith(double_zero_);
  int pos=0;
  for (int i=0;i<basic_number;i++) {
    int temp_row=(*temp_LP.row)[i];
    if (temp_row>=m) {
      current_status=infeasible;
      break;
    } else {
      int ro=(*row)[pos]=temp_row;
      (*y_trans)(0,ro)=(*temp_LP.y_trans)(0,temp_row);
      pos++;
    }
  }
#ifdef DEBUG
  cout << "after setting initial feasible guess" << endl;
  printOn(cout);
#endif
  delete temp_A;
  delete temp_b;
  return current_status;
}

template<> void LinearProgram<double>::computeSolution() {
//TRACER_CALL(t,"LinearProgram::computeSolution");
#ifdef DEBUG
//cout << "basic_number = " << basic_number << endl;
//cout << "row = \n" << *row << endl;
//cout << "column = \n" << *column << endl;
//cout << "basic_matrix :\n" << *basic_matrix << endl;
#endif
  Matrix<double> *rhs_basic=OPERATOR_NEW Matrix<double>(basic_number,1);
  Matrix<double> *soln_basic=OPERATOR_NEW Matrix<double>(basic_number,1);
  for (int i=0;i<basic_number;i++) {
    (*rhs_basic)(i,0)=b((*row)[i],0);
  }
#ifdef DEBUG
//cout << "rhs_basic = \n" << *rhs_basic << endl;
#endif
  char trans='n';
  double alpha=double_one_;
  double beta=double_zero_;
  int inc_rhs=1;
  int inc_soln=1;
  F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
    basic_inverse->addr(),&basic_number,rhs_basic->addr(),&inc_rhs,
    &beta,soln_basic->addr(),&inc_soln);
#ifdef DEBUG
//cout << "soln_basic = \n" << *soln_basic << endl;
#endif
  x->fillWith(double_zero_);
  for (int j=0;j<basic_number;j++) {
    (*x)((*column)[j],0)=(*soln_basic)(j,0);
  }
#ifdef DEBUG
//cout << "x = \n" << *x << endl;
#endif

  { for (int j=0;j<basic_number;j++) {
    (*rhs_basic)(j,0)=c_trans(0,(*column)[j]);
  } }
  trans='T';
  F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
    basic_inverse->addr(),&basic_number,rhs_basic->addr(),&inc_rhs,
    &beta,soln_basic->addr(),&inc_soln);

  y_trans->fillWith(double_zero_);
  { for (int i=0;i<basic_number;i++) {
    (*y_trans)(0,(*row)[i])=(*soln_basic)(i,0);
  } }
  delete soln_basic;
  delete rhs_basic;
#ifdef DEBUG
//cout << "y_trans = \n" << *y_trans << endl;
#endif
}

template<> int LinearProgram<double>::basicPrimalPivot(Matrix<double> &h,
double &ratio) {
//TRACER_CALL(t,"LinearProgram::basicPrimalPivot");
  int i_basic=basic_number;
  ratio=huge_;
  if (basic_number>0) {
    Matrix<double> rhs;
    rhs.copy(h);
    char trans='n';
    double alpha=double_one_;
    double beta=double_zero_;
    int inc_rhs=1;
    int inc_h=1;
    F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
      basic_inverse->addr(),&basic_number,rhs.addr(),&inc_rhs,
      &beta,h.addr(),&inc_h);
#ifdef DEBUG
    cout << "after solve, h = \n" << h << endl;
#endif
    for (int i=0;i<basic_number;i++) {
      double hi=h(i,0);
      if (hi<=double_zero_) continue;
      double xi=(*x)((*column)[i],0);
      if (i_basic>=basic_number || hi*ratio<xi) {
	i_basic=i;
	ratio=xi/hi;
      }
    }
  }
  return i_basic;
}

template<> void LinearProgram<double>::add(int i_non_basic,
int j_non_basic) {
//TRACER_CALL(t,"LinearProgram::add");
#ifdef DEBUG
  cout << "adding row " << (*row)[i_non_basic] << endl;
  cout << "adding column " << (*column)[j_non_basic] << endl;
#endif
  int rowi=(*row)[i_non_basic];
  (*row)[i_non_basic]=(*row)[basic_number];
  (*row)[basic_number]=rowi;
  int colj=(*column)[j_non_basic];
  (*column)[j_non_basic]=(*column)[basic_number];
  (*column)[basic_number]=colj;
  SquareMatrix<double> *new_basic_inverse=
    OPERATOR_NEW SquareMatrix<double>(basic_number+1,basic_number+1);
#ifdef DEBUG
  cout << "basic matrix:" << endl;
  { for (int i=0;i<=basic_number;i++) {
    for (int j=0;j<=basic_number;j++) {
      cout << A((*row)[i],(*column)[j]) << " ";
    }
    cout << endl;
  } }
  cout << "basic_inverse:\n" << *basic_inverse << endl;
#endif

  Matrix<double> *temp=OPERATOR_NEW Matrix<double>(basic_number,1);
  char trans='n';
  double alpha=double_one_;
  double beta=double_zero_;
  int inc_rhs=1;
  int inc_soln=1;
  if (basic_number>0) {
    for (int i=0;i<basic_number;i++) {
      (*temp)(i,0)=-A((*row)[i],(*column)[basic_number]);
    }
    F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
      basic_inverse->addr(),&basic_number,temp->addr(),&inc_rhs,
      &beta,new_basic_inverse->addr(0,basic_number),&inc_soln);
#ifdef DEBUG
    cout << "-basic_inverse * new column:" << endl;
    { for (int i=0;i<basic_number;i++) {
      cout << (*new_basic_inverse)(i,basic_number) << " ";
    } }
    cout << endl;
#endif
  }
  double diag=A((*row)[basic_number],(*column)[basic_number]);
  for (int k=0;k<basic_number;k++) {
    diag += A((*row)[basic_number],(*column)[k]) 
          * (*new_basic_inverse)(k,basic_number);
  }
  diag=double_one_/diag;
  (*new_basic_inverse)(basic_number,basic_number)=diag;
#ifdef DEBUG
  cout << "new_basic_inverse(basic_number,basic_number) = " << diag 
       << endl;
#endif
  if (basic_number>0) {
    for (int k=0;k<basic_number;k++) {
      (*new_basic_inverse)(k,basic_number) *= diag;
    }
#ifdef DEBUG
    cout << "last column of new_basic_inverse:" << endl;
    { for (int i=0;i<=basic_number;i++) {
      cout << (*new_basic_inverse)(i,basic_number) << " ";
    } }
    cout << endl;
#endif

    for (int j=0;j<basic_number;j++) {
      (*temp)(j,0)=-A((*row)[basic_number],(*column)[j]);
    }
    trans='T';
    inc_soln=new_basic_inverse->size(0);
    F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
      basic_inverse->addr(),&basic_number,temp->addr(),&inc_rhs,
      &beta,new_basic_inverse->addr(basic_number,0),&inc_soln);
#ifdef DEBUG
    cout << "-basic_inverse * new row:" << endl;
    { for (int i=0;i<basic_number;i++) {
      cout << (*new_basic_inverse)(basic_number,i) << " ";
    } }
    cout << endl;
#endif
    { for (int j=0;j<basic_number;j++) {
      for (int i=0;i<basic_number;i++) {
	(*new_basic_inverse)(i,j)=(*basic_inverse)(i,j)
	  +(*new_basic_inverse)(i,basic_number)
	  *(*new_basic_inverse)(basic_number,j);
      }
      (*new_basic_inverse)(basic_number,j) *= diag;
    } }
  }
  basic_number++;
  delete temp;
  delete basic_inverse;
  basic_inverse=new_basic_inverse;
#ifdef DEBUG
  cout << "new basic_inverse:\n" << *basic_inverse << endl;
#endif
}

template<> void  LinearProgram<double>::switchRows(int i_basic,
int i_non_basic) {
//TRACER_CALL(t,"LinearProgram::switchRows");
#ifdef DEBUG
  cout << "switching basic row " << (*row)[i_basic] 
       << " with non-basic row " << (*row)[i_non_basic] << endl;
#endif
//update inverse of active constraint matrix
  Matrix<double> *temp_rhs=OPERATOR_NEW Matrix<double>(1,basic_number); 
  Matrix<double> *temp_soln=OPERATOR_NEW Matrix<double>(1,basic_number); 
  for (int j=0;j<basic_number;j++) {
    (*temp_rhs)(0,j)=A((*row)[i_non_basic],(*column)[j]);
  }
  char trans='T';
  double alpha=double_one_;
  double beta=double_zero_;
  int inc_rhs=temp_rhs->size(0);
  int inc_soln=temp_soln->size(0);
  F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
    basic_inverse->addr(),&basic_number,temp_rhs->addr(),&inc_rhs,
    &beta,temp_soln->addr(),&inc_soln);
  double ti=double_one_/(*temp_soln)(0,i_basic);
  int incA=1;
  { for (int j=0;j<basic_number;j++) {
    if (j!=i_basic) {
      alpha=-(*temp_soln)(0,j)*ti;
      F77NAME(daxpy)(&basic_number,&alpha,
        basic_inverse->addr(0,i_basic),&incA,
        basic_inverse->addr(0,j),&incA);
    }
  } }
  F77NAME(dscal)(&basic_number,&ti,basic_inverse->addr(0,i_basic),
    &incA);
  delete temp_rhs;
  delete temp_soln;
#ifdef DEBUG
//cout << "\tupdated basic_inverse = " << endl;
//cout << *basic_inverse << endl;
#endif

//switch rows
  int rowi=(*row)[i_non_basic];
  (*row)[i_non_basic]=(*row)[i_basic];
  (*row)[i_basic]=rowi;
#ifdef DEBUG
//cout << "\tupdated row = " << endl;
//cout << *row << endl;
#endif
}

template<> void  LinearProgram<double>::switchColumns(int j_basic,
int j_non_basic) {
//TRACER_CALL(t,"LinearProgram::switchColumns");
#ifdef DEBUG
  cout << "switching basic column " << (*column)[j_basic] 
       << " with non-basic column " << (*column)[j_non_basic] << endl;
#endif
//update inverse of active constraint matrix
  Matrix<double> *temp_rhs=OPERATOR_NEW Matrix<double>(basic_number,1); 
  Matrix<double> *temp_soln=OPERATOR_NEW Matrix<double>(basic_number,1); 
  for (int i=0;i<basic_number;i++) {
    (*temp_rhs)(i,0)=A((*row)[i],(*column)[j_non_basic]);
  }
  char trans='n';
  double alpha=double_one_;
  double beta=double_zero_;
  int inc_rhs=1;
  int inc_soln=1;
  F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
    basic_inverse->addr(),&basic_number,temp_rhs->addr(),&inc_rhs,
    &beta,temp_soln->addr(),&inc_soln);
  double tj=double_one_/(*temp_soln)(j_basic,0);
  { for (int i=0;i<basic_number;i++) {
    if (i!=j_basic) {
      alpha=-(*temp_soln)(i,0)*tj;
      F77NAME(daxpy)(&basic_number,&alpha,
        basic_inverse->addr(j_basic,0),&basic_number,
        basic_inverse->addr(i,0),&basic_number);
    }
  } }
  F77NAME(dscal)(&basic_number,&tj,basic_inverse->addr(j_basic,0),
    &basic_number);
  delete temp_rhs;
  delete temp_soln;
#ifdef DEBUG
//cout << "\tupdated basic_inverse = " << endl;
//cout << *basic_inverse << endl;
#endif

//switch columns
  int colj=(*column)[j_non_basic];
  (*column)[j_non_basic]=(*column)[j_basic];
  (*column)[j_basic]=colj;
#ifdef DEBUG
//cout << "\tupdated column = " << endl;
//cout << *column << endl;
#endif
}
template<> int LinearProgram<double>::basicDualPivot(Matrix<double> &h,
double &ratio) {
//TRACER_CALL(t,"LinearProgram::basicDualPivot");
  int j_basic=basic_number;
  ratio=huge_;
  if (basic_number>0) {
    Matrix<double> rhs;
    rhs.copy(h);
    char trans='T';
    double alpha=double_one_;
    double beta=double_zero_;
    int inc_rhs=1;
    int inc_h=1;
    F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
      basic_inverse->addr(),&basic_number,rhs.addr(),&inc_rhs,
      &beta,h.addr(),&inc_h);
#ifdef DEBUG
    cout << "after solve, h = \n" << h << endl;
#endif
    for (int j=0;j<basic_number;j++) {
      double hj=h(j,0);
      if (hj<=double_zero_) continue;
      double yj=(*y_trans)(0,(*row)[j]);
      if (j_basic>=basic_number || hj*ratio<yj) {
	j_basic=j;
	ratio=yj/hj;
#ifdef DEBUG
	cout << "j,ratio = " << j << " " << ratio << endl;
#endif
      }
    }
  }
  return j_basic;
}

void testLinearProgram(double scalar) {
  scalar=double_one_;
//
//problem 6, Winston p. 286
//solves dual
//adds one variable
//y_trans=[1,0]
//Matrix<double> A(2,2); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,2); 
//c(0,0)=1.; c(0,1)= 2.;
//A(0,0)=1.; A(0,1)=-1.; b(0,0)= 3.;
//A(1,0)=1.; A(1,1)= 1.; b(1,0)= 1.;

//problem 6, Winston p. 286
//solves primal
//adds one variable
//x=[1,0]
//Matrix<double> A(2,2); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,2); 
//c(0,0)=-3.; c(0,1)=-1.;
//A(0,0)=-1.; A(0,1)=-1.; b(0,0)=-1.;
//A(1,0)= 1.; A(1,1)=-1.; b(1,0)=-2.;

//Leather Limited problem, Winston p. 121
//solves dual
//adds row 0 and col 1; adds row 1 and col 0
//y_trans = [20,20]
//Matrix<double> A(2,2); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,2); 
//c(0,0)=40.; c(0,1)=60.;
//A(0,0)= 1.; A(0,1)= 2.; b(0,0)= 4.;
//A(1,0)= 1.; A(1,1)= 1.; b(1,0)= 3.;

//Leather Limited problem, Winston p. 121
//solves primal
//adds row 1 and col 0; adds row 0 and col 1
//x = [20,20]
//Matrix<double> A(2,2); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,2); 
//c(0,0)=-4.; c(0,1)=-3.;
//A(0,0)=-1.; A(0,1)=-1.; b(0,0)=-40.;
//A(1,0)=-2.; A(1,1)=-1.; b(1,0)=-60.;

//Giapetto problem, Winston p. 45
//solves dual
//adds row 0 & col 2; adds row 1 and col 0; switch cols 1 and 2
//y_trans = [20,60]
//Matrix<double> A(2,3); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,3); 
//c(0,0)=100.; c(0,1)=80.; c(0,2)=40.;
//A(0,0)=  2.; A(0,1)= 1.; A(0,2)= 1.; b(0,0)= 3.;
//A(1,0)=  1.; A(1,1)= 1.; A(1,2)= 0.; b(1,0)= 2.;

//Giapetto problem, Winston p. 45
//solves primal
//adds row 2 & col 0; adds row 0 and col 1; switch rows 1 and 2
//x = [20,60]
//Matrix<double> A(3,2); 
//Matrix<double> b(3,1); 
//Matrix<double> c(1,2); 
//c(0,0)=-3.; c(0,1)=-2.;
//A(0,0)=-2.; A(0,1)=-1.; b(0,0)=-100.;
//A(1,0)=-1.; A(1,1)=-1.; b(1,0)= -80.;
//A(2,0)=-1.; A(2,1)= 0.; b(2,0)= -40.;

//Dorian problem, Winston p. 58
//solves dual
//adds row 0 & col 0; adds row 1 and col 1
//y_trans = [5,7.5]
//Matrix<double> A(2,2); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,2); 
//c(0,0)=50.; c(0,1)=100.;
//A(0,0)= 7.; A(0,1)=  2.; b(0,0)= 28.;
//A(1,0)= 2.; A(1,1)= 12.; b(1,0)= 24.;

//Auto company, Winston p. 61
//multiple solutions
//solves dual
//adds row 0 and col 0
//y_trans = [40,0]
//Matrix<double> A(2,2); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,2); 
//c(0,0)=120.; c(0,1)=50.;
//A(0,0)=  3.; A(0,1)= 1.; b(0,0)= 3.;
//A(1,0)=  2.; A(1,1)= 1.; b(1,0)= 2.;
//

//
//Problem 2, Winston p. 198
//solves primal
//adds row 0 and col 0
//x = [40,0]
//Matrix<double> A(2,2); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,2); 
//c(0,0)=-4.; c(0,1)= 1.;
//A(0,0)=-3.; A(0,1)=-1.; b(0,0)=-6.;
//A(1,0)= 1.; A(1,1)=-2.; b(1,0)= 0.;

//Auto company, Winston p. 63
//infeasible
//solves dual
//adds row 0 and col 0
//Matrix<double> A(2,4); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,4); 
//c(0,0)=120.; c(0,1)=50.; c(0,2)=-30.; c(0,3)=-20.;
//A(0,0)=  3.; A(0,1)= 1.; A(0,2)=  1.; A(0,3)=  0.; b(0,0)= 3.;
//A(1,0)=  2.; A(1,1)= 1.; A(1,2)=  0.; A(1,3)=  1.; b(1,0)= 2.;

//Winston p. 64
//unbounded
//Matrix<double> A(2,2); 
//Matrix<double> b(2,1); 
//Matrix<double> c(1,2); 
//c(0,0)= 1.; c(0,1)=-6.;
//A(0,0)= 1.; A(0,1)= 2.; b(0,0)= 2.;
//A(1,0)=-1.; A(1,1)= 1.; b(1,0)=-1.;

//Dakotta problem, Winston p. 280
//solves dual
//adds row 0 and col 2; adds row2 and col 1
//y_trans = [2,0,8]
//Matrix<double> A(3,3); 
//Matrix<double> b(3,1); 
//Matrix<double> c(1,3); 
//c(0,0)= 48.; c(0,1)=20. ; c(0,2)= 8. ;
//A(0,0)=  8.; A(0,1)= 4. ; A(0,2)= 2. ; b(0,0)= 60.;
//A(1,0)=  6.; A(1,1)= 2. ; A(1,2)= 1.5; b(1,0)= 30.;
//A(2,0)=  1.; A(2,1)= 1.5; A(2,2)= 0.5; b(2,0)= 20.;

//Diet problem, Winston p. 68
//solves dual
//adds row 0 and col 1; adds row 2 and col 2; switches row 0 with row 1
//y_trans=[0,2.5,7.5,0]
//Matrix<double> A(4,4); 
//Matrix<double> b(4,1); 
//Matrix<double> c(1,4); 
//c(0,0)= 50.; c(0,1)= 20.; c(0,2)= 30.; c(0,3)= 80.;
//A(0,0)=400.; A(0,1)=200.; A(0,2)=150.; A(0,3)=500.; b(0,0)= 500.;
//A(1,0)=  3.; A(1,1)=  2.; A(1,2)=  0.; A(1,3)=  0.; b(1,0)=   6.;
//A(2,0)=  2.; A(2,1)=  2.; A(2,2)=  4.; A(2,3)=  4.; b(2,0)=  10.;
//A(3,0)=  2.; A(3,1)=  4.; A(3,2)=  1.; A(3,3)=  5.; b(3,0)=   8.;
//

//Post Office, Winston p. 71
//solves dual
//adds row 3 and col 0
//adds row 5 and col 1
//adds row 0 and col 3
//adds row 2 and col 5
//adds row 6 and col 2
//y_trans = [1/3,0,1/3,1/3,0,1/3,0]
  Matrix<double> A(7,7); 
  Matrix<double> b(7,1); 
  Matrix<double> c(1,7); 
  for (int j=0;j<7;j++) {
    c(0,j)=1.;
    for (int k=0;k<5;k++) A((j+k)%7,j)=1.;
    A((j+5)%7,j)=A((j+6)%7,j)=0.;
  }
  b(0,0)=17.; b(1,0)=13.; b(2,0)=15.; b(3,0)=19.; b(4,0)=14.;
  b(5,0)=16.; b(6,0)=11.;

//int m=A.size(0),n=A.size(1);
  LinearProgram<double> LP(A,b,c);
  LP.findBasicFeasibleGuess();

  while (LP.currentStatus()==primal_feasible || 
  LP.currentStatus()==dual_feasible) {
    LP.simplexStep();
//  cout << LP << endl;
    cout << "\nstatus = " << LP.currentStatus()
         << endl;
    cout << "current primal value = " << LP.currentValue() <<endl;
    cout << "current primal solution = " << LP.currentPrimalSolution() 
         << endl;
    cout << "current dual solution = " << LP.currentDualSolution() 
         << endl;
    char dummy[80];
    cout << "hit RETURN to continue" << endl;
    cin.getline(dummy,80);
  }

  cout << "\n\nafter iteration, status = " << LP.currentStatus()
       << endl;
  cout << "final primal value = " << LP.currentValue() << endl;
  cout << "final primal solution = " << LP.currentPrimalSolution() 
       << endl;
  cout << "final dual solution = " << LP.currentDualSolution() << endl;
}
*/
