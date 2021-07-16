#include <complex>
#include <float.h>
#include <limits>
#include <math.h>
#include "Vector.H"
float float_mone_=-1.;
float float_one_=1.;
float float_zero_=0.;
float float_undefined_=numeric_limits<float>::infinity();
complex<float> complex_float_mone_=complex<float>(-1.,0.);
complex<float> complex_float_one_=complex<float>(1.,0.);
complex<float> complex_float_zero_=complex<float>(0.,0.);
complex<float> complex_float_undefined_=
  complex<float>(float_undefined_,float_undefined_);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  void F77NAME(saxpy)(const int &n,const float &a,const float *x,
    const int &incx,float *y,const int &incy);
  void F77NAME(slabad)(float &small,float &large);
  float F77NAME(slamch)(const char&);
  float F77NAME(dzasum)(const int&,const complex<float>*,const int&);
  float F77NAME(scnrm2)(const int &n,complex<float> *x,
    const int &incx);
  int F77NAME(ilaenv)(const int &ispec,const char *name,const char *opts,
    const int &n1,const int &n2,const int &n3,const int &n4);
  int F77NAME(ilatrans)(const char &trans);
  int F77NAME(ilauplo)(const char &trans);
  void F77NAME(caxpy)(const int &n,const complex<float> &a,
    const complex<float> *x,const int &incx,complex<float> *y,
    const int &incy);
  void F77NAME(ccopy)(const int &n,const complex<float> *sx,
    const int &incx,complex<float> *sy,const int &incy);
  complex<float> F77NAME(cdotc)(const int &n,const complex<float> *sx,
    const int &incx,const complex<float> *sy,const int &incy);
  complex<float> F77NAME(cdotu)(const int &n,const complex<float> *sx,
    const int &incx,const complex<float> *sy,const int &incy);
  void F77NAME(csrscl)(const int &n,const float &sa,complex<float> *sx,
    const int &incx);
  void F77NAME(cgbcon)(const char &norm,const int &n,const int &kl,
    const int &ku,const complex<float> *AB,const int &ldab,int *ipiv,
    const float &anorm,float &rcond,complex<float> *work,
    float *rwork,int &info);
  void F77NAME(cgbconnp)(const char &norm,const int &n,const int &kl,
    const int &ku,const complex<float> *AB,const int &ldab,
    const float &anorm,float &rcond,complex<float> *work,
    float *rwork,int &info);
  void F77NAME(cgbequ)(const int &m,const int &n,const int &kl,
    const int &ku,const complex<float> *AB,const int &ldab,float *r,
    float *c,float &rowcnd,float &colcnd,float &amax,int &info);
  void F77NAME(cgbamv)(const char &trans,const int &m,const int &n,
    const int &kl,const int &ku,const float &alpha,
    const complex<float> *A,const int &lda,const complex<float> *x,
    const int &incx,const float &beta,float *y,
    const int &incy);
  void F77NAME(cgbmv)(const char &trans,const int &m,const int &n,
    const int &kl,const int &ku,const complex<float> &alpha,
    const complex<float> *A,const int &lda,const complex<float> *x,
    const int &incx,const complex<float> &beta,complex<float> *y,
    const int &incy);
  void F77NAME(cgbtf2np)(const int &m,const int &n,const int &kl,
    const int &ku,complex<float> *AB,const int &ldab,int &info);
  void F77NAME(cgbtrf)(const int &m,const int &n,const int &kl,
    const int &ku,complex<float> *AB,const int &ldab,int *ipiv,
    int &info);
  void F77NAME(cgbtrs)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const complex<float> *AB,
    const int &ldab,const int *ipiv,complex<float> *b,const int &ldb,
    int &info);
  void F77NAME(cgbtrsnp)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const complex<float> *AB,
    const int &ldab,complex<float> *b,const int &ldb,int &info);
  void F77NAME(csscal)(const int &n,const float &a,
    complex<float> *x,const int &incx);
  void F77NAME(cgeamv)(const int &trans,const int &m,const int &n,
    float &alpha,const complex<float> *A,const int &lda,
    const complex<float> *x,const int &incx,const float &beta,
    float *y,const int &incy);
  void F77NAME(cgecon)(const char &norm,const int &n,
    const complex<float> *A,const int &lda,const float &anorm,
    float &rcond,complex<float> *work,float *rwork,int &info);
  void F77NAME(cgeequ)(const int &m,const int &n,
    const complex<float> *A,const int &lda,float *r,float *c,
    float &rowcnd,float &rolcnd,float &amax,int &info);
  void F77NAME(cgeev)(const char &jobvl,const char &jobvr,const int &n,
    complex<float> *A,const int &lda,complex<float> *w,
    complex<float> *vl,const int &ldvl,complex<float> *vr,
    const int &ldvr,complex<float> *work,const int &lwork,
    float *rwork,int &info);
  void F77NAME(cgels)(const char &trans,const int &m,const int &n,
    const int &nrhs,complex<float> *a,const int &lda,complex<float> *b,
    const int &ldb,complex<float> *work,const int &lwork,
    const int &info);
  void F77NAME(cgemm)(const char &transa,const char &transb,
    const int &m,const int &n,const int &k, 
    const complex<float> &alpha,const complex<float> *A,
    const int &lda,const complex<float> *B,const int &ldb,
    const complex<float> &beta,complex<float> *C,const int &ldc);
  void F77NAME(cgemv)(const char &trans,const int &m,const int &n,
    const complex<float> &alpha,const complex<float> *a,
    const int &lda,const complex<float> *x,const int &incx,
    const complex<float> &beta,complex<float> *y,const int &incy);
  void F77NAME(cgelqf)(const int &m,const int &n,complex<float> *A,
    const int &lda,complex<float> *tau,complex<float> *work,
    int &lwork,int &info);
  void F77NAME(cgeqp3)(const int &m,const int &n,complex<float> *A,
    const int &lda,int *jpvt,complex<float> *tau,complex<float> *work,
    const int &lwork,float *rwork,int &info);
  void F77NAME(cgeqrf)(const int &m,const int &n,complex<float> *A,
    const int &lda,complex<float> *tau,complex<float> *work,
    int &lwork,int &info);
  void F77NAME(cgerc)(const int&,const int&,const complex<float>&,
    const complex<float>*,const int&,const complex<float>*,const int&,
    complex<float>*,const int&);
  void F77NAME(cgeru)(const int&,const int&,const complex<float>&,
    const complex<float>*,const int&,const complex<float>*,const int&,
    complex<float>*,const int&);
  void F77NAME(cgesv)(int &n,int &nrhs,complex<float> *a,int &lda,
    int *ipiv,complex<float> *b, int &ldb, int &info);
  void F77NAME(cgesvd)(const char &jobu,const char &jobvt,const int &m, 
    const int &n,const complex<float> *A,const int &lda,
    float *S,complex<float> *U,const int &ldu,complex<float> *Vt,
    const int &ldvt,complex<float> *work,const int &lwork,
    float *rwork,int &info);
  void F77NAME(cgesc2)(const int &n,const complex<float> *A,
    const int &lda,complex<float> *rhs,const int *ipiv,const int *jpiv,
    complex<float> &scale);
  void F77NAME(cgetc2)(const int &n,complex<float> *A,const int &lda,
    int *ipiv,int *jpiv,int &info);
  void F77NAME(cgetrf)(int &m,int &n,complex<float> *a,int &lda,
    int *ipiv,int &info);
  void F77NAME(cgetrs)(const char &trans,const int &n,const int &nrhs,
    complex<float> *a,const int &lda,int *ipiv,complex<float> *b,
    const int &ldb,int &info);
  void F77NAME(cgetri)(int &n,complex<float> *a,int &lda,int *ipiv,
    complex<float> *work,int &lwork,int &info);
  void F77NAME(cgtcon)(const char &norm,const int &n,
    const complex<float> *L,const complex<float> *D,
    const complex<float> *U,const complex<float> *U2,
    const int *ipiv,const float &anorm,float &rcond,
    complex<float> *work,int &info);
  void F77NAME(cgtconnp)(const char &norm,const int &n,
    const complex<float> *L,const complex<float> *D,
    const complex<float> *U,const float &anorm,float &rcond,
    complex<float> *work,int &info);
  void F77NAME(cgtmv)(const int &n,const complex<float> &alpha,
    const complex<float> *L,const complex<float> *D,
    const complex<float> *U,const complex<float> *x,const int &incx,
    const complex<float> &beta,complex<float> *y,const int &incy);
  void F77NAME(cgtrfs)(const char &trans,const int &n,const int &nrhs,
    const complex<float> *dl,const complex<float> *d,
    const complex<float> *du,const complex<float> *dlf,
    const complex<float> *df,const complex<float> *duf,
    const complex<float> *du2,const int *ipiv,const complex<float> *B,
    const int &ldb,complex<float> *X,const int &ldx,float *ferr,
    float *berr,complex<float> *work,float *rwork,int &info);
  void F77NAME(cgtrfsnp)(const char &trans,const int &n,const int &nrhs,
    const complex<float> *dl,const complex<float> *d,
    const complex<float> *du,const complex<float> *dlf,
    const complex<float> *df,const complex<float> *duf,
    const complex<float> *B,const int &ldb,complex<float> *X,
    const int &ldx,float *ferr,float *berr,complex<float> *work,
    float *rwork,int &info);
  void F77NAME(cgtrfsr)(const char &trans,const int &n,const int &nrhs,
    const complex<float> *dl,const complex<float> *d,
    const complex<float> *du,const complex<float> *dlf,
    const complex<float> *df,const complex<float> *duf,
    const complex<float> *du2,const int *ipiv,const complex<float> *B,
    const int &ldb,complex<float> *X,const int &ldx,float *ferr,
    float *berr,complex<float> *work,float *rwork,int &info);
  void F77NAME(cgtrfsrnp)(const char &trans,const int &n,const int &nrhs,
    const complex<float> *dl,const complex<float> *d,
    const complex<float> *du,const complex<float> *dlf,
    const complex<float> *df,const complex<float> *duf,
    const complex<float> *B,const int &ldb,complex<float> *X,
    const int &ldx,float *ferr,float *berr,complex<float> *work,
    float *rwork,int &info);
  void F77NAME(chtmv)(const int &n,const complex<float> &alpha,
    const complex<float> *L,const float *D,const complex<float> *x,
    const int &incx,const complex<float> &beta,complex<float> *y,
    const int &incy);
  void F77NAME(zsteqr)(const char &compz,const int &n,float *D,
    complex<float> *E,complex<float> *Z,const int &ldz,float *work,
    int &info);
  void F77NAME(cgtsv)(const int &n,const int &nrhs,complex<float> *dl,
    complex<float> *d,complex<float> *du,complex<float> *b,
    const int &ldb,int &info);
  void F77NAME(cgtsvnp)(const int &n,const int &nrhs,complex<float> *dl,
    complex<float> *d,complex<float> *du,complex<float> *b,
    const int &ldb,int &info);
  void F77NAME(cgtsvr)(const int &n,const int &nrhs,complex<float> *dl,
    complex<float> *d,complex<float> *du,complex<float> *b,
    const int &ldb,int &info);
  void F77NAME(cgtsvrnp)(const int &n,const int &nrhs,complex<float> *dl,
    complex<float> *d,complex<float> *du,complex<float> *b,
    const int &ldb,int &info);
  void F77NAME(cgttrf)(const int &n,complex<float> *L,
    complex<float> *D,complex<float> *U,complex<float> *U2,int *ipiv,
    int &info);
  void F77NAME(cgttrfnp)(const int &n,complex<float> *L,
    complex<float> *D,complex<float> *U,int &info);
  void F77NAME(cgttrs)(const char &trans,const int &n,const int &nrhs,
    const complex<float> *L,const complex<float> *D,
    const complex<float> *U,const complex<float> *U2,
    const int *ipiv,complex<float> *B,const int &ldb,int &info);
  void F77NAME(chbev)(const char &jobz,const char &uplo,const int &n,
    const int &kd,complex<float> *AB,const int &ldab,float *w,
    complex<float> *Z,const int &ldz,complex<float> *work,float *rwork,
    int &info);
  void F77NAME(chbamv)(const char &uplo,const int &n,const int &k,
    const float &alpha,const complex<float> *A,const int &lda,
    const complex<float> *x,const int &incx,const float &beta,
    float *y,const int &incy);
  void F77NAME(chbmv)(const char &uplo,const int &n,const int &k,
    const complex<float> &alpha,const complex<float> *A,const int &lda,
    const complex<float> *x,const int &incx,const complex<float> &beta,
    complex<float> *y,const int &incy);
  void F77NAME(checon)(const char &uplo,const int &n,
    const complex<float> *A,const int &lda,const int *ipiv,
    const float &anorm,float &rcond,complex<float> *work,int &info);
  void F77NAME(cheequb)(const char &uplo,const int &n,
    const complex<float> *A,const int &lda,float *s,float &scond,
    float &amax,complex<float> *work,int &info);
  void F77NAME(cheev)(const char &jobz,const char &uplo,const int &n,
    complex<float> *A,const int &lda,float *w,complex<float> *work,
    const int &lwork,float *rwork,int &info);
  void F77NAME(chemm)(const char &side,const char &uplo,const int &m,
    const int &n,const complex<float> &alpha,const complex<float> *A,
    const int &lda,const complex<float> *B,const int &ldb,
    const complex<float> &beta,complex<float> *C,const int &ldc);
  void F77NAME(chemv)(const char &uplo,const int &n,
    const complex<float> &alpha,const complex<float> *A,
    const int &lda,const complex<float> *x,const int &incx,
    const complex<float> &beta,complex<float> *y,const int &incy);
  void F77NAME(cher)(const char &uplo,const int &n,
    const complex<float> &alpha,const complex<float> *X,
    const int &incx,complex<float> *A,const int &lda);
  void F77NAME(cher2)(const char &uplo,const int &n,
    const complex<float> &alpha,const complex<float> *X,
    const int &incx,const complex<float> *y,const int &ldy,
    complex<float> *A,const int &lda);
  void F77NAME(cher2k)(const char &uplo,const char &trans,const int &n,
    const int &k,const complex<float> &alpha,const complex<float> *A,
    const int &lda,const complex<float> *B,const int &ldb,
    const complex<float> &beta,complex<float> *C,const int &ldc);
  void F77NAME(cherk)(const char &uplo,const char &trans,const int &n,
    const int &k,const complex<float> &alpha,const complex<float> *A,
    const int &lda,const complex<float> &beta,complex<float> *C,
    const int &ldc);
  void F77NAME(chetrf)(const char &uplo,const int &n,complex<float> *a,
    const int &lda,int *ipiv,complex<float> *work,const int &lwork,
    int &info);
  void F77NAME(chetri)(const char &uplo,const int &n,complex<float> *a,
    const int &lda,int *ipiv,complex<float> *work,int &info);
  void F77NAME(chetrs)(const char &uplo,const int &n,const int &nrhs,
    complex<float> *a,const int &lda,int *ipiv,complex<float> *b,
    const int &ldb,int &info);
  void F77NAME(chseqr)(const char &job,const char &compz,const int &n,
    const int &ilo,const int &ihi,complex<float> *H,const int &ldh,
    complex<float> *w,complex<float> *Z,const int &ldz,
    complex<float> *work,const int &lwork,int &info);
  void F77NAME(clacn2)(const int &n,complex<float> *v,
    complex<float> *x,float &est,int &kase,int *isave);
  void F77NAME(clacpy)(const char &uplo,const int &m,const int &n,
    const complex<float> *A,const int &lda,complex<float> *B,
    const int &ldb);
  void F77NAME(cladiv)(const complex<float>&,const complex<float>&,
    const complex<float>&,const complex<float>&,complex<float>&,
    complex<float>&);
  void F77NAME(claic1)(const int &job,const int &j,
    const complex<float> *x,const float &sest,const complex<float> *w,
    const complex<float> &gamma,float &sestpr,complex<float> &s,
    complex<float> &c);
  float F77NAME(clangb)(const char &norm,const int &n,const int &kl,
    const int &ku,const complex<float> *AB,const int &ldab,
    float *work);
  float F77NAME(clange)(const char &norm,const int &m,const int &n,
    const complex<float> *A,const int &lda,float *work);
  float F77NAME(clanhb)(const char &norm,const char &uplo,const int &n,
    const int &k,const complex<float> *AB,const int &ldab,float *work);
  float F77NAME(clanhs)(const char &norm,const int &n,
    const complex<float> *A,const int &lda,float *work);
  float F77NAME(clangt)(const char &norm,const int &n,
    const complex<float> *L,const complex<float> *D,
    const complex<float> *U);
  float F77NAME(clanhe)(const char &norm,const char &uplo,const int &n,
    const complex<float> *A,const int &lda,float *work);
  float F77NAME(clanht)(const char &norm,const int &n,
    const float *D,const complex<float> *E);
  float F77NAME(clansb)(const char &norm,const char &uplo,
    const int &n,const int &k,const complex<float> *AB,const int &ldab,
    float *work);
  float F77NAME(clantb)(const char &norm,const char &uplo,
    const char &diag,const int &n,const int &k,const complex<float> *AB,
    const int &ldab,float *work);
  float F77NAME(clantr)(const char &norm,const char &uplo,
    const char &diag,const int &m,const int &n,const complex<float> *A,
    const int &lda,float *work);
  void F77NAME(claqgb)(const int &m,const int &n,const int &kl,
    const int &ku,complex<float> *AB,const int &ldab,float *r,
    float *c,const float &rowcnd,const float &colcnd,
    const float &amax,char &equed);
  void F77NAME(claqge)(const int &m,const int &n,complex<float> *A,
    const int &lda,const float *r,const float *c,const float &rowcnd,
    const float &colcnd,const float &amax,char &equed);
  void F77NAME(claqhe)(const char &uplo,const int &n,complex<float> *A,
    const int &lda,const float *s,const float &scond,const float &amax,
    char &equed);
  void F77NAME(claqsb)(const char &uplo,const int &n,const int &kd,
    complex<float> *AB,const int &ldab,const float *s,
    const float &scond,const float &amax,char &equed);
  void F77NAME(claqsy)(const char &uplo,const int &n,complex<float> *A,
    const int &lda,float *s,const float &scond,const float &amax,
    char &equed);
  void F77NAME(clargv)(const int &n,complex<float> *x,const int &incx,
    complex<float> *y,const int &incy,float *c,
    const int &incc);
  void F77NAME(clartv)(const int &n,complex<float> *x,const int &incx,
    complex<float> *y,const int &incy,const float *c,
    const complex<float> *s,const int &incc);
  void F77NAME(clascl)(const char &type,const int &kl,const int &ku,
    const float &cfrom,const float &cto,const int &m,const int &n,
    complex<float> *A,const int &lda,int &info);
  void F77NAME(claset)(const char &uplo,const int &m,const int &n,
    const complex<float> &alpha,const complex<float> &beta,
    complex<float> *A,const int &lda);
  void F77NAME(claswp)(const int &n,complex<float> *A,const int &lda,
    const int &k1,const int &k2,const int *ipiv,const int &incx);
  void F77_NAME(cla_heamv)(const int &uplo,const int &n,
    const float &alpha,const complex<float> *A,const int &lda,
    const complex<float> *x,const int &incx,const float &beta,
    float *y,const int &incy);
  float F77_NAME(cla_porpvgrw)(const char &uplo,const int &ncols,
    const complex<float> *A,const int &lda,const complex<float> *AF,
    const int &ldaf,complex<float> *work);
  float F77_NAME(cla_syrpvgrw)(const char &uplo,const int &n,
    const int &info,const complex<float> *A,const int &lda,
    const complex<float> *AF,const int &ldaf,const int *ipiv,
    complex<float> *work);
  void F77NAME(corgqr)(int &m,int &n,int &k,complex<float> *a,int &lda,
    const complex<float> *tau,complex<float> *work,int &lwork,
    int &info);
  void F77NAME(cpbcon)(const char &uplo,const int &n,const int &kd,
    const complex<float> *AB,const int &ldab,const float &anorm,
    float &rcond,complex<float> *work,float *rwork,int &info);
  void F77NAME(cpbequ)(const char &uplo,const int &n,const int &kd,
    const complex<float> *AB,const int &ldab,float *s,float &scond,
    float &amax,int &info);
  void F77NAME(cpbsv)(const char &uplo,const int &n,const int &kd,
    const int &nrhs,complex<float> *AB,const int &ldab,
    complex<float> *B,const int &ldb,int &info);
  void F77NAME(cpbtrf)(const char &uplo,const int &n,const int &kd,
    complex<float> *AB,const int &ldab,int &info);
  void F77NAME(cpbtrs)(const char &uplo,int &n,const int &kd,
    const int &nrhs,const complex<float> *AB,const int &ldab,
    complex<float> *b,const int &ldb,int &info);
  void F77NAME(cpocon)(const char &uplo,const int &n,
    const complex<float> *a,const int &lda,const complex<float> &anorm,
    float &rcond,complex<float> *work,float *rwork,int &info);
  void F77NAME(cpoequb)(const int &n,const complex<float> *A,
    const int &lda,float *s,float &scond,float &amax,int &info);
  void F77NAME(cpotrf)(const char &uplo,const int &n,complex<float> *A,
    const int &lda,int &info);
  void F77NAME(cpotri)(const char &uplo,const int &n,complex<float> *A,
    const int &lda,int &info);
  void F77NAME(cpotrs)(const char &uplo,const int &n,const int &nrhs,
    complex<float> *A,const int &lda,complex<float> *B,const int &ldb,
    int &info);
  void F77NAME(cptcon)(const int &n,const float *d,
    const complex<float> *e,const float &anorm,float &rcond,
    float *rwork,int &info);
  void F77NAME(cptrfs)(const char &uplo,const int &n,const int &nrhs,
    const float *d,const complex<float> *E,const float *df,
    const complex<float> *ef,const complex<float> *B,const int &ldb,
    complex<float> *X,const int &ldx,float *ferr,float *berr,
    complex<float> *work,float *rwork,int &info);
  void F77NAME(cptrfsr)(const char &uplo,const int &n,const int &nrhs,
    const float *d,const complex<float> *E,const float *df,
    const complex<float> *ef,const complex<float> *B,const int &ldb,
    complex<float> *X,const int &ldx,float *ferr,float *berr,
    complex<float> *work,float *rwork,int &info);
  void F77NAME(cptsv)(const int &n,const int &nrhs,float *d,
    complex<float> *e,complex<float> *b,const int &ldb,int &info);
  void F77NAME(cpttrf)(const int &n,float *d,complex<float> *e,
    int &info);
  void F77NAME(cpttrs)(const char &uplo,const int &n,const int &nrhs,
    float *d,complex<float> *e,complex<float> *B,const int &ldb,
    int &info);
  void F77NAME(crot)(const int &n,complex<float> *sx,const int &incx,
    complex<float> *sy,const int &incy,const float &c,
    const complex<float> &s);
  void F77NAME(crotg)(complex<float> &sa,complex<float> &sb,
    float &c,complex<float> &s);
  void F77NAME(cscal)(const int &n,const complex<float> &a,
    complex<float> *x,const int &incx);
  void F77NAME(cswap)(const int &n,complex<float> *sx,const int &incx,
    complex<float> *sy,const int &incy);
  void F77NAME(ctbsv)(const char &uplo,const char &trans,const char &diag,
    const int &n,const int &k,const complex<float> *A,const int &lda,
    complex<float> *x,const int &incx);
  void F77NAME(ctrcon)(const char &norm,const char &uplo,
    const char &diag,const int &n,const complex<float> *A,const int &lda,
    float &rcond,complex<float> *work,float *rwork,int &info);
  void F77NAME(ctrevc)(const char &side,const char &howmny,
    const bool *select,const int &n,complex<float> *T,const int &ldt,
    complex<float> *VL,const int &ldvl,complex<float> *VR,
    const int &ldvr,const int &mm,int &m,complex<float> *work,
    float *rwork,int &info);
  void F77NAME(ctrmv)(const char &uplo,const char &trans,
    const char &diag,const int &n,const complex<float> *A,const int &lda,
    complex<float> *x,const int &incx);
  void F77NAME(ctrsv)(const char&,const char&,const char&,const int &n,
    const complex<float> *A,const int &lda,complex<float> *x,
    const int &incx);
  void F77NAME(ctrmm)(const char &side,const char &uplo,
    const char &transa,const char &diag,const int &m,const int &n,
    const complex<float> &alpha,const complex<float> *A,const int &lda,
    complex<float> *B,const int &ldb);
  void F77NAME(ctrsm)(const char &side,const char &uplo,
    const char &transa,const char &diag,const int &m,const int &n,
    const complex<float> &alpha,const complex<float> *a,const int &lda,
    complex<float> *b,const int &ldb);
  void F77NAME(ctrtri)(const char &uplo,const char &diag,const int &n,
    complex<float> *A,const int &lda,int &info);
  void F77NAME(ctzrzf)(const int &m,const int &n,complex<float> *A,
    const int &lda,complex<float> *tau,complex<float> *work,
    const int &lwork,int &info);
  void F77NAME(cunglq)(const int &m,const int &n,const int &k,
    complex<float> *A,const int &lda,const complex<float> *tau,
    complex<float> *work,const int &lwork,int &info);
  void F77NAME(cungqr)(const int &m,const int &n,const int &k,
    complex<float> *A,const int &lda,const complex<float> *tau,
    complex<float> *work,const int &lwork,int &info);
  void F77NAME(cunmlq)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,complex<float> *A,const int &lda,
    complex<float> *tau,complex<float> *C,const int &ldc,
    complex<float> *work,int &lwork,int &info);
  void F77NAME(cunmqr)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,complex<float> *A,const int &lda,
    complex<float> *tau,complex<float> *C,const int &ldc,
    complex<float> *work,int &lwork,int &info);
  void F77NAME(cunmrz)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,const int &l,const complex<float> *A,
    const int &lda,const complex<float> *tau,complex<float> *C,
    const int &ldc,complex<float> *work,const int &lwork,int &info);
  int F77NAME(ilaclc)(const int &m,const int &n,const complex<float> *A,
    const int &lda);
  int F77NAME(ilaclr)(const int &m,const int &n,const complex<float> *A,
    const int &lda);
  int F77NAME(icamax)(const int &n,const complex<float> *x,
    const int &incx);
  int F77NAME(icmin)(const int &n,const complex<float> *a,
    const int &inca);
  int F77NAME(icsumn)(const int &n,const complex<float> *a,
    const int &inca);
  int F77NAME(icsump)(const int &n,const complex<float> *a,
    const int &inca);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Vector.C"
#ifdef DEBUG
template<> const complex<float>
  Vector<float,complex<float> >::undefined_(
    numeric_limits<float>::infinity(),numeric_limits<float>::infinity());
#endif
template<> const complex<float>
  Vector<float,complex<float> >::zero_(float_zero_,float_zero_);
template<> const complex<float>
  Vector<float,complex<float> >::one_(float_one_,float_zero_);
template<> const complex<float>
  Vector<float,complex<float> >::mone_(float_mone_,float_zero_);

template<> int Vector<float,complex<float> >::amax() const { 
  return F77NAME(icamax)(sz,data,1)-1;
}
template<> float Vector<float,complex<float> >::asum() const {
  return F77NAME(dzasum)(sz,data,1);
}
template<> void Vector<float,complex<float> >::axpy(complex<float> a,
const Vector<float,complex<float> > &x) {
//TRACER_CALL(t,"Vector<float,complex<float> >::axpy");
#ifdef DEBUG
//cout << "\ta = " << a << endl;
//cout << "\tsz,x.sz = " << sz << " " << x.sz << endl;
//printOn(cout);
//x.printOn(cout);
#endif
  F77NAME(caxpy)(min(sz,x.sz),a,x.data,1,data,1);
#ifdef DEBUG
//printOn(cout);
#endif
}
template<> complex<float> Vector<float,complex<float> >::dot(
const Vector<float,complex<float> > &x) const {
  OBSOLETE(0);
}
template<> complex<float> Vector<float,complex<float> >::dotc(
const Vector<float,complex<float> > &x) const {
  return F77NAME(cdotc)(min(sz,x.sz),x.data,1,data,1);
}
template<> complex<float> Vector<float,complex<float> >::dotu(
const Vector<float,complex<float> > &x) const {
  return F77NAME(cdotu)(min(sz,x.sz),x.data,1,data,1);
}
template<> float Vector<float,complex<float> >::nrm2() const {
  return F77NAME(scnrm2)(sz,data,1);
}
template<> void Vector<float,complex<float> >::rot(
Vector<float,complex<float> > &x,float c,float s) {
  F77NAME(crot)(min(sz,x.sz),x.data,1,data,1,c,s);
}
template<> void Vector<float,complex<float> >::scal(complex<float> a)
{
  F77NAME(cscal)(sz,a,data,1);
}
template<> void Vector<float,complex<float> >::swap(
Vector<float,complex<float> > &x) {
  F77NAME(cswap)(min(sz,x.sz),x.data,1,data,1);
}
template<> void Vector<float,complex<float> >::largv(
Vector<float,complex<float> > &y,Vector<float,float> &c) {
  int n=size();
  CHECK_SAME(n,y.size());
  CHECK_SAME(n,c.size());
  F77NAME(clargv)(n,data,1,y.data,1,c.addr(),1);
}
template<> void Vector<float,complex<float> >::lartv(
Vector<float,complex<float> > &y,
const Vector<float,float> &c,
const Vector<float,complex<float> > &s) {
  int n=size();
  CHECK_SAME(n,y.size());
  CHECK_SAME(n,c.size());
  CHECK_SAME(n,s.size());
  F77NAME(clartv)(n,addr(),1,y.addr(),1,c.addr(),s.addr(),1);
}

template class Vector<float,complex<float> >; 
template void testVector(const float&,const complex<float>&);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Matrix.C"
#ifdef DEBUG
template<> const complex<float>
  Matrix<float,complex<float> >::undefined_(
    numeric_limits<float>::infinity(),numeric_limits<float>::infinity());
#endif
template<> const complex<float>
  Matrix<float,complex<float> >::zero_(float_zero_,float_zero_);
template<> const complex<float>
  Matrix<float,complex<float> >::one_(float_one_,float_zero_);
template<> const complex<float>
  Matrix<float,complex<float> >::mone_(float_mone_,float_zero_);

/*
template<> Matrix<float,complex<float> >*
Matrix<float,complex<float> >::transpose() const {
  int m=size(0),n=size(1);
  Matrix<float,complex<float> > *T=
    OPERATOR_NEW Matrix<float,complex<float> >(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(m,addr(0,j),1,T->addr(j,0),n);
  }
  return T;
}

template<> Matrix<float,complex<float> >*
Matrix<float,complex<float> >::conjugateTranspose() const {
  int m=size(0),n=size(1);
  Matrix<float,complex<float> > *T=
    OPERATOR_NEW Matrix<float,complex<float> >(n,m);
  for (int j=0;j<n;j++) {
    const complex<float> *col_j=addr(0,j);
    complex<float> *T_row_j=T->addr(j,0);
    for (int i=0;i<m;i++,col_j++,T_row_j+=n) {
      *T_row_j=conj(*col_j); 
    }
  }
  return T;
}
*/

template<> void Matrix<float,complex<float> >::interchangeColumns(
int i,int j) {
  int m=size(0),n=size(1);
  CHECK_BOUNDS(i,0,n)
  CHECK_BOUNDS(j,0,n)
  F77NAME(cswap)(m,addr(0,i),1,addr(0,j),1);
}

template<> void Matrix<float,complex<float> >::interchangeRows(
int i,int j) {
  int m=size(0),n=size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,m)
  F77NAME(cswap)(n,addr(i,0),m,addr(j,0),m);
}

template<> void Matrix<float,complex<float> >::gemv(
complex<float> alpha,const Vector<float,complex<float> > &x,
complex<float> beta,Vector<float,complex<float> > &y,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(x.size(),m);
    CHECK_SAME(y.size(),n);
    F77NAME(cgemv)(trans,m,n,alpha,addr(),m,x.addr(),1,beta,y.addr(),1);
  } else {
    CHECK_SAME(x.size(),n);
    CHECK_SAME(y.size(),m);
    F77NAME(cgemv)(trans,m,n,alpha,addr(),m,x.addr(),1,beta,y.addr(),1);
  }
}

template<> void Matrix<float,complex<float> >::ger(
complex<float> alpha,const Vector<float,complex<float> > &x,
const Vector<float,complex<float> > &y) {
  OBSOLETE("inappropriate for this class");
}

template<> void Matrix<float,complex<float> >::gerc(
complex<float> alpha,const Vector<float,complex<float> > &x,
const Vector<float,complex<float> > &y) {
  int m=size(0),n=size(1);
  CHECK_SAME(x.size(),m);
  CHECK_SAME(y.size(),n);
  F77NAME(cgerc)(m,n,alpha,x.addr(),1,y.addr(),1,addr(),m);
}

template<> void Matrix<float,complex<float> >::geru(
complex<float> alpha,const Vector<float,complex<float> > &x,
const Vector<float,complex<float> > &y) {
  int m=size(0),n=size(1);
  CHECK_SAME(x.size(),m);
  CHECK_SAME(y.size(),n);
  F77NAME(cgeru)(m,n,alpha,x.addr(),1,y.addr(),1,addr(),m);
}

template<> void Matrix<float,complex<float> >::gemm(
complex<float> alpha,const Matrix<float,complex<float> > &A,
const Matrix<float,complex<float> > &B,complex<float> beta,
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
  F77NAME(cgemm)(transa,transb,m,n,k,alpha,A.addr(),A.size(0),
    B.addr(),B.size(0),beta,addr(),m);
}

template<> float Matrix<float,complex<float> >::equilibrate(
Vector<float,float> &r,Vector<float,float> &c,
float &rowcnd,float &colcnd) const {
  int m=size(0),n=size(1);
  float amax=numeric_limits<float>::infinity();
  int info;
  F77NAME(cgeequ)(m,n,addr(),m,r.addr(),c.addr(),rowcnd,colcnd,
    amax,info);
  CHECK_SAME(info,0);
  return amax;
}

template<> void Matrix<float,complex<float> >::copyFrom(char uplo,
int m,int n,const Matrix<float,complex<float> > &A) {
  int s0=size(0),as0=A.size(0);
  m=min(m,min(s0,as0));
  n=min(n,min(size(1),A.size(1)));
  if (uplo=='A' || uplo=='a') {
    F77NAME(clacpy)(uplo,m,n,A.addr(),as0,addr(),s0);
//  otherwise, slacpy only copies a triangular part:
  } else if (uplo=='L' || uplo=='l') {
    for (int j=0;j<n;j++) {
      if (j<m) {
        F77NAME(ccopy)(m-j,A.addr(j,j),1,addr(j,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+1,m),A.addr(0,j),1,addr(0,j),1);
    }
  }
}

template<> void Matrix<float,complex<float> >::scale(char type,
int kl,int ku,float denominator,float numerator) {
  int m=size(0),info;
  F77NAME(clascl)(type,kl,ku,denominator,numerator,m,size(1),addr(),
    m,info);
  CHECK_SAME(info,0);
}

template<> void Matrix<float,complex<float> >::set(char uplo,
complex<float> offdiag,complex<float> diag) {
  int m=size(0);
  F77NAME(claset)(uplo,m,size(1),offdiag,diag,addr(),m);
}

template<> int Matrix<float,complex<float> >::lastNonzeroColumn()
const {
  int m=size(0);
  return F77NAME(ilaclc)(m,size(1),addr(),m);
}

template<> int Matrix<float,complex<float> >::lastNonzeroRow() const {
  int m=size(0);
  return F77NAME(ilaclr)(m,size(1),addr(),m);
}

template<> float Matrix<float,complex<float> >::normFrobenius()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(clange)('F',m,size(1),addr(),m,work);
}

template<> float Matrix<float,complex<float> >::normInfinity() const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float val=F77NAME(clange)('I',m,size(1),addr(),m,work);
  delete [] work;
  return val;
}

template<> float Matrix<float,complex<float> >::normMaxEntry() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clange)('M',m,size(1),addr(),m,work);
}

template<> float Matrix<float,complex<float> >::normOne() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clange)('1',m,size(1),addr(),m,work);
}

template<> Matrix<float,complex<float> >*
Matrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(M.size(0),k);
  Matrix<float,complex<float> > *P=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(cgemm)('N','N',m,n,k,complex_float_one_,addr(),m,M.addr(),k,
    complex_float_zero_,P->addr(),m);
  return P;
}

template<> Vector<float,complex<float> >*
Matrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(v.size(),n);
  Vector<float,complex<float> > *P=
    OPERATOR_NEW Vector<float,complex<float> >(m);
  F77NAME(cgemv)('N',m,n,complex_float_one_,addr(),m,v.addr(),1,
    complex_float_zero_,P->addr(),1);
  return P;
}

template<> void Matrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  int m=size(0), n=size(1);
  Matrix<float,complex<float> > AF(*this);
  int info;
  if (m==n) {
    CHECK_SAME(m,b.size())
    CHECK_SAME(n,x.size())
    x.copy(b);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(cgetrf)(m,m,AF.addr(),m,ipiv,info);
    if (info==0) {
      F77NAME(cgetrs)(trans,m,1,AF.addr(),m,ipiv,x.addr(),m,info);
    }
    CHECK_SAME(info,0)
    delete [] ipiv;
  } else {
    complex<float> w=complex_float_undefined_;
    int lwork=-1;
    bool transposed=(trans!='N' && trans!='n');
    if (m>n) {
      if (transposed) {
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())
        x.copyFrom(n,b);
        F77NAME(cgels)(trans,m,n,1,AF.addr(),m,x.addr(),m,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cgels)(trans,m,n,1,AF.addr(),m,x.addr(),m,work,
          lwork,info);
        CHECK_SAME(info,0)
        delete [] work;
      } else {
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        Vector<float,complex<float> > xtmp(m);
        xtmp.copy(b);
        F77NAME(cgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),m,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),m,work,
          lwork,info);
        CHECK_SAME(info,0)

        x.copyFrom(n,xtmp);
        delete [] work;
      }
    } else {
      if (transposed) {
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())
        Vector<float,complex<float> > xtmp(n);
        xtmp.copy(b);
        F77NAME(cgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),n,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),n,work,
          lwork,info);
        CHECK_SAME(info,0)

        x.copyFrom(m,xtmp);
        delete [] work;
      } else {
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        x.copyFrom(m,b);
        F77NAME(cgels)(trans,m,n,1,AF.addr(),m,x.addr(),n,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cgels)(trans,m,n,1,AF.addr(),m,x.addr(),n,work,
          lwork,info);
        CHECK_SAME(info,0)
        delete [] work;
      }
    }
  }
}

template<> void Matrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  int m=size(0), n=size(1);
  Matrix<float,complex<float> > AF(*this);
  int info;
  bool left_side=(side=='L' || side=='l');
  bool transposed=(trans!='N' && trans!='n');
  if (m==n) {
    X.copy(B);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(cgetrf)(m,m,AF.addr(),m,ipiv,info);
    if (info==0) {
      if (left_side) {
        int k=B.size(1);
        CHECK_SAME(k,X.size(1))
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        F77NAME(cgetrs)(trans,m,k,AF.addr(),m,ipiv,X.addr(),m,info);
      } else {
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        CHECK_SAME(m,B.size(1))
        CHECK_SAME(n,X.size(1))
        if (transposed) {
          for (int j=0;j<m-1;j++) {
            if (ipiv[j]-1!=j) {
              F77NAME(cswap)(k,X.addr(0,j),1,X.addr(0,ipiv[j]-1),1);
            }
          }
          F77NAME(ctrsm)('R','L','T','U',k,m,complex_float_one_,
            AF.addr(),m,X.addr(),k);
          F77NAME(ctrsm)('R','U','T','N',k,m,complex_float_one_,
            AF.addr(),m,X.addr(),k);
        } else {
          F77NAME(ctrsm)('R','U','N','N',k,m,complex_float_one_,
            AF.addr(),m,X.addr(),k);
          F77NAME(ctrsm)('R','L','N','U',k,m,complex_float_one_,
            AF.addr(),m,X.addr(),k);
          for (int j=m-2;j>=0;j--) {
            if (ipiv[j]-1!=j) {
              F77NAME(cswap)(k,X.addr(0,j),1,X.addr(0,ipiv[j]-1),1);
            }
          }
        }
      }
    }
    CHECK_SAME(info,0)
    delete [] ipiv; ipiv=0;
  } else {
    complex<float> w=complex_float_undefined_;
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
          F77NAME(cgels)(trans,m,n,k,AF.addr(),m,X.addr(),m,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w.real());
          complex<float> *work=
            OPERATOR_NEW_BRACKET(complex<float>,lwork);
          F77NAME(cgels)(trans,m,n,k,AF.addr(),m,X.addr(),m,work,
            lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        } else {
          CHECK_SAME(m,B.size(0))
          CHECK_SAME(n,X.size(0))
          Matrix<float,complex<float> > Xtmp(B);
          F77NAME(cgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),m,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w.real());
          complex<float> *work=
            OPERATOR_NEW_BRACKET(complex<float>,lwork);
          F77NAME(cgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),m,work,
            lwork,info);
          CHECK_SAME(info,0)

          X.copyFrom('A',n,k,Xtmp);
          delete [] work; work=0;
        }
      } else { // m>n, right side
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        complex<float> *tau=OPERATOR_NEW_BRACKET(complex<float>,n);
        F77NAME(cgeqrf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cgeqrf)(m,n,AF.addr(),m,tau,work,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        if (transposed) {
          CHECK_SAME(m,B.size(1))
          CHECK_SAME(n,X.size(1))
          Matrix<float,complex<float> > Xtmp(B);

          lwork=-1;
          F77NAME(cunmqr)('R','N',k,m,n,AF.addr(),m,tau,Xtmp.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w.real());
          work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
          F77NAME(cunmqr)('R','N',k,m,n,AF.addr(),m,tau,Xtmp.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
          X.copyFrom('A',k,n,Xtmp);

          F77NAME(ctrsm)('R','U','T','N',k,n,complex_float_one_,
            AF.addr(),m,X.addr(),k);
        } else {
          CHECK_SAME(n,B.size(1))
          CHECK_SAME(m,X.size(1))
          X=zero_;
          X.copyFrom('A',k,n,B);
          F77NAME(ctrsm)('R','U','N','N',k,n,complex_float_one_,
            AF.addr(),m,X.addr(),k);
          lwork=-1;
          F77NAME(cunmqr)('R','C',k,m,n,AF.addr(),m,tau,X.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w.real());
          work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
          F77NAME(cunmqr)('R','C',k,m,n,AF.addr(),m,tau,X.addr(),k,
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
          Matrix<float,complex<float> > Xtmp(B);
          F77NAME(cgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),n,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w.real());
          complex<float> *work=
            OPERATOR_NEW_BRACKET(complex<float>,lwork);
          F77NAME(cgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),n,work,
            lwork,info);
          CHECK_SAME(info,0)

          X.copyFrom('A',m,k,Xtmp);
          delete [] work; work=0;
        } else {
          CHECK_SAME(m,B.size(0))
          CHECK_SAME(n,X.size(0))
          X.copyFrom('A',m,k,B);
          F77NAME(cgels)(trans,m,n,k,AF.addr(),m,X.addr(),n,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w.real());
          complex<float> *work=
            OPERATOR_NEW_BRACKET(complex<float>,lwork);
          F77NAME(cgels)(trans,m,n,k,AF.addr(),m,X.addr(),n,work,
            lwork,info);
          CHECK_SAME(info,0)
          delete [] work;
        }
      } else { // m<n, right side
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        complex<float> *tau=OPERATOR_NEW_BRACKET(complex<float>,n);
        F77NAME(cgelqf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)
        lwork=static_cast<int>(w.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cgelqf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        if (transposed) {
          CHECK_SAME(m,B.size(1))
          CHECK_SAME(n,X.size(1))
          X=zero_;
          X.copyFrom('A',k,m,B);
          F77NAME(ctrsm)('R','L','T','N',k,m,complex_float_one_,
            AF.addr(),m,X.addr(),k);
          lwork=-1;
          F77NAME(cunmlq)('R','N',k,n,m,AF.addr(),m,tau,X.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w.real());
          work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
          F77NAME(cunmlq)('R','N',k,n,m,AF.addr(),m,tau,X.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        } else {
          CHECK_SAME(n,B.size(1))
          CHECK_SAME(m,X.size(1))
          Matrix<float,complex<float> > Xtmp(B);
          lwork=-1;
          F77NAME(cunmlq)('R','C',k,n,m,AF.addr(),m,tau,Xtmp.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w.real());
          work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
          F77NAME(cunmlq)('R','C',k,n,m,AF.addr(),m,tau,Xtmp.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
          X.copyFrom('A',k,m,Xtmp);

          F77NAME(ctrsm)('R','L','N','N',k,m,complex_float_one_,
            AF.addr(),m,X.addr(),k);
        }
        delete [] tau; tau=0;
      }
    }
  }
}

template class Matrix<float,complex<float> >; 
template void testMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SquareMatrix.C"

template<> SquareMatrix<float,complex<float> >*
SquareMatrix<float,complex<float> >::operator*(
const SquareMatrix<float,complex<float> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *P=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  F77NAME(cgemm)('N','N',n,n,n,complex_float_one_,addr(),n,S.addr(),n,
    complex_float_zero_,P->addr(),n);
  return P;
}

template<> Matrix<float,complex<float> >*
SquareMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *P=(m==n ?
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n) :
    OPERATOR_NEW Matrix<float,complex<float> >(m,n));
  F77NAME(cgemm)('N','N',m,n,m,complex_float_one_,addr(),m,M.addr(),m,
    complex_float_zero_,P->addr(),m);
  return P;
}

/*
template<> SquareMatrix<float,complex<float> >*
SquareMatrix<float,complex<float> >::transpose() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *T=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n,addr(0,j),1,T->addr(j,0),n);
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >*
SquareMatrix<float,complex<float> >::conjugateTranspose() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *T=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    const complex<float> *col_j=addr(0,j);
    complex<float> *T_row_j=T->addr(j,0);
    for (int i=0;i<n;i++,col_j++,T_row_j+=n) {
      *T_row_j=conj(*col_j);
    }
  }
  return T;
}
*/

template<> float 
SquareMatrix<float,complex<float> >::reciprocalConditionNumber(
char norm) const {
  int n=size(0),info;
  float rcond;
  float *rwork=OPERATOR_NEW_BRACKET(float,2*n);
  float anorm=F77NAME(clange)(norm,n,n,addr(),n,rwork);

  SquareMatrix<float,complex<float> > *AF=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(cgetrf)(n,n,AF->addr(0,0),n,ipiv,info);
  CHECK_SAME(info,0)
  delete [] ipiv; ipiv=0;

  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(cgecon)(norm,n,AF->addr(),n,anorm,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] iwork; iwork=0;
  delete [] rwork; rwork=0;
  delete[] work; work=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SquareMatrix<float,complex<float> >*
SquareMatrix<float,complex<float> >::inverse() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *Ainv=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  *Ainv = *this;
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(cgetrf)(n,n,Ainv->addr(0,0),n,ipiv,info);
  CHECK_SAME(info,0)

  complex<float> w(complex_float_undefined_);
  int lwork=-1;
  F77NAME(cgetri)(n,Ainv->addr(0,0),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(cgetri)(n,Ainv->addr(0,0),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)

  delete [] ipiv;
  delete [] work;
  return Ainv;
}
*/

template<> Vector<float,complex<float> >*
SquareMatrix<float,complex<float> >::eigenvalues(
SquareMatrix<float,complex<float> > *&V,
SquareMatrix<float,complex<float> > *&U) const {
  int n=size(0);
  if (V!=0) CHECK_SAME(n,V->size(0));
  if (U!=0) CHECK_SAME(n,U->size(0));
  char jobvl=(V==0 ? 'N' : 'V');
  char jobvr=(U==0 ? 'N' : 'V');
  complex<float> *vl=(V==0 ? 0 : V->addr());
  complex<float> *vr=(U==0 ? 0 : U->addr());
  SquareMatrix<float,complex<float> > AF(n);
  AF.copy(*this);
  Vector<float,complex<float> > *lambda =
    OPERATOR_NEW Vector<float,complex<float> >(n);
  complex<float> w;
  float *rwork=OPERATOR_NEW_BRACKET(float,2*n);
  int lwork=-1,info;
  F77NAME(cgeev)(jobvl,jobvr,n,AF.addr(),n,lambda->addr(),vl,n,vr,n,&w,
    lwork,rwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(cgeev)(jobvl,jobvr,n,AF.addr(),n,lambda->addr(),vl,n,vr,n,
    work,lwork,rwork,info);
  CHECK_TEST(info==0);
  delete [] work;
  delete [] rwork;
  return lambda;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const SquareMatrix<float,complex<float> > &S) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<float,complex<float> > *P=(m==n ?
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n) :
    OPERATOR_NEW Matrix<float,complex<float> >(m,n));
  F77NAME(cgemm)('N','N',m,n,n,complex_float_one_,M.addr(),m,S.addr(),n,
    complex_float_zero_,P->addr(),m);
  return P;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "TrapezoidalMatrix.C"
template<> complex<float>
  TrapezoidalMatrix<float,complex<float> >::outofbounds_=
  complex<float>(float_zero_,float_zero_);
template<> complex<float>
  TrapezoidalMatrix<float,complex<float> >::safety_=
  complex<float>(float_zero_,float_zero_);

template class TrapezoidalMatrix<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> Matrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::makeMatrix() const {
  int m=size(0),n=size(1);
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(this)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(m-j,addr(j,j),1,M->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*M)(j,j)=complex_float_one_;
      if (j+1<m) F77NAME(ccopy)(m-j-1,addr(j+1,j),1,M->addr(j+1,j),1);
    }
  }
  return M;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >& 
LowerTrapezoidalMatrix<float,complex<float> >::operator+=(
const LowerTrapezoidalMatrix<float,complex<float> > &L) {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0))
  CHECK_SAME(n,L.size(1))
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m-j,one_,L.addr(j,j),1,addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(caxpy)(m-j-1,one_,L.addr(j+1,j),1,addr(j+1,j),1);
      }
      (*this)(j,j)+=complex_float_one_;
    }
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >&
LowerTrapezoidalMatrix<float,complex<float> >::operator-=(
const LowerTrapezoidalMatrix<float,complex<float> > &L) {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0))
  CHECK_SAME(n,L.size(1))
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m-j,mone_,L.addr(j,j),1,addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(caxpy)(m-j-1,mone_,L.addr(j+1,j),1,addr(j+1,j),1);
      }
      (*this)(j,j)-=complex_float_one_;
    }
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >& 
LowerTrapezoidalMatrix<float,complex<float> >::operator*=(
complex<float> d) {
  int m=size(0),n=size(1);
  for (int j=0;j<n;j++) {
    F77NAME(cscal)(m-j,d,addr(j,j),1);
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >&
LowerTrapezoidalMatrix<float,complex<float> >::operator/=(
complex<float> d) {
  int m=size(0),n=size(1);
  complex<float> dinv=one_/d;
  for (int j=0;j<n;j++) {
    F77NAME(cscal)(m-j,dinv,addr(j,j),1);
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<float,complex<float> > *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(m,n);
  S->copy(*this);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m-j,complex_float_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m-j-1,complex_float_one_,L.addr(j+1,j),1,
        S->addr(j+1,j),1);
      (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *S=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(m-j,complex_float_one_,addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<float,complex<float> > *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(m,n);
  S->copy(*this);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m-j,complex_float_mone_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m-j-1,complex_float_mone_,L.addr(j+1,j),1,
        S->addr(j+1,j),1);
      (*S)(j,j)-=complex_float_one_;
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *D=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(claset)('U',m,n,complex_float_zero_,complex_float_zero_,
    D->addr(),m);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(m-j,addr(j,j),1,D->addr(j,j),1);
    F77NAME(caxpy)(m,complex_float_mone_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >* 
LowerTrapezoidalMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
// to compute the jth column of the product
//   [ L_11   0  ] [ 0 ] = [    0   ]
//   [ L_21 L_22 ] [ m ] = [ L_22 m ]
//   [ L_31 L_32 ]       = [ L_32 m ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  LowerTrapezoidalMatrix<float,complex<float> > *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(m,n);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(k-j,L.addr(j,j),1,P->addr(j,j),1);
      F77NAME(ctrmv)('L','N','N',k-j,addr(j,j),m,P->addr(j,j),1);
      if (m>k) {
        F77NAME(cgemv)('N',m-k,k-j,complex_float_one_,addr(k,j),m,
          L.addr(j,j),1,complex_float_zero_,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*P)(j,j)=(*this)(j,j);
      if (j<k-1) {
        F77NAME(ccopy)(k-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
        F77NAME(ctrmv)('L','N','N',k-j-1,addr(j+1,j+1),m,
          P->addr(j+1,j),1);
        F77NAME(caxpy)(k-j-1,complex_float_one_,addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
      if (m>k) {
        if (j<k-1) {
          F77NAME(cgemv)('N',m-k,k-j-1,complex_float_one_,addr(k,j+1),m,
            L.addr(j+1,j),1,complex_float_zero_,P->addr(k,j),1);
        }
        F77NAME(caxpy)(m-k,complex_float_one_,addr(k,j),1,
          P->addr(k,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<float,complex<float> > *P=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(clacpy)('A',k,n,M.addr(),k,P->addr(),m);
  F77NAME(ctrmm)('L','L','N','N',k,n,complex_float_one_,addr(),m,
    P->addr(),m);
  if (m>k) {
    F77NAME(cgemm)('N','N',m-k,n,k,complex_float_one_,addr(k,0),m,
      M.addr(),k,complex_float_zero_,P->addr(k,0),m);
  }
  return P;
}

template<> Vector<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<float,complex<float> > *p=
    OPERATOR_NEW Vector<float,complex<float> >(m);
  F77NAME(ccopy)(n,v.addr(),1,p->addr(),1);
  F77NAME(ctrmv)('L','N','N',n,addr(),m,p->addr(),1);
  if (m>n) {
    F77NAME(cgemv)('N',m-n,n,complex_float_one_,addr(n,0),m,v.addr(),1,
      complex_float_zero_,p->addr(n),1);
  }
  return p;
}

/*
template<> UpperTrapezoidalMatrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::transpose() const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<float,complex<float> > *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(m-j,addr(j,j),1,U->addr(j,j),n);
  }
  return U;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::conjugateTranspose()
const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<float,complex<float> > *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(n,m);
  for (int j=0;j<n;j++) {
    const complex<float> *col_j=addr(j,j);
    complex<float> *U_row_j=U->addr(j,j);
    for (int i=j;i<m;i++,col_j++,U_row_j+=n) {
      *U_row_j=conj(*col_j);
    }
  }
  return U;
}
*/

template<> Vector<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::trmv(
const Vector<float,complex<float> > &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<float,complex<float> > *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<float,complex<float> >(n);
    F77NAME(ccopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(ctrmv)('L',trans,'N',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(cgemv)(trans,m-n,n,one_,addr(n,0),m,x.addr(n),1,one_,
        p->addr(),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<float,complex<float> >(m);
    F77NAME(ccopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(ctrmv)('L','N','N',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(cgemv)('N',m-n,n,one_,addr(n,0),m,x.addr(),1,zero_,
        p->addr(n),1);
    }
  }
  return p;
}

template<> Matrix<float,complex<float> >*
LowerTrapezoidalMatrix<float,complex<float> >::trmm(
const Matrix<float,complex<float> > &X,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<float,complex<float> > *P=0;
  if (side=='L' || side=='l') {
    int k=X.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,X.size(0));
      P=OPERATOR_NEW Matrix<float,complex<float> >(n,k);
      F77NAME(clacpy)('A',n,k,X.addr(),m,P->addr(),n);
      F77NAME(ctrmm)('L','L',trans,'N',n,k,one_,addr(),m,P->addr(),n);
      if (m>n) {
        F77NAME(cgemm)(trans,'N',n,k,m-n,one_,addr(n,0),m,X.addr(n,0),m,
          one_,P->addr(),n);
      }
    } else {
      CHECK_SAME(n,X.size(0));
      P=OPERATOR_NEW Matrix<float,complex<float> >(m,k);
      F77NAME(clacpy)('A',n,k,X.addr(),n,P->addr(),m);
      F77NAME(ctrmm)('L','L','N','N',n,k,one_,addr(),m,P->addr(),m);
      if (m>n) {
        F77NAME(cgemm)('N','N',m-n,k,n,one_,addr(n,0),m,X.addr(),n,
          zero_,P->addr(n,0),m);
      }
    }
  } else {
    int k=X.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,X.size(1));
      P=OPERATOR_NEW Matrix<float,complex<float> >(k,m);
      F77NAME(clacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(ctrmm)('R','L',trans,'N',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(cgemm)('N',trans,k,m-n,n,one_,X.addr(),k,addr(n,0),m,
          zero_,P->addr(0,n),k);
      }
    } else {
      CHECK_SAME(m,X.size(1));
      P=OPERATOR_NEW Matrix<float,complex<float> >(k,n);
      F77NAME(clacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(ctrmm)('R','L','N','N',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(cgemm)('N','N',k,n,m-n,one_,X.addr(0,n),k,addr(n,0),m,
          one_,P->addr(),k);
      }
    }
  }
  return P;
}

template<> void LowerTrapezoidalMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(ccopy)(n,b.addr(),1,x.addr(),1);
    if (m>n) { // use trailing entries of x as free variables
      F77NAME(cgemv)(trans,m-n,n,mone_,addr(n,0),m,x.addr(n),1,one_,
        x.addr(),1);
    }
    F77NAME(ctrsv)('L',trans,'N',n,addr(),m,x.addr(),1);
  } else { // calling routine will have to check consistency conditions
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(ccopy)(n,b.addr(),1,x.addr(),1);
    F77NAME(ctrsv)('L','N','N',n,addr(),m,x.addr(),1);
  }
}

template<> void LowerTrapezoidalMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(clacpy)('A',n,k,B.addr(),n,X.addr(),n);
      if (m>n) { // use these entries of X as free variables
        F77NAME(cgemm)(trans,'N',n,k,m-n,mone_,addr(n,0),m,X.addr(n,0),m,
          one_,X.addr(),m);
      }
      F77NAME(ctrsm)('L','L',trans,'N',n,k,one_,addr(),m,X.addr(),m);
    } else { // calling routine will have to check consistency conditions
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(clacpy)('A',n,k,B.addr(),m,X.addr(),n);
      F77NAME(ctrsm)('L','L','N','N',n,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(clacpy)('A',k,n,B.addr(),k,X.addr(),k);
      F77NAME(ctrsm)('R','L',trans,'N',k,n,one_,addr(),m,X.addr(),k);
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(clacpy)('A',k,n,B.addr(),k,X.addr(),k);
      if (m>n) { // use these entries of X as free variables
        F77NAME(cgemm)('N','N',k,n,m-n,mone_,X.addr(0,n),n,addr(n,0),m,
          one_,X.addr(),k);
      }
      F77NAME(ctrsm)('R','L','N','N',k,n,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> float
LowerTrapezoidalMatrix<float,complex<float> >::normFrobenius() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('F','L','N',m,size(1),addr(),m,work);
}

template<> float
LowerTrapezoidalMatrix<float,complex<float> >::normInfinity() const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float result=F77NAME(clantr)('I','L','N',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> float
LowerTrapezoidalMatrix<float,complex<float> >::normMaxEntry() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('M','L','N',m,size(1),addr(),m,work);
}

template<> float
LowerTrapezoidalMatrix<float,complex<float> >::normOne() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('O','L','N',m,size(1),addr(),m,work);
}

template<> Matrix<float,complex<float> >* operator+(
const Matrix<float,complex<float> > &M,
const LowerTrapezoidalMatrix<float,complex<float> > &L) {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *S=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  S->copy(M);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m-j,complex_float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(caxpy)(m-j-1,complex_float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const LowerTrapezoidalMatrix<float,complex<float> > &L) {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *D=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  D->copy(M);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m-j,complex_float_mone_,L.addr(j,j),1,
        D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(caxpy)(m-j-1,complex_float_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
      (*D)(j,j)-=complex_float_one_;
    }
  }
  return D;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const LowerTrapezoidalMatrix<float,complex<float> > &L) {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  char diagL=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0 ?
    'N' : 'U');
  int m=M.size(0),k=L.size(0),n=L.size(1);
  CHECK_SAME(k,M.size(1));
  Matrix<float,complex<float> > *P=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(clacpy)('A',m,n,M.addr(),m,P->addr(),m);
  F77NAME(ctrmm)('R','L','N',diagL,m,n,complex_float_one_,L.addr(),k,
    P->addr(),m);
  if (k>n) {
    F77NAME(cgemm)('N','N',m,n,k-n,complex_float_one_,M.addr(0,n),m,
      L.addr(n,0),k,complex_float_one_,P->addr(),m);
  }
  return P;
}

template class LowerTrapezoidalMatrix<float,complex<float> >;
template void testLowerTrapezoidalMatrix(float,complex<float> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float
LowerTriangularMatrix<float,complex<float> >::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info;
  F77NAME(ctrcon)(norm,'L','N',n,addr(),n,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork;
  delete [] work;
  return rcond;
}

/*
template<> LowerTriangularMatrix<float,complex<float> >*
LowerTriangularMatrix<float,complex<float> >::inverse() const {
  LowerTriangularMatrix<float,complex<float> > *L=
    OPERATOR_NEW LowerTriangularMatrix<float,complex<float> >(*this);
  int n=size(0);
  int info;
  F77NAME(ctrtri)('L','N',n,L->addr(),n,info);
  CHECK_TEST(info==0);
  return L;
}
*/

template class LowerTriangularMatrix<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void
UnitLowerTrapezoidalMatrix<float,complex<float> >::copyFrom(int m,int n,
const Matrix<float,complex<float> > &L) {
  m=min(m,min(size(0),L.size(0)));
  n=min(n,min(size(1),L.size(1)));
  for (int j=0;j<n;j++) {
    if (j+1<m) {
      F77NAME(ccopy)(m-j-1,L.addr(j+1,j),1,addr(j+1,j),1);
    }
  }
}

/*
template<> UnitUpperTrapezoidalMatrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::transpose() const {
  int m=size(0),n=size(1);
  UnitUpperTrapezoidalMatrix<float,complex<float> > *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<float,complex<float> >(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(m-j-1,addr(j+1,j),1,U->addr(j,j+1),n);
  }
  return U;
}

template<> UnitUpperTrapezoidalMatrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::conjugateTranspose()
const {
  int m=size(0),n=size(1);
  UnitUpperTrapezoidalMatrix<float,complex<float> > *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<float,complex<float> >(n,m);
  for (int j=0;j<n;j++) {
    const complex<float> *col_j=addr(j+1,j);
    complex<float> *U_row_j=U->addr(j,j+1);
    for (int i=j+1;i<m;i++,col_j++,U_row_j+=n) {
      *U_row_j=conj(*col_j);
    }
  }
  return U;
}
*/

template<> LowerTrapezoidalMatrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<float,complex<float> > *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(m,n);
  if (L_non_unit) {
    S->copy(L);
    for (int j=0;j<n;j++) {
      (*S)(j,j)+=complex_float_one_;
      F77NAME(caxpy)(m-j-1,complex_float_one_,addr(j+1,j),1,
        S->addr(j+1,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=2.;
      if (j<m-1) {
        F77NAME(ccopy)(m-j-1,addr(j+1,j),1,S->addr(j+1,j),1);
        F77NAME(caxpy)(m-j-1,complex_float_one_,addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *S=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    (*S)(j,j)+=complex_float_one_;
    F77NAME(caxpy)(m-j-1,complex_float_one_,addr(j+1,j),1,
      S->addr(j+1,j),1);
  }
  return S;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<float,complex<float> > *dif=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(m,n);
  for (int j=0;j<n;j++) {
    (*dif)(j,j)=
      (L_non_unit ? complex_float_one_-L(j,j) : complex_float_zero_);
    if (j<m-1) {
      F77NAME(ccopy)(m-j-1,addr(j+1,j),1,dif->addr(j+1,j),1);
      F77NAME(caxpy)(m-j-1,complex_float_mone_,L.addr(j+1,j),1,
        dif->addr(j+1,j),1);
    }
  }
  return dif;
}

template<> Matrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *dif=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(claset)('U',m,n,complex_float_zero_,complex_float_zero_,
    dif->addr(),m);
  for (int j=0;j<n;j++) {
    (*dif)(j,j)=complex_float_one_-M(j,j);
    F77NAME(ccopy)(m-j-1,addr(j+1,j),1,dif->addr(j+1,j),1);
    F77NAME(caxpy)(m-j-1,complex_float_mone_,M.addr(j+1,j),1,
      dif->addr(j+1,j),1);
  }
  return dif;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator*(
complex<float> d) const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<float,complex<float> > *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(m,n);
  for (int j=0;j<n;j++) {
    (*P)(j,j)=d;
    if (j<m-1) {
      F77NAME(ccopy)(m-j-1,addr(j+1,j),1,P->addr(j+1,j),1);
      F77NAME(cscal)(m-j-1,d,P->addr(j+1,j),1);
    }
  }
  return P;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator/(
complex<float> d) const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<float,complex<float> > *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(m,n);
  complex<float> dinv=complex_float_one_/d;
  for (int j=0;j<n;j++) {
    (*P)(j,j)=dinv;
    if (j<m-1) {
      F77NAME(ccopy)(m-j-1,addr(j+1,j),1,P->addr(j+1,j),1);
      F77NAME(cscal)(m-j-1,dinv,P->addr(j+1,j),1);
    }
  }
  return P;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
// to compute the jth column of the product
//   [ L_11   0  ] [ 0 ] = [    0   ]
//   [ L_21 L_22 ] [ m ] = [ L_22 m ]
//   [ L_31 L_32 ]       = [ L_32 m ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  LowerTrapezoidalMatrix<float,complex<float> > *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(m,n);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(k-j,L.addr(j,j),1,P->addr(j,j),1);
      F77NAME(ctrmv)('L','N','U',k-j,addr(j,j),m,P->addr(j,j),1);
      if (m>k) {
        F77NAME(cgemv)('N',m-k,k-j,complex_float_one_,addr(k,j),m,
          L.addr(j,j),1,complex_float_zero_,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*P)(j,j)=complex_float_one_;
      if (j<k-1) {
        F77NAME(ccopy)(k-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
        F77NAME(ctrmv)('L','N','U',k-j-1,addr(j+1,j+1),m,
          P->addr(j+1,j),1);
        F77NAME(caxpy)(k-j-1,complex_float_one_,addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
      if (m>k) {
        if (j<k-1) {
          F77NAME(cgemv)('N',m-k,k-j-1,complex_float_one_,addr(k,j+1),m,
            L.addr(j+1,j),1,complex_float_zero_,P->addr(k,j),1);
        }
        F77NAME(caxpy)(m-k,complex_float_one_,addr(k,j),1,
          P->addr(k,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<float,complex<float> > *P=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n,zero_);
  F77NAME(clacpy)('A',k,n,M.addr(),k,P->addr(),m);
  F77NAME(ctrmm)('L','L','N','U',k,n,complex_float_one_,addr(),m,
    P->addr(),m);
  if (m>k) {
    F77NAME(cgemm)('N','N',m-k,n,k,complex_float_one_,addr(k,0),m,
      M.addr(),k,complex_float_one_,P->addr(k,0),m);
  }
  return P;
}

template<> Vector<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<float,complex<float> > *p=
    OPERATOR_NEW Vector<float,complex<float> >(m,zero_);
  F77NAME(ccopy)(n,v.addr(),1,p->addr(),1);
  F77NAME(ctrmv)('L','N','U',n,addr(),m,p->addr(),1);
  if (m>n) {
    F77NAME(cgemv)('N',m-n,n,complex_float_one_,addr(n,0),m,
      v.addr(),1,complex_float_one_,p->addr(n),m);
  }
  return p;
}

template<> Vector<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::trmv(
const Vector<float,complex<float> > &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<float,complex<float> > *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<float,complex<float> >(n);
    F77NAME(ccopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(ctrmv)('L',trans,'U',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(cgemv)(trans,m-n,n,one_,addr(n,0),m,x.addr(n),1,one_,
        p->addr(),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<float,complex<float> >(m);
    F77NAME(ccopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(ctrmv)('L','N','U',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(cgemv)('N',m-n,n,one_,addr(n,0),m,x.addr(),1,zero_,
        p->addr(n),1);
    }
  }
  return p;
}

template<> Matrix<float,complex<float> >*
UnitLowerTrapezoidalMatrix<float,complex<float> >::trmm(
const Matrix<float,complex<float> > &X,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<float,complex<float> > *P=0;
  if (side=='L' || side=='l') {
    int k=X.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,X.size(0));
      P=OPERATOR_NEW Matrix<float,complex<float> >(n,k);
      F77NAME(clacpy)('A',n,k,X.addr(),m,P->addr(),n);
      F77NAME(ctrmm)('L','L',trans,'U',n,k,one_,addr(),m,P->addr(),n);
      if (m>n) {
        F77NAME(cgemm)(trans,'N',n,k,m-n,one_,addr(n,0),m,X.addr(n,0),m,
          one_,P->addr(),n);
      }
    } else {
      CHECK_SAME(n,X.size(0));
      P=OPERATOR_NEW Matrix<float,complex<float> >(m,k);
      F77NAME(clacpy)('A',n,k,X.addr(),n,P->addr(),m);
      F77NAME(ctrmm)('L','L','N','U',n,k,one_,addr(),m,P->addr(),m);
      if (m>n) {
        F77NAME(cgemm)('N','N',m-n,k,n,one_,addr(n,0),m,X.addr(),n,
          zero_,P->addr(n,0),m);
      }
    }
  } else {
    int k=X.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,X.size(1));
      P=OPERATOR_NEW Matrix<float,complex<float> >(k,m);
      F77NAME(clacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(ctrmm)('R','L',trans,'U',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(cgemm)('N',trans,k,m-n,n,one_,X.addr(),k,addr(n,0),m,
          zero_,P->addr(0,n),k);
      }
    } else {
      CHECK_SAME(m,X.size(1));
      P=OPERATOR_NEW Matrix<float,complex<float> >(k,n);
      F77NAME(clacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(ctrmm)('R','L','N','U',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(cgemm)('N','N',k,n,m-n,one_,X.addr(0,n),k,addr(n,0),m,
          one_,P->addr(),k);
      }
    }
  }
  return P;
}

template<> void
UnitLowerTrapezoidalMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(ccopy)(n,b.addr(),1,x.addr(),1);
    if (m>n) { // use trailing entries of x as free variables
      F77NAME(cgemv)(trans,m-n,n,mone_,addr(n,0),m,x.addr(n),1,one_,
        x.addr(),1);
    }
    F77NAME(ctrsv)('L',trans,'U',n,addr(),m,x.addr(),1);
  } else { // calling routine will have to check consistency conditions
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(ccopy)(n,b.addr(),1,x.addr(),1);
    F77NAME(ctrsv)('L','N','U',n,addr(),m,x.addr(),1);
  }
}

template<> void
UnitLowerTrapezoidalMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(clacpy)('A',n,k,B.addr(),n,X.addr(),n);
      if (m>n) { // use these entries of X as free variables
        F77NAME(cgemm)(trans,'N',n,k,m-n,mone_,addr(n,0),m,X.addr(n,0),m,
          one_,X.addr(),m);
      }
      F77NAME(ctrsm)('L','L',trans,'U',n,k,one_,addr(),m,X.addr(),m);
    } else { // calling routine will have to check consistency conditions
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(clacpy)('A',n,k,B.addr(),m,X.addr(),n);
      F77NAME(ctrsm)('L','L','N','U',n,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(clacpy)('A',k,n,B.addr(),k,X.addr(),k);
      F77NAME(ctrsm)('R','L',trans,'U',k,n,one_,addr(),m,X.addr(),k);
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(clacpy)('A',k,n,B.addr(),k,X.addr(),k);
      if (m>n) { // use these entries of X as free variables
        F77NAME(cgemm)('N','N',k,n,m-n,mone_,X.addr(0,n),n,addr(n,0),m,
          one_,X.addr(),k);
      }
      F77NAME(ctrsm)('R','L','N','U',k,n,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> float
UnitLowerTrapezoidalMatrix<float,complex<float> >::normFrobenius()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('F','L','U',m,size(1),addr(),m,work);
}

template<> float
UnitLowerTrapezoidalMatrix<float,complex<float> >::normInfinity()
const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float result=F77NAME(clantr)('I','L','U',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> float
UnitLowerTrapezoidalMatrix<float,complex<float> >::normMaxEntry()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('M','L','U',m,size(1),addr(),m,work);
}

template<> float
UnitLowerTrapezoidalMatrix<float,complex<float> >::normOne() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('O','L','U',m,size(1),addr(),m,work);
}

template class UnitLowerTrapezoidalMatrix<float,complex<float> >;
template void testUnitLowerTrapezoidalMatrix(float,complex<float> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float UnitLowerTriangularMatrix<float,complex<float> >
::reciprocalConditionNumber(char norm) const {
  int n=size(0);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info;
  F77NAME(ctrcon)(norm,'L','U',n,addr(),n,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork;
  delete [] work;
  return rcond;
}

/*
template<> UnitLowerTriangularMatrix<float,complex<float> >*
UnitLowerTriangularMatrix<float,complex<float> >::inverse() const {
  UnitLowerTriangularMatrix<float,complex<float> > *result=
    OPERATOR_NEW UnitLowerTriangularMatrix<float,complex<float> >(*this);
  int n=size(0);
  int info;
  F77NAME(ctrtri)('L','U',n,result->addr(),n,info);
  CHECK_TEST(info==0);
  return result;
}
*/

template class UnitLowerTriangularMatrix<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> Matrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::makeMatrix() const {
  int m=size(0),n=size(1);
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(this)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+1,m),addr(0,j),1,M->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j,m),addr(0,j),1,M->addr(0,j),1);
      if (j<m) (*M)(j,j)=complex_float_one_;
    }
  }
  return M;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >& 
UpperTrapezoidalMatrix<float,complex<float> >::operator+=(
const UpperTrapezoidalMatrix<float,complex<float> > &M) {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0))
  CHECK_SAME(n,M.size(1))
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(j+1,m),one_,M.addr(0,j),1,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >&
UpperTrapezoidalMatrix<float,complex<float> >::operator-=(
const UpperTrapezoidalMatrix<float,complex<float> > &M) {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0))
  CHECK_SAME(n,M.size(1))
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(j+1,m),mone_,M.addr(0,j),1,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >& 
UpperTrapezoidalMatrix<float,complex<float> >::operator*=(
complex<float> d) {
  int m=size(0),n=size(1);
  for (int j=0;j<n;j++) {
    F77NAME(cscal)(min(j+1,m),d,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >&
UpperTrapezoidalMatrix<float,complex<float> >::operator/=(
complex<float> d) {
  int m=size(0),n=size(1);
  complex<float> dinv=complex_float_one_/d;
  for (int j=0;j<n;j++) {
    F77NAME(cscal)(min(j+1,m),dinv,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<float,complex<float> > *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(m,n);
  if (U_non_unit) {
    S->copy(U);
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(min(j+1,m),complex_float_one_,addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      if (j>0) {
        F77NAME(caxpy)(min(j,m),complex_float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(m,zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(j+1,addr(0,j),1,S->addr(0,j),1);
      F77NAME(caxpy)(m-j,complex_float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(j+1,addr(0,j),1,S->addr(0,j),1);
      if (j<m-1) {
        F77NAME(caxpy)(m-j-1,complex_float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *S=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(m,j+1),complex_float_one_,addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<float,complex<float> > *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(m,n);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      F77NAME(caxpy)(min(j+1,m),complex_float_mone_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      if (j>0) {
        F77NAME(caxpy)(min(j,m),complex_float_mone_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)-=complex_float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *D=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(j+1,addr(0,j),1,D->addr(0,j),1);
      F77NAME(caxpy)(m-j,complex_float_mone_,L.addr(j,j),1,
        D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(j+1,addr(0,j),1,D->addr(0,j),1);
      if (j<m-1) {
        F77NAME(caxpy)(m-j-1,complex_float_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
      (*D)(j,j)-=complex_float_one_;
    }
  }
  return D;
}

template<> Matrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *D=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(j+1,addr(0,j),1,D->addr(0,j),1);
    F77NAME(caxpy)(m,complex_float_mone_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >* 
UpperTrapezoidalMatrix<float,complex<float> >::operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
// [ U_1 , U_2 ] [ V_11 , V_12 , V_13 ]
//               [   0  , V_22 , V_23 ]
//   = [ U_1 V_11 , U_1 V_12 + U_2 V_22 , U_1 V_13 + U_2 V_23 ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  int m=size(0),k=size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  UpperTrapezoidalMatrix<float,complex<float> > *P=OPERATOR_NEW
    UpperTrapezoidalMatrix<float,complex<float> >(m,n,
      complex_float_zero_);
  if (U_non_unit) {
    for (int j=0;j<m;j++) {
      F77NAME(ccopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ctrmv)('U','N','N',j+1,addr(),m,P->addr(0,j),1);
    }
    for (int j=m;j<k;j++) {
      F77NAME(ccopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ctrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
      F77NAME(cgemv)('N',m,j-m+1,complex_float_one_,addr(0,m),m,
        U.addr(m,j),1,complex_float_one_,P->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(ccopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(ctrmv)('U','N','N',j,addr(),m,P->addr(0,j),1);
        F77NAME(caxpy)(j,complex_float_one_,addr(0,j),1,P->addr(0,j),1);
      }
      (*P)(j,j)=(*this)(j,j);
    }
    for (int j=m;j<k;j++) {
      F77NAME(ccopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ctrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
      if (j>m) {
        F77NAME(cgemv)('N',m,j-m,complex_float_one_,addr(0,m),m,
          U.addr(m,j),1,complex_float_one_,P->addr(0,j),1);
      }
      F77NAME(caxpy)(m,complex_float_one_,addr(0,j),1,P->addr(0,j),1);
    }
  }
  for (int j=k;j<n;j++) {
    F77NAME(ccopy)(m,U.addr(0,j),1,P->addr(0,j),1);
    F77NAME(ctrmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
    F77NAME(cgemv)('N',m,k-m,complex_float_one_,addr(0,m),m,
      U.addr(m,j),1,complex_float_one_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  if (m>=n) {
//  note that
//  [ U_11 U_12 U_13 ] [ L_1 ] = [ U_11 L_1 + U_12 L_2 + U_13 L_3 ]
//  [      U_22 U_23 ] [ L_2 ] = [            U_22 L_2 + U_23 L_3 ]
//                     [ L_3 ]
    if (L_non_unit) { // U_11 L_1:
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(cgemv)('N',j,n-j,complex_float_one_,addr(0,j),m,
            L.addr(j,j),1,complex_float_zero_,P->addr(0,j),1);
        }
        F77NAME(ccopy)(n-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(ctrmv)('U','N','N',n-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(caxpy)(j,complex_float_one_,addr(0,j),1,
            P->addr(0,j),1);
          if (j<n-1) {
            F77NAME(cgemv)('N',j,n-j-1,complex_float_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,complex_float_zero_,P->addr(0,j),1);
          }
        }
        if (j<n-1) {
          F77NAME(ccopy)(n-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(ctrmv)('U','N','N',n-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(cdotu)(n-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (n<k) { // U_12 L_2 + U_13 L_3:
      F77NAME(cgemm)('N','N',n,n,k-n,float_one_,addr(0,n),m,
        L.addr(n,0),k,float_one_,P->addr(0,0),m);
    }
    if (n<m) {
      F77NAME(clacpy)('A',m-n,n,L.addr(n,0),k,P->addr(n,0),m);
      F77NAME(ctrmm)('L','U','N','N',m-n,n,complex_float_one_,
        addr(n,n),m,P->addr(n,0),m); // U_22 L_2
      if (m<k) {
        F77NAME(cgemm)('N','N',m-n,n,k-m,complex_float_one_,addr(n,m),m,
          L.addr(m,0),k,complex_float_one_,P->addr(n,0),m); // U_23 L_3
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
          F77NAME(cgemv)('N',j,m-j,complex_float_one_,addr(0,j),m,
            L.addr(j,j),1,complex_float_zero_,P->addr(0,j),1);
        }
        F77NAME(ccopy)(m-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(ctrmv)('U','N','N',m-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(caxpy)(j,complex_float_one_,addr(0,j),1,
            P->addr(0,j),1);
          if (j<m-1) {
            F77NAME(cgemv)('N',j,m-j-1,complex_float_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,complex_float_zero_,P->addr(0,j),1);
          }
        }
        if (j<m-1) {
          F77NAME(ccopy)(m-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(ctrmv)('U','N','N',m-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(cdotu)(m-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (m<k) { // U_2 L_21 + U_3 L_31 
      F77NAME(cgemm)('N','N',m,m,k-m,complex_float_one_,addr(0,m),m,
        L.addr(m,0),k,complex_float_one_,P->addr(0,0),m);
    }
    char diagL=(L_non_unit ? 'N' : 'U');
    F77NAME(clacpy)('A',m,n-m,addr(0,m),m,P->addr(0,m),m);
    F77NAME(ctrmm)('R','L','N',diagL,m,n-m,complex_float_one_,
      L.addr(m,m),k,P->addr(0,m),m); // U_2 L_22
    if (n<k) {
      F77NAME(cgemm)('N','N',m,n-m,k-n,complex_float_one_,addr(0,n),m,
        L.addr(n,m),k,complex_float_one_,P->addr(0,m),m); // U_3 L_32
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
// [ U_1 U_2 ] [ M_1 ] = U_1 M_1 + U_2 M_2
//             [ M_2 ]
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<float,complex<float> > *P=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(clacpy)('A',m,n,M.addr(),k,P->addr(),m);
  F77NAME(ctrmm)('L','U','N','N',m,n,complex_float_one_,addr(),m,
    P->addr(),m); // U_1 M_1
  if (k>m) {
    F77NAME(cgemm)('N','N',m,n,k-m,complex_float_one_,addr(0,m),m,
      M.addr(m,0),k,complex_float_one_,P->addr(),m); // U_2 M_2
  }
  return P;
}

template<> Vector<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<float,complex<float> > *p=
    OPERATOR_NEW Vector<float,complex<float> >(m);
  F77NAME(ccopy)(m,v.addr(),1,p->addr(),1);
  F77NAME(ctrmv)('U','N','N',m,addr(),m,p->addr(),1);
  F77NAME(cgemv)('N',m,n-m,complex_float_one_,addr(0,m),m,v.addr(m),1,
    complex_float_one_,p->addr(),1);
  return p;
}

/*
template<> LowerTrapezoidalMatrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::transpose() const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<float,complex<float> > *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(n,m);
  for (int i=0;i<m;i++) {
    F77NAME(ccopy)(n-i,addr(i,i),m,L->addr(i,i),1);
  }
  return L;
}

template<> LowerTrapezoidalMatrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::conjugateTranspose()
const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<float,complex<float> > *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,complex<float> >(n,m);
  for (int i=0;i<m;i++) {
    const complex<float> *row_i=addr(i,i);
    complex<float> *L_col_i=L->addr(i,i);
    for (int j=i;j<n;j++,row_i+=m,L_col_i++) {
      *L_col_i=conj(*row_i);
    }
  }
  return L;
}
*/

template<> Vector<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::trmv(
const Vector<float,complex<float> > &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<float,complex<float> > *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<float,complex<float> >(n);
    F77NAME(ccopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(ctrmv)('U',trans,'N',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(cgemv)(trans,m,n-m,one_,addr(0,m),m,x.addr(),1,zero_,
        p->addr(m),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<float,complex<float> >(m);
    F77NAME(ccopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(ctrmv)('U','N','N',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(cgemv)('N',m,n-m,one_,addr(0,m),m,x.addr(m),1,one_,
        p->addr(),1);
    }
  }
  return p;
}

template<> Matrix<float,complex<float> >*
UpperTrapezoidalMatrix<float,complex<float> >::trmm(
const Matrix<float,complex<float> > &M,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<float,complex<float> > *P=0;
  if (side=='L' || side=='l') {
    int k=M.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,M.size(0));
      P=OPERATOR_NEW Matrix<float,complex<float> >(n,k);
      F77NAME(clacpy)('A',m,k,M.addr(),m,P->addr(),n);
      F77NAME(ctrmm)('L','U',trans,'N',m,k,one_,addr(),m,P->addr(),n);
      if (n>m) {
        F77NAME(cgemm)(trans,'N',n-m,k,m,one_,addr(0,m),m,M.addr(),m,
          zero_,P->addr(m,0),n);
      }
    } else {
      CHECK_SAME(n,M.size(0));
      P=OPERATOR_NEW Matrix<float,complex<float> >(m,k);
      F77NAME(clacpy)('A',m,k,M.addr(),n,P->addr(),m);
      F77NAME(ctrmm)('L','U','N','N',m,k,one_,addr(),m,P->addr(),m);
      if (n>m) {
        F77NAME(cgemm)('N','N',m,k,n-m,one_,addr(0,m),m,M.addr(m,0),n,
          one_,P->addr(),m);
      }
    }
  } else {
    int k=M.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,M.size(1));
      P=OPERATOR_NEW Matrix<float,complex<float> >(k,m);
      F77NAME(clacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(ctrmm)('R','U',trans,'N',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(cgemm)('N',trans,k,m,n-m,one_,M.addr(0,m),k,addr(0,m),m,
          one_,P->addr(),k);
      }
    } else {
      CHECK_SAME(m,M.size(1));
      P=OPERATOR_NEW Matrix<float,complex<float> >(k,n);
      F77NAME(clacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(ctrmm)('R','U','N','N',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(cgemm)('N','N',k,n-m,m,one_,M.addr(),k,addr(0,m),m,
          zero_,P->addr(0,m),k);
      }
    }
  }
  return P;
}

template<> float
UpperTrapezoidalMatrix<float,complex<float> >::normFrobenius() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('F','U','N',m,size(1),addr(),m,work);
}

template<> float
UpperTrapezoidalMatrix<float,complex<float> >::normInfinity() const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float result=F77NAME(clantr)('I','U','N',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> float
UpperTrapezoidalMatrix<float,complex<float> >::normMaxEntry() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('M','U','N',m,size(1),addr(),m,work);
}

template<> float
UpperTrapezoidalMatrix<float,complex<float> >::normOne() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('O','U','N',m,size(1),addr(),m,work);
}

template<> void UpperTrapezoidalMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  char diag=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(this) ==0 ?
    'N' : 'U');
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    // calling routine will have to check consistency conditions
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(ccopy)(m,b.addr(),1,x.addr(),1);
    F77NAME(ctrsv)('U',trans,'N',m,addr(),m,x.addr(),1);
  } else {
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(ccopy)(m,b.addr(),1,x.addr(),1);
    if (n>m) { // use trailing entries of x as free variables
      F77NAME(cgemv)('N',m,n-m,mone_,addr(0,m),m,x.addr(m),1,one_,
        x.addr(),1);
    }
    F77NAME(ctrsv)('U','N','N',m,addr(),m,x.addr(),1);
  }
}

template<> void UpperTrapezoidalMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  char diag=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(this) ==0 ?
    'N' : 'U');
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(clacpy)('A',m,k,B.addr(),n,X.addr(),m);
      F77NAME(ctrsm)('L','U',trans,'N',m,k,one_,addr(),m,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(clacpy)('A',m,k,B.addr(),m,X.addr(),n);
      if (n>m) { // use these entries of X as free variables
        F77NAME(cgemm)('N','N',m,k,n-m,mone_,addr(0,m),m,X.addr(m,0),n,
          one_,X.addr(),n);
      }
      F77NAME(ctrsm)('L','U','N','N',m,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(clacpy)('A',k,m,B.addr(),k,X.addr(),k);
      if (n>m) { // use these entries of X as free variables
        F77NAME(cgemm)('N',trans,k,m,n-m,mone_,X.addr(0,m),k,addr(0,m),m,
          one_,X.addr(),k);
      }
      F77NAME(ctrsm)('R','U',trans,'N',k,m,one_,addr(),m,X.addr(),k);
    } else { // calling routine will have to check consistency
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(clacpy)('A',k,m,B.addr(),k,X.addr(),k);
      F77NAME(ctrsm)('R','U','N','N',k,m,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> SquareMatrix<float,complex<float> >* operator+(
const UnitLowerTrapezoidalMatrix<float,complex<float> > &L,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    (*S)(j,j)=complex_float_one_;
    if (j<m-1) {
      F77NAME(ccopy)(m-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
    }
    if (U_non_unit) {
      F77NAME(caxpy)(j+1,complex_float_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    if (L_non_unit) {
      F77NAME(ccopy)(m-j,L.addr(j,j),1,S->addr(j,j),1);
    } else {
      (*S)(j,j)=complex_float_one_;
      if (j<m-1) {
        F77NAME(ccopy)(m-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
    }
    if (U_non_unit) {
      F77NAME(caxpy)(j+1,complex_float_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >* operator+(
const Matrix<float,complex<float> > &M,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(not_trapezoidal);
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  Matrix<float,complex<float> > *S=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  S->copy(M);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(min(j+1,m),complex_float_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j,m),complex_float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const UnitLowerTrapezoidalMatrix<float,complex<float> > &L,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  SquareMatrix<float,complex<float> > *D=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    (*D)(j,j)=complex_float_one_;
    if (j<m-1) {
      F77NAME(ccopy)(m-j-1,L.addr(j+1,j),1,D->addr(j+1,j),1);
    }
    if (U_non_unit) {
      F77NAME(caxpy)(j+1,complex_float_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_mone_,U.addr(0,j),1,
          D->addr(0,j),1);
      }
      (*D)(j,j)-=complex_float_one_;
    }
  }
  return D;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  SquareMatrix<float,complex<float> > *D=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    if (L_non_unit) {
      F77NAME(ccopy)(m-j,L.addr(j,j),1,D->addr(j,j),1);
    } else {
      (*D)(j,j)=complex_float_one_;
      if (j<m-1) {
        F77NAME(ccopy)(m-j-1,L.addr(j+1,j),1,D->addr(j+1,j),1);
      }
    }
    if (U_non_unit) {
      F77NAME(caxpy)(j+1,complex_float_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_mone_,U.addr(0,j),1,
          D->addr(0,j),1);
      }
      (*D)(j,j)-=complex_float_one_;
    }
  }
  return D;
}

template<> Matrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(not_trapezoidal);
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  Matrix<float,complex<float> > *D=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  D->copy(M);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(min(j+1,m),complex_float_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j,m),complex_float_mone_,U.addr(0,j),1,
          D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)-=complex_float_one_;
    }
  }
  return D;
}

template<> Matrix<float,complex<float> >* operator*(
const UnitLowerTrapezoidalMatrix<float,complex<float> > &L,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
//  Note that
//  [ L_1 ] [ U_1 U_2 ] = [ L_1 U_1 , L_1 U_2 ]
//  [ L_2 ]             = [ L_2 U_1 , L_2 U_2 ]
//  and that
//  [ L_11      ] [ u ] = [ L_11 u ]
//  [ L_21 L_22 ] [   ] = [ L_21 u ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U) ==0);
  int m=L.size(0),k=L.size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<k;j++) { // L_1 U_1
    if (U_non_unit) {
      F77NAME(ccopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ctrmv)('L','N','U',j+1,L.addr(),m,P->addr(0,j),1);
      if (j+1<k) {
        F77NAME(cgemv)('N',k-j-1,j+1,complex_float_one_,L.addr(j+1,0),m,
          U.addr(0,j),1,complex_float_zero_,P->addr(j+1,j),1);
      }
    } else {
      if (j>0) {
        F77NAME(ccopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(ctrmv)('L','N','U',j,L.addr(),m,P->addr(0,j),1);
        (*P)(j,j)=F77NAME(cdotu)(j,L.addr(j,0),m,U.addr(0,j),1);
        if (j+1<k) {
          F77NAME(cgemv)('N',k-j-1,j,complex_float_one_,L.addr(j+1,0),m,
            U.addr(0,j),1,complex_float_zero_,P->addr(j+1,j),1);
        }
      }
      (*P)(j,j)+=L(j,j);
      if (j+1<k) {
        F77NAME(caxpy)(k-j-1,complex_float_one_,L.addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
    }
  }
  if (n>k) { // L_1 U_2
    F77NAME(clacpy)('A',k,n-k,U.addr(0,k),k,P->addr(0,k),m);
    F77NAME(ctrmm)('L','L','N','U',k,n-k,complex_float_one_,L.addr(),m,
      P->addr(0,k),m);
  }
  if (m>k) {
    char diagU=(U_non_unit ? 'N' : 'U');
    F77NAME(clacpy)('A',m-k,k,L.addr(k,0),m,P->addr(k,0),m);
    F77NAME(ctrmm)('R','U','N',diagU,m-k,k,complex_float_one_,U.addr(),k,
      P->addr(k,0),m); // L_2 U_1
    if (n>k) { // L_2 U_2
      F77NAME(cgemm)('N','N',m-k,n-k,k,complex_float_one_,L.addr(k,0),m,
        U.addr(0,k),k,complex_float_zero_,P->addr(k,k),m);
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >* operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
//  Note that
//  [ L_1 ] [ U_1 U_2 ] = [ L_1 U_1 , L_1 U_2 ]
//  [ L_2 ]             = [ L_2 U_1 , L_2 U_2 ]
//  and that
//  [ L_11      ] [ u ] = [ L_11 u ]
//  [ L_21 L_22 ] [   ] = [ L_21 u ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  char diagL=(L_non_unit ? 'N' : 'U');
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  int m=L.size(0),k=L.size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<k;j++) { // L_1 U_1
    if (U_non_unit) {
      F77NAME(ccopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ctrmv)('L','N',diagL,j+1,L.addr(),m,P->addr(0,j),1);
      if (j+1<k) {
        F77NAME(cgemv)('N',k-j-1,j+1,complex_float_one_,L.addr(j+1,0),m,
          U.addr(0,j),1,complex_float_zero_,P->addr(j+1,0),1);
      }
    } else {
      if (j>0) {
        F77NAME(ccopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(ctrmv)('L','N',diagL,j,L.addr(),m,P->addr(0,j),1);
        (*P)(j,j)=F77NAME(cdotu)(j,L.addr(j,0),m,U.addr(0,j),1);
        if (j+1<k) {
          F77NAME(cgemv)('N',k-j-1,j,complex_float_one_,L.addr(j+1,0),m,
            U.addr(0,j),1,complex_float_zero_,P->addr(j+1,j),1);
        }
      }
      (*P)(j,j)+=(L_non_unit ? L(j,j) : complex_float_one_);
      if (j+1<k) {
        F77NAME(caxpy)(k-j-1,complex_float_one_,L.addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
    }
  }
  if (n>k) { // L_1 U_2
    F77NAME(clacpy)('A',k,n-k,U.addr(0,k),k,P->addr(0,k),m);
    F77NAME(ctrmm)('L','L','N',diagL,k,n-k,complex_float_one_,L.addr(),m,
      P->addr(0,k),m);
  }
  if (m>k) {
    char diagU=(U_non_unit ? 'N' : 'U');
    F77NAME(clacpy)('A',m-k,k,L.addr(k,0),m,P->addr(k,0),m);
    F77NAME(ctrmm)('R','U','N',diagU,m-k,k,complex_float_one_,U.addr(),k,
      P->addr(k,0),m); // L_2 U_1
    if (n>k) { // L_2 U_2
      F77NAME(cgemm)('N','N',m-k,n-k,k,complex_float_one_,L.addr(k,0),m,
        U.addr(0,k),k,complex_float_zero_,P->addr(k,k),m); 
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const UpperTrapezoidalMatrix<float,complex<float> > &U) {
// Note that
// M [ U_1 , U_2 ] = [ M U_1 , M U_2 ]
  char diagU=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U) ==0 ?
    'N' : 'U');
  int m=M.size(0),k=U.size(0),n=U.size(1);
  CHECK_SAME(k,M.size(1));
  Matrix<float,complex<float> > *P=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(clacpy)('A',m,k,M.addr(),m,P->addr(),m);
  F77NAME(ctrmm)('R','U','N',diagU,m,k,complex_float_one_,U.addr(),k,
    P->addr(),m); // M U_1
  if (n>k) { // M U_2
    F77NAME(cgemm)('N','N',m,n-k,k,complex_float_one_,M.addr(),m,
      U.addr(0,k),k,complex_float_zero_,P->addr(0,k),m);
  }
  return P;
}

template class UpperTrapezoidalMatrix<float,complex<float> >;
template void testUpperTrapezoidalMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float
UpperTriangularMatrix<float,complex<float> >::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info;
  F77NAME(ctrcon)(norm,'U','N',n,addr(),n,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork;
  delete [] work;
  return rcond;
}

/*
template<> UpperTriangularMatrix<float,complex<float> >*
UpperTriangularMatrix<float,complex<float> >::inverse() const {
  UpperTriangularMatrix<float,complex<float> > *I=
    OPERATOR_NEW UpperTriangularMatrix<float,complex<float> >(*this);
  int n=size(0);
  int info;
  F77NAME(ctrtri)('U','N',n,I->addr(),n,info);
  CHECK_TEST(info==0);
  return I;
}
*/

template class UpperTriangularMatrix<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void
UnitUpperTrapezoidalMatrix<float,complex<float> >::copyFrom(int m,int n,
const Matrix<float,complex<float> > &U) {
  m=min(m,min(size(0),U.size(0)));
  n=min(n,min(size(1),U.size(1)));
  for (int j=1;j<n;j++) {
    F77NAME(ccopy)(min(j,m),U.addr(0,j),1,addr(0,j),1);
  }
}

template<> UpperTrapezoidalMatrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<float,complex<float> > *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(m,n);
  if (U_non_unit) {
    S->copy(U);
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j,m),complex_float_one_,addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=complex_float_one_;
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(ccopy)(min(j,m),addr(0,j),1,S->addr(0,j),1);
        F77NAME(caxpy)(min(j,m),complex_float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)=2.;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(m,complex_float_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(ccopy)(j,addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)=complex_float_one_;
      F77NAME(caxpy)(m-j,complex_float_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(ccopy)(j,addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)=2.;
      if (m<j-1) {
        F77NAME(caxpy)(m-j-1,complex_float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *S=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(caxpy)(min(j,m),complex_float_one_,addr(0,j),1,
        S->addr(0,j),1);
    }
    if (j<m) (*S)(j,j)+=complex_float_one_;
  }
  return S;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<float,complex<float> > *D=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(m,n);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(ccopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)=complex_float_one_;
      F77NAME(caxpy)(min(j+1,m),complex_float_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(ccopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
        F77NAME(caxpy)(min(j,m),complex_float_mone_,addr(0,j),1,
          D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)=complex_float_zero_;
    }
  }
  return D;
}

template<> SquareMatrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *D=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(m,complex_float_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(ccopy)(j,addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)=complex_float_one_;
      F77NAME(caxpy)(m-j,complex_float_mone_,L.addr(j,j),1,
        D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(ccopy)(j,addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)=complex_float_zero_;
      if (m<j-1) {
        F77NAME(caxpy)(m-j-1,complex_float_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
    }
  }
  return D;
}

template<> Matrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *D=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(ccopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
    }
    if (j<m) (*D)(j,j)=complex_float_one_;
    F77NAME(caxpy)(m,complex_float_one_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator*(
complex<float> d) const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<float,complex<float> > *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(m,n);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(ccopy)(min(j,m),addr(0,j),1,P->addr(0,j),1);
      F77NAME(cscal)(min(j,m),d,P->addr(0,j),1);
    }
    if (j<m) (*P)(j,j)=d;
  }
  return P;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator/(
complex<float> d) const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<float,complex<float> > *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(m,n);
  complex<float> dinv=complex_float_one_/d;
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(ccopy)(min(j,m),addr(0,j),1,P->addr(0,j),1);
      F77NAME(cscal)(min(j,m),dinv,P->addr(0,j),1);
    }
    if (j<m) (*P)(j,j)=dinv;
  }
  return P;
}

template<> UpperTrapezoidalMatrix<float,complex<float> >* 
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
// [ U_1 , U_2 ] [ V_11 , V_12 , V_13 ]
//               [   0  , V_22 , V_23 ]
//   = [ U_1 V_11 , U_1 V_12 + U_2 V_22 , U_1 V_13 + U_2 V_23 ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  int m=size(0),k=size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  UpperTrapezoidalMatrix<float,complex<float> > *P=OPERATOR_NEW
    UpperTrapezoidalMatrix<float,complex<float> >(m,n,
      complex_float_zero_);
  if (U_non_unit) {
    for (int j=0;j<m;j++) {
      F77NAME(ccopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ctrmv)('U','N','U',j+1,addr(),m,P->addr(0,j),1);
    }
    for (int j=m;j<k;j++) {
      F77NAME(ccopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ctrmv)('U','N','U',m,addr(),m,P->addr(0,j),1);
      F77NAME(cgemv)('N',m,j-m+1,complex_float_one_,addr(0,m),m,
        U.addr(m,j),1,complex_float_one_,P->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(ccopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(ctrmv)('U','N','U',j,addr(),m,P->addr(0,j),1);
        F77NAME(caxpy)(j,complex_float_one_,addr(0,j),1,P->addr(0,j),1);
      }
      (*P)(j,j)=complex_float_one_;
    }
    for (int j=m;j<k;j++) {
      F77NAME(ccopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(ctrmv)('U','N','U',m,addr(),m,P->addr(0,j),1);
      if (j>m) {
        F77NAME(cgemv)('N',m,j-m,complex_float_one_,addr(0,m),m,
          U.addr(m,j),1,complex_float_one_,P->addr(0,j),1);
      }
      F77NAME(caxpy)(m,complex_float_one_,addr(0,j),1,P->addr(0,j),1);
    }
  }
  for (int j=k;j<n;j++) {
    F77NAME(ccopy)(m,U.addr(0,j),1,P->addr(0,j),1);
    F77NAME(ctrmv)('U','N','U',m,addr(),m,P->addr(0,j),1);
    F77NAME(cgemv)('N',m,k-m,complex_float_one_,addr(0,m),m,
      U.addr(m,j),1,complex_float_one_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  if (m>=n) {
//  note that
//  [ U_11 U_12 U_13 ] [ L_1 ] = [ U_11 L_1 + U_12 L_2 + U_13 L_3 ]
//  [      U_22 U_23 ] [ L_2 ] = [            U_22 L_2 + U_23 L_3 ]
//                     [ L_3 ]
    if (L_non_unit) { // U_11 L_1:
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(cgemv)('N',j,n-j,complex_float_one_,addr(0,j),m,
            L.addr(j,j),1,complex_float_zero_,P->addr(0,j),1);
        }
        F77NAME(ccopy)(n-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(ctrmv)('U','N','U',n-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(caxpy)(j,complex_float_one_,addr(0,j),1,
            P->addr(0,j),1);
          if (j<n-1) {
            F77NAME(cgemv)('N',j,n-j-1,complex_float_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,complex_float_zero_,P->addr(0,j),1);
          }
        }
        if (j<n-1) {
          F77NAME(ccopy)(n-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(ctrmv)('U','N','U',n-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(cdotu)(n-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=complex_float_one_;
      }
    }
    if (n<k) { // U_12 L_2 + U_13 L_3:
      F77NAME(cgemm)('N','N',n,n,k-n,float_one_,addr(0,n),m,
        L.addr(n,0),k,float_one_,P->addr(0,0),m);
    }
    if (n<m) {
      F77NAME(clacpy)('A',m-n,n,L.addr(n,0),k,P->addr(n,0),m);
      F77NAME(ctrmm)('L','U','N','U',m-n,n,complex_float_one_,
        addr(n,n),m,P->addr(n,0),m); // U_22 L_2
      if (m<k) {
        F77NAME(cgemm)('N','N',m-n,n,k-m,complex_float_one_,addr(n,m),m,
          L.addr(m,0),k,complex_float_one_,P->addr(n,0),m); // U_23 L_3
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
          F77NAME(cgemv)('N',j,m-j,complex_float_one_,addr(0,j),m,
            L.addr(j,j),1,complex_float_zero_,P->addr(0,j),1);
        }
        F77NAME(ccopy)(m-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(ctrmv)('U','N','U',m-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(caxpy)(j,complex_float_one_,addr(0,j),1,
            P->addr(0,j),1);
          if (j<m-1) {
            F77NAME(cgemv)('N',j,m-j-1,complex_float_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,complex_float_zero_,P->addr(0,j),1);
          }
        }
        if (j<m-1) {
          F77NAME(ccopy)(m-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(ctrmv)('U','N','U',m-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(cdotu)(m-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (m<k) { // U_2 L_21 + U_3 L_31 
      F77NAME(cgemm)('N','N',m,m,k-m,complex_float_one_,addr(0,m),m,
        L.addr(m,0),k,complex_float_one_,P->addr(0,0),m);
    }
    char diagL=(L_non_unit ? 'N' : 'U');
    F77NAME(clacpy)('A',m,n-m,addr(0,m),m,P->addr(0,m),m);
    F77NAME(ctrmm)('R','L','U',diagL,m,n-m,complex_float_one_,
      L.addr(m,m),k,P->addr(0,m),m); // U_2 L_22
    if (n<k) {
      F77NAME(cgemm)('N','N',m,n-m,k-n,complex_float_one_,addr(0,n),m,
        L.addr(n,m),k,complex_float_one_,P->addr(0,m),m); // U_3 L_32
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
// [ U_1 U_2 ] [ M_1 ] = U_1 M_1 + U_2 M_2
//             [ M_2 ]
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<float,complex<float> > *P=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(clacpy)('A',m,n,M.addr(),k,P->addr(),m);
  F77NAME(ctrmm)('L','U','N','U',m,n,complex_float_one_,addr(),m,
    P->addr(),m); // U_1 M_1
  if (k>m) {
    F77NAME(cgemm)('N','N',m,n,k-m,complex_float_one_,addr(0,m),m,
      M.addr(m,0),k,complex_float_one_,P->addr(),m); // U_2 M_2
  }
  return P;
}

template<> Vector<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<float,complex<float> > *p=
    OPERATOR_NEW Vector<float,complex<float> >(m);
  F77NAME(ccopy)(m,v.addr(),1,p->addr(),1);
  F77NAME(ctrmv)('U','N','U',m,addr(),m,p->addr(),1);
  F77NAME(cgemv)('N',m,n-m,complex_float_one_,addr(0,m),m,v.addr(m),1,
    complex_float_one_,p->addr(),1);
  return p;
}

/*
template<> UnitLowerTrapezoidalMatrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::transpose() const {
  int m=size(0),n=size(1);
  UnitLowerTrapezoidalMatrix<float,complex<float> > *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<float,complex<float> >(n,m);
  for (int i=0;i<m;i++) {
    F77NAME(ccopy)(n-i-1,addr(i,i+1),m,L->addr(i+1,i),1);
  }
  return L;
}

template<> UnitLowerTrapezoidalMatrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::conjugateTranspose()
const {
  int m=size(0),n=size(1);
  UnitLowerTrapezoidalMatrix<float,complex<float> > *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<float,complex<float> >(n,m);
  for (int i=0;i<m;i++) {
    const complex<float> *row_i=addr(i,i+1);
    complex<float> *L_col_i=L->addr(i+1,i);
    for (int j=i+1;j<n;j++,row_i+=m,L_col_i++) {
      *L_col_i=conj(*row_i);
    }
  }
  return L;
}
*/

template<> Vector<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::trmv(
const Vector<float,complex<float> > &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<float,complex<float> > *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<float,complex<float> >(n);
    F77NAME(ccopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(ctrmv)('U',trans,'U',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(cgemv)(trans,m,n-m,one_,addr(0,m),m,x.addr(),1,zero_,
        p->addr(m),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<float,complex<float> >(m);
    F77NAME(ccopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(ctrmv)('U','N','U',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(cgemv)('N',m,n-m,one_,addr(0,m),m,x.addr(m),1,one_,
        p->addr(),1);
    }
  }
  return p;
}

template<> Matrix<float,complex<float> >*
UnitUpperTrapezoidalMatrix<float,complex<float> >::trmm(
const Matrix<float,complex<float> > &M,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<float,complex<float> > *P=0;
  if (side=='L' || side=='l') {
    int k=M.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,M.size(0));
      P=OPERATOR_NEW Matrix<float,complex<float> >(n,k);
      F77NAME(clacpy)('A',m,k,M.addr(),m,P->addr(),n);
      F77NAME(ctrmm)('L','U',trans,'U',m,k,one_,addr(),m,P->addr(),n);
      if (n>m) {
        F77NAME(cgemm)(trans,'N',n-m,k,m,one_,addr(0,m),m,M.addr(),m,
          zero_,P->addr(m,0),n);
      }
    } else {
      CHECK_SAME(n,M.size(0));
      P=OPERATOR_NEW Matrix<float,complex<float> >(m,k);
      F77NAME(clacpy)('A',m,k,M.addr(),n,P->addr(),m);
      F77NAME(ctrmm)('L','U','N','U',m,k,one_,addr(),m,P->addr(),m);
      if (n>m) {
        F77NAME(cgemm)('N','N',m,k,n-m,one_,addr(0,m),m,M.addr(m,0),n,
          one_,P->addr(),m);
      }
    }
  } else {
    int k=M.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,M.size(1));
      P=OPERATOR_NEW Matrix<float,complex<float> >(k,m);
      F77NAME(clacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(ctrmm)('R','U',trans,'U',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(cgemm)('N',trans,k,m,n-m,one_,M.addr(0,m),k,addr(0,m),m,
          one_,P->addr(),k);
      }
    } else {
      CHECK_SAME(m,M.size(1));
      P=OPERATOR_NEW Matrix<float,complex<float> >(k,n);
      F77NAME(clacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(ctrmm)('R','U','N','U',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(cgemm)('N','N',k,n-m,m,one_,M.addr(),k,addr(0,m),m,
          zero_,P->addr(0,m),k);
      }
    }
  }
  return P;
}

template<> float
UnitUpperTrapezoidalMatrix<float,complex<float> >::normFrobenius()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('F','U','U',m,size(1),addr(),m,work);
}

template<> float
UnitUpperTrapezoidalMatrix<float,complex<float> >::normInfinity()
const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float result=F77NAME(clantr)('I','U','U',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> float
UnitUpperTrapezoidalMatrix<float,complex<float> >::normMaxEntry()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('M','U','U',m,size(1),addr(),m,work);
}

template<> float
UnitUpperTrapezoidalMatrix<float,complex<float> >::normOne() const {
  int m=size(0);
  float *work=0;
  return F77NAME(clantr)('O','U','U',m,size(1),addr(),m,work);
}

template<> void
UnitUpperTrapezoidalMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    // calling routine will have to check consistency conditions
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(ccopy)(m,b.addr(),1,x.addr(),1);
    F77NAME(ctrsv)('U',trans,'U',m,addr(),m,x.addr(),1);
  } else {
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(ccopy)(m,b.addr(),1,x.addr(),1);
    if (n>m) { // use trailing entries of x as free variables
      F77NAME(cgemv)('N',m,n-m,mone_,addr(0,m),m,x.addr(m),1,one_,
        x.addr(),1);
    }
    F77NAME(ctrsv)('U','N','U',m,addr(),m,x.addr(),1);
  }
}

template<> void
UnitUpperTrapezoidalMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  bool not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,complex<float> >*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(clacpy)('A',m,k,B.addr(),n,X.addr(),m);
      F77NAME(ctrsm)('L','U',trans,'U',m,k,one_,addr(),m,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(clacpy)('A',m,k,B.addr(),m,X.addr(),n);
      if (n>m) { // use these entries of X as free variables
        F77NAME(cgemm)('N','N',m,k,n-m,mone_,addr(0,m),m,X.addr(m,0),n,
          one_,X.addr(),n);
      }
      F77NAME(ctrsm)('L','U','N','U',m,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(clacpy)('A',k,m,B.addr(),k,X.addr(),k);
      if (n>m) { // use these entries of X as free variables
        F77NAME(cgemm)('N',trans,k,m,n-m,mone_,X.addr(0,m),k,addr(0,m),m,
          one_,X.addr(),k);
      }
      F77NAME(ctrsm)('R','U',trans,'U',k,m,one_,addr(),m,X.addr(),k);
    } else { // calling routine will have to check consistency
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(clacpy)('A',k,m,B.addr(),k,X.addr(),k);
      F77NAME(ctrsm)('R','U','N','U',k,m,one_,addr(),m,X.addr(),k);
    }
  }
}

template class UnitUpperTrapezoidalMatrix<float,complex<float> >;
template void testUnitUpperTrapezoidalMatrix(float,complex<float> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float
UnitUpperTriangularMatrix<float,complex<float> >::
reciprocalConditionNumber(char norm) const {
  int n=size(0);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info;
  F77NAME(ctrcon)(norm,'U','U',n,addr(),n,rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork;
  delete [] work;
  return rcond;
}

/*
template<> UnitUpperTriangularMatrix<float,complex<float> >*
UnitUpperTriangularMatrix<float,complex<float> >::inverse() const {
  UnitUpperTriangularMatrix<float,complex<float> > *I=
    OPERATOR_NEW UnitUpperTriangularMatrix<float,complex<float> >(*this);
  int n=size(0);
  int info;
  F77NAME(ctrtri)('U','U',n,I->addr(),n,info);
  CHECK_TEST(info==0);
  return I;
}
*/

template class UnitUpperTriangularMatrix<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "OrthogonalMatrix.C"

/*
template<> OrthogonalMatrix<float,complex<float> >*
OrthogonalMatrix<float,complex<float> >::transpose() const {
  int m=size(0),n=size(1);
  OrthogonalMatrix<float,complex<float> > *Q=
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(m,addr(0,j),1,Q->addr(j,0),n);
  }
  return Q;
}

template<> OrthogonalMatrix<float,complex<float> >*
OrthogonalMatrix<float,complex<float> >::conjugateTranspose() const {
  int m=size(0),n=size(1);
  OrthogonalMatrix<float,complex<float> > *Q=
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(n,m);
  for (int j=0;j<n;j++) {
    const complex<float> *col_j=addr(0,j);
    complex<float> *Q_row_j=Q->addr(j,0);
    for (int i=0;i<m;i++,col_j++,Q_row_j+=n) {
      *Q_row_j=conj(*col_j);
    }
  }
  return Q;
}
*/

/*
template<> void OrthogonalMatrix<float,complex<float> >::solve(
const UpperTrapezoidalMatrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  int m=size(0), n=size(1);
  bool B_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&B)==0);
  X=complex_float_zero_;
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<=min(j,n-1);i++) {
            F77NAME(caxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<min(j,n);i++) {
            F77NAME(caxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
          if (j<n) {
            F77NAME(caxpy)(m,complex_float_one_,addr(0,j),1,
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
            X(i,j)=F77NAME(cdotc)(min(j+1,m),addr(0,i),1,B.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(cdotc)(min(j,m),addr(0,i),1,B.addr(0,j),1);
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
            X(i,j)=F77NAME(cdotu)(m-i,B.addr(i,i),k,addr(i,j),1);
          }
        }
      } else {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=(*this)(i,j)
              +F77NAME(cdotu)(m-i-1,B.addr(i,i+1),k,addr(i+1,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      if (B_non_unit) {
        for (int i=0;i<k;i++) {
          for (int j=i;j<n;j++) {
            F77NAME(caxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      } else {
        for (int i=0;i<k;i++) {
          F77NAME(caxpy)(m,float_one_,addr(0,i),1,X.addr(i,0),k);
          for (int j=i+1;j<n;j++) {
            F77NAME(caxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      }
    }
  }
}
*/

/*
template<> void OrthogonalMatrix<float,complex<float> >::solve(
const LowerTrapezoidalMatrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  int m=size(0), n=size(1);
  bool B_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&B)==0);
  X=complex_float_zero_;
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=j;i<n;i++) {
            F77NAME(caxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          F77NAME(caxpy)(m,complex_float_one_,addr(0,j),1,X.addr(0,j),1);
          for (int i=j+1;i<n;i++) {
            F77NAME(caxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(cdotc)(m-j,addr(j,i),1,B.addr(j,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=(*this)(j,i)
              +F77NAME(cdotc)(m-j-1,addr(j+1,i),1,B.addr(j+1,j),1);
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
            X(i,j)=F77NAME(cdotu)(min(i+1,m),B.addr(i,0),k,addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=F77NAME(cdotu)(min(i,m),B.addr(i,0),k,addr(0,j),1);
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
            F77NAME(caxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      } else {
        for (int i=0;i<k;i++) {
          for (int j=0;j<min(i,n);j++) {
            F77NAME(caxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
          if (i<n) {
            F77NAME(caxpy)(m,complex_float_one_,addr(0,i),1,
              X.addr(i,0),k);
          }
        }
      }
    }
  }
}
*/

template<> void OrthogonalMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  int m=size(0), n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans=='C' || trans=='c') { // X=Q*B is smallest solution
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(cgemm)('N','N',m,k,n,float_one_,addr(),m,B.addr(),n,
        complex_float_zero_,X.addr(),m);
    } else if (trans=='T' || trans=='t') {
    // X=conj(Q)*B=conj(Q*conj(B)) is smallest
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      X=complex_float_zero_;
      for (int kk=0;kk<k;kk++) {
        for (int j=0;j<n;j++) {
          F77NAME(caxpy)(m,conj(B(j,kk)),addr(0,j),1,X.addr(0,kk),1);
        }
        complex<float> *Xik=X.addr(0,kk);
        for (int i=0;i<m;i++,Xik++) *Xik=conj(*Xik);
      }
    } else { // X=Q^H*B minimizes residual
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
//    F77NAME(cgemm)('C','N',n,k,m,float_one_,addr(),m,B.addr(),m,
//      complex_float_zero_,X.addr(),n);
      Vector<float,complex<float> > *r=
        OPERATOR_NEW Vector<float,complex<float> >(m);
      for (int l=0;l<k;l++) {
        F77NAME(ccopy)(m,B.addr(0,l),1,r->addr(),1);
        for (int j=0;j<n;j++) {
          complex<float> &Xjl=X(j,l);
          Xjl=F77NAME(cdotc)(m,addr(0,j),1,r->addr(),1);
          F77NAME(caxpy)(m,-Xjl,addr(0,j),1,r->addr(),1);
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
//    F77NAME(cgemm)('N','N',k,n,m,float_one_,B.addr(),k,addr(),m,
//      complex_float_zero_,X.addr(),k);
      Vector<float,complex<float> > *r=
        OPERATOR_NEW Vector<float,complex<float> >(m);
      for (int l=0;l<k;l++) {
        for (int i=0;i<m;i++) (*r)[i]=conj(B(l,i));
        for (int j=0;j<n;j++) {
          complex<float> Xlj=F77NAME(cdotc)(m,addr(0,j),1,r->addr(),1);
          F77NAME(caxpy)(m,-Xlj,addr(0,j),1,r->addr(),1);
          X(l,j)=conj(Xlj);
        }
      }
      delete r; r=0;
    } else if (trans=='T' || trans=='t') {
    // X=B*conj(Q) minimizes resid
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
//    X=complex_float_zero_;
//    for (int j=0;j<n;j++) {
//      for (int i=0;i<m;i++) {
//        F77NAME(caxpy)(k,conj((*this)(i,j)),B.addr(0,i),1,
//          X.addr(0,j),1);
//      }
//    }
      Vector<float,complex<float> > *r=
        OPERATOR_NEW Vector<float,complex<float> >(m);
      for (int l=0;l<k;l++) {
        F77NAME(ccopy)(m,B.addr(l,0),k,r->addr(),1);
        for (int j=0;j<n;j++) {
          complex<float> &Xlj=X(l,j);
          Xlj=F77NAME(cdotc)(m,addr(0,j),1,r->addr(),1);
          F77NAME(caxpy)(m,-Xlj,addr(0,j),1,r->addr(),1);
        }
      }
      delete r; r=0;
    } else { // X=B*Q^H is smallest solution
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(cgemm)('N','C',k,m,n,float_one_,B.addr(),k,addr(),m,
        complex_float_zero_,X.addr(),k);
    }
  }
}

template<> void OrthogonalMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  int m=size(0), n=size(1);
  if (trans=='C' || trans=='c') { // x=Q*b is smallest solution
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(cgemv)('N',m,n,float_one_,addr(),m,b.addr(),1,
      complex_float_zero_,x.addr(),1);
  } else if (trans=='T' || trans=='t') {
  // x=conj(Q)*b=conj(Q*conj(b)) is smallest
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    x=complex_float_zero_;
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(m,conj(b[j]),addr(0,j),1,x.addr(),1);
    }
    complex<float> *xi=x.addr();
    for (int i=0;i<m;i++,xi++) *xi=conj(*xi);
  } else { // x=Q^H*b minimizes residual
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(cgemv)('C',m,n,float_one_,addr(),m,b.addr(),1,
      complex_float_zero_,x.addr(),1);
  }
}

template class OrthogonalMatrix<float,complex<float> >;
template void testOrthogonalMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SymmetricMatrix.C"

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,addr(j,j),1,M->addr(j,j),1);
    if (j+1<n) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Mji=M->addr(j,j+1);
      for (int i=j+1;i<n;i++,Aij++,Mji+=n) *Mji=conj(*Aij);
    }
  }
  return M;
}

template<> void SymmetricMatrix<float,complex<float> >::fillWith(
complex<float> d) {
  Matrix<float,complex<float> >::set('L',d,d.real());
}

template<> complex<float>
SymmetricMatrix<float,complex<float> >::operator()(int i,int j) const {
  return (j<=i ? *(this->addr(i,j)) : conj(*(this->addr(j,i))) );
}

template<> SymmetricMatrix<float,complex<float> >&
SymmetricMatrix<float,complex<float> >::operator+=(
const SymmetricMatrix<float,complex<float> > &S) {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0))
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_one_,S.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<float,complex<float> >&
SymmetricMatrix<float,complex<float> >::operator-=(
const SymmetricMatrix<float,complex<float> > &S) {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0))
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_mone_,S.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<float,complex<float> >&
SymmetricMatrix<float,complex<float> >::operator*=(
complex<float> scalar) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    F77NAME(csscal)(n-j,scalar.real(),addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<float,complex<float> >&
SymmetricMatrix<float,complex<float> >::operator/=(
complex<float> scalar) {
  int n=this->size(0);
  CHECK_NONZERO(abs(scalar))
  for (int j=0;j<n;j++) {
    F77NAME(csscal)(n-j,float_one_/scalar.real(),addr(j,j),1);
  }
  return *this;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  int n=this->size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *col_j=addr(j+1,j);
      complex<float> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(min(j+1,n),one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j,n),one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      if (j<n) (*S)(j,j)+=one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int n=this->size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *col_j=addr(j+1,j);
      complex<float> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(n-j,one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(caxpy)(n-j-1,one_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      (*S)(j,j)+=one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator+(
const SquareMatrix<float,complex<float> > &S) const {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *T=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  T->copy(S);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,one_,addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *col_j=addr(j+1,j);
      complex<float> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,T_row_j+=n) *T_row_j+=conj(*col_j);
    }
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  int n=this->size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,one_,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *col_j=addr(j+1,j);
      complex<float> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j+=conj(*col_j);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  int n=this->size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *col_j=addr(j+1,j);
      complex<float> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(min(j+1,n),mone_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j,n),mone_,U.addr(0,j),1,S->addr(0,j),1);
      }
      if (j<n) (*S)(j,j)+=mone_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const UnitUpperTrapezoidalMatrix<float,complex<float> > &U,
const SymmetricMatrix<float,complex<float> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,complex<float> > *T=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(ccopy)(min(j,n),U.addr(0,j),1,T->addr(0,j),1);
    }
    if (j<n) (*T)(j,j)=complex_float_one_;
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *S_col_j=S.addr(j+1,j);
      complex<float> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const SymmetricMatrix<float,complex<float> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,complex<float> > *T=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+1,n),U.addr(0,j),1,T->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(ccopy)(min(j,n),U.addr(0,j),1,T->addr(0,j),1);
      }
      if (j<n) (*T)(j,j)=complex_float_one_;
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *S_col_j=S.addr(j+1,j);
      complex<float> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int n=this->size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *col_j=addr(j+1,j);
      complex<float> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(n-j,mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(caxpy)(n-j-1,mone_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      (*S)(j,j)+=mone_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const UnitLowerTrapezoidalMatrix<float,complex<float> > &L,
const SymmetricMatrix<float,complex<float> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *T=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    if (j<n-1) {
      F77NAME(ccopy)(n-j-1,L.addr(j+1,j),1,T->addr(j+1,j),1);
    }
    (*T)(j,j)=complex_float_one_;
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *S_col_j=S.addr(j+1,j);
      complex<float> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const SymmetricMatrix<float,complex<float> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *T=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(n-j,L.addr(j,j),1,T->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(ccopy)(n-j-1,L.addr(j+1,j),1,T->addr(j+1,j),1);
      }
      (*T)(j,j)=complex_float_one_;
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *S_col_j=S.addr(j+1,j);
      complex<float> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator-(
const SquareMatrix<float,complex<float> > &S) const {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *T=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *col_j=addr(j+1,j);
      complex<float> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,T_row_j+=n) *T_row_j=conj(*col_j);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n,mone_,S.addr(0,j),1,T->addr(0,j),1);
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const SquareMatrix<float,complex<float> > &S,
const SymmetricMatrix<float,complex<float> > &SS) {
  int n=SS.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *T=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  T->copy(S);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_mone_,SS.addr(j,j),1,
      T->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *S_col_j=S.addr(j+1,j);
      complex<float> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  int n=this->size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *col_j=addr(j+1,j);
      complex<float> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_row_j+=n) *S_row_j=conj(*col_j);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n,mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const SymmetricMatrix<float,complex<float> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *T=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  T->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *S_col_j=S.addr(j+1,j);
      complex<float> *T_row_j=T->addr(j,j+1);
      for (int i=j+1;i<n;i++,S_col_j++,T_row_j+=n) {
        *T_row_j-=conj(*S_col_j);
      }
    }
  }
  return T;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator*(complex<float> d)
const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    const complex<float> *col_j=addr(j,j);
    complex<float> *S_col_j=S->addr(j,j);
    *S_col_j=(*col_j)*d;
    if (j<n-1) {
      S_col_j++;
      col_j++;
      complex<float> *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_col_j++,S_row_j+=n) {
        *S_col_j=(*col_j)*d;
        *S_row_j=conj(*col_j)*d;
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator/(complex<float> d)
const {
  CHECK_NONZERO(abs(d));
  return operator*(complex_float_one_/d);
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator*(
const SymmetricMatrix<float,complex<float> > &S) const {
// compute by bordering: note that
// [ sigma s^H ] [ tau t^H ] = [ sigma tau + s^H t , sigma t^H + s^H T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^H +  S  T ]
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *T=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-1;k>=0;k--) {
    if (k<n-1) {
      F77NAME(chemv)('L',n-k-1,complex_float_one_,addr(k+1,k+1),n,
        S.addr(k+1,k),1,complex_float_zero_,T->addr(k+1,k),1); // S t
      complex<float> *tt=t+(k+1);
      F77NAME(chemv)('L',n-k-1,complex_float_one_,S.addr(k+1,k+1),n,
        addr(k+1,k),1,complex_float_zero_,tt,1); // s^H T = ( T s )^H
      complex<float> *tk=T->addr(k,k+1);
      for (int i=k+1;i<n;i++,tk+=n,tt++) *tk=conj(*tt);
      (*T)(k,k)=F77NAME(cdotc)(n-k-1,addr(k+1,k),1,S.addr(k+1,k),1);//s^Ht
    }
    F77NAME(cgerc)(n-k,n-k,complex_float_one_,addr(k,k),1,S.addr(k,k),1,
      T->addr(k,k),n);
  }
  delete [] t; t=0;
  return T;
}

template<> Matrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
// compute by bordering: note that
// S [ U_1 , U_2 ] = [ S U_1 , S U_2 ]
// and that
// [ sigma s^H ] [ upsilon u^T ] = [ sigma upsilon , sigma u^T + s^H U ]
// [   s    S  ] [          U  ] [ [   s   upsilon ,   s   u^T +  S  U ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  char diag=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0
    ? 'N' : 'U');
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,m);
  for (int k=m-1;k>=0;k--) { // S U_1
    if (k<m-1) {
      F77NAME(ccopy)(m-k-1,addr(k+1,k),1,t+(k+1),1);
      complex<float> *tt=t+(k+1);
      F77NAME(ctrmv)('U','C',diag,m-k-1,U.addr(k+1,k+1),m,tt,1);
      complex<float> *mk=M->addr(k,k+1);
      for (int j=k+1;j<m;j++,mk+=m,tt++) *mk=conj(*tt); // s^H U
    }
    if (diag=='N') {
      F77NAME(cgeru)(m-k,m-k,complex_float_one_,addr(k,k),1,
        U.addr(k,k),m,M->addr(k,k),m);
    } else {
      F77NAME(caxpy)(m-k,complex_float_one_,addr(k,k),1,M->addr(k,k),1);
      F77NAME(cgeru)(m-k,m-k-1,complex_float_one_,addr(k,k),1,
        U.addr(k,k+1),m,M->addr(k,k+1),m);
    }
  }
  if (n>m) { // S U_2
    F77NAME(chemm)('L','L',m,n-m,complex_float_one_,addr(),m,
      U.addr(0,m),m,complex_float_zero_,M->addr(0,m),m);
  }
  delete [] t; t=0;
  return M;
}

template<> Matrix<float,complex<float> >* operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const SymmetricMatrix<float,complex<float> > &S) {
//compute by bordering: note that
// [ U_1 U_2 ] [ S_11 S_21^H ]
//             [ S_21  S_22  ]
//   = [ U_1 S_11 + U_2 S_21 , U_1 S_21^H + U_2 S_22 ]
// and that
// [ U    u    ] [  S  bar(s)] = [ U S +   u    s^T , U bar(s) + u sigma ]
// [   upsilon ] [ s^T sigma ] = [      upsilon s^T ,    upsilon   sigma ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  char diag=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0
    ? 'N' : 'U');
  for (int k=0;k<m;k++) { // U_1 S_11
    if (k>0) {
      complex<float> *mk=M->addr(0,k);
      const complex<float> *sk=S.addr(k,0);
      for (int i=0;i<k;i++,mk++,sk+=n) *mk=conj(*sk);
      F77NAME(ctrmv)('U','N',diag,k,U.addr(),m,M->addr(0,k),1); //U bar(s)
    }
    if (diag=='N') {
      F77NAME(cgeru)(k+1,k+1,complex_float_one_,U.addr(0,k),1,
        S.addr(k,0),n,M->addr(),m);
    } else {
      F77NAME(cgeru)(k,k+1,complex_float_one_,U.addr(0,k),1,
        S.addr(k,0),n,M->addr(),m);
      F77NAME(caxpy)(k+1,complex_float_one_,S.addr(k,0),n,
        M->addr(k,0),m);
    }
  }
  if (n>m) {
    F77NAME(cgemm)('N','N',m,m,n-m,complex_float_one_,U.addr(0,m),m,
      S.addr(m,0),n,complex_float_one_,M->addr(),m); // U_2 S_21
    for (int j=m;j<n;j++) {
      complex<float> *mj=M->addr(0,j);
      const complex<float> *sj=S.addr(j,0);
      for (int i=0;i<m;i++,mj++,sj+=n) *mj=conj(*sj);
    }
    F77NAME(ctrmm)('L','U','N',diag,m,n-m,complex_float_one_,U.addr(),m,
      M->addr(0,m),m); // U_1 S_21^T
    F77NAME(chemm)('R','L',m,n-m,complex_float_one_,S.addr(m,m),n,
      U.addr(0,m),m,complex_float_one_,M->addr(0,m),m); // U_2 S_22
  }
  return M;
}

template<> Matrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
// compute by bordering: note that
// [ S_11 S_21^H ] [ L_1 ] = [ S_11 L_1 + S_21^H L_2 ]
// [ S_21  S_22  ] [ L_2 ] = [ S_21 L_1 +  S_22  L_2 ]
// and that
// [  S  bar(s)] [   L          ] = [  S  L +   s   ell^T ,bar(s) lambda ]
// [ s^T sigma ] [ ell^T lambda ] = [ s^T L + sigma ell^T , sigma lambda 
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(1));
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0
    ? 'N' : 'U');
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,m);
  for (int k=0;k<n;k++) { // S_11 L_1
    if (k>0) {
      F77NAME(ccopy)(k,addr(k,0),m,M->addr(k,0),m);
      F77NAME(ctrmv)('L','T',diag,k,L.addr(),m,M->addr(k,0),m); // s^T L
    }
    const complex<float> *sk=addr(k,0);
    complex<float> *tt=t;
    for (int i=0;i<=k;i++,sk+=m,tt++) *tt=conj(*sk);
    if (diag=='N') {
      F77NAME(cgeru)(k+1,k+1,complex_float_one_,t,1,L.addr(k,0),m,
        M->addr(),m);
    } else {
      F77NAME(cgeru)(k+1,k,complex_float_one_,t,1,L.addr(k,0),m,
        M->addr(),m);
      F77NAME(caxpy)(k+1,complex_float_one_,t,1,M->addr(0,k),1);
    }
  }
  if (m>n) {
    F77NAME(cgemm)('C','N',n,n,m-n,complex_float_one_,addr(n,0),m,
      L.addr(n,0),m,complex_float_one_,M->addr(),m); // S_21^H L_2
    F77NAME(clacpy)('A',m-n,n,addr(n,0),m,M->addr(n,0),m);
    F77NAME(ctrmm)('R','L','N',diag,m-n,n,complex_float_one_,L.addr(),m,
      M->addr(n,0),m); // S_21 L_1
    F77NAME(chemm)('L','L',m-n,n,complex_float_one_,addr(n,n),m,
      L.addr(n,0),m,complex_float_one_,M->addr(n,0),m); // S_22 L_2
  }
  delete [] t; t=0;
  return M;
}

template<> Matrix<float,complex<float> >* operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const SymmetricMatrix<float,complex<float> > &S) {
// compute by bordering: note that
// [ L_1 ] S = [ L_1 S ]
// [ L_2 ]   = [ L_2 S ]
// and that
// [ lambda   ] [ sigma s^H ] = [ lambda sigma       , lambda s^H       ]
// [  ell   L ] [   s    S  ] = [  ell   sigma + L s ,   ell  s^H + L S ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0
    ? 'N' : 'U');
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-1;k>=0;k--) { // L_1 S
    if (k<n-1) {
      F77NAME(ccopy)(n-k-1,S.addr(k+1,k),1,M->addr(k+1,k),1);
      F77NAME(ctrmv)('L','N',diag,n-k-1,L.addr(k+1,k+1),m,
        M->addr(k+1,k),1); // L s
    }
    complex<float> *tk=t+k;
    F77NAME(ccopy)(n-k,S.addr(k,k),1,tk,1);
    if (diag=='N') {
      F77NAME(cgerc)(n-k,n-k,complex_float_one_,L.addr(k,k),1,
        tk,1,M->addr(k,k),m);
    } else {
      F77NAME(cgerc)(n-k-1,n-k,complex_float_one_,L.addr(k+1,k),1,
        tk,1,M->addr(k+1,k),m);
      for (int i=k;i<n;i++) t[i]=conj(t[i]);
      F77NAME(caxpy)(n-k,complex_float_one_,tk,1,M->addr(k,k),m);
    }
  }
  if (m>n) { // L_2 S
    F77NAME(chemm)('R','L',m-n,n,complex_float_one_,S.addr(),n,
      L.addr(n,0),m,complex_float_zero_,M->addr(n,0),m);
  } 
  delete [] t; t=0;
  return M;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator*(
const SquareMatrix<float,complex<float> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *T=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  F77NAME(chemm)('L','L',n,n,complex_float_one_,addr(),n,S.addr(),n,
    complex_float_zero_,T->addr(),n);
  return T;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SquareMatrix<float,complex<float> > &A,
const SymmetricMatrix<float,complex<float> > &B) {
// compute by bordering: note that
// [ mu v^T ] [ sigma s^H ] = [ mu sigma + v^T s , mu s^H + v^T S ]
// [  m  M  ] [   s    S  ] = [  m sigma +  M  s ,  m s^H +  M  S ]
  int n=A.size(0);
  CHECK_SAME(n,B.size(0));
  SquareMatrix<float,complex<float> > *T=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  complex<float> *v=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-1;k>=0;k--) {
    if (k+1<n) {
      F77NAME(cgemv)('N',n-k-1,n-k-1,float_one_,A.addr(k+1,k+1),n,
        B.addr(k+1,k),1,float_zero_,T->addr(k+1,k),1); // M s
      complex<float> *ti=t;
      const complex<float> *Aik=A.addr(k,k+1);
      for (int i=k+1;i<n;i++,ti++,Aik+=n) *ti=conj(*Aik);
      F77NAME(chemv)('L',n-k-1,float_one_,B.addr(k+1,k+1),n,
        t,1,float_zero_,v,1);
      complex<float> *vi=v;
      complex<float> *Tki=T->addr(k,k+1);
      for (int i=k+1;i<n;i++,vi++,Tki+=n) *Tki=conj(*vi);
      //v^T S = (S bar(v))^H
      (*T)(k,k)=
        F77NAME(cdotu)(n-k-1,A.addr(k,k+1),n,B.addr(k+1,k),1); // v^T s
    }
    F77NAME(cgerc)(n-k,n-k,float_one_,A.addr(k,k),1,B.addr(k,k),1,
      T->addr(k,k),n);
  }
  delete [] t; t=0;
  delete [] v; v=0;
  return T;
}

template<> Matrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  int m=size(0),n=M.size(1);
  CHECK_SAME(m,M.size(0));
  Matrix<float,complex<float> > *T=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(chemm)('L','L',m,n,complex_float_one_,addr(),m,M.addr(),m,
    complex_float_zero_,T->addr(),m);
  return T;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const SymmetricMatrix<float,complex<float> > &S) {
  int m=M.size(0),n=S.size(0);
  CHECK_SAME(n,M.size(1));
  Matrix<float,complex<float> > *T=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  F77NAME(chemm)('R','L',m,n,complex_float_one_,S.addr(),n,M.addr(),m,
    complex_float_zero_,T->addr(),m);
  return T;
}

template<> Vector<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<float,complex<float> > *w=
    OPERATOR_NEW Vector<float,complex<float> >(n);
  F77NAME(chemv)('L',n,complex_float_one_,addr(),n,v.addr(),1,
    complex_float_zero_,w->addr(),1);
  return w;
}

// y := A * x * alpha + y * beta
template<> void SymmetricMatrix<float,complex<float> >::symv(
complex<float> alpha,const Vector<float,complex<float> > &x,
complex<float> beta,Vector<float,complex<float> > &y) const {
  int n=this->size(0);
  F77NAME(chemv)('L',n,alpha,addr(),n,x.addr(),1,beta,y.addr(),1);
}

// A += x * alpha * x^T
template<> void SymmetricMatrix<float,complex<float> >::syr(
complex<float> alpha,const Vector<float,complex<float> > &x) {
  int n=this->size(0);
  F77NAME(cher)('L',n,alpha,x.addr(),1,addr(),n);
}

// A += x * alpha * y^T + y * alpha * x^T
template<> void SymmetricMatrix<float,complex<float> >::syr2(
complex<float> alpha,const Vector<float,complex<float> > &x,
const Vector<float,complex<float> > &y) {
  int n=this->size(0);
  F77NAME(cher2)('L',n,alpha,x.addr(),1,y.addr(),1,addr(),n);
}

// C := A * alpha * B + C * beta
template<> void SymmetricMatrix<float,complex<float> >::symm(
complex<float> alpha,const Matrix<float,complex<float> > &B,
complex<float> beta,Matrix<float,complex<float> > &C,char side) const {
  int m=C.size(0),n=C.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,size(0));
    F77NAME(chemm)('L','L',m,n,alpha,addr(),m,B.addr(),m,beta,
      C.addr(),m);
  } else {
    CHECK_SAME(n,size(0));
    F77NAME(chemm)('R','L',m,n,alpha,addr(),n,B.addr(),m,beta,
      C.addr(),m);
  }
}

// C := A * alpha * A^T + C * beta
template<> void SymmetricMatrix<float,complex<float> >::syrk(
complex<float> alpha,const Matrix<float,complex<float> > &A,
complex<float> beta,char transa) {
  int m=A.size(0),n=A.size(1);
  if (transa!='N' && transa!='n') {
    CHECK_SAME(n,size(0));
    F77NAME(cherk)('L',transa,n,m,alpha,A.addr(),m,beta,addr(),n);
  } else {
    CHECK_SAME(m,size(0));
    F77NAME(cherk)('L',transa,m,n,alpha,A.addr(),m,beta,addr(),m);
  }
}

// C := A * alpha * B^T + B * alpha + A^T + C * beta
template<> void SymmetricMatrix<float,complex<float> >::syr2k(
complex<float> alpha,const Matrix<float,complex<float> > &A,
const Matrix<float,complex<float> > &B,complex<float> beta,
char transab) {
  int m=A.size(0),n=A.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (transab!='N' && transab!='n') {
    CHECK_SAME(n,size(0));
    F77NAME(cher2k)('L',transab,n,m,alpha,A.addr(),m,B.addr(),m,beta,
      addr(),n);
  } else {
    CHECK_SAME(m,size(0));
    F77NAME(cher2k)('L',transab,m,n,alpha,A.addr(),m,B.addr(),m,beta,
      addr(),m);
  }
}

/*
// y := abs(A) * abs(x) * alpha + abs(y) * beta
template<> void SymmetricMatrix<float,complex<float> >::syamv(float alpha,
const Vector<float,complex<float> > &x,float beta,Vector<float,complex<float> > &y)
const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,y.size());
  F77_NAME(sla_syamv)('L',n,alpha,addr(),n,x.addr(),1,beta,y.addr(),1);
}
*/

template<> float SymmetricMatrix<float,complex<float> >::equilibrate(
Vector<float,float> &s,float &scond) const {
  int n=size(0);
  CHECK_SAME(n,s.size());
  float amax=numeric_limits<float>::infinity();
  int info;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,3*n);
  F77NAME(cheequb)('L',n,addr(),n,s.addr(),scond,amax,work,info);
  CHECK_SAME(info,0);
  delete [] work;
  return amax;
}

template<> float
SymmetricMatrix<float,complex<float> >::normFrobenius() const {
  int n=size(0);
  float *work=0;
  return F77NAME(clanhe)('F','L',n,addr(),n,work);
}

template<> float
SymmetricMatrix<float,complex<float> >::normInfinity() const {
  int n=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,n);
  float val=F77NAME(clanhe)('I','L',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> float
SymmetricMatrix<float,complex<float> >::normMaxEntry() const {
  int n=size(0);
  float *work=0;
  return F77NAME(clanhe)('M','L',n,addr(),n,work);
}

template<> float
SymmetricMatrix<float,complex<float> >::normOne() const {
  int n=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,n);
  float val=F77NAME(clanhe)('O','L',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> float
SymmetricMatrix<float,complex<float> >::reciprocalConditionNumber()
const {
  int n=size(0);

  SymmetricMatrix<float,complex<float> > *AF=
    OPERATOR_NEW SymmetricMatrix<float,complex<float> >(*this);
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  complex<float> w=complex_float_undefined_;
  F77NAME(chetrf)('L',n,AF->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0);

  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(chetrf)('L',n,AF->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;

  float *rwork=OPERATOR_NEW_BRACKET(float,2*n);
  float anorm=F77NAME(clanhe)('O','L',n,addr(),n,rwork);
  delete rwork; rwork=0;

  float rcond;
  work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  F77NAME(checon)('L',n,AF->addr(),n,ipiv,anorm,rcond,work,info);
  delete [] ipiv; ipiv=0;
  delete [] work; work=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SymmetricMatrix<float,complex<float> >*
SymmetricMatrix<float,complex<float> >::inverse() const {
  int n=size(0);
  SymmetricMatrix<float,complex<float> > *Ainv=
    OPERATOR_NEW SymmetricMatrix<float,complex<float> >(*this);
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  complex<float> w=undefined_;
  F77NAME(chetrf)('L',n,Ainv->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0);

  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(chetrf)('L',n,Ainv->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;

  work=OPERATOR_NEW_BRACKET(complex<float>,n);
  F77NAME(chetri)('L',n,Ainv->addr(),n,ipiv,work,info);
  CHECK_SAME(info,0)

  delete [] ipiv;
  delete [] work;
  return Ainv;
}
*/

template<> Vector<float,float>*
SymmetricMatrix<float,complex<float> >::eigenvalues(
OrthogonalMatrix<float,complex<float> > *&Q) const {
  int n=size(0);
  if (Q!=0) CHECK_SAME(n,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Q->Matrix<float,complex<float> >::copyFrom('L',n,n,*this);
  Vector<float,float> *lambda=OPERATOR_NEW Vector<float,float>(n);
  complex<float> w;
  int lwork=-1,info;
  float *rwork=OPERATOR_NEW_BRACKET(float,max(1,3*n-2));
  F77NAME(cheev)(jobz,'L',n,Q->addr(),n,lambda->addr(),&w,lwork,rwork,
    info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(cheev)(jobz,'L',n,Q->addr(),n,lambda->addr(),work,lwork,rwork,
    info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete [] rwork; rwork=0;

  return lambda;
}

template<> void SymmetricMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char) const {
  int n=size(0);
  SymmetricMatrix<float,complex<float> > AF(*this);
  int info;
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  x.copy(b);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  complex<float> w=complex_float_undefined_;
  F77NAME(chetrf)('L',n,AF.addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(chetrf)('L',n,AF.addr(),n,ipiv,work,lwork,info);
  delete [] work; work=0;

  if (info==0) {
    F77NAME(chetrs)('L',n,1,AF.addr(),n,ipiv,x.addr(),n,info);
  }
  CHECK_SAME(info,0)
  delete [] ipiv; ipiv=0;
}

template<> void SymmetricMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char) const {
  int n=size(0);
  SymmetricMatrix<float,complex<float> > AF(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  complex<float> w=complex_float_undefined_;
  int info;
  F77NAME(chetrf)('L',n,AF.addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(chetrf)('L',n,AF.addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;

  bool left_side=(side=='L' || side=='l');
  if (left_side) {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1))
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    X.copy(B);
    F77NAME(chetrs)('L',n,1,AF.addr(),n,ipiv,X.addr(),n,info);
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0))
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
    for (int i=0;i<k;i++) {
      complex<float> *tt=t;
      const complex<float> *bi=B.addr(i,0);
      for (int j=0;j<n;j++,tt++,bi+=k) *tt=conj(*bi);
      F77NAME(chetrs)('L',n,1,AF.addr(),n,ipiv,t,n,info);
      CHECK_SAME(info,0)
      tt=t;
      complex<float> *xi=X.addr(i,0);
      for (int j=0;j<n;j++,tt++,xi+=k) *xi=conj(*tt);
    }
    delete [] t; t=0;
  }
  delete [] ipiv; ipiv=0;
}

template class SymmetricMatrix<float,complex<float> >;
template void testSymmetricMatrix(float,complex<float> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveMatrix<float,complex<float> >&
SymmetricPositiveMatrix<float,complex<float> >::operator*=(
complex<float> scalar) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    F77NAME(csscal)(n-j,abs(scalar),addr(j,j),1);
  }
  return *this;
}

template<> SymmetricPositiveMatrix<float,complex<float> >&
SymmetricPositiveMatrix<float,complex<float> >::operator/=(
complex<float> scalar) {
  int n=this->size(0);
  CHECK_NONZERO(abs(scalar))
  for (int j=0;j<n;j++) {
    F77NAME(csscal)(n-j,float_one_/abs(scalar),addr(j,j),1);
  }
  return *this;
}

template<> float
SymmetricPositiveMatrix<float,complex<float> >::equilibrate(
Vector<float,float> &s,float &scond) const {
  int n=size(0);
  CHECK_SAME(n,s.size());
  float amax=numeric_limits<float>::infinity();
  int info;
  F77NAME(cpoequb)(n,addr(),n,s.addr(),scond,amax,info);
  CHECK_SAME(info,0);
  return amax;
}

template<> float
SymmetricPositiveMatrix<float,complex<float> >::reciprocalConditionNumber(
) const {
  int n=size(0);
  SymmetricPositiveMatrix<float,complex<float> > *AF=
    OPERATOR_NEW SymmetricPositiveMatrix<float,complex<float> >(*this);
  int info;
  F77NAME(cpotrf)('L',n,AF->addr(),n,info);
  CHECK_SAME(info,0)

  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  float anorm=F77NAME(clanhe)('O','L',n,addr(),n,rwork);
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float rcond;
  F77NAME(cpocon)('L',n,addr(),n,anorm,rcond,work,rwork,info);
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SymmetricPositiveMatrix<float,complex<float> >*
SymmetricPositiveMatrix<float,complex<float> >::inverse() const {
  int n=size(0);
  SymmetricPositiveMatrix<float,complex<float> > *Ainv=
    OPERATOR_NEW SymmetricPositiveMatrix<float,complex<float> >(*this);
  int info;
  F77NAME(cpotrf)('L',n,Ainv->addr(),n,info);
  CHECK_SAME(info,0);

  F77NAME(cpotri)('L',n,Ainv->addr(),n,info);
  CHECK_SAME(info,0)
  return Ainv;
}
*/

template<> void SymmetricPositiveMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char) const {
  int n=size(0);
  SymmetricPositiveMatrix<float,complex<float> > AF(*this);
  int info;
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  x.copy(b);
  F77NAME(cpotrf)('L',n,AF.addr(),n,info);
  CHECK_SAME(info,0)

  F77NAME(cpotrs)('L',n,1,AF.addr(),n,x.addr(),n,info);
  CHECK_SAME(info,0)
}

template<> void SymmetricPositiveMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char) const {
  int n=size(0);
  SymmetricPositiveMatrix<float,complex<float> > AF(*this);
  int info;
  F77NAME(cpotrf)('L',n,AF.addr(),n,info);
  CHECK_SAME(info,0)

  bool left_side=(side=='L' || side=='l');
  if (left_side) {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1))
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    X.copy(B);
    F77NAME(cpotrs)('L',n,1,AF.addr(),n,X.addr(),n,info);
    CHECK_SAME(info,0)
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0))
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
    for (int i=0;i<k;i++) {
      complex<float> *tt=t;
      const complex<float> *bi=B.addr(i,0);
      for (int j=0;j<n;j++,tt++,bi+=k) *tt=conj(*bi);
      F77NAME(cpotrs)('L',n,1,AF.addr(),n,t,n,info);
      CHECK_SAME(info,0)
      tt=t;
      complex<float> *xi=X.addr(i,0);
      for (int j=0;j<n;j++,tt++,xi+=k) *xi=conj(*tt);
    }
    delete [] t; t=0;
  }
}

template class SymmetricPositiveMatrix<float,complex<float> >;
template void testSymmetricPositiveMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "BandMatrix.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> complex<float>
  TridiagonalMatrix<float,complex<float> >::safety_ =
  complex_float_zero_;
template<> const complex<float>
  TridiagonalMatrix<float,complex<float> >::outofbounds_ =
  complex_float_zero_;
template<> const complex<float>
  TridiagonalMatrix<float,complex<float> >::undefined_ =
  complex<float>(numeric_limits<float>::infinity(),
    numeric_limits<float>::infinity());

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  F77NAME(ccopy)(n,D->addr(),1,M->addr(0,0),n+1);
  F77NAME(ccopy)(n-1,L->addr(),1,M->addr(1,0),n+1);
  F77NAME(ccopy)(n-1,U->addr(),1,M->addr(0,1),n+1);
  return M;
}

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator+(
const SymmetricMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int k=0;k<n;k++) {
    F77NAME(ccopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) F77NAME(ccopy)(n-k-1,M.addr(k+1,k),1,S->addr(k,k+1),n);
  }
  F77NAME(caxpy)(n,complex_float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(ccopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)=complex_float_one_;
      if (k<n-1) F77NAME(ccopy)(n-k-1,M.addr(k+1,k),1,S->addr(k+1,k),1);
    }
  }
  F77NAME(caxpy)(n,complex_float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<float,complex<float> > *S=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(ccopy)(k+1,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(ccopy)(k,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)=complex_float_one_;
    }
  }
  F77NAME(caxpy)(n,complex_float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  S->copy(M);
  F77NAME(caxpy)(n,complex_float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator-(
const SymmetricMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  F77NAME(ccopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(ccopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    F77NAME(caxpy)(n-k,float_mone_,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) {
      F77NAME(caxpy)(n-k-1,float_mone_,M.addr(k+1,k),1,S->addr(k,k+1),n);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const SymmetricMatrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) { 
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  for (int k=0;k<n;k++) {
    F77NAME(ccopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) F77NAME(ccopy)(n-k-1,M.addr(k+1,k),1,S->addr(k,k+1),n);
  }
  F77NAME(caxpy)(n,complex_float_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  F77NAME(ccopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(ccopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(caxpy)(n-k,float_mone_,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)-=float_one_;
      if (k<n-1) {
        F77NAME(caxpy)(n-k-1,float_mone_,M.addr(k+1,k),1,
          S->addr(k+1,k),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) { 
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(ccopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)=float_one_;
      if (k<n-1) {
        F77NAME(ccopy)(n-k-1,M.addr(k+1,k),1,S->addr(k+1,k),1);
      }
    }
  }
  F77NAME(caxpy)(n,complex_float_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<float,complex<float> > *S=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  F77NAME(ccopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(ccopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(caxpy)(k+1,float_mone_,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(caxpy)(k,float_mone_,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)-=float_one_;
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >* operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<float,complex<float> > *S=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(ccopy)(k+1,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(ccopy)(k,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)=float_one_;
    }
  }
  F77NAME(caxpy)(n,complex_float_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  F77NAME(ccopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(ccopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  *S-=M;
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  S->copy(M);
  F77NAME(caxpy)(n,complex_float_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator*(
const SymmetricMatrix<float,complex<float> > &M) const {
// compute by bordering: note that
// [    tau     upsilon e_0^T ] [ sigma s^H ]
// [ e_0 lambda       T       ] [   s    S  ]
//   = [ tau sigma + upsilon e_0^T s , tau s^H + upsilon e_0^T S ]
//   = [ e_0 lambda sigma +      T s , e_0 lambda s^H + T S      ]
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  (*S)(n-2,n-2)=(*D)[n-2]*M(n-2,n-2)+(*U)[n-2]*M(n-1,n-2);
  (*S)(n-1,n-2)=(*L)[n-2]*M(n-2,n-2)+(*D)[n-1]*M(n-1,n-2);
  (*S)(n-2,n-1)=(*D)[n-2]*conj(M(n-1,n-2))+(*U)[n-2]*M(n-1,n-1);
  (*S)(n-1,n-1)=(*L)[n-2]*conj(M(n-1,n-2))+(*D)[n-1]*M(n-1,n-1);
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-3;k>=0;k--) {
    F77NAME(cgtmv)(n-k-1,complex_float_one_,L->addr(k+1),D->addr(k+1),
      U->addr(k+1),M.addr(k+1,k),1,complex_float_zero_,S->addr(k+1,k),1);
      // T s
    complex<float> *ti=t;
    const complex<float> *Mik=M.addr(k+1,k+1);
    for (int i=k+1;i<n;i++,ti++,Mik++) *ti=conj(*Mik);
    F77NAME(caxpy)(n-k-1,(*U)[k],t,1,S->addr(k,k+1),n); // upsilon e_0^T S
    (*S)(k,k)=(*U)[k]*M(k+1,k); // upsilon e_0^T s
    ti=t;
    Mik=M.addr(k,k);
    for (int i=k;i<n;i++,ti++,Mik++) *ti=conj(*Mik);
    F77NAME(caxpy)(n-k,(*D)[k],t,1,S->addr(k,k),n); // tau [ sigma , s^H ]
    F77NAME(caxpy)(n-k,(*L)[k],t,1,S->addr(k+1,k),n);
      // e_0 lambda [ sigma, s^H ]
  }
  delete [] t; t=0;
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SymmetricMatrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) {
// compute by bordering: note that
// [ sigma s^H ] [    tau     upsilon e_0^T ]
// [   s    S  ] [ e_0 lambda       T       ]
//   = [ sigma tau + s^H e_0 lambda , sigma upsilon e_0^T + s^H T ]
//   = [     s tau +   S e_0 lambda ,     s upsilon e_0^T +   S T ]
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  (*S)(n-2,n-2)=M(n-2,n-2)*T(n-2,n-2)+conj(M(n-1,n-2))*T(n-1,n-2);
  (*S)(n-1,n-2)=M(n-1,n-2)*T(n-2,n-2)+     M(n-1,n-1) *T(n-1,n-2);
  (*S)(n-2,n-1)=M(n-2,n-2)*T(n-2,n-1)+conj(M(n-1,n-2))*T(n-1,n-1);
  (*S)(n-1,n-1)=M(n-1,n-2)*T(n-2,n-1)+     M(n-1,n-1) *T(n-1,n-1);
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-3;k>=0;k--) {
    complex<float> *tj=t;
    const complex<float> *Mjk=M.addr(k+1,k);
    for (int j=k+1;j<n;j++,tj++,Mjk++) *tj=conj(*Mjk);
    F77NAME(cgtmv)(n-k-1,float_one_,T.addr(k+1,k+2),
      T.addr(k+1,k+1),T.addr(k+2,k+1),t,1,float_zero_,S->addr(k,k+1),n);
      // s^H T = ( T^T bar(s) )^T
    F77NAME(caxpy)(n-k-1,T(k+1,k),M.addr(k+1,k+1),1,S->addr(k+1,k),1);
      // S e_0 lambda
    (*S)(k,k)=conj(M(k+1,k))*T(k+1,k); // s^H e_0 lambda
    F77NAME(caxpy)(n-k,T(k,k),M.addr(k,k),1,S->addr(k,k),1);
      // [ sigma ] tau
      // [   s   ]
    F77NAME(caxpy)(n-k,T(k,k+1),M.addr(k,k),1,S->addr(k,k+1),1);
      // [ sigma ] upsilon
      // [   s   ]
  }
  delete [] t; t=0;
  return S;
}

template<> Matrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &M) const {
// to compute jth column of product:
//   [      T_11        e_j tau_12 e_0^T ] [  0  ]
//   [ e_0 tau_21 e_j^t      T_22        ] [ ell ]
//     = [   e_j tau_12 e_0^T ell ]
//     = [        T_22 ell        ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *S=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  if (M_non_unit) {
    F77NAME(cgtmv)(m,complex_float_one_,L->addr(),D->addr(),U->addr(),
      M.addr(),1,complex_float_zero_,S->addr(),1);
    for (int j=1;j<n;j++) {
      (*S)(j-1,j)=(*U)[j-1]*M(j,j); // e_j upsilon e_0^T ell
      if (j<m-1) { // T_22 ell
        F77NAME(cgtmv)(m-j,complex_float_one_,L->addr(j),D->addr(j),
          U->addr(j),M.addr(j,j),1,complex_float_zero_,S->addr(j,j),1);
      } else (*S)(j,j)=(*D)[j]*M(j,j);
    }
  } else {
// note that
// [     tau    , upsilon e_0^T ] [ 1 ] = [ tau + upsilon e_0^T ell ]
// [ e_0 lambda ,        T      ] [ell] = [ e_0 lambda      + T ell ]
    (*S)(0,0)=(*D)[0]+(*U)[0]*M(1,0);
    (*S)(1,0)=(*L)[0];
    F77NAME(cgtmv)(m-1,complex_float_one_,L->addr(1),D->addr(1),
      U->addr(1),M.addr(1,0),1,complex_float_one_,S->addr(1,0),1);
    for (int j=1;j<n;j++) {
      (*S)(j-1,j)=(*U)[j-1]; // e_j upsilon e_0^T ell
      (*S)(j,j)=(*D)[j]; // tau
      if (j<m-1) {
        (*S)(j,j)+=(*U)[j]*M(j+1,j); // upsilon e_0^T ell
        (*S)(j+1,j)=(*L)[j]; // e_0 lambda
      }
      if (j<m-2) { // T ell
        F77NAME(cgtmv)(m-j-1,complex_float_one_,L->addr(j+1),
          D->addr(j+1),U->addr(j+1),M.addr(j+1,j),1,
          complex_float_one_,S->addr(j+1,j),1);
      } else if (j<m-1) (*S)(j+1,j)+=(*D)[j+1]*M(j+1,j);
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >* operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) {
// compute by columns: note that
// L T e_j = L e_{j-1} T_{j-1,j} + L e_j T_{j,j} + L e_{j+1{ T){j+1,j}
//         = L e_{j-1} U[j-1] + L e_j D[j] + L e_{j+1} L[j]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(m-j+1,T(j-1,j),M.addr(j-1,j-1),1,
          S->addr(j-1,j),1);
      }
      F77NAME(caxpy)(m-j,T(j,j),M.addr(j,j),1,S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(caxpy)(m-j-1,T(j+1,j),M.addr(j+1,j+1),1,S->addr(j+1,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        (*S)(j-1,j)=T(j-1,j);
        F77NAME(caxpy)(m-j,T(j-1,j),M.addr(j,j-1),1,S->addr(j,j),1);
      }
      (*S)(j,j)+=T(j,j);
      if (j<m-1) {
        F77NAME(caxpy)(m-j-1,T(j,j),M.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      if (j<n-1) {
        (*S)(j+1,j)+=T(j+1,j);
        if (j<m-2) {
          F77NAME(caxpy)(m-j-2,T(j+1,j),M.addr(j+2,j+1),1,
            S->addr(j+2,j),1);
        }
      }
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &M) const {
// compute by columns: note that
// T [ U_1 , U_2 ] = [ T U_1 , T U_2 ]
// and that
// [      T_11       ,e_j upsilon e_0^T ][ u ] = [             T_11 u ]
// [ e_0 lambda e_j^T,      T_22        ][ 0 ] = [ e_0 lambda e_j^T u ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  if (M_non_unit) {
    (*S)(0,0)=(*D)[0]*M(0,0);
    (*S)(1,0)=(*L)[0]*M(0,0);
    for (int j=1;j<m-1;j++) { // T U_1
      F77NAME(cgtmv)(j+1,complex_float_one_,L->addr(),D->addr(),
        U->addr(),M.addr(0,j),1,complex_float_zero_,S->addr(0,j),1);
      (*S)(j+1,j)=(*L)[j]*M(j,j);
    }
    for (int j=m-1;j<n;j++) { // T U_2
      F77NAME(cgtmv)(m,complex_float_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,complex_float_zero_,S->addr(0,j),1);
    }
  } else {
    (*S)(0,0)=(*D)[0];
    (*S)(1,0)=(*L)[0];
    for (int j=1;j<m;j++) {
      F77NAME(cgtmv)(j,complex_float_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,complex_float_zero_,S->addr(0,j),1);
      (*S)(j-1,j)+=(*U)[j-1];
      (*S)(j,j)=(*L)[j-1]*M(j-1,j)+(*D)[j];
      if (j<m-1) (*S)(j+1,j)=(*L)[j];
    }
    for (int j=m;j<n;j++) {
      F77NAME(cgtmv)(m,complex_float_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,complex_float_zero_,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >* operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) {
// compute by columns: note that
// U T e_j = U e_{j-1} T_{j-1,j} + U e_j T_{j,j} + U e_{j+1{ T){j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j,m),T(j-1,j),M.addr(0,j-1),1,S->addr(0,j),1);
      }
      F77NAME(caxpy)(min(j+1,m),T(j,j),M.addr(0,j),1,S->addr(0,j),1);
      if (j<n-1) {
        F77NAME(caxpy)(min(j+2,m),T(j+1,j),M.addr(0,j+1),1,
          S->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j-1,m),T(j-1,j),M.addr(0,j-1),1,
          S->addr(0,j),1);
        if (j<=m) (*S)(j-1,j)=T(j-1,j);
      }
      F77NAME(caxpy)(min(j,m),T(j,j),M.addr(0,j),1,S->addr(0,j),1);
      if (j<m) (*S)(j,j)+=T(j,j);
      if (j<n-1) {
        F77NAME(caxpy)(min(j+1,m),T(j+1,j),M.addr(0,j+1),1,
          S->addr(0,j),1);
        if (j<m-1) (*S)(j+1,j)+=T(j+1,j);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator*(
const SquareMatrix<float,complex<float> > &M) const {
  int m=M.size(0);
  CHECK_SAME(m,size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(m,complex_float_zero_);
  for (int j=0;j<m;j++) {
    F77NAME(cgtmv)(m,float_one_,L->addr(),D->addr(),U->addr(),
      M.addr(0,j),1,float_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SquareMatrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) {
  int m=M.size(0);
  CHECK_SAME(m,T.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(m,complex_float_zero_);
  for (int j=0;j<m;j++) {
    for (int k=max(0,j-1);k<=min(m-1,j+1);k++) {
      F77NAME(caxpy)(m,T(k,j),M.addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(cgtmv)(m,float_one_,L->addr(),D->addr(),U->addr(),
      M.addr(0,j),1,float_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const TridiagonalMatrix<float,complex<float> > &T) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
      F77NAME(caxpy)(m,T(k,j),M.addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Vector<float,complex<float> >*
TridiagonalMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<float,complex<float> > *w=OPERATOR_NEW
    Vector<float,complex<float> >(n,complex_float_zero_);
  F77NAME(cgtmv)(n,complex_float_one_,L->addr(),D->addr(),U->addr(),
    v.addr(),1,complex_float_zero_,w->addr(),1);
  return w;
}

template<> float
TridiagonalMatrix<float,complex<float> >::normFrobenius() const {
  return F77NAME(clangt)('F',dim,L->addr(),D->addr(),U->addr());
}

template<>
float TridiagonalMatrix<float,complex<float> >::normInfinity() const {
  return F77NAME(clangt)('I',dim,L->addr(),D->addr(),U->addr());
}

template<>
float TridiagonalMatrix<float,complex<float> >::normMaxEntry() const {
  return F77NAME(clangt)('M',dim,L->addr(),D->addr(),U->addr());
}

template<>
float TridiagonalMatrix<float,complex<float> >::normOne() const {
  return F77NAME(clangt)('O',dim,L->addr(),D->addr(),U->addr());
}

template<> float
TridiagonalMatrix<float,complex<float> >::reciprocalConditionNumber(
char norm) const {
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *DF=
    OPERATOR_NEW Vector<float,complex<float> >(dim);
  Vector<float,complex<float> > *UF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *UF2=
    OPERATOR_NEW Vector<float,complex<float> >(max(0,dim-2));
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(cgttrf)(dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),ipiv,
    info);
  float rcond=numeric_limits<float>::infinity();
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*dim);
  float anorm=F77NAME(clangt)(norm,dim,L->addr(),D->addr(),U->addr());
  F77NAME(cgtcon)(norm,dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
    ipiv,anorm,rcond,work,info);
  delete [] work; work=0;
  delete [] ipiv; ipiv=0;
  delete UF2; UF2=0;
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
  return rcond;
}

template<> void TridiagonalMatrix<float,complex<float> >::gtmv(
complex<float> alpha,const Vector<float,complex<float> > &x,
complex<float> beta,Vector<float,complex<float> > &b,char trans)
const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  if (trans=='N' || trans=='n') {
    F77NAME(cgtmv)(n,alpha,L->addr(),D->addr(),U->addr(),
      x.addr(),1,beta,b.addr(),1);
  } else if (trans=='T' || trans=='t') {
    F77NAME(cgtmv)(n,alpha,U->addr(),D->addr(),L->addr(),
      x.addr(),1,beta,b.addr(),1);
  } else {
    if (abs(beta)==float_zero_) b=complex_float_zero_;
    else b*=beta;
    if (abs(alpha)==float_zero_) return;
    if (abs(x[0])>float_zero_) {
      complex<float> temp=x[0]*alpha;
      b[0]+=conj((*D)[0])*temp;
      b[1]+=conj((*U)[0])*temp;
    }
    for (int i=1;i<n-1;i++) {
      if (abs(x[i])>float_zero_) {
        complex<float> temp=x[i]*alpha;
        b[i-1]+=conj((*L)[i-1])*temp;
        b[i]+=conj((*D)[i])*temp;
        b[i+1]+=conj((*U)[i])*temp;
      }
    }
    if (abs(x[n-1])>float_zero_) {
      complex<float> temp=x[n-1]*alpha;
      b[n-2]+=conj((*L)[n-2])*temp;
      b[n-1]+=conj((*D)[n-1])*temp;
    }
  }
}

template<> void TridiagonalMatrix<float,complex<float> >::gtmm(
complex<float> alpha,const Matrix<float,complex<float> > &X,
complex<float> beta,Matrix<float,complex<float> > &B,char side,
char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,size(0));
    if (trans=='N' || trans=='n') { // A X alpha + B beta
      for (int j=0;j<n;j++) {
        F77NAME(cgtmv)(m,alpha,L->addr(),D->addr(),U->addr(),
          X.addr(0,j),1,beta,B.addr(0,j),1);
      }
    } else if (trans=='T' || trans=='t') { // A^T X alpha + B beta
      for (int j=0;j<n;j++) {
        F77NAME(cgtmv)(m,alpha,U->addr(),D->addr(),L->addr(),
          X.addr(0,j),1,beta,B.addr(0,j),1);
      }
    } else { // A^H X alpha + B beta
      if (abs(beta)==float_zero_) B=complex_float_zero_;
      else B*=beta;
      if (abs(alpha)==float_zero_) return;
      for (int j=0;j<n;j++) {
        if (abs(X(0,j))>float_zero_) {
          complex<float> temp=X(0,j)*alpha;
          B(0,j)+=conj((*D)[0])*temp;
          B(1,j)+=conj((*U)[0])*temp;
        }
        for (int i=1;i<m-1;i++) {
          if (abs(X(i,j))>float_zero_) {
            complex<float> temp=X(i,j)*alpha;
            B(i-1,j)+=conj((*L)[i-1])*temp;
            B(i,j)+=conj((*D)[i])*temp;
            B(i+1,j)+=conj((*U)[i])*temp;
          }
        }
        if (abs(X(m-1,j))>float_zero_) {
          complex<float> temp=X(m-1,j)*alpha;
          B(m-2,j)+=conj((*L)[m-2])*temp;
          B(m-1,j)+=conj((*D)[m-1])*temp;
        }
      }
    }
  } else {
    CHECK_SAME(n,size(0));
    if (abs(beta)==float_zero_) B=float_zero_;
    else if (beta!=float_one_) B*=beta;
    if (trans=='N' || trans=='n') { // X A alpha + B beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(caxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else if (trans=='T' || trans=='t') { // X A^T alpha + B beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(caxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else { // X A^H alpha + B beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(caxpy)(m,conj((*this)(j,k))*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

template<> void TridiagonalMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *DF=
    OPERATOR_NEW Vector<float,complex<float> >(dim);
  Vector<float,complex<float> > *UF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  x.copy(b);
  int info;
  if (trans!='N' && trans!='n') {
    F77NAME(cgtsv)(dim,1,UF->addr(),DF->addr(),LF->addr(),x.addr(),dim,
      info);
  } else {
    F77NAME(cgtsv)(dim,1,LF->addr(),DF->addr(),UF->addr(),x.addr(),dim,
      info);
  }
  CHECK_SAME(info,0);
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
}

template<> void TridiagonalMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *DF=
    OPERATOR_NEW Vector<float,complex<float> >(dim);
  Vector<float,complex<float> > *UF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  int info;
  X.copy(B);
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    if (trans!='N' && trans!='n') {
      F77NAME(cgtsv)(m,n,UF->addr(),DF->addr(),LF->addr(),X.addr(),m,
        info);
    } else {
      F77NAME(cgtsv)(m,n,LF->addr(),DF->addr(),UF->addr(),X.addr(),m,
        info);
    }
  } else {
    CHECK_SAME(dim,n);
    Vector<float,complex<float> > *UF2=
      OPERATOR_NEW Vector<float,complex<float> >(dim-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
    F77NAME(cgttrf)(n,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
      ipiv,info);
    CHECK_TEST(info==0);
    complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
    char rtrans=(trans!='N' && trans!='n' ? 'N' : 'T');
    for (int i=0;i<m;i++) {
      F77NAME(ccopy)(n,B.addr(i,0),m,t,1);
      F77NAME(cgttrs)(rtrans,n,1,LF->addr(),DF->addr(),UF->addr(),
        UF2->addr(),ipiv,t,n,info);
      F77NAME(ccopy)(n,t,1,X.addr(i,0),m);
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

template class TridiagonalMatrix<float,complex<float> >;
template void testTridiagonalMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> complex<float>
  SymmetricTridiagonalMatrix<float,complex<float> >::safety_=
  complex_float_zero_;
template<> const complex<float>
  SymmetricTridiagonalMatrix<float,complex<float> >::outofbounds_=
  complex_float_zero_;
template<> const complex<float>
  SymmetricTridiagonalMatrix<float,complex<float> >::undefined_=
  complex<float>(numeric_limits<float>::infinity(),
    numeric_limits<float>::infinity());

template<> SymmetricTridiagonalMatrix<float,complex<float> >::
SymmetricTridiagonalMatrix(int n,complex<float> d) : dim(n) {
  L=OPERATOR_NEW Vector<float,complex<float> >(n-1,d);
  D=OPERATOR_NEW Vector<float,float>(n,d.real());
}

template<> SquareMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  F77NAME(ccopy)(n-1,L->addr(),1,M->addr(1,0),n+1);
  const float *Aii=D->addr();
  complex<float> *Mii=M->addr(0,0);
  for (int i=0;i<dim;i++,Aii++,Mii+=n+1) *Mii=*Aii;
  const complex<float> *Aip1i=L->addr();
  complex<float> *Miip1=M->addr(0,1);
  for (int i=0;i<dim-1;i++,Aip1i++,Miip1+=n+1) *Miip1=conj(*Aip1i);
  return M;
}

template<> void
SymmetricTridiagonalMatrix<float,complex<float> >::fillWith(
complex<float> scalar) {
  *L=scalar; *D=scalar.real();
}

template<> complex<float>
SymmetricTridiagonalMatrix<float,complex<float> >::
upperDiagonalValue(int i) const {
  return conj((*L)[i]);
}

template<> SymmetricTridiagonalMatrix<float,complex<float> >&
SymmetricTridiagonalMatrix<float,complex<float> >::operator=(
complex<float> scalar) {
  *L=scalar; *D=scalar.real(); return *this;
}

template<> SymmetricTridiagonalMatrix<float,complex<float> >&
SymmetricTridiagonalMatrix<float,complex<float> >::operator*=(
complex<float> scalar) {
  int n=this->size(0);
  (*D)*=scalar.real();
  (*L)*=scalar;
  return *this;
}

template<> SymmetricTridiagonalMatrix<float,complex<float> >&
SymmetricTridiagonalMatrix<float,complex<float> >::operator/=(
complex<float> scalar) {
  int n=this->size(0);
  CHECK_NONZERO(scalar.real())
  (*D)/=scalar.real();
  (*L)/=scalar;
  return *this;
}

template<> TridiagonalMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator*(
complex<float> d) const {
  int n=size(0);
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  float *di=D->addr();
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) *Sii=*di;
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),1);
  complex<float> *li=L->addr();
  complex<float> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1=conj(*li);
  *S*=d;
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator/(
complex<float> d) const {
  int n=size(0);
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  float *di=D->addr();
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) *Sii=*di;
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),1);
  complex<float> *li=L->addr();
  complex<float> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1=conj(*li);
  *S/=d;
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator+(
const TridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  S->copy(T);
  float *di=D->addr();
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) *Sii+=*di;
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,S->addr(1,0),1);
  complex<float> *li=L->addr();
  complex<float> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1+=conj(*li);
  return S;
}

template<> SymmetricMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator+(
const SymmetricMatrix<float,complex<float> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricMatrix<float,complex<float> > *T=
    OPERATOR_NEW SymmetricMatrix<float,complex<float> >(n);
  T->copy(S);
  float *di=D->addr();
  complex<float> *Tii=T->addr(0,0);
  for (int i=0;i<n;i++,di++,Tii+=n+1) {
    *Tii+=complex<float>(*di,float_zero_);
  }
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,T->addr(1,0),n+1);
  return T;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0); 
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    complex<float> *Sjj=S->addr(0,0);
    for (int j=0;j<n;j++,Sjj+=n+1) *Sjj=float_one_;
    for (int j=0;j<n-1;j++) {
      F77NAME(ccopy)(n-j-1,M.addr(j+1,j),1,S->addr(j+1,j),1);
    }
  }
  float *di=D->addr();
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii+=n+1) {
    *Sii+=complex<float>(*di,float_zero_);
  }
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,S->addr(1,0),n+1);
  complex<float> *lj=L->addr();
  complex<float> *Sjjp1=S->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Sjjp1+=n+1) *Sjjp1+=conj(*lj);
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
   UpperHessenbergMatrix<float,complex<float> >(n,complex_float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0); 
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(j+1,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=1;j<n;j++) F77NAME(ccopy)(j,U.addr(0,j),1,H->addr(0,j),1);
    complex<float> *Hjj=H->addr(0,0);
    for (int j=0;j<n;j++,Hjj+=n+1) *Hjj=complex_float_one_;
  }
  float *di=D->addr();
  complex<float> *Hii=H->addr(0,0);
  for (int i=0;i<n;i++,di++,Hii+=n+1) {
    *Hii+=complex<float>(*di,float_zero_);
  }
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,H->addr(1,0),n+1);
  complex<float> *lj=L->addr();
  complex<float> *Hjjp1=H->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Hjjp1+=n+1) *Hjjp1+=conj(*lj);
  return H;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  S->copy(M);
  float *di=D->addr();
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii+=*di;
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,S->addr(1,0),n+1);
  complex<float> *li=L->addr();
  complex<float> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1+=n+1) *Siip1+=conj(*li);
  return S;
}

template<> SymmetricTridiagonalMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW SymmetricTridiagonalMatrix<float,complex<float> >(n);
  S->copy(*this);
  F77NAME(saxpy)(n,float_mone_,T.D->addr(),1,S->D->addr(),1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.L->addr(),1,S->L->addr(),1);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator-(
const TridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  float *di=D->addr();
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) *Sii=*di;
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),1);
  complex<float> *li=L->addr();
  complex<float> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1=conj(*li);
  F77NAME(caxpy)(n,complex_float_mone_,T.addr(0,0),1,S->addr(0,0),1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(1,0),1,S->addr(1,0),1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(0,1),1,S->addr(0,1),1);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >* operator-(
const TridiagonalMatrix<float,complex<float> > &T,
const SymmetricTridiagonalMatrix<float,complex<float> > &St) {
  int n=St.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  S->copy(T);
  float *di=St.diagonalAddr(0);
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii++) {
    *Sii-=complex<float>(*di,float_zero_);
  }
  F77NAME(caxpy)(n-1,complex_float_mone_,St.lowerDiagonalAddr(0),1,
    S->addr(1,0),1);
  complex<float> *li=St.lowerDiagonalAddr(0);
  complex<float> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1++) *Siip1-=conj(*li);
  return S;
}

template<> SymmetricMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator-(
const SymmetricMatrix<float,complex<float> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricMatrix<float,complex<float> > *T=OPERATOR_NEW
    SymmetricMatrix<float,complex<float> >(n,complex_float_zero_);
  float *di=D->addr();
  complex<float> *Tii=T->addr(0,0);
  for (int i=0;i<n;i++,di++,Tii+=n+1) {
    *Tii=complex<float>(*di,float_zero_);
  }
  F77NAME(ccopy)(n-1,L->addr(),1,T->addr(1,0),n+1);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,complex_float_mone_,S.addr(j,j),1,T->addr(j,j),1);
  }
  return T;
}

template<> SymmetricMatrix<float,complex<float> >* operator-(
const SymmetricMatrix<float,complex<float> > &S,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,complex<float> > *M=
    OPERATOR_NEW SymmetricMatrix<float,complex<float> >(n);
  M->copy(S);
  float *di=T.diagonalAddr(0);
  complex<float> *Mii=M->addr(0,0);
  for (int i=0;i<n;i++,di++,Mii+=n+1) {
    *Mii-=complex<float>(*di,float_zero_);
  }
  F77NAME(caxpy)(n-1,complex_float_mone_,T.lowerDiagonalAddr(0),1,
    M->addr(1,0),n+1);
  return M;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0); 
  float *di=D->addr();
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di;
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  complex<float> *li=L->addr();
  complex<float> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1+=n+1) *Siip1=conj(*li);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(n-j,complex_float_mone_,M.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    Sii=S->addr(0,0);
    for (int i=0;i<n;i++,Sii+=n+1) *Sii-=float_one_;
    for (int j=0;j<n-1;j++) {
      F77NAME(caxpy)(n-j-1,complex_float_mone_,M.addr(j+1,j),1,
        S->addr(j+1,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &M,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0); 
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    complex<float> *Sjj=S->addr(0,0);
    for (int j=0;j<n;j++,Sjj+=n+1) *Sjj=complex_float_one_;
    for (int j=0;j<n-1;j++) {
      F77NAME(ccopy)(n-j-1,M.addr(j+1,j),1,
        S->addr(j+1,j),1);
    }
  }
  float *dj=T.diagonalAddr(0);
  complex<float> *Sjj=S->addr(0,0);
  for (int j=0;j<n;j++,dj++,Sjj+=n+1) *Sjj-=*dj;
  F77NAME(caxpy)(n-1,complex_float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),n+1);
  complex<float> *lj=T.lowerDiagonalAddr(0);
  complex<float> *Sjjp1=S->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Sjjp1+=n+1) *Sjjp1-=conj(*lj);
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
   UpperHessenbergMatrix<float,complex<float> >(n,complex_float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0); 
  float *dj=D->addr();
  complex<float> *Hjj=H->addr(0,0);
  for (int j=0;j<n;j++,dj++,Hjj+=n+1) *Hjj=*dj;
  F77NAME(ccopy)(n-1,L->addr(),1,H->addr(1,0),n+1);
  complex<float> *lj=L->addr(0);
  complex<float> *Hjjp1=H->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Hjjp1+=n+1) *Hjjp1=conj(*lj);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(j+1,complex_float_mone_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  } else {
    Hjj=H->addr(0,0);
    for (int j=0;j<n;j++,Hjj+=n+1) *Hjj-=float_one_;
    for (int j=1;j<n;j++) {
      F77NAME(caxpy)(j,complex_float_mone_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  }
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >* operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
   UpperHessenbergMatrix<float,complex<float> >(n,complex_float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0); 
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(j+1,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    complex<float> *Hjj=H->addr(0,0);
    for (int j=0;j<n;j++,Hjj+=n+1) *Hjj=float_one_;
    for (int j=1;j<n;j++) F77NAME(ccopy)(j,U.addr(0,j),1,H->addr(0,j),1);
  }
  float *dj=T.diagonalAddr(0);
  complex<float> *Hjj=H->addr(0,0);
  for (int j=0;j<n;j++,dj++,Hjj+=n+1) {
    *Hjj-=complex<float>(*dj,float_zero_);
  }
  F77NAME(caxpy)(n-1,complex_float_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  complex<float> *lj=T.lowerDiagonalAddr(0);
  complex<float> *Hjjp1=H->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Hjjp1+=n+1) *Hjjp1-=conj(*lj);
  return H;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  float *dj=D->addr();
  complex<float> *Sjj=S->addr(0,0);
  for (int j=0;j<n;j++,dj++,Sjj+=n+1) *Sjj=*dj;
  F77NAME(ccopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  complex<float> *lj=L->addr();
  complex<float> *Sjjp1=S->addr(0,1);
  for (int j=0;j<n-1;j++,lj++,Sjjp1+=n+1) *Sjjp1=conj(*lj);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n,complex_float_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  S->copy(M);
  float *di=T.diagonalAddr(0);
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii-=*di;
  F77NAME(caxpy)(n-1,complex_float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),n+1);
  complex<float> *li=T.lowerDiagonalAddr(0);
  complex<float> *Siip1=S->addr(0,1);
  for (int i=0;i<n-1;i++,li++,Siip1+=n+1) *Siip1-=conj(*li);
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator*(
const SymmetricMatrix<float,complex<float> > &M) const {
// compute by bordering: note that
// [     tau    , bar(lambda) e_0^T ] [ sigma S^H ]
// [ e_0 lambda ,           T       ] [   S    S  ]
//   = [ tau sigma + bar(lambda) e_0^T s , tau s^H + bar(lambda) e_0^T S ]
//   = [ e_0 lambda sigma +          T s , e_0 lambda s^H +          T S ]
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  (*S)(n-2,n-2)=(*D)[n-2]*     M(n-2,n-2) +conj((*L)[n-2])*M(n-1,n-2);
  (*S)(n-1,n-2)=(*L)[n-2]*     M(n-2,n-2) +     (*D)[n-1] *M(n-1,n-2);
  (*S)(n-2,n-1)=(*D)[n-2]*conj(M(n-1,n-2))+conj((*L)[n-2])*M(n-1,n-1);
  (*S)(n-1,n-1)=(*L)[n-2]*conj(M(n-1,n-2))+     (*D)[n-1] *M(n-1,n-1);
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-3;k>=0;k--) {
    F77NAME(chtmv)(n-k-1,float_one_,L->addr(k+1),D->addr(k+1),
      M.addr(k+1,k),1,float_zero_,S->addr(k+1,k),1); // T s
    complex<float> *ti=t;
    const complex<float> *Mik=M.addr(k+1,k+1);
    for (int i=k+1;i<n;i++,ti++,Mik++) *ti=conj(*Mik);
    F77NAME(caxpy)(n-k-1,conj((*L)[k]),t,1,S->addr(k,k+1),n);
      // bar(lambda) e_0^T S = bar(lambda) ( S e_0 )^H
    (*S)(k,k)=conj((*L)[k])*M(k+1,k); // bar(lambda) e_0^T s
    ti=t;
    Mik=M.addr(k,k);
    for (int i=k;i<n;i++,ti++,Mik++) *ti=conj(*Mik);
    F77NAME(caxpy)(n-k,(*D)[k],t,1,S->addr(k,k),n);
      // tau [ sigma , s^H ]
    F77NAME(caxpy)(n-k,(*L)[k],t,1,S->addr(k+1,k),n);
      // e_0 lamba [ sigma , s^H ]
  }
  delete [] t; t=0;
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SymmetricMatrix<float,complex<float> > &M,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
// compute by bordering: note that
// [ sigma s^H ] [     tau    , bar(lambda) e_0^T ]
// [   s    S  ] [ e_0 lambda ,           T       ]
//   = [ sigma tau + s^H e_0 lambda , sigma bar(lambda) e_0^T + s^H T ]
//   = [   s   tau +  S  e_o lambda ,   s   bar(lambda) e_0^T +  S  T ]
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  (*S)(n-2,n-2)=M(n-2,n-2)*T.diagonalValue(n-2)
               +conj(M(n-1,n-2))*T.lowerDiagonalValue(n-2);
  (*S)(n-1,n-2)=M(n-1,n-2)*T.diagonalValue(n-2)
               +     M(n-1,n-1) *T.lowerDiagonalValue(n-2);
  (*S)(n-2,n-1)=M(n-2,n-2)*T.upperDiagonalValue(n-2)
               +conj(M(n-1,n-2))*T.diagonalValue(n-1);
  (*S)(n-1,n-1)=M(n-1,n-2)*T.upperDiagonalValue(n-2)
               +     M(n-1,n-1) *T.diagonalValue(n-1);
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-3;k>=0;k--) {
    F77NAME(chtmv)(n-k-1,complex_float_one_,T.lowerDiagonalAddr(k+1),
      T.diagonalAddr(k+1),M.addr(k+1,k),1,complex_float_zero_,t,1);
    complex<float> *tj=t;
    complex<float> *Skj=S->addr(k,k+1);
    for (int j=k+1;j<n;j++,tj++,Skj+=n) *Skj=conj(*tj);
      // s^H T = ( T s )^H
    F77NAME(caxpy)(n-k-1,T.lowerDiagonalValue(k),M.addr(k+1,k+1),1,
      S->addr(k+1,k),1); // S e_0 lambda
    (*S)(k,k)=conj(M(k+1,k))*T.lowerDiagonalValue(k); // s^H e_0 lambda
    F77NAME(caxpy)(n-k,T.diagonalValue(k),M.addr(k,k),1,S->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(caxpy)(n-k,T.upperDiagonalValue(k),M.addr(k,k),1,
      S->addr(k,k+1),1);
      // [sigma ] bar(lambda)
      // [  s   ]
  }
  delete [] t; t=0;
  return S;
}

template<> Matrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &M) const {
// compute by columns: note that
// [     S_11        , e_j bar(sigma) e_0^T ] [ 0 ]
// [ e_0 sigma e_j^T ,         S_22         ] [ell]
//   = [ e_j bar(sigma) e_0^T ell ]
//   = [ S_22 ell ]
  int m=M.size(0),n=M.size(1);;
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) (*S)(j-1,j)=upperDiagonalValue(j-1)*M(j,j);
        // sigma e_0^T ell
      if (j+1<m) { // S_22 ell
        F77NAME(chtmv)(m-j,complex_float_one_,L->addr(j),D->addr(j),
          M.addr(j,j),1,complex_float_zero_,S->addr(j,j),1);
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
          F77NAME(chtmv)(m-j-1,complex_float_one_,L->addr(j+1),
            D->addr(j+1),M.addr(j+1,j),1,complex_float_zero_,
            S->addr(j+1,j),1);
        } else (*S)(j+1,j)+=diagonalValue(j+1)*M(j+1,j);
      }
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >* operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &M,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
// compute by columns: note that
// L T e_j = L e_{j-1} T_{j-1,j} + L e_j T_{j,j} + L e_{j+1} T_{j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(m-j+1,T.upperDiagonalValue(j-1),M.addr(j-1,j-1),1,
          S->addr(j-1,j),1);
      }
      F77NAME(caxpy)(m-j,T.diagonalValue(j),M.addr(j,j),1,S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(caxpy)(m-j-1,T.lowerDiagonalValue(j),M.addr(j+1,j+1),1,
          S->addr(j+1,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        (*S)(j-1,j)=T.upperDiagonalValue(j-1);
        F77NAME(caxpy)(m-j,T.upperDiagonalValue(j-1),M.addr(j,j-1),1,
          S->addr(j,j),1);
      }
      (*S)(j,j)+=T.diagonalValue(j);
      if (j<m-1) {
        F77NAME(caxpy)(m-j-1,T.diagonalValue(j),M.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      if (j<n-1) {
        (*S)(j+1,j)+=T.lowerDiagonalValue(j);
        if (j<m-2) {
          F77NAME(caxpy)(m-j-2,T.lowerDiagonalValue(j),M.addr(j+2,j+1),1,
            S->addr(j+2,j),1);
        }
      }
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &M) const {
// compute by columns: note that
// [       S_11      , e_k bar(sigma) e_0^T ] [ u ]= [ S_11 u ]
// [ e_0 sigma e_k^T ,         S_22         ] [   ]= [ e_0 sigma e_k^T u ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<m-1;j++) {
      F77NAME(chtmv)(j+1,complex_float_one_,L->addr(),D->addr(),
        M.addr(0,j),1,complex_float_zero_,S->addr(0,j),1);
      (*S)(j+1,j)=lowerDiagonalValue(j)*M(j,j);
    }
    for (int j=m-1;j<n;j++) {
      F77NAME(chtmv)(m,complex_float_one_,L->addr(),D->addr(),
        M.addr(0,j),1,complex_float_zero_,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(chtmv)(j,complex_float_one_,L->addr(),D->addr(),
          M.addr(0,j),1,complex_float_zero_,S->addr(0,j),1);
        (*S)(j-1,j)+=upperDiagonalValue(j-1);
        (*S)(j,j)=lowerDiagonalValue(j-1)*M(j-1,j);
      }
      (*S)(j,j)+=diagonalValue(j);
      if (j+1<m) (*S)(j+1,j)=lowerDiagonalValue(j);
    }
    for (int j=m;j<n;j++) {
      F77NAME(chtmv)(m,complex_float_one_,L->addr(),D->addr(),
        M.addr(0,j),1,complex_float_zero_,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >* operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &M,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
// compute by columns: note that
// U T e_j = U e_{j-1} T_{j-1,j} + U e_j T_{j,j} + U e_{j+1} T_{j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j,m),T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
          S->addr(0,j),1);
      }
      F77NAME(caxpy)(min(j+1,m),T.diagonalValue(j),M.addr(0,j),1,
        S->addr(0,j),1);
      if (j<n-1) {
        F77NAME(caxpy)(min(j+2,m),T.lowerDiagonalValue(j),M.addr(0,j+1),1,
          S->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(min(j-1,m),T.upperDiagonalValue(j-1),
          M.addr(0,j-1),1,S->addr(0,j),1);
        if (j<=m) (*S)(j-1,j)=T.upperDiagonalValue(j-1);
      }
      F77NAME(caxpy)(min(j,m),T.diagonalValue(j),M.addr(0,j),1,
        S->addr(0,j),1);
      if (j<m) (*S)(j,j)+=T.diagonalValue(j);
      if (j<n-1) {
        F77NAME(caxpy)(min(j+1,m),T.lowerDiagonalValue(j),M.addr(0,j+1),1,
          S->addr(0,j),1);
        if (j<m-1) (*S)(j+1,j)+=T.lowerDiagonalValue(j);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator*(
const SquareMatrix<float,complex<float> > &M) const {
  int n=M.size(0);
  CHECK_SAME(n,size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(chtmv)(n,float_one_,L->addr(),D->addr(),
      M.addr(0,j),1,float_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SquareMatrix<float,complex<float> > &M,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(caxpy)(n,T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(caxpy)(n,T.diagonalValue(j),M.addr(0,j),1,S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(caxpy)(n,T.lowerDiagonalValue(j),M.addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(chtmv)(m,float_one_,L->addr(),D->addr(),
      M.addr(0,j),1,float_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,complex<float> > *S=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(caxpy)(m,T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(caxpy)(m,T.diagonalValue(j),M.addr(0,j),1,S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(caxpy)(m,T.lowerDiagonalValue(j),M.addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> Vector<float,complex<float> >*
SymmetricTridiagonalMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<float,complex<float> > *w=OPERATOR_NEW
    Vector<float,complex<float> >(n,complex_float_zero_);
  F77NAME(chtmv)(n,complex_float_one_,L->addr(),D->addr(),
    v.addr(),1,float_zero_,w->addr(),1);
  return w;
}

template<> float
SymmetricTridiagonalMatrix<float,complex<float> >::normFrobenius()
const {
  return F77NAME(clanht)('F',dim,D->addr(),L->addr());
}

template<> float
SymmetricTridiagonalMatrix<float,complex<float> >::normInfinity()
const {
  return F77NAME(clanht)('I',dim,D->addr(),L->addr());
}

template<> float
SymmetricTridiagonalMatrix<float,complex<float> >::normMaxEntry()
const {
  return F77NAME(clanht)('M',dim,D->addr(),L->addr());
}

template<> float
SymmetricTridiagonalMatrix<float,complex<float> >::normOne() const {
  return F77NAME(clanht)('O',dim,D->addr(),L->addr());
}

template<> float SymmetricTridiagonalMatrix<float,complex<float> >
::reciprocalConditionNumber(char norm) const {
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *UF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *UF2=
    OPERATOR_NEW Vector<float,complex<float> >(max(0,dim-2));

  complex<float> *zd=OPERATOR_NEW_BRACKET(complex<float>,dim);
  const float *di=diagonalAddr(0);
  complex<float> *zdi=zd;
  for (int i=0;i<dim;i++,di++,zdi++) {
    *zdi=complex<float>(*di,float_zero_);
  }
  LF->copy(*L);
  UF->copy(*L);
  float anorm=F77NAME(clangt)(norm,dim,L->addr(),zd,L->addr());

  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(cgttrf)(dim,LF->addr(),zd,UF->addr(),UF2->addr(),ipiv,info);
  float rcond=numeric_limits<float>::infinity();
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*dim);
  F77NAME(cgtcon)(norm,dim,LF->addr(),zd,UF->addr(),UF2->addr(),
    ipiv,anorm,rcond,work,info);
  delete [] work; work=0;
  delete [] ipiv; ipiv=0;
  delete UF2; UF2=0;
  delete UF; UF=0;
  delete [] zd; zd=0;
  delete LF; LF=0;
  return rcond;
}

template<> void SymmetricTridiagonalMatrix<float,complex<float> >::stmv(
complex<float> alpha,const Vector<float,complex<float> > &x,
complex<float> beta,Vector<float,complex<float> > &b) const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  F77NAME(chtmv)(n,alpha,L->addr(),D->addr(),x.addr(),1,beta,b.addr(),1);
}

template<> void SymmetricTridiagonalMatrix<float,complex<float> >::stmm(
complex<float> alpha,const Matrix<float,complex<float> > &X,
complex<float> beta,Matrix<float,complex<float> > &B,char side) const {
  if (side=='L' || side=='l') {
    int m=size(0),n=X.size(1);
    CHECK_SAME(n,B.size(1));
    CHECK_SAME(m,X.size(0));
    CHECK_SAME(m,B.size(0));
    for (int j=0;j<n;j++) {
      F77NAME(chtmv)(m,alpha,L->addr(),D->addr(),X.addr(0,j),1,beta,
        B.addr(0,j),1);
    }
  } else {
    int m=X.size(0),n=X.size(1);
    CHECK_SAME(n,size(0));
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,B.size(1));
    if (abs(beta)==float_zero_) B=complex_float_zero_;
    else if (beta!=float_one_) B*=beta;
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(caxpy)(m,upperDiagonalValue(j-1)*alpha,X.addr(0,j-1),1,
          B.addr(0,j),1);
      }
      F77NAME(caxpy)(m,diagonalValue(j)*alpha,X.addr(0,j),1,
        B.addr(0,j),1);
      if (j<n-1) {
        F77NAME(caxpy)(m,lowerDiagonalValue(j)*alpha,X.addr(0,j+1),1,
          B.addr(0,j),1);
      }
    }
  }
}

template<> void
SymmetricTridiagonalMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,Vector<float,
complex<float> > &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *UF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *DF=
    OPERATOR_NEW Vector<float,complex<float> >(dim);
  float *di=D->addr();
  complex<float> *dfi=DF->addr();
  for (int i=0;i<dim;i++,di++,dfi++) {
    *dfi=complex<float>(*di,float_zero_);
  }
  const complex<float> *li=L->addr();
  complex<float> *ufi=UF->addr();
  for (int i=0;i<dim-1;i++,li++,ufi++) *ufi=conj(*li);
  LF->copy(*L);
  x.copy(b);
  int info;
  F77NAME(cgtsv)(dim,1,LF->addr(),DF->addr(),UF->addr(),x.addr(),dim,
    info);
  CHECK_SAME(info,0);
  delete DF; DF=0;
  delete LF; LF=0;
  delete UF; UF=0;
}

template<> void
SymmetricTridiagonalMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *UF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,complex<float> > *DF=
    OPERATOR_NEW Vector<float,complex<float> >(dim);
  float *di=D->addr();
  complex<float> *dfi=DF->addr();
  for (int i=0;i<dim;i++,di++,dfi++) {
    *dfi=complex<float>(*di,float_zero_);
  }
  const complex<float> *li=L->addr();
  complex<float> *ufi=UF->addr();
  for (int i=0;i<dim-1;i++,li++,ufi++) *ufi=conj(*li);
  LF->copy(*L);
  int info;
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    X.copy(B);
    F77NAME(cgtsv)(m,n,LF->addr(),DF->addr(),UF->addr(),X.addr(),m,info);
    CHECK_SAME(info,0);
  } else {
    CHECK_SAME(dim,n);
    complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
    Vector<float,complex<float> > *UF2=
      OPERATOR_NEW Vector<float,complex<float> >(dim-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
    F77NAME(cgttrf)(n,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),ipiv,
      info);
    CHECK_TEST(info==0);
    for (int i=0;i<m;i++) {
      F77NAME(ccopy)(n,B.addr(i,0),m,t,1);
      F77NAME(cgttrs)('T',n,1,LF->addr(),DF->addr(),UF->addr(),
        UF2->addr(),ipiv,t,n,info);
      CHECK_TEST(info==0);
      F77NAME(ccopy)(n,t,1,X.addr(i,0),m);
    }
    delete [] ipiv; ipiv=0;
    delete UF2; UF2=0;
    delete [] t; t=0;
  }
  delete DF; DF=0;
  delete LF; LF=0;
  delete UF; UF=0;
}

template class SymmetricTridiagonalMatrix<float,complex<float> >;
template void testSymmetricTridiagonalMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveMatrix<float,complex<float> >*
SymmetricPositiveTridiagonalMatrix<float,complex<float> >::operator+(
const SymmetricPositiveMatrix<float,complex<float> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricPositiveMatrix<float,complex<float> > *T=
    OPERATOR_NEW SymmetricPositiveMatrix<float,complex<float> >(n);
  T->copy(S);
  float *di=D->addr();
  complex<float> *Tii=T->addr(0,0);
  for (int i=0;i<n;i++,di++,Tii+=n+1) {
    *Tii+=complex<float>(*di,float_zero_);
  }
  F77NAME(caxpy)(n-1,complex_float_one_,L->addr(),1,T->addr(1,0),n+1);
  return T;
}

template<> float
SymmetricPositiveTridiagonalMatrix<float,complex<float> >
::reciprocalConditionNumber(char norm) const {
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  LF->copy(*L);
  DF->copy(*D);
  int info;
  F77NAME(cpttrf)(dim,DF->addr(),LF->addr(),info);
  CHECK_TEST(info==0);
  float rcond=numeric_limits<float>::infinity();
  float *work=OPERATOR_NEW_BRACKET(float,2*dim);
  float anorm=normOne();
  F77NAME(cptcon)(dim,DF->addr(),LF->addr(),anorm,rcond,work,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete DF; DF=0;
  delete LF; LF=0;
  return rcond;
}

/*
template<> Vector<float,float>*
SymmetricPositiveTridiagonalMatrix<float,float>::eigenvalues(
OrthogonalMatrix<float,float> *&Q) const {
  if (Q!=0) CHECK_SAME(dim,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Vector<float,float> *lambda =OPERATOR_NEW Vector<float,float>(dim);
  lambda->copy(*D);
  Vector<float,float> *L_copy =OPERATOR_NEW Vector<float,float>(dim-1);
  L_copy->copy(*L);
  float *work=OPERATOR_NEW_BRACKET(float,2*dim-2);
  int info;
  float *qa=( Q==0 ? 0 : Q->addr() );
  F77NAME(dstev)(jobz,dim,lambda->addr(),L_copy->addr(),qa,dim,work,info);
  delete [] work; work=0;
  delete L_copy; L_copy=0;
  return lambda;
}
*/

template<> void
SymmetricPositiveTridiagonalMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  LF->copy(*L);
  DF->copy(*D);
  x.copy(b);
  int info;
  F77NAME(cptsv)(dim,1,DF->addr(),LF->addr(),x.addr(),dim,info);
  CHECK_SAME(info,0);
  delete DF; DF=0;
  delete LF; LF=0;
}

template<> void
SymmetricPositiveTridiagonalMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<float,complex<float> > *LF=
    OPERATOR_NEW Vector<float,complex<float> >(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  LF->copy(*L);
  DF->copy(*D);
  int info;
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    X.copy(B);
    F77NAME(cptsv)(m,n,DF->addr(),LF->addr(),X.addr(),m,info);
    CHECK_SAME(info,0);
  } else {
    CHECK_SAME(dim,n);
    complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
    F77NAME(cpttrf)(n,DF->addr(),LF->addr(),info);
    CHECK_SAME(info,0);
    for (int i=0;i<m;i++) {
      const complex<float> *Bij=B.addr(i,0);
      complex<float> *tj=t;
      for (int j=0;j<n;j++,Bij++,tj++) *tj=conj(*Bij);
      F77NAME(cpttrs)('L',n,1,DF->addr(),LF->addr(),t,n,info);
      CHECK_SAME(info,0);
      complex<float> *Xij=X.addr(i,0);
      tj=t;
      for (int j=0;j<n;j++,Xij++,tj++) *Xij=conj(*tj);
    }
    delete [] t; t=0;
  }
  delete DF; DF=0;
  delete LF; LF=0;
}

template class
  SymmetricPositiveTridiagonalMatrix<float,complex<float> >;
template void
  testSymmetricPositiveTridiagonalMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SquareMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  F77NAME(ccopy)(n,D->addr(),1,M->addr(0,0),n+1);
  return M;
}

template<> SymmetricTridiagonalMatrix<float,complex<float> >*
operator+(const DiagonalMatrix<float,float> &A,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<float,complex<float> >(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_one_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >*
operator+(const DiagonalMatrix<float,complex<float> > &A,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  F77NAME(ccopy)(n-1,T.lowerDiagonalAddr(0),1,S->lowerDiagonalAddr(),1);
  const complex<float> *di=A.addr();
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Sii=S->diagonalAddr();
  for (int i=0;i<n;i++,di++,Tii++,Sii++) *Sii=*di+*Tii;
  const complex<float> *Tip1i=T.lowerDiagonalAddr(0);
  complex<float> *Siip1=S->upperDiagonalAddr();
  for (int i=0;i<n-1;i++,Tip1i++,Siip1++) *Siip1=conj(*Tip1i);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator+(
const TridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=OPERATOR_NEW
    TridiagonalMatrix<float,complex<float> >(n);
  S->copy(T);
  F77NAME(caxpy)(n,complex_float_one_,addr(),1,S->addr(0,0),1);
  return S;
}

template<> SymmetricMatrix<float,complex<float> >* operator+(
const DiagonalMatrix<float,float> &A,
const SymmetricMatrix<float,complex<float> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,complex<float> > *S=
    OPERATOR_NEW SymmetricMatrix<float,complex<float> >(n);
  S->copy(T);
  const float *di=A.addr();
  complex<float> *Sii=S->addr();
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii+=*di;
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator+(
const DiagonalMatrix<float,complex<float> > &A,
const SymmetricMatrix<float,complex<float> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,T.addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=A[j];
    if (j>0) {
      const complex<float> *Tji=T.addr(j,0);
      complex<float> *Sij=S->addr(0,j);
      for (int i=0;i<j;i++,Tji+=n,Sij++) *Sij=conj(*Tji);
    }
  }
  return S;
}

template<> UpperTriangularMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<float,complex<float> > *S=OPERATOR_NEW
    UpperTriangularMatrix<float,complex<float> >(n);
  S->copy(U);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    F77NAME(caxpy)(n,complex_float_one_,addr(),1,S->addr(),n+1);
  } else {
    const complex<float> *di=D->addr();
    complex<float> *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di+float_one_;
  }
  return S;
}

template<> LowerTriangularMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<float,complex<float> > *S=OPERATOR_NEW
    LowerTriangularMatrix<float,complex<float> >(n);
  S->copy(L);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    F77NAME(caxpy)(n,float_one_,addr(),1,S->addr(),n+1);
  } else {
    const complex<float> *di=D->addr();
    complex<float> *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di+float_one_;
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n);
  S->copy(M);
  F77NAME(caxpy)(n,complex_float_one_,addr(),1,S->addr(),n+1);
  return S;
}

template<> SymmetricTridiagonalMatrix<float,complex<float> >* operator-(
const DiagonalMatrix<float,float> &A,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<float,complex<float> >(n);
  S->copy(T);
  *S*=float_mone_;
  F77NAME(saxpy)(n,float_one_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >* operator-(
const DiagonalMatrix<float,complex<float> > &A,
const SymmetricTridiagonalMatrix<float,complex<float> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  const complex<float> *di=A.addr();
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Sii=S->diagonalAddr();
  for (int i=0;i<n;i++,di++,Tii++,Sii++) *Sii=*di-*Tii;
  const complex<float> *Tip1i=T.lowerDiagonalAddr(0);
  complex<float> *Siip1=S->upperDiagonalAddr();
  complex<float> *Sip1i=S->lowerDiagonalAddr();
  for (int i=0;i<n-1;i++,Tip1i++,Siip1++,Sip1i++) {
    *Sip1i=-*Tip1i;
    *Siip1=-conj(*Tip1i);
  }
  return S;
}

template<> SymmetricTridiagonalMatrix<float,complex<float> >* operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<float,complex<float> >(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_mone_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >* operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const DiagonalMatrix<float,complex<float> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  F77NAME(ccopy)(n-1,T.lowerDiagonalAddr(0),1,S->lowerDiagonalAddr(),1);
  const complex<float> *di=A.addr();
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Sii=S->diagonalAddr();
  for (int i=0;i<n;i++,di++,Tii++,Sii++) *Sii=*Tii-*di;
  const complex<float> *Tip1i=T.lowerDiagonalAddr(0);
  complex<float> *Siip1=S->upperDiagonalAddr();
  for (int i=0;i<n-1;i++,Tip1i++,Siip1++) *Siip1=conj(*Tip1i);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator-(
const TridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=OPERATOR_NEW
    TridiagonalMatrix<float,complex<float> >(n);
  S->copy(T);
  *S*=float_mone_;
  F77NAME(caxpy)(n,complex_float_one_,addr(),1,S->addr(0,0),1);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >* operator-(
const TridiagonalMatrix<float,complex<float> > &T,
const DiagonalMatrix<float,complex<float> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=OPERATOR_NEW
    TridiagonalMatrix<float,complex<float> >(n);
  S->copy(T);
  F77NAME(caxpy)(n,complex_float_mone_,A.addr(),1,S->addr(0,0),1);
  return S;
}

template<> SymmetricMatrix<float,complex<float> >* operator-(
const DiagonalMatrix<float,float> &A,
const SymmetricMatrix<float,complex<float> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricMatrix<float,complex<float> >(n);
  S->copy(T);
  *S*=complex_float_mone_;
  const float *di=A.addr();
  complex<float> *Sii=S->addr();
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii+=*di;
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const DiagonalMatrix<float,complex<float> > &A,
const SymmetricMatrix<float,complex<float> > &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    (*S)(j,j)=A[j]-T(j,j);
    if (j+1<n) {
      const complex<float> *Tij=T.addr(j+1,j);
      complex<float> *Sij=S->addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Tij++,Sij++,Sji+=n) {
        *Sij=-*Tij;
        *Sji=-conj(*Tij);
      }
    }
  }
  return S;
}

template<> SymmetricMatrix<float,complex<float> >* operator-(
const SymmetricMatrix<float,complex<float> > &T,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricMatrix<float,complex<float> >(n);
  S->copy(T);
  const float *di=A.addr();
  complex<float> *Sii=S->addr();
  for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii-=*di;
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const SymmetricMatrix<float,complex<float> > &T,
const DiagonalMatrix<float,complex<float> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    (*S)(j,j)=T(j,j)-A[j];
    if (j+1<n) {
      const complex<float> *Tij=T.addr(j+1,j);
      complex<float> *Sij=S->addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Tij++,Sij++,Sji+=n) {
        *Sij=*Tij;
        *Sji=conj(*Tij);
      }
    }
  }
  return S;
}

template<> UpperTriangularMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<float,complex<float> > *S=OPERATOR_NEW
    UpperTriangularMatrix<float,complex<float> >(n);
  S->copy(U);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U) ==0) {
    *S*=float_mone_;
    F77NAME(caxpy)(n,complex_float_one_,addr(),1,S->addr(),n+1);
  } else {
    const complex<float> *dj=D->addr();
    for (int j=0;j<n;j++,dj++) {
      complex<float> *Sij=S->addr(0,j);
      const complex<float> *Uij=U.addr(0,j);
      for (int i=0;i<j;i++,Sij++,Uij++) *Sij=-*Uij;
      *Sij=*dj-float_one_;
    }
  }
  return S;
}

template<> UpperTriangularMatrix<float,complex<float> >* operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const DiagonalMatrix<float,complex<float> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<float,complex<float> > *S=OPERATOR_NEW
    UpperTriangularMatrix<float,complex<float> >(n);
  S->copy(U);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    F77NAME(caxpy)(n,complex_float_mone_,A.addr(),1,S->addr(),n+1);
  } else {
    const complex<float> *di=A.addr();
    complex<float> *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=float_one_-*di;
  }
  return S;
}

template<> LowerTriangularMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<float,complex<float> > *S=OPERATOR_NEW
    LowerTriangularMatrix<float,complex<float> >(n);
  S->copy(L);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    *S*=complex_float_mone_;
    F77NAME(caxpy)(n,complex_float_one_,addr(),1,S->addr(),n+1);
  } else {
    const complex<float> *dj=D->addr();
    for (int j=0;j<n;j++,dj++) {
      complex<float> *Sij=S->addr(j,j);
      *Sij=*dj-float_one_;
      if (j+1<n) {
        Sij++;
        const complex<float> *Lij=L.addr(j+1,j);
        for (int i=j+1;i<n;i++,Sij++,Lij++) *Sij=-*Lij;
      }
    }
  }
  return S;
}

template<> LowerTriangularMatrix<float,complex<float> >* operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const DiagonalMatrix<float,complex<float> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<float,complex<float> > *S=
    OPERATOR_NEW LowerTriangularMatrix<float,complex<float> >(n);
  S->copy(L);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    F77NAME(caxpy)(n,complex_float_mone_,A.addr(),1,S->addr(),n+1);
  } else {
    const complex<float> *di=A.addr();
    complex<float> *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=float_one_-*di;
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  S->copy(M);
  *S*=complex_float_mone_;
  F77NAME(caxpy)(n,complex_float_one_,addr(),1,S->addr(),n+1);
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const DiagonalMatrix<float,complex<float> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  S->copy(M);
  F77NAME(caxpy)(n,complex_float_mone_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >* 
DiagonalMatrix<float,complex<float> >::operator*(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  for (int i=0;i<n;i++) {
    complex<float> di=(*D)[i];
    if (i>0) (*S)(i,i-1)=di*T.lowerDiagonalValue(i-1);
    (*S)(i,i)=di*T.diagonalValue(i);
    if (i<n-1) (*S)(i,i+1)=di*conj(T.lowerDiagonalValue(i));
  }
  return S;
}

template<> TridiagonalMatrix<float,complex<float> >* operator*(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const DiagonalMatrix<float,complex<float> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,complex<float> > *S=
    OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    complex<float> dj=A[j];
    if (j>0) (*S)(j-1,j)=conj(T.lowerDiagonalValue(j-1))*dj;
    (*S)(j,j)=T.diagonalValue(j)*dj;
    if (j<n-1) (*S)(j+1,j)=T.lowerDiagonalValue(j)*dj;
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
DiagonalMatrix<float,complex<float> >::operator*(
const SymmetricMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(dim);
  for (int i=0;i<n;i++) {
    complex<float> di=D->operator[](i);
    if (i>0) {
      const complex<float> *Tij=T.addr(i,0);
      complex<float> *Sij=S->addr(i,0);
      for (int j=0;j<i;j++,Tij+=n,Sij+=n) *Sij=di*(*Tij);
    }
    const complex<float> *Tji=T.addr(i,i);
    complex<float> *Sij=S->addr(i,i);
    for (int j=i;j<n;j++,Tji++,Sij+=n) *Sij=di*conj(*Tji);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SymmetricMatrix<float,complex<float> > &T,
const DiagonalMatrix<float,complex<float> > &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    complex<float> dj=A[j];
    if (j>0) {
      const complex<float> *Tji=T.addr(j,0);
      complex<float> *Sij=S->addr(0,j);
      for (int i=0;i<j;i++,Tji+=n,Sij++) *Sij=conj(*Tji)*dj;
    }
    const complex<float> *Tij=T.addr(j,j);
    complex<float> *Sij=S->addr(j,j);
    for (int i=j;i<n;i++,Tij++,Sij++) *Sij=(*Tij)*dj;
  }
  return S;
}

template class DiagonalMatrix<float,complex<float> >;
template void testDiagonalMatrix(float,complex<float> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const complex<float>
  UpperHessenbergMatrix<float,complex<float> >::outofbounds_ =
  complex_float_zero_;
template<> const complex<float> UpperHessenbergMatrix<float,
  complex<float> >::undefined_(numeric_limits<float>::infinity(),
    numeric_limits<float>::infinity());
template<> complex<float>
  UpperHessenbergMatrix<float,complex<float> >::safety_ =
  complex_float_zero_;

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(n,j+2),addr(0,j),1,M->addr(0,j),1);
  }
  return M;
}

template<> UpperHessenbergMatrix<float,complex<float> >&
UpperHessenbergMatrix<float,complex<float> >::operator+=(
const UpperHessenbergMatrix<float,complex<float> > &H) {
  int n=size(0);
  CHECK_SAME(n,H.size(0));
  complex<float> *colj=addr();
  const complex<float> *Hcolj=H.addr();
  for (int j=0;j<n;j++,colj+=n,Hcolj+=n) {
    F77NAME(caxpy)(min(n,j+2),complex_float_one_,Hcolj,1,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<float,complex<float> >&
UpperHessenbergMatrix<float,complex<float> >::operator-=(
const UpperHessenbergMatrix<float,complex<float> > &H) {
  int n=size(0);
  CHECK_SAME(n,H.size(0));
  complex<float> *colj=addr();
  const complex<float> *Hcolj=H.addr();
  for (int j=0;j<n;j++,colj+=n,Hcolj+=n) {
    F77NAME(caxpy)(min(n,j+2),complex_float_mone_,Hcolj,1,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<float,complex<float> >&
UpperHessenbergMatrix<float,complex<float> >::operator*=(
complex<float> scalar) {
  int n=size(0);
  complex<float> *colj=addr();
  for (int j=0;j<n;j++,colj+=n) {
    F77NAME(cscal)(min(n,j+2),scalar,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<float,complex<float> >&
UpperHessenbergMatrix<float,complex<float> >::operator/=(
complex<float> scalar) {
  int n=size(0);
  complex<float> *colj=addr();
  complex<float> s=complex_float_one_/scalar;
  for (int j=0;j<n;j++,colj+=n) {
    F77NAME(cscal)(min(n,j+2),s,colj,1);
  }
  return *this;
}

template<> void UpperHessenbergMatrix<float,complex<float> >::copy(
const Matrix<float,complex<float> > &M) {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(j+2,n),M.addr(0,j),1,addr(0,j),1);
  }
}

/*
template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::transpose() const {
  int n=size(0);
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(dcopy)(min(j+2,n),addr(0,j),1,S->addr(j,0),n);
  }
  return S;
}
*/

template<> void UpperHessenbergMatrix<float,complex<float> >::copyFrom(
int m,const SquareMatrix<float,complex<float> > &S) {
  m=min(m,size(0));
  m=min(m,S.size(0));
  for (int j=0;j<m;j++) {
    F77NAME(ccopy)(min(j+2,m),S.addr(0,j),1,addr(0,j),1);
  }
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator+(
const DiagonalMatrix<float,complex<float> > &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
      complex_float_zero_);
  H->copy(*this);
  F77NAME(caxpy)(n,complex_float_one_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator+(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(*this);
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Hii=H->addr(0,0);
  for (int i=0;i<n;i++,Tii++,Hii+=n+1) *Hii+=*Tii;
  F77NAME(caxpy)(n-1,complex_float_one_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  const complex<float> *Tip1i=T.lowerDiagonalAddr(0);
  complex<float> *Hiip1=H->addr(0,1);
  for (int i=0;i<n-1;i++,Tip1i++,Hiip1+=n+1) *Hiip1+=conj(*Tip1i);
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator+(
const TridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(*this);
  F77NAME(caxpy)(n,complex_float_one_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator+(
const SymmetricMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      const complex<float> *Mij=M.addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Mij++,Sji+=n) *Sji=conj(*Mij);
    }
    F77NAME(caxpy)(min(j+2,n),float_one_,addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      F77NAME(caxpy)(n-j,complex_float_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)+=complex_float_one_;
      if (j<n-1) {
        F77NAME(caxpy)(n-j-1,complex_float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(*this);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(j+1,complex_float_one_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)+=complex_float_one_;
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_one_,U.addr(0,j),1,
          H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(j+2,n),complex_float_one_,addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator-(
const DiagonalMatrix<float,complex<float> > &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(*this);
  F77NAME(caxpy)(n,complex_float_mone_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >* operator-(
const DiagonalMatrix<float,complex<float> > &D,
const UpperHessenbergMatrix<float,complex<float> > &U) {
  int n=U.size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(U);
  *H*=complex_float_mone_;
  F77NAME(caxpy)(n,complex_float_one_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(*this);
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Hii=H->addr(0,0);
  for (int i=0;i<n;i++,Tii++,Hii+=n+1) *Hii-=*Tii;
  F77NAME(caxpy)(n-1,complex_float_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  const complex<float> *Tip1i=T.lowerDiagonalAddr(0);
  complex<float> *Hiip1=H->addr(0,1);
  for (int i=0;i<n-1;i++,Tip1i++,Hiip1+=n+1) *Hiip1-=conj(*Tip1i);
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >* operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const UpperHessenbergMatrix<float,complex<float> > &U) {
  int n=U.size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(U);
  *H*=complex_float_mone_;
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Hii=H->addr(0,0);
  for (int i=0;i<n;i++,Tii++,Hii+=n+1) *Hii+=*Tii;
  F77NAME(caxpy)(n-1,complex_float_one_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  const complex<float> *Tip1i=T.lowerDiagonalAddr(0);
  complex<float> *Hiip1=H->addr(0,1);
  for (int i=0;i<n-1;i++,Tip1i++,Hiip1+=n+1) *Hiip1+=conj(*Tip1i);
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator-(
const TridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(*this);
  F77NAME(caxpy)(n,complex_float_mone_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_mone_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >* operator-(
const TridiagonalMatrix<float,complex<float> > &T,
const UpperHessenbergMatrix<float,complex<float> > &U) {
  int n=U.size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(U);
  *H*=complex_float_mone_;
  F77NAME(caxpy)(n,complex_float_one_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(caxpy)(n-1,complex_float_one_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator-(
const SymmetricMatrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(n-j,float_mone_,M.addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      const complex<float> *Mij=M.addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Mij++,Sji+=n) *Sji-=conj(*Mij);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const SymmetricMatrix<float,complex<float> > &M,
const UpperHessenbergMatrix<float,complex<float> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      const complex<float> *Mij=M.addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,Mij++,Sji+=n) *Sji=conj(*Mij);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(j+2,n),complex_float_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      F77NAME(caxpy)(n-j,complex_float_mone_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)-=complex_float_one_;
      if (j<n-1) {
        F77NAME(caxpy)(n-j-1,complex_float_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const UpperHessenbergMatrix<float,complex<float> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      F77NAME(caxpy)(min(j+2,n),complex_float_mone_,H.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=complex_float_one_;
      if (j<n-1) {
        F77NAME(ccopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      F77NAME(caxpy)(min(j+2,n),complex_float_mone_,H.addr(0,j),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,
    complex_float_zero_);
  H->copy(*this);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(j+1,complex_float_mone_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)-=complex_float_one_;
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_mone_,U.addr(0,j),1,
          H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> UpperHessenbergMatrix<float,complex<float> >* operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const UpperHessenbergMatrix<float,complex<float> > &M) {
  int n=M.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,complex<float> > *H=OPERATOR_NEW
    UpperHessenbergMatrix<float,complex<float> >(n,float_zero_);
  H->copy(M);
  *H*=complex_float_mone_;
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(j+1,complex_float_one_,U.addr(0,j),1,
        H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)+=complex_float_one_;
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_one_,U.addr(0,j),1,
          H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
    F77NAME(caxpy)(n,complex_float_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const UpperHessenbergMatrix<float,complex<float> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n,M.addr(0,j),1,S->addr(0,j),1);
    F77NAME(caxpy)(min(j+2,n),complex_float_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const UpperHessenbergMatrix<float,complex<float> > &H2) const {
  int n=size(0);
  CHECK_SAME(n,H2.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<=min(j+1,n-1);k++) {
      F77NAME(caxpy)(min(k+2,n),H2(k,j),addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const DiagonalMatrix<float,complex<float> > &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,complex<float> > *S=
    OPERATOR_NEW UpperHessenbergMatrix<float,complex<float> >(n);
  S->copy(*this);
  for (int j=0;j<n;j++) {
    F77NAME(cscal)(min(j+2,n),D[j],S->addr(0,j),1);
  }
  return S;
}

template<> UpperHessenbergMatrix<float,complex<float> >* operator*(
const DiagonalMatrix<float,complex<float> > &D,
const UpperHessenbergMatrix<float,complex<float> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,complex<float> > *S=
    OPERATOR_NEW UpperHessenbergMatrix<float,complex<float> >(n);
  S->copy(H);
  for (int i=0;i<n;i++) {
    int j=max(0,i-1);
    F77NAME(cscal)(n-j,D[i],S->addr(i,j),n);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
//compute by columns: note that
// H T e_j = H e_{j-1} T_{j-1,j} + H e_j T_{j,j} + H e_{j+1} T_{j+1,j}
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(caxpy)(min(j+1,n),T.upperDiagonalValue(j-1),addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(caxpy)(min(j+2,n),T.diagonalValue(j),addr(0,j),1,
      S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(caxpy)(min(j+3,n),T.lowerDiagonalValue(j),addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const UpperHessenbergMatrix<float,complex<float> > &H) {
// compute by rows: note that
// e_i^T T H
//   = T_{i,i-1} e_{i-1}^T H + T_{i,i} e_i^T H + T_{i,i+1} e_{i+1}^T H
  int n=H.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int i=0;i<n;i++) {
    if (i>0) {
      int j=max(0,i-2);
      F77NAME(caxpy)(n-j,T.lowerDiagonalValue(i-1),H.addr(i-1,j),n,
        S->addr(i,j),n);
    }
    int j=max(0,i-1);
    F77NAME(caxpy)(n-j,T.diagonalValue(i),H.addr(i,j),n,S->addr(i,j),n);
    if (i<n-1) {
      F77NAME(caxpy)(n-i,T.upperDiagonalValue(i),H.addr(i+1,i),n,
        S->addr(i,i),n);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const TridiagonalMatrix<float,complex<float> > &T) const {
//compute by columns: note that
// H T e_j = H e_{j-1} T_{j-1,j} + H e_j T_{j,j} + H e_{j+1} T_{j+1,j}
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(j-1,0);k<=min(j+1,n-1);k++) {
      F77NAME(caxpy)(min(k+2,n),T(k,j),addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const TridiagonalMatrix<float,complex<float> > &T,
const UpperHessenbergMatrix<float,complex<float> > &H) {
// compute by rows: note that
// e_i^T T H
//   = T_{i,i-1} e_{i-1}^T H + T_{i,i} e_i^T H + T_{i,i+1} e_{i+1}^T H
  int n=H.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int i=0;i<n;i++) {
    for (int k=max(i-1,0);k<=min(i+1,n-1);k++) {
      int j=max(0,k-1);
      F77NAME(caxpy)(min(n-k+1,n),T(i,k),H.addr(k,j),n,S->addr(i,j),n);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const SymmetricMatrix<float,complex<float> > &S) const {
// compute by bordering: note that
// [     eta_11 h^T ] [ sigma s^H ]
// [ e_0 eta_21  H  ] [   s    S  ]
//   = [     eta_11 sigma + h^T s ,     eta_11 s^H + h^T S ]
//   = [ e_0 eta_21 sigma +  H  s , e_0 eta_21 s^H +  H  S ]
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  (*M)(n-1,n-1)=(*this)(n-1,n-1)*S(n-1,n-1);
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  complex<float> *v=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-2;k>=0;k--) {
    for (int j=k+1;j<n;j++) { // H s
      F77NAME(caxpy)(min(n-k-1,j-k+1),S(j,k),addr(k+1,j),1,
        M->addr(k+1,k),1);
    }
    const complex<float> *hkj=addr(k,k+1);
    complex<float> *tj=t;
    for (int j=k+1;j<n;j++,hkj+=n,tj++) *tj=conj(*hkj);
    F77NAME(chemv)('L',n-k-1,float_one_,S.addr(k+1,k+1),n,
      t,1,float_zero_,v,1);
    complex<float> *vj=v;
    complex<float> *Mkj=M->addr(k,k+1);
    for (int j=k+1;j<n;j++,vj++,Mkj+=n) *Mkj=conj(*vj);
      //h^T S=(S bar(h))^H
    (*M)(k,k)=F77NAME(cdotu)(n-k-1,addr(k,k+1),n,S.addr(k+1,k),1);//h^T s
    F77NAME(cgerc)(2,n-k,complex_float_one_,addr(k,k),1,S.addr(k,k),1,
      M->addr(k,k),n);
  }
  delete t; t=0;
  delete v; v=0;
  return M;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SymmetricMatrix<float,complex<float> > &S,
const UpperHessenbergMatrix<float,complex<float> > &H) {
// compute by bordering: note that
// [ sigma s^H ] [     eta_11 h^T ]
// [   s    S  ] [ e_0 eta_21  H  ]
//   = [ sigma eta_11 + s^H e_0 eta_21 , sigma h^T + s^H H ]
//   = [   s   eta_11 +  S  e_0 eta_21 ,   s   h^T +  S  H ]
  int n=H.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  (*M)(n-1,n-1)=S(n-1,n-1)*H(n-1,n-1);
  for (int k=n-2;k>=0;k--) { // s^T H
    F77NAME(caxpy)(n-k-1,H(k+1,k),S.addr(k+1,k+1),1,M->addr(k+1,k),1);
      // S e_0 eta_21
    for (int j=k+1;j<n;j++) { // s^H H
      (*M)(k,j)=
        F77NAME(cdotc)(min(n-k-1,j-k+1),S.addr(k+1,k),1,H.addr(k+1,j),1);
    }
    (*M)(k,k)=S(k,k+1)*H(k+1,k); // s^H e_0 eta_21
    F77NAME(cgeru)(n-k,n-k,complex_float_one_,S.addr(k,k),1,
      H.addr(k,k),n,M->addr(k,k),n);
  }
  return M;
}

template<> Matrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *M=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<m;k++) {
        F77NAME(caxpy)(min(k+2,m),L(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(min(j+2,m),addr(0,j),1,M->addr(0,j),1);
      for (int k=j+1;k<m;k++) {
        F77NAME(caxpy)(min(k+2,m),L(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<float,complex<float> >* operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const UpperHessenbergMatrix<float,complex<float> > &H) {
// compute by columns: note that
//   [ L_11      ] [ h ] = [ L_11 h ]
//   [ L_21 L_22 ] [   ] = [ L_21 h ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,H.size(0));
  Matrix<float,complex<float> > *M=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n,float_zero_);
  char diag=(dynamic_cast<const UnitLowerTrapezoidalMatrix<float,
    complex<float> >*>(&L)==0 ? 'N' : 'U');
  for (int j=0;j<n;j++) {
    int k=min(j+2,n);
    F77NAME(ccopy)(k,H.addr(0,j),1,M->addr(0,j),1);
    F77NAME(ctrmv)('L','N',diag,k,L.addr(),m,M->addr(0,j),1); // L_11 h
    if (k<m) {
      F77NAME(cgemv)('N',m-k,k,complex_float_one_,L.addr(k,0),m,
        H.addr(0,j),1,complex_float_zero_,M->addr(k,j),1); // L_21 h
    }
  }
  return M;
}

template<> Matrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
// compute by columns: note that
//   H [ U_1 , U_2 ] = [ H U_1 , H U_2 ]
// and that
//   [      H_11     H_12 ] [ u ] = [     H_11      u ]
//   [ e_0 eta e_k^T H_22 ] [   ] = [ e_0 eta e_k^T u ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *M=
    OPERATOR_NEW Matrix<float,complex<float> >(m,n,float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0);
  if (U_non_unit) { // H U_1
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,m-1);k++) {
        F77NAME(caxpy)(min(k+2,m),U(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,m);k++) {
        F77NAME(caxpy)(min(k+2,m),U(k,j),addr(0,k),1,M->addr(0,j),1);
      }
      if (j<m) {
        F77NAME(caxpy)(min(j+2,m),float_one_,addr(0,j),1,M->addr(0,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<float,complex<float> >* operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const UpperHessenbergMatrix<float,complex<float> > &H) {
// compute by columns: note that
//   [ U_1 , U_2 ] [     H_11          , H_12 ]
//                 [ e_0 eta e_{m-1}^T , H_22 ]
//   = [ U_1 H_11 + U_2 e_0 eta e_{m-1}^T , U_1 H_12 + U_2 H_22 ]
// and that
//   [ U_11 U_12 ] [ h ] = [ U_11 h ]
//   [      U_22 ] [   ] = [        ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,H.size(0));
  char diag=(dynamic_cast<const UnitUpperTrapezoidalMatrix<float,
    complex<float> >*>(&U)==0 ? 'N' : 'U');
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<m;j++) { // U_1 H_11
    int k=min(j+2,m);
    F77NAME(ccopy)(k,H.addr(0,j),1,M->addr(0,j),1);
    F77NAME(ctrmv)('U','N',diag,k,U.addr(),m,M->addr(0,j),1); // U_11 h
  }
  if (m<n) {
    F77NAME(caxpy)(m,H(m,m-1),U.addr(0,m),1,M->addr(0,m-1),1);
      // U_2 e_0 eta e_{m-1}^T
    F77NAME(clacpy)('A',m,n-m,H.addr(0,m),m,M->addr(0,m),m);
    F77NAME(ctrmm)('L','U','N',diag,m,n-m,complex_float_one_,U.addr(),m,
      M->addr(0,m),m); // U_1 H_12
    F77NAME(cgemm)('N','N',m,n-m,n-m,complex_float_one_,U.addr(0,m),m,
      H.addr(m,m),n,complex_float_one_,M->addr(0,m),m); // U_2 H_22
  }
  return M;
}

template<> SquareMatrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const SquareMatrix<float,complex<float> > &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<n;k++) {
      F77NAME(caxpy)(min(k+2,n),S(k,j),addr(0,k),1,M->addr(0,j),1);
    }
  }
  return M;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SquareMatrix<float,complex<float> > &S,
const UpperHessenbergMatrix<float,complex<float> > &H) {
  int n=H.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *M=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    F77NAME(cgemv)('N',n,min(j+2,n),complex_float_one_,S.addr(),n,
      H.addr(0,j),1,complex_float_zero_,M->addr(0,j),1);
  }
  return M;
}

template<> Matrix<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,complex<float> > *R=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<m;k++) {
      F77NAME(caxpy)(min(m,k+2),M(k,j),addr(0,k),1,R->addr(0,j),1);
    }
  }
  return R;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const UpperHessenbergMatrix<float,complex<float> > &H) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,H.size(0));
  Matrix<float,complex<float> > *R=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<=min(n-1,j+1);k++) {
      F77NAME(caxpy)(m,H(k,j),M.addr(0,k),1,R->addr(0,j),1);
    }
  }
  return R;
}

template<> Vector<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<float,complex<float> > *r=
    OPERATOR_NEW Vector<float,complex<float> >(n,complex_float_zero_);
  for (int k=0;k<n;k++) {
    F77NAME(caxpy)(min(n,k+2),v[k],addr(0,k),1,r->addr(),1);
  }
  return r;
}

template<> void UpperHessenbergMatrix<float,complex<float> >::uhmv(
complex<float> alpha,const Vector<float,complex<float> > &x,
complex<float> beta,Vector<float,complex<float> > &b,char trans)
const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  if (abs(beta)==float_zero_) b=complex_float_zero_;
  else b*=beta;
  if (trans=='N' || trans=='n') { // b=H*x*alpha+b*beta
    complex<float> *xj=x.addr();
    for (int j=0;j<n;j++,xj++) {
      F77NAME(caxpy)(min(n,j+2),(*xj)*alpha,addr(0,j),1,b.addr(),1); 
    }
  } else if (trans=='T' || trans=='t') { // b=H^T*x*alpha+b*beta
    complex<float> *bj=b.addr();
    for (int j=0;j<n;j++,bj++) {
      *bj+=F77NAME(cdotu)(min(n,j+2),addr(0,j),1,x.addr(),1)*alpha;
    }
  } else { // b=H^H*x*alpha+b*beta
    complex<float> *bj=b.addr();
    for (int j=0;j<n;j++,bj++) {
      *bj+=F77NAME(cdotc)(min(n,j+2),addr(0,j),1,x.addr(),1)*alpha;
    }
  }
}

template<> void UpperHessenbergMatrix<float,complex<float> >::uhmm(
complex<float> alpha,const Matrix<float,complex<float> > &X,
complex<float> beta,Matrix<float,complex<float> > &B,char side,
char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (abs(beta)==float_zero_) B=complex_float_zero_;
  else B*=beta;
  if (side=='L' || side=='l') {
    CHECK_SAME(m,size(0));
    if (trans=='N' || trans=='n') { // B=H*X*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int k=0;k<m;k++) {
          F77NAME(caxpy)(min(m,k+2),X(k,j)*alpha,addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else if (trans=='T' || trans=='t') { // B=H^T*X*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int i=0;i<m;i++) {
          B(i,j)+=F77NAME(cdotu)(min(m,i+2),addr(0,i),1,X.addr(0,j),1)
                *alpha;
        }
      }
    } else { // B=H^H*X*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int i=0;i<m;i++) {
          B(i,j)+=F77NAME(cdotc)(min(m,i+2),addr(0,i),1,X.addr(0,j),1)
                *alpha;
        }
      }
    }
  } else {
    CHECK_SAME(n,size(0));
    if (trans=='N' || trans=='n') { // B=X*H*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int k=0;k<=min(n-1,j+1);k++) {
          F77NAME(caxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else if (trans=='T' || trans=='t') { // // B=X*H^T*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<n;k++) {
          F77NAME(caxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else { // B=X*H^H*alpha+B*beta
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<n;k++) {
          F77NAME(caxpy)(m,conj((*this)(j,k))*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

template<> float
UpperHessenbergMatrix<float,complex<float> >::normFrobenius() const {
  int n=size(0);
  float *work=0;
  return F77NAME(clanhs)('F',n,addr(),n,work);
}

template<> float
UpperHessenbergMatrix<float,complex<float> >::normInfinity() const {
  int n=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,n);
  float val=F77NAME(clanhs)('I',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> float
UpperHessenbergMatrix<float,complex<float> >::normMaxEntry() const {
  int n=size(0);
  float *work=0;
  return F77NAME(clanhs)('M',n,addr(),n,work);
}

template<> float
UpperHessenbergMatrix<float,complex<float> >::normOne() const {
  int n=size(0);
  float *work=0;
  return F77NAME(clanhs)('O',n,addr(),n,work);
}

template<> Vector<float,complex<float> >*
UpperHessenbergMatrix<float,complex<float> >::eigenvalues(
SquareMatrix<float,complex<float> > *&V,
SquareMatrix<float,complex<float> > *&U) const {
  int n=size(0);
  if (V!=0) CHECK_SAME(n,V->size(0));
  if (U!=0) CHECK_SAME(n,U->size(0));
  UpperHessenbergMatrix<float,complex<float> > *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,complex<float> >(*this);
  char job=(V==0 && U==0 ? 'E' : 'S');
  char compz=(V==0 && U==0 ? 'N' : 'I');
  Vector<float,complex<float> > *lambda =
    OPERATOR_NEW Vector<float,complex<float> >(n);
  OrthogonalMatrix<float,complex<float> > *Z=(V==0 && U==0 ? 0 :
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(n,n));
  complex<float> *za=(Z==0 ? 0 : Z->addr());
  complex<float> ww(complex_float_undefined_);
  int lwork=-1;
  int info;
  F77NAME(chseqr)(job,compz,n,1,n,H->addr(),n,lambda->addr(),za,n,&ww,
    lwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(ww.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(chseqr)(job,compz,n,1,n,H->addr(),n,lambda->addr(),za,n,work,
    lwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;

  if (Z!=0) {
    char side=(V==0 ? 'R' : (U==0 ? 'L' : 'B') );
    bool *select=0;
    complex<float> *vla=0;
//  SquareMatrix<float,complex<float> > *Vl=0;
    if (V!=0) {
//    Vl=OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
      V->copy(*Z);
      vla=V->addr();
    }
    complex<float> *vra=0;
//  SquareMatrix<float,complex<float> > *Vr=0;
    if (U!=0) {
//    Vr=OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
      U->copy(*Z);
      vra=U->addr();
    }
    work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
    int mout=-1;
    float *rwork=OPERATOR_NEW_BRACKET(float,n);
    F77NAME(ctrevc)(side,'B',select,n,H->addr(),n,vla,n,vra,n,n,mout,work,
      rwork,info);
    CHECK_TEST(info==0);
    if (V!=0) {
      for (int j=0;j<n;j++) {
        float val=F77NAME(scnrm2)(n,V->addr(0,j),1);
        F77NAME(csscal)(n,1./val,V->addr(0,j),1);
      }
    }
    if (U!=0) {
      for (int j=0;j<n;j++) {
        float val=F77NAME(scnrm2)(n,U->addr(0,j),1);
        F77NAME(csscal)(n,1./val,U->addr(0,j),1);
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

template<> int UpperHessenbergMatrix<float,complex<float> >::factor(
int *ipiv) {
  int n=size(0);
  for (int k=0;k<n;k++) ipiv[k]=k;
  int info=0;
  for (int k=0;k<n-1;k++) {
    if (abs((*this)(k,k))<abs((*this)(k+1,k)) ) {
      ipiv[k]=k+1;
      F77NAME(cswap)(n-k,addr(k,k),n,addr(k+1,k),n);
    }
    if (abs((*this)(k,k))>float_zero_) {
      (*this)(k+1,k)/=(*this)(k,k);
      F77NAME(caxpy)(n-k-1,-(*this)(k+1,k),addr(k,k+1),n,
        addr(k+1,k+1),n);
    } else {
      info=k+1;
      break;
    }
  }
  return info;
}

template<> void UpperHessenbergMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  int n=size(0);
  CHECK_SAME(n,b.size());
  CHECK_SAME(n,x.size());
  x.copy(b);
  UpperHessenbergMatrix<float,complex<float> > *HF=
    OPERATOR_NEW UpperHessenbergMatrix<float,complex<float> >(*this);
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
      x[j]=(x[j]-F77NAME(cdotu)(j,HF->addr(0,j),1,x.addr(),1))
          /(*HF)(j,j);
    }
    for (int i=n-2;i>=0;i--) { // back-solve L^T Q z = y
      int ip=ipiv[i];
      complex<float> temp=x[i]-(*HF)(i+1,i)*x[i+1];
      x[i]=x[ip];
      x[ip]=temp;
    }
    if (trans=='C' || trans=='c') {
      for (int j=0;j<n;j++) x[j]=conj(x[j]); // x = conj(y)
    }
  } else { // H = Q^T L U and H x = b ==> Q^T L U x = b
    for (int i=0;i<n-1;i++) { // forward-solve Q^T L  = b
      int ip=ipiv[i];
      complex<float> temp=x[2*i-ip+1]-(*HF)(i+1,i)*x[ip];
      x[i]=x[ip];
      x[i+1]=temp;
    }
    for (int j=n-1;j>0;j--) { // back-solve U x = y
      x[j]/=(*HF)(j,j);
      F77NAME(caxpy)(j,-x[j],HF->addr(0,j),1,x.addr(),1);
    }
    x[0]/=(*HF)(0,0);
  }
  delete ipiv; ipiv=0;
  delete HF; HF=0;
}

template<> void UpperHessenbergMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  int n=size(0);
  UpperHessenbergMatrix<float,complex<float> > *HF=
    OPERATOR_NEW UpperHessenbergMatrix<float,complex<float> >(*this);
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
        complex<float> *X0k=X.addr(0,k);
        complex<float> *Xjk=X0k;
        if (trans=='C' || trans=='c') {
          for (int j=0;j<n;j++,Xjk+=nrhs) *Xjk=conj(*Xjk);
        }
        *X0k/=(*HF)(0,0); // forward-solve U^T Y = conj(B)
        Xjk=X0k+nrhs;
        for (int j=1;j<n;j++,Xjk+=nrhs) {
          *Xjk=(*Xjk-F77NAME(cdotu)(j,HF->addr(0,j),1,X0k,1))/(*HF)(j,j);
        }
        for (int i=n-2;i>=0;i--) { // back-solve L^T Q Z = Y
          int ip=ipiv[i];
          complex<float> temp=X(i,k)-(*HF)(i+1,i)*X(i+1,k);
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
          complex<float> temp=X(2*i-ip+1,k)-(*HF)(i+1,i)*X(ip,k);
          X(i,k)=X(ip,k);
          X(i+1,k)=temp;
        }
        for (int j=n-1;j>0;j--) { // back-solve U X = Y
          X(j,k)/=(*HF)(j,j);
          F77NAME(caxpy)(j,-X(j,k),HF->addr(0,j),1,X.addr(0,k),1);
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
        complex<float> *Xk0=X.addr(k,0);
        complex<float> *Xki=Xk0;
        if (trans=='C' || trans=='c') {
          for (int i=0;i<n;i++,Xki+=nrhs) *Xki=conj(*Xki);
        }
        for (int i=0;i<n-1;i++) { // forward-solve Q^T L Y = B^H
          int ip=ipiv[i];
          complex<float> temp=X(k,2*i-ip+1)-(*HF)(i+1,i)*X(k,ip);
          X(k,i)=X(k,ip);
          X(k,i+1)=temp;
        }
        complex<float> *Xkj=X.addr(k,n-1);
        for (int j=n-1;j>0;j--,Xkj-=nrhs) { // back-solve U Z = Y
          *Xkj/=(*HF)(j,j);
          F77NAME(caxpy)(j,-(*Xkj),HF->addr(0,j),1,Xk0,nrhs);
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
        complex<float> *Xk0=X.addr(k,0);
        complex<float> *Xkj=Xk0;
        *Xkj/=(*HF)(0,0); // forward-solve U^T Y = B^T
        for (int j=1;j<n;j++,Xkj+=nrhs) {
          *Xkj=(*Xkj-F77NAME(cdotu)(j,HF->addr(0,j),1,Xk0,nrhs))
              /(*HF)(j,j);
        }
        for (int i=n-2;i>=0;i--) { // back-solve L^T Q X^T = Y
          int ip=ipiv[i];
          complex<float> temp=X(k,i)-(*HF)(i+1,i)*X(k,i+1);
          X(k,i)=X(k,ip);
          X(k,ip)=temp;
        }
      }
    }
  }
  delete ipiv; ipiv=0;
  delete HF; HF=0;
}

template class UpperHessenbergMatrix<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const complex<float>
  BandMatrix<float,complex<float> >::outofbounds_=complex_float_zero_;
template<> complex<float> BandMatrix<float,complex<float> >::safety_ =
  complex_float_zero_;

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::makeMatrix() const {
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      M->addr(ibeg,j),1);
  }
  return M;
}

template<> BandMatrix<float,complex<float> >&
BandMatrix<float,complex<float> >::operator+=(
const BandMatrix<float,complex<float> > &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  CHECK_SAME(nsup,B.nsup);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,complex_float_one_,
      B.addr(ibeg,j),1,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<float,complex<float> >&
BandMatrix<float,complex<float> >::operator-=(
const BandMatrix<float,complex<float> > &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  CHECK_SAME(nsup,B.nsup);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<float,complex<float> >&
BandMatrix<float,complex<float> >::operator*=(complex<float> d) {
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(cscal)(min(dim-1,j+nsub)-ibeg+1,d,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<float,complex<float> >&
BandMatrix<float,complex<float> >::operator/=(complex<float> d) {
  CHECK_TEST(abs(d)>float_zero_);
  complex<float> dinv=complex_float_one_/d;
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(cscal)(min(dim-1,j+nsub)-ibeg+1,dinv,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator+(
const BandMatrix<float,complex<float> > &B) const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,B.nsub),
    max(nsup,B.nsup),complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    ibeg=max(0,j-B.nsup);
    F77NAME(caxpy)(min(dim-1,j+B.nsub)-ibeg+1,complex_float_one_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* 
BandMatrix<float,complex<float> >::operator+(
const UpperHessenbergMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    F77NAME(caxpy)(min(dim,j+2),complex_float_one_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator+(
const DiagonalMatrix<float,complex<float> > &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<float,complex<float> > *S=
    OPERATOR_NEW BandMatrix<float,complex<float> >(*this);
  F77NAME(caxpy)(dim,complex_float_one_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator+(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,1),max(nsup,1),
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(caxpy)(dim-1,complex_float_one_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),nt);
  complex<float> *Sjj=S->addr(0,0);
  const float *Tjj=T.diagonalAddr(0);
  for (int j=0;j<dim;j++,Sjj+=S->nt,Tjj++) *Sjj+=*Tjj;
  complex<float> *Sjjp1=S->addr(0,1);
  const complex<float> *Tjm1j=T.lowerDiagonalAddr(0);
  for (int j=0;j<dim-1;j++,Sjjp1+=S->nt,Tjm1j++) *Sjjp1+=conj(*Tjm1j);
  return S;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator+(
const TridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,1),max(nsup,1),
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(caxpy)(dim,complex_float_one_,T.diagonalAddr(),1,
    S->addr(0,0),S->nt);
  F77NAME(caxpy)(dim-1,complex_float_one_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S->nt);
  F77NAME(caxpy)(dim-1,complex_float_one_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S->nt);
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator+(
const SymmetricMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<dim;j++) {
    F77NAME(caxpy)(dim-j,complex_float_one_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(caxpy)(dim-j-1,complex_float_one_,H.addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(caxpy)(j+1,complex_float_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(caxpy)(dim-j,complex_float_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      (*S)(j,j)+=complex_float_one_;
      if (j+1<dim) {
        F77NAME(caxpy)(dim-j-1,complex_float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  S->copy(M);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,complex_float_one_,
      addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const BandMatrix<float,complex<float> > &B) const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,B.nsub),
    max(nsup,B.nsup),complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    ibeg=max(0,j-B.nsup);
    F77NAME(caxpy)(min(B.dim-1,j+B.nsub)-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const UpperHessenbergMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    F77NAME(caxpy)(min(dim,j+2),complex_float_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const UpperHessenbergMatrix<float,complex<float> > &H,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(n,j+2),H.addr(0,j),1,S->addr(0,j),1);
    int ibeg=max(0,j-B.supDiags());
    F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const DiagonalMatrix<float,complex<float> > &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<float,complex<float> > *S=
    OPERATOR_NEW BandMatrix<float,complex<float> >(*this);
  F77NAME(caxpy)(dim,complex_float_mone_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> BandMatrix<float,complex<float> >* operator-(
const DiagonalMatrix<float,complex<float> > &D,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,B.subDiags(),B.supDiags(),
    complex_float_zero_);
  F77NAME(ccopy)(n,D.addr(),1,S->addr(0,0),S->bands());
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,1),max(nsup,1),
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(caxpy)(dim-1,complex_float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  complex<float> *Sjj=S->addr(0,0);
  const float *Tjj=T.diagonalAddr(0);
  for (int j=0;j<dim;j++,Sjj+=S->nt,Tjj++) *Sjj-=*Tjj;
  complex<float> *Sjjp1=S->addr(0,1);
  const complex<float> *Tjp1j=T.lowerDiagonalAddr(0);
  for (int j=0;j<dim-1;j++,Sjjp1+=S->nt,Tjp1j++) *Sjjp1+=conj(*Tjp1j);
  return S;
}

template<> BandMatrix<float,complex<float> >* operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,max(B.subDiags(),1),
    max(B.supDiags(),1),complex_float_zero_);
  int nb=S->bands();
  F77NAME(ccopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(1,0),nb);
  complex<float> *Sjj=S->addr(0,0);
  const float *Tjj=T.diagonalAddr(0);
  for (int j=0;j<n;j++,Sjj+=nb,Tjj++) *Sjj=*Tjj;
  complex<float> *Sjjp1=S->addr(0,1);
  const complex<float> *Tjp1j=T.lowerDiagonalAddr(0);
  for (int j=0;j<n-1;j++,Sjjp1+=nb,Tjp1j++) *Sjjp1=conj(*Tjp1j);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const TridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,1),max(nsup,1),
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(caxpy)(dim,complex_float_mone_,T.diagonalAddr(),1,
    S->addr(0,0),nt);
  F77NAME(caxpy)(dim-1,complex_float_mone_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),nt);
  F77NAME(caxpy)(dim-1,complex_float_mone_,T.upperDiagonalAddr(),1,
    S->addr(0,1),nt);
  return S;
}

template<> BandMatrix<float,complex<float> >* operator-(
const TridiagonalMatrix<float,complex<float> > &T,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,max(B.subDiags(),1),
    max(B.supDiags(),1),complex_float_zero_);
  F77NAME(ccopy)(n,T.diagonalAddr(),1,S->addr(0,0),S->bands());
  F77NAME(ccopy)(n-1,T.lowerDiagonalAddr(),1,S->addr(1,0),S->bands());
  F77NAME(ccopy)(n-1,T.upperDiagonalAddr(),1,S->addr(0,1),S->bands());
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const SymmetricMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<dim;j++) {
    F77NAME(caxpy)(dim-j,complex_float_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim) {
      complex<float> *Sji=S->addr(j,j+1);
      const complex<float> *Hij=H.addr(j+1,j);
      for (int i=j+1;i<dim;i++,Sji+=dim,Hij++) *Sji-=conj(*Hij);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const SymmetricMatrix<float,complex<float> > &H,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,H.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      complex<float> *Sji=S->addr(j,j+1);
      const complex<float> *Hij=H.addr(j+1,j);
      for (int i=j+1;i<n;i++,Sji+=n,Hij++) *Sji=conj(*Hij);
    }
  }
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(caxpy)(j+1,complex_float_mone_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_mone_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)-=complex_float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(j+1,U.addr(0,j),1,S->addr(0,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-B.supDiags());
      if (j>0) F77NAME(ccopy)(j,U.addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)=complex_float_one_;
      F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(caxpy)(dim-j,complex_float_mone_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      (*S)(j,j)-=complex_float_one_;
      if (j+1<dim) {
        F77NAME(caxpy)(dim-j-1,complex_float_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=complex_float_one_;
      if (j+1<n) F77NAME(ccopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  S->operator-=(M);
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  S->copy(M);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(caxpy)(min(n-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const BandMatrix<float,complex<float> > &B) const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,nsub+B.nsub,nsup+B.nsup,
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=max(0,j-B.nsup);k<=min(dim-1,j+B.nsub);k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,B(k,j),
        addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const UpperHessenbergMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=0;k<=min(dim-1,j+1);k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,H(k,j),
        addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const UpperHessenbergMatrix<float,complex<float> > &H,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(caxpy)(min(n,k+2),B(k,j),H.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const DiagonalMatrix<float,complex<float> > &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<float,complex<float> > *P=
    OPERATOR_NEW BandMatrix<float,complex<float> >(*this);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(cscal)(min(dim-1,j+nsub)-ibeg+1,D[j],
      P->addr(ibeg,j),1);
  }
  return P;
}

template<> BandMatrix<float,complex<float> >* operator*(
const DiagonalMatrix<float,complex<float> > &D,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,B.subDiags(),B.supDiags());
  P->copy(B);
  int stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    int jbeg=max(0,i-B.subDiags());
    int jend=min(n-1,i+B.supDiags());
    F77NAME(cscal)(min(n-1,i+B.supDiags())-jbeg+1,D[i],
      P->addr(i,jbeg),stride);
  }
  return P;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,nsub+1,nsup+1,
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    if (j>0) {
      int ibeg=max(0,j-1-nsup);
      F77NAME(caxpy)(min(dim-1,nsub+j-1)-ibeg+1,
        T.upperDiagonalValue(j-1),addr(ibeg,j-1),1,P->addr(ibeg,j),1);
    }
    int ibeg=max(0,j-nsup);
    F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,T.diagonalValue(j),
      addr(ibeg,j),1,P->addr(ibeg,j),1);
    if (j<dim-1) {
      int ibeg=max(0,j+1-nsup);
      F77NAME(caxpy)(min(dim-1,j+1+nsub)-ibeg+1,
        T.lowerDiagonalValue(j),addr(ibeg,j+1),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> BandMatrix<float,complex<float> >* operator*(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,B.subDiags()+1,B.supDiags()+1,
    complex_float_zero_);
  int B_stride=B.bands()-1;
  int P_stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    if (i>0) {
      int jbeg=max(0,i-1-B.subDiags());
      F77NAME(caxpy)(min(n-1,i-1+B.supDiags())-jbeg+1,
        T.lowerDiagonalValue(i-1),B.addr(i-1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
    int jbeg=max(0,i-B.subDiags());
    F77NAME(caxpy)(min(n-1,i+B.supDiags())-jbeg+1,T.diagonalValue(i),
      B.addr(i,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    if (i<n-1) {
      int jbeg=max(0,i+1-B.subDiags());
      F77NAME(caxpy)(min(n-1,i+1+B.supDiags())-jbeg+1,
        T.upperDiagonalValue(i),B.addr(i+1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
  }
  return P;
}

template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const TridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,nsub+1,nsup+1,
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    if (j>0) {
      int ibeg=max(0,j-1-nsup);
      F77NAME(caxpy)(min(dim-1,nsub+j-1)-ibeg+1,T(j-1,j),
        addr(ibeg,j-1),1,P->addr(ibeg,j),1);
    }
    int ibeg=max(0,j-nsup);
    F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,T(j,j),addr(ibeg,j),1,
      P->addr(ibeg,j),1);
    if (j<dim-1) {
      int ibeg=max(0,j+1-nsup);
      F77NAME(caxpy)(min(dim-1,j+1+nsub)-ibeg+1,T(j+1,j),
        addr(ibeg,j+1),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> BandMatrix<float,complex<float> >* operator*(
const TridiagonalMatrix<float,complex<float> > &T,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,B.subDiags()+1,B.supDiags()+1,
    complex_float_zero_);
  int B_stride=B.bands()-1;
  int P_stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    if (i>0) {
      int jbeg=max(0,i-1-B.subDiags());
      F77NAME(caxpy)(min(n-1,i-1+B.supDiags())-jbeg+1,T(i,i-1),
        B.addr(i-1,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    }
    int jbeg=max(0,i-B.subDiags());
    F77NAME(caxpy)(min(n-1,i+B.supDiags())-jbeg+1,T(i,i),
      B.addr(i,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    if (i<n-1) {
      int jbeg=max(0,i+1-B.subDiags());
      F77NAME(caxpy)(min(n-1,i+1+B.supDiags())-jbeg+1,
        T(i,i+1),B.addr(i+1,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    }
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const SymmetricMatrix<float,complex<float> > &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=0;k<dim;k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,S(k,j),addr(ibeg,k),1,
        P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SymmetricMatrix<float,complex<float> > &S,
const BandMatrix<float,complex<float> > &B) {
//compute by bordering: note that
//  [ sigma s^H ] [ beta c^T ]
//  [   s    S  ] [   b   B  ]
//  = [ sigma beta + s^H b , sigma c^T + s^H B ]
//  = [   s   beta +  S  b ,   s  c^T +   S  B ]
  int n=B.size(0),nsub=B.subDiags(),nsup=B.supDiags();
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  complex<float> *t=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int k=n-1;k>=0;k--) {
    int iend=min(n-1,k+nsub);
    if (k<n-1) {
      // s^H B = ( B^H s )^H
      F77NAME(cgbmv)('C',n-k-1,n-k-1,nsub,nsup,
        complex_float_one_,B.addr(k+1-nsup,k+1),B.bands(),
        S.addr(k+1,k),1,complex_float_zero_,t+k+1,1);
      for (int j=k+1;j<n;j++) (*P)(k,j)=conj(t[j]);
      // S b: note that
      // [ S_11 , S_21^H ] [ b ] = [ S_11 b ]
      // [ S_21 ,  S_22  ] [ 0 ] = [ S_21 b ]
      if (nsub>0) { // S_11 b
        F77NAME(chemv)('L',min(n-k-1,nsub),complex_float_one_,
          S.addr(k+1,k+1),n,B.addr(k+1,k),1,complex_float_zero_,
          P->addr(k+1,k),1);
        if (nsub<n-k-1) { // S_21 b
          F77NAME(cgemv)('N',n-k-1-nsub,nsub,complex_float_one_,
            S.addr(k+1+nsub,k+1),n,B.addr(k+1,k),1,complex_float_zero_,
            P->addr(k+1+nsub,k),1);
        }
      }
    }
    if (iend>k) { // s^H b
      (*P)(k,k)=F77NAME(cdotc)(iend-k,S.addr(k+1,k),1,B.addr(k+1,k),1);
    }
    int jend=min(n-1,k+nsup);
    F77NAME(cgeru)(n-k,jend-k+1,float_one_,S.addr(k,k),1,
      B.addr(k,k),B.bands()-1,P->addr(k,k),n);
  }
  delete [] t; t=0;
  return P;
}

template<> Matrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  CHECK_SAME(dim,U.size(0));
  int n=U.size(1);
  Matrix<float,complex<float> > *M=
    OPERATOR_NEW Matrix<float,complex<float> >(dim,U.size(1),
    complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,dim-1);k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,dim);k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
      if (j<dim) {
        int ibeg=max(0,j-nsup);
        F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,float_one_,
          addr(ibeg,j),1,M->addr(ibeg,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<float,complex<float> >* operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const BandMatrix<float,complex<float> > &B) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(caxpy)(min(k+1,m),B(k,j),U.addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(caxpy)(min(k,m),B(k,j),U.addr(0,k),1,M->addr(0,j),1);
        if (k<m) (*M)(k,j)+=B(k,j);
      }
    }
  }
  return M;
}

template<> Matrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  CHECK_SAME(dim,L.size(0));
  int n=L.size(1);
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(dim,L.size(1),complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<dim;k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,complex_float_one_,
        addr(ibeg,j),1,M->addr(ibeg,j),1);
      for (int k=j+1;k<dim;k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<float,complex<float> >* operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const BandMatrix<float,complex<float> > &B) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,complex<float> > *M=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(caxpy)(m-k,B(k,j),L.addr(k,k),1,M->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        (*M)(k,j)+=B(k,j);
        if (k+1<m) {
          F77NAME(caxpy)(m-k-1,B(k,j),L.addr(k+1,k),1,M->addr(k+1,j),1);
        }
      }
    }
  }
  return M;
}

template<> SquareMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const SquareMatrix<float,complex<float> > &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(cgbmv)('N',dim,dim,nsub,nsup,complex_float_one_,addr(),nt,
      S.addr(0,j),1,complex_float_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SquareMatrix<float,complex<float> > &S,
const BandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(caxpy)(n,B(k,j),S.addr(0,j),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  CHECK_SAME(dim,M.size(0));
  int n=M.size(1);
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(dim,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(cgbmv)('N',dim,dim,nsub,nsup,complex_float_one_,addr(),nt,
      M.addr(0,j),1,complex_float_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const BandMatrix<float,complex<float> > &B) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(caxpy)(m,B(k,j),M.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Vector<float,complex<float> >*
BandMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  CHECK_SAME(dim,v.size());
  Vector<float,complex<float> > *w=
    OPERATOR_NEW Vector<float,complex<float> >(dim,float_zero_);
  F77NAME(cgbmv)('N',dim,dim,nsub,nsup,complex_float_one_,addr(),nt,
    v.addr(),1,complex_float_zero_,w->addr(),1);
  return w;
}

template<> void BandMatrix<float,complex<float> >::gbmv(
complex<float> alpha,const Vector<float,complex<float> > &x,
complex<float> beta,Vector<float,complex<float> > &b,char trans)
const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  F77NAME(cgbmv)(trans,dim,dim,nsub,nsup,alpha,addr(),nt,x.addr(),1,
    beta,b.addr(),1);
}

template<> void BandMatrix<float,complex<float> >::gbmm(
complex<float> alpha,const Matrix<float,complex<float> > &X,
complex<float> beta,Matrix<float,complex<float> > &B,char side,
char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') { // B=op(A)*X*alpha+B*beta
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      F77NAME(cgbmv)(trans,dim,dim,nsub,nsup,alpha,addr(),nt,
        X.addr(0,j),1,beta,B.addr(0,j),1);
    }
  } else { // B=X*op(A)*alpha+B*beta
    CHECK_SAME(n,dim);
    if (abs(beta)==float_zero_) B=float_zero_;
    else B*=beta;
    if (trans=='N' || trans=='n') {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsup);k<=min(dim-1,j+nsub);k++) {
          F77NAME(caxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else if (trans=='C' || trans=='c') {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsub);k<=min(dim-1,j+nsup);k++) {
          F77NAME(caxpy)(m,conj((*this)(j,k))*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsub);k<=min(dim-1,j+nsup);k++) {
          F77NAME(caxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

/*
template<> BandMatrix<float,complex<float> >*
BandMatrix<float,complex<float> >::transpose() const {
  BandMatrix<float,complex<float> > *X=
    OPERATOR_NEW BandMatrix<float,complex<float> >(dim,nsup,nsub);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsub);
    F77NAME(dcopy)(min(dim-1,j+nsup)-i+1,addr(j,i),nsub+nsup,
      X->addr(i,j),1);
  }
  return X;
}
*/

template<> float BandMatrix<float,complex<float> >::equilibrate(
Vector<float,float> &r,Vector<float,float> &c,float &rowcnd,
float &colcnd) const {
  CHECK_SAME(dim,r.size());
  CHECK_SAME(dim,c.size());
  float amax;
  int info;
  F77NAME(cgbequ)(dim,dim,nsub,nsup,addr(),nt,r.addr(),c.addr(),
    rowcnd,colcnd,amax,info);
  CHECK_TEST(info==0);
  return amax;
}

template<> float BandMatrix<float,complex<float> >::normFrobenius()
const {
  float *work=0;
  return F77NAME(clangb)('F',dim,nsub,nsup,addr(),nt,work);
}

template<> float BandMatrix<float,complex<float> >::normInfinity()
const {
  float *work=OPERATOR_NEW_BRACKET(float,dim);
  float val=F77NAME(clangb)('I',dim,nsub,nsup,addr(),nt,work);
  delete work;
  return val;
}

template<> float BandMatrix<float,complex<float> >::normMaxEntry()
const {
  float *work=0;
  return F77NAME(clangb)('M',dim,nsub,nsup,addr(),nt,work);
}

template<> float BandMatrix<float,complex<float> >::normOne() const {
  float *work=0;
  return F77NAME(clangb)('O',dim,nsub,nsup,addr(),nt,work);
}

template<> float
BandMatrix<float,complex<float> >::reciprocalConditionNumber(char norm)
const {
  BandMatrix<float,complex<float> > *BF=
    OPERATOR_NEW BandMatrix<float,complex<float> >(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(cgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);

  float anorm=(norm=='I' || norm=='i' ? normInfinity() : normOne());
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*dim);
  float *rwork=OPERATOR_NEW_BRACKET(float,dim);
  F77NAME(cgbcon)(norm,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,anorm,
    rcond,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete [] rwork; rwork=0;
  delete [] ipiv; ipiv=0;
  delete BF; BF=0;
  return rcond;
}

template<> void BandMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char trans) const {
  CHECK_SAME(dim,b.size())
  CHECK_SAME(dim,x.size())
  BandMatrix<float,complex<float> > *BF=
    OPERATOR_NEW BandMatrix<float,complex<float> >(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(cgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);

  x.copy(b);
  F77NAME(cgbtrs)(trans,dim,nsub,nsup,1,BF->addr(),BF->nt,ipiv,
    x.addr(),dim,info);
  CHECK_TEST(info==0);
  delete [] ipiv;
  delete BF; BF=0;
}

template<> void BandMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char trans) const {
  BandMatrix<float,complex<float> > *BF=
    OPERATOR_NEW BandMatrix<float,complex<float> >(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(ccopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(cgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(dim,B.size(0))
    CHECK_SAME(dim,X.size(0))
    CHECK_SAME(nrhs,X.size(1))
    X.copy(B);
    F77NAME(cgbtrs)(trans,dim,nsub,nsup,nrhs,BF->addr(),BF->nt,ipiv,
      X.addr(),dim,info);
    CHECK_TEST(info==0);
  } else { // cgbtrs
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
          complex<float> *Xij=X.addr(0,j);
          for (int i=0;i<nrhs;i++,Xij++) *Xij=conj(*Xij);
        }
      }
      if (nsub>0) { // solve Q^T L Y^T = B^T: cgbtrs
        for (int j=0;j<dim-1;j++) {
          int lm=min(nsub,dim-j-1);
          int jp=ipiv[j]-1;
          if (jp!=j) F77NAME(cswap)(nrhs,X.addr(0,jp),1,X.addr(0,j),1);
          F77NAME(cgeru)(nrhs,lm,float_mone_,BF->addr(j+1,j),1,
            X.addr(0,j),1,X.addr(0,j+1),1);
        }
      }
      for (int i=0;i<nrhs;i++) { // solve U X^T = Y^T: ctbsv loop 20
        for (int j=dim-1;j>=0;j--) {
          if (abs(X(i,j))>float_zero_) {
            X(i,j)/=(*BF)(j,j);
            int ii=max(0,j-nsub-nsup);
            F77NAME(caxpy)(j-ii,-X(i,j),BF->addr(ii,j),1,
              X.addr(i,ii),nrhs);
          }
        }
      }
      if (trans=='C' || trans=='c') {
        for (int j=0;j<dim;j++) {
          complex<float> *Xij=X.addr(0,j);
          for (int i=0;i<nrhs;i++,Xij++) *Xij=conj(*Xij);
        }
      }
    } else {
    // A = Q^T L U and B = X A = X Q^T L U ==> U^T L^T Q X^T = B^T
      for (int i=0;i<nrhs;i++) { // solve U^T Y^T=B^T: ctbsv loop 100
        for (int j=0;j<dim;j++) {
          int ii=max(0,j-nsub-nsup);
          X(i,j)=(X(i,j)
            -F77NAME(cdotu)(j-i,BF->addr(i,j),1,X.addr(i,ii),1))
            /(*BF)(j,j);
        }
      }
      if (nsub>0) { // solve L^T Q X^T = Y^T: cgbtrs loop 40
        for (int j=dim-2;j>=0;j--) {
          int lm=min(nsub,dim-1-j);
          F77NAME(cgemv)('T',nrhs,lm,float_mone_,X.addr(0,j+1),nrhs,
            BF->addr(j+1,j),1,float_one_,X.addr(0,j),nrhs);
          int jp=ipiv[j]-1;
          if (jp!=j) {
            F77NAME(cswap)(nrhs,X.addr(0,jp),1,X.addr(0,j),1);
          }
        }
      }
    }
  }

  delete [] ipiv;
  delete BF; BF=0;
}

template class BandMatrix<float,complex<float> >;
//template void testBandMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const complex<float>
  SymmetricBandMatrix<float,complex<float> >::outofbounds_
  = complex_float_zero_;
template<> complex<float>
  SymmetricBandMatrix<float,complex<float> >::safety_ =
  complex_float_zero_;

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::makeMatrix() const {
  SquareMatrix<float,complex<float> > *M=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(dim-j,nsub+1),addr(j,j),1,M->addr(j,j),1);
    if (j+1<dim) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Mji=M->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Mji+=dim) {
        *Mji=conj(*Aij);
      }
    }
  }
  return M;
}

template<> void SymmetricBandMatrix<float,complex<float> >::fillWith(
complex<float> d) {
  complex<float> *colj=addr();
  for (int j=0;j<dim;j++,colj+=nt) {
    colj[0]=d.real();
    for (int i=1;i<min(dim-j,nt);i++) colj[i]=d;
  }
}

template<> complex<float>
SymmetricBandMatrix<float,complex<float> >::operator()(int i,int j)
const {
  if (i-j>=0 && i-j<=nsub) return *AB->addr(i-j,j);
  if (j-i>0 && j-i<=nsub) return conj(*AB->addr(j-i,i));
  return outofbounds_;
}

template<>
SymmetricBandMatrix<float,complex<float> >::SymmetricBandMatrix(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) : nsub(1),
nt(2) {
  dim=T.size(0);
  AB=OPERATOR_NEW Matrix<float,complex<float> >(2,dim);
  complex<float> *dii=addr(0,0);
  float *Tii=T.diagonalAddr(0);
  for (int i=0;i<dim;i++,dii+=2,Tii++) *dii=*Tii;
  F77NAME(ccopy)(dim-1,T.lowerDiagonalAddr(0),1,addr(1,0),2);
  (*AB)(1,dim-1)=complex_float_undefined_;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::makeBandMatrix() const {
  BandMatrix<float,complex<float> > *B=
    OPERATOR_NEW BandMatrix<float,complex<float> >(dim,nsub,nsub);
  int stride=B->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,B->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Bji=B->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Bji+=stride) {
        *Bji=conj(*Aij);
      }
    }
  }
  return B;
}

template<> SymmetricMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::makeSymmetricMatrix()
const {
  SymmetricMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >&
SymmetricBandMatrix<float,complex<float> >::operator+=(
const SymmetricBandMatrix<float,complex<float> > &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  for (int j=0;j<dim;j++) {
    F77NAME(caxpy)(min(nt,dim-j),complex_float_one_,B.addr(j,j),1,
      addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<float,complex<float> >&
SymmetricBandMatrix<float,complex<float> >::operator-=(
const SymmetricBandMatrix<float,complex<float> > &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  for (int j=0;j<dim;j++) {
    F77NAME(caxpy)(min(nt,dim-j),complex_float_mone_,B.addr(j,j),1,
      addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<float,complex<float> >&
SymmetricBandMatrix<float,complex<float> >::operator*=(
complex<float> d) {
  for (int j=0;j<dim;j++) {
    F77NAME(cscal)(min(nt,dim-j),d,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<float,complex<float> >&
SymmetricBandMatrix<float,complex<float> >::operator/=(
complex<float> d) {
  CHECK_TEST(abs(d)>float_zero_);
  complex<float> dinv=complex_float_one_/d;
  for (int j=0;j<dim;j++) {
    F77NAME(cscal)(min(nt,dim-j),dinv,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const SymmetricBandMatrix<float,complex<float> > &B) const {
  CHECK_SAME(dim,B.dim);
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(dim,max(nsub,B.nsub),
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(caxpy)(min(B.nt,dim-j),complex_float_one_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const BandMatrix<float,complex<float> > &B) const {
  CHECK_SAME(dim,B.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,B.subDiags()),
    max(nsub,B.supDiags()),complex_float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=stride) {
        *Sji=conj(*Aij);
      }
    }
    int ibeg=max(0,j-B.supDiags());
    F77NAME(caxpy)(min(dim-1,j+B.subDiags())-ibeg+1,complex_float_one_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const UpperHessenbergMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
        *Sji=conj(*Aij);
      }
    }
    F77NAME(caxpy)(min(dim,j+2),complex_float_one_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const DiagonalMatrix<float,complex<float> > &D) const {
  CHECK_SAME(dim,D.size(0));
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(*this);
  F77NAME(caxpy)(dim,complex_float_one_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(dim,max(nsub,1),
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  F77NAME(caxpy)(dim-1,complex_float_one_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<dim;i++,Tii++,Sii+=S->nt) *Sii+=*Tii;
  return S;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const TridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,1),max(nsub,1),
    complex_float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=stride) {
        *Sji=conj(*Aij);
      }
    }
  }
  F77NAME(caxpy)(dim,complex_float_one_,T.diagonalAddr(),1,
    S->addr(0,0),S->bands());
  F77NAME(caxpy)(dim-1,complex_float_one_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S->bands());
  F77NAME(caxpy)(dim-1,complex_float_one_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S->bands());
  return S;
}

template<> SymmetricMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const SymmetricMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SymmetricMatrix<float,complex<float> > *S=makeSymmetricMatrix();
  for (int j=0;j<dim;j++) {
    F77NAME(caxpy)(dim-j,complex_float_one_,H.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<float> *Aij=addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      F77NAME(caxpy)(j+1,complex_float_one_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<float> *Aij=addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)+=complex_float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<float> *Aij=addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      F77NAME(caxpy)(dim-j,complex_float_one_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<float> *Aij=addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      (*S)(j,j)+=complex_float_one_;
      if (j+1<dim) {
        F77NAME(caxpy)(dim-j-1,complex_float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator+(
const Matrix<float,complex<float> > &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  S->copy(M);
  for (int j=0;j<dim;j++) {
    F77NAME(caxpy)(min(nt,dim-j),complex_float_one_,addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
        *Sji+=conj(*Aij);
      }
    }
  }
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const SymmetricBandMatrix<float,complex<float> > &B) const {
  CHECK_SAME(dim,B.dim);
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(dim,max(nsub,B.nsub),
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(caxpy)(min(B.nt,dim-j),complex_float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const BandMatrix<float,complex<float> > &B) const {
  CHECK_SAME(dim,B.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,B.subDiags()),
    max(nsub,B.supDiags()),complex_float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=stride) {
        *Sji=conj(*Aij);
      }
    }
    int ibeg=max(0,j-B.supDiags());
    F77NAME(caxpy)(min(dim-1,j+B.subDiags())-ibeg+1,complex_float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >* operator-(
const BandMatrix<float,complex<float> > &B,
const SymmetricBandMatrix<float,complex<float> > &H) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,max(H.subDiags(),B.subDiags()),
    max(H.subDiags(),B.supDiags()),complex_float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(ccopy)(min(n-1,j+B.subDiags())-ibeg+1,B.addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(H.bands(),n-j),complex_float_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<n && H.subDiags()>0) {
      const complex<float> *Hij=H.addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(n-1,j+H.subDiags());i++,Hij++,Sji+=stride) {
        *Sji-=conj(*Hij);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const UpperHessenbergMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
        *Sji=conj(*Aij);
      }
    }
    F77NAME(caxpy)(min(dim,j+2),complex_float_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const UpperHessenbergMatrix<float,complex<float> > &H,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(n,j+2),H.addr(0,j),1,S->addr(0,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<n && B.subDiags()>0) {
      const complex<float> *Bij=B.addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
        *Sji-=conj(*Bij);
      }
    }
  }
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const DiagonalMatrix<float,complex<float> > &D) const {
  CHECK_SAME(dim,D.size(0));
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(*this);
  F77NAME(caxpy)(dim,complex_float_mone_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >* operator-(
const DiagonalMatrix<float,complex<float> > &D,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(n,B.subDiags(),
    complex_float_zero_);
  int nt=B.bands();
  F77NAME(ccopy)(n,D.addr(),1,S->addr(0,0),nt);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(nt,n-j),complex_float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(dim,max(nsub,1),
    complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  F77NAME(caxpy)(dim-1,complex_float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<dim;i++,Tii++,Sii+=S->nt) *Sii-=*Tii;
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >* operator-(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(n,max(B.subDiags(),1),
    complex_float_zero_);
  F77NAME(ccopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(1,0),S->bands());
  const float *Tii=T.diagonalAddr(0);
  complex<float> *Sii=S->addr(0,0);
  for (int i=0;i<n;i++,Tii++,Sii+=S->bands()) *Sii=*Tii;
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const TridiagonalMatrix<float,complex<float> > &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,max(nsub,1),max(nsub,1),
    complex_float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=stride) {
        *Sji=conj(*Aij);
      }
    }
  }
  int S_nt=S->bands();
  F77NAME(caxpy)(dim,complex_float_mone_,T.diagonalAddr(),1,
    S->addr(0,0),S_nt);
  F77NAME(caxpy)(dim-1,complex_float_mone_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S_nt);
  F77NAME(caxpy)(dim-1,complex_float_mone_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S_nt);
  return S;
}

template<> BandMatrix<float,complex<float> >* operator-(
const TridiagonalMatrix<float,complex<float> > &T,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,complex<float> > *S=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,max(B.subDiags(),1),
    max(B.subDiags(),1),complex_float_zero_);
  int S_nt=S->bands();
  F77NAME(ccopy)(n,T.diagonalAddr(),1,S->addr(0,0),S_nt);
  F77NAME(ccopy)(n-1,T.lowerDiagonalAddr(),1,S->addr(1,0),S_nt);
  F77NAME(ccopy)(n-1,T.upperDiagonalAddr(),1,S->addr(0,1),S_nt);
  int stride=S->bands()-1;
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<n && B.subDiags()>0) {
      const complex<float> *Bij=B.addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=stride) {
        *Sji-=conj(*Bij);
      }
    }
  }
  return S;
}

template<> SymmetricMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const SymmetricMatrix<float,complex<float> > &H) const {
  CHECK_SAME(dim,H.size(0));
  SymmetricMatrix<float,complex<float> > *S=makeSymmetricMatrix();
  for (int j=0;j<dim;j++) {
    F77NAME(caxpy)(dim-j,complex_float_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricMatrix<float,complex<float> >* operator-(
const SymmetricMatrix<float,complex<float> > &H,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SymmetricMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(n-j,H.addr(j,j),1,S->addr(j,j),1);
    F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<float> *Aij=addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      F77NAME(caxpy)(j+1,complex_float_mone_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<float> *Aij=addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      if (j>0) {
        F77NAME(caxpy)(j,complex_float_mone_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      (*S)(j,j)-=complex_float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(j+1,U.addr(0,j),1,S->addr(0,j),1);
    }
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,
        B.addr(j,j),1,S->addr(j,j),1);
      if (j+1<n && B.subDiags()>0) {
        const complex<float> *Bij=B.addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
          *Sji-=conj(*Bij);
        }
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) F77NAME(ccopy)(j,U.addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)=complex_float_one_;
    }
    for (int j=0;j<n;j++) {
      F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,
        B.addr(j,j),1,S->addr(j,j),1);
      if (j+1<n && B.subDiags()>0) {
        const complex<float> *Bij=B.addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
          *Sji-=conj(*Bij);
        }
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<float> *Aij=addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      F77NAME(caxpy)(dim-j,complex_float_mone_,L.addr(j,j),1,
        S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j+1<dim && nsub>0) {
        const complex<float> *Aij=addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
          *Sji=conj(*Aij);
        }
      }
      (*S)(j,j)-=complex_float_one_;
      if (j+1<dim) {
        F77NAME(caxpy)(dim-j-1,complex_float_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(ccopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,
        B.addr(j,j),1,S->addr(j,j),1);
      if (j+1<n && B.subDiags()>0) {
        const complex<float> *Bij=B.addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
          *Sji-=conj(*Bij);
        }
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=complex_float_one_;
      if (j+1<n) {
        F77NAME(ccopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,
        B.addr(j,j),1,S->addr(j,j),1);
      if (j+1<n && B.subDiags()>0) {
        const complex<float> *Bij=B.addr(j+1,j);
        complex<float> *Sji=S->addr(j,j+1);
        for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
          *Sji-=conj(*Bij);
        }
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator-(
const Matrix<float,complex<float> > &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(ccopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim && nsub>0) {
      const complex<float> *Aij=addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(dim-1,j+nsub);i++,Aij++,Sji+=dim) {
        *Sji=conj(*Aij);
      }
    }
    F77NAME(caxpy)(dim,complex_float_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,complex<float> >* operator-(
const Matrix<float,complex<float> > &M,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,complex<float> > *S=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(min(B.bands(),n-j),complex_float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<n && B.subDiags()>0) {
      const complex<float> *Bij=B.addr(j+1,j);
      complex<float> *Sji=S->addr(j,j+1);
      for (int i=j+1;i<=min(n-1,j+B.subDiags());i++,Bij++,Sji+=n) {
        *Sji-=conj(*Bij);
      }
    }
  }
  return S;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const SymmetricBandMatrix<float,complex<float> > &B) const {
// compute by bordering: note that
// [ sigma s^H ] [ tau t^H ] = [ sigma tau + s^H t , sigma t^H + s^H T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^H +  S  T ]
  CHECK_SAME(dim,B.dim);
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,nsub+B.nsub,nsub+B.nsub,
    complex_float_zero_);
  int nb=min(nsub,B.nsub);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*B(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(dim-1,k+B.nsub);kk++) { // S t
        int ibeg=max(k+1,kk-B.nsub);
        if (kk>ibeg) {
          complex<float> Bkkk=B(kk,k);
          const complex<float> *Akki=addr(kk,ibeg);
          complex<float> *Pik=P->addr(ibeg,k);
          for (int i=ibeg;i<kk;i++,Akki+=nsub,Pik++) {
            *Pik+=conj(*Akki)*Bkkk;
          }
        }
        ibeg=max(kk,ibeg);
        F77NAME(caxpy)(min(kk+nsub,dim-1)-ibeg+1,B(kk,k),
          addr(ibeg,kk),1,P->addr(ibeg,k),1);
      }
      for (int j=k+1;j<=min(k+nsub+B.nsub,dim-1);j++) { // s^H T
        int kbeg=max(k+1,j-B.nsub);
        int kend=min(j-1,k+nsub);
        if (kbeg<=kend) {
          (*P)(k,j)=conj(F77NAME(cdotu)(kend-kbeg+1,addr(kbeg,k),1,
            B.addr(j,kbeg),B.nt-1));
        }
        kbeg=max(kbeg,j);
        kend=min(dim-1,min(k+nsub,j+B.nsub));
        if (kend>=kbeg) {
          (*P)(k,j)+=
            F77NAME(cdotc)(kend-kbeg+1,addr(kbeg,k),1,B.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^H t
        F77NAME(cdotc)(min(nb,dim-k-1),addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(cgerc)(min(dim-k,nt),min(dim-k,B.nt),complex_float_one_,
      addr(k,k),1,B.addr(k,k),1,P->addr(k,k),P->bands()-1);
  }
  return P;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const BandMatrix<float,complex<float> > &B) const {
// compute by bordering: note that
// [ sigma s^H ] [ beta c^T ] = [ sigma beta + s^H b , sigma c^T + s^H B ]
// [   s    S  ] [  b    B  ] = [   s   beta +  S  b ,   s   c^T +  S  B ]
  CHECK_SAME(dim,B.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,nsub+B.subDiags(),
    nsub+B.supDiags(),complex_float_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*B(dim-1,dim-1);
  int nb=min(nsub,B.subDiags());
  for (int k=dim-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(dim-1,k+B.subDiags());kk++) { // S b
        int ibeg=max(k+1,kk-B.supDiags());
        if (kk>ibeg) {
          complex<float> Bkkk=B(kk,k);
          const complex<float> *Akki=addr(kk,ibeg);
          complex<float> *Pik=P->addr(ibeg,k);
          for (int i=ibeg;i<kk;i++,Akki+=nsub,Pik++) {
            *Pik+=conj(*Akki)*Bkkk;
          }
        }
        ibeg=max(kk,ibeg);
        F77NAME(caxpy)(min(kk+nsub,dim-1)-ibeg+1,B(kk,k),
          addr(ibeg,kk),1,P->addr(ibeg,k),1);
      }
      for (int j=k+1;j<=min(k+P->supDiags(),dim-1);j++) { // s^H B
        int kbeg=max(k+1,j-B.supDiags());
        int kend=min(dim-1,min(k+nsub,j+B.subDiags()));
        if (kbeg<=kend) {
          (*P)(k,j)=F77NAME(cdotc)(kend-kbeg+1,addr(kbeg,k),1,
            B.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^H b
        F77NAME(cdotc)(min(nb,dim-k-1),addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(cgeru)(min(dim-k,nt),min(dim-k,B.supDiags()+1),
      complex_float_one_,addr(k,k),1,B.addr(k,k),1,
      P->addr(k,k),P->bands()-1);
  }
  return P;
}

template<> BandMatrix<float,complex<float> >* operator*(
const BandMatrix<float,complex<float> > &B,
const SymmetricBandMatrix<float,complex<float> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,B.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,S.subDiags()+B.subDiags(),
      S.subDiags()+B.supDiags(),complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-S.subDiags());k<=min(n-1,j+S.subDiags());k++) {
      int ibeg=max(0,k-B.supDiags());
      F77NAME(caxpy)(min(n-1,k+B.subDiags())-ibeg+1,S(k,j),
        B.addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const UpperHessenbergMatrix<float,complex<float> > &H) const {
// compute by bordering: note that
// [ sigma s^H ] [     eta_11 h^T ]
// [   s    S  ] [ e_0 eta_21  H  ]
// = [ sigma eta_11 + s^H e_0 eta_21 , sigma h^T + s^H H ]
// = [   s   eta_11 +  S  e_0 eta_21 ,   s   h^T +  S  H ]
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*H(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nsub>0) {
      F77NAME(caxpy)(min(dim-k-1,nsub+1),H(k+1,k),addr(k+1,k+1),1,
        P->addr(k+1,k),1); // S  e_0 eta_21
      for (int j=k+1;j<dim;j++) { // s^H H
        (*P)(k,j)=
          F77NAME(cdotc)(min(dim-k-1,min(j-k+1,nsub)),addr(k+1,k),1,
            H.addr(k+1,j),1);
      }
      (*P)(k,k)=conj((*this)(k+1,k))*H(k+1,k); // s^H e_0 eta_21
    }
    F77NAME(cgeru)(min(dim-k,nsub+1),dim-k,complex_float_one_,
      addr(k,k),1,H.addr(k,k),dim,P->addr(k,k),dim);
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const UpperHessenbergMatrix<float,complex<float> > &H,
const SymmetricBandMatrix<float,complex<float> > &S) {
  int n=S.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-S.subDiags());k<=min(n-1,j+S.subDiags());k++) {
      F77NAME(caxpy)(min(n,k+2),S(k,j),H.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const DiagonalMatrix<float,complex<float> > &D) const {
  BandMatrix<float,complex<float> > *P=makeBandMatrix();
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsub);
    F77NAME(cscal)(min(dim-1,j+nsub)-ibeg+1,D[j],P->addr(ibeg,j),1);
  }
  return P;
}

template<> BandMatrix<float,complex<float> >* operator*(
const DiagonalMatrix<float,complex<float> > &D,
const SymmetricBandMatrix<float,complex<float> > &S) {
  int n=S.size(0);
  BandMatrix<float,complex<float> > *P=S.makeBandMatrix();
  int stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    int jbeg=max(0,i-S.subDiags());
    F77NAME(cscal)(min(n-1,i+S.subDiags())-jbeg+1,D[i],
      P->addr(i,jbeg),stride);
  }
  return P;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
// compute by bordering: note that
// [ sigma s^H ] [     tau    , conj(lambda) e_0^T ]
// [   s    S  ] [ e_0 lambda ,                 T  ]
//   = [ sigma tau + s^H e_0 lambda , sigma conj(lambda) e_0^T + s^H T ]
//   = [   s   tau +  S  e_o lambda ,   s   conj(lambda) e_0^T +  S  T ]
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,nsub+1,nsub+1,
    complex_float_zero_);
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
    F77NAME(caxpy)(min(dim-k-1,nsub+1),T.lowerDiagonalValue(k),
      addr(k+1,k+1),1,P->addr(k+1,k),1); // S e_0 lambda
    F77NAME(caxpy)(min(dim-k,nsub+1),T.diagonalValue(k),addr(k,k),1,
      P->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(caxpy)(min(dim-k,nsub+1),T.upperDiagonalValue(k),addr(k,k),1,
      P->addr(k,k+1),1);
     // [ sigma ] conj(lambda)
     // [   s   ]
  }
  return P;
}

template<> BandMatrix<float,complex<float> >* operator*(
const SymmetricTridiagonalMatrix<float,complex<float> > &T,
const SymmetricBandMatrix<float,complex<float> > &S) {
// compute by bordering: note that
// [     tau    , conj(lambda) e_0^T ] [ sigma s^H ]
// [ e_0 lambda ,            T       ] [   s    S  ]
// = [ tau sigma + conj(lambda) e_0^T s, tau s^H  + conj(lambda) e_0^T S ]
// = [ e_0 lambda sigma +          T s , e_0 lambda s^H +            T S ]
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=S.subDiags();
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,nb+1,nb+1,complex_float_zero_);
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
    complex<float> Tkkp1=T.upperDiagonalValue(k);
    const complex<float> *Sjkp1=S.addr(k+1,k+1);
    complex<float> *Pkj=P->addr(k,k+1);
    for (int j=k+1;j<=min(k+1+nb,n-1);j++,Sjkp1++,Pkj+=stride) {
      *Pkj+=Tkkp1*conj(*Sjkp1);
    }
      // conj(lambda) e_0^T S = conj(lambda) ( S e_0 )^H
    float Tkk=T.diagonalValue(k);
    const complex<float> *Sjk=S.addr(k,k);
    Pkj=P->addr(k,k);
    for (int j=k;j<=min(n-1,k+nb);j++,Sjk++,Pkj+=stride) {
      *Pkj+=Tkk*conj(*Sjk);
    }
      // tau [ sigma , s^H ]
    complex<float> Tkp1k=T.lowerDiagonalValue(k);
    Sjk=S.addr(k,k);
    complex<float> *Pkp1j=P->addr(k+1,k);
    for (int j=k;j<=min(n-1,k+nb);j++,Sjk++,Pkp1j+=stride) {
      *Pkp1j+=Tkp1k*conj(*Sjk);
    }
      // e_0 lambda [ sigma , s^H ]
  }
  return P;
}

template<> BandMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const TridiagonalMatrix<float,complex<float> > &T) const {
// compute by bordering: note that
// [ sigma s^H ] [    tau     upsilon e_0^T ]
// [   s    S  ] [ e_0 lambda       T       ]
//   = [ sigma tau + s^H e_0 lambda , sigma upsilon e_0^T + s^H T ]
//   = [     s tau +   S e_0 lambda ,     s upsilon e_0^T +   S T ]
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(dim,nsub+1,nsub+1,
    complex_float_zero_);
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
    F77NAME(caxpy)(min(dim-k-1,nsub+1),T(k+1,k),addr(k+1,k+1),1,
      P->addr(k+1,k),1); // S e_0 lambda
    F77NAME(caxpy)(min(dim-k,nsub+1),T(k,k),addr(k,k),1,P->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(caxpy)(min(dim-k,nsub+1),T(k,k+1),addr(k,k),1,
      P->addr(k,k+1),1);
      // [ sigma ] lambda
      // [   s   ]
  }
  return P;
}

template<> BandMatrix<float,complex<float> >* operator*(
const TridiagonalMatrix<float,complex<float> > &T,
const SymmetricBandMatrix<float,complex<float> > &S) {
// compute by bordering: note that
// [     tau    , upsilon e_0^T ] [ sigma s^H ]
// [ e_0 lambda ,       T       ] [   s    S  ]
//   = [ tau sigma + upsilon e_0^T s , tau s^H        + upsilon e_0^T S ]
//   = [ e_0 lambda sigma +      T s , e_0 lambda s^H +             T S ]
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=S.subDiags();
  BandMatrix<float,complex<float> > *P=OPERATOR_NEW
    BandMatrix<float,complex<float> >(n,nb+1,nb+1,complex_float_zero_);
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
    complex<float> Tkkp1=T(k,k+1);
    const complex<float> *Sjkp1=S.addr(k+1,k+1);
    complex<float> *Pkj=P->addr(k,k+1);
    for (int j=k+1;j<=min(n-1,k+nb+1);j++,Sjkp1++,Pkj+=stride) {
      *Pkj+=Tkkp1*conj(*Sjkp1);
    }
      // upsilon e_0^T S = upsilon ( S e_0 )^T
    complex<float> Tkk=T(k,k);
    const complex<float> *Sjk=S.addr(k,k);
    Pkj=P->addr(k,k);
    for (int j=k;j<=min(n-1,k+nb);j++,Sjk++,Pkj+=stride) {
      *Pkj+=Tkk*conj(*Sjk);
    }
      // tau [ sigma , s^H ]
    complex<float> Tkp1k=T(k+1,k);
    Sjk=S.addr(k,k);
    complex<float> *Pkp1j=P->addr(k+1,k);
    for (int j=k;j<=min(n-1,k+nb);j++,Sjk++,Pkp1j+=stride) {
      *Pkp1j+=Tkp1k*conj(*Sjk);
    }
      // e_0 lambda [ sigma , s^H ]
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const SymmetricMatrix<float,complex<float> > &S) const {
// compute by bordering: note that
// [ sigma s^H ] [ tau t^H ] = [ sigma tau + s^H t , sigma t^H + s^H T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^H +  S  T ]
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(dim,complex_float_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*S(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nsub>0) {
      F77NAME(chbmv)('L',dim-k-1,nsub,complex_float_one_,
        addr(k+1,k+1),nt,S.addr(k+1,k),1,complex_float_zero_,
        P->addr(k+1,k),1); // S t
      for (int j=k+1;j<dim;j++) { // s^H T
        int kend=min(j-1,k+nsub);
        if (k<kend) {
          (*P)(k,j)=
            conj(F77NAME(cdotu)(kend-k,addr(k+1,k),1,S.addr(j,k+1),dim));
        }
        int kbeg=max(k+1,j);
        kend=min(dim-1,k+nsub);
        if (kbeg<=kend) {
          (*P)(k,j)+=
            F77NAME(cdotc)(kend-kbeg+1,addr(kbeg,k),1,S.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^H t
        F77NAME(cdotc)(min(nsub,dim-k-1),addr(k+1,k),1,S.addr(k+1,k),1);
    }
    F77NAME(cgerc)(min(dim-k,nt),dim-k,complex_float_one_,addr(k,k),1,
      S.addr(k,k),1,P->addr(k,k),dim);
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SymmetricMatrix<float,complex<float> > &S,
const SymmetricBandMatrix<float,complex<float> > &B) {
// compute by bordering: note that
// [ sigma s^H ] [ tau t^H ] = [ sigma tau + s^H t , sigma t^H + s^H T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^H +  S  T ]
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  (*P)(n-1,n-1)=S(n-1,n-1)*B(n-1,n-1);
  int nb=B.subDiags();
  for (int k=n-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(n-1,k+nb);kk++) { // S t 
        int ibeg=k+1;
        complex<float> Bkkk=B(kk,k);
        if (kk>ibeg) {
          const complex<float> *Skki=S.addr(kk,ibeg);
          complex<float> *Pik=P->addr(ibeg,k);
          for (int i=ibeg;i<kk;i++,Skki+=n,Pik++) {
            *Pik+=conj(*Skki)*Bkkk;
          }
        }
        F77NAME(caxpy)(n-kk,Bkkk,S.addr(kk,kk),1,P->addr(kk,k),1);
      }
      for (int j=k+1;j<n;j++) { // s^H T
        int kbeg=max(k+1,j-nb);
        if (kbeg<j) {
          (*P)(k,j)=conj(F77NAME(cdotu)(j-kbeg,S.addr(kbeg,k),1,
            B.addr(j,kbeg),B.bands()-1));
        }
        kbeg=max(kbeg,j);
        (*P)(k,j)+=F77NAME(cdotc)(min(n-1,j+nb)-kbeg+1,
          S.addr(kbeg,k),1,B.addr(kbeg,j),1);
      }
      (*P)(k,k)= // s^H t
        F77NAME(cdotc)(min(nb,n-k-1),S.addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(cgerc)(n-k,min(n-k,B.bands()),complex_float_one_,
      S.addr(k,k),1,B.addr(k,k),1,P->addr(k,k),n);
  }
  return P;
}

template<> Matrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U) const {
  CHECK_SAME(dim,U.size(0));
  int n=U.size(1);
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(dim,n,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,dim-1);k++) {
        int ibeg=max(0,k-nsub);
        complex<float> Ukj=U(k,j);
        if (ibeg<k) {
          const complex<float> *Aki=addr(k,ibeg);
          complex<float> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<k;i++,Aki+=nsub,Pij++) {
            *Pij+=conj(*Aki)*Ukj;
          }
        }
        ibeg=max(ibeg,k);
        F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,Ukj,addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,dim);k++) {
        int ibeg=max(0,k-nsub);
        complex<float> Ukj=U(k,j);
        if (ibeg<k) {
          const complex<float> *Aki=addr(k,ibeg);
          complex<float> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<k;i++,Aki+=nsub,Pij++) {
            *Pij+=conj(*Aki)*Ukj;
          }
        }
        ibeg=max(ibeg,k);
        F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,Ukj,addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
      if (j<dim) {
        int ibeg=max(0,j-nsub);
        if (ibeg<j) {
          const complex<float> *Aji=addr(j,ibeg);
          complex<float> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<j;i++,Aji+=nsub,Pij++) {
            *Pij+=conj(*Aji);
          }
        }
        ibeg=max(ibeg,j);
        F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,complex_float_one_,
          addr(ibeg,j),1,P->addr(ibeg,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >* operator*(
const UpperTrapezoidalMatrix<float,complex<float> > &U,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,complex<float> >*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(caxpy)(min(k+1,m),B(k,j),U.addr(0,k),1,P->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(caxpy)(min(k,m),B(k,j),U.addr(0,k),1,P->addr(0,j),1);
        if (k<m) (*P)(k,j)+=B(k,j);
      }
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L) const {
  CHECK_SAME(dim,L.size(0));
  int n=L.size(1);
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(dim,n,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<dim;k++) {
        int ibeg=max(0,k-nsub);
        complex<float> Lkj=L(k,j);
        if (ibeg<k) {
          const complex<float> *Aki=addr(k,ibeg);
          complex<float> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<k;i++,Aki+=nsub,Pij++) {
            *Pij+=conj(*Aki)*Lkj;
          }
        }
        ibeg=max(ibeg,k);
        F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,Lkj,addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-nsub);
      if (ibeg<j) {
        const complex<float> *Aji=addr(j,ibeg);
        complex<float> *Pij=P->addr(ibeg,j);
        for (int i=ibeg;i<j;i++,Aji+=nsub,Pij++) {
          *Pij+=conj(*Aji);
        }
      }
      ibeg=max(ibeg,j);
      F77NAME(caxpy)(min(dim-1,j+nsub)-ibeg+1,complex_float_one_,
        addr(ibeg,j),1,P->addr(ibeg,j),1);
      for (int k=j+1;k<dim;k++) {
        int ibeg=max(0,k-nsub);
        complex<float> Lkj=L(k,j);
        if (ibeg<k) {
          const complex<float> *Aki=addr(k,ibeg);
          complex<float> *Pij=P->addr(ibeg,j);
          for (int i=ibeg;i<k;i++,Aki+=nsub,Pij++) {
            *Pij+=conj(*Aki)*Lkj;
          }
        }
        ibeg=max(ibeg,k);
        F77NAME(caxpy)(min(dim-1,k+nsub)-ibeg+1,Lkj,addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >* operator*(
const LowerTrapezoidalMatrix<float,complex<float> > &L,
const SymmetricBandMatrix<float,complex<float> > &B) {
int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,complex<float> >*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(caxpy)(m-k,B(k,j),L.addr(k,k),1,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        (*P)(k,j)+=B(k,j);
        if (k+1<m) {
          F77NAME(caxpy)(m-k-1,B(k,j),L.addr(k+1,k),1,P->addr(k+1,j),1);
        }
      }
    }
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const SquareMatrix<float,complex<float> > &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<float,complex<float> > *P=
    OPERATOR_NEW SquareMatrix<float,complex<float> >(dim);
  for (int j=0;j<dim;j++) {
    F77NAME(chbmv)('L',dim,nsub,complex_float_one_,addr(),nt,
      S.addr(0,j),1,complex_float_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> SquareMatrix<float,complex<float> >* operator*(
const SquareMatrix<float,complex<float> > &S,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,complex<float> > *P=OPERATOR_NEW
    SquareMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(caxpy)(n,B(k,j),S.addr(0,j),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Matrix<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const Matrix<float,complex<float> > &M) const {
  CHECK_SAME(dim,M.size(0));
  int n=M.size(1);
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(dim,n);
  for (int j=0;j<dim;j++) {
    F77NAME(chbmv)('L',dim,nsub,complex_float_one_,addr(),nt,
      M.addr(0,j),1,complex_float_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<float,complex<float> >* operator*(
const Matrix<float,complex<float> > &M,
const SymmetricBandMatrix<float,complex<float> > &B) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,complex<float> > *P=OPERATOR_NEW
    Matrix<float,complex<float> >(m,n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(caxpy)(m,B(k,j),M.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Vector<float,complex<float> >*
SymmetricBandMatrix<float,complex<float> >::operator*(
const Vector<float,complex<float> > &v) const {
  CHECK_SAME(dim,v.size());
  Vector<float,complex<float> > *p=
    OPERATOR_NEW Vector<float,complex<float> >(dim);
  F77NAME(chbmv)('L',dim,nsub,complex_float_one_,addr(),nt,v.addr(),1,
    complex_float_zero_,p->addr(),1);
  return p;
}

template<> void SymmetricBandMatrix<float,complex<float> >::sbmv(
complex<float> alpha,const Vector<float,complex<float> > &x,
complex<float> beta,Vector<float,complex<float> > &b) const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  F77NAME(chbmv)('L',dim,nsub,alpha,addr(),nt,x.addr(),1,beta,
    b.addr(),1);
}

template<> void SymmetricBandMatrix<float,complex<float> >::sbmm(
complex<float> alpha,const Matrix<float,complex<float> > &X,
complex<float> beta,Matrix<float,complex<float> > &B,char side)
const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      F77NAME(chbmv)('L',dim,nsub,alpha,addr(),nt,X.addr(0,j),1,beta,
        B.addr(0,j),1);
    }
  } else {
    CHECK_SAME(n,dim);
    if (abs(beta)==complex_float_zero_) B=complex_float_zero_;
    else B*=beta;
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-nsub);k<=min(dim-1,j+nsub);k++) {
        F77NAME(caxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
          B.addr(0,j),1);
      }
    }
  }
}

template<> float
SymmetricBandMatrix<float,complex<float> >::normFrobenius() const {
  float *work=0;
  return F77NAME(clanhb)('F','L',dim,nsub,addr(),nt,work);
}

template<> float
SymmetricBandMatrix<float,complex<float> >::normInfinity() const {
  float *work=OPERATOR_NEW_BRACKET(float,dim);;
  float val=F77NAME(clanhb)('I','L',dim,nsub,addr(),nt,work);
  delete [] work; work=0;
  return val;
}

template<> float
SymmetricBandMatrix<float,complex<float> >::normMaxEntry() const {
  float *work=0;
  return F77NAME(clanhb)('M','L',dim,nsub,addr(),nt,work);
}

template<> float
SymmetricBandMatrix<float,complex<float> >::normOne() const {
  float *work=OPERATOR_NEW_BRACKET(float,dim);;
  float val=F77NAME(clanhb)('O','L',dim,nsub,addr(),nt,work);
  delete [] work; work=0;
  return val;
}

template<> Vector<float,float>* 
SymmetricBandMatrix<float,complex<float> >::eigenvalues(
OrthogonalMatrix<float,complex<float> > *&Q) const {
  if (Q!=0) CHECK_SAME(dim,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Vector<float,float> *lambda=OPERATOR_NEW Vector<float,float>(dim);
  SymmetricBandMatrix<float,complex<float> > *copy=
    OPERATOR_NEW SymmetricBandMatrix<float,complex<float> >(*this);
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,dim);
  float *rwork=OPERATOR_NEW_BRACKET(float,max(1,3*dim-2));
  int info;
  complex<float> *qa=( Q==0 ? 0 : Q->addr() );
  F77NAME(chbev)(jobz,'L',dim,nsub,copy->addr(),copy->bands(),
    lambda->addr(),qa,dim,work,rwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete [] rwork; rwork=0;
  delete copy; copy=0;
  return lambda;
}

template class SymmetricBandMatrix<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveBandMatrix<float,complex<float> >*
SymmetricPositiveBandMatrix<float,complex<float> >::operator+(
const SymmetricPositiveBandMatrix<float,complex<float> > &B) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,B.size(0));
  SymmetricPositiveBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricPositiveBandMatrix<float,complex<float> >(n,
    max(ns,B.subDiags()),complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(caxpy)(min(B.bands(),n-j),complex_float_one_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricPositiveBandMatrix<float,complex<float> >*
SymmetricPositiveBandMatrix<float,complex<float> >::operator+(
const SymmetricPositiveTridiagonalMatrix<float,complex<float> > &T)
const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricPositiveBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricPositiveBandMatrix<float,complex<float> >(n,max(ns,1),
    complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=T(j,j);
    if (j+1<n) (*S)(j+1,j)+=T(j+1,j);
  }
  return S;
}

template<> SymmetricBandMatrix<float,complex<float> >*
SymmetricPositiveBandMatrix<float,complex<float> >::operator+(
const SymmetricTridiagonalMatrix<float,complex<float> > &T) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricBandMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricBandMatrix<float,complex<float> >(n,max(ns,1),
    complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=T(j,j);
    if (j+1<n) (*S)(j+1,j)+=T(j+1,j);
  }
  return S;
}

template<> SymmetricPositiveMatrix<float,complex<float> >*
SymmetricPositiveBandMatrix<float,complex<float> >::operator+(
const SymmetricPositiveMatrix<float,complex<float> > &M) const {
  int n=size(0);
  int nb=bands();
  CHECK_SAME(n,M.size(0));
  SymmetricPositiveMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricPositiveMatrix<float,complex<float> >(n,
    complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(caxpy)(n-j,complex_float_one_,M.addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricMatrix<float,complex<float> >*
SymmetricPositiveBandMatrix<float,complex<float> >::operator+(
const SymmetricMatrix<float,complex<float> > &T) const {
  int n=size(0);
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,complex<float> > *S=OPERATOR_NEW
    SymmetricMatrix<float,complex<float> >(n,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(caxpy)(n-j,complex_float_one_,T.addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> float
SymmetricPositiveBandMatrix<float,complex<float> >::equilibrate(
Vector<float,float> &s,float &scond) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,s.size());
  float amax;
  int info;
  F77NAME(cpbequ)('L',n,ns,addr(),nb,s.addr(),scond,amax,info);
  CHECK_TEST(info==0);
  return amax;
}

template<> float
SymmetricPositiveBandMatrix<float,complex<float> >
::reciprocalConditionNumber() const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  SymmetricPositiveBandMatrix<float,complex<float> > *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<float,complex<float> >(
    *this);
  int info;
  F77NAME(cpbtrf)('L',n,ns,BF->addr(),nb,info);
  CHECK_TEST(info==0);

  float anorm=normOne();
  float rcond=numeric_limits<float>::infinity();
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  F77NAME(cpbcon)('L',n,ns,BF->addr(),nb,anorm,rcond,work,rwork,info);
  delete [] work; work=0;
  delete [] rwork; rwork=0;
  delete BF; BF=0;
  return rcond;
}

template<> void
SymmetricPositiveBandMatrix<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,char) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  SymmetricPositiveBandMatrix<float,complex<float> > *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<float,complex<float> >(
    *this);
  x.copy(b);
  int info;
  F77NAME(cpbsv)('L',n,ns,1,BF->addr(),nb,x.addr(),n,info);
  CHECK_TEST(info==0);
  delete BF; BF=0;
}

template<> void
SymmetricPositiveBandMatrix<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,char side,char) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  SymmetricPositiveBandMatrix<float,complex<float> > *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<float,complex<float> >(
    *this);
  int info;
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    CHECK_SAME(nrhs,X.size(1))
    X.copy(B);
    F77NAME(cpbsv)('L',n,ns,1,BF->addr(),nb,X.addr(),n,info);
    CHECK_TEST(info==0);
  } else {
    int nrhs=B.size(0);
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    CHECK_SAME(nrhs,X.size(0))
    X.copy(B);
    F77NAME(cpbtrf)('L',n,ns,BF->addr(),nb,info);
    CHECK_TEST(info==0);
    for (int j=0;j<nrhs;j++) { // dpbtrs loop 20:
      F77NAME(ctbsv)('L','N','N',n,ns,BF->addr(),nb,X.addr(0,j),nrhs);
      F77NAME(ctbsv)('L','T','N',n,ns,BF->addr(),nb,X.addr(0,j),nrhs);
    }
  }
  delete BF; BF=0;
}

template class SymmetricPositiveBandMatrix<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SpecializedMatrix.C"
template class CompanionMatrix<float,complex<float> >;
template class HadamardMatrix<float,complex<float> >;
template class HankelMatrix<float,complex<float> >;
template class HilbertMatrix<float,complex<float> >;
template class KahanMatrix<float,complex<float> >;
template class LauchliMatrix<float,complex<float> >;
template class PascalMatrix<float,complex<float> >;
template class RosserMatrix<float,complex<float> >;
template class ToeplitzMatrix<float,complex<float> >;
template class SymmetricToeplitzMatrix<float,complex<float> >;
template class VandermondeMatrix<float,complex<float> >;
template class WilkinsonMatrix<float,complex<float> >;
template void testSpecializedMatrix(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "GaussianFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> GaussianFactorization<float,complex<float>,
SquareMatrix<float,complex<float> > >::GaussianFactorization(
const SquareMatrix<float,complex<float> >& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),piv_op(po),equ_op(eo),
equed('N'),colcnd(numeric_limits<float>::infinity()),
rowcnd(numeric_limits<float>::infinity()) {
  int n=A.size(0);
  LU=OPERATOR_NEW SquareMatrix<float,complex<float> >(n);
  LU->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    r=OPERATOR_NEW Vector<float,float>(n);
    c=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    F77NAME(cgeequ)(n,n,LU->addr(),n,r->addr(),c->addr(),rowcnd,colcnd,
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
      F77NAME(claqge)(n,n,LU->addr(),n,r->addr(),c->addr(),
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
        if (abs((*LU)(k,k))==float_zero_) {
          info=k+1;
          break;
        }
	complex<float> alpha=complex_float_one_/(*LU)(k,k);
	int nrows=n-k-1;
	if (nrows>0) {
	  F77NAME(cscal)(nrows,alpha,LU->addr(k+1,k),1);
	  for (int j=k+1;j<n;j++) {
	    alpha=-(*LU)(k,j);
	    F77NAME(caxpy)(nrows,alpha,LU->addr(k+1,k),1,
			   LU->addr(k+1,j),1);
	  }
	}
      }
      break;
    }
    case Factorization::PIVOT_ROWS_AND_COLUMNS: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      jpiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(cgetc2)(n,LU->addr(),n,ipiv,jpiv,info);
      break;
    }
    default: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(cgetrf)(n,n,LU->addr(),n,ipiv,info);
      break;
    }
  }
  float *work=0;
  if (info>0) { // see dgesvx
    rpvgrw=F77NAME(clantr)('M','U','N',info,info,LU->addr(),n,work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
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
    rpvgrw=F77NAME(clantr)('M','U','N',n,n,LU->addr(),n,work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
  }
}

template<> void GaussianFactorization<float,complex<float>,
SquareMatrix<float,complex<float> > >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,Factorization::TRANSPOSE_OPTION to) {
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
      const float *ri=r->addr();
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
    if (ipiv!=0) F77NAME(claswp)(1,x.addr(),n,1,n-1,ipiv,1);
    F77NAME(ctrsm)('L','L','N','U',n,1,complex_float_one_,LU->addr(),n,
      x.addr(),n);
    F77NAME(ctrsm)('L','U','N','N',n,1,complex_float_one_,LU->addr(),n,
      x.addr(),n);
    if (jpiv!=0) F77NAME(claswp)(1,x.addr(),n,1,n-1,jpiv,-1);
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ci=c->addr();
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
  } else {
    // A^T X = B ==> U^T L^T Q^T R^{-1} X = P^T C B
    // A^H X = B ==> U^H L^H Q^T R^{-1} X = P^T C B
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ci=c->addr();
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
    if (jpiv!=0) F77NAME(claswp)(1,x.addr(),n,1,n-1,jpiv,1);
    if (to==Factorization::TRANSPOSE) {
      F77NAME(ctrsm)('L','U','T','N',n,1,complex_float_one_,LU->addr(),n,
        x.addr(),n);
      F77NAME(ctrsm)('L','L','T','U',n,1,complex_float_one_,LU->addr(),n,
        x.addr(),n);
    } else {
      F77NAME(ctrsm)('L','U','C','N',n,1,complex_float_one_,LU->addr(),n,
        x.addr(),n);
      F77NAME(ctrsm)('L','L','C','U',n,1,complex_float_one_,LU->addr(),n,
        x.addr(),n);
    }
    if (ipiv!=0) F77NAME(claswp)(1,x.addr(),n,1,n-1,ipiv,-1);
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ri=r->addr();
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
  }
}

template<> void GaussianFactorization<float,complex<float>,
SquareMatrix<float,complex<float> > >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Factorization::TRANSPOSE_OPTION to,
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
          const float *ri=r->addr();
          complex<float> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
      if (ipiv!=0) F77NAME(claswp)(n,X.addr(),m,1,m-1,ipiv,1);
      F77NAME(ctrsm)('L','L','N','U',m,n,complex_float_one_,
                     LU->addr(),m,X.addr(),m);
      F77NAME(ctrsm)('L','U','N','N',m,n,complex_float_one_,
                     LU->addr(),m,X.addr(),m);
      if (jpiv!=0) F77NAME(claswp)(n,X.addr(),m,1,m-1,jpiv,-1);
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ci=c->addr();
          complex<float> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
    } else {
      // A^T X = B ==> U^T L^T Q^T R^{-1} X = P^T C B
      // A^H X = B ==> U^H L^H Q^T R^{-1} X = P^T C B
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ci=c->addr();
          complex<float> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
      if (jpiv!=0) F77NAME(claswp)(n,X.addr(),m,1,m-1,jpiv,1);
      if (to==Factorization::TRANSPOSE) {
        F77NAME(ctrsm)('L','U','T','N',m,n,complex_float_one_,
                       LU->addr(),m,X.addr(),m);
        F77NAME(ctrsm)('L','L','T','U',m,n,complex_float_one_,
                       LU->addr(),m,X.addr(),m);
      } else {
        F77NAME(ctrsm)('L','U','C','N',m,n,complex_float_one_,
                       LU->addr(),m,X.addr(),m);
        F77NAME(ctrsm)('L','L','C','U',m,n,complex_float_one_,
                       LU->addr(),m,X.addr(),m);
      }
      if (ipiv!=0) F77NAME(claswp)(n,X.addr(),m,1,m-1,ipiv,-1);
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ri=r->addr();
          complex<float> *Xij=X.addr(0,j);
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
        for (int j=0;j<n;j++) F77NAME(csscal)(m,(*c)[j],X.addr(0,j),1);
      }
      if (jpiv!=0) {
        for (int j=n-2;j>=0;j--) {
          cout << "jpiv[" << j << "] = " << jpiv[j] << endl;
          if (j!=jpiv[j]-1) {
            F77NAME(cswap)(m,X.addr(0,j),1,X.addr(0,jpiv[j]-1),1);
          }
        }
      }
      F77NAME(ctrsm)('R','U','N','N',m,n,complex_float_one_,
                     LU->addr(),n,X.addr(),m);
      F77NAME(ctrsm)('R','L','N','U',m,n,complex_float_one_,
                     LU->addr(),n,X.addr(),m);
      if (ipiv!=0) {
        for (int i=0;i<n-1;i++) {
          if (i!=ipiv[i]-1) {
            F77NAME(cswap)(m,X.addr(0,i),1,X.addr(0,ipiv[i]-1),1);
          }
        }
      }
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(csscal)(m,(*r)[j],X.addr(0,j),1);
      }
    } else {
      // X A^T = B ==> X C^{-1} P U^T L^T = B R Q
      // X A^H = B ==> X C^{-1} P U^H L^H = B R Q
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(csscal)(m,(*r)[j],X.addr(0,j),1);
      }
      if (ipiv!=0) {
        for (int i=n-2;i>=0;i--) {
          if (i!=ipiv[i]-1) {
            F77NAME(cswap)(m,X.addr(0,i),1,X.addr(0,ipiv[i]-1),1);
          }
        }
      }
      if (to==Factorization::TRANSPOSE) {
        F77NAME(ctrsm)('R','L','T','U',m,n,complex_float_one_,
                       LU->addr(),n,X.addr(),m);
        F77NAME(ctrsm)('R','U','T','N',m,n,complex_float_one_,
                       LU->addr(),n,X.addr(),m);
      } else {
        F77NAME(ctrsm)('R','L','C','U',m,n,complex_float_one_,
                       LU->addr(),n,X.addr(),m);
        F77NAME(ctrsm)('R','U','C','N',m,n,complex_float_one_,
                       LU->addr(),n,X.addr(),m);
      }
      if (jpiv!=0) {
        for (int j=0;j<n-1;j++) {
          if (j!=jpiv[j]-1) {
            F77NAME(cswap)(m,X.addr(0,j),1,X.addr(0,jpiv[j]-1),1);
          }
        }
      }
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(csscal)(m,(*c)[j],X.addr(0,j),1);
      }
    }
  }
}

template<> float GaussianFactorization<float,complex<float>,
SquareMatrix<float,complex<float> > >::reciprocalConditionNumber(
Factorization::CONDITION_NUMBER_NORM cnn) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  float anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,2*n);
  int info=0;
  F77NAME(cgecon)(norm,n,LU->addr(),n,anorm,rcond,work,rwork,info);
  CHECK_SAME(info,0)
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void GaussianFactorization<float,complex<float>,
SquareMatrix<float,complex<float> > >::improve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float &berr,float &ferr,
Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(n+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > residual(n);
  Vector<float,float> work(n);
  Vector<float,complex<float> > v(n);
  char trans='N';
  if (to==Factorization::TRANSPOSE) trans='T';
  else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(cgemv)(trans,n,n,complex_float_mone_,A_original->addr(),n,
      x.addr(),1,complex_float_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(cgeamv)(trans,n,n,float_one_,A_original->addr(),n,
      x.addr(),1,float_one_,work.addr(),1);

    berr=float_zero_;
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
    F77NAME(caxpy)(n,complex_float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  complex<float> *residuali=residual.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(clacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual,to);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual,to);
    }
  }
  int i=F77NAME(icamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void GaussianFactorization<float,complex<float>,
SquareMatrix<float,complex<float> > >::improve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Vector<float,float> &berr,
Vector<float,float> &ferr,Factorization::TRANSPOSE_OPTION to,
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
  float nz=static_cast<float>(k+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > x(k);
  Vector<float,complex<float> > rhs(k);
  Vector<float,complex<float> > residual(k);
  Vector<float,float> work(k);
  Vector<float,complex<float> > v(k);
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
      F77NAME(ccopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(ccopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        const complex<float> *Xji=X.addr(j,0);
        complex<float> *xi=x.addr();
        const complex<float> *Bji=B.addr(j,0);
        complex<float> *rhsi=rhs.addr();
        for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
          *xi=conj(*Xji);
          *rhsi=conj(*Bji);
        }
      } else {
        F77NAME(ccopy)(k,X.addr(j,0),m,x.addr(),1);
        F77NAME(ccopy)(k,B.addr(j,0),m,rhs.addr(),1);
      }
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(cgemv)(trans,k,k,complex_float_mone_,A_original->addr(),k,
        x.addr(),1,complex_float_one_,residual.addr(),float_one_);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77NAME(cgeamv)(trans,k,k,float_one_,A_original->addr(),k,
        x.addr(),1,float_one_,work.addr(),1);

      float &s=berr[j];
      s=float_zero_;
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
      F77NAME(caxpy)(k,complex_float_one_,residual.addr(),1,x.addr(),1);
      lstres=s;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(ccopy)(k,x.addr(),1,X.addr(0,j),1);
      } else F77NAME(ccopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<float> *Xji=X.addr(j,0);
        const complex<float> *xi=x.addr();
        for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
      } else {
        F77NAME(ccopy)(k,x.addr(),1,X.addr(j,0),m);
      }
    }

    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(clacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual,lto);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual,lto);
      }
    }
    int i=F77NAME(icamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr[j]/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> GaussianFactorization<float,complex<float>,
TridiagonalMatrix<float,complex<float> > >::GaussianFactorization(
const TridiagonalMatrix<float,complex<float> >& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),
piv_op(Factorization::PIVOT_ROWS),equ_op(Factorization::NO_EQUILIBRATION),
equed('N'),colcnd(numeric_limits<float>::infinity()),
rowcnd(numeric_limits<float>::infinity()) {
  int n=A.size(0);
  anormi=A.normInfinity();
  anormm=A.normMaxEntry();
  anormo=A.normOne();
  rpvgrw=float_one_; // not computed
}

template<> void GaussianFactorization<float,complex<float>,
TridiagonalMatrix<float,complex<float> > >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,Factorization::TRANSPOSE_OPTION to) {
//constructor did not factor
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  LU=OPERATOR_NEW TridiagonalMatrix<float,complex<float> >(n);
  LU->copy(*A_original);
  if (&x!=&b) x.copy(b);
  int info;
  if (to==Factorization::NO_TRANSPOSE) {
    if (piv_op==Factorization::PIVOT_ROWS) {
      F77NAME(cgtsv)(n,1,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),x.addr(),n,info);
    } else {
      F77NAME(cgtsvnp)(n,1,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),x.addr(),n,info);
    }
  } else {
    if (to==Factorization::CONJUGATE_TRANSPOSE) {
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,xi++) *xi=conj(*xi);
    }
    if (piv_op==Factorization::PIVOT_ROWS) {
      F77NAME(cgtsv)(n,1,LU->upperDiagonalAddr(),LU->diagonalAddr(),
        LU->lowerDiagonalAddr(),x.addr(),n,info);
    } else {
      F77NAME(cgtsvnp)(n,1,LU->upperDiagonalAddr(),LU->diagonalAddr(),
        LU->lowerDiagonalAddr(),x.addr(),n,info);
    }
    if (to==Factorization::CONJUGATE_TRANSPOSE) {
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,xi++) *xi=conj(*xi);
    }
  }
  CHECK_TEST(info==0);
  delete LU; LU=0;
}

template<> void GaussianFactorization<float,complex<float>,
TridiagonalMatrix<float,complex<float> > >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Factorization::TRANSPOSE_OPTION to,
Factorization::SIDE_OPTION so) {
//constructor did not factor
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  LU=OPERATOR_NEW
    TridiagonalMatrix<float,complex<float> >(A_original->size(0));
  LU->copy(*A_original);
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) { // A X = B
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(cgtsv)(m,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(cgtsvnp)(m,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      }
    } else { // A^T X = B, or A^H X = B <==> A^T conj(X) = conj(B)
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<float> *Xij=X.addr();
        for (int ij=0;ij<m*n;ij++,Xij++) *Xij=conj(*Xij);
      }
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(cgtsv)(m,n,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(cgtsvnp)(m,n,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      }
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<float> *Xij=X.addr();
        for (int ij=0;ij<m*n;ij++,Xij++) *Xij=conj(*Xij);
      }
    }
  } else {
    CHECK_SAME(n,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) { // X A = B ==> A^T X^T = B^T
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(cgtsvr)(n,m,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(cgtsvrnp)(n,m,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      }
    } else {
      // X A^T = B ==> A X^T = B^T
      // X A^H = B ==> A X^H = B^H
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<float> *Xij=X.addr();
        for (int ij=0;ij<m*n;ij++,Xij++) *Xij=conj(*Xij);
      }
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(cgtsvr)(n,m,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(cgtsvrnp)(n,m,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      }
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<float> *Xij=X.addr();
        for (int ij=0;ij<m*n;ij++,Xij++) *Xij=conj(*Xij);
      }
    }
  }
  CHECK_TEST(info==0);
  delete LU; LU=0;
}

template<> float GaussianFactorization<float,complex<float>,
TridiagonalMatrix<float,complex<float> > >::reciprocalConditionNumber(
Factorization::CONDITION_NUMBER_NORM cnn) {
  LU=OPERATOR_NEW
    TridiagonalMatrix<float,complex<float> >(A_original->size(0));
  LU->copy(*A_original);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  float anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  int info=0;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<float,complex<float> > *u2=
      OPERATOR_NEW Vector<float,complex<float> >(n-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(cgttrf)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    F77NAME(cgtcon)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,info);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(cgttrfnp)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    F77NAME(cgtconnp)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,info);
  }
  CHECK_SAME(info,0)
  delete LU; LU=0;
  delete [] work; work=0;
  return rcond;
}

template<> void GaussianFactorization<float,complex<float>,
TridiagonalMatrix<float,complex<float> > >::improve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float &berr,float &ferr,
Factorization::TRANSPOSE_OPTION to) {
  LU=OPERATOR_NEW
    TridiagonalMatrix<float,complex<float> >(A_original->size(0));
  LU->copy(*A_original);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  char norm=(to==Factorization::NO_TRANSPOSE ? 'O' : 'I');
  float anorm=(to==Factorization::NO_TRANSPOSE ? anormo : anormi);
  char trans='N';
  if (to==Factorization::TRANSPOSE) trans='T';
  else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info;
  float rcond;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<float,complex<float> > *u2=
      OPERATOR_NEW Vector<float,complex<float> >(n-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(cgttrf)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    CHECK_TEST(info==0);
    F77NAME(cgtcon)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,info);
    CHECK_TEST(info==0);
    F77NAME(cgtrfs)(trans,n,1,A_original->lowerDiagonalAddr(),
      A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
      LU->lowerDiagonalAddr(),LU->diagonalAddr(),LU->upperDiagonalAddr(),
      u2->addr(),ipiv,b.addr(),n,x.addr(),n,&ferr,&berr,work,rwork,info);
    CHECK_TEST(info==0);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(cgttrfnp)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    CHECK_TEST(info==0);
    F77NAME(cgtconnp)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,info);
    CHECK_TEST(info==0);
    F77NAME(cgtrfsnp)(trans,n,1,A_original->lowerDiagonalAddr(),
      A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
      LU->lowerDiagonalAddr(),LU->diagonalAddr(),LU->upperDiagonalAddr(),
      b.addr(),n,x.addr(),n,&ferr,&berr,work,rwork,info);
    CHECK_TEST(info==0);
  }
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  delete LU; LU=0;
}

template<> void GaussianFactorization<float,complex<float>,
TridiagonalMatrix<float,complex<float> > >::improve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Vector<float,float> &berr,
Vector<float,float> &ferr,Factorization::TRANSPOSE_OPTION to,
Factorization::SIDE_OPTION so) {
  LU=OPERATOR_NEW
    TridiagonalMatrix<float,complex<float> >(A_original->size(0));
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
  float anorm=(to==Factorization::NO_TRANSPOSE ? anormo : anormi);
  char trans='N';
  if (to==Factorization::TRANSPOSE) trans='T';
  else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*k);
  float *rwork=OPERATOR_NEW_BRACKET(float,k);
  int info;
  float rcond;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<float,complex<float> > *u2=
      OPERATOR_NEW Vector<float,complex<float> >(k-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,k);
    F77NAME(cgttrf)(k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    CHECK_TEST(info==0);
    F77NAME(cgtcon)(norm,k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,info);
    CHECK_TEST(info==0);
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(cgtrfs)(trans,k,n,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),u2->addr(),ipiv,B.addr(),m,X.addr(),m,
        ferr.addr(),berr.addr(),work,rwork,info);
    } else {
      F77NAME(cgtrfsr)(trans,k,m,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),u2->addr(),ipiv,B.addr(),m,X.addr(),m,
        ferr.addr(),berr.addr(),work,rwork,info);
    }
    CHECK_TEST(info==0);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(cgttrfnp)(k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    CHECK_TEST(info==0);
    F77NAME(cgtconnp)(norm,k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,info);
    CHECK_TEST(info==0);
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(cgtrfsnp)(trans,k,n,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),B.addr(),m,X.addr(),m,ferr.addr(),
        berr.addr(),work,rwork,info);
    } else {
      F77NAME(cgtrfsrnp)(trans,k,m,A_original->lowerDiagonalAddr(),
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

template<> GaussianFactorization<float,complex<float>,
BandMatrix<float,complex<float> > >::GaussianFactorization(
const BandMatrix<float,complex<float> >& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),piv_op(po),equ_op(eo),
equed('N'),colcnd(numeric_limits<float>::infinity()),
rowcnd(numeric_limits<float>::infinity()) {
  int n=A.size(0),nsub=A.subDiags(),nsup=A.supDiags();
  LU=OPERATOR_NEW BandMatrix<float,complex<float> >(n,nsub,
    (po==Factorization::NO_PIVOTING ? nsup : nsub+nsup),float_zero_);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-nsup);
    int iend=min(n-1,j+nsub);
    F77NAME(ccopy)(iend-ibeg+1,A.addr(ibeg,j),1,LU->addr(ibeg,j),1);
  }

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    r=OPERATOR_NEW Vector<float,float>(n);
    c=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    F77NAME(cgbequ)(n,n,LU->subDiags(),LU->supDiags(),LU->addr(),
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
      F77NAME(claqgb)(n,n,nsub,nsub+nsup,LU->addr(),LU->bands(),
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
      F77NAME(cgbtf2np)(n,n,nsub,nsup,LU->addr(),LU->bands(),info);
      break;
    }
    default: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(cgbtrf)(n,n,nsub,nsup,LU->addr(),LU->bands(),ipiv,info);
    }
  }
  float *work=0;
  if (info>0) { // see dgbsvx
    rpvgrw=F77NAME(clantb)('M','U','N',info,min(info-1,LU->supDiags()),
      LU->addr(),LU->bands(),work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
    cerr << "zero pivot in GaussianFactorization::GaussianFactorization"
       << "\n zero pivot number,reciprocal pivot growth = " << info
       << " " << rpvgrw << endl;
    if (LU) delete LU; LU=0;
    if (r) delete r; r=0;
    if (c) delete c; c=0;
    if (ipiv) delete [] ipiv;ipiv=0;
    A_original=0;
  } else {
    rpvgrw=F77NAME(clantb)('M','U','N',n,LU->supDiags(),LU->addr(),
      LU->bands(),work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
  }
}

template<> void GaussianFactorization<float,complex<float>,
BandMatrix<float,complex<float> > >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,Factorization::TRANSPOSE_OPTION to) {
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
      const float *ri=r->addr();
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
    if (ipiv!=0) {
      F77NAME(cgbtrs)('N',n,A_original->subDiags(),A_original->supDiags(),
        1,LU->addr(),LU->bands(),ipiv,x.addr(),n,info);
    } else {
      F77NAME(cgbtrsnp)('N',n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),x.addr(),n,info);
    }
    CHECK_TEST(info==0);
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ci=c->addr();
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
  } else {
    // A^T X = B ==> U^T L^T Q^T R^{-1} X = C B
    // A^H X = B ==> U^H L^H Q^T R^{-1} X = C B
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ci=c->addr();
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
    char trans=(to==Factorization::TRANSPOSE ? 'T' : 'C');
    if (ipiv!=0) {
      F77NAME(cgbtrs)(trans,n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
        x.addr(),n,info);
    } else {
      F77NAME(cgbtrsnp)(trans,n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),x.addr(),n,info);
    }
    CHECK_TEST(info==0);
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ri=r->addr();
      complex<float> *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
  }
}

template<> void GaussianFactorization<float,complex<float>,
BandMatrix<float,complex<float> > >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Factorization::TRANSPOSE_OPTION to,
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
          const float *ri=r->addr();
          complex<float> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
      if (ipiv!=0) {
        F77NAME(cgbtrs)('N',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),ipiv,
          X.addr(),m,info);
      } else {
        F77NAME(cgbtrsnp)('N',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),X.addr(),m,
          info);
      }
      CHECK_TEST(info==0);
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ci=c->addr();
          complex<float> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
    } else {
      // A^T X = B ==> U^T L^T Q^T R^{-1} X = C B
      // A^H X = B ==> U^H L^H Q^T R^{-1} X = C B
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ci=c->addr();
          complex<float> *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
      char trans=(to==Factorization::TRANSPOSE ? 'T' : 'C');
      if (ipiv!=0) {
        F77NAME(cgbtrs)(trans,m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),ipiv,
          X.addr(),m,info);
      } else {
        F77NAME(cgbtrsnp)(trans,m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),X.addr(),m,
          info);
      }
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ri=r->addr();
          complex<float> *Xij=X.addr(0,j);
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
        for (int j=0;j<n;j++) F77NAME(csscal)(m,(*c)[j],X.addr(0,j),1);
      }
      Vector<float,complex<float> > *x=
        OPERATOR_NEW Vector<float,complex<float> >(n);
      for (int i=0;i<m;i++) {
        F77NAME(ccopy)(n,X.addr(i,0),m,x->addr(),1);
        if (ipiv!=0) {
          F77NAME(cgbtrs)('T',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
            x->addr(),n,info);
        } else {
          F77NAME(cgbtrsnp)('T',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),x->addr(),n,
            info);
        }
        F77NAME(ccopy)(n,x->addr(),1,X.addr(i,0),m);
      }
      delete x; x=0;
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(csscal)(m,(*r)[j],X.addr(0,j),1);
      }
    } else {
      // X A^T = B ==> A X^T = B^T
      //           ==> L U C^{-1} X^T = Q^T ( B R )^T
      // X A^H = B ==> A X^H = B^H
      //           ==> L U C^{-1} X^H = Q^T ( B R )^H
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(csscal)(m,(*r)[j],X.addr(0,j),1);
      }
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<float> *Xij=X.addr();
        for (int i=0;i<m*n;i++,Xij++) *Xij=conj(*Xij);
      }
      Vector<float,complex<float> > *x=
        OPERATOR_NEW Vector<float,complex<float> >(n);
      for (int i=0;i<m;i++) {
        F77NAME(ccopy)(n,X.addr(i,0),m,x->addr(),1);
        if (ipiv!=0) {
          F77NAME(cgbtrs)('N',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
            x->addr(),n,info);
        } else {
          F77NAME(cgbtrsnp)('N',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),x->addr(),n,
            info);
        }
        F77NAME(ccopy)(n,x->addr(),1,X.addr(i,0),m);
      }
      delete x; x=0;
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(csscal)(m,(*c)[j],X.addr(0,j),1);
      }
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<float> *Xij=X.addr();
        for (int i=0;i<m*n;i++,Xij++) *Xij=conj(*Xij);
      }
    }
  }
}

template<> float GaussianFactorization<float,complex<float>,
BandMatrix<float,complex<float> > >::reciprocalConditionNumber(
Factorization::CONDITION_NUMBER_NORM cnn) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  float anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info=0;
  if (piv_op==Factorization::PIVOT_ROWS) {
    F77NAME(cgbcon)(norm,n,A_original->subDiags(),A_original->supDiags(),
      LU->addr(),LU->bands(),ipiv,anorm,rcond,work,rwork,info);
  } else {
    F77NAME(cgbconnp)(norm,n,A_original->subDiags(),
      A_original->supDiags(),LU->addr(),LU->bands(),anorm,rcond,work,
      rwork,info);
  }
  CHECK_SAME(info,0)
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void GaussianFactorization<float,complex<float>,
BandMatrix<float,complex<float> > >::improve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float &berr,float &ferr,
Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(A_original->bands()+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > residual(n);
  Vector<float,float> work(n);
  char trans='N';
  if (to==Factorization::TRANSPOSE) trans='T';
  else if (to==Factorization::CONJUGATE_TRANSPOSE) trans='C';
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(cgbmv)(trans,n,n,A_original->subDiags(),
      A_original->supDiags(),float_mone_,A_original->addr(),
      A_original->bands(),x.addr(),1,float_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(cgbamv)(trans,n,n,A_original->subDiags(),
      A_original->supDiags(),float_one_,A_original->addr(),
      A_original->bands(),x.addr(),1,float_one_,work.addr(),1);

    berr=float_zero_;
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
    F77NAME(caxpy)(n,float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  Vector<float,complex<float> > v(n);
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(clacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual,to);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual,to);
    }
  }
  int i=F77NAME(icamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void GaussianFactorization<float,complex<float>,
BandMatrix<float,complex<float> > >::improve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Vector<float,float> &berr,
Vector<float,float> &ferr,Factorization::TRANSPOSE_OPTION to,
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
  float nz=static_cast<float>(A_original->bands()+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > x(k);
  Vector<float,complex<float> > rhs(k);
  Vector<float,complex<float> > residual(k);
  Vector<float,float> work(k);
  Vector<float,complex<float> > v(k);
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
      F77NAME(ccopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(ccopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        const complex<float> *Xji=X.addr(j,0);
        complex<float> *xi=x.addr();
        const complex<float> *Bji=B.addr(j,0);
        complex<float> *rhsi=rhs.addr();
        for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
          *xi=conj(*Xji);
          *rhsi=conj(*Bji);
        }
      } else {
        F77NAME(ccopy)(k,X.addr(j,0),m,x.addr(),1);
        F77NAME(ccopy)(k,B.addr(j,0),m,rhs.addr(),1);
      }
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(cgbmv)(trans,k,k,A_original->subDiags(),
        A_original->supDiags(),float_mone_,
        A_original->addr(),A_original->bands(),x.addr(),1,float_one_,
        residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77NAME(cgbamv)(trans,k,k,A_original->subDiags(),
        A_original->supDiags(),float_one_,A_original->addr(),
        A_original->bands(),x.addr(),1,float_one_,work.addr(),1);
      float &s=berr[j];
      s=float_zero_;
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
      F77NAME(caxpy)(k,complex_float_one_,residual.addr(),1,x.addr(),1);
      lstres=s;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(ccopy)(k,x.addr(),1,X.addr(0,j),1);
      } else F77NAME(ccopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      if (to==Factorization::CONJUGATE_TRANSPOSE) {
        complex<float> *Xji=X.addr(j,0);
        const complex<float> *xi=x.addr();
        for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
      } else {
        F77NAME(ccopy)(k,x.addr(),1,X.addr(j,0),m);
      }
    }

    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(clacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual,lto);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual,lto);
      }
    }
    int i=F77NAME(icamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template class GaussianFactorization<float,complex<float>,
  SquareMatrix<float,complex<float> > >;
template void testGaussianFactorization(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "CholeskyFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> CholeskyFactorization<float,complex<float>,
SymmetricPositiveMatrix<float,complex<float> > >::CholeskyFactorization(
const SymmetricPositiveMatrix<float,complex<float> >& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equed('N'),scond(numeric_limits<float>::infinity()) {
  int n=A.size(0);
  L=OPERATOR_NEW SymmetricPositiveMatrix<float,complex<float> >(n);
  L->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,3*n);
    int info;
    F77NAME(cheequb)('L',n,L->addr(),n,s->addr(),scond,amax,work,info);
    delete [] work; work=0;
    equed='N';
    F77NAME(claqhe)('L',n,L->addr(),n,s->addr(),scond,amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=L->normInfinity();
  anormm=L->normMaxEntry();
  anormo=L->normOne();

  int info=0;
  F77NAME(cpotrf)('L',n,L->addr(),n,info);
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,4*n);
  if (info!=0) { // see dposvxx
    rpvgrw=F77_NAME(cla_porpvgrw)('L',info,A_original->addr(),n,
      L->addr(),n, work);
    delete L; L=0;
    if (s!=0) delete s; s=0;
  } else {
    rpvgrw=F77_NAME(cla_porpvgrw)('L',n,A_original->addr(),n,
      L->addr(),n, work);
  }
  delete [] work; work=0;
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveMatrix<float,complex<float> > >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x) {
  CHECK_POINTER(L);
//constructor factored S A S = M D M^T
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,L->size(0));
  if (&x!=&b) x.copy(b);
//A x = b ==> L L^T S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    complex<float> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  int info;
  F77NAME(cpotrs)('L',n,1,L->addr(),n,x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    complex<float> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveMatrix<float,complex<float> > >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Factorization::SIDE_OPTION so) {
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
        const float *si=s->addr();
        complex<float> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    int info;
    F77NAME(cpotrs)('L',m,n,L->addr(),m,X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const float *si=s->addr();
        complex<float> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,L->size(0));
//  X A = B ==> X S^{-1} L L^T = B S
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(csscal)(m,(*s)[j],X.addr(0,j),1);
    }
    F77NAME(ctrsm)('R','L','C','N',m,n,complex_float_one_,L->addr(),n,
      X.addr(),m);
    F77NAME(ctrsm)('R','L','N','N',m,n,complex_float_one_,L->addr(),n,
      X.addr(),m);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(csscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> float CholeskyFactorization<float,complex<float>,
SymmetricPositiveMatrix<float,complex<float> > >::
reciprocalConditionNumber() {
  CHECK_POINTER(L);
  int n=L->size(0);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info=0;
  F77NAME(cpocon)('L',n,L->addr(),n,anormo,rcond,work,rwork,info);
  CHECK_SAME(info,0)
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveMatrix<float,complex<float> > >::improve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float &berr,float &ferr) {
  CHECK_POINTER(L);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(n+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > residual(n);
  Vector<float,float> work(n);
  Vector<float,complex<float> > v(n);
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(chemv)('L',n,complex_float_mone_,A_original->addr(),n,
      x.addr(),1,complex_float_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77_NAME(cla_heamv)(F77NAME(ilauplo)('L'),n,float_one_,
      A_original->addr(),n,x.addr(),1,float_one_,work.addr(),1);

    berr=float_zero_;
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
    F77NAME(caxpy)(n,complex_float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  complex<float> *residuali=residual.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(clacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(icamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isave; isave=0;
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveMatrix<float,complex<float> > >::improve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Vector<float,float> &berr,
Vector<float,float> &ferr,Factorization::SIDE_OPTION so) {
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
  float nz=static_cast<float>(k+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > x(k);
  Vector<float,complex<float> > rhs(k);
  Vector<float,complex<float> > residual(k);
  Vector<float,float> work(k);
  Vector<float,complex<float> > v(k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(ccopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(ccopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      const complex<float> *Xji=X.addr(j,0);
      complex<float> *xi=x.addr();
      const complex<float> *Bji=B.addr(j,0);
      complex<float> *rhsi=rhs.addr();
      for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
        *xi=conj(*Xji);
        *rhsi=conj(*Bji);
      }
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(chemv)('L',k,complex_float_mone_,A_original->addr(),k,
        x.addr(),1,complex_float_one_,residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77_NAME(cla_heamv)(F77NAME(ilauplo)('L'),k,float_one_,
        A_original->addr(),k,x.addr(),1,float_one_,work.addr(),1);
      float &berrj=berr[j];
      berrj=float_zero_;
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
      F77NAME(caxpy)(k,complex_float_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(ccopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      complex<float> *Xji=X.addr(j,0);
      const complex<float> *xi=x.addr();
      for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
    }

    complex<float> *residuali=residual.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(clacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(icamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr[j]/=lstres;
  }
  delete [] isave; isave=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> CholeskyFactorization<float,complex<float>,
SymmetricPositiveTridiagonalMatrix<float,complex<float> > >::
CholeskyFactorization(
const SymmetricPositiveTridiagonalMatrix<float,complex<float> >& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equ_op(Factorization::NO_EQUILIBRATION),equed('N'),
scond(numeric_limits<float>::infinity()) {
  int n=A.size(0);
  anormi=A.normInfinity();
  anormm=A.normMaxEntry();
  anormo=A.normOne();
  rpvgrw=float_one_; // not computed
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveTridiagonalMatrix<float,complex<float> > >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x) {
//constructor did not factor
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<float,complex<float> >(n);
  L->copy(*A_original);
  if (&x!=&b) x.copy(b);
  int info;
  F77NAME(cptsv)(n,1,L->diagonalAddr(0),L->lowerDiagonalAddr(0),
    x.addr(),n,info);
  CHECK_TEST(info==0);
  delete L; L=0;
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveTridiagonalMatrix<float,complex<float> > >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Factorization::SIDE_OPTION so) {
//constructor did not factor
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<float,complex<float> >(
    A_original->size(0));
  L->copy(*A_original);
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,L->size(0));
    F77NAME(cptsv)(m,n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),
      X.addr(),m,info);
  } else {
    CHECK_SAME(n,L->size(0));
    F77NAME(cpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
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

template<> float CholeskyFactorization<float,complex<float>,
SymmetricPositiveTridiagonalMatrix<float,complex<float> > >::
reciprocalConditionNumber() {
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<float,complex<float> >(
    A_original->size(0));
  L->copy(*A_original);
  int n=L->size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,n);
  int info=0;
  F77NAME(cpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  F77NAME(cptcon)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,work,info);
  CHECK_SAME(info,0)
  delete L; L=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveTridiagonalMatrix<float,complex<float> > >::improve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float &berr,float &ferr) {
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<float,complex<float> >(
    A_original->size(0));
  L->copy(*A_original);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info;
  float rcond;
  F77NAME(cpttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  CHECK_TEST(info==0);
  F77NAME(cptcon)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,rwork,info);
  CHECK_TEST(info==0);
  F77NAME(cptrfs)('L',n,1,A_original->diagonalAddr(0),
    A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
    L->lowerDiagonalAddr(0),b.addr(),n,x.addr(),n,&ferr,&berr,
    work,rwork,info);
  CHECK_TEST(info==0);
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  delete L; L=0;
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveTridiagonalMatrix<float,complex<float> > >::improve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Vector<float,float> &berr,
Vector<float,float> &ferr,Factorization::SIDE_OPTION so) {
  L=OPERATOR_NEW
    SymmetricPositiveTridiagonalMatrix<float,complex<float> >(
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
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,k);
  float *rwork=OPERATOR_NEW_BRACKET(float,k);
  int info;
  float rcond;
  F77NAME(cpttrf)(k,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  CHECK_TEST(info==0);
  F77NAME(cptcon)(k,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,rwork,info);
  CHECK_TEST(info==0);
  if (so==Factorization::LEFT_SIDE) {
    F77NAME(cptrfs)('L',k,n,A_original->diagonalAddr(0),
      A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
      L->lowerDiagonalAddr(0),B.addr(),m,X.addr(),m,ferr.addr(),
      berr.addr(),work,rwork,info);
  } else {
    F77NAME(cptrfsr)('L',k,m,A_original->diagonalAddr(0),
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

template<> CholeskyFactorization<float,complex<float>,
SymmetricPositiveBandMatrix<float,complex<float> > >::
CholeskyFactorization(
const SymmetricPositiveBandMatrix<float,complex<float> >& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equ_op(eo),equed('N'),scond(numeric_limits<float>::infinity()) {
  int n=A.size(0),nsub=A.subDiags();
  L=OPERATOR_NEW SymmetricPositiveBandMatrix<float,complex<float> >(n,
    nsub,complex_float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(ccopy)(min(n-j,nsub+1),A.addr(j,j),1,L->addr(j,j),1);
  }

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    F77NAME(cpbequ)('L',n,L->subDiags(),L->addr(),L->bands(),s->addr(),
      scond,amax,info);
    equed='N';
    F77NAME(claqsb)('L',n,nsub,L->addr(),L->bands(),s->addr(),scond,
      amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=L->normInfinity();
  anormm=L->normMaxEntry();
  anormo=L->normOne();

  int info=0;
  F77NAME(cpbtrf)('L',n,nsub,L->addr(),L->bands(),info);
  float *work=0;
  if (info>0) { // see dgbsvx
    rpvgrw=F77NAME(clansb)('M','L',info,min(info-1,L->subDiags()),
      L->addr(),L->bands(),work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
    cerr << "zero pivot in CholeskyFactorization::CholeskyFactorization"
       << "\n zero pivot number,reciprocal pivot growth = " << info
       << " " << rpvgrw << endl;
    if (L) delete L; L=0;
    if (s) delete s; s=0;
    A_original=0;
  } else {
    rpvgrw=F77NAME(clansb)('M','L',n,L->subDiags(),L->addr(),
      L->bands(),work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
  }
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveBandMatrix<float,complex<float> > >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^H
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  if (&x!=&b) x.copy(b);
  int info;
//A x = b ==> L L^H S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    complex<float> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  {
  F77NAME(cpbtrs)('L',n,L->subDiags(),1,L->addr(),L->bands(),
    x.addr(),n,info);
  }
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    complex<float> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveBandMatrix<float,complex<float> > >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Factorization::SIDE_OPTION so) {
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
        const float *si=s->addr();
        complex<float> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    F77NAME(cpbtrs)('L',m,L->subDiags(),n,L->addr(),L->bands(),
      X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const float *si=s->addr();
        complex<float> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,L->size(0));
//  X A = B ==> X S^{-1} L L^H = B S
//          ==> L L^H S^{-1} X^H = S B^H
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(csscal)(m,(*s)[j],X.addr(0,j),1);
    }
    Vector<float,complex<float> > *x=
      OPERATOR_NEW Vector<float,complex<float> >(n);
    for (int i=0;i<m;i++) {
      complex<float> *Xij=X.addr(i,0);
      complex<float> *xj=x->addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *xj=conj(*Xij);
      F77NAME(cpbtrs)('L',n,L->subDiags(),1,L->addr(),L->bands(),
        x->addr(),n,info);
      xj=x->addr();
      Xij=X.addr(i,0);
      for (int j=0;j<n;j++,Xij+=m,xj++) *Xij=conj(*xj);
    }
    delete x; x=0;
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(csscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> float CholeskyFactorization<float,complex<float>,
SymmetricPositiveBandMatrix<float,complex<float> > >::
reciprocalConditionNumber() {
  CHECK_POINTER(L);
  int n=L->size(0);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  float *rwork=OPERATOR_NEW_BRACKET(float,n);
  int info=0;
  F77NAME(cpbcon)('L',n,L->subDiags(),L->addr(),L->bands(),anormo,
    rcond,work,rwork,info);
  CHECK_SAME(info,0)
  delete [] rwork; rwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveBandMatrix<float,complex<float> > >::improve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float &berr,float &ferr) {
  CHECK_POINTER(L);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(A_original->bands()+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > residual(n);
  Vector<float,float> work(n);
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(chbmv)('L',n,A_original->subDiags(),complex_float_mone_,
      A_original->addr(),A_original->bands(),x.addr(),1,
      complex_float_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(chbamv)('L',n,A_original->subDiags(),float_one_,
      A_original->addr(),A_original->bands(),x.addr(),1,
      float_one_,work.addr(),1);

    berr=float_zero_;
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
    F77NAME(caxpy)(n,complex_float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  Vector<float,complex<float> > v(n);
  int kase=0;
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(clacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(icamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isave; isave=0;
}

template<> void CholeskyFactorization<float,complex<float>,
SymmetricPositiveBandMatrix<float,complex<float> > >::improve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,
Vector<float,float> &berr,
Vector<float,float> &ferr,Factorization::SIDE_OPTION so) {
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
  float nz=static_cast<float>(A_original->bands()+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > x(k);
  Vector<float,complex<float> > rhs(k);
  Vector<float,complex<float> > residual(k);
  Vector<float,float> work(k);
  Vector<float,complex<float> > v(k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(ccopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(ccopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      const complex<float> *Xji=X.addr(j,0);
      complex<float> *xi=x.addr();
      const complex<float> *Bji=B.addr(j,0);
      complex<float> *rhsi=rhs.addr();
      for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
        *xi=conj(*Xji);
        *rhsi=conj(*Bji);
      }
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(chbmv)('L',k,A_original->subDiags(),complex_float_mone_,
        A_original->addr(),A_original->bands(),x.addr(),1,
        complex_float_one_,residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77NAME(chbamv)('L',k,A_original->subDiags(),float_one_,
        A_original->addr(),A_original->bands(),x.addr(),1,
        float_one_,work.addr(),1);
      float &berrj=berr[j];
      berrj=float_zero_;
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
      F77NAME(caxpy)(k,complex_float_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(ccopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      complex<float> *Xji=X.addr(j,0);
      const complex<float> *xi=x.addr();
      for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
    }

    complex<float> *residuali=residual.addr();
    float *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(clacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(icamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr/=lstres;
  }
  delete [] isave; isave=0;
}

template class CholeskyFactorization<float,complex<float>,
  SymmetricPositiveMatrix<float,complex<float> > >;
template void testCholeskyFactorization(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "MDMtFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> MDMtFactorization<float,complex<float> >::MDMtFactorization(
const SymmetricMatrix<float,complex<float> >& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),MD(0),ipiv(0),s(0),
equed('N'),scond(numeric_limits<float>::infinity()) {
  int n=A.size(0);
  MD=OPERATOR_NEW SymmetricMatrix<float,complex<float> >(n);
  MD->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,3*n);
    F77NAME(cheequb)('L',n,MD->addr(),n,s->addr(),scond,amax,work,info);
    delete [] work; work=0;
    equed='N';
    F77NAME(claqsy)('L',n,MD->addr(),n,s->addr(),scond,amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=MD->normInfinity();
  anormm=MD->normMaxEntry();
  anormo=MD->normOne();

  int info=0;
  ipiv=OPERATOR_NEW_BRACKET(int,n);
  complex<float> workq=complex_float_zero_;
  int lwork=-1;
  F77NAME(chetrf)('L',n,MD->addr(),n,ipiv,&workq,lwork,info);
  lwork=max(2*n,static_cast<int>(workq.real()));
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(chetrf)('L',n,MD->addr(),n,ipiv,work,lwork,info);
  if (info!=0) { // see dposvxx
    rpvgrw=F77_NAME(cla_syrpvgrw)('L',n,info,A_original->addr(),n,
      MD->addr(),n,ipiv,work);
    delete MD; MD=0;
    if (s!=0) delete s; s=0;
    if (ipiv!=0) delete ipiv; ipiv=0;
  } else {
    rpvgrw=F77_NAME(cla_syrpvgrw)('L',n,n,A_original->addr(),n,
      MD->addr(),n,ipiv,work);
  }
  delete [] work; work=0;
}

template<> void MDMtFactorization<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x) {
  CHECK_POINTER(MD);
//constructor factored S A S = M D M^H
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,MD->size(0));
  if (&x!=&b) x.copy(b);
//A x = b ==> M D M^H S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    complex<float> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  int info;
  F77NAME(chetrs)('L',n,1,MD->addr(),n,ipiv,x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    complex<float> *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void MDMtFactorization<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Factorization::SIDE_OPTION so) {
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
        const float *si=s->addr();
        complex<float> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    int info;
    F77NAME(chetrs)('L',m,n,MD->addr(),m,ipiv,X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const float *si=s->addr();
        complex<float> *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,MD->size(0));
//  X A = B ==> X S^{-1} M D M^T = B S
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(csscal)(m,(*s)[j],X.addr(0,j),1);
    }
    Vector<float,complex<float> > x(n);
    for (int i=0;i<m;i++) {
      complex<float> *Xij=X.addr(i,0);
      complex<float> *xj=x.addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *xj=conj(*Xij);
      int info;
      F77NAME(chetrs)('L',n,1,MD->addr(),n,ipiv,x.addr(),n,info);
      CHECK_TEST(info==0);
      Xij=X.addr(i,0);
      xj=x.addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *Xij=conj(*xj);
    }
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(csscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> float
MDMtFactorization<float,complex<float> >::reciprocalConditionNumber() {
  CHECK_POINTER(MD);
  int n=MD->size(0);
  float rcond;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,2*n);
  int info=0;
  F77NAME(checon)('L',n,MD->addr(),n,ipiv,anormo,rcond,work,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  return rcond;
}

template<> void MDMtFactorization<float,complex<float> >::improve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float &berr,float &ferr) {
  CHECK_POINTER(MD);
  int n=MD->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(n+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > residual(n);
  Vector<float,float> work(n);
  Vector<float,complex<float> > v(n);
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(chemv)('L',n,complex_float_mone_,A_original->addr(),n,
      x.addr(),1,complex_float_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77_NAME(cla_heamv)(F77NAME(ilauplo)('L'),n,float_one_,
      A_original->addr(),n,x.addr(),1,float_one_,work.addr(),1);

    berr=float_zero_;
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
    F77NAME(caxpy)(n,complex_float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  complex<float> *residuali=residual.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(clacn2)(n,v.addr(),residual.addr(),ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(icamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isave; isave=0;
}

template<> void MDMtFactorization<float,complex<float> >::improve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Vector<float,float> &berr,
Vector<float,float> &ferr,Factorization::SIDE_OPTION so) {
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
  float nz=static_cast<float>(k+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,complex<float> > x(k);
  Vector<float,complex<float> > rhs(k);
  Vector<float,complex<float> > residual(k);
  Vector<float,float> work(k);
  Vector<float,complex<float> > v(k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(ccopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(ccopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      const complex<float> *Xji=X.addr(j,0);
      complex<float> *xi=x.addr();
      const complex<float> *Bji=B.addr(j,0);
      complex<float> *rhsi=rhs.addr();
      for (int i=0;i<k;i++,Xji+=m,Bji+=m,xi++,rhsi++) {
        *xi=conj(*Xji);
        *rhsi=conj(*Bji);
      }
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(chemv)('L',k,complex_float_mone_,A_original->addr(),k,
        x.addr(),1,complex_float_one_,residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77_NAME(cla_heamv)(F77NAME(ilauplo)('L'),k,float_one_,
        A_original->addr(),k,x.addr(),1,float_one_,work.addr(),1);
      float &berrj=berr[j];
      berrj=float_zero_;
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
      F77NAME(caxpy)(k,complex_float_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(ccopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      complex<float> *Xji=X.addr(j,0);
      const complex<float> *xi=x.addr();
      for (int i=0;i<k;i++,Xji+=m,xi++) *Xji=conj(*xi);
    }

    complex<float> *residuali=residual.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(clacn2)(k,v.addr(),residual.addr(),ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(icamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr[j]/=lstres;
  }
  delete [] isave; isave=0;
}

template class MDMtFactorization<float,complex<float> >;
template void testMDMtFactorization(float,complex<float>);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "HouseholderQRFactorization.C"
template<> HouseholderQRFactorization<float,complex<float> >
::HouseholderQRFactorization(const Matrix<float,complex<float> > &A,
Factorization::PIVOT_OPTION po) :
QR(0),tau(0),jpvt(0),piv_op(po),A_original(&A),iascl(0),ascl(float_one_),
anrm(numeric_limits<float>::infinity()) {
//TRACER_CALL(t,"HouseholderQRFactorization::HouseholderQRFactorization");
  int m=A.size(0),n=A.size(1);
  QR=OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  QR->copy(A);
  float dw=numeric_limits<float>::infinity();
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);
  anrm=F77NAME(clange)('M',m,n,A.addr(),m,&dw);
  if (anrm>float_zero_ && anrm<smlnum) {
    ascl=smlnum;
    iascl=1;
  } else if (anrm>bignum) {
    ascl=bignum;
    iascl=2;
  } else if (anrm==float_zero_) {
    delete QR; QR=0;
    return;
  }
  int info=0;
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,anrm,ascl,m,n,QR->addr(),m,info);
    CHECK_SAME(info,0);
  }

  if (piv_op==Factorization::NO_PIVOTING) {
    complex<float> zw=complex_float_undefined_;
    int lwork=-1;
    if (m>=n) {
      tau=OPERATOR_NEW Vector<float,complex<float> >(n);
      F77NAME(cgeqrf)(m,n,QR->addr(),m,tau->addr(),&zw,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cgeqrf)(m,n,QR->addr(),m,tau->addr(),work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
    } else {
      tau=OPERATOR_NEW Vector<float,complex<float> >(m);
      F77NAME(cgelqf)(m,n,QR->addr(),m,tau->addr(),&zw,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cgelqf)(m,n,QR->addr(),m,tau->addr(),work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
    }
  } else {
    int mn=min(m,n);
    int nb=F77NAME(ilaenv)(1,"CGEQRF"," ",m,n,-1,1);
    int lwkmin=mn+max(2*mn,n+1);
    int lwork=max(lwkmin,mn+2*n+nb*(n+1));
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);

    tau=OPERATOR_NEW Vector<float,complex<float> >(mn);
    jpvt=OPERATOR_NEW_BRACKET(int,n);
    float *rwork=OPERATOR_NEW_BRACKET(float,2*n);
    F77NAME(cgeqp3)(m,n,QR->addr(),m,jpvt,tau->addr(),work,lwork,
      rwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    delete [] rwork; rwork=0;
  }
}

template<> OrthogonalMatrix<float,complex<float> >* 
HouseholderQRFactorization<float,complex<float> >::orthogonalPart()
const {
//TRACER_CALL(t,"HouseholderQRFactorization::orthogonalPart");
  int m=QR->size(0),n=QR->size(1);
  complex<float> w=complex_float_undefined_;
  if (m>=n) {
    OrthogonalMatrix<float,complex<float> > *Q=
      OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(m,m);
    Q->copyFrom('A',m,n,*QR);

    int lwork=-1;
    int info;
    F77NAME(cungqr)(m,m,n,Q->addr(),m,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w.real());
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cungqr)(m,m,n,Q->addr(),m,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    return Q;
  } else {
    OrthogonalMatrix<float,complex<float> > *Q=
      OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(n,n);
    Q->copyFrom('A',m,n,*QR);

    int lwork=-1;
    int info;
    F77NAME(cunglq)(n,n,m,Q->addr(),n,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w.real());
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunglq)(n,n,m,Q->addr(),n,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    return Q;
  }
}

template<> float
HouseholderQRFactorization<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,Factorization::TRANSPOSE_OPTION tr)
const {
//TRACER_CALL(t,"HouseholderQRFactorization::solve");
  int m=QR->size(0),n=QR->size(1);
  float dw=numeric_limits<float>::infinity();
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(clange)('M',brow,1,b.addr(),b.size(),&dw);
  int ibscl=0;
  float bscl=float_one_;
  if (bnrm>float_zero_ && bnrm<smlnum) {
    ibscl=1;
    bscl=smlnum;
  } else if (bnrm > bignum) {
    ibscl=2;
    bscl=bignum;
  }

  int info;
  complex<float> zw=complex_float_undefined_;
  int lwork=-1;
  int scllen=0;
  float residual_norm=numeric_limits<float>::infinity();
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
        Vector<float,complex<float> > xtmp(m);
        xtmp.copy(b);
        if (ibscl!=0) {
          F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      Q' b = [ y \\ z ]:
        F77NAME(cunmqr)('L','C',m,1,n,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,&zw,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cunmqr)('L','C',m,1,n,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,work,lwork,info);
        CHECK_SAME(info,0)

//      solve R x = y:
        F77NAME(ctrsv)('U','N','N',n,QR->addr(),m,xtmp.addr(),1);
        delete [] work; work=0;
        x.copyFrom(n,xtmp);
        residual_norm=F77NAME(scnrm2)(m-n,xtmp.addr(n),1);
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

        x=complex_float_zero_;
        x.copyFrom(m,b);
        if (ibscl!=0) {
          F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,x.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      solve R' v = b
        F77NAME(ctrsv)('U','C','N',n,QR->addr(),m,x.addr(),1);
//      Q [v\\0]
        F77NAME(cunmqr)('L','N',m,1,n,QR->addr(),m,tau->addr(),
                        x.addr(),m,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cunmqr)('L','N',m,1,n,QR->addr(),m,tau->addr(),
                        x.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=float_zero_;
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
        x=complex_float_zero_;
        x.copyFrom(m,b);
        if (ibscl!=0) {
          F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,x.addr(),m,info);
          CHECK_SAME(info,0)
        }
        F77NAME(ctrsv)('L','N','N',m,QR->addr(),m,x.addr(),1);
        F77NAME(cunmlq)('L','C',n,1,m,QR->addr(),m,tau->addr(),
                        x.addr(),n,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cunmlq)('L','C',n,1,m,QR->addr(),m,tau->addr(),
                        x.addr(),n,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=float_zero_;
        scllen=n;
      } else { // min || b - A' x ||
//      if A=[L,0]Q, r=b-A'x and 0=Ar then
//        [v]=Qr=Qb-[L']x=[y]-[L'x] and 0=[L,0]Qr=Lv ==> v=0, w=z
//        [w]       [0 ]  [z] [ 0 ]
//      so compute [y]=Qb, solve L'x=y
//                 [z]
        CHECK_SAME(n,b.size());
        CHECK_SAME(m,x.size());
        Vector<float,complex<float> > xtmp(n);
        xtmp.copy(b);
        if (ibscl!=0) {
          F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
          CHECK_SAME(info,0)
        }
        F77NAME(cunmlq)('L','N',m,1,m,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cunmlq)('L','N',m,1,m,QR->addr(),m,tau->addr(),
                        xtmp.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        F77NAME(ctrsv)('L','C','N',m,QR->addr(),m,xtmp.addr(),1);
        x.copyFrom(m,xtmp);
        residual_norm=F77NAME(scnrm2)(n-m,xtmp.addr(m),1);
        scllen=m;
      }
    }
  } else {
    int mn=min(m,n);
    if (tr==Factorization::NO_TRANSPOSE) {// min || x || s.t. A x = b
//    TRACER_CALL(t,"HouseholderQRFactorization::solve");
      CHECK_SAME(m,b.size())
      CHECK_SAME(n,x.size())
#ifdef DEBUG
//    cout << "\tA_original = " << endl;
//    A_original->printOn(cout);
//    cout << "\tb = " << endl;
//    b.printOn(cout);
//    printOn(cout);
#endif
      Vector<float,complex<float> > xtmp(m);
      xtmp.copy(b);
      if (ibscl!=0) {
        F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
        CHECK_SAME(info,0)
      }
#ifdef DEBUG
//    cout << "\tscaled b = " << endl;
//    xtmp.printOn(cout);
#endif
//    Q' b = [ y \\ z ]:
      F77NAME(cunmqr)('L','C',m,1,mn,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,&zw,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cunmqr)('L','C',m,1,mn,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,work,lwork,info);
      CHECK_SAME(info,0)
#ifdef DEBUG
//    cout << "\tQ' b = " << endl;
//    xtmp.printOn(cout);
#endif
      delete [] work; work=0;
//    solve R x' = y:
      F77NAME(ctrsv)('U','N','N',mn,QR->addr(),m,xtmp.addr(),1);
#ifdef DEBUG
//    cout << "\tx' = " << endl;
//    xtmp.printOn(cout);
#endif
//    permute
      for (int i=0;i<n;i++) x[jpvt[i]-1]=xtmp[i];
      residual_norm=F77NAME(scnrm2)(m-mn,xtmp.addr(mn),1);
#ifdef DEBUG
//    cout << "\tx = " << endl;
//    x.printOn(cout);
#endif
    } else {
      CHECK_SAME(n,b.size())
      CHECK_SAME(m,x.size())
      Vector<float,complex<float> > xtmp(n);
//    P'b:
      for (int i=0;i<n;i++) xtmp[i]=b[jpvt[i]-1];
      if (ibscl!=0) {
        F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    solve R' v = P' b
      F77NAME(ctrsv)('U','C','N',mn,QR->addr(),m,xtmp.addr(),1);
      x=complex_float_zero_;
      x.copyFrom(mn,xtmp);
//    Q [v\\0]
      lwork=-1;
      F77NAME(cunmqr)('L','N',m,1,mn,QR->addr(),m,tau->addr(),
                      x.addr(),m,&zw,lwork,info );
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cunmqr)('L','N',m,1,mn,QR->addr(),m,tau->addr(),
                      x.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      if (mn<n) {
        residual_norm=F77NAME(scnrm2)(n-mn,xtmp.addr(mn),1);
      } else residual_norm=float_zero_;
    }
  }
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(clascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  return residual_norm;
}

template<> void
HouseholderQRFactorization<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,
Vector<float,complex<float> > &residual_norm,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::solve");
  int m=QR->size(0),n=QR->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  CHECK_SAME(k,residual_norm.size())
  float dw=numeric_limits<float>::infinity();
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(clange)('M',brow,k,B.addr(),B.size(0),&dw);
  int ibscl=0;
  float bscl=float_one_;
  if (bnrm>float_zero_ && bnrm<smlnum) {
    ibscl=1;
    bscl=smlnum;
  } else if (bnrm > bignum) {
    ibscl=2;
    bscl=bignum;
  }

  int info;
  complex<float> zw=complex_float_undefined_;
  int lwork=-1;
  int scllen=0;
  if (piv_op==Factorization::NO_PIVOTING) {
    if (m>=n) {
      if (tr==Factorization::NO_TRANSPOSE) { // min || B - A x ||
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        Matrix<float,complex<float> > Xtmp(m,k);
        Xtmp.copy(B);
        if (ibscl!=0) {
          F77NAME(clascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      Q' B = [ Y \\ Z ]:
        F77NAME(cunmqr)('L','C',m,k,n,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,&zw,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cunmqr)('L','C',m,k,n,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,work,lwork,info);
        CHECK_SAME(info,0)

//      solve R X = Y:
        F77NAME(ctrsm)('L','U','N','N',n,k,complex_float_one_,
          QR->addr(),m,Xtmp.addr(),m);
        delete [] work; work=0;
        X.copyFrom('A',n,k,Xtmp);
        for (int j=0;j<k;j++) {
          residual_norm[j]=F77NAME(scnrm2)(m-n,Xtmp.addr(n,j),1);
        }
        scllen=n;
      } else { // min || x || s.t. A^t x = B
        CHECK_SAME(n,B.size(0))
        CHECK_SAME(m,X.size(0))

        X=complex_float_zero_;
        X.copyFrom('A',m,k,B);
        if (ibscl!=0) {
          F77NAME(clascl)('G',0,0,bnrm,bscl,brow,k,X.addr(),m,info);
          CHECK_SAME(info,0)
        }
        F77NAME(ctrsm)('L','U','C','N',n,k,complex_float_one_,
          QR->addr(),m,X.addr(),m);
        F77NAME(cunmqr)('L','N',m,k,n,QR->addr(),m,tau->addr(),
                        X.addr(),m,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cunmqr)('L','N',m,k,n,QR->addr(),m,tau->addr(),
                        X.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=float_zero_;
        scllen=m;
      }
    } else {
      if (tr==Factorization::NO_TRANSPOSE) { // min || x || s.t. A x = B
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        X=complex_float_zero_;
        X.copyFrom('A',m,k,B);
        if (ibscl!=0) {
          F77NAME(clascl)('G',0,0,bnrm,bscl,brow,k,X.addr(),m,info);
          CHECK_SAME(info,0)
        }
//      solve L X = B:
        F77NAME(ctrsm)('L','L','N','N',m,k,complex_float_one_,
          QR->addr(),m,X.addr(),n);
//      Q' V
        F77NAME(cunmlq)('L','C',n,k,m,QR->addr(),m,tau->addr(),
                        X.addr(),n,&zw,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cunmlq)('L','C',n,k,m,QR->addr(),m,tau->addr(),
                        X.addr(),n,work,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        residual_norm=float_zero_;
        scllen=n;
      } else { // // min || B - A' X ||
        CHECK_SAME(n,B.size(0))
        CHECK_SAME(m,X.size(0))

        Matrix<float,complex<float> > Xtmp(n,k);
        Xtmp.copy(B);
        if (ibscl!=0) {
          F77NAME(clascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
          CHECK_SAME(info,0)
        }
        F77NAME(cunmlq)('L','N',m,k,m,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,&zw,lwork,info );
        CHECK_SAME(info,0)

        lwork=static_cast<int>(zw.real());
        complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
        F77NAME(cunmlq)('L','N',m,k,m,QR->addr(),m,tau->addr(),
                        Xtmp.addr(),m,work,lwork,info );
        CHECK_SAME(info,0)
        delete [] work; work=0;

        F77NAME(ctrsm)('L','L','C','N',m,k,complex_float_one_,
          QR->addr(),m,Xtmp.addr(),n);
        X.copyFrom('A',m,k,Xtmp);
        residual_norm=float_zero_;
        for (int j=0;j<k;j++) {
          residual_norm[j]=F77NAME(scnrm2)(n-m,Xtmp.addr(m,j),1);
        }
        scllen=m;
      }
    }
  } else {
    int mn=min(m,n);
    if (tr==Factorization::NO_TRANSPOSE) {
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<float,complex<float> > Xtmp(m,k);
      Xtmp.copy(B);
      if (ibscl!=0) {
        F77NAME(clascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    Q' B = [ Y \\ Z ]:
      F77NAME(cunmqr)('L','C',m,k,mn,QR->addr(),m,tau->addr(),
                      Xtmp.addr(),m,&zw,lwork,info);
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cunmqr)('L','C',m,k,mn,QR->addr(),m,tau->addr(),
                      Xtmp.addr(),m,work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
//    solve R V = Y:
      F77NAME(ctrsv)('U','N','N',mn,QR->addr(),m,Xtmp.addr(),1);
//    P V:
      for (int j=0;j<k;j++) {
        for (int i=0;i<n;i++) X(jpvt[i]-1,j)=Xtmp(i,j);
        residual_norm[j]=
          F77NAME(scnrm2)(m-mn,Xtmp.addr(mn,j),1);
      }
    } else {
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<float,complex<float> > Xtmp(n,k);
//    P'B:
      for (int j=0;j<k;j++) {
        for (int i=0;i<n;i++) Xtmp(i,j)=B(jpvt[i]-1,j);
      }
      if (ibscl!=0) {
        F77NAME(clascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    solve R' V = P' B
      F77NAME(ctrsv)('U','C','N',mn,QR->addr(),m,Xtmp.addr(),1);
      X=complex_float_zero_;
      X.copyFrom('A',mn,k,Xtmp);
//    Q [V\\0]
      lwork=-1;
      F77NAME(cunmqr)('L','N',m,k,mn,QR->addr(),m,tau->addr(),
                      X.addr(),m,&zw,lwork,info );
      CHECK_SAME(info,0)
      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cunmqr)('L','N',m,k,mn,QR->addr(),m,tau->addr(),
                      X.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      if (mn<n) {
        for (int j=0;j<k;j++) {
          residual_norm[j]=
            F77NAME(scnrm2)(n-mn,Xtmp.addr(mn,j),1);
        }
      } else residual_norm=float_zero_;
    }
  }
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(clascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
}

template<> void
HouseholderQRFactorization<float,complex<float> >::improve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::improve");
  CHECK_TEST(piv_op==Factorization::NO_PIVOTING)
  int m=QR->size(0),n=QR->size(1);
  float dw=numeric_limits<float>::infinity();
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(clange)('M',brow,1,b.addr(),b.size(),&dw);
  int ibscl=0;
  float bscl=float_one_;
  if (bnrm>float_zero_ && bnrm<smlnum) {
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
    F77NAME(clascl)('G',0,0,bnrm,bscl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,ascl,anrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  complex<float> zw=complex_float_undefined_;
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
      Vector<float,complex<float> > r(m);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    r=b-Ax
      F77NAME(cgemv)('N',m,n,complex_float_mone_,A_original->addr(),m,
        x.addr(),1,complex_float_one_,r.addr(),1);
      Vector<float,complex<float> > c(n,complex_float_zero_);
//    -dc=A'r
      F77NAME(cgemv)('C',m,n,complex_float_one_,A_original->addr(),m,
        r.addr(),1,complex_float_zero_,c.addr(),1);
//    solve R'(-dv) =-dc:
      F77NAME(ctrsv)('U','C','N',n,QR->addr(),m,c.addr(),1);
//    solve Rdx=-v:
      F77NAME(ctrsv)('U','N','N',n,QR->addr(),m,c.addr(),1);
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
      Vector<float,complex<float> > r(n);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    r=b-A'x
      F77NAME(cgemv)('C',m,n,complex_float_mone_,A_original->addr(),m,
        x.addr(),1,complex_float_one_,r.addr(),1);
//    solve R'dv=r:
      Vector<float,complex<float> > c(m,complex_float_zero_);
      c.copyFrom(n,r);
      F77NAME(ctrsv)('U','C','N',n,QR->addr(),m,c.addr(),1);
//    dx=Q[dv\\0]
      F77NAME(cunmqr)('L','C',m,1,n,QR->addr(),m,tau->addr(),
                      c.addr(),m,&zw,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cunmqr)('L','C',m,1,n,QR->addr(),m,tau->addr(),
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
      Vector<float,complex<float> > v(n,complex_float_zero_);
      v.copyFrom(m,b);
      if (ibscl!=0) {
        F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,v.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    r=b-Ax
      F77NAME(cgemv)('N',m,n,complex_float_mone_,A_original->addr(),m,
        x.addr(),1,complex_float_one_,v.addr(),1);
//    solve Ldv=r:
      F77NAME(ctrsv)('L','N','N',m,QR->addr(),m,v.addr(),1);
//    dx=A'[dv\\dg]:
      F77NAME(cunmlq)('L','C',n,1,m,QR->addr(),m,tau->addr(),
                      v.addr(),n,&zw,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cunmlq)('L','C',n,1,m,QR->addr(),m,tau->addr(),
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
      Vector<float,complex<float> > r(n);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    r=b-A'x
      F77NAME(cgemv)('C',m,n,complex_float_mone_,A_original->addr(),m,
        x.addr(),1,complex_float_one_,r.addr(),1);
      Vector<float,complex<float> > c(m,complex_float_zero_);
//    -dc=Ar
      F77NAME(cgemv)('N',m,n,complex_float_one_,A_original->addr(),m,
        r.addr(),1,complex_float_zero_,c.addr(),1);
//    solve L(-dv)=-dc:
      F77NAME(ctrsv)('L','N','N',m,QR->addr(),m,c.addr(),1);
//    solve L'dx=-dv:
      F77NAME(ctrsv)('L','C','N',m,QR->addr(),m,c.addr(),1);
      x+=c;
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(clascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
}

template<> void
HouseholderQRFactorization<float,complex<float> >
::improve(const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Factorization::TRANSPOSE_OPTION tr)
const {
//TRACER_CALL(t,"HouseholderQRFactorization::improve");
  CHECK_TEST(piv_op==Factorization::NO_PIVOTING)
  int m=QR->size(0),n=QR->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  float dw=numeric_limits<float>::infinity();
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);
  
  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(clange)('M',brow,k,B.addr(),B.size(0),&dw);
  int ibscl=0;
  float bscl=float_one_;
  if (bnrm>float_zero_ && bnrm<smlnum) {
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
    F77NAME(clascl)('G',0,0,bnrm,bscl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,ascl,anrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  complex<float> zw=complex_float_undefined_;
  if (m>=n) {
    if (tr==Factorization::NO_TRANSPOSE) { // min || B - A x ||
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<float,complex<float> > R(m,k);
      R.copy(B);
//    R=B-AX
      F77NAME(cgemm)('N','N',m,k,n,complex_float_mone_,
        A_original->addr(),m,X.addr(),n,complex_float_one_,R.addr(),m);
      Matrix<float,complex<float> > C(n,k,complex_float_zero_);
//    -dC=A'R
      F77NAME(cgemm)('C','N',n,k,m,complex_float_one_,
        A_original->addr(),m,R.addr(),m,complex_float_zero_,C.addr(),n);
//    solve R'(-dV)=-dC:
      F77NAME(ctrsm)('L','U','C','N',n,k,complex_float_one_,QR->addr(),m,
        C.addr(),n);
//    solve RdX=-V:
      F77NAME(ctrsm)('L','U','N','N',n,k,complex_float_one_,QR->addr(),m,
        C.addr(),n);
      X+=C;
      scllen=n;
    } else { // min || X || s.t. A^t X = B
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<float,complex<float> > R(n,k);
      R.copy(B);
//    R=B-A'X
      F77NAME(cgemm)('C','N',n,k,m,complex_float_mone_,
        A_original->addr(),m,X.addr(),m,complex_float_one_,R.addr(),n);
//    solve R'dV=R:
      Matrix<float,complex<float> > C(m,k,complex_float_zero_);
      C.copyFrom('A',n,k,R);
      F77NAME(ctrsm)('L','U','C','N',n,k,complex_float_one_,QR->addr(),m,
        C.addr(),n);
//    dX=Q[dV\\0]
      F77NAME(cunmqr)('L','C',m,k,n,QR->addr(),m,tau->addr(),
                      C.addr(),m,&zw,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cunmqr)('L','C',m,k,n,QR->addr(),m,tau->addr(),
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
      Matrix<float,complex<float> > V(n,k,complex_float_zero_);
      V.copyFrom('A',m,k,B);
//    R=B-AX
      F77NAME(cgemm)('N','N',m,k,n,complex_float_mone_,
        A_original->addr(),m,X.addr(),n,complex_float_one_,V.addr(),m);
//    solve LdV=R:
      F77NAME(ctrsm)('L','L','N','N',m,k,complex_float_one_,QR->addr(),m,
        V.addr(),n);
//    dX=Q'[dV\\0]:
      F77NAME(cunmlq)('L','C',n,k,m,QR->addr(),m,tau->addr(),
                      V.addr(),n,&zw,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(zw.real());
      complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
      F77NAME(cunmlq)('L','C',n,k,m,QR->addr(),m,tau->addr(),
                      V.addr(),n,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      X+=V;
      scllen=n;
    } else { // min || B - A' X ||
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<float,complex<float> > R(n,k);
      R.copy(B);
//    R=B-A'X
      F77NAME(cgemm)('C','N',n,k,m,complex_float_mone_,
        A_original->addr(),m,X.addr(),m,complex_float_one_,R.addr(),n);
      Matrix<float,complex<float> > C(m,k,complex_float_zero_);
//    (-dC)=A R
      F77NAME(cgemm)('N','N',m,k,n,complex_float_one_,
        A_original->addr(),m,R.addr(),n,complex_float_zero_,C.addr(),m);
//    L(-dV)=-dC
      F77NAME(ctrsm)('L','L','N','N',m,k,complex_float_one_,QR->addr(),m,
        C.addr(),m);
//    L'dX=-dV
      F77NAME(ctrsm)('L','L','C','N',m,k,complex_float_one_,QR->addr(),m,
        C.addr(),m);
      X+=C;
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(clascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
}

template class HouseholderQRFactorization<float,complex<float> >;
//template void testHouseholderQRFactorization(float,complex<float> );
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "CompleteOrthogonalDecomposition.C"
template<> CompleteOrthogonalDecomposition<float,complex<float> >::
CompleteOrthogonalDecomposition(const Matrix<float,complex<float> > &A,
float rcond) : the_rank(0),URV(0),utau(0),vtau(0),jpvt(0),iascl(0),
ascl(float_one_),anrm(numeric_limits<float>::infinity()),A_original(&A) {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::CompleteOrthogonalDecomposition");
  int m=A.size(0),n=A.size(1);
  URV=OPERATOR_NEW Matrix<float,complex<float> >(m,n);
  URV->copy(A);
  float dw=numeric_limits<float>::infinity();
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);
  anrm=F77NAME(clange)('M',m,n,A.addr(),m,&dw);
  if (anrm>float_zero_ && anrm<smlnum) {
    ascl=smlnum;
    iascl=1;
  } else if (anrm>bignum) {
    ascl=bignum;
    iascl=2;
  } else if (anrm==float_zero_) {
    delete URV; URV=0;
    return;
  }
  int info=0;
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,anrm,ascl,m,n,URV->addr(),m,info);
    CHECK_SAME(info,0);
  }

  int mn=min(m,n);
  int nb1=F77NAME(ilaenv)(1,"DGEQRF"," ",m,n,-1,1);
  int nb2=F77NAME(ilaenv)(1,"DGERQF"," ",m,n,-1,1);
  int nb=max(nb1,nb2);
  int lwkmin=mn+max(2*mn,n+1);
  int lwork=max(lwkmin,mn+2*n+nb*(n+1));
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);

  utau=OPERATOR_NEW Vector<float,complex<float> >(mn);
  jpvt=OPERATOR_NEW_BRACKET(int,n);
  float *rwork=OPERATOR_NEW_BRACKET(float,2*n);
  F77NAME(cgeqp3)(m,n,URV->addr(),m,jpvt,utau->addr(),work,lwork,
    rwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  delete [] rwork; rwork=0;
  float smax=abs((*URV)(0,0));
  if (smax<=float_zero_) return;
  float smin=smax;

  the_rank=1;
  Vector<float,complex<float> > xmin(mn);
  xmin[0]=complex_float_one_;
  Vector<float,complex<float> > xmax(mn);
  xmax[0]=complex_float_one_;
  while (the_rank<mn) {
    int i=the_rank+1;
    float sminpr=numeric_limits<float>::infinity();
    complex<float> s1=complex_float_undefined_,
      c1=complex_float_undefined_;
//                 (job,j,x,sest,w,gamma,sestpr,s,c)
    F77NAME(claic1)(2,the_rank,xmin.addr(),smin,URV->addr(0,i-1),
      (*URV)(i-1,i-1),sminpr,s1,c1);
    float smaxpr=numeric_limits<float>::infinity();
    complex<float> s2=complex_float_undefined_,
      c2=complex_float_undefined_;
    F77NAME(claic1)(1,the_rank,xmax.addr(),smax,URV->addr(0,i-1),
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
    vtau=OPERATOR_NEW Vector<float,complex<float> >(the_rank);
    complex<float> zw=complex_float_undefined_;
    int lwork=-1;
    F77NAME(ctzrzf)(the_rank,n,URV->addr(),m,vtau->addr(),&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(ctzrzf)(the_rank,n,URV->addr(),m,vtau->addr(),work,lwork,
      info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
  }
}

template<> OrthogonalMatrix<float,complex<float> >* 
CompleteOrthogonalDecomposition<float,complex<float> >
::leftOrthogonalPart() const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::leftOrthogonalPart");
  int m=URV->size(0),n=URV->size(1);
  OrthogonalMatrix<float,complex<float> > *Q=
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(m,m);
  Q->copyFrom('A',m,the_rank,*URV);

  complex<float> w=complex_float_undefined_;
  int lwork=-1;
  int info;
  F77NAME(cungqr)(m,m,the_rank,Q->addr(),m,utau->addr(),&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(cungqr)(m,m,the_rank,Q->addr(),m,utau->addr(),work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  return Q;
}

template<> OrthogonalMatrix<float,complex<float> >* 
CompleteOrthogonalDecomposition<float,complex<float> >
::rightOrthogonalPartTransposed() const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::rightOrthogonalPartTransposed");
  int m=URV->size(0),n=URV->size(1);
  OrthogonalMatrix<float,complex<float> > *Q=
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(n,n);

  int lwork=-1;
  int info;
  complex<float> zw=complex_float_undefined_;
  F77NAME(cunmrz)('L','C',n,n,the_rank,n-the_rank,URV->addr(),m,
    vtau->addr(),Q->addr(),n,&zw,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(zw.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(cunmrz)('L','C',n,n,the_rank,n-the_rank,URV->addr(),m,
    vtau->addr(),Q->addr(),n,work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  return Q;
}

template<> float
CompleteOrthogonalDecomposition<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,Factorization::TRANSPOSE_OPTION tr)
const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve");
  int m=URV->size(0),n=URV->size(1);
  float dw=numeric_limits<float>::infinity();
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(clange)('M',brow,1,b.addr(),b.size(),&dw);
  int ibscl=0;
  float bscl=float_one_;
  if (bnrm>float_zero_ && bnrm<smlnum) {
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
  float residual_norm=numeric_limits<float>::infinity();
  complex<float> zw=complex_float_undefined_;
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
    Vector<float,complex<float> > xtmp(m);
    xtmp.copy(b);
    if (ibscl!=0) {
      F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
      CHECK_SAME(info,0)
    }
//  U' b = [ y \\ z ]:
    F77NAME(cunmqr)('L','C',m,1,mn,URV->addr(),m,utau->addr(),
                    xtmp.addr(),m,&zw,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(zw.real());
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunmqr)('L','C',m,1,mn,URV->addr(),m,utau->addr(),
                    xtmp.addr(),m,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R v = y:
    F77NAME(ctrsv)('U','N','N',the_rank,URV->addr(),m,xtmp.addr(),1);
    Vector<float,complex<float> > px(n,complex_float_zero_);
    px.copyFrom(the_rank,xtmp);
    lwork=-1;
//  P'x=V[v\\0]:
    F77NAME(cunmrz)('L','C',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),px.addr(),n,&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunmrz)('L','C',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),px.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  permute
    for (int i=0;i<n;i++) x[jpvt[i]-1]=px[i];
    residual_norm=F77NAME(scnrm2)(m-the_rank,xtmp.addr(the_rank),1);
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
    Vector<float,complex<float> > xtmp(n);
//  P'b:
    for (int i=0;i<n;i++) xtmp[i]=b[jpvt[i]-1];
    if (ibscl!=0) {
      F77NAME(clascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
      CHECK_SAME(info,0)
    }
    lwork=-1;
//  [y\\z]=V'(P'b):
    F77NAME(cunmrz)('L','N',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),xtmp.addr(),n,&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunmrz)('L','N',n,1,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),xtmp.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R' v = b
    F77NAME(ctrsv)('U','C','N',the_rank,URV->addr(),m,xtmp.addr(),1);

    x=complex_float_zero_;
    x.copyFrom(the_rank,xtmp);
//  U [v\\0]
    lwork=-1;
    F77NAME(cunmqr)('L','N',m,1,mn,URV->addr(),m,utau->addr(),
                    x.addr(),m,&zw,lwork,info );
    CHECK_SAME(info,0)

    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunmqr)('L','N',m,1,mn,URV->addr(),m,utau->addr(),
                    x.addr(),m,work,lwork,info );
    CHECK_SAME(info,0)
    delete [] work; work=0;
    residual_norm=F77NAME(scnrm2)(n-the_rank,xtmp.addr(the_rank),1);
  }
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(clascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
  return residual_norm;
}

template<> void
CompleteOrthogonalDecomposition<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,Vector<float,float> &residual_norm,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve");
  int m=URV->size(0),n=URV->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  CHECK_SAME(k,residual_norm.size())
  float dw=numeric_limits<float>::infinity();
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(clange)('M',brow,k,B.addr(),B.size(0),&dw);
  int ibscl=0;
  float bscl=float_one_;
  if (bnrm>float_zero_ && bnrm<smlnum) {
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
  residual_norm=numeric_limits<float>::infinity();
  complex<float> zw=complex_float_undefined_;
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
    Matrix<float,complex<float> > Xtmp(m,k);
    Xtmp.copy(B);
    if (ibscl!=0) {
      F77NAME(clascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
      CHECK_SAME(info,0)
    }
//  U' B = [ Y \\ Z ]:
    F77NAME(cunmqr)('L','C',m,k,mn,URV->addr(),m,utau->addr(),
                    Xtmp.addr(),m,&zw,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(zw.real());
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunmqr)('L','C',m,k,mn,URV->addr(),m,utau->addr(),
                    Xtmp.addr(),m,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R V = Y:
    F77NAME(ctrsm)('L','U','N','N',the_rank,k,complex_float_one_,
      URV->addr(),m,Xtmp.addr(),m);
    Matrix<float,complex<float> > PX(n,k,complex_float_zero_);
    PX.copyFrom('A',the_rank,k,Xtmp);
    lwork=-1;
//  P'X=V[V\\0]:
    F77NAME(cunmrz)('L','C',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),PX.addr(),n,&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunmrz)('L','C',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),PX.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  permute
    for (int j=0;j<k;j++) {
      for (int i=0;i<n;i++) X(jpvt[i]-1,j)=PX(i,j);
      residual_norm[j]=
        F77NAME(scnrm2)(m-the_rank,Xtmp.addr(the_rank,j),1);
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
    Matrix<float,complex<float> > Xtmp(n,k);
//  P'B:
    for (int j=0;j<k;j++) {
      for (int i=0;i<n;i++) Xtmp(i,j)=B(jpvt[i]-1,j);
    }
    if (ibscl!=0) {
      F77NAME(clascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
      CHECK_SAME(info,0)
    }
    lwork=-1;
//  [Y\\Z]=V'(P'B):
    F77NAME(cunmrz)('L','N',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),Xtmp.addr(),n,&zw,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(zw.real());
    complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunmrz)('L','N',n,k,the_rank,n-the_rank,URV->addr(),m,
      vtau->addr(),Xtmp.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R' V = B
    F77NAME(ctrsm)('L','U','C','N',the_rank,k,complex_float_one_,
      URV->addr(),m,Xtmp.addr(),n);

    X=complex_float_zero_;
    X.copyFrom('A',the_rank,k,Xtmp);
//  U [V\\0]
    lwork=-1;
    F77NAME(cunmqr)('L','N',m,k,mn,URV->addr(),m,utau->addr(),
                    X.addr(),m,&zw,lwork,info );
    CHECK_SAME(info,0)

    lwork=static_cast<int>(zw.real());
    work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
    F77NAME(cunmqr)('L','N',m,k,mn,URV->addr(),m,utau->addr(),
                    X.addr(),m,work,lwork,info );
    CHECK_SAME(info,0)
    delete [] work; work=0;
    for (int j=0;j<k;j++) {
      residual_norm[j]=
        F77NAME(scnrm2)(n-the_rank,Xtmp.addr(the_rank,j),1);
    }
  }
  if (iascl!=0) {
    F77NAME(clascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(clascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
}

template class CompleteOrthogonalDecomposition<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "GramSchmidtQRFactorization.C"
template<> float
  GramSchmidtQRFactorization<float,complex<float> >::tol=
  float_one_/sqrt(2.);

template<> void
GramSchmidtQRFactorization<float,complex<float> >::reorthogonalize(
int j) {
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,R->size(0))
  CHECK_BOUNDS(j,0,n)
  for (int k=0;k<j;k++) {
    complex<float> dot=F77NAME(cdotc)(m,Q->addr(0,k),1,Q->addr(0,j),1)
      /(*R)(k,k);
    (*R)(k,j)+=dot;
    F77NAME(caxpy)(m,-dot,Q->addr(0,k),1,Q->addr(0,j),1);
  }
  (*R)(j,j)=F77NAME(cdotc)(m,Q->addr(0,j),1,Q->addr(0,j),1);
}

template<> GramSchmidtQRFactorization<float,complex<float> >
::GramSchmidtQRFactorization(const Matrix<float,complex<float> > &A) :
A_original(&A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::GramSchmidtQRFactorization");
  int m=A.size(0),n=A.size(1);
  assert(m>=n);
  Q=OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(m,n);
  R=OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(n,n);
  Q->copy(A);

  float *norm=OPERATOR_NEW_BRACKET(float,n);
  for (int j=0;j<n;j++) {
    norm[j]=real(F77NAME(cdotc)(m,Q->addr(0,j),1,Q->addr(0,j),1));
  }
  float tol2=tol*tol;
  for (int k=0;k<n;k++) {
    complex<float> length=
      F77NAME(cdotc)(m,Q->addr(0,k),1,Q->addr(0,k),1);
    float rlength=real(length);
    (*R)(k,k)=rlength;
    if (rlength<tol*norm[k]) reorthogonalize(k);
    if (rlength>float_zero_ &&k<n-1) {
      F77NAME(cgemv)('C',m,n-k-1,complex_float_one_/(*R)(k,k),
        Q->addr(0,k+1),m,Q->addr(0,k),1,complex_float_zero_,
        R->addr(k,k+1),n);
      for (int j=k+1;j<n;j++) (*R)(k,j)=conj((*R)(k,j));
      F77NAME(cgeru)(m,n-k-1,complex_float_mone_,Q->addr(0,k),1,
        R->addr(k,k+1),n,Q->addr(0,k+1),m);
    }
  }
  delete [] norm;
}

template<> void
GramSchmidtQRFactorization<float,complex<float> >::solveOverdetermined(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,
Vector<float,complex<float> > &residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,b.size());
  CHECK_SAME(n,x.size());
  CHECK_SAME(m,residual.size());
  residual.copy(b);
  for (int j=0;j<n;j++) {
    complex<float> dot=
      -F77NAME(cdotc)(m,Q->addr(0,j),1,residual.addr(),1)/(*R)(j,j);
    x[j]=-dot;
    F77NAME(caxpy)(m,dot,Q->addr(0,j),1,residual.addr(),1);
  }
  F77NAME(ctrsv)('U','N','U',n,R->addr(),n,x.addr(),1);
}

template<> void
GramSchmidtQRFactorization<float,complex<float> >::solveOverdetermined(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,
Matrix<float,complex<float> > &Residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,X.size(0));
  CHECK_SAME(k,X.size(1));
  CHECK_SAME(m,Residual.size(0));
  CHECK_SAME(k,Residual.size(1));
  Residual.copy(B);
  for (int j=0;j<n;j++) {
    F77NAME(cgemv)('C',m,k,complex_float_one_/(*R)(j,j),
      Residual.addr(),m,Q->addr(0,j),1,
      complex_float_zero_,X.addr(j,0),n);
    for (int kk=0;kk<k;kk++) X(j,kk)=conj(X(j,kk));
    F77NAME(cgeru)(m,k,complex_float_mone_,Q->addr(0,j),1,X.addr(j,0),n,
      Residual.addr(),m);
  }
  F77NAME(ctrsm)('L','U','N','U',n,k,complex_float_one_,R->addr(),n,
    X.addr(),n);
}

template<> void
GramSchmidtQRFactorization<float,complex<float> >::solveUnderdetermined(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveUnderdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,b.size());
  CHECK_SAME(m,x.size());

  Vector<float,complex<float> > *v=
    OPERATOR_NEW Vector<float,complex<float> >(n);
  v->copy(b);

  F77NAME(ctrsv)('U','C','U',n,R->addr(),n,v->addr(),1);
  x=complex_float_zero_;
  for (int j=0;j<n;j++) {
    complex<float> omega=((*v)[j]
      -F77NAME(cdotc)(m,Q->addr(0,j),1,x.addr(),1)/(*R)(j,j))/(*R)(j,j);
    F77NAME(caxpy)(m,omega,Q->addr(0,j),1,x.addr(),1);
  }
  delete v; v=0;
}

template<> void
GramSchmidtQRFactorization<float,complex<float> >::solveUnderdetermined(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveUnderdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(n,B.size(0));
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(k,X.size(1));

  Matrix<float,complex<float> > *V=
    OPERATOR_NEW Matrix<float,complex<float> >(n,k);
  V->copy(B);
  F77NAME(ctrsm)('L','U','C','U',n,k,complex_float_one_,R->addr(),n,
    V->addr(),n);
  X=complex_float_zero_;
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,k);
  for (int j=0;j<n;j++) {
    F77NAME(cgemv)('C',m,k,complex_float_one_,X.addr(),m,Q->addr(0,j),1,
      complex_float_zero_,work,1);
    for (int kk=0;kk<k;kk++) {
      work[kk]=((*V)(j,kk)-conj(work[kk]))/(*R)(j,j);
    }
    F77NAME(cgeru)(m,k,complex_float_one_,Q->addr(0,j),1,work,1,
      X.addr(),m);
  }
  delete [] work; work=0;
  delete V; V=0;
}

template<> void GramSchmidtQRFactorization<float,complex<float> >
::improveOverdetermined(const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,
Vector<float,complex<float> > &residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveOverdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,b.size());
  CHECK_SAME(n,x.size());
  CHECK_SAME(m,residual.size());
  Vector<float,complex<float> > *delta_b=
    OPERATOR_NEW Vector<float,complex<float> >(m);
  Vector<float,complex<float> > *delta_z=
    OPERATOR_NEW Vector<float,complex<float> >(n,0.);
  Vector<float,complex<float> > *delta_y=
    OPERATOR_NEW Vector<float,complex<float> >(n,0.);
  F77NAME(cgemv)('N',m,n,complex_float_mone_,A_original->addr(),m,
    x.addr(),1,complex_float_zero_,delta_b->addr(),1);
  F77NAME(cgemv)('C',m,n,complex_float_mone_,A_original->addr(),m,
    residual.addr(),1,complex_float_zero_,delta_z->addr(),1);
  for (int i=0;i<m;i++) (*delta_b)[i]+=b[i]-residual[i];
  F77NAME(ctrsv)('U','C','U',n,R->addr(),n,delta_z->addr(),1);
  for (int j=0;j<n;j++) {
    (*delta_y)[j]=(F77NAME(cdotc)(m,Q->addr(0,j),1,delta_b->addr(),1)
                  -(*delta_z)[j])/(*R)(j,j);
    F77NAME(caxpy)(m,-(*delta_y)[j],Q->addr(0,j),1,delta_b->addr(),1);
  }
  F77NAME(ctrsv)('U','N','U',n,R->addr(),n,delta_y->addr(),1);
  x+=(*delta_y);
  residual+=(*delta_b);
  delete delta_y; delta_y=0;
  delete delta_z; delta_z=0;
  delete delta_b; delta_b=0;
}

template<> void GramSchmidtQRFactorization<float,complex<float> >
::improveOverdetermined(const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,
Matrix<float,complex<float> > &Residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,X.size(0));
  CHECK_SAME(k,X.size(1));
  CHECK_SAME(m,Residual.size(0));
  CHECK_SAME(k,Residual.size(1));
  Matrix<float,complex<float> > *delta_B=
    OPERATOR_NEW Matrix<float,complex<float> >(m,k);
  Matrix<float,complex<float> > *delta_Z=
    OPERATOR_NEW Matrix<float,complex<float> >(n,k,0.);
  Matrix<float,complex<float> > *delta_Y=
    OPERATOR_NEW Matrix<float,complex<float> >(n,k,0.);
  F77NAME(cgemm)('N','N',m,k,n,complex_float_mone_,A_original->addr(),m,
    X.addr(),n,complex_float_zero_,delta_B->addr(),m);
  F77NAME(cgemm)('C','N',n,k,m,complex_float_mone_,A_original->addr(),m,
    Residual.addr(),m,complex_float_zero_,delta_Z->addr(),n);
  for (int j=0;j<k;j++) {
    for (int i=0;i<m;i++) (*delta_B)(i,j)+=B(i,j)-Residual(i,j);
  }
  F77NAME(ctrsm)('L','U','C','U',n,k,complex_float_one_,R->addr(),n,
    delta_Z->addr(),n);
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,k);
  for (int j=0;j<n;j++) {
    F77NAME(cgemv)('C',m,k,complex_float_one_,delta_B->addr(),m,
      Q->addr(0,j),1,complex_float_zero_,work,1);
    for (int kk=0;kk<k;kk++) {
      (*delta_Y)(j,kk)=(conj(work[kk])-(*delta_Z)(j,kk))/(*R)(j,j);
    }
    F77NAME(cgeru)(m,k,complex_float_mone_,Q->addr(0,j),1,
      delta_Y->addr(j,0),n,delta_B->addr(),m);
  }
  delete [] work; work=0;
  F77NAME(ctrsm)('L','U','N','U',n,k,complex_float_one_,R->addr(),n,
    delta_Y->addr(),n);
  X+=(*delta_Y);
  Residual+=(*delta_B);
  delete delta_Y; delta_Y=0;
  delete delta_Z; delta_Z=0;
  delete delta_B; delta_B=0;
}

template<> void GramSchmidtQRFactorization<float,complex<float> >
::improveUnderdetermined(const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveUnderdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,b.size());
  CHECK_SAME(m,x.size());

  Vector<float,complex<float> > *y=
    OPERATOR_NEW Vector<float,complex<float> >(n);
  Vector<float,complex<float> > *delta_x=
    OPERATOR_NEW Vector<float,complex<float> >(m);
  Vector<float,complex<float> > *delta_by=
    OPERATOR_NEW Vector<float,complex<float> >(n);
  y->copy(b);
  delta_x->copy(x);
  delta_by->copy(b);
  F77NAME(ctrsv)('U','C','U',n,R->addr(),n,y->addr(),1);
  for (int j=0;j<n;j++) (*y)[j]/=(*R)(j,j);
  F77NAME(cgemv)('N',m,n,complex_float_one_,Q->addr(),m,y->addr(),1,
    complex_float_mone_,delta_x->addr(),1);
  F77NAME(cgemv)('C',m,n,complex_float_mone_,A_original->addr(),m,
    x.addr(),1,complex_float_one_,delta_by->addr(),1);
  F77NAME(ctrsv)('U','C','U',n,R->addr(),n,delta_by->addr(),1);
  for (int j=0;j<n;j++) {
    complex<float> omega=
      (F77NAME(cdotc)(m,Q->addr(0,j),1,delta_x->addr(),1)-(*delta_by)[j])
      /(*R)(j,j);
    F77NAME(caxpy)(m,-omega,Q->addr(0,j),1,delta_x->addr(),1);
  }
  x+=(*delta_x);
  delete delta_by; delta_by=0;
  delete delta_x; delta_x=0;
  delete y; y=0;
}

template<> void GramSchmidtQRFactorization<float,complex<float> >
::improveUnderdetermined(const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveUnderdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(n,B.size(0));
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(k,X.size(1));

  Matrix<float,complex<float> > *Y=
    OPERATOR_NEW Matrix<float,complex<float> >(n,k);
  Matrix<float,complex<float> > *delta_X=
    OPERATOR_NEW Matrix<float,complex<float> >(m,k);
  Matrix<float,complex<float> > *delta_BY=
    OPERATOR_NEW Matrix<float,complex<float> >(n,k);
  Y->copy(B);
  delta_X->copy(X);
  delta_BY->copy(B);

  F77NAME(ctrsm)('L','U','C','U',n,k,complex_float_one_,R->addr(),n,
    Y->addr(),n);
  for (int j=0;j<n;j++) {
    F77NAME(csscal)(k,float_one_/real((*R)(j,j)),Y->addr(j,0),n);
  }
  F77NAME(cgemm)('N','N',m,k,n,complex_float_one_,Q->addr(),m,
    Y->addr(),n,complex_float_mone_,delta_X->addr(),m);
  F77NAME(cgemm)('C','N',n,k,m,complex_float_mone_,A_original->addr(),m,
    X.addr(),m,complex_float_one_,delta_BY->addr(),n);
  F77NAME(ctrsm)('L','U','C','U',n,k,complex_float_one_,R->addr(),n,
    delta_BY->addr(),n);
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,k);
  for (int j=0;j<n;j++) {
    F77NAME(cgemv)('C',m,k,complex_float_one_,delta_X->addr(),m,
      Q->addr(0,j),1,complex_float_zero_,work,1);
    for (int kk=0;kk<k;kk++) {
      work[kk]=(conj(work[kk])-(*delta_BY)(j,kk))/(*R)(j,j);
    }
    F77NAME(cgeru)(m,k,complex_float_mone_,Q->addr(0,j),1,work,1,
      delta_X->addr(),m);
  }
  delete [] work; work=0;
  X+=(*delta_X);
  delete delta_BY; delta_BY=0;
  delete delta_X; delta_X=0;
  delete Y; Y=0;
}

template<> void
GramSchmidtQRFactorization<float,complex<float> >::dropColumn(int j) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropColumn");
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(j,0,n)
  for (int k=j+1;k<n;k++) {
    float c;
    complex<float> s;
    F77NAME(crotg)(R->operator()(k,k),R->operator()(j,k),c,s);
    int ncols=n-k-1;
    if (ncols>0) {
      F77NAME(crot)(ncols,R->addr(k,k+1),n,R->addr(j,k+1),n,c,s);
    }
    F77NAME(crot)(m,Q->addr(0,k),1,Q->addr(0,j),1,c,s);
  }

  for (int k=j+1;k<n;k++) {
    for (int i=0;i<j;i++) (*R)(i,k-1)=(*R)(i,k);
    for (int i=j;i<k;i++) (*R)(i,k-1)=(*R)(i+1,k);
    F77NAME(ccopy)(m,Q->addr(0,k),1,Q->addr(0,k-1),1);
  }

  OrthogonalMatrix<float,complex<float> > *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(m,n-1);
  new_Q->copyFrom('A',m,n-1,*Q);
  delete Q;
  Q=new_Q;
  UpperTrapezoidalMatrix<float,complex<float> > *new_R=OPERATOR_NEW
    UpperTrapezoidalMatrix<float,complex<float> >(n-1,n-1);
  new_R->copyFrom(n-1,n-1,*R);
  delete R;
  R=new_R;
}

template<> void
GramSchmidtQRFactorization<float,complex<float> >::addColumn(int j,
const Matrix<float,complex<float> > &A) {
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,A.size(0))

  OrthogonalMatrix<float,complex<float> > *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(m,n+1);
  new_Q->copyFrom('A',m,n,*Q);
  delete Q;
  Q=new_Q;
  F77NAME(ccopy)(m,A.addr(0,j),1,Q->addr(0,n),1);

  UpperTrapezoidalMatrix<float,complex<float> > *new_R=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,complex<float> >(n+1,n+1);
  new_R->copyFrom(n,n,*R);
  delete R;
  R=new_R;

  float norm=F77NAME(scnrm2)(m,Q->addr(0,n),1);
  for (j=0;j<n;j++) {
    complex<float>  dot=
      -F77NAME(cdotc)(m,Q->addr(0,j),1,Q->addr(0,n),1);
    (*R)(j,n)=-dot;
    F77NAME(caxpy)(m,dot,Q->addr(0,j),1,Q->addr(0,n),1);
  }
  float nrm=F77NAME(scnrm2)(m,Q->addr(0,n),1);
  (*R)(n,n)=nrm;
  if (nrm<tol*norm) reorthogonalize(n);
  complex<float> a=complex_float_one_/(*R)(n,n);
  F77NAME(cscal)(m,a,Q->addr(0,n),1);
}

template<> void
GramSchmidtQRFactorization<float,complex<float> >::dropRow(int i) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropRow");
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(i,0,m)
  for (int k=i+1;k<m;k++) {
    F77NAME(cswap)(n,Q->addr(k-1,0),m,Q->addr(k,0),m);
  }

  complex<float> *work_column=OPERATOR_NEW_BRACKET(complex<float>,m);
  for (int i=0;i<m;i++) work_column[i]=complex_float_zero_;
  int Mmin1=m-1;
  for (int j=0;j<n;j++) {
    F77NAME(caxpy)(Mmin1,conj((*Q)(m-1,j)),Q->addr(0,j),1,work_column,1);
  }
  F77NAME(csscal)(Mmin1,float_mone_,work_column,1);

  float nrm=sqrt(float_one_-pow(F77NAME(scnrm2)(n,Q->addr(m-1,0),m),2));
  work_column[m-1]=nrm;
  float alpha=float_one_/nrm;
  F77NAME(csscal)(Mmin1,alpha,work_column,1);

  complex<float> *work_row=OPERATOR_NEW_BRACKET(complex<float>,n);
  for (int j=0;j<n;j++) work_row[j]=complex_float_zero_;
  for (int j=n-1;j>=0;j--) {
    float c;
    complex<float> s;
    F77NAME(crotg)(work_column[Mmin1],Q->operator()(Mmin1,j),c,s);
    F77NAME(crot)(Mmin1,work_column,1,Q->addr(0,j),1,c,s);
    int ncols=n-j;
    F77NAME(crot)(ncols,work_row+j,1,R->addr(j,j),n,c,s);
  }

  OrthogonalMatrix<float,complex<float> > *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(Mmin1,n);
  new_Q->copyFrom('A',Mmin1,n,*Q);
  delete Q;
  Q=new_Q;

  delete [] work_row;
  delete [] work_column;
}

template<> void
GramSchmidtQRFactorization<float,complex<float> >::addRow(int i,
const Matrix<float,complex<float> > &A) {
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,A.size(1))

  OrthogonalMatrix<float,complex<float> > *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(m+1,n);
  new_Q->copyFrom('A',m,n,*Q);
  delete Q;
  Q=new_Q;
  for (int j=0;j<n;j++) (*Q)(m,j)=complex_float_zero_;

  complex<float> *work_row=OPERATOR_NEW_BRACKET(complex<float>,n);
  int incA=A.size(0);
  F77NAME(ccopy)(n,A.addr(i,0),incA,work_row,1);

  complex<float> *work_column=OPERATOR_NEW_BRACKET(complex<float>,m+1);
  for (int k=0;k<m;k++) work_column[k]=complex_float_zero_;
  work_column[m]=complex_float_one_;

  int Mp1=m+1;
  for (int j=0;j<n;j++) {
    float c;
    complex<float> s;
    F77NAME(crotg)(R->operator()(j,j),work_row[j],c,s);
    int ncols=n-j-1;
    if (ncols>0) {
      F77NAME(crot)(ncols,R->addr(j,j+1),n,work_row+j+1,1,c,s);
    }
    F77NAME(crot)(Mp1,Q->addr(0,j),1,work_column,1,c,s);
  }
  delete [] work_column;
  delete [] work_row;
}

template class GramSchmidtQRFactorization<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SingularValueDecomposition.C"
template<> SingularValueDecomposition<float,complex<float> >::
SingularValueDecomposition(const Matrix<float,complex<float> > &A) :
U(0),Vtranspose(0),s(0),A_original(&A) {
//TRACER_CALL(t,"SingularValueDecomposition::SingularValueDecomposition");
  int m=A.size(0),n=A.size(1);
  int mn=min(m,n);
  U=OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(m,mn);
  Vtranspose=OPERATOR_NEW OrthogonalMatrix<float,complex<float> >(mn,n);
  s=OPERATOR_NEW Vector<float,float>(mn);
  Matrix<float,complex<float> > Acopy(m,n);
  Acopy.copy(A);

  int info=0;
  complex<float> w=complex_float_undefined_;
  int lwork=-1;
  float *rwork=OPERATOR_NEW_BRACKET(float,5*mn);
  F77NAME(cgesvd)('S','S',m,n,Acopy.addr(),m,s->addr(),U->addr(),m,
    Vtranspose->addr(),mn,&w,lwork,rwork,info);
  lwork=static_cast<int>(w.real());
  complex<float> *work=OPERATOR_NEW_BRACKET(complex<float>,lwork);
  F77NAME(cgesvd)('S','S',m,n,Acopy.addr(),m,s->addr(),U->addr(),m,
    Vtranspose->addr(),mn,work,lwork,rwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;
  delete [] rwork; rwork=0;
}

template<> void
SingularValueDecomposition<float,complex<float> >::solve(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float rcond,
Factorization::TRANSPOSE_OPTION to) const {
//TRACER_CALL(t,"SingularValueDecomposition::solve");
  CHECK_TEST(rcond<float_one_);
  int m=A_original->size(0),n=A_original->size(1);
  int mn=min(m,n);
  if (rcond<float_zero_) rcond=F77NAME(slamch)('P');
  float thr=max(rcond*(*s)[0],F77NAME(slamch)('S'));
  int rank=0;
  for (int i=0;i<mn;i++) {
    if ((*s)[i]>thr) rank++;
    else break;
  }
  complex<float> *y=OPERATOR_NEW_BRACKET(complex<float>,rank);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(cgemv)('C',m,rank,float_one_,U->addr(),m,b.addr(),1,
      float_zero_,y,1);
    for (int i=0;i<rank;i++) y[i]/=(*s)[i];
    F77NAME(cgemv)('C',rank,n,float_one_,Vtranspose->addr(),mn,y,1,
      float_zero_,x.addr(),1);
  } else {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(cgemv)('N',rank,n,float_one_,Vtranspose->addr(),mn,
      b.addr(),1,float_zero_,y,1);
    for (int i=0;i<rank;i++) y[i]/=(*s)[i];
    F77NAME(cgemv)('N',m,rank,float_one_,U->addr(),m,y,1,float_zero_,
      x.addr(),1);
  }
  delete [] y; y=0;
}

template<> void
SingularValueDecomposition<float,complex<float> >::regularize(
const Vector<float,complex<float> > &b,
Vector<float,complex<float> > &x,float ridge,
Factorization::TRANSPOSE_OPTION to) const {
  CHECK_TEST(ridge>=float_zero_);
  int m=A_original->size(0),n=A_original->size(1);
  int mn=min(m,n);
  float sr=sqrt(ridge);
  complex<float> *y=OPERATOR_NEW_BRACKET(complex<float>,mn);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(cgemv)('C',m,mn,float_one_,U->addr(),m,b.addr(),1,
      float_zero_,y,1);
    for (int i=0;i<mn;i++) {
      float si=(*s)[i];
      float scale=max(si,sr);
      si/=scale;
      float ri=sr/scale;
      y[i]*=si/(scale*(si*si+ri*ri));
    }
    F77NAME(cgemv)('C',mn,n,float_one_,Vtranspose->addr(),mn,y,1,
      float_zero_,x.addr(),1);
  } else {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(cgemv)('N',mn,n,float_one_,Vtranspose->addr(),mn,b.addr(),1,
      float_zero_,y,1);
    for (int i=0;i<mn;i++) {
      float si=(*s)[i];
      float scale=max(abs(si),sr);
      si/=scale;
      float ri=sr/scale;
      y[i]*=si/(scale*(si*si+ri*ri));
    }
    F77NAME(cgemv)('N',m,mn,float_one_,U->addr(),m,y,1,float_zero_,
      x.addr(),1);
  }
  delete [] y; y=0;
}

template<> void
SingularValueDecomposition<float,complex<float> >::solve(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,float rcond,
Factorization::TRANSPOSE_OPTION to) const {
//TRACER_CALL(t,"SingularValueDecomposition::solve");
  CHECK_TEST(rcond<float_one_);
  int m=A_original->size(0),n=A_original->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1));
  int mn=min(m,n);
  if (rcond<float_zero_) rcond=F77NAME(slamch)('P');
  float thr=max(rcond*(*s)[0],F77NAME(slamch)('S'));
  int rank=0;
  for (int i=0;i<mn;i++) {
    if ((*s)[i]>thr) rank++;
    else break;
  }
  complex<float> *Y=OPERATOR_NEW_BRACKET(complex<float>,rank*k);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,X.size(0));
    F77NAME(cgemm)('C','N',rank,k,m,float_one_,U->addr(),m,B.addr(),m,
      float_zero_,Y,rank);
    for (int i=0;i<rank;i++) F77NAME(csrscl)(k,(*s)[i],Y+i,rank);
    F77NAME(cgemm)('C','N',n,k,rank,float_one_,Vtranspose->addr(),mn,
      Y,rank,float_zero_,X.addr(),n);
  } else {
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(m,X.size(0));
    F77NAME(cgemm)('N','N',rank,k,n,float_one_,Vtranspose->addr(),mn,
      B.addr(),n,float_zero_,Y,rank);
    for (int i=0;i<rank;i++) F77NAME(csrscl)(k,(*s)[i],Y+i,rank);
    F77NAME(cgemm)('N','N',m,k,rank,float_one_,U->addr(),m,Y,rank,
      float_zero_,X.addr(),m);
  }
  delete [] Y; Y=0;
}

template<> void
SingularValueDecomposition<float,complex<float> >::regularize(
const Matrix<float,complex<float> > &B,
Matrix<float,complex<float> > &X,float ridge,
Factorization::TRANSPOSE_OPTION to) const {
  CHECK_TEST(ridge>=float_zero_);
  int m=A_original->size(0),n=A_original->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1));
  int mn=min(m,n);
  float sr=sqrt(ridge);
  complex<float> *Y=OPERATOR_NEW_BRACKET(complex<float>,mn*k);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,X.size(0));
    F77NAME(cgemm)('C','N',mn,k,m,float_one_,U->addr(),m,B.addr(),m,
      float_zero_,Y,mn);
    for (int i=0;i<mn;i++) {
      float si=(*s)[i];
      float scale=max(abs(si),sr);
      si/=scale;
      float ri=sr/scale;
      F77NAME(csscal)(k,si/(scale*(si*si+ri*ri)),Y+i,mn);
    }
    F77NAME(cgemm)('C','N',n,k,mn,float_one_,Vtranspose->addr(),mn,Y,mn,
      float_zero_,X.addr(),n);
  } else {
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(m,X.size(0));
    F77NAME(cgemm)('N','N',mn,k,n,float_one_,Vtranspose->addr(),mn,
      B.addr(),n,float_zero_,Y,mn);
    for (int i=0;i<mn;i++) {
      float si=(*s)[i];
      float scale=max(abs(si),sr);
      si/=scale;
      float ri=sr/scale;
      F77NAME(csscal)(k,si/(scale*(si*si+ri*ri)),Y+i,mn);
    }
    F77NAME(cgemm)('N','N',m,k,mn,float_one_,U->addr(),m,Y,mn,
      float_zero_,X.addr(),m);
  }
  delete [] Y; Y=0;
}

template class SingularValueDecomposition<float,complex<float> >;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> CholeskyFactorization<float>::CholeskyFactorization(
const SymmetricPositiveMatrix<float> &A) : anorm(float_zero_) {
//TRACER_CALL(t,"CholeskyFactorization::CF");
  int n=A.size(0);
  int inc=1;
  for (int j=0;j<n;j++) {
    int k=n-j;
    float col_sum=F77NAME(dasum)(&k,A.addr(j,j),&inc);
    col_sum+=F77NAME(dasum)(&j,A.addr(j,0),&n);
    if (col_sum>anorm) anorm=col_sum;
  }
  SymmetricPositiveMatrix<float> *Acopy=
    OPERATOR_NEW SymmetricPositiveMatrix<float>(A);
  int info;
  char uplo='L';
  F77NAME(dpotrf)(&uplo,&n,Acopy->addr(),&n,&info);
  CHECK_SAME(int,info,0)
  L=OPERATOR_NEW LowerTrapezoidalMatrix<float>(n,n);
  for (int j=0;j<n;j++) {
    for (int i=j;i<n;i++) (*L)(i,j)=(*Acopy)(i,j);
  }
  delete Acopy;
}

template<> Matrix<float>* CholeskyFactorization<float>::solveAXeqB(
const Matrix<float> &B) const {
//TRACER_CALL(t,"CholeskyFactorization::solveAXeqB");
  int n=L->size(0); 
  CHECK_SAME(int,n,B.size(0))
  int k=B.size(1);

  Matrix<float> *X=OPERATOR_NEW Matrix<float>(B);
  char uplo='L';
  int info;
  F77NAME(dpotrs)(&uplo,&n,&k,L->addr(),&n,X->addr(),&n,&info);
  CHECK_SAME(int,info,0)
  return X;
}

template<> float CholeskyFactorization<float>::conditionNumber() const {
  int n=L->size(0);
  float rcond=float_zero_;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  char uplo='L';
  F77NAME(dpocon)(&uplo,&n,L->addr(),&n,&anorm,&rcond,work,iwork,
                  &info);
  CHECK_SAME(int,info,0)
  delete [] iwork;
  delete [] work;
  if (rcond>float_zero_) return float_one_/rcond;
  return numeric_limits<float>::infinity();
}

#ifdef __GNUC__
  template class CholeskyFactorization<float>;
  template ostream& operator<<(ostream&,
                               const CholeskyFactorization<float>&);
  template void testCholeskyFactorization(float);
#endif
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> MDMtFactorization<float>::MDMtFactorization(
const SymmetricMatrix<float> &A) : anorm(float_zero_) {
//TRACER_CALL(t,"MDMtFactorization::MDMtF");
  int n=A.size(0);
  ipiv=OPERATOR_NEW_BRACKET(int,n);
  int inc=1;
  for (int j=0;j<n;j++) {
    int k=n-j;
    float col_sum=F77NAME(dasum)(&k,A.addr(j,j),&inc);
    col_sum+=F77NAME(dasum)(&j,A.addr(j,0),&n);
    if (col_sum>anorm) anorm=col_sum;
  }
  SymmetricMatrix<float> *Acopy=OPERATOR_NEW SymmetricMatrix<float>(A);
  int info;
  char uplo='L';
  int nb=LaEnvBlockSize("DSYTRF",A);
  int lwork=n*nb;
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(dsytrf)(&uplo,&n,Acopy->addr(),&n,ipiv->addr(),work,&lwork,
                  &info);
  CHECK_SAME(int,info,0)
  delete [] work;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<float>(n,n);
  for (int j=0;j<n;j++) {
    for (int i=j;i<n;i++) (*L)(i,j)=(*Acopy)(i,j);
  }
  delete Acopy;
}

template<> Matrix<float>* MDMtFactorization<float>::solveAXeqB(
const Matrix<float> &B) const {
//TRACER_CALL(t,"MDMtFactorization::solveAXeqB");
  int n=L->size(0);
  CHECK_SAME(int,n,B.size(0))
  int k=B.size(1);

  Matrix<float> *X=OPERATOR_NEW Matrix<float>(B);
  char uplo='L';
  int info;
  F77NAME(dsytrs)(&uplo,&n,&k,L->addr(),&n,ipiv->addr(),X->addr(),&n,
                  &info);
  CHECK_SAME(int,info,0)
  return X;
}

template<> float MDMtFactorization<float>::conditionNumber() const {
  int n=L->size(0);
  float rcond=float_zero_;
  float *work=OPERATOR_NEW_BRACKET(float,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  char uplo='L';
  F77NAME(dsycon)(&uplo,&n,L->addr(),&n,ipiv->addr(),&anorm,&rcond,work,
                  iwork,&info);
  CHECK_SAME(int,info,0)
  delete [] iwork;
  delete [] work;
  if (rcond>float_zero_) return float_one_/rcond;
  return numeric_limits<float>::infinity();
}

#ifdef __GNUC__
  template class MDMtFactorization<float>;
  template ostream& operator<<(ostream&,
                               const MDMtFactorization<float>&);
  template void testMDMtFactorization(float);
#endif
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> SingularValueDecomposition<float>::SingularValueDecomposition(
const Matrix<float> &A) {
//TRACER_CALL(t,"SingularValueDecomposition::SVD");
  int m=A.size(0),n=A.size(1),minmn=min(m,n);
  Sigma=OPERATOR_NEW Vector<float>(minmn);
  int info=0;
  float *work=0;
  int lwork=5*max(m,n);
  if (m>=n) {
    U=OPERATOR_NEW Matrix<float>(A);
    V_transpose=OPERATOR_NEW Matrix<float>(n,n);
    char jobu='O';
    char jobvt='S';
//  int lwork=max(1,max( 3*n + m, 5*n - 4));
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(dgesvd)(&jobu,&jobvt,&m,&n,U->addr(),&m,Sigma->addr(),
      0,&m,V_transpose->addr(),&n,work,&lwork,&info);
  } else {
    U=OPERATOR_NEW Matrix<float>(m,m);
    V_transpose=OPERATOR_NEW Matrix<float>(A);
    char jobu='S';
    char jobvt='O';
//  int lwork=max(1,max( 3*m + n, 5*m - 4));
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(dgesvd)(&jobu,&jobvt,&m,&n,V_transpose->addr(),&m,
      Sigma->addr(),U->addr(),&m,0,&m,work,&lwork,&info);
  }
  CHECK_SAME(int,info,0)
  if (work) delete work;
}

template<> float SingularValueDecomposition<float>::conditionNumber(
float cutoff) const {
  CHECK_POSITIVE(cutoff)
  int r=rank(cutoff);
  if (r<1) return numeric_limits<float>::infinity();
  if (r==1) return float_one_;
  return (*Sigma)[0]/(*Sigma)[r-1];
}

template<> Matrix<float>* SingularValueDecomposition<float>::solveAXeqB(
const Matrix<float> &B,float cutoff,transpose_option option) const {
//TRACER_CALL(t,"SingularValueDecomposition::solveAXeqB");
  int m=U->size(0),n=V_transpose->size(1),minmn=Sigma->size();
  int K=B.size(1);
  int r=rank(cutoff);
  float alpha=float_one_;
  float beta=float_zero_;
  char trans='N';
  switch (option) {
    case no_matrix_transpose: {
      CHECK_SAME(int,m,B.size(0))
      Matrix<float> C(r,K);
      char transa='T';
      F77NAME(dgemm)(&transa,&trans,&r,&K,&m,&alpha,U->addr(),&m,
        B.addr(),&m,&beta,C.addr(),&r);
      for (int i=0;i<r;i++) {
        float temp=1./(*Sigma)[i];
        F77NAME(dscal)(&K,&temp,C.addr(i,0),&r);
      }
      Matrix<float> *X=OPERATOR_NEW Matrix<float>(n,K);
      F77NAME(dgemm)(&transa,&trans,&n,&K,&r,&alpha,
        V_transpose->addr(),&minmn,C.addr(),&r,&beta,X->addr(),&n);
      return X;
      break;
    }
    case matrix_transpose: {
      CHECK_SAME(int,n,B.size(0))
      Matrix<float> C(r,K);
      F77NAME(dgemm)(&trans,&trans,&r,&K,&n,&alpha,V_transpose->addr(),
        &n,B.addr(),&n,&beta,C.addr(),&r);
      for (int i=0;i<r;i++) {
        float temp=1./(*Sigma)[i];
        F77NAME(dscal)(&K,&temp,C.addr(i,0),&r);
      }
      Matrix<float> *X=OPERATOR_NEW Matrix<float>(m,K);
      F77NAME(dgemm)(&trans,&trans,&m,&K,&r,&alpha,U->addr(),&m,
        C.addr(),&r,&beta,X->addr(),&m);
      return X;
      break;
    }
    default:
      return 0;
  }
}

#ifdef __GNUC__
  template class SingularValueDecomposition<float>;
  template ostream& operator<<(ostream&,const SingularValueDecomposition<float>&);
  template void testSingularValueDecomposition(float);
#endif
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//#include "LinearProgram.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
#ifdef __GNUC__
  template class VirtualLinearProgram<float>;
  template ostream& operator<<(ostream&,
                               const VirtualLinearProgram<float>&);
#endif
template<> float VirtualLinearProgram<float>::zero_ = float_zero_;
template<> float VirtualLinearProgram<float>::one_ = float_one_;
template<> float VirtualLinearProgram<float>::huge_ = numeric_limits<float>::infinity();

template<> float VirtualLinearProgram<float>::currentValue() const {
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
  template class SFLinearProgram<float>;
  template ostream& operator<<(ostream&,const SFLinearProgram<float>&);
#endif

template<> SFLinearProgram<float>::SFLinearProgram(
const Matrix<float> &Ain,const Matrix<float> &bin,
const Matrix<float> &cin) : VirtualLinearProgram<float>(Ain,bin,cin) {
  int m=A.size(0),n=A.size(1);
  assert(m<=n);

  int incb=1;
  int imin=F77NAME(idmin)(&m,b.addr(),&incb)-1;
//if (b(imin,0)<zero_) throw bad_constraint_vector;
  assert(b(imin,0)>=zero_);

  r_trans = OPERATOR_NEW Matrix<float>(1,n-m);
  h=OPERATOR_NEW Matrix<float>(m,1);
  column=OPERATOR_NEW_BRACKET(int,n);
  for (int j=0;j<n;j++) (*column)[j]=j;
  new_basic=n;
}

template<> status_option 
SFLinearProgram<float> ::findBasicFeasibleGuess() {
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
  SquareMatrix<float> *active=OPERATOR_NEW SquareMatrix<float>(m,m);
  A.copyInto(*active);
  basic_inverse=inverse(*active);
  delete active;

  Matrix<float> *x_basic = productOf(*basic_inverse,b);
  int incx=1;
  int imin=F77NAME(idmin)(&m,x_basic->addr(),&incx)-1;
  float xmin=(*x_basic)(imin,0);
  x->fillWith(float_zero_);
  x_basic->copyInto(*x);
  delete x_basic;

  if (xmin<float_zero_) {
//  augment the program and find a basic feasible solution
    Matrix<float> *augmented_A=
      OPERATOR_NEW Matrix<float>(m,m+n,float(float_zero_));
    Matrix<float> *augmented_c_trans=
      OPERATOR_NEW Matrix<float>(1,m+n,float(float_zero_));
    int incaA=1;
    for (int j=0;j<m;j++) (*augmented_A)(j,j)=float_one_;
    { for (int j=0;j<n;j++) {
      F77NAME(dcopy)(&m,A.addr(0,j),&incA,augmented_A->addr(0,j+m),
		     &incaA);
      (*augmented_c_trans)(0,j+m)=c_trans(0,j);
    } }

//  initialize and solve the augmented program
    SFLinearProgram<float> LP(*augmented_A,b,*augmented_c_trans);
    LP.basic_inverse=OPERATOR_NEW SquareMatrix<float>(m,m);
    augmented_A->copyInto(*LP.basic_inverse);
    LP.x->fillWith(float_zero_);
    b.copyInto(*(LP.x));
    LP.y_trans->fillWith(float_zero_);
    c_trans.copyInto(*(LP.r_trans));
    int incr=LP.r_trans->size(0);
    LP.new_basic=F77NAME(idmin)(&n,r_trans->addr(),&incr)-1;
    LP.current_status=primal_feasible;
    while (LP.simplexStep() == primal_feasible) {;}

//  if no augmented variables are basic, store basic variables here
//  we could have trouble if an augmented variable were basic and zero
    x->fillWith(float_zero_);
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
    float yj=float_zero_;
    for (int k=0;k<m;k++) {
      yj += c_trans(0,(*column)[k]) * (*basic_inverse)(k,j);
    }
    (*y_trans)(0,j)=yj;
  }

//find cost reduction vector
  int jmin=largestCostReduction();
  current_status=
    ((*r_trans)(0,jmin)<float_zero_ ? primal_feasible : optimal);
  new_basic=jmin+m;

  return current_status;
}

template<> int SFLinearProgram<float>::largestCostReduction() {
//TRACER_CALL(t,"SFLinearProgram::largestCostReduction");
  int m=A.size(0),n=A.size(1);
  int jmin=0; float rmin=float_zero_;
  int incy=y_trans->size(0);
  int incA=1;
  for (int j=0;j<n-m;j++) {
    int colj=(*column)[j+m];
    float rj=c_trans(0,colj)
          -F77NAME(ddot)(&m,y_trans->addr(),&incy,A.addr(0,colj),&incA);
    if (rj < rmin) { jmin=j; rmin=rj; }
    (*r_trans)(0,j)=rj;
  }
  int NminM=n-m;
  int incr=r_trans->size(0);
  return F77NAME(idmin)(&NminM,r_trans->addr(),&incr)-1;
}

template<> status_option SFLinearProgram<float>::simplexStep() {
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
  float epsmin=float_zero_;
  char trans='n';
  float alpha=float_one_;
  float beta=float_zero_;
  int inch=1;
  int incA=1;
  int col_new_basic=(*column)[new_basic];
  F77NAME(dgemv)(&trans,&m,&m,&alpha,basic_inverse->addr(),&m,
		 A.addr(0,col_new_basic),&incA,&beta,h->addr(),&inch);
  for (int i=0;i<m;i++) {
    float hi=(*h)(i,0);
    int coli=(*column)[i];
    if (hi>float_zero_) {
      float xi=(*x)(coli,0);
      if (imin>=n || xi < epsmin*hi) { 
	imin=i; colmin=coli; epsmin=xi/hi;
      } else if (xi<=float_zero_ && coli<colmin) { 
	imin=i; colmin=coli; epsmin=float_zero_; 
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
  (*x)((*column)[imin],0)=float_zero_;
  (*x)(col_new_basic,0)=epsmin;
#ifdef DEBUG
//cout << "\tupdated x" << endl;
//cout << *x << endl;
#endif

//update dual solution
  float hmin=(*h)(imin,0);
  alpha=(*r_trans)(0,new_basic-m)/hmin;
  int incy=1;
  F77NAME(saxpy)(&m,&alpha,basic_inverse->addr(imin,0),&m,
		 y_trans->addr(),&incy);
#ifdef DEBUG
//cout << "\tupdated y_trans" << endl;
//cout << *y_trans << endl;
#endif

//update inverse of active constraint matrix
  { for (int i=0;i<m;i++) {
    if (i!=imin) {
      alpha=-(*h)(i,0)/hmin;
      F77NAME(saxpy)(&m,&alpha,basic_inverse->addr(imin,0),&m,
		     basic_inverse->addr(i,0),&m);
    }
  } }
  alpha=float_one_/hmin;
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
    ((*r_trans)(0,jmin)<float_zero_ ? primal_feasible : optimal);
  new_basic=jmin+m;
#ifdef DEBUG
//cout << "\tupdated r_trans = " << endl;
//cout << *r_trans << endl;
//cout << "new_basic = " << (*column)[new_basic] << endl;
#endif

  return current_status;
}

template<> void SFLinearProgram<float>::costBounds(
Matrix<float> &lower_trans,Matrix<float> &upper_trans) const {
  int m=A.size(0),n=A.size(1);
  CHECK_SAME(int,n,lower_trans.size(1))
  CHECK_SAME(int,n,upper_trans.size(1))

  int inc0=1;
  int inc1=A.size(0);
//perturbation in basic variable
  for (int i=0;i<m;i++) {
    float lo=-huge_;
    float hi=huge_;
    for (int j=m;j<n;j++) {
      int colj=(*column)[j];
      float denom=F77NAME(ddot)(&m,basic_inverse->addr(i,0),&inc1,
				A.addr(0,colj),&inc0);
      float rj=(*r_trans)(0,j-m);
      if (denom<float_zero_ && denom*lo>rj) lo=rj/denom;
      else if (denom>float_zero_ && denom*hi>rj) hi=rj/denom;
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

template<> float SFLinearProgram<float>::costSensitivity(int j,
Matrix<float> &dual_derivative) const {
  int m=A.size(0);
  CHECK_SAME(int,m,dual_derivative.size(1))
  for (int i=0;i<m;i++) {
    if ((*column)[i]==j) {
      for (int k=0;k<m;k++) dual_derivative(0,k)=(*basic_inverse)(i,k);
      return (*x)(j,0);
    }
  }
  dual_derivative.fillWith(float_zero_);
  return float_zero_;
}

template<> void SFLinearProgram<float>::arrayBounds(int i,int j,
float &lower,float &upper) const {
//TRACER_CALL(t,"SFLinearProgram::arrayBounds");
  int m=A.size(0),n=A.size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,n)
  char safe='S';
  float sfmin=F77NAME(slamch)(&safe);
  upper=float_one_/sfmin;
  lower=-upper;
  int incA=1;

//perturbation in basic column of A
  float yi=(*y_trans)(0,i);
  float xj=(*x)(j,0);
  for (int jj=0;jj<m;jj++) {
    if ((*column)[jj]==j) {
      float ainv_ji=(*basic_inverse)(jj,i);
      if (fabs(ainv_ji)>sfmin) {
	float bound=-float_one_/ainv_ji;
	if (ainv_ji>float_zero_) lower=bound;
	else upper=bound;
      }
      for (int k=0;k<m;k++) {
        if (k!=jj) {
	  float xk=(*x)((*column)[k],0);
	  float denom=(*basic_inverse)(k,i)*xj-ainv_ji*xk;
	  if (fabs(denom)>xk*sfmin) {
	    float bound=xk/denom;
	    if (denom>float_zero_ && bound<upper) upper=bound;
	    if (denom<float_zero_ && bound>lower) lower=bound;
	  }
	}
      }
      for (int kk=m;kk<n;kk++) {
	float rk=(*r_trans)(0,kk-m);
	float term=F77NAME(ddot)(&m,basic_inverse->addr(jj,0),&m,
			   A.addr(0,(*column)[kk]),&incA);
	float denom=yi*term-ainv_ji*rk;
	if (fabs(denom)>rk*sfmin) {
	  float bound=rk/denom;
	  if (denom>float_zero_ && bound<upper) upper=bound;
	  if (denom<float_zero_ && bound>lower) lower=bound;
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
      float rj=(*r_trans)(0,jj-m);
      if (fabs(yi)>rj*sfmin) {
	float bound=A(i,j)+rj/yi;
	if (yi>float_zero_ && bound<upper) upper=bound;
	if (yi<float_zero_ && bound>lower) bound=lower;
      }
    }
  } }
}

template<> void testSFLinearProgram(float scalar) {
  scalar=float_one_;
//Giapetto problem, p. 45 of Winston 
  Matrix<float> A(3,5); 
  Matrix<float> b(3,1); 
  Matrix<float> c(1,5); 
  c(0,0)=0.; c(0,1)=0.; c(0,2)=0.; c(0,3)=-3.; c(0,4)=-2.;
  A(0,0)=1.; A(0,1)=0.; A(0,2)=0.; A(0,3)= 2.; A(0,4)= 1.; b(0,0)=100.;
  A(1,0)=0.; A(1,1)=1.; A(1,2)=0.; A(1,3)= 1.; A(1,4)= 1.; b(1,0)= 80.;
  A(2,0)=0.; A(2,1)=0.; A(2,2)=1.; A(2,3)= 1.; A(2,4)= 0.; b(2,0)= 40.;

//
//Dorian problem, p. 130 of Winston 
//Matrix<float> A(4,7); 
//Matrix<float> b(4,1); 
//Matrix<float> c(1,7); 
//c(0,0)=0.; c(0,1)=0.; c(0,2)=0.; c(0,3)=0.; c(0,4)=-60.; c(0,5)=-30. ; c(0,6)=-20. ;
//A(0,0)=1.; A(0,1)=0.; A(0,2)=0.; A(0,3)=0.; A(0,4)=  8.; A(0,5)=  6. ; A(0,6)= 1. ; b(0,0)=48.;
//A(1,0)=0.; A(1,1)=1.; A(1,2)=0.; A(1,3)=0.; A(1,4)=  4.; A(1,5)=  2. ; A(1,6)= 1.5; b(1,0)=20.;
//A(2,0)=0.; A(2,1)=0.; A(2,2)=1.; A(2,3)=0.; A(2,4)=  2.; A(2,5)=  1.5; A(2,6)= 0.5; b(2,0)= 8.;
//A(3,0)=0.; A(3,1)=0.; A(3,2)=0.; A(3,3)=1.; A(3,4)=  0.; A(3,5)=  1. ; A(3,6)= 0. ; b(3,0)= 5.;

//first two unknowns not basic, feasible
//Matrix<float> A(2,5); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,5); 
//c(0,0)= 5.; c(0,1)=0.; c(0,2)=0.; c(0,3)= 3.; c(0,4)=-2.;
//A(0,0)=-6.; A(0,1)=0.; A(0,2)=1.; A(0,3)=-2.; A(0,4)= 2.; b(0,0)= 6.;
//A(1,0)=-3.; A(1,1)=1.; A(1,2)=0.; A(1,3)= 6.; A(1,4)= 3.; b(1,0)=15.;

//problem unbounded, p. 146
//Matrix<float> A(2,6); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,6); 
//c(0,0)=0.; c(0,1)=0.; c(0,2)=-36.; c(0,3)=-30.; c(0,4)= 3.; c(0,5)= 4.;
//A(0,0)=1.; A(0,1)=0.; A(0,2)=  1.; A(0,3)=  1.; A(0,4)=-1.; A(0,5)= 0.; b(0,0)= 5.;
//A(1,0)=0.; A(1,1)=1.; A(1,2)=  6.; A(1,3)=  5.; A(1,4)= 0.; A(1,5)=-1.; b(1,0)=10.;

//basic feasible solution is degenerate, p. 159
//Matrix<float> A(2,4);
//Matrix<float> b(2,1); 
//Matrix<float> c(1,4); 
//c(0,0)=0.; c(0,1)=0.; c(0,2)=-5.; c(0,3)=-2.;
//A(0,0)=1.; A(0,1)=0.; A(0,2)= 1.; A(0,3)= 1.; b(0,0)=6.;
//A(1,0)=0.; A(1,1)=1.; A(1,2)= 1.; A(1,3)=-1.; b(1,0)=0.;

//problem cycles without lexico-graphic ordering, p. 161
//Matrix<float> A(2,6);
//Matrix<float> b(2,1); 
//Matrix<float> c(1,6); 
//c(0,0)=0.; c(0,1)=0.; c(0,2)=-2.  ; c(0,3)=-3.; c(0,4)= 1.   ; c(0,5)=12.;
//A(0,0)=1.; A(0,1)=0.; A(0,2)=-2.  ; A(0,3)=-9.; A(0,4)= 1.   ; A(0,5)= 9.; b(0,0)=0.;
//A(1,0)=0.; A(1,1)=1.; A(1,2)=1./3.; A(1,3)= 1.; A(1,4)=-1./3.; A(1,5)=-2.; b(1,0)=0.;
//

  int m=A.size(0),n=A.size(1);
  SFLinearProgram<float> LP(A,b,c);
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

  Matrix<float> *c_lower=OPERATOR_NEW Matrix<float>(1,n);
  Matrix<float> *c_upper=OPERATOR_NEW Matrix<float>(1,n);
  LP.costBounds(*c_lower,*c_upper);
  cout << "lower bounds on cost vector:" << endl;
  cout << *c_lower << endl;
  cout << "upper bounds on cost vector:" << endl;
  cout << *c_upper << endl;

  Matrix<float> *dual_derivative=OPERATOR_NEW Matrix<float>(1,m);
  for (int j=0;j<n;j++) {
    float dzdc=LP.costSensitivity(j,*dual_derivative);
    cout << "\nj,dzdc = " << j << " " << dzdc << endl;
    cout << "dual_derivative = " << *dual_derivative << endl;
  }

  Matrix<float> *b_lower=OPERATOR_NEW Matrix<float>(m,1);
  Matrix<float> *b_upper=OPERATOR_NEW Matrix<float>(m,1);
  LP.constraintBounds(*b_lower,*b_upper);
  cout << "lower bounds on constraint vector:" << endl;
  cout << *b_lower << endl;
  cout << "upper bounds on constraint vector:" << endl;
  cout << *b_upper << endl;

  Matrix<float> *primal_derivative=OPERATOR_NEW Matrix<float>(n,1);
  for (int i=0;i<m;i++) {
    float dzdb=LP.constraintSensitivity(i,*primal_derivative);
    cout << "\ni,dzdb = " << i << " " << dzdb << endl;
    cout << "primal_derivative = " << *primal_derivative << endl;
  }

  float A_lower,A_upper;
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
      float dzda=
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
  template class LinearProgram<float>;
//template LinearProgram<float>::LinearProgram(const Matrix<float>&,const Matrix<float>&,const Matrix<float>&);
//template float LinearProgram<float>::currentValue() const;
//template status_option LinearProgram<float>::simplexStep();
//template void LinearProgram<float>::costBounds(Matrix<float>&,
//   Matrix<float>&) const;
//template void LinearProgram<float>::constraintBounds(Matrix<float>&,
//   Matrix<float>&) const;
//template void LinearProgram<float>::printOn(ostream&) const;
  template ostream& operator<<(ostream&,const LinearProgram<float>&);
#endif

template<> status_option LinearProgram<float>::findBasicFeasibleGuess() {
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
LinearProgram<float>::findPrimalBasicFeasibleGuess() {
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
  Matrix<float> *temp_A=OPERATOR_NEW Matrix<float>(m,n+m-ib,float_zero_);
  A.copyInto(*temp_A);
  Matrix<float> *temp_c_trans=OPERATOR_NEW Matrix<float>(1,n+m-ib,float_zero_);
  int pos_b=n;
  for (int i=0;i<m;i++) {
    if (b(i,0)>float_zero_) {
      (*temp_A)(i,pos_b)=float_one_;
      (*temp_c_trans)(0,pos_b)=float_one_;
      pos_b++;
    }
  }
  LinearProgram temp_LP(*temp_A,b,*temp_c_trans);
  temp_LP.basic_number=pos_b-n;
  delete temp_LP.basic_inverse;
  temp_LP.basic_inverse=OPERATOR_NEW SquareMatrix<float>(temp_LP.basic_number,0);
  basic_inverse->fillWith(float_zero_);
  for (int k=0;k<temp_LP.basic_number;k++) {
    (*temp_LP.basic_inverse)(k,k)=float_one_;
  }
  temp_LP.x->fillWith(float_zero_);
  pos_b=n;
  { for (int i=0;i<m;i++) {
    if (b(i,0)>float_zero_) (*temp_LP.x)(pos_b++,0)=b(i,0);
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
  x->fillWith(float_zero_);
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
LinearProgram<float>::findDualBasicFeasibleGuess() {
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
  Matrix<float> *temp_A=OPERATOR_NEW Matrix<float>(m+n-ic,n,float_zero_);
  A.copyInto(*temp_A);
  Matrix<float> *temp_b=OPERATOR_NEW Matrix<float>(m+n-ic,1,float_zero_);
  int neg_c=m;
  for (int j=0;j<n;j++) {
    if (c_trans(0,j)<float_zero_) {
      (*temp_A)(neg_c,j)=-float_one_;
      (*temp_b)(neg_c,0)=float_one_;
      neg_c++;
    }
  }
  LinearProgram temp_LP(*temp_A,*temp_b,c_trans);
  temp_LP.basic_number=neg_c;
  delete temp_LP.basic_inverse;
  temp_LP.basic_inverse=OPERATOR_NEW SquareMatrix<float>(temp_LP.basic_number,0);
  basic_inverse->fillWith(float_zero_);
  for (int k=0;k<temp_LP.basic_number;k++) {
    (*temp_LP.basic_inverse)(k,k)=float_one_;
  }
  temp_LP.y_trans->fillWith(float_zero_);
  neg_c=m;
  { for (int j=0;j<n;j++) {
    if (c_trans(0,j)<float_zero_) (*temp_LP.y_trans)(0,neg_c++)=-c_trans(0,j);
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
  y_trans->fillWith(float_zero_);
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

template<> void LinearProgram<float>::computeSolution() {
//TRACER_CALL(t,"LinearProgram::computeSolution");
#ifdef DEBUG
//cout << "basic_number = " << basic_number << endl;
//cout << "row = \n" << *row << endl;
//cout << "column = \n" << *column << endl;
//cout << "basic_matrix :\n" << *basic_matrix << endl;
#endif
  Matrix<float> *rhs_basic=OPERATOR_NEW Matrix<float>(basic_number,1);
  Matrix<float> *soln_basic=OPERATOR_NEW Matrix<float>(basic_number,1);
  for (int i=0;i<basic_number;i++) {
    (*rhs_basic)(i,0)=b((*row)[i],0);
  }
#ifdef DEBUG
//cout << "rhs_basic = \n" << *rhs_basic << endl;
#endif
  char trans='n';
  float alpha=float_one_;
  float beta=float_zero_;
  int inc_rhs=1;
  int inc_soln=1;
  F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
    basic_inverse->addr(),&basic_number,rhs_basic->addr(),&inc_rhs,
    &beta,soln_basic->addr(),&inc_soln);
#ifdef DEBUG
//cout << "soln_basic = \n" << *soln_basic << endl;
#endif
  x->fillWith(float_zero_);
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

  y_trans->fillWith(float_zero_);
  { for (int i=0;i<basic_number;i++) {
    (*y_trans)(0,(*row)[i])=(*soln_basic)(i,0);
  } }
  delete soln_basic;
  delete rhs_basic;
#ifdef DEBUG
//cout << "y_trans = \n" << *y_trans << endl;
#endif
}

template<> int LinearProgram<float>::basicPrimalPivot(Matrix<float> &h,
float &ratio) {
//TRACER_CALL(t,"LinearProgram::basicPrimalPivot");
  int i_basic=basic_number;
  ratio=huge_;
  if (basic_number>0) {
    Matrix<float> rhs;
    rhs.copy(h);
    char trans='n';
    float alpha=float_one_;
    float beta=float_zero_;
    int inc_rhs=1;
    int inc_h=1;
    F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
      basic_inverse->addr(),&basic_number,rhs.addr(),&inc_rhs,
      &beta,h.addr(),&inc_h);
#ifdef DEBUG
    cout << "after solve, h = \n" << h << endl;
#endif
    for (int i=0;i<basic_number;i++) {
      float hi=h(i,0);
      if (hi<=float_zero_) continue;
      float xi=(*x)((*column)[i],0);
      if (i_basic>=basic_number || hi*ratio<xi) {
	i_basic=i;
	ratio=xi/hi;
      }
    }
  }
  return i_basic;
}

template<> void LinearProgram<float>::add(int i_non_basic,
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
  SquareMatrix<float> *new_basic_inverse=
    OPERATOR_NEW SquareMatrix<float>(basic_number+1,basic_number+1);
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

  Matrix<float> *temp=OPERATOR_NEW Matrix<float>(basic_number,1);
  char trans='n';
  float alpha=float_one_;
  float beta=float_zero_;
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
  float diag=A((*row)[basic_number],(*column)[basic_number]);
  for (int k=0;k<basic_number;k++) {
    diag += A((*row)[basic_number],(*column)[k]) 
          * (*new_basic_inverse)(k,basic_number);
  }
  diag=float_one_/diag;
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

template<> void  LinearProgram<float>::switchRows(int i_basic,
int i_non_basic) {
//TRACER_CALL(t,"LinearProgram::switchRows");
#ifdef DEBUG
  cout << "switching basic row " << (*row)[i_basic] 
       << " with non-basic row " << (*row)[i_non_basic] << endl;
#endif
//update inverse of active constraint matrix
  Matrix<float> *temp_rhs=OPERATOR_NEW Matrix<float>(1,basic_number); 
  Matrix<float> *temp_soln=OPERATOR_NEW Matrix<float>(1,basic_number); 
  for (int j=0;j<basic_number;j++) {
    (*temp_rhs)(0,j)=A((*row)[i_non_basic],(*column)[j]);
  }
  char trans='T';
  float alpha=float_one_;
  float beta=float_zero_;
  int inc_rhs=temp_rhs->size(0);
  int inc_soln=temp_soln->size(0);
  F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
    basic_inverse->addr(),&basic_number,temp_rhs->addr(),&inc_rhs,
    &beta,temp_soln->addr(),&inc_soln);
  float ti=float_one_/(*temp_soln)(0,i_basic);
  int incA=1;
  { for (int j=0;j<basic_number;j++) {
    if (j!=i_basic) {
      alpha=-(*temp_soln)(0,j)*ti;
      F77NAME(saxpy)(&basic_number,&alpha,
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

template<> void  LinearProgram<float>::switchColumns(int j_basic,
int j_non_basic) {
//TRACER_CALL(t,"LinearProgram::switchColumns");
#ifdef DEBUG
  cout << "switching basic column " << (*column)[j_basic] 
       << " with non-basic column " << (*column)[j_non_basic] << endl;
#endif
//update inverse of active constraint matrix
  Matrix<float> *temp_rhs=OPERATOR_NEW Matrix<float>(basic_number,1); 
  Matrix<float> *temp_soln=OPERATOR_NEW Matrix<float>(basic_number,1); 
  for (int i=0;i<basic_number;i++) {
    (*temp_rhs)(i,0)=A((*row)[i],(*column)[j_non_basic]);
  }
  char trans='n';
  float alpha=float_one_;
  float beta=float_zero_;
  int inc_rhs=1;
  int inc_soln=1;
  F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
    basic_inverse->addr(),&basic_number,temp_rhs->addr(),&inc_rhs,
    &beta,temp_soln->addr(),&inc_soln);
  float tj=float_one_/(*temp_soln)(j_basic,0);
  { for (int i=0;i<basic_number;i++) {
    if (i!=j_basic) {
      alpha=-(*temp_soln)(i,0)*tj;
      F77NAME(saxpy)(&basic_number,&alpha,
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
template<> int LinearProgram<float>::basicDualPivot(Matrix<float> &h,
float &ratio) {
//TRACER_CALL(t,"LinearProgram::basicDualPivot");
  int j_basic=basic_number;
  ratio=huge_;
  if (basic_number>0) {
    Matrix<float> rhs;
    rhs.copy(h);
    char trans='T';
    float alpha=float_one_;
    float beta=float_zero_;
    int inc_rhs=1;
    int inc_h=1;
    F77NAME(dgemv)(&trans,&basic_number,&basic_number,&alpha,
      basic_inverse->addr(),&basic_number,rhs.addr(),&inc_rhs,
      &beta,h.addr(),&inc_h);
#ifdef DEBUG
    cout << "after solve, h = \n" << h << endl;
#endif
    for (int j=0;j<basic_number;j++) {
      float hj=h(j,0);
      if (hj<=float_zero_) continue;
      float yj=(*y_trans)(0,(*row)[j]);
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

void testLinearProgram(float scalar) {
  scalar=float_one_;
//
//problem 6, Winston p. 286
//solves dual
//adds one variable
//y_trans=[1,0]
//Matrix<float> A(2,2); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,2); 
//c(0,0)=1.; c(0,1)= 2.;
//A(0,0)=1.; A(0,1)=-1.; b(0,0)= 3.;
//A(1,0)=1.; A(1,1)= 1.; b(1,0)= 1.;

//problem 6, Winston p. 286
//solves primal
//adds one variable
//x=[1,0]
//Matrix<float> A(2,2); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,2); 
//c(0,0)=-3.; c(0,1)=-1.;
//A(0,0)=-1.; A(0,1)=-1.; b(0,0)=-1.;
//A(1,0)= 1.; A(1,1)=-1.; b(1,0)=-2.;

//Leather Limited problem, Winston p. 121
//solves dual
//adds row 0 and col 1; adds row 1 and col 0
//y_trans = [20,20]
//Matrix<float> A(2,2); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,2); 
//c(0,0)=40.; c(0,1)=60.;
//A(0,0)= 1.; A(0,1)= 2.; b(0,0)= 4.;
//A(1,0)= 1.; A(1,1)= 1.; b(1,0)= 3.;

//Leather Limited problem, Winston p. 121
//solves primal
//adds row 1 and col 0; adds row 0 and col 1
//x = [20,20]
//Matrix<float> A(2,2); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,2); 
//c(0,0)=-4.; c(0,1)=-3.;
//A(0,0)=-1.; A(0,1)=-1.; b(0,0)=-40.;
//A(1,0)=-2.; A(1,1)=-1.; b(1,0)=-60.;

//Giapetto problem, Winston p. 45
//solves dual
//adds row 0 & col 2; adds row 1 and col 0; switch cols 1 and 2
//y_trans = [20,60]
//Matrix<float> A(2,3); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,3); 
//c(0,0)=100.; c(0,1)=80.; c(0,2)=40.;
//A(0,0)=  2.; A(0,1)= 1.; A(0,2)= 1.; b(0,0)= 3.;
//A(1,0)=  1.; A(1,1)= 1.; A(1,2)= 0.; b(1,0)= 2.;

//Giapetto problem, Winston p. 45
//solves primal
//adds row 2 & col 0; adds row 0 and col 1; switch rows 1 and 2
//x = [20,60]
//Matrix<float> A(3,2); 
//Matrix<float> b(3,1); 
//Matrix<float> c(1,2); 
//c(0,0)=-3.; c(0,1)=-2.;
//A(0,0)=-2.; A(0,1)=-1.; b(0,0)=-100.;
//A(1,0)=-1.; A(1,1)=-1.; b(1,0)= -80.;
//A(2,0)=-1.; A(2,1)= 0.; b(2,0)= -40.;

//Dorian problem, Winston p. 58
//solves dual
//adds row 0 & col 0; adds row 1 and col 1
//y_trans = [5,7.5]
//Matrix<float> A(2,2); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,2); 
//c(0,0)=50.; c(0,1)=100.;
//A(0,0)= 7.; A(0,1)=  2.; b(0,0)= 28.;
//A(1,0)= 2.; A(1,1)= 12.; b(1,0)= 24.;

//Auto company, Winston p. 61
//multiple solutions
//solves dual
//adds row 0 and col 0
//y_trans = [40,0]
//Matrix<float> A(2,2); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,2); 
//c(0,0)=120.; c(0,1)=50.;
//A(0,0)=  3.; A(0,1)= 1.; b(0,0)= 3.;
//A(1,0)=  2.; A(1,1)= 1.; b(1,0)= 2.;
//

//
//Problem 2, Winston p. 198
//solves primal
//adds row 0 and col 0
//x = [40,0]
//Matrix<float> A(2,2); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,2); 
//c(0,0)=-4.; c(0,1)= 1.;
//A(0,0)=-3.; A(0,1)=-1.; b(0,0)=-6.;
//A(1,0)= 1.; A(1,1)=-2.; b(1,0)= 0.;

//Auto company, Winston p. 63
//infeasible
//solves dual
//adds row 0 and col 0
//Matrix<float> A(2,4); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,4); 
//c(0,0)=120.; c(0,1)=50.; c(0,2)=-30.; c(0,3)=-20.;
//A(0,0)=  3.; A(0,1)= 1.; A(0,2)=  1.; A(0,3)=  0.; b(0,0)= 3.;
//A(1,0)=  2.; A(1,1)= 1.; A(1,2)=  0.; A(1,3)=  1.; b(1,0)= 2.;

//Winston p. 64
//unbounded
//Matrix<float> A(2,2); 
//Matrix<float> b(2,1); 
//Matrix<float> c(1,2); 
//c(0,0)= 1.; c(0,1)=-6.;
//A(0,0)= 1.; A(0,1)= 2.; b(0,0)= 2.;
//A(1,0)=-1.; A(1,1)= 1.; b(1,0)=-1.;

//Dakotta problem, Winston p. 280
//solves dual
//adds row 0 and col 2; adds row2 and col 1
//y_trans = [2,0,8]
//Matrix<float> A(3,3); 
//Matrix<float> b(3,1); 
//Matrix<float> c(1,3); 
//c(0,0)= 48.; c(0,1)=20. ; c(0,2)= 8. ;
//A(0,0)=  8.; A(0,1)= 4. ; A(0,2)= 2. ; b(0,0)= 60.;
//A(1,0)=  6.; A(1,1)= 2. ; A(1,2)= 1.5; b(1,0)= 30.;
//A(2,0)=  1.; A(2,1)= 1.5; A(2,2)= 0.5; b(2,0)= 20.;

//Diet problem, Winston p. 68
//solves dual
//adds row 0 and col 1; adds row 2 and col 2; switches row 0 with row 1
//y_trans=[0,2.5,7.5,0]
//Matrix<float> A(4,4); 
//Matrix<float> b(4,1); 
//Matrix<float> c(1,4); 
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
  Matrix<float> A(7,7); 
  Matrix<float> b(7,1); 
  Matrix<float> c(1,7); 
  for (int j=0;j<7;j++) {
    c(0,j)=1.;
    for (int k=0;k<5;k++) A((j+k)%7,j)=1.;
    A((j+5)%7,j)=A((j+6)%7,j)=0.;
  }
  b(0,0)=17.; b(1,0)=13.; b(2,0)=15.; b(3,0)=19.; b(4,0)=14.;
  b(5,0)=16.; b(6,0)=11.;

//int m=A.size(0),n=A.size(1);
  LinearProgram<float> LP(A,b,c);
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
