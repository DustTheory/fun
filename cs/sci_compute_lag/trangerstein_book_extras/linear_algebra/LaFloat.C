#include <complex>
#include <float.h>
#include <limits>
#include <math.h>
#include <string.h>
#include "Tracer.H"
#include "Vector.H"
const float float_half_=0.5;
const float float_mone_=-1.;
const float float_one_=1.;
const float float_undefined_=numeric_limits<float>::infinity();
const float float_zero_=0.;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  int F77NAME(isamax)(const int &n,const float *x,const int &incx);
  int F77NAME(ismax)(const int &n,const float *a,const int &inca);
  int F77NAME(ismin)(const int &n,const float *a,const int &inca);
  int F77NAME(issumn)(const int &n,const float *a,const int &inca);
  int F77NAME(issump)(const int &n,const float *a,const int &inca);
  float F77NAME(sasum)(const int&,const float*,const int&);
  void F77NAME(saxpy)(const int &n,const float &a,const float *x,
    const int &incx,float *y,const int &incy);
  void F77NAME(scopy)(const int &n,const float *sx,const int &incx,
    float *sy,const int &incy);
  float F77NAME(sdot)(const int &n,const float *sx,const int &incx,
    const float *sy,const int &incy);
  void F77NAME(sgbcon)(const char &norm,const int &n,const int &kl,
    const int &ku,const float *AB,const int &ldab,const int *ipiv,
    const float &anorm,float &rcond,float *work,int *iwork,int &info);
  void F77NAME(sgbconnp)(const char &norm,const int &n,const int &kl,
    const int &ku,const float *AB,const int &ldab,const float &anorm,
    float &rcond,float *work,int *iwork,int &info);
  void F77NAME(sgbequ)(const int &m,const int &n,const int &kl,
    const int &ku,const float *AB,const int &ldab,float *r,float *c,
    float &rowcnd,float &colcnd,float &amax,int &info);
  void F77NAME(sgbmv)(const char &trans,const int &m,const int &n,
    const int &kl,const int &ku,const float &alpha,const float *A,
    const int &lda,const float *x,const int &incx,const float &beta,
    float *y,const int &incy);
  void F77NAME(sgbtf2)(const int &m,const int &n,const int &kl,
    const int &ku,float *AB,const int &ldab,int *ipiv,int &info);
  void F77NAME(sgbtf2np)(const int &m,const int &n,const int &kl,
    const int &ku,float *AB,const int &ldab,int &info);
  void F77NAME(sgbtrf)(const int &m,const int &n,const int &kl,
    const int &ku,float *AB,const int &ldab,int *ipiv,int &info);
  void F77NAME(sgbtrs)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const float *AB,const int &ldab,
    const int *ipiv,float *B,const int &ldb,int &info);
  void F77NAME(sgbtrsnp)(const char &trans,const int &n,const int &kl,
    const int &ku,const int &nrhs,const float *AB,const int &ldab,
    float *B,const int &ldb,int &info);
  void F77NAME(sgecon)(const char &trans,const int &n,const float *a,
    const int &lda,const float &anorm,
    float &rcond,float *work,int *iwork,int &info);
  float F77NAME(sgeequ)(const int &m,const int &n,const float *A,
    const int &lda,float *r,float *c,float &rowcnd,float &rolcnd,
    float &amax,int &info);
  void F77NAME(sgeev)(const char &jobvl,const char &jobvr,const int &n,
    float *A,const int &lda,float *wr,float *wi,float *vl,
    const int &ldvl,float *vr,const int &ldvr,float *work,
    const int &lwork,int &info);
  void F77NAME(sgels)(const char &trans,const int &m,const int &n,
    const int &nrhs,float *a,const int &lda,float *b,const int &ldb,
    float *work,const int &lwork,const int &info);
  void F77NAME(sgemm)(const char &transa,const char &transb,
    const int &m,const int &n,const int &k, 
    const float &alpha,const float *A,
    const int &lda,const float *B,const int &ldb,
    const float &beta,float *C,const int &ldc);
  void F77NAME(sgemv)(const char &trans,const int &m,const int &n,
    const float &alpha,const float *a,const int &lda,const float *x,
    const int &incx,const float &beta,float *y,const int &incy);
  void F77NAME(sgelqf)(const int &m,const int &n,float *A,const int &lda,
    float *tau,float *work,int &lwork,int &info);
  void F77NAME(sgeqp3)(const int &m,const int &n,float *A,const int &lda,
    int *jpvt,float *tau,float *work,const int &lwork,int &info);
  void F77NAME(sgeqrf)(const int &m,const int &n,float *A,const int &lda,
    float *tau,float *work,int &lwork,int &info);
  void F77NAME(sger)(const int&,const int&,const float&,const float*,
    const int&,const float*,const int&,float*,const int&);
  void F77NAME(sgesv)(int &n,int &nrhs,float *a,int &lda,int *ipiv, 
    float *b, int &ldb, int &info);
  void F77NAME(sgesvd)(const char &jobu,const char &jobvt,const int &m, 
    const int &n,const float *A,const int &lda,
    float *S,float *U,const int &ldu,float *Vt,
    const int &ldvt,float *work,const int &lwork, 
    int &info);
  void F77NAME(sgesc2)(const int &n,const float *A,const int &lda,
    float *rhs,const int *ipiv,const int *jpiv,float &scale);
  void F77NAME(sgetc2)(const int &n,float *A,const int &lda,
    int *ipiv,int *jpiv,int &info);
  void F77NAME(sgetrf)(int &m,int &n,float *a,int &lda,int *ipiv,
    int &info);
  void F77NAME(sgetrs)(const char &trans,const int &n,const int &nrhs,
    float *a,const int &lda,int *ipiv,float *b,const int &ldb,
    int &info);
  void F77NAME(sgetri)(int &n,float *a,int &lda,int *ipiv,float *work,
    int &lwork,int &info);
  void F77NAME(sgtcon)(const char &norm,const int &n,const float *dl,
    const float *d,const float *du,const float *du2,const int *ipiv,
    const float &anorm,float &rcond,float *work,int *iwork,int &info);
  void F77NAME(sgtconnp)(const char &norm,const int &n,const float *dl,
    const float *d,const float *du,const float &anorm,float &rcond,
    float *work,int *iwork,int &info);
  void F77NAME(sgtmv)(const int &n,const float &alpha,const float *dl,
    const float *d,const float *du,const float *x,const int &incx,
    const float &beta,float *y,const int &incy);
  void F77NAME(sgtrfs)(const char &trans,const int &n,const int &nrhs,
    const float *dl,const float *d,const float *du,const float *dlf,
    const float *df,const float *duf,const float *du2,const int *ipiv,
    const float *B,const int &ldb,float *X,const int &ldx,float *ferr,
    float *berr,float *work,int *iwork,int &info);
  void F77NAME(sgtrfsnp)(const char &trans,const int &n,const int &nrhs,
    const float *dl,const float *d,const float *du,const float *dlf,
    const float *df,const float *duf,const float *B,const int &ldb,
    float *X,const int &ldx,float *ferr,float *berr,float *work,
    int *iwork,int &info);
  void F77NAME(sgtrfsr)(const char &trans,const int &n,const int &nrhs,
    const float *dl,const float *d,const float *du,const float *dlf,
    const float *df,const float *duf,const float *du2,const int *ipiv,
    const float *B,const int &ldb,float *X,const int &ldx,float *ferr,
    float *berr,float *work,int *iwork,int &info);
  void F77NAME(sgtrfsrnp)(const char &trans,const int &n,const int &nrhs,
    const float *dl,const float *d,const float *du,const float *dlf,
    const float *df,const float *duf,const float *B,const int &ldb,
    float *X,const int &ldx,float *ferr,float *berr,float *work,
    int *iwork,int &info);
  void F77NAME(sgtsv)(const int &n,const int &nrhs,const float *dl,
    const float *d,const float *du,float *b,const int &ldb,int &info);
  void F77NAME(sgtsvnp)(const int &n,const int &nrhs,const float *dl,
    const float *d,const float *du,float *b,const int &ldb,int &info);
  void F77NAME(sgtsvr)(const int &n,const int &nrhs,const float *dl,
    const float *d,const float *du,float *b,const int &ldb,int &info);
  void F77NAME(sgtsvrnp)(const int &n,const int &nrhs,const float *dl,
    const float *d,const float *du,float *b,const int &ldb,int &info);
  void F77NAME(sgttrf)(const int &n,float *dl,float *d,float *du,
    float *du2,int *ipiv,int &info);
  void F77NAME(sgttrfnp)(const int &n,float *dl,float *d,float *du,
    int &info);
  void F77NAME(sgttrs)(const char &trans,const int &n,const int &nrhs,
    float *dl,float *d,float *du,float *du2,int *ipiv,
    float *B,const int &ldb,int &info);
  void F77NAME(shseqr)(const char &job,const char &compz,const int &n,
    const int &ilo,const int &ihi,float *H,const int &ldh,float *wr,
    float *wi,float *Z,const int &ldz,float *work,const int &lwork,
    int &info);
  void F77NAME(slabad)(float &small,float &large);
  void F77NAME(slacn2)(const int &n,float *v,float *x,int *isgn,
    float &est,int &kase,int *isave);
  void F77NAME(slacpy)(const char &uplo,const int &m,const int &n,
    const float *A,const int &lda,float *B,const int &ldb);
  void F77NAME(sladiv)(const float&,const float&,const float&,
    const float&,float&,float&);
  void F77NAME(sgbamv)(const char &trans,const int &m,const int &n,
    const int &kl,const int &ku,const float &alpha,const float *A,
    const int &lda,const float *x,const int &incx,const float &beta,
    float *y,const int &incy);
  void F77NAME(slaic1)(const int &job,const int &j,const float *x,
    const float &sest,const float *w,const float &gamma,
    float &sestpr,float &s,float &c);
  float F77NAME(slamch)(const char&);
  float F77NAME(slangb)(const char &norm,const int &n,const int &kl,
    const int &lu,const float *AB,const int &ldab,float *work);
  float F77NAME(slange)(const char &norm,const int &m,const int &n,
    const float *A,const int &lda,float *work);
  float F77NAME(slangt)(const char &norm,const int &n,const float *dl,
    const float *d,const float *du);
  float F77NAME(slanhs)(const char &norm,const int &n,const float *A,
    const int &lda,float *work);
  float F77NAME(slansb)(const char &norm,const char &uplo,const int &n,
    const int &k,const float *AB,const int &ldab,float *work);
  float F77NAME(slanst)(const char &norm,const int &n,const float *d,
    const float *e);
  void F77NAME(slapmt)(const bool &forwrd,const int &m,const int &n,
    float *X,const int &ldx,int *k);
  void F77NAME(slaqgb)(const int &m,const int &n,const int &kl,
    const int &ku,float *AB,const int &ldab,const float *r,
    const float *c,const float &rowcnd,const float &colcnd,
    const float &amax,char &equed);
  void F77NAME(slaqsb)(const char &uplo,const int &n,const int &kd,
    float *AB,const int &ldab,const float *s,const float &scond,
    const float &amax,char &equed);
  void F77NAME(slaqsy)(const char &uplo,const int &n,float *A,
    const int &lda,const float *s,const float &scond,
    const float &amax,char &equed);
  void F77NAME(slaswp)(const int &n,float *A,const int &lda,
    const int &k1,const int &k2,const int *ipiv,const int &incx);
  float F77NAME(slansy)(const char &norm,const char &uplo,const int &n,
    const float *A,const int &lda,float *work);
  float F77NAME(slantb)(const char &norm,const char &uplo,
    const char &diag,const int &n,const int &k,const float *AB,
    const int &ldab,float *work);
  float F77NAME(slantr)(const char &norm,const char &uplo,
    const char &diag,const int &m,const int &n,const float *A,
    const int &lda,float *work);
  void F77NAME(slaqge)(const int &M,const int &n,float *A,const int &lda,
    const float *r,const float *c,const float &rowcnd,
    const float &colcnd,const float &amax,char &equed);
  void F77NAME(slargv)(const int &n,float *x,const int &incx,
    float *y,const int &incy,float *c,const int &incc);
  void F77NAME(slartv)(const int &n,float *x,const int &incx,
    float *y,const int &incy,const float *c,const float *s,
    const int &incc);
  void F77NAME(slascl)(const char &type,const int &kl,const int &ku,
    const float &cfrom,const float &cto,const int &m,const int &n,
    float *A,const int &lda,int &info);
  void F77NAME(slaset)(const char &uplo,const int &m,const int &n,
    const float &alpha,const float &beta,float *A,const int &lda);
  void F77_NAME(sla_geamv)(const int &trans,const int &m,const int &n,
    const float &alpha,const float *a,const int &lda,const float *x,
    const int &incx,const float &beta,float *y,const int &incy);
  float F77_NAME(sla_porpvgrw)(const char &uplo,const int &ncols,
    const float *A,const int &lda,const float *AF,const int &ldaf,
    float *work);
  void F77_NAME(sla_syamv)(const int &uplo,const int &n,
    const float &alpha,const float *A,const int &lda,
    const float *x,const int &ldx,const float &beta,
    float *y,const int &incy);
  float F77_NAME(sla_syrpvgrw)(const char &uplo,const int &n,
    const int &info,const float *A,const int &lda,const float *AF,
    const int &ldaf,const int *ipiv,float *work);
  float F77NAME(snrm2)(const int &n,float *x,const int &incx);
  void F77NAME(sorglq)(int &m,int &n,int &k,float *a,int &lda,
    const float *tau,float *work,int &lwork,
    int &info);
  void F77NAME(sorgqr)(const int &m,const int &n,const int &k,float *A,
    const int &lda,const float *tau,float *work,const int &lwork,
    int &info);
  void F77NAME(sormlq)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,float *A,const int &lda,float *tau,
    float *C,const int &ldc,float *work,int &lwork,int &info);
  void F77NAME(sormqr)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,float *A,const int &lda,float *tau,
    float *C,const int &ldc,float *work,int &lwork,int &info);
  void F77NAME(sormrz)(const char &side,const char &trans,const int &m,
    const int &n,const int &k,const int &l,const float *A,const int &lda,
    const float *tau,float *C,const int &ldc,float *work,
    const int &lwork,int &info);
  void F77NAME(spbcon)(const char &uplo,const int &n,const int &kd,
    float *AB,const int &ldab,const float &anorm,float &rcond,
    float *work,int *iwork,int &info);
  void F77NAME(spbequ)(const char &uplo,const int &n,const int &kd,
    const float *AB,const int &ldab,float *s,float &scond,
    float &amax,int &info);
  void F77NAME(spbsv)(const char &uplo,const int &n,const int &kd,
    const int &nrhs,float *AB,const int &ldab,float *B,const int &ldb,
    int &info);
  void F77NAME(spbtrf)(const char &uplo,const int &n,const int &kd,
    float *AB,const int &ldab,int &info);
  void F77NAME(spbtrs)(const char &uplo,const int &n,const int &kd,
    const int &nrhs,const float *AB,const int &ldab,float *B,
    const int &ldb,int &info);
  void F77NAME(spocon)(const char &uplo,const int &n,const float *a,
    const int &lda,const float &anorm,float &rcond,float *work,
    int *iwork,int &info);
//void F77NAME(spoequ)(const int &n,const float *A,const int &lda,
//  float *s,float &scond,float &amax,int &info);
  void F77NAME(spoequb)(const int &n,const float *A,const int &lda,
    float *s,float &scond,float &amax,int &info);
  void F77NAME(spotrf)(const char &uplo,const int &n,float *A,
    const int &lda,int &info);
  void F77NAME(spotri)(const char &uplo,const int &n,float *A,
    const int &lda,int &info);
  void F77NAME(spotrs)(const char &uplo,const int &n,const int &nrhs,
    const float *A,const int &lda,float *B,const int &ldb,int &info);
  void F77NAME(sptcon)(const int &n,const float *d,const float *e,
    const float &anorm,float &rcond,float *work,int &info);
  void F77NAME(sptsv)(const int &n,const int &nrhs,float *d,float *e,
    float *b,const int &ldb,int &info);
  void F77NAME(sptrfs)(const int &n,const int &nrhs,const float *d,
    const float *e,const float *df,const float *ef,const float *B,
    const int &ldb,float *X,const int &ldx,float *ferr,float *berr,
    float *work,int &info);
  void F77NAME(sptrfsr)(const int &n,const int &nrhs,const float *d,
    const float *e,const float *df,const float *ef,const float *B,
    const int &ldb,float *X,const int &ldx,float *ferr,float *berr,
    float *work,int &info);
  void F77NAME(spttrf)(const int &n,float *d,float *e,int &info);
  void F77NAME(spttrs)(const int &n,const int &nrhs,const float *d,
    const float *e, float *B,const int &ldb,int &info);
  void F77NAME(srot)(const int &n,float *sx,const int &incx,float *sy,
    const int &incy,const float &c,const float &s);
  void F77NAME(srotg)(float &sa,float &sb,float &c,float &s);
  void F77NAME(srscl)(const int &n,const float &sa,float *sx,
    const int &incx);
  void F77NAME(ssbev)(const char &jobz,const char &uplo,const int &n,
    const int &kd,float *AB,const int &ldab,float *w,float *Z,
    const int &ldz,float *work,int &info);
  void F77NAME(ssbamv)(const char &uplo,const int &n,const int &k,
    const float &alpha,const float *A,const int &lda,const float *x,
    const int &incx,const float &beta,float *y,const int &incy);
  void F77NAME(ssbmv)(const char &uplo,const int &n,const int &k,
    const float &alpha,const float *A,const int &lda,const float *x,
    const int &incx,const float &beta,float *y,const int &incy);
  void F77NAME(sscal)(const int &n,const float &a,float *x,
    const int &incx);
  void F77NAME(ssteqr)(const char &jobz,const int &n,float *d,float *e,
    float *Z,const int &ldz,float *work,int &info);
  void F77NAME(sstmv)(const int &n,const float &alpha,const float *dl,
    const float *d,const float *x,const int &incx,const float &beta,
    float *y,const int &incy);
  void F77NAME(sswap)(const int &n,float *sx,const int &incx,float *sy,
    const int &incy);
  void F77NAME(ssycon)(const char &uplo,const int &n,const float *a,
    const int &lda,const int *ipiv,const float &anorm,float &rcond,
    float *work,int *iwork,int &info);
  void F77NAME(ssyequb)(const char &uplo,const int &n,const float *A,
    const int &lda,float *s,float &scond,float &amax,float *work,
    int &info);
  void F77NAME(ssyev)(const char &jobz,const char &uplo,const int &n,
    float *A,const int &lda,float *W,float *work,int &lwork,
    int &info);
  void F77NAME(ssymm)(const char &side,const char &uplo,const int &m,
    const int &n,const float &alpha,const float *A,const int &lda,
    const float *B,const int &ldb,const float &beta,float *C,
    const int &ldc);
  void F77NAME(ssymv)(const char &uplo,const int &n,const float &alpha,
    const float *A,const int &lda,const float *x,const int &incx,
    const float &beta,float *y,const int &incy);
  void F77NAME(ssyr)(const char &uplo,const int &n,const float &alpha,
    const float *x,const int &incx,float *A,const int &lda);
  void F77NAME(ssyrk)(const char &uplo,const char &trans,const int &n,
    const int &k,const float &alpha,const float *A,const int &lda,
    const float &beta,float *C,const int &ldc);
  void F77NAME(ssyr2k)(const char &uplo,const char &trans,const int &n,
    const int &k,const float &alpha,const float *A,const int &lda,
    const float *B,const int &ldb,const float &beta,float *C,
    const int &ldc);
  void F77NAME(ssyr2)(const char &uplo,const int &n,const float &alpha,
    const float *x,const int &incx,const float *y,const int &incy,
    float *A,const int &lda);
  void F77NAME(ssytrf)(const char &uplo,const int &n,float *a,
    const int &lda,int *ipiv,float *work,const int &lwork,int &info);
  void F77NAME(ssytri)(const char &uplo,const int &n,float *a,
    const int &lda,int *ipiv,float *work,int &info);
  void F77NAME(ssytrs)(const char &uplo,const int &n,const int &nrhs,
    float *a,const int &lda,int *ipiv,float *b,const int &ldb,
    int &info);
  void F77NAME(stbsv)(const char &uplo,const char &trans,const char &diag,
    const int &n,const int &k,const float *A,const int &lda,float *x,
    const int &incx);
  void F77NAME(strcon)(const char &norm,const char &uplo,
    const char &diag,const int &n,const float *A,const int &lda,
    float &rcond,float *work,int *iwork,int &info);
  void F77NAME(strevc)(const char &side,const char &howmny,bool *select,
    const int &n,const float *T,const int &ldt,float *Vl,
    const int &ldvl,float *Vr,const int &ldvr,const int &mm,
    int &m,float *work,int &info);
  void F77NAME(strmv)(const char &uplo,const char &trans,
    const char &diag,const int &n,const float *A,const int &lda,
    float *x,const int &incx);
  void F77NAME(strsv)(const char&,const char&,const char&,const int &n,
    const float *A,const int &lda,float *x,const int &incx);
  void F77NAME(strmm)(const char &side,const char &uplo,
    const char &transa,const char &diag,const int &m,const int &n,
    const float &alpha,const float *A,const int &lda,float *B,
    const int &ldb);
  void F77NAME(strsm)(const char &side,const char &uplo,
    const char &transa,const char &diag,const int &m,const int &n,
    const float &alpha,const float *a,const int &lda,float *b,
    const int &ldb);
  void F77NAME(strtri)(const char &uplo,const char &diag,const int &n,
    float *A,const int &lda,int &info);
  void F77NAME(stzrzf)(const int &m,const int &n,float *A,const int &lda,
    float *tau,float *work,const int &lwork,int &info);
  float F77NAME(scnrm2)(const int &n,const complex<float> *x,
    const int &incx);
  int F77NAME(ilaslc)(const int &m,const int &n,const float *A,
    const int &lda);
  int F77NAME(ilaslr)(const int &m,const int &n,const float *A,
    const int &lda);
  int F77NAME(ilaenv)(const int &ispec,const char *name,const char *opts,
    const int &n1,const int &n2,const int &n3,const int &n4);
  int F77NAME(ilatrans)(const char &trans);
  int F77NAME(ilauplo)(const char &trans);
  void F77NAME(csscal)(const int &n,const float &da,
    complex<float> *zx,const int &incx);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Vector.C"
#ifdef DEBUG
template<> const float Vector<float,float>::undefined_ = HUGE_VAL;
#endif
template<> const float Vector<float,float>::zero_ = float_zero_;
template<> const float Vector<float,float>::one_ = float_one_;
template<> const float Vector<float,float>::mone_ = float_mone_;

template<> int Vector<float,float>::amax() const { 
  return F77NAME(isamax)(sz,data,1)-1;
}
template<> float Vector<float,float>::asum() const {
  return F77NAME(sasum)(sz,data,1);
}
template<> void Vector<float,float>::axpy(float a,
const Vector<float,float> &x) {
  F77NAME(saxpy)(min(sz,x.sz),a,x.data,1,data,1);
}
template<> float Vector<float,float>::dot(
const Vector<float,float> &x) const {
  return F77NAME(sdot)(min(sz,x.sz),x.data,1,data,1);
}
template<> float Vector<float,float>::dotc(
const Vector<float,float> &x) const {
  OBSOLETE(0);
}
template<> float Vector<float,float>::dotu(
const Vector<float,float> &x) const {
  OBSOLETE(0);
}
template<> float Vector<float,float>::nrm2() const {
  return F77NAME(snrm2)(sz,data,1);
}
template<> void Vector<float,float>::rot(Vector<float,float> &x,
float c,float s) {
  F77NAME(srot)(min(sz,x.sz),x.data,1,data,1,c,s);
}
template<> void Vector<float,float>::scal(float a) {
  F77NAME(sscal)(sz,a,data,1);
}
template<> void Vector<float,float>::swap(Vector<float,float> &x) {
  F77NAME(sswap)(min(sz,x.sz),x.data,1,data,1);
}
template<> void Vector<float,float>::largv(Vector<float,float> &y,
Vector<float,float> &c) {
  int n=size();
  CHECK_SAME(n,y.size());
  CHECK_SAME(n,c.size());
  F77NAME(slargv)(n,data,1,y.data,1,c.data,1);
}
template<> void Vector<float,float>::lartv(Vector<float,float> &y,
const Vector<float,float> &c,const Vector<float,float> &s) {
  int n=size();
  CHECK_SAME(n,y.size());
  CHECK_SAME(n,c.size());
  CHECK_SAME(n,s.size());
  F77NAME(slartv)(n,addr(),1,y.addr(),1,c.addr(),s.addr(),1);
}

template class Vector<float,float>; 
template void testVector(const float&,const float&);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "Matrix.C"

#ifdef DEBUG
template<> const float Matrix<float,float>::undefined_ = HUGE_VAL;
#endif
template<> const float Matrix<float,float>::zero_ = float_zero_;
template<> const float Matrix<float,float>::one_ = float_one_;
template<> const float Matrix<float,float>::mone_ = float_mone_;

/*
template<> Matrix<float,float>* Matrix<float,float>::transpose()
const {
  int m=size(0),n=size(1);
  Matrix<float,float> *X=OPERATOR_NEW Matrix<float,float>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(m,addr(0,j),1,X->addr(j,0),n);
  }
  return X;
}

template<> Matrix<float,float>*
Matrix<float,float>::conjugateTranspose() const {
  return transpose();
}
*/

template<> void Matrix<float,float>::interchangeColumns(int i,int j) {
  int m=size(0),n=size(1);
  CHECK_BOUNDS(i,0,n)
  CHECK_BOUNDS(j,0,n)
  F77NAME(sswap)(m,addr(0,i),1,addr(0,j),1);
}

template<> void Matrix<float,float>::interchangeRows(int i,int j) {
  int m=size(0),n=size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,m)
  F77NAME(sswap)(n,addr(i,0),m,addr(j,0),m);
}

template<> void Matrix<float,float>::gemv(float alpha,
const Vector<float,float> &x,float beta,Vector<float,float> &y,
char trans) const {
//TRACER_CALL(t,"Matrix::gemv");
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(x.size(),m);
    CHECK_SAME(y.size(),n);
    F77NAME(sgemv)('T',m,n,alpha,addr(),m,x.addr(),1,beta,y.addr(),1);
  } else {
    CHECK_SAME(x.size(),n);
    CHECK_SAME(y.size(),m);
    F77NAME(sgemv)('N',m,n,alpha,addr(),m,x.addr(),1,beta,y.addr(),1);
  }
}

template<> void Matrix<float,float>::ger(float alpha,
const Vector<float,float> &x,const Vector<float,float> &y) {
  int m=size(0),n=size(1);
  CHECK_SAME(x.size(),m);
  CHECK_SAME(y.size(),n);
  F77NAME(sger)(m,n,alpha,x.addr(),1,y.addr(),1,addr(),m);
}

template<> void Matrix<float,float>::gerc(float alpha,
const Vector<float,float> &x,const Vector<float,float> &y) {
  OBSOLETE("inappropriate for this class");
}

template<> void Matrix<float,float>::geru(float alpha,
const Vector<float,float> &x,const Vector<float,float> &y) {
  OBSOLETE("inappropriate for this class");
}

template<> void Matrix<float,float>::gemm(float alpha,
const Matrix<float,float> &A,const Matrix<float,float> &B,
float beta,char transa,char transb) {
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
  F77NAME(sgemm)(transa,transb,m,n,k,alpha,A.addr(),A.size(0),
    B.addr(),B.size(0),beta,addr(),m);
}

template<> float Matrix<float,float>::equilibrate(
Vector<float,float> &r,Vector<float,float> &c,float &rowcnd,
float &colcnd) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,r.size());
  CHECK_SAME(n,c.size());
  float amax=float_undefined_;
  int info;
  F77NAME(sgeequ)(m,n,addr(),m,r.addr(),c.addr(),rowcnd,colcnd,
    amax,info);
  CHECK_SAME(info,0);
  return amax;
}

template<> void Matrix<float,float>::copyFrom(char uplo,int m,int n,
const Matrix<float,float> &A) {
  int s0=size(0),as0=A.size(0);
  m=min(m,min(s0,as0));
  n=min(n,min(size(1),A.size(1)));
  if (uplo=='A' || uplo=='a') {
    F77NAME(slacpy)(uplo,m,n,A.addr(),as0,addr(),s0);
//  otherwise, slacpy only copies a triangular part:
  } else if (uplo=='L' || uplo=='l') {
    for (int j=0;j<n;j++) {
      if (j<m) {
        F77NAME(scopy)(m-j,A.addr(j,j),1,addr(j,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+1,m),A.addr(0,j),1,addr(0,j),1);
    }
  }
}

template<> void Matrix<float,float>::scale(char type,int kl,int ku,
float denominator,float numerator) {
  int m=size(0),info;
  F77NAME(slascl)(type,kl,ku,denominator,numerator,m,size(1),addr(),
    m,info);
  CHECK_SAME(info,0);
}

template<> void Matrix<float,float>::set(char uplo,float offdiag,
float diag) {
  int m=size(0);
  F77NAME(slaset)(uplo,m,size(1),offdiag,diag,addr(),m);
}

template<> int Matrix<float,float>::lastNonzeroColumn() const {
  int m=size(0);
  return F77NAME(ilaslc)(m,size(1),addr(),m);
}

template<> int Matrix<float,float>::lastNonzeroRow() const {
  int m=size(0);
  return F77NAME(ilaslr)(m,size(1),addr(),m);
}

template<> float Matrix<float,float>::normFrobenius() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slange)('F',m,size(1),addr(),m,work);
}

template<> float Matrix<float,float>::normInfinity() const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float val=F77NAME(slange)('I',m,size(1),addr(),m,work);
  delete [] work; work=0;
  return val;
}

template<> float Matrix<float,float>::normMaxEntry() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slange)('M',m,size(1),addr(),m,work);
}

template<> float Matrix<float,float>::normOne() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slange)('1',m,size(1),addr(),m,work);
}

template<> Matrix<float,float>* Matrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(M.size(0),k);
  Matrix<float,float> *P=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(sgemm)('N','N',m,n,k,float_one_,addr(),m,M.addr(),k,
    float_zero_,P->addr(),m);
  return P;
}

template<> Vector<float,float>* Matrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(v.size(),n);
  Vector<float,float> *P=OPERATOR_NEW Vector<float,float>(m);
  F77NAME(sgemv)('N',m,n,float_one_,addr(),m,v.addr(),1,float_zero_,
    P->addr(),1);
  return P;
}

template<> void Matrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  int m=size(0), n=size(1);
  Matrix<float,float> AF(*this);
  int info;
  if (m==n) {
    CHECK_SAME(m,b.size())
    CHECK_SAME(n,x.size())
    x.copy(b);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(sgetrf)(m,m,AF.addr(),m,ipiv,info);
    if (info==0) {
      F77NAME(sgetrs)(trans,m,1,AF.addr(),m,ipiv,x.addr(),m,info);
    }
    CHECK_SAME(info,0)
    delete [] ipiv;
  } else {
    float w=HUGE_VAL;
    int lwork=-1;
    bool transposed=(trans!='N' && trans!='n');
    if (m>n) {
      if (transposed) {
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())
        x.copyFrom(n,b);
        F77NAME(sgels)(trans,m,n,1,AF.addr(),m,x.addr(),m,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        float *work=OPERATOR_NEW_BRACKET(float,lwork);
        F77NAME(sgels)(trans,m,n,1,AF.addr(),m,x.addr(),m,work,
          lwork,info);
        CHECK_SAME(info,0)
        delete [] work;
      } else {
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        Vector<float,float> xtmp(m);
        xtmp.copy(b);
        F77NAME(sgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),m,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        float *work=OPERATOR_NEW_BRACKET(float,lwork);
        F77NAME(sgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),m,work,
          lwork,info);
        CHECK_SAME(info,0)

        x.copyFrom(n,xtmp);
        delete [] work;
      }
    } else {
      if (transposed) {
        CHECK_SAME(n,b.size())
        CHECK_SAME(m,x.size())
        Vector<float,float> xtmp(n);
        xtmp.copy(b);
        F77NAME(sgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),n,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        float *work=OPERATOR_NEW_BRACKET(float,lwork);
        F77NAME(sgels)(trans,m,n,1,AF.addr(),m,xtmp.addr(),n,work,
          lwork,info);
        CHECK_SAME(info,0)

        x.copyFrom(m,xtmp);
        delete [] work;
      } else {
        CHECK_SAME(m,b.size())
        CHECK_SAME(n,x.size())
        x.copyFrom(m,b);
        F77NAME(sgels)(trans,m,n,1,AF.addr(),m,x.addr(),n,
          &w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        float *work=OPERATOR_NEW_BRACKET(float,lwork);
        F77NAME(sgels)(trans,m,n,1,AF.addr(),m,x.addr(),n,work,
          lwork,info);
        CHECK_SAME(info,0)
        delete [] work;
      }
    }
  }
}

template<> void Matrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  int m=size(0), n=size(1);
  Matrix<float,float> AF(*this);
  int info;
  bool left_side=(side=='L' || side=='l');
  bool transposed=(trans!='N' && trans!='n');
  if (m==n) {
    X.copy(B);
    int *ipiv=OPERATOR_NEW_BRACKET(int,m);
    F77NAME(sgetrf)(m,m,AF.addr(),m,ipiv,info);
    if (info==0) {
      if (left_side) {
        int k=B.size(1);
        CHECK_SAME(k,X.size(1))
        CHECK_SAME(m,B.size(0))
        CHECK_SAME(n,X.size(0))
        F77NAME(sgetrs)(trans,m,k,AF.addr(),m,ipiv,X.addr(),m,info);
      } else {
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        CHECK_SAME(m,B.size(1))
        CHECK_SAME(n,X.size(1))
        if (transposed) {
          for (int j=0;j<m-1;j++) {
            if (ipiv[j]-1!=j) {
              F77NAME(sswap)(k,X.addr(0,j),1,X.addr(0,ipiv[j]-1),1);
            }
          }
          F77NAME(strsm)('R','L','T','U',k,m,float_one_,AF.addr(),m,
            X.addr(),k);
          F77NAME(strsm)('R','U','T','N',k,m,float_one_,AF.addr(),m,
            X.addr(),k);
        } else {
          F77NAME(strsm)('R','U','N','N',k,m,float_one_,AF.addr(),m,
            X.addr(),k);
          F77NAME(strsm)('R','L','N','U',k,m,float_one_,AF.addr(),m,
            X.addr(),k);
          for (int j=m-2;j>=0;j--) {
            if (ipiv[j]-1!=j) {
              F77NAME(sswap)(k,X.addr(0,j),1,X.addr(0,ipiv[j]-1),1);
            }
          }
        }
      }
    }
    CHECK_SAME(info,0)
    delete [] ipiv; ipiv=0;
  } else {
    float w=HUGE_VAL;
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
          F77NAME(sgels)(trans,m,n,k,AF.addr(),m,X.addr(),m,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w);
          float *work=OPERATOR_NEW_BRACKET(float,lwork);
          F77NAME(sgels)(trans,m,n,k,AF.addr(),m,X.addr(),m,work,
            lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        } else {
          CHECK_SAME(m,B.size(0))
          CHECK_SAME(n,X.size(0))
          Matrix<float,float> Xtmp(B);
          F77NAME(sgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),m,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w);
          float *work=OPERATOR_NEW_BRACKET(float,lwork);
          F77NAME(sgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),m,work,
            lwork,info);
          CHECK_SAME(info,0)

          X.copyFrom('A',n,k,Xtmp);
          delete [] work; work=0;
        }
      } else { // m>n, right side
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        float *tau=OPERATOR_NEW_BRACKET(float,n);
        F77NAME(sgeqrf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)

        lwork=static_cast<int>(w);
        float *work=OPERATOR_NEW_BRACKET(float,lwork);
        F77NAME(sgeqrf)(m,n,AF.addr(),m,tau,work,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        if (transposed) {
          CHECK_SAME(m,B.size(1))
          CHECK_SAME(n,X.size(1))
          Matrix<float,float> Xtmp(B);

          lwork=-1;
          F77NAME(sormqr)('R','N',k,m,n,AF.addr(),m,tau,Xtmp.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w);
          work=OPERATOR_NEW_BRACKET(float,lwork);
          F77NAME(sormqr)('R','N',k,m,n,AF.addr(),m,tau,Xtmp.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
          X.copyFrom('A',k,n,Xtmp);

          F77NAME(strsm)('R','U','T','N',k,n,float_one_,AF.addr(),m,
            X.addr(),k);
        } else {
          CHECK_SAME(n,B.size(1))
          CHECK_SAME(m,X.size(1))
          X=zero_;
          X.copyFrom('A',k,n,B);
          F77NAME(strsm)('R','U','N','N',k,n,float_one_,AF.addr(),m,
            X.addr(),k);
          lwork=-1;
          F77NAME(sormqr)('R','T',k,m,n,AF.addr(),m,tau,X.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w);
          work=OPERATOR_NEW_BRACKET(float,lwork);
          F77NAME(sormqr)('R','T',k,m,n,AF.addr(),m,tau,X.addr(),k,
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
          Matrix<float,float> Xtmp(B);
          F77NAME(sgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),n,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w);
          float *work=OPERATOR_NEW_BRACKET(float,lwork);
          F77NAME(sgels)(trans,m,n,k,AF.addr(),m,Xtmp.addr(),n,work,
            lwork,info);
          CHECK_SAME(info,0)

          X.copyFrom('A',m,k,Xtmp);
          delete [] work; work=0;
        } else {
          CHECK_SAME(m,B.size(0))
          CHECK_SAME(n,X.size(0))
          X.copyFrom('A',m,k,B);
          F77NAME(sgels)(trans,m,n,k,AF.addr(),m,X.addr(),n,
            &w,lwork,info);
          CHECK_SAME(info,0)

          lwork=static_cast<int>(w);
          float *work=OPERATOR_NEW_BRACKET(float,lwork);
          F77NAME(sgels)(trans,m,n,k,AF.addr(),m,X.addr(),n,work,
            lwork,info);
          CHECK_SAME(info,0)
          delete [] work;
        }
      } else { // m<n, right side
        int k=B.size(0);
        CHECK_SAME(k,X.size(0))
        float *tau=OPERATOR_NEW_BRACKET(float,n);
        F77NAME(sgelqf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)
        lwork=static_cast<int>(w);
        float *work=OPERATOR_NEW_BRACKET(float,lwork);
        F77NAME(sgelqf)(m,n,AF.addr(),m,tau,&w,lwork,info);
        CHECK_SAME(info,0)
        delete [] work; work=0;
        if (transposed) {
          CHECK_SAME(m,B.size(1))
          CHECK_SAME(n,X.size(1))
          X=zero_;
          X.copyFrom('A',k,m,B);
          F77NAME(strsm)('R','L','T','N',k,m,float_one_,AF.addr(),m,
            X.addr(),k);
          lwork=-1;
          F77NAME(sormlq)('R','N',k,n,m,AF.addr(),m,tau,X.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w);
          work=OPERATOR_NEW_BRACKET(float,lwork);
          F77NAME(sormlq)('R','N',k,n,m,AF.addr(),m,tau,X.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
        } else {
          CHECK_SAME(n,B.size(1))
          CHECK_SAME(m,X.size(1))
          Matrix<float,float> Xtmp(B);
          lwork=-1;
          F77NAME(sormlq)('R','T',k,n,m,AF.addr(),m,tau,Xtmp.addr(),k,
            &w,lwork,info);
          CHECK_SAME(info,0)
          lwork=static_cast<int>(w);
          work=OPERATOR_NEW_BRACKET(float,lwork);
          F77NAME(sormlq)('R','T',k,n,m,AF.addr(),m,tau,Xtmp.addr(),k,
            work,lwork,info);
          CHECK_SAME(info,0)
          delete [] work; work=0;
          X.copyFrom('A',k,m,Xtmp);

          F77NAME(strsm)('R','L','N','N',k,m,float_one_,AF.addr(),m,
            X.addr(),k);
        }
        delete [] tau; tau=0;
      }
    }
  }
}

template class Matrix<float,float>; 
template void testMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SquareMatrix.C"

template<> SquareMatrix<float,float>*
SquareMatrix<float,float>::operator*(
const SquareMatrix<float,float> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  F77NAME(sgemm)('N','N',n,n,n,float_one_,addr(),n,S.addr(),n,
    float_zero_,P->addr(),n);
  return P;
}

template<> Matrix<float,float>*
SquareMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *P=(m==n ?
    OPERATOR_NEW SquareMatrix<float,float>(n) :
    OPERATOR_NEW Matrix<float,float>(m,n));
  F77NAME(sgemm)('N','N',m,n,m,float_one_,addr(),m,M.addr(),m,
    float_zero_,P->addr(),m);
  return P;
}

/*
template<> SquareMatrix<float,float>*
SquareMatrix<float,float>::transpose() const {
  int n=size(0);
  SquareMatrix<float,float> *X=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n,addr(0,j),1,X->addr(j,0),n);
  }
  return X;
}

template<> SquareMatrix<float,float>*
SquareMatrix<float,float>::conjugateTranspose() const {
  return transpose();
}
*/

template<> float 
SquareMatrix<float,float>::reciprocalConditionNumber(char norm) const {
  int n=size(0),info;
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,4*n);
  float anorm=F77NAME(slange)(norm,n,n,addr(),n,work);

  SquareMatrix<float,float> *AF=
    OPERATOR_NEW SquareMatrix<float,float>(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(sgetrf)(n,n,AF->addr(0,0),n,ipiv,info);
  CHECK_SAME(info,0)
  delete [] ipiv; ipiv=0;

  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(sgecon)(norm,n,AF->addr(),n,anorm,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete[] work;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SquareMatrix<float,float>*
SquareMatrix<float,float>::inverse() const {
  int n=size(0);
  SquareMatrix<float,float> *Ainv=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  *Ainv = *this;
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(sgetrf)(n,n,Ainv->addr(0,0),n,ipiv,info);
  CHECK_SAME(info,0)

  float w=HUGE_VAL;
  int lwork=-1;
  F77NAME(sgetri)(n,Ainv->addr(0,0),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(sgetri)(n,Ainv->addr(0,0),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)

  delete [] ipiv;
  delete [] work;
  return Ainv;
}
*/

template<> Vector<float,complex<float> >*
SquareMatrix<float,float>::eigenvalues(
SquareMatrix<float,complex<float> > *&V,
SquareMatrix<float,complex<float> > *&U) const {
  int n=size(0);
  if (V!=0) CHECK_SAME(n,V->size(0));
  if (U!=0) CHECK_SAME(n,U->size(0));
  char jobvl=(V==0 ? 'N' : 'V');
  char jobvr=(U==0 ? 'N' : 'V');
  float *wr=OPERATOR_NEW_BRACKET(float,n);
  float *wi=OPERATOR_NEW_BRACKET(float,n);
  float *vl=(V==0 ? 0 : OPERATOR_NEW_BRACKET(float,n*n));
  float *vr=(U==0 ? 0 : OPERATOR_NEW_BRACKET(float,n*n));
  SquareMatrix<float,float> AF(n);
  AF.copy(*this);
  float w;
  int lwork=-1,info;
  F77NAME(sgeev)(jobvl,jobvr,n,AF.addr(),n,wr,wi,vl,n,vr,n,&w,
    lwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(sgeev)(jobvl,jobvr,n,AF.addr(),n,wr,wi,vl,n,vr,n,work,
    lwork,info);
  CHECK_TEST(info==0);

  Vector<float,complex<float> > *lambda =
    OPERATOR_NEW Vector<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    complex<float> &ev=(*lambda)[j];
    ev.real()=wr[j];
    ev.imag()=wi[j];
  }
  for (int j=0;j<n;) {
    int jn=j*n;
    if (wi[j]>zero_) {
      int jp1n=jn+n;
      if (V!=0) {
        for (int i=0;i<n;i++) {
          complex<float> &z=V->operator()(i,j);
          z.real()=vl[i+jn];
          z.imag()=vl[i+jp1n];
          complex<float> &zz=V->operator()(i,j+1);
          zz.real()=vl[i+jn];
          zz.imag()=-vl[i+jp1n];
        }
      }
      if (U!=0) {
        for (int i=0;i<n;i++) {
          complex<float> &z=U->operator()(i,j);
          z.real()=vr[i+jn];
          z.imag()=vr[i+jp1n];
          complex<float> &zz=U->operator()(i,j+1);
          zz.real()=vr[i+jn];
          zz.imag()=-vr[i+jp1n];
        }
      }
      j+=2;
    } else {
      if (V!=0) {
        for (int i=0;i<n;i++) {
          complex<float> &z=V->operator()(i,j);
          z.real()=vl[i+jn];
          z.imag()=zero_;
        }
      }
      if (U!=0) {
        for (int i=0;i<n;i++) {
          complex<float> &z=U->operator()(i,j);
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

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,const SquareMatrix<float,float> &S) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<float,float> *P=(m==n ?
    OPERATOR_NEW SquareMatrix<float,float>(n) :
    OPERATOR_NEW Matrix<float,float>(m,n));
  F77NAME(sgemm)('N','N',m,n,n,float_one_,M.addr(),m,S.addr(),n,
    float_zero_,P->addr(),m);
  return P;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "TrapezoidalMatrix.C"
template<> float TrapezoidalMatrix<float,float>::outofbounds_ =
  float_zero_;
template<> float TrapezoidalMatrix<float,float>::safety_ =
  float_zero_;

template class TrapezoidalMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> Matrix<float,float>*
LowerTrapezoidalMatrix<float,float>::makeMatrix() const {
  int m=size(0),n=size(1);
  Matrix<float,float> *M=OPERATOR_NEW
    Matrix<float,float>(m,n,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(this)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(m-j,addr(j,j),1,M->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*M)(j,j)=float_one_;
      if (j+1<m) F77NAME(scopy)(m-j-1,addr(j+1,j),1,M->addr(j+1,j),1);
    }
  }
  return M;
}

template<> LowerTrapezoidalMatrix<float,float>& 
LowerTrapezoidalMatrix<float,float>::operator+=(
const LowerTrapezoidalMatrix<float,float> &L) {
  bool L_non_unit=
    (dynamic_cast<const UnitLowerTrapezoidalMatrix<float,float>*>(&L)
    ==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0))
  CHECK_SAME(n,L.size(1))
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(m-j,one_,L.addr(j,j),1,addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(saxpy)(m-j-1,one_,L.addr(j+1,j),1,addr(j+1,j),1);
      }
      (*this)(j,j)+=float_one_;
    }
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<float,float>&
LowerTrapezoidalMatrix<float,float>::operator-=(
const LowerTrapezoidalMatrix<float,float> &L) {
  bool L_non_unit=
    (dynamic_cast<const UnitLowerTrapezoidalMatrix<float,float>*>(&L)
    ==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0))
  CHECK_SAME(n,L.size(1))
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(m-j,mone_,L.addr(j,j),1,addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(saxpy)(m-j-1,mone_,L.addr(j+1,j),1,addr(j+1,j),1);
      }
      (*this)(j,j)-=float_one_;
    }
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<float,float>& 
LowerTrapezoidalMatrix<float,float>::operator*=(float d) {
  int m=size(0),n=size(1);
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(m-j,d,addr(j,j),1);
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<float,float>&
LowerTrapezoidalMatrix<float,float>::operator/=(float d) {
  int m=size(0),n=size(1);
  float dinv=float_one_/d;
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(m-j,dinv,addr(j,j),1);
  }
  return *this;
}

template<> LowerTrapezoidalMatrix<float,float>*
LowerTrapezoidalMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<float,float> *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(m,n);
  S->copy(*this);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(m-j,float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(m-j-1,float_one_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> Matrix<float,float>*
LowerTrapezoidalMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *S=OPERATOR_NEW Matrix<float,float>(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(m-j,float_one_,addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> LowerTrapezoidalMatrix<float,float>*
LowerTrapezoidalMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  bool L_non_unit=
   (dynamic_cast<const UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<float,float> *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(m,n);
  S->copy(*this);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(m-j,float_mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(m-j-1,float_mone_,L.addr(j+1,j),1,
        S->addr(j+1,j),1);
      (*S)(j,j)-=float_one_;
    }
  }
  return S;
}

template<> Matrix<float,float>*
LowerTrapezoidalMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *D=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(slaset)('U',m,n,float_zero_,float_zero_,D->addr(),m);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(m-j,addr(j,j),1,D->addr(j,j),1);
    F77NAME(saxpy)(m,float_mone_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> LowerTrapezoidalMatrix<float,float>* 
LowerTrapezoidalMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &L) const {
// to compute the jth column of the product
//   [ L_11   0  ] [ 0 ] = [    0   ]
//   [ L_21 L_22 ] [ m ] = [ L_22 m ]
//   [ L_31 L_32 ]       = [ L_32 m ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  LowerTrapezoidalMatrix<float,float> *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(m,n);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(k-j,L.addr(j,j),1,P->addr(j,j),1);
      F77NAME(strmv)('L','N','N',k-j,addr(j,j),m,P->addr(j,j),1);
      if (m>k) {
        F77NAME(sgemv)('N',m-k,k-j,float_one_,addr(k,j),m,L.addr(j,j),1,
          float_zero_,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*P)(j,j)=(*this)(j,j);
      if (j<k-1) {
        F77NAME(scopy)(k-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
        F77NAME(strmv)('L','N','N',k-j-1,addr(j+1,j+1),m,
          P->addr(j+1,j),1);
        F77NAME(saxpy)(k-j-1,float_one_,addr(j+1,j),1,P->addr(j+1,j),1);
      }
      if (m>k) {
        if (j<k-1) {
          F77NAME(sgemv)('N',m-k,k-j-1,float_one_,addr(k,j+1),m,
            L.addr(j+1,j),1,float_zero_,P->addr(k,j),1);
        }
        F77NAME(saxpy)(m-k,float_one_,addr(k,j),1,P->addr(k,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<float,float>*
LowerTrapezoidalMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<float,float> *P=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(slacpy)('A',k,n,M.addr(),k,P->addr(),m);
  F77NAME(strmm)('L','L','N','N',k,n,float_one_,addr(),m,
    P->addr(),m);
  if (m>k) {
    F77NAME(sgemm)('N','N',m-k,n,k,float_one_,addr(k,0),m,
      M.addr(),k,float_zero_,P->addr(k,0),m);
  }
  return P;
}

template<> Vector<float,float>*
LowerTrapezoidalMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<float,float> *p=OPERATOR_NEW Vector<float,float>(m);
  F77NAME(scopy)(n,v.addr(),1,p->addr(),1);
  F77NAME(strmv)('L','N','N',n,addr(),m,p->addr(),1);
  if (m>n) {
    F77NAME(sgemv)('N',m-n,n,float_one_,addr(n,0),m,v.addr(),1,
      float_zero_,p->addr(n),1);
  }
  return p;
}

/*
template<> UpperTrapezoidalMatrix<float,float>*
LowerTrapezoidalMatrix<float,float>::transpose() const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<float,float> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(m-j,addr(j,j),1,U->addr(j,j),n);
  }
  return U;
}

template<> UpperTrapezoidalMatrix<float,float>*
LowerTrapezoidalMatrix<float,float>::conjugateTranspose() const {
  return transpose();
}
*/

template<> Vector<float,float>*
LowerTrapezoidalMatrix<float,float>::trmv(
const Vector<float,float> &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<float,float> *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<float,float>(n);
    F77NAME(scopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(strmv)('L','T','N',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(sgemv)('T',m-n,n,one_,addr(n,0),m,x.addr(n),1,one_,
        p->addr(),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<float,float>(m);
    F77NAME(scopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(strmv)('L','N','N',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(sgemv)('N',m-n,n,one_,addr(n,0),m,x.addr(),1,zero_,
        p->addr(n),1);
    }
  }
  return p;
}

template<> Matrix<float,float>*
LowerTrapezoidalMatrix<float,float>::trmm(
const Matrix<float,float> &X,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<float,float> *P=0;
  if (side=='L' || side=='l') {
    int k=X.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,X.size(0));
      P=OPERATOR_NEW Matrix<float,float>(n,k);
      F77NAME(slacpy)('A',n,k,X.addr(),m,P->addr(),n);
      F77NAME(strmm)('L','L','T','N',n,k,one_,addr(),m,P->addr(),n);
      if (m>n) {
        F77NAME(sgemm)('T','N',n,k,m-n,one_,addr(n,0),m,X.addr(n,0),m,
          one_,P->addr(),n);
      }
    } else {
      CHECK_SAME(n,X.size(0));
      P=OPERATOR_NEW Matrix<float,float>(m,k);
      F77NAME(slacpy)('A',n,k,X.addr(),n,P->addr(),m);
      F77NAME(strmm)('L','L','N','N',n,k,one_,addr(),m,P->addr(),m);
      if (m>n) {
        F77NAME(sgemm)('N','N',m-n,k,n,one_,addr(n,0),m,X.addr(),n,
          zero_,P->addr(n,0),m);
      }
    }
  } else {
    int k=X.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,X.size(1));
      P=OPERATOR_NEW Matrix<float,float>(k,m);
      F77NAME(slacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(strmm)('R','L','T','N',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(sgemm)('N','T',k,m-n,n,one_,X.addr(),k,addr(n,0),m,
          zero_,P->addr(0,n),k);
      }
    } else {
      CHECK_SAME(m,X.size(1));
      P=OPERATOR_NEW Matrix<float,float>(k,n);
      F77NAME(slacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(strmm)('R','L','N','N',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(sgemm)('N','N',k,n,m-n,one_,X.addr(0,n),k,addr(n,0),m,
          one_,P->addr(),k);
      }
    }
  }
  return P;
}

template<> void LowerTrapezoidalMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(scopy)(n,b.addr(),1,x.addr(),1);
    if (m>n) { // use trailing entries of x as free variables
      F77NAME(sgemv)('T',m-n,n,mone_,addr(n,0),m,x.addr(n),1,one_,
        x.addr(),1);
    }
    F77NAME(strsv)('L','T','N',n,addr(),m,x.addr(),1);
  } else { // calling routine will have to check consistency conditions
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(scopy)(n,b.addr(),1,x.addr(),1);
    F77NAME(strsv)('L','N','N',n,addr(),m,x.addr(),1);
  }
}

template<> void LowerTrapezoidalMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(slacpy)('A',n,k,B.addr(),n,X.addr(),n);
      if (m>n) { // use these entries of X as free variables
        F77NAME(sgemm)('T','N',n,k,m-n,mone_,addr(n,0),m,X.addr(n,0),m,
          one_,X.addr(),m);
      }
      F77NAME(strsm)('L','L','T','N',n,k,one_,addr(),m,X.addr(),m);
    } else { // calling routine will have to check consistency conditions
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(slacpy)('A',n,k,B.addr(),m,X.addr(),n);
      F77NAME(strsm)('L','L','N','N',n,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(slacpy)('A',k,n,B.addr(),k,X.addr(),k);
      F77NAME(strsm)('R','L','T','N',k,n,one_,addr(),m,X.addr(),k);
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(slacpy)('A',k,n,B.addr(),k,X.addr(),k);
      if (m>n) { // use these entries of X as free variables
        F77NAME(sgemm)('N','N',k,n,m-n,mone_,X.addr(0,n),n,addr(n,0),m,
          one_,X.addr(),k);
      }
      F77NAME(strsm)('R','L','N','N',k,n,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> float LowerTrapezoidalMatrix<float,float>::normFrobenius()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('F','L','N',m,size(1),addr(),m,work);
}

template<> float LowerTrapezoidalMatrix<float,float>::normInfinity()
const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float result=F77NAME(slantr)('I','L','N',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> float LowerTrapezoidalMatrix<float,float>::normMaxEntry()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('M','L','N',m,size(1),addr(),m,work);
}

template<> float LowerTrapezoidalMatrix<float,float>::normOne()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('O','L','N',m,size(1),addr(),m,work);
}

template<> Matrix<float,float>* operator+(
const Matrix<float,float> &M,
const LowerTrapezoidalMatrix<float,float> &L) {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *S=OPERATOR_NEW Matrix<float,float>(m,n);
  S->copy(M);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(m-j,float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(saxpy)(m-j-1,float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> Matrix<float,float>* operator-(
const Matrix<float,float> &M,
const LowerTrapezoidalMatrix<float,float> &L) {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *D=OPERATOR_NEW Matrix<float,float>(m,n);
  D->copy(M);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(m-j,float_mone_,L.addr(j,j),1,D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<m-1) {
        F77NAME(saxpy)(m-j-1,float_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
      (*D)(j,j)-=float_one_;
    }
  }
  return D;
}

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,
const LowerTrapezoidalMatrix<float,float> &L) {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  char diagL=
    (dynamic_cast<const UnitLowerTrapezoidalMatrix<float,float>*>(&L)
    ==0 ? 'N' : 'U');
  int m=M.size(0),k=L.size(0),n=L.size(1);
  CHECK_SAME(k,M.size(1));
  Matrix<float,float> *P=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(slacpy)('A',m,n,M.addr(),m,P->addr(),m);
  F77NAME(strmm)('R','L','N',diagL,m,n,float_one_,L.addr(),k,
    P->addr(),m);
  if (k>n) {
    F77NAME(sgemm)('N','N',m,n,k-n,float_one_,M.addr(0,n),m,
      L.addr(n,0),k,float_one_,P->addr(),m);
  }
  return P;
}

template class LowerTrapezoidalMatrix<float,float>;
template void testLowerTrapezoidalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float
LowerTriangularMatrix<float,float>::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  F77NAME(strcon)(norm,'L','N',n,addr(),n,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete [] work;
  return rcond;
}

/*
template<> LowerTriangularMatrix<float,float>*
LowerTriangularMatrix<float,float>::inverse() const {
  LowerTriangularMatrix<float,float> *L=
    OPERATOR_NEW LowerTriangularMatrix<float,float>(*this);
  int n=size(0);
  int info;
  F77NAME(strtri)('L','N',n,L->addr(),n,info);
  CHECK_TEST(info==0);
  return L;
}
*/

template class LowerTriangularMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void
UnitLowerTrapezoidalMatrix<float,float>::copyFrom(int m,int n,
const Matrix<float,float> &L) {
  m=min(m,min(size(0),L.size(0)));
  n=min(n,min(size(1),L.size(1)));
  for (int j=0;j<n;j++) {
    if (j+1<m) {
      F77NAME(scopy)(m-j-1,L.addr(j+1,j),1,addr(j+1,j),1);
    }
  }
}

/*
template<> UnitUpperTrapezoidalMatrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::transpose() const {
  int m=size(0),n=size(1);
  UnitUpperTrapezoidalMatrix<float,float> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<float,float>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(m-j-1,addr(j+1,j),1,U->addr(j,j+1),n);
  }
  return U;
}

template<> UnitUpperTrapezoidalMatrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::conjugateTranspose() const {
  return transpose();
}
*/

template<> LowerTrapezoidalMatrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<float,float> *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(m,n);
  if (L_non_unit) {
    S->copy(L);
    for (int j=0;j<n;j++) {
      (*S)(j,j)+=float_one_;
      F77NAME(saxpy)(m-j-1,float_one_,addr(j+1,j),1,S->addr(j+1,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=2.;
      if (j<m-1) {
        F77NAME(scopy)(m-j-1,addr(j+1,j),1,S->addr(j+1,j),1);
        F77NAME(saxpy)(m-j-1,float_one_,addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> Matrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *S=OPERATOR_NEW Matrix<float,float>(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    (*S)(j,j)+=float_one_;
    F77NAME(saxpy)(m-j-1,float_one_,addr(j+1,j),1,S->addr(j+1,j),1);
  }
  return S;
}

template<> LowerTrapezoidalMatrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTrapezoidalMatrix<float,float> *dif=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(m,n);
  for (int j=0;j<n;j++) {
    (*dif)(j,j)=(L_non_unit ? float_one_-L(j,j) : float_zero_);
    if (j<m-1) {
      F77NAME(scopy)(m-j-1,addr(j+1,j),1,dif->addr(j+1,j),1);
      F77NAME(saxpy)(m-j-1,float_mone_,L.addr(j+1,j),1,
        dif->addr(j+1,j),1);
    }
  }
  return dif;
}

template<> Matrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *dif=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(slaset)('U',m,n,float_zero_,float_zero_,dif->addr(),m);
  for (int j=0;j<n;j++) {
    (*dif)(j,j)=float_one_-M(j,j);
    F77NAME(scopy)(m-j-1,addr(j+1,j),1,dif->addr(j+1,j),1);
    F77NAME(saxpy)(m-j-1,float_mone_,M.addr(j+1,j),1,dif->addr(j+1,j),1);
  }
  return dif;
}

template<> LowerTrapezoidalMatrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator*(float d) const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<float,float> *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(m,n);
  for (int j=0;j<n;j++) {
    (*P)(j,j)=d;
    if (j<m-1) {
      F77NAME(scopy)(m-j-1,addr(j+1,j),1,P->addr(j+1,j),1);
      F77NAME(sscal)(m-j-1,d,P->addr(j+1,j),1);
    }
  }
  return P;
}

template<> LowerTrapezoidalMatrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator/(float d) const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<float,float> *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(m,n);
  float dinv=float_one_/d;
  for (int j=0;j<n;j++) {
    (*P)(j,j)=dinv;
    if (j<m-1) {
      F77NAME(scopy)(m-j-1,addr(j+1,j),1,P->addr(j+1,j),1);
      F77NAME(sscal)(m-j-1,dinv,P->addr(j+1,j),1);
    }
  }
  return P;
}

template<> LowerTrapezoidalMatrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &L) const {
// to compute the jth column of the product
//   [ L_11   0  ] [ 0 ] = [    0   ]
//   [ L_21 L_22 ] [ m ] = [ L_22 m ]
//   [ L_31 L_32 ]       = [ L_32 m ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  LowerTrapezoidalMatrix<float,float> *P=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(m,n);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(k-j,L.addr(j,j),1,P->addr(j,j),1);
      F77NAME(strmv)('L','N','U',k-j,addr(j,j),m,P->addr(j,j),1);
      if (m>k) {
        F77NAME(sgemv)('N',m-k,k-j,float_one_,addr(k,j),m,L.addr(j,j),1,
          float_zero_,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*P)(j,j)=float_one_;
      if (j<k-1) {
        F77NAME(scopy)(k-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
        F77NAME(strmv)('L','N','U',k-j-1,addr(j+1,j+1),m,
          P->addr(j+1,j),1);
        F77NAME(saxpy)(k-j-1,float_one_,addr(j+1,j),1,P->addr(j+1,j),1);
      }
      if (m>k) {
        if (j<k-1) {
          F77NAME(sgemv)('N',m-k,k-j-1,float_one_,addr(k,j+1),m,
            L.addr(j+1,j),1,float_zero_,P->addr(k,j),1);
        }
        F77NAME(saxpy)(m-k,float_one_,addr(k,j),1,P->addr(k,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  F77NAME(slacpy)('A',k,n,M.addr(),k,P->addr(),m);
  F77NAME(strmm)('L','L','N','U',k,n,float_one_,addr(),m,
    P->addr(),m);
  if (m>k) {
    F77NAME(sgemm)('N','N',m-k,n,k,float_one_,addr(k,0),m,
      M.addr(),k,float_one_,P->addr(k,0),m);
  }
  return P;
}

template<> Vector<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<float,float> *p=
    OPERATOR_NEW Vector<float,float>(m,float_zero_);
  F77NAME(scopy)(n,v.addr(),1,p->addr(),1);
  F77NAME(strmv)('L','N','U',n,addr(),m,p->addr(),1);
  if (m>n) {
    F77NAME(sgemv)('N',m-n,n,float_one_,addr(n,0),m,
      v.addr(),1,float_one_,p->addr(n),m);
  }
  return p;
}

template<> Vector<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::trmv(
const Vector<float,float> &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<float,float> *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<float,float>(n);
    F77NAME(scopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(strmv)('L','T','U',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(sgemv)('T',m-n,n,one_,addr(n,0),m,x.addr(n),1,one_,
        p->addr(),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<float,float>(m);
    F77NAME(scopy)(n,x.addr(),1,p->addr(),1);
    F77NAME(strmv)('L','N','U',n,addr(),m,p->addr(),1);
    if (m>n) {
      F77NAME(sgemv)('N',m-n,n,one_,addr(n,0),m,x.addr(),1,zero_,
        p->addr(n),1);
    }
  }
  return p;
}

template<> Matrix<float,float>*
UnitLowerTrapezoidalMatrix<float,float>::trmm(
const Matrix<float,float> &X,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<float,float> *P=0;
  if (side=='L' || side=='l') {
    int k=X.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,X.size(0));
      P=OPERATOR_NEW Matrix<float,float>(n,k);
      F77NAME(slacpy)('A',n,k,X.addr(),m,P->addr(),n);
      F77NAME(strmm)('L','L','T','U',n,k,one_,addr(),m,P->addr(),n);
      if (m>n) {
        F77NAME(sgemm)('T','N',n,k,m-n,one_,addr(n,0),m,X.addr(n,0),m,
          one_,P->addr(),n);
      }
    } else {
      CHECK_SAME(n,X.size(0));
      P=OPERATOR_NEW Matrix<float,float>(m,k);
      F77NAME(slacpy)('A',n,k,X.addr(),n,P->addr(),m);
      F77NAME(strmm)('L','L','N','U',n,k,one_,addr(),m,P->addr(),m);
      if (m>n) {
        F77NAME(sgemm)('N','N',m-n,k,n,one_,addr(n,0),m,X.addr(),n,
          zero_,P->addr(n,0),m);
      }
    }
  } else {
    int k=X.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,X.size(1));
      P=OPERATOR_NEW Matrix<float,float>(k,m);
      F77NAME(slacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(strmm)('R','L','T','U',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(sgemm)('N','T',k,m-n,n,one_,X.addr(),k,addr(n,0),m,
          zero_,P->addr(0,n),k);
      }
    } else {
      CHECK_SAME(m,X.size(1));
      P=OPERATOR_NEW Matrix<float,float>(k,n);
      F77NAME(slacpy)('A',k,n,X.addr(),k,P->addr(),k);
      F77NAME(strmm)('R','L','N','U',k,n,one_,addr(),m,P->addr(),k);
      if (m>n) {
        F77NAME(sgemm)('N','N',k,n,m-n,one_,X.addr(0,n),k,addr(n,0),m,
          one_,P->addr(),k);
      }
    }
  }
  return P;
}

template<> void UnitLowerTrapezoidalMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(scopy)(n,b.addr(),1,x.addr(),1);
    if (m>n) { // use trailing entries of x as free variables
      F77NAME(sgemv)('T',m-n,n,mone_,addr(n,0),m,x.addr(n),1,one_,
        x.addr(),1);
    }
    F77NAME(strsv)('L','T','U',n,addr(),m,x.addr(),1);
  } else { // calling routine will have to check consistency conditions
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(scopy)(n,b.addr(),1,x.addr(),1);
    F77NAME(strsv)('L','N','U',n,addr(),m,x.addr(),1);
  }
}

template<> void UnitLowerTrapezoidalMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(slacpy)('A',n,k,B.addr(),n,X.addr(),n);
      if (m>n) { // use these entries of X as free variables
        F77NAME(sgemm)('T','N',n,k,m-n,mone_,addr(n,0),m,X.addr(n,0),m,
          one_,X.addr(),m);
      }
      F77NAME(strsm)('L','L','T','U',n,k,one_,addr(),m,X.addr(),m);
    } else { // calling routine will have to check consistency conditions
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(slacpy)('A',n,k,B.addr(),m,X.addr(),n);
      F77NAME(strsm)('L','L','N','U',n,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(slacpy)('A',k,n,B.addr(),k,X.addr(),k);
      F77NAME(strsm)('R','L','T','U',k,n,one_,addr(),m,X.addr(),k);
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(slacpy)('A',k,n,B.addr(),k,X.addr(),k);
      if (m>n) { // use these entries of X as free variables
        F77NAME(sgemm)('N','N',k,n,m-n,mone_,X.addr(0,n),n,addr(n,0),m,
          one_,X.addr(),k);
      }
      F77NAME(strsm)('R','L','N','U',k,n,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> float
UnitLowerTrapezoidalMatrix<float,float>::normFrobenius() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('F','L','U',m,size(1),addr(),m,work);
}

template<> float
UnitLowerTrapezoidalMatrix<float,float>::normInfinity() const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float result=F77NAME(slantr)('I','L','U',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> float
UnitLowerTrapezoidalMatrix<float,float>::normMaxEntry() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('M','L','U',m,size(1),addr(),m,work);
}

template<> float
UnitLowerTrapezoidalMatrix<float,float>::normOne() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('O','L','U',m,size(1),addr(),m,work);
}

template class UnitLowerTrapezoidalMatrix<float,float>;
template void testUnitLowerTrapezoidalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float
UnitLowerTriangularMatrix<float,float>::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  F77NAME(strcon)(norm,'L','U',n,addr(),n,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete [] work;
  return rcond;
}

/*
template<> UnitLowerTriangularMatrix<float,float>*
UnitLowerTriangularMatrix<float,float>::inverse() const {
  UnitLowerTriangularMatrix<float,float> *result=
    OPERATOR_NEW UnitLowerTriangularMatrix<float,float>(*this);
  int n=size(0);
  int info;
  F77NAME(strtri)('L','U',n,result->addr(),n,info);
  CHECK_TEST(info==0);
  return result;
}
*/

template class UnitLowerTriangularMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> Matrix<float,float>*
UpperTrapezoidalMatrix<float,float>::makeMatrix() const {
  int m=size(0),n=size(1);
  Matrix<float,float> *M=OPERATOR_NEW
    Matrix<float,float>(m,n,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(this)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+1,m),addr(0,j),1,M->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j,m),addr(0,j),1,M->addr(0,j),1);
      if (j<m) (*M)(j,j)=float_one_;
    }
  }
  return M;
}

template<> UpperTrapezoidalMatrix<float,float>& 
UpperTrapezoidalMatrix<float,float>::operator+=(
const UpperTrapezoidalMatrix<float,float> &M) {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0))
  CHECK_SAME(n,M.size(1))
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(j+1,m),one_,M.addr(0,j),1,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<float,float>&
UpperTrapezoidalMatrix<float,float>::operator-=(
const UpperTrapezoidalMatrix<float,float> &M) {
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0))
  CHECK_SAME(n,M.size(1))
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(j+1,m),mone_,M.addr(0,j),1,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<float,float>& 
UpperTrapezoidalMatrix<float,float>::operator*=(float d) {
  int m=size(0),n=size(1);
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(min(j+1,m),d,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<float,float>&
UpperTrapezoidalMatrix<float,float>::operator/=(float d) {
  int m=size(0),n=size(1);
  float dinv=float_one_/d;
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(min(j+1,m),dinv,addr(0,j),1);
  }
  return *this;
}

template<> UpperTrapezoidalMatrix<float,float>*
UpperTrapezoidalMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<float,float> *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(m,n);
  if (U_non_unit) {
    S->copy(U);
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(min(j+1,m),float_one_,addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      if (j>0) {
        F77NAME(saxpy)(min(j,m),float_one_,addr(0,j),1,S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
UpperTrapezoidalMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(m,zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(j+1,addr(0,j),1,S->addr(0,j),1);
      F77NAME(saxpy)(m-j,float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(j+1,addr(0,j),1,S->addr(0,j),1);
      if (j<m-1) {
        F77NAME(saxpy)(m-j-1,float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> Matrix<float,float>*
UpperTrapezoidalMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *S=OPERATOR_NEW Matrix<float,float>(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(m,j+1),float_one_,addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> UpperTrapezoidalMatrix<float,float>*
UpperTrapezoidalMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<float,float> *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(m,n);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      F77NAME(saxpy)(min(j+1,m),float_mone_,U.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+1,m),addr(0,j),1,S->addr(0,j),1);
      if (j>0) {
        F77NAME(saxpy)(min(j,m),float_mone_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)-=float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
UpperTrapezoidalMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *D=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(j+1,addr(0,j),1,D->addr(0,j),1);
      F77NAME(saxpy)(m-j,float_mone_,L.addr(j,j),1,D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(j+1,addr(0,j),1,D->addr(0,j),1);
      if (j<m-1) {
        F77NAME(saxpy)(m-j-1,float_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
      (*D)(j,j)-=float_one_;
    }
  }
  return D;
}

template<> Matrix<float,float>*
UpperTrapezoidalMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *D=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(j+1,addr(0,j),1,D->addr(0,j),1);
    F77NAME(saxpy)(m,float_mone_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> UpperTrapezoidalMatrix<float,float>* 
UpperTrapezoidalMatrix<float,float>::operator*(
const UpperTrapezoidalMatrix<float,float> &U) const {
// [ U_1 , U_2 ] [ V_11 , V_12 , V_13 ]
//               [   0  , V_22 , V_23 ]
//   = [ U_1 V_11 , U_1 V_12 + U_2 V_22 , U_1 V_13 + U_2 V_23 ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  int m=size(0),k=size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  UpperTrapezoidalMatrix<float,float> *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(m,n,float_zero_);
  if (U_non_unit) {
    for (int j=0;j<m;j++) {
      F77NAME(scopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(strmv)('U','N','N',j+1,addr(),m,P->addr(0,j),1);
    }
    for (int j=m;j<k;j++) {
      F77NAME(scopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(strmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
      F77NAME(sgemv)('N',m,j-m+1,float_one_,addr(0,m),m,U.addr(m,j),1,
        float_one_,P->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(scopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(strmv)('U','N','N',j,addr(),m,P->addr(0,j),1);
        F77NAME(saxpy)(j,float_one_,addr(0,j),1,P->addr(0,j),1);
      }
      (*P)(j,j)=(*this)(j,j);
    }
    for (int j=m;j<k;j++) {
      F77NAME(scopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(strmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
      if (j>m) {
        F77NAME(sgemv)('N',m,j-m,float_one_,addr(0,m),m,U.addr(m,j),1,
          float_one_,P->addr(0,j),1);
      }
      F77NAME(saxpy)(m,float_one_,addr(0,j),1,P->addr(0,j),1);
    }
  }
  for (int j=k;j<n;j++) {
    F77NAME(scopy)(m,U.addr(0,j),1,P->addr(0,j),1);
    F77NAME(strmv)('U','N','N',m,addr(),m,P->addr(0,j),1);
    F77NAME(sgemv)('N',m,k-m,float_one_,addr(0,m),m,U.addr(m,j),1,
      float_one_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<float,float>*
UpperTrapezoidalMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  if (m>=n) {
//  note that
//  [ U_11 U_12 U_13 ] [ L_1 ] = [ U_11 L_1 + U_12 L_2 + U_13 L_3 ]
//  [      U_22 U_23 ] [ L_2 ] = [            U_22 L_2 + U_23 L_3 ]
//                     [ L_3 ]
    if (L_non_unit) { // U_11 L_1:
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(sgemv)('N',j,n-j,float_one_,addr(0,j),m,L.addr(j,j),1,
            float_zero_,P->addr(0,j),1);
        }
        F77NAME(scopy)(n-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(strmv)('U','N','N',n-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(saxpy)(j,float_one_,addr(0,j),1,P->addr(0,j),1);
          if (j<n-1) {
            F77NAME(sgemv)('N',j,n-j-1,float_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,float_zero_,P->addr(0,j),1);
          }
        }
        if (j<n-1) {
          F77NAME(scopy)(n-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(strmv)('U','N','N',n-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(sdot)(n-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (n<k) { // U_12 L_2 + U_13 L_3:
      F77NAME(sgemm)('N','N',n,n,k-n,float_one_,addr(0,n),m,
        L.addr(n,0),k,float_one_,P->addr(0,0),m);
    }
    if (n<m) {
      F77NAME(slacpy)('A',m-n,n,L.addr(n,0),k,P->addr(n,0),m);
      F77NAME(strmm)('L','U','N','N',m-n,n,float_one_,addr(n,n),m,
        P->addr(n,0),m); // U_22 L_2
      if (m<k) {
        F77NAME(sgemm)('N','N',m-n,n,k-m,float_one_,addr(n,m),m,
          L.addr(m,0),k,float_one_,P->addr(n,0),m); // U_23 L_3
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
          F77NAME(sgemv)('N',j,m-j,float_one_,addr(0,j),m,L.addr(j,j),1,
            float_zero_,P->addr(0,j),1);
        }
        F77NAME(scopy)(m-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(strmv)('U','N','N',m-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(saxpy)(j,float_one_,addr(0,j),1,P->addr(0,j),1);
          if (j<m-1) {
            F77NAME(sgemv)('N',j,m-j-1,float_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,float_zero_,P->addr(0,j),1);
          }
        }
        if (j<m-1) {
          F77NAME(scopy)(m-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(strmv)('U','N','N',m-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(sdot)(m-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=(*this)(j,j);
      }
    }
    if (m<k) { // U_2 L_21 + U_3 L_31 
      F77NAME(sgemm)('N','N',m,m,k-m,float_one_,addr(0,m),m,
        L.addr(m,0),k,float_one_,P->addr(0,0),m);
    }
    char diagL=(L_non_unit ? 'N' : 'U');
    F77NAME(slacpy)('A',m,n-m,addr(0,m),m,P->addr(0,m),m);
    F77NAME(strmm)('R','L','N',diagL,m,n-m,float_one_,L.addr(m,m),k,
      P->addr(0,m),m); // U_2 L_22
    if (n<k) {
      F77NAME(sgemm)('N','N',m,n-m,k-n,float_one_,addr(0,n),m,
        L.addr(n,m),k,float_one_,P->addr(0,m),m); // U_3 L_32
    }
  }
  return P;
}

template<> Matrix<float,float>*
UpperTrapezoidalMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
// [ U_1 U_2 ] [ M_1 ] = U_1 M_1 + U_2 M_2
//             [ M_2 ]
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<float,float> *P=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(slacpy)('A',m,n,M.addr(),k,P->addr(),m);
  F77NAME(strmm)('L','U','N','N',m,n,float_one_,addr(),m,
    P->addr(),m); // U_1 M_1
  if (k>m) {
    F77NAME(sgemm)('N','N',m,n,k-m,float_one_,addr(0,m),m,
      M.addr(m,0),k,float_one_,P->addr(),m); // U_2 M_2
  }
  return P;
}

template<> Vector<float,float>*
UpperTrapezoidalMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<float,float> *p=OPERATOR_NEW Vector<float,float>(m);
  F77NAME(scopy)(m,v.addr(),1,p->addr(),1);
  F77NAME(strmv)('U','N','N',m,addr(),m,p->addr(),1);
  F77NAME(sgemv)('N',m,n-m,float_one_,addr(0,m),m,v.addr(m),1,
    float_one_,p->addr(),1);
  return p;
}

/*
template<> LowerTrapezoidalMatrix<float,float>*
UpperTrapezoidalMatrix<float,float>::transpose() const {
  int m=size(0),n=size(1);
  LowerTrapezoidalMatrix<float,float> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(j+1,m),addr(0,j),1,L->addr(j,0),n);
  }
  return L;
}

template<> LowerTrapezoidalMatrix<float,float>*
UpperTrapezoidalMatrix<float,float>::conjugateTranspose() const {
  return transpose();
}
*/

template<> Vector<float,float>*
UpperTrapezoidalMatrix<float,float>::trmv(
const Vector<float,float> &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<float,float> *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<float,float>(n);
    F77NAME(scopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(strmv)('U','T','N',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(sgemv)('T',m,n-m,one_,addr(0,m),m,x.addr(),1,zero_,
        p->addr(m),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<float,float>(m);
    F77NAME(scopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(strmv)('U','N','N',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(sgemv)('N',m,n-m,one_,addr(0,m),m,x.addr(m),1,one_,
        p->addr(),1);
    }
  }
  return p;
}

template<> Matrix<float,float>*
UpperTrapezoidalMatrix<float,float>::trmm(
const Matrix<float,float> &M,char side,char trans) const {
//TRACER_CALL(tr,"UpperTrapezoidalMatrix::trmm");
  int m=size(0),n=size(1);
  Matrix<float,float> *P=0;
  if (side=='L' || side=='l') {
    int k=M.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,M.size(0));
      P=OPERATOR_NEW Matrix<float,float>(n,k);
      F77NAME(slacpy)('A',m,k,M.addr(),m,P->addr(),n);
      F77NAME(strmm)('L','U','T','N',m,k,one_,addr(),m,P->addr(),n);
      if (n>m) {
        F77NAME(sgemm)('T','N',n-m,k,m,one_,addr(0,m),m,M.addr(),m,
          zero_,P->addr(m,0),n);
      }
    } else {
      CHECK_SAME(n,M.size(0));
      P=OPERATOR_NEW Matrix<float,float>(m,k);
      F77NAME(slacpy)('A',m,k,M.addr(),n,P->addr(),m);
      F77NAME(strmm)('L','U','N','N',m,k,one_,addr(),m,P->addr(),m);
      if (n>m) {
        F77NAME(sgemm)('N','N',m,k,n-m,one_,addr(0,m),m,M.addr(m,0),n,
          one_,P->addr(),m);
      }
    }
  } else {
    int k=M.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,M.size(1));
      P=OPERATOR_NEW Matrix<float,float>(k,m);
      F77NAME(slacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(strmm)('R','U','T','N',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(sgemm)('N','T',k,m,n-m,one_,M.addr(0,m),k,addr(0,m),m,
          one_,P->addr(),k);
      }
    } else {
      CHECK_SAME(m,M.size(1));
      P=OPERATOR_NEW Matrix<float,float>(k,n);
      F77NAME(slacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(strmm)('R','U','N','N',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(sgemm)('N','N',k,n-m,m,one_,M.addr(),k,addr(0,m),m,
          zero_,P->addr(0,m),k);
      }
    }
  }
  return P;
}

template<> float UpperTrapezoidalMatrix<float,float>::normFrobenius()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('F','U','N',m,size(1),addr(),m,work);
}

template<> float UpperTrapezoidalMatrix<float,float>::normInfinity()
const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float result=F77NAME(slantr)('I','U','N',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> float UpperTrapezoidalMatrix<float,float>::normMaxEntry()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('M','U','N',m,size(1),addr(),m,work);
}

template<> float UpperTrapezoidalMatrix<float,float>::normOne()
const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('O','U','N',m,size(1),addr(),m,work);
}

template<> void UpperTrapezoidalMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    // calling routine will have to check consistency conditions
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(scopy)(m,b.addr(),1,x.addr(),1);
    F77NAME(strsv)('U','T','N',m,addr(),m,x.addr(),1);
  } else {
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(scopy)(m,b.addr(),1,x.addr(),1);
    if (n>m) { // use trailing entries of x as free variables
      F77NAME(sgemv)('N',m,n-m,mone_,addr(0,m),m,x.addr(m),1,one_,
        x.addr(),1);
    }
    F77NAME(strsv)('U','N','N',m,addr(),m,x.addr(),1);
  }
}

template<> void UpperTrapezoidalMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(slacpy)('A',m,k,B.addr(),n,X.addr(),m);
      F77NAME(strsm)('L','U','T','N',m,k,one_,addr(),m,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(slacpy)('A',m,k,B.addr(),m,X.addr(),n);
      if (n>m) { // use these entries of X as free variables
        F77NAME(sgemm)('N','N',m,k,n-m,mone_,addr(0,m),m,X.addr(m,0),n,
          one_,X.addr(),n);
      }
      F77NAME(strsm)('L','U','N','N',m,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(slacpy)('A',k,m,B.addr(),k,X.addr(),k);
      if (n>m) { // use these entries of X as free variables
        F77NAME(sgemm)('N','T',k,m,n-m,mone_,X.addr(0,m),k,addr(0,m),m,
          one_,X.addr(),k);
      }
      F77NAME(strsm)('R','U','T','N',k,m,one_,addr(),m,X.addr(),k);
    } else { // calling routine will have to check consistency
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(slacpy)('A',k,m,B.addr(),k,X.addr(),k);
      F77NAME(strsm)('R','U','N','N',k,m,one_,addr(),m,X.addr(),k);
    }
  }
}

template<> SquareMatrix<float,float>* operator+(
const UnitLowerTrapezoidalMatrix<float,float> &L,
const UpperTrapezoidalMatrix<float,float> &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    (*S)(j,j)=float_one_;
    if (j<m-1) {
      F77NAME(scopy)(m-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
    }
    if (U_non_unit) {
      F77NAME(saxpy)(j+1,float_one_,U.addr(0,j),1,S->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(saxpy)(j,float_one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator+(
const LowerTrapezoidalMatrix<float,float> &L,
const UpperTrapezoidalMatrix<float,float> &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    if (L_non_unit) {
      F77NAME(scopy)(m-j,L.addr(j,j),1,S->addr(j,j),1);
    } else {
      (*S)(j,j)=float_one_;
      if (j<m-1) {
        F77NAME(scopy)(m-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
    }
    if (U_non_unit) {
      F77NAME(saxpy)(j+1,float_one_,U.addr(0,j),1,S->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(saxpy)(j,float_one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> Matrix<float,float>* operator+(
const Matrix<float,float> &M,
const UpperTrapezoidalMatrix<float,float> &U) {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(not_trapezoidal);
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  Matrix<float,float> *S=OPERATOR_NEW Matrix<float,float>(m,n);
  S->copy(M);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(min(j+1,m),float_one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j,m),float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const UnitLowerTrapezoidalMatrix<float,float> &L,
const UpperTrapezoidalMatrix<float,float> &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  SquareMatrix<float,float> *D=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    (*D)(j,j)=float_one_;
    if (j<m-1) {
      F77NAME(scopy)(m-j-1,L.addr(j+1,j),1,D->addr(j+1,j),1);
    }
    if (U_non_unit) {
      F77NAME(saxpy)(j+1,float_mone_,U.addr(0,j),1,D->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(saxpy)(j,float_mone_,U.addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)-=float_one_;
    }
  }
  return D;
}

template<> SquareMatrix<float,float>* operator-(
const LowerTrapezoidalMatrix<float,float> &L,
const UpperTrapezoidalMatrix<float,float> &U) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  SquareMatrix<float,float> *D=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    if (L_non_unit) {
      F77NAME(scopy)(m-j,L.addr(j,j),1,D->addr(j,j),1);
    } else {
      (*D)(j,j)=float_one_;
      if (j<m-1) {
        F77NAME(scopy)(m-j-1,L.addr(j+1,j),1,D->addr(j+1,j),1);
      }
    }
    if (U_non_unit) {
      F77NAME(saxpy)(j+1,float_mone_,U.addr(0,j),1,D->addr(0,j),1);
    } else {
      if (j>0) {
        F77NAME(saxpy)(j,float_mone_,U.addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)-=float_one_;
    }
  }
  return D;
}

template<> Matrix<float,float>* operator-(
const Matrix<float,float> &M,
const UpperTrapezoidalMatrix<float,float> &U) {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(not_trapezoidal);
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  Matrix<float,float> *D=OPERATOR_NEW Matrix<float,float>(m,n);
  D->copy(M);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(min(j+1,m),float_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j,m),float_mone_,U.addr(0,j),1,
          D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)-=float_one_;
    }
  }
  return D;
}

template<> Matrix<float,float>* operator*(
const UnitLowerTrapezoidalMatrix<float,float> &L,
const UpperTrapezoidalMatrix<float,float> &U) {
//  Note that
//  [ L_1 ] [ U_1 U_2 ] = [ L_1 U_1 , L_1 U_2 ]
//  [ L_2 ]             = [ L_2 U_1 , L_2 U_2 ]
//  and that
//  [ L_11      ] [ u ] = [ L_11 u ]
//  [ L_21 L_22 ] [   ] = [ L_21 u ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U) ==0);
  int m=L.size(0),k=L.size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<k;j++) { // L_1 U_1
    if (U_non_unit) {
      F77NAME(scopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(strmv)('L','N','U',j+1,L.addr(),m,P->addr(0,j),1);
      if (j+1<k) {
        F77NAME(sgemv)('N',k-j-1,j+1,float_one_,L.addr(j+1,0),m,
          U.addr(0,j),1,float_zero_,P->addr(j+1,j),1);
      }
    } else {
      if (j>0) {
        F77NAME(scopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(strmv)('L','N','U',j,L.addr(),m,P->addr(0,j),1);
        (*P)(j,j)=F77NAME(sdot)(j,L.addr(j,0),m,U.addr(0,j),1);
        if (j+1<k) {
          F77NAME(sgemv)('N',k-j-1,j,float_one_,L.addr(j+1,0),m,
            U.addr(0,j),1,float_zero_,P->addr(j+1,j),1);
        }
      }
      (*P)(j,j)+=float_one_;
      if (j+1<k) {
        F77NAME(saxpy)(k-j-1,float_one_,L.addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
    }
  }
  if (n>k) { // L_1 U_2
    F77NAME(slacpy)('A',k,n-k,U.addr(0,k),k,P->addr(0,k),m);
    F77NAME(strmm)('L','L','N','U',k,n-k,float_one_,L.addr(),m,
      P->addr(0,k),m);
  }
  if (m>k) {
    char diagU=(U_non_unit ? 'N' : 'U');
    F77NAME(slacpy)('A',m-k,k,L.addr(k,0),m,P->addr(k,0),m);
    F77NAME(strmm)('R','U','N',diagU,m-k,k,float_one_,U.addr(),k,
      P->addr(k,0),m); // L_2 U_1
    if (n>k) { // L_2 U_2
      F77NAME(sgemm)('N','N',m-k,n-k,k,float_one_,L.addr(k,0),m,
        U.addr(0,k),k,float_zero_,P->addr(k,k),m);
    }
  }
  return P;
}

template<> Matrix<float,float>* operator*(
const LowerTrapezoidalMatrix<float,float> &L,
const UpperTrapezoidalMatrix<float,float> &U) {
//  Note that
//  [ L_1 ] [ U_1 U_2 ] = [ L_1 U_1 , L_1 U_2 ]
//  [ L_2 ]             = [ L_2 U_1 , L_2 U_2 ]
//  and that
//  [ L_11      ] [ u ] = [ L_11 u ]
//  [ L_21 L_22 ] [   ] = [ L_21 u ]
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  char diagL=(L_non_unit ? 'N' : 'U');
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  int m=L.size(0),k=L.size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<k;j++) { // L_1 U_1
    if (U_non_unit) {
      F77NAME(scopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(strmv)('L','N',diagL,j+1,L.addr(),m,P->addr(0,j),1);
      if (j+1<k) {
        F77NAME(sgemv)('N',k-j-1,j+1,float_one_,L.addr(j+1,0),m,
          U.addr(0,j),1,float_zero_,P->addr(j+1,0),1);
      }
    } else {
      if (j>0) {
        F77NAME(scopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(strmv)('L','N',diagL,j,L.addr(),m,P->addr(0,j),1);
        (*P)(j,j)=F77NAME(sdot)(j,L.addr(j,0),m,U.addr(0,j),1);
        if (j+1<k) {
          F77NAME(sgemv)('N',k-j-1,j,float_one_,L.addr(j+1,0),m,
            U.addr(0,j),1,float_zero_,P->addr(j+1,j),1);
        }
      }
      (*P)(j,j)+=(L_non_unit ? L(j,j) : float_one_);
      if (j+1<k) {
        F77NAME(saxpy)(k-j-1,float_one_,L.addr(j+1,j),1,
          P->addr(j+1,j),1);
      }
    }
  }
  if (n>k) { // L_1 U_2
    F77NAME(slacpy)('A',k,n-k,U.addr(0,k),k,P->addr(0,k),m);
    F77NAME(strmm)('L','L','N',diagL,k,n-k,float_one_,L.addr(),m,
      P->addr(0,k),m);
  }
  if (m>k) {
    char diagU=(U_non_unit ? 'N' : 'U');
    F77NAME(slacpy)('A',m-k,k,L.addr(k,0),m,P->addr(k,0),m);
    F77NAME(strmm)('R','U','N',diagU,m-k,k,float_one_,U.addr(),k,
      P->addr(k,0),m); // L_2 U_1
    if (n>k) { // L_2 U_2
      F77NAME(sgemm)('N','N',m-k,n-k,k,float_one_,L.addr(k,0),m,
        U.addr(0,k),k,float_zero_,P->addr(k,k),m); 
    }
  }
  return P;
}

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,
const UpperTrapezoidalMatrix<float,float> &U) {
// Note that
// M [ U_1 , U_2 ] = [ M U_1 , M U_2 ]
  char diagU=
    (dynamic_cast<const UnitUpperTrapezoidalMatrix<float,float>*>(&U)
    ==0 ? 'N' : 'U');
  int m=M.size(0),k=U.size(0),n=U.size(1);
  CHECK_SAME(k,M.size(1));
  Matrix<float,float> *P=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(slacpy)('A',m,k,M.addr(),m,P->addr(),m);
  F77NAME(strmm)('R','U','N',diagU,m,k,float_one_,U.addr(),k,
    P->addr(),m); // M U_1
  if (n>k) { // M U_2
    F77NAME(sgemm)('N','N',m,n-k,k,float_one_,M.addr(),m,U.addr(0,k),k,
      float_zero_,P->addr(0,k),m);
  }
  return P;
}

template class UpperTrapezoidalMatrix<float,float>;
template void testUpperTrapezoidalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float
UpperTriangularMatrix<float,float>::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  F77NAME(strcon)(norm,'U','N',n,addr(),n,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete [] work;
  return rcond;
}

/*
template<> UpperTriangularMatrix<float,float>*
UpperTriangularMatrix<float,float>::inverse() const {
  UpperTriangularMatrix<float,float> *I=
    OPERATOR_NEW UpperTriangularMatrix<float,float>(*this);
  int n=size(0);
  int info;
  F77NAME(strtri)('U','N',n,I->addr(),n,info);
  CHECK_TEST(info==0);
  return I;
}
*/

template class UpperTriangularMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void
UnitUpperTrapezoidalMatrix<float,float>::copyFrom(int m,int n,
const Matrix<float,float> &U) {
  m=min(m,min(size(0),U.size(0)));
  n=min(n,min(size(1),U.size(1)));
  for (int j=1;j<n;j++) {
    F77NAME(scopy)(min(j,m),U.addr(0,j),1,addr(0,j),1);
  }
}

template<> UpperTrapezoidalMatrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<float,float> *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(m,n);
  if (U_non_unit) {
    S->copy(U);
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j,m),float_one_,addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)+=float_one_;
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(min(j,m),addr(0,j),1,S->addr(0,j),1);
        F77NAME(saxpy)(min(j,m),float_one_,U.addr(0,j),1,
          S->addr(0,j),1);
      }
      if (j<m) (*S)(j,j)=2.;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(m,float_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(j,addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)=float_one_;
      F77NAME(saxpy)(m-j,float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(j,addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)=2.;
      if (m<j-1) {
        F77NAME(saxpy)(m-j-1,float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> Matrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  bool M_not_trapezoidal=(dynamic_cast<const
    TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *S=OPERATOR_NEW Matrix<float,float>(m,n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(saxpy)(min(j,m),float_one_,addr(0,j),1,S->addr(0,j),1);
    }
    if (j<m) (*S)(j,j)+=float_one_;
  }
  return S;
}

template<> UpperTrapezoidalMatrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &U) const {
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTrapezoidalMatrix<float,float> *D=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(m,n);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)=float_one_;
      F77NAME(saxpy)(min(j+1,m),float_mone_,U.addr(0,j),1,
        D->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
        F77NAME(saxpy)(min(j,m),float_mone_,addr(0,j),1,D->addr(0,j),1);
      }
      if (j<m) (*D)(j,j)=float_zero_;
    }
  }
  return D;
}

template<> SquareMatrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  int m=size(0),n=size(1);
  CHECK_SAME(m,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *D=
    OPERATOR_NEW SquareMatrix<float,float>(m,float_zero_);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(j,addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)=float_one_;
      F77NAME(saxpy)(m-j,float_mone_,L.addr(j,j),1,D->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(j,addr(0,j),1,D->addr(0,j),1);
      }
      (*D)(j,j)=float_zero_;
      if (m<j-1) {
        F77NAME(saxpy)(m-j-1,float_mone_,L.addr(j+1,j),1,
          D->addr(j+1,j),1);
      }
    }
  }
  return D;
}

template<> Matrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  bool M_not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&M)==0);
  CHECK_TEST(M_not_trapezoidal);
  int m=size(0),n=size(1);
  CHECK_SAME(m,M.size(0));
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *D=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(scopy)(min(j,m),addr(0,j),1,D->addr(0,j),1);
    }
    if (j<m) (*D)(j,j)=float_one_;
    F77NAME(saxpy)(m,float_one_,M.addr(0,j),1,D->addr(0,j),1);
  }
  return D;
}

template<> UpperTrapezoidalMatrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator*(float d) const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<float,float> *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(m,n);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(scopy)(min(j,m),addr(0,j),1,P->addr(0,j),1);
      F77NAME(sscal)(min(j,m),d,P->addr(0,j),1);
    }
    if (j<m) (*P)(j,j)=d;
  }
  return P;
}

template<> UpperTrapezoidalMatrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator/(float d) const {
  int m=size(0),n=size(1);
  UpperTrapezoidalMatrix<float,float> *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(m,n);
  float dinv=float_one_/d;
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(scopy)(min(j,m),addr(0,j),1,P->addr(0,j),1);
      F77NAME(sscal)(min(j,m),dinv,P->addr(0,j),1);
    }
    if (j<m) (*P)(j,j)=dinv;
  }
  return P;
}

template<> UpperTrapezoidalMatrix<float,float>* 
UnitUpperTrapezoidalMatrix<float,float>::operator*(
const UpperTrapezoidalMatrix<float,float> &U) const {
// [ U_1 , U_2 ] [ V_11 , V_12 , V_13 ]
//               [   0  , V_22 , V_23 ]
//   = [ U_1 V_11 , U_1 V_12 + U_2 V_22 , U_1 V_13 + U_2 V_23 ]
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  int m=size(0),k=size(1),n=U.size(1);
  CHECK_SAME(k,U.size(0));
  UpperTrapezoidalMatrix<float,float> *P=
    OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(m,n,float_zero_);
  if (U_non_unit) {
    for (int j=0;j<m;j++) { // U_1 V_11
      F77NAME(scopy)(j+1,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(strmv)('U','N','U',j+1,addr(),m,P->addr(0,j),1);
    }
    for (int j=m;j<k;j++) {
//    the next two lines could have been slacpy & strmm, outside j loop
      F77NAME(scopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(strmv)('U','N','U',m,addr(),m,P->addr(0,j),1); // U_1 V_12
//    the next two lines could have been strmm('R',...)
//    but we could not replace both xxxmv with xxxmm
      F77NAME(sgemv)('N',m,j-m+1,float_one_,addr(0,m),m,U.addr(m,j),1,
        float_one_,P->addr(0,j),1); // U_2 V_22
    }
  } else {
    for (int j=0;j<m;j++) { // U_1 V_11
      if (j>0) {
        F77NAME(scopy)(j,U.addr(0,j),1,P->addr(0,j),1);
        F77NAME(strmv)('U','N','U',j,addr(),m,P->addr(0,j),1);
        F77NAME(saxpy)(j,float_one_,addr(0,j),1,P->addr(0,j),1);
      }
      (*P)(j,j)=float_one_;
    }
    for (int j=m;j<k;j++) {
      F77NAME(scopy)(m,U.addr(0,j),1,P->addr(0,j),1);
      F77NAME(strmv)('U','N','U',m,addr(),m,P->addr(0,j),1); // U_1 V_12
      if (j>m) { // U_2 V_22
        F77NAME(sgemv)('N',m,j-m,float_one_,addr(0,m),m,U.addr(m,j),1,
          float_one_,P->addr(0,j),1);
      }
      F77NAME(saxpy)(m,float_one_,addr(0,j),1,P->addr(0,j),1);
    }
  }
  for (int j=k;j<n;j++) {
    F77NAME(scopy)(m,U.addr(0,j),1,P->addr(0,j),1);
    F77NAME(strmv)('U','N','U',m,addr(),m,P->addr(0,j),1); // U_1 V_13
    F77NAME(sgemv)('N',m,k-m,float_one_,addr(0,m),m,U.addr(m,j),1,
      float_one_,P->addr(0,j),1); // U_2 V_23
  }
  return P;
}

template<> Matrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int m=size(0),k=size(1),n=L.size(1);
  CHECK_SAME(k,L.size(0));
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  if (m>=n) {
//  note that
//  [ U_11 U_12 U_13 ] [ L_1 ] = [ U_11 L_1 + U_12 L_2 + U_13 L_3 ]
//  [      U_22 U_23 ] [ L_2 ] = [            U_22 L_2 + U_23 L_3 ]
//                     [ L_3 ]
    if (L_non_unit) { // U_11 L_1:
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(sgemv)('N',j,n-j,float_one_,addr(0,j),m,L.addr(j,j),1,
            float_zero_,P->addr(0,j),1);
        }
        F77NAME(scopy)(n-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(strmv)('U','N','U',n-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        if (j>0) {
          F77NAME(saxpy)(j,float_one_,addr(0,j),1,P->addr(0,j),1);
          if (j<n-1) {
            F77NAME(sgemv)('N',j,n-j-1,float_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,float_zero_,P->addr(0,j),1);
          }
        }
        if (j<n-1) {
          F77NAME(scopy)(n-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(strmv)('U','N','U',n-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(sdot)(n-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=float_one_;
      }
    }
    if (n<k) { // U_12 L_2 + U_13 L_3:
      F77NAME(sgemm)('N','N',n,n,k-n,float_one_,addr(0,n),m,
        L.addr(n,0),k,float_one_,P->addr(0,0),m);
    }
    if (n<m) {
      F77NAME(slacpy)('A',m-n,n,L.addr(n,0),k,P->addr(n,0),m);
      F77NAME(strmm)('L','U','N','U',m-n,n,float_one_,addr(n,n),m,
        P->addr(n,0),m); // U_22 L_2
      if (m<k) {
        F77NAME(sgemm)('N','N',m-n,n,k-m,float_one_,addr(n,m),m,
          L.addr(m,0),k,float_one_,P->addr(n,0),m); // U_23 L_3
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
          F77NAME(sgemv)('N',j,m-j,float_one_,addr(0,j),m,L.addr(j,j),1,
            float_zero_,P->addr(0,j),1);
        }
        F77NAME(scopy)(m-j,L.addr(j,j),1,P->addr(j,j),1);
        F77NAME(strmv)('U','N','U',m-j,addr(j,j),m,P->addr(j,j),1);
      }
    } else {
      for (int j=0;j<m;j++) {
        if (j>0) {
          F77NAME(saxpy)(j,float_one_,addr(0,j),1,P->addr(0,j),1);
          if (j<m-1) {
            F77NAME(sgemv)('N',j,m-j-1,float_one_,addr(0,j+1),m,
              L.addr(j+1,j),1,float_zero_,P->addr(0,j),1);
          }
        }
        if (j<m-1) {
          F77NAME(scopy)(m-j-1,L.addr(j+1,j),1,P->addr(j+1,j),1);
          F77NAME(strmv)('U','N','U',m-j-1,addr(j+1,j+1),m,
            P->addr(j+1,j),1);
          (*P)(j,j)=F77NAME(sdot)(m-j-1,addr(j,j+1),m,L.addr(j+1,j),1);
        }
        (*P)(j,j)+=float_one_;
      }
    }
    if (m<k) { // U_2 L_21 + U_3 L_31 
      F77NAME(sgemm)('N','N',m,m,k-m,float_one_,addr(0,m),m,
        L.addr(m,0),k,float_one_,P->addr(0,0),m);
    }
    char diagL=(L_non_unit ? 'N' : 'U');
    F77NAME(slacpy)('A',m,n-m,addr(0,m),m,P->addr(0,m),m);
    F77NAME(strmm)('R','L','U',diagL,m,n-m,float_one_,L.addr(m,m),k,
      P->addr(0,m),m); // U_2 L_22
    if (n<k) {
      F77NAME(sgemm)('N','N',m,n-m,k-n,float_one_,addr(0,n),m,
        L.addr(n,m),k,float_one_,P->addr(0,m),m); // U_3 L_32
    }
  }
  return P;
}

template<> Matrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
// [ U_1 U_2 ] [ M_1 ] = U_1 M_1 + U_2 M_2
//             [ M_2 ]
  int m=size(0),k=size(1),n=M.size(1);
  CHECK_SAME(k,M.size(0));
  Matrix<float,float> *P=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(slacpy)('A',m,n,M.addr(),k,P->addr(),m);
  F77NAME(strmm)('L','U','N','U',m,n,float_one_,addr(),m,
    P->addr(),m);
  if (k>m) {
    F77NAME(sgemm)('N','N',m,n,k-m,float_one_,addr(0,m),m,
      M.addr(m,0),k,float_one_,P->addr(),m);
  }
  return P;
}

template<> Vector<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int m=size(0),n=size(1);
  CHECK_SAME(n,v.size());
  Vector<float,float> *p=OPERATOR_NEW Vector<float,float>(m);
  F77NAME(scopy)(m,v.addr(),1,p->addr(),1);
  F77NAME(strmv)('U','N','U',m,addr(),m,p->addr(),1);
  F77NAME(sgemv)('N',m,n-m,float_one_,addr(0,m),m,v.addr(m),1,
    float_one_,p->addr(),1);
  return p;
}

/*
template<> UnitLowerTrapezoidalMatrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::transpose() const {
  int m=size(0),n=size(1);
  UnitLowerTrapezoidalMatrix<float,float> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<float,float>(n,m);
  for (int i=0;i<m;i++) {
    F77NAME(scopy)(n-i-1,addr(i,i+1),m,L->addr(i+1,i),1);
  }
  return L;
}

template<> UnitLowerTrapezoidalMatrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::conjugateTranspose() const {
  return transpose();
}
*/

template<> Vector<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::trmv(
const Vector<float,float> &x,char trans) const {
  int m=size(0),n=size(1);
  Vector<float,float> *p=0;
  if (trans!='N' && trans!='n') {
    CHECK_SAME(m,x.size());
    p=OPERATOR_NEW Vector<float,float>(n);
    F77NAME(scopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(strmv)('U','T','U',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(sgemv)('T',m,n-m,one_,addr(0,m),m,x.addr(),1,zero_,
        p->addr(m),1);
    }
  } else {
    CHECK_SAME(n,x.size());
    p=OPERATOR_NEW Vector<float,float>(m);
    F77NAME(scopy)(m,x.addr(),1,p->addr(),1);
    F77NAME(strmv)('U','N','U',m,addr(),m,p->addr(),1);
    if (n>m) {
      F77NAME(sgemv)('N',m,n-m,one_,addr(0,m),m,x.addr(m),1,one_,
        p->addr(),1);
    }
  }
  return p;
}

template<> Matrix<float,float>*
UnitUpperTrapezoidalMatrix<float,float>::trmm(
const Matrix<float,float> &M,char side,char trans) const {
  int m=size(0),n=size(1);
  Matrix<float,float> *P=0;
  if (side=='L' || side=='l') {
    int k=M.size(1);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,M.size(0));
      P=OPERATOR_NEW Matrix<float,float>(n,k);
      F77NAME(slacpy)('A',m,k,M.addr(),m,P->addr(),n);
      F77NAME(strmm)('L','U','T','U',m,k,one_,addr(),m,P->addr(),n);
      if (n>m) {
        F77NAME(sgemm)('T','N',n-m,k,m,one_,addr(0,m),m,M.addr(),m,
          zero_,P->addr(m,0),n);
      }
    } else {
      CHECK_SAME(n,M.size(0));
      P=OPERATOR_NEW Matrix<float,float>(m,k);
      F77NAME(slacpy)('A',m,k,M.addr(),n,P->addr(),m);
      F77NAME(strmm)('L','U','N','U',m,k,one_,addr(),m,P->addr(),m);
      if (n>m) {
        F77NAME(sgemm)('N','N',m,k,n-m,one_,addr(0,m),m,M.addr(m,0),n,
          one_,P->addr(),m);
      }
    }
  } else {
    int k=M.size(0);
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,M.size(1));
      P=OPERATOR_NEW Matrix<float,float>(k,m);
      F77NAME(slacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(strmm)('R','U','T','U',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(sgemm)('N','T',k,m,n-m,one_,M.addr(0,m),k,addr(0,m),m,
          one_,P->addr(),k);
      }
    } else {
      CHECK_SAME(m,M.size(1));
      P=OPERATOR_NEW Matrix<float,float>(k,n);
      F77NAME(slacpy)('A',k,m,M.addr(),k,P->addr(),k);
      F77NAME(strmm)('R','U','N','U',k,m,one_,addr(),m,P->addr(),k);
      if (n>m) {
        F77NAME(sgemm)('N','N',k,n-m,m,one_,M.addr(),k,addr(0,m),m,
          zero_,P->addr(0,m),k);
      }
    }
  }
  return P;
}

template<> float
UnitUpperTrapezoidalMatrix<float,float>::normFrobenius() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('F','U','U',m,size(1),addr(),m,work);
}

template<> float
UnitUpperTrapezoidalMatrix<float,float>::normInfinity() const {
  int m=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,m);
  float result=F77NAME(slantr)('I','U','U',m,size(1),addr(),m,work);
  delete [] work;
  return result;
}

template<> float
UnitUpperTrapezoidalMatrix<float,float>::normMaxEntry() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('M','U','U',m,size(1),addr(),m,work);
}

template<> float
UnitUpperTrapezoidalMatrix<float,float>::normOne() const {
  int m=size(0);
  float *work=0;
  return F77NAME(slantr)('O','U','U',m,size(1),addr(),m,work);
}

template<> void UnitUpperTrapezoidalMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  int m=size(0),n=size(1);
  if (trans!='N' && trans!='n') {
    // calling routine will have to check consistency conditions
    CHECK_SAME(m,x.size());
    CHECK_SAME(n,b.size());
    F77NAME(scopy)(m,b.addr(),1,x.addr(),1);
    F77NAME(strsv)('U','T','U',m,addr(),m,x.addr(),1);
  } else {
    CHECK_SAME(n,x.size());
    CHECK_SAME(m,b.size());
    F77NAME(scopy)(m,b.addr(),1,x.addr(),1);
    if (n>m) { // use trailing entries of x as free variables
      F77NAME(sgemv)('N',m,n-m,mone_,addr(0,m),m,x.addr(m),1,one_,
        x.addr(),1);
    }
    F77NAME(strsv)('U','N','U',m,addr(),m,x.addr(),1);
  }
}

template<> void UnitUpperTrapezoidalMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  bool not_trapezoidal=
    (dynamic_cast<const TrapezoidalMatrix<float,float>*>(&B)==0);
  CHECK_TEST(not_trapezoidal);
  int m=size(0),n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      // calling routine will have to check consistency
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(slacpy)('A',m,k,B.addr(),n,X.addr(),m);
      F77NAME(strsm)('L','U','T','U',m,k,one_,addr(),m,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      F77NAME(slacpy)('A',m,k,B.addr(),m,X.addr(),n);
      if (n>m) { // use these entries of X as free variables
        F77NAME(sgemm)('N','N',m,k,n-m,mone_,addr(0,m),m,X.addr(m,0),n,
          one_,X.addr(),n);
      }
      F77NAME(strsm)('L','U','N','U',m,k,one_,addr(),m,X.addr(),n);
    }
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(m,B.size(1));
      CHECK_SAME(n,X.size(1));
      F77NAME(slacpy)('A',k,m,B.addr(),k,X.addr(),k);
      if (n>m) { // use these entries of X as free variables
        F77NAME(sgemm)('N','T',k,m,n-m,mone_,X.addr(0,m),k,addr(0,m),m,
          one_,X.addr(),k);
      }
      F77NAME(strsm)('R','U','T','U',k,m,one_,addr(),m,X.addr(),k);
    } else { // calling routine will have to check consistency
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(slacpy)('A',k,m,B.addr(),k,X.addr(),k);
      F77NAME(strsm)('R','U','N','U',k,m,one_,addr(),m,X.addr(),k);
    }
  }
}

template class UnitUpperTrapezoidalMatrix<float,float>;
template void testUnitUpperTrapezoidalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float
UnitUpperTriangularMatrix<float,float>::reciprocalConditionNumber(
char norm) const {
  int n=size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  F77NAME(strcon)(norm,'U','U',n,addr(),n,rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] iwork;
  delete [] work;
  return rcond;
}

/*
template<> UnitUpperTriangularMatrix<float,float>*
UnitUpperTriangularMatrix<float,float>::inverse() const {
  UnitUpperTriangularMatrix<float,float> *I=
    OPERATOR_NEW UnitUpperTriangularMatrix<float,float>(*this);
  int n=size(0);
  int info;
  F77NAME(strtri)('U','U',n,I->addr(),n,info);
  CHECK_TEST(info==0);
  return I;
}
*/

template class UnitUpperTriangularMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "OrthogonalMatrix.C"

/*
template<> OrthogonalMatrix<float,float>*
OrthogonalMatrix<float,float>::transpose() const {
  int m=size(0),n=size(1);
  OrthogonalMatrix<float,float> *Q=
    OPERATOR_NEW OrthogonalMatrix<float,float>(n,m);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(m,addr(0,j),1,Q->addr(j,0),n);
  }
  return Q;
}

template<> OrthogonalMatrix<float,float>*
OrthogonalMatrix<float,float>::conjugateTranspose() const {
  return transpose();
}
*/

/*
template<> void OrthogonalMatrix<float,float>::solve(
const UpperTrapezoidalMatrix<float,float> &B,
Matrix<float,float> &X,char side,char trans) const {
  int m=size(0), n=size(1);
  bool B_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&B)==0);
  X=float_zero_;
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<=min(j,n-1);i++) {
            F77NAME(saxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<min(j,n);i++) {
            F77NAME(saxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
          if (j<n) {
            F77NAME(saxpy)(m,float_one_,addr(0,j),1,X.addr(0,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(sdot)(min(j+1,m),addr(0,i),1,B.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(sdot)(min(j,m),addr(0,i),1,B.addr(0,j),1);
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
            X(i,j)=F77NAME(sdot)(m-i,B.addr(i,i),k,addr(i,j),1);
          }
        }
      } else {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=(*this)(i,j)
              +F77NAME(sdot)(m-i-1,B.addr(i,i+1),k,addr(i+1,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      if (B_non_unit) {
        for (int i=0;i<k;i++) {
          for (int j=i;j<n;j++) {
            F77NAME(saxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      } else {
        for (int i=0;i<k;i++) {
          F77NAME(saxpy)(m,float_one_,addr(0,i),1,X.addr(i,0),k);
          for (int j=i+1;j<n;j++) {
            F77NAME(saxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      }
    }
  }
}
*/

/*
template<> void OrthogonalMatrix<float,float>::solve(
const LowerTrapezoidalMatrix<float,float> &B,
Matrix<float,float> &X,char side,char trans) const {
  int m=size(0), n=size(1);
  bool B_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&B)==0);
  X=float_zero_;
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=j;i<n;i++) {
            F77NAME(saxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          F77NAME(saxpy)(m,float_one_,addr(0,j),1,X.addr(0,j),1);
          for (int i=j+1;i<n;i++) {
            F77NAME(saxpy)(m,B(i,j),addr(0,i),1,X.addr(0,j),1);
          }
        }
      }
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
      if (B_non_unit) {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=F77NAME(sdot)(m-j,addr(j,i),1,B.addr(j,j),1);
          }
        }
      } else {
        for (int j=0;j<k;j++) {
          for (int i=0;i<n;i++) {
            X(i,j)=(*this)(j,i)
              +F77NAME(sdot)(m-j-1,addr(j+1,i),1,B.addr(j+1,j),1);
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
            X(i,j)=F77NAME(sdot)(min(i+1,m),B.addr(i,0),k,addr(0,j),1);
          }
        }
      } else {
        for (int j=0;j<n;j++) {
          for (int i=0;i<k;i++) {
            X(i,j)=F77NAME(sdot)(min(i,m),B.addr(i,0),k,addr(0,j),1);
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
            F77NAME(saxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
        }
      } else {
        for (int i=0;i<k;i++) {
          for (int j=0;j<min(i,n);j++) {
            F77NAME(saxpy)(m,B(i,j),addr(0,j),1,X.addr(i,0),k);
          }
          if (i<n) {
            F77NAME(saxpy)(m,float_one_,addr(0,i),1,X.addr(i,0),k);
          }
        }
      }
    }
  }
}
*/

template<> void OrthogonalMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  int m=size(0), n=size(1);
  if (side=='L' || side=='l') {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1));
    if (trans!='N' && trans!='n') {
      CHECK_SAME(n,B.size(0));
      CHECK_SAME(m,X.size(0));
      F77NAME(sgemm)('N','N',m,k,n,float_one_,addr(),m,B.addr(),n,
        float_zero_,X.addr(),m);
    } else {
      CHECK_SAME(m,B.size(0));
      CHECK_SAME(n,X.size(0));
//    F77NAME(sgemm)('T','N',n,k,m,float_one_,addr(),m,B.addr(),m,
//      float_zero_,X.addr(),n);
      Vector<float,float> *r=OPERATOR_NEW Vector<float,float>(m);
      for (int l=0;l<k;l++) {
        F77NAME(scopy)(m,B.addr(0,l),1,r->addr(),1);
        for (int j=0;j<n;j++) {
          float &Xjl=X(j,l);
          Xjl=F77NAME(sdot)(m,r->addr(),1,addr(0,j),1); 
          F77NAME(saxpy)(m,-Xjl,addr(0,j),1,r->addr(),1);
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
//    F77NAME(sgemm)('N','N',k,n,m,float_one_,B.addr(),k,addr(),m,
//      float_zero_,X.addr(),k);
      Vector<float,float> *r=OPERATOR_NEW Vector<float,float>(m);
      for (int l=0;l<k;l++) {
        F77NAME(scopy)(m,B.addr(l,0),k,r->addr(),1);
        for (int j=0;j<n;j++) {
          float &Xlj=X(l,j);
          Xlj=F77NAME(sdot)(m,r->addr(),1,addr(0,j),1);
          F77NAME(saxpy)(m,-Xlj,addr(0,j),1,r->addr(),1);
        }
      }
      delete r; r=0;
    } else {
      CHECK_SAME(n,B.size(1));
      CHECK_SAME(m,X.size(1));
      F77NAME(sgemm)('N','T',k,m,n,float_one_,B.addr(),k,addr(),m,
        float_zero_,X.addr(),k);
    }
  }
}

template<> void OrthogonalMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  int m=size(0), n=size(1);
  if (trans!='N' && trans!='n') {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(sgemv)('N',m,n,float_one_,addr(),m,b.addr(),1,float_zero_,
      x.addr(),1);
  } else {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
//  F77NAME(sgemv)('T',m,n,float_one_,addr(),m,b.addr(),1,float_zero_,
//    x.addr(),1);
    Vector<float,float> *r=OPERATOR_NEW Vector<float,float>(m);
    r->copy(b);
    for (int j=0;j<n;j++) {
      x[j]=F77NAME(sdot)(m,r->addr(),1,addr(0,j),1); 
      F77NAME(saxpy)(m,-x[j],addr(0,j),1,r->addr(),1);
    }
    delete r; r=0;
  }
}

template class OrthogonalMatrix<float,float>;
template void testOrthogonalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SymmetricMatrix.C"

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,float> *M=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,addr(j,j),1,M->addr(j,j),1);
    if (j+1<n) {
      F77NAME(scopy)(n-j-1,addr(j+1,j),1,M->addr(j,j+1),n);
    }
  }
  return M;
}

template<> void SymmetricMatrix<float,float>::fillWith(float d) {
  Matrix<float,float>::set('L',d,d);
}

template<> float
SymmetricMatrix<float,float>::operator()(int i,int j) const {
  return (j<=i ? *(this->addr(i,j)) : *(this->addr(j,i)) );
}

template<> SymmetricMatrix<float,float>&
SymmetricMatrix<float,float>::operator+=(
const SymmetricMatrix<float,float> &S) {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0))
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_one_,S.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<float,float>&
SymmetricMatrix<float,float>::operator-=(
const SymmetricMatrix<float,float> &S) {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0))
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,S.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<float,float>&
SymmetricMatrix<float,float>::operator*=(float scalar) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(n-j,scalar,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricMatrix<float,float>&
SymmetricMatrix<float,float>::operator/=(float scalar) {
  int n=this->size(0);
  CHECK_NONZERO(scalar)
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(n-j,float_one_/scalar,addr(j,j),1);
  }
  return *this;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &U) const {
  int n=this->size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(scopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(min(j+1,n),one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j,n),one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      if (j<n) (*S)(j,j)+=one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int n=this->size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(scopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(n-j,one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(saxpy)(n-j-1,one_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      (*S)(j,j)+=one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator+(
const SquareMatrix<float,float> &S) const {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  T->copy(S);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,one_,addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n-j-1,one_,addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  int n=this->size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,one_,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n-j-1,one_,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &U) const {
  int n=this->size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(scopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(min(j+1,n),mone_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j,n),mone_,U.addr(0,j),1,S->addr(0,j),1);
      }
      if (j<n) (*S)(j,j)+=mone_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const UnitUpperTrapezoidalMatrix<float,float> &U,
const SymmetricMatrix<float,float> &S) {
  int n=S.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(scopy)(min(j,n),U.addr(0,j),1,T->addr(0,j),1);
    }
    if (j<n) (*T)(j,j)=float_one_;
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n-j-1,float_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<float,float>* operator-(
const UpperTrapezoidalMatrix<float,float> &U,
const SymmetricMatrix<float,float> &S) {
  int n=S.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+1,n),U.addr(0,j),1,T->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(min(j,n),U.addr(0,j),1,T->addr(0,j),1);
      }
      if (j<n) (*T)(j,j)=float_one_;
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n-j-1,float_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int n=this->size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(scopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(n-j,mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(saxpy)(n-j-1,mone_,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      (*S)(j,j)+=mone_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const UnitLowerTrapezoidalMatrix<float,float> &L,
const SymmetricMatrix<float,float> &S) {
  int n=S.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    if (j<n-1) {
      F77NAME(scopy)(n-j-1,L.addr(j+1,j),1,T->addr(j+1,j),1);
    }
    (*T)(j,j)=float_one_;
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n-j-1,float_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<float,float>* operator-(
const LowerTrapezoidalMatrix<float,float> &L,
const SymmetricMatrix<float,float> &S) {
  int n=S.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(n-j,L.addr(j,j),1,T->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j<n-1) {
        F77NAME(scopy)(n-j-1,L.addr(j+1,j),1,T->addr(j+1,j),1);
      }
      (*T)(j,j)=float_one_;
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n-j-1,float_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator-(
const SquareMatrix<float,float> &S) const {
  int n=this->size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(scopy)(n-j-1,addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n,mone_,S.addr(0,j),1,T->addr(0,j),1);
  }
  return T;
}

template<> SquareMatrix<float,float>* operator-(
const SquareMatrix<float,float> &S,
const SymmetricMatrix<float,float> &SS) {
  int n=SS.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  T->copy(S);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,SS.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n-j-1,float_mone_,SS.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  int n=this->size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(scopy)(n-j-1,addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n,mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const Matrix<float,float> &M,
const SymmetricMatrix<float,float> &S) {
  int n=S.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  T->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,S.addr(j,j),1,T->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n-j-1,float_mone_,S.addr(j+1,j),1,T->addr(j,j+1),n);
    }
  }
  return T;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator*(float d) const {
  int n=size(0);
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) { 
    const float *col_j=addr(j,j);
    float *S_col_j=S->addr(j,j);
    *S_col_j=(*col_j)*d;
    if (j<n-1) {
      S_col_j++;
      col_j++;
      float *S_row_j=S->addr(j,j+1);
      for (int i=j+1;i<n;i++,col_j++,S_col_j++,S_row_j+=n) {
        *S_row_j=*S_col_j=(*col_j)*d;
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator/(float d) const {
  CHECK_NONZERO(d);
  return operator*(float_one_/d);
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator*(
const SymmetricMatrix<float,float> &S) const {
// compute by bordering: note that
// [ sigma s^T ] [ tau t^T ] = [ sigma tau + s^T t , sigma t^T + s^T T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^T +  S  T ]
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int k=n-1;k>=0;k--) {
    if (k<n-1) {
      F77NAME(ssymv)('L',n-k-1,float_one_,addr(k+1,k+1),n,
        S.addr(k+1,k),1,float_zero_,T->addr(k+1,k),1); // S t
      F77NAME(ssymv)('L',n-k-1,float_one_,S.addr(k+1,k+1),n,
        addr(k+1,k),n,float_zero_,T->addr(k,k+1),n); // s^T T
      (*T)(k,k)=F77NAME(sdot)(n-k-1,addr(k+1,k),1,S.addr(k+1,k),1);//s^T t
    }
    F77NAME(sger)(n-k,n-k,float_one_,addr(k,k),1,S.addr(k,k),1,
      T->addr(k,k),n);
  }
  return T;
}

template<> Matrix<float,float>*
SymmetricMatrix<float,float>::operator*(
const UpperTrapezoidalMatrix<float,float> &U) const {
// compute by bordering: note that
// S [ U_1 , U_2 ] = [ S U_1 , S U_2 ]
// and that
// [ sigma s^T ] [ upsilon u^T ] = [ sigma upsilon , sigma u^T + s^T U ]
// [   s    S  ] [          U  ] [ [   s   upsilon ,   s   u^T +  S  U ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  char diag=(
    dynamic_cast<const UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0
    ? 'N' : 'U');
  for (int k=m-1;k>=0;k--) { // S U_1
    if (k<m-1) { // s^T U
      F77NAME(scopy)(m-k-1,addr(k+1,k),1,M->addr(k,k+1),m);
      F77NAME(strmv)('U','T',diag,m-k-1,U.addr(k+1,k+1),m,
        M->addr(k,k+1),m);
    }
    if (diag=='N') {
      F77NAME(sger)(m-k,m-k,float_one_,addr(k,k),1,U.addr(k,k),m,
        M->addr(k,k),m);
    } else {
      F77NAME(saxpy)(m-k,float_one_,addr(k,k),1,M->addr(k,k),1);
      F77NAME(sger)(m-k,m-k-1,float_one_,addr(k,k),1,U.addr(k,k+1),m,
        M->addr(k,k+1),m);
    }
  }
  if (n>m) { // S U_2
    F77NAME(ssymm)('L','L',m,n-m,float_one_,addr(),m,U.addr(0,m),m,
      float_zero_,M->addr(0,m),m);
  }
  return M;
// compute by columns: note that
// S [ U_1 , U_2 ] = [ S U_1 , S U_2 ]
// and that
// [ S_11 S_21^T ] [ u ] = [ S_11 u ]
// [ S_21  S_22  ] [   ] = [ S_21 u ]
//bool U_non_unit=(dynamic_cast<const
//  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
//if (U_non_unit) { // S U_1
//  for (int j=0;j<m;j++) {
//    F77NAME(ssymv)('L',j+1,float_one_,addr(),m,U.addr(0,j),1,
//      float_zero_,M->addr(0,j),1); // S_11 u
//    if (j+1<m) {
//      F77NAME(sgemv)('N',m-j-1,j+1,float_one_,addr(j+1,0),m,
//        U.addr(0,j),1,float_zero_,M->addr(j+1,j),1); // S_21 u
//    }
//  }
//} else {
//  for (int j=0;j<m;j++) {
//    if (j>0) {
//      F77NAME(ssymv)('L',j,float_one_,addr(),m,U.addr(0,j),1,
//        float_zero_,M->addr(0,j),1);
//      F77NAME(saxpy)(j,float_one_,addr(j,0),m,M->addr(0,j),1);
//      (*M)(j,j)=F77NAME(sdot)(j,addr(j,0),m,U.addr(0,j),1);
//    }
//    (*M)(j,j)+=(*this)(j,j);
//    if (j+1<m) {
//      if (j>0) {
//        F77NAME(sgemv)('N',m-j-1,j,float_one_,addr(j+1,0),m,
//          U.addr(0,j),1,float_zero_,M->addr(j+1,j),1);
//      }
//      F77NAME(saxpy)(m-j-1,float_one_,addr(j+1,j),1,M->addr(j+1,j),1);
//    }
//  }
//}
//if (n>m) { // S U_2
//  F77NAME(ssymm)('L','L',m,n-m,float_one_,addr(),m,U.addr(0,m),m,
//    float_zero_,M->addr(0,m),m);
//}
}

template<> Matrix<float,float>* operator*(
const UpperTrapezoidalMatrix<float,float> &U,
const SymmetricMatrix<float,float> &S) {
//compute by bordering: note that
// [ U_1 U_2 ] [ S_11 S_21^T ]
//             [ S_21  S_22  ]
//   = [ U_1 S_11 + U_2 S_21 , U_1 S_21^T + U_2 S_22 ]
// and that
// [ U    u    ] [  S    s   ] = [ U S +    u    s^T , U s + u    sigma ]
// [   upsilon ] [ s^T sigma ] = [       upsilon s^T ,    upsilon sigma ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  char diag=(
    dynamic_cast<const UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0
    ? 'N' : 'U');
  for (int k=0;k<m;k++) { // U_1 S_11
    if (k>0) {
      F77NAME(scopy)(k,S.addr(k,0),n,M->addr(0,k),1);
      F77NAME(strmv)('U','N',diag,k,U.addr(),m,M->addr(0,k),1); // U s
    }
    if (diag=='N') {
      F77NAME(sger)(k+1,k+1,float_one_,U.addr(0,k),1,S.addr(k,0),n,
        M->addr(),m);
    } else {
      F77NAME(sger)(k,k+1,float_one_,U.addr(0,k),1,S.addr(k,0),n,
        M->addr(),m);
      F77NAME(saxpy)(k+1,float_one_,S.addr(k,0),n,M->addr(k,0),m);
    }
  }
  if (n>m) {
    F77NAME(sgemm)('N','N',m,m,n-m,float_one_,U.addr(0,m),m,
      S.addr(m,0),n,float_one_,M->addr(),m); // U_2 S_21
    for (int j=m;j<n;j++) {
      F77NAME(scopy)(m,S.addr(j,0),n,M->addr(0,j),1);
    }
    F77NAME(strmm)('L','U','N',diag,m,n-m,float_one_,U.addr(),m,
      M->addr(0,m),m); // U_1 S_21^T
    F77NAME(ssymm)('R','L',m,n-m,float_one_,S.addr(m,m),n,
      U.addr(0,m),m,float_one_,M->addr(0,m),m); // U_2 S_22
  }
  return M;
}

template<> Matrix<float,float>*
SymmetricMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &L) const {
// compute by bordering: note that
// [ S_11 S_21^T ] [ L_1 ] = [ S_11 L_1 + S_21^T L_2 ]
// [ S_21  S_22  ] [ L_2 ] = [ S_21 L_1 +  S_22  L_2 ]
// and that
// [  S    s   ] [   L          ] = [  S  L +   s   ell^T ,   s   lambda ]
// [ s^T sigma ] [ ell^T lambda ] = [ s^T L + sigma ell^T , sigma lambda ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(1));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0 ? 'N' : 'U');
  for (int k=0;k<n;k++) { // S_11 L_1
    if (k>0) {
      F77NAME(scopy)(k,addr(k,0),m,M->addr(k,0),m);
      F77NAME(strmv)('L','T',diag,k,L.addr(),m,M->addr(k,0),m); // s^T L
    }
    if (diag=='N') {
      F77NAME(sger)(k+1,k+1,float_one_,addr(k,0),m,L.addr(k,0),m,
        M->addr(),m);
    } else {
      F77NAME(sger)(k+1,k,float_one_,addr(k,0),m,L.addr(k,0),m,
        M->addr(),m);
      F77NAME(saxpy)(k+1,float_one_,addr(k,0),m,M->addr(0,k),1);
    }
  }
  if (m>n) {
    F77NAME(sgemm)('T','N',n,n,m-n,float_one_,addr(n,0),m,
      L.addr(n,0),m,float_one_,M->addr(),m); // S_21^T L_2
    F77NAME(slacpy)('A',m-n,n,addr(n,0),m,M->addr(n,0),m);
    F77NAME(strmm)('R','L','N',diag,m-n,n,float_one_,L.addr(),m,
      M->addr(n,0),m); // S_21 L_1
    F77NAME(ssymm)('L','L',m-n,n,float_one_,addr(n,n),m,L.addr(n,0),m,
      float_one_,M->addr(n,0),m); // S_22 L_2
  }
  return M;
// compute by columns: note that
// [ S_11 S_21^T ] [     ] = [ S_21^T ell ]
// [ S_21  S_22  ] [ ell ] = [  S_22  ell ]
//bool L_non_unit=(dynamic_cast<const
//  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
//if (L_non_unit) {
//  for (int j=0;j<n;j++) {
//    if (j>0) {
//      F77NAME(sgemv)('T',m-j,j,float_one_,addr(j,0),m,L.addr(j,j),1,
//        float_zero_,M->addr(0,j),1);
//    }
//    F77NAME(ssymv)('L',m-j,float_one_,addr(j,j),m,L.addr(j,j),1,
//      float_zero_,M->addr(j,j),1);
//  }
//} else {
//  for (int j=0;j<n;j++) {
//    if (j>0) {
//      if (j+1<m) {
//        F77NAME(sgemv)('T',m-j-1,j,float_one_,addr(j+1,0),m,
//          L.addr(j+1,j),1,float_zero_,M->addr(0,j),1);
//      }
//      F77NAME(saxpy)(j,float_one_,addr(j,0),m,M->addr(0,j),1);
//    }
//    (*M)(j,j)=F77NAME(sdot)(m-j-1,addr(j+1,j),1,L.addr(j+1,j),1);
//    F77NAME(ssymv)('L',m-j-1,float_one_,addr(j+1,j+1),m,
//      L.addr(j+1,j),1,float_zero_,M->addr(j+1,j),1);
//    (*M)(j,j)+=(*this)(j,j);
//    F77NAME(saxpy)(m-j-1,float_one_,addr(j+1,j),1,M->addr(j+1,j),1);
//  }
//}
//return M;
}

template<> Matrix<float,float>* operator*(
const LowerTrapezoidalMatrix<float,float> &L,
const SymmetricMatrix<float,float> &S) {
// compute by bordering: note that
// [ L_1 ] S = [ L_1 S ]
// [ L_2 ]   = [ L_2 S ]
// and that
// [ lambda   ] [ sigma s^T ] = [ lambda sigma       , lambda s^T       ]
// [  ell   L ] [   s    S  ] = [  ell   sigma + L s ,   ell  s^T + L S ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,S.size(0));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0 ? 'N' : 'U');
  for (int k=n-1;k>=0;k--) { // L_1 S
    if (k<n-1) {
      F77NAME(scopy)(n-k-1,S.addr(k+1,k),1,M->addr(k+1,k),1);
      F77NAME(strmv)('L','N',diag,n-k-1,L.addr(k+1,k+1),m,
        M->addr(k+1,k),1); // L s
    }
    if (diag=='N') {
      F77NAME(sger)(n-k,n-k,float_one_,L.addr(k,k),1,S.addr(k,k),1,
        M->addr(k,k),m);
    } else {
      F77NAME(saxpy)(n-k,float_one_,S.addr(k,k),1,M->addr(k,k),m);
      F77NAME(sger)(n-k-1,n-k,float_one_,L.addr(k+1,k),1,S.addr(k,k),1,
        M->addr(k+1,k),m);
    }
  }
  if (m>n) { // L_2 S
    F77NAME(ssymm)('R','L',m-n,n,float_one_,S.addr(),n,L.addr(n,0),m,
      float_zero_,M->addr(n,0),m);
  }
  return M;
}

template<> SquareMatrix<float,float>*
SymmetricMatrix<float,float>::operator*(
const SquareMatrix<float,float> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  F77NAME(ssymm)('L','L',n,n,float_one_,addr(),n,S.addr(),n,
    float_zero_,T->addr(),n);
  return T;
}

template<> SquareMatrix<float,float>* operator*(
const SquareMatrix<float,float> &A,
const SymmetricMatrix<float,float> &B) {
// compute by bordering: note that
// [  M   v ] [  S    s   ] = [  M  S +  v s^T ,  M  s +  v sigma ]
// [ m^T mu ] [ s^T sigma ] = [ m^T S + mu s^T , m^T s + mu sigma ]
  int n=A.size(0);
  CHECK_SAME(n,B.size(0));
  SquareMatrix<float,float> *T=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
//for (int k=0;k<n;k++) {
//  if (k>0) {
//    F77NAME(sgemv)('N',k,k,float_one_,A.addr(),n,B.addr(k,0),n,
//      float_zero_,T->addr(0,k),1); // M s
//    F77NAME(ssymv)('L',k,float_one_,B.addr(),n,A.addr(k,0),n,
//      float_zero_,T->addr(k,0),n); // m^T S
//    (*T)(k,k)=F77NAME(sdot)(k,A.addr(k,0),n,B.addr(k,0),n); // m^T s
//  }
//  F77NAME(sger)(k+1,k+1,float_one_,A.addr(0,k),1,B.addr(k,0),n,
//    T->addr(),n);
//}
//return T;
// compute by bordering: note that
// [ mu v^T ] [ sigma s^T ] = [ mu sigma + v^T s , mu s^T + v^T S ]
// [  m  M  ] [   s    S  ] = [  m sigma +  M  s ,  m s^T +  M  S ]
  for (int k=n-1;k>=0;k--) {
    if (k+1<n) {
      F77NAME(sgemv)('N',n-k-1,n-k-1,float_one_,A.addr(k+1,k+1),n,
        B.addr(k+1,k),1,float_zero_,T->addr(k+1,k),1); // M s
      F77NAME(ssymv)('L',n-k-1,float_one_,B.addr(k+1,k+1),n,
        A.addr(k,k+1),n,float_zero_,T->addr(k,k+1),n);//v^T S = (S v)^T
      (*T)(k,k)=
        F77NAME(sdot)(n-k-1,A.addr(k,k+1),n,B.addr(k+1,k),1); // v^T s
    }
    F77NAME(sger)(n-k,n-k,float_one_,A.addr(k,k),1,B.addr(k,k),1,
      T->addr(k,k),n);
  }
  return T;
}

template<> Matrix<float,float>*
SymmetricMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  int m=size(0),n=M.size(1);
  CHECK_SAME(m,M.size(0));
  Matrix<float,float> *T=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(ssymm)('L','L',m,n,float_one_,addr(),m,M.addr(),m,
    float_zero_,T->addr(),m);
  return T;
}

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,const SymmetricMatrix<float,float> &S) {
  int m=M.size(0),n=S.size(0);
  CHECK_SAME(n,M.size(1));
  Matrix<float,float> *T=OPERATOR_NEW Matrix<float,float>(m,n);
  F77NAME(ssymm)('R','L',m,n,float_one_,S.addr(),n,M.addr(),m,
    float_zero_,T->addr(),m);
  return T;
}

template<> Vector<float,float>*
SymmetricMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<float,float> *w=OPERATOR_NEW Vector<float,float>(n);
  F77NAME(ssymv)('L',n,float_one_,addr(),n,v.addr(),1,float_zero_,
    w->addr(),1);
  return w;
}

// y := A * x * alpha + y * beta
template<> void SymmetricMatrix<float,float>::symv(float alpha,
const Vector<float,float> &x,float beta,Vector<float,float> &y)
const {
  int n=this->size(0);
  F77NAME(ssymv)('L',n,alpha,addr(),n,x.addr(),1,beta,y.addr(),1);
}

// A += x * alpha * x^T
template<> void SymmetricMatrix<float,float>::syr(float alpha,
const Vector<float,float> &x) {
  int n=this->size(0);
  F77NAME(ssyr)('L',n,alpha,x.addr(),1,addr(),n);
}

// A += x * alpha * y^T + y * alpha * x^T
template<> void SymmetricMatrix<float,float>::syr2(float alpha,
const Vector<float,float> &x,const Vector<float,float> &y) {
  int n=this->size(0);
  F77NAME(ssyr2)('L',n,alpha,x.addr(),1,y.addr(),1,addr(),n);
}

// C := A * alpha * B + C * beta
template<> void SymmetricMatrix<float,float>::symm(float alpha,
const Matrix<float,float> &B,float beta,Matrix<float,float> &C,
char side) const {
  int m=C.size(0),n=C.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,size(0));
    F77NAME(ssymm)('L','L',m,n,alpha,addr(),m,B.addr(),m,beta,
      C.addr(),m);
  } else {
    CHECK_SAME(n,size(0));
    F77NAME(ssymm)('R','L',m,n,alpha,addr(),n,B.addr(),m,beta,
      C.addr(),m);
  }
}

// C := A * alpha * A^T + C * beta
template<> void SymmetricMatrix<float,float>::syrk(float alpha,
const Matrix<float,float> &A,float beta,char transa) {
  int m=A.size(0),n=A.size(1);
  if (transa!='N' && transa!='n') {
    CHECK_SAME(n,size(0));
    F77NAME(ssyrk)('L',transa,n,m,alpha,A.addr(),m,beta,addr(),n);
  } else {
    CHECK_SAME(m,size(0));
    F77NAME(ssyrk)('L',transa,m,n,alpha,A.addr(),m,beta,addr(),m);
  }
}

// C := A * alpha * B^T + B * alpha + A^T + C * beta
template<> void SymmetricMatrix<float,float>::syr2k(float alpha,
const Matrix<float,float> &A,const Matrix<float,float> &B,
float beta,char transab) {
  int m=A.size(0),n=A.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (transab!='N' && transab!='n') {
    CHECK_SAME(n,size(0));
    F77NAME(ssyr2k)('L',transab,n,m,alpha,A.addr(),m,B.addr(),m,beta,
      addr(),n);
  } else {
    CHECK_SAME(m,size(0));
    F77NAME(ssyr2k)('L',transab,m,n,alpha,A.addr(),m,B.addr(),m,beta,
      addr(),m);
  }
}

/*
// y := abs(A) * abs(x) * alpha + abs(y) * beta
template<> void SymmetricMatrix<float,float>::syamv(float alpha,
const Vector<float,float> &x,float beta,Vector<float,float> &y)
const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,y.size());
  F77_NAME(sla_syamv)('L',n,alpha,addr(),n,x.addr(),1,beta,y.addr(),1);
}
*/

template<> float SymmetricMatrix<float,float>::equilibrate(
Vector<float,float> &s,float &scond) const {
  int n=size(0);
  CHECK_SAME(n,s.size());
  float amax=float_undefined_;
  int info;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  F77NAME(ssyequb)('L',n,addr(),n,s.addr(),scond,amax,work,info);
  CHECK_SAME(info,0);
  delete [] work;
  return amax;
}

template<> float SymmetricMatrix<float,float>::normFrobenius() const {
  int n=size(0);
  float *work=0;
  return F77NAME(slansy)('F','L',n,addr(),n,work);
}

template<> float SymmetricMatrix<float,float>::normInfinity() const {
  int n=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,n);
  float val=F77NAME(slansy)('I','L',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> float SymmetricMatrix<float,float>::normMaxEntry() const {
  int n=size(0);
  float *work=0;
  return F77NAME(slansy)('M','L',n,addr(),n,work);
}

template<> float SymmetricMatrix<float,float>::normOne() const {
  int n=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,n);
  float val=F77NAME(slansy)('O','L',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> float
SymmetricMatrix<float,float>::reciprocalConditionNumber() const {
  int n=size(0);

  SymmetricMatrix<float,float> *AF=
    OPERATOR_NEW SymmetricMatrix<float,float>(*this);
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  float w=HUGE_VAL;
  F77NAME(ssytrf)('L',n,AF->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0);

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(ssytrf)('L',n,AF->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;

  work=OPERATOR_NEW_BRACKET(float,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  float rcond;
  float anorm=F77NAME(slansy)('O','L',n,addr(),n,work);
  F77NAME(ssycon)('L',n,AF->addr(),n,ipiv,anorm,rcond,work,iwork,info);
  delete [] ipiv; ipiv=0;
  delete [] work; work=0;
  delete [] iwork; iwork=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SymmetricMatrix<float,float>*
SymmetricMatrix<float,float>::inverse() const {
  int n=size(0);
  SymmetricMatrix<float,float> *Ainv=
    OPERATOR_NEW SymmetricMatrix<float,float>(*this);
  int info;
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  float w=HUGE_VAL;
  F77NAME(ssytrf)('L',n,Ainv->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0);

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(ssytrf)('L',n,Ainv->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;

  work=OPERATOR_NEW_BRACKET(float,n);
  F77NAME(ssytri)('L',n,Ainv->addr(),n,ipiv,work,info);
  CHECK_SAME(info,0)

  delete [] ipiv;
  delete [] work;
  return Ainv;
}
*/

template<> Vector<float,float>*
SymmetricMatrix<float,float>::eigenvalues(
OrthogonalMatrix<float,float> *&Q) const {
  int n=size(0);
  if (Q!=0) CHECK_SAME(n,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Q->Matrix<float,float>::copyFrom('L',n,n,*this);
  Vector<float,float> *lambda =OPERATOR_NEW Vector<float,float>(n);
  float w;
  int lwork=-1,info;
  F77NAME(ssyev)(jobz,'L',n,Q->addr(),n,lambda->addr(),&w,lwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(ssyev)(jobz,'L',n,Q->addr(),n,lambda->addr(),work,lwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;

  return lambda;
}

template<> void SymmetricMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char) const {
  int n=size(0);
  SymmetricMatrix<float,float> AF(*this);
  int info;
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  x.copy(b);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  float w=HUGE_VAL;
  F77NAME(ssytrf)('L',n,AF.addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(ssytrf)('L',n,AF.addr(),n,ipiv,work,lwork,info);
  delete [] work; work=0;

  if (info==0) {
    F77NAME(ssytrs)('L',n,1,AF.addr(),n,ipiv,x.addr(),n,info);
  }
  CHECK_SAME(info,0)
  delete [] ipiv; ipiv=0;
}

template<> void SymmetricMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,char)
const {
  int n=size(0);
  SymmetricMatrix<float,float> AF(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int lwork=-1;
  float w=HUGE_VAL;
  int info;
  F77NAME(ssytrf)('L',n,AF.addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(ssytrf)('L',n,AF.addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;

  bool left_side=(side=='L' || side=='l');
  if (left_side) {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1))
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    X.copy(B);
    F77NAME(ssytrs)('L',n,1,AF.addr(),n,ipiv,X.addr(),n,info);
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0))
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    float *t=OPERATOR_NEW_BRACKET(float,n);
    for (int i=0;i<k;i++) {
      F77NAME(scopy)(n,B.addr(i,0),k,t,1);
      F77NAME(ssytrs)('L',n,1,AF.addr(),n,ipiv,t,n,info);
      F77NAME(scopy)(n,t,1,X.addr(i,0),k);
      CHECK_SAME(info,0)
    }
    delete [] t; t=0;
  }
  delete [] ipiv; ipiv=0;
}

template class SymmetricMatrix<float,float>;
template void testSymmetricMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveMatrix<float,float>&
SymmetricPositiveMatrix<float,float>::operator*=(float scalar) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(n-j,abs(scalar),addr(j,j),1);
  }
  return *this;
}

template<> SymmetricPositiveMatrix<float,float>&
SymmetricPositiveMatrix<float,float>::operator/=(float scalar) {
  int n=this->size(0);
  CHECK_NONZERO(scalar)
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(n-j,float_one_/abs(scalar),addr(j,j),1);
  }
  return *this;
}

template<> float SymmetricPositiveMatrix<float,float>::equilibrate(
Vector<float,float> &s,float &scond) const {
  int n=size(0);
  CHECK_SAME(n,s.size());
  float amax=float_undefined_;
  int info;
  F77NAME(spoequb)(n,addr(),n,s.addr(),scond,amax,info);
  CHECK_SAME(info,0);
  return amax;
}

template<> float
SymmetricPositiveMatrix<float,float>::reciprocalConditionNumber(
) const {
  int n=size(0);
  SymmetricPositiveMatrix<float,float> *AF=
  OPERATOR_NEW SymmetricPositiveMatrix<float,float>(*this);
  int info;
  F77NAME(spotrf)('L',n,AF->addr(),n,info);
  CHECK_SAME(info,0)

  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  float rcond;
  float anorm=F77NAME(slansy)('O','L',n,addr(),n,work);
  F77NAME(spocon)('L',n,addr(),n,anorm,rcond,work,iwork,info);
  delete [] work; work=0;
  delete [] iwork; iwork=0;
  delete AF; AF=0;
  return rcond;
}

/*
template<> SymmetricPositiveMatrix<float,float>*
SymmetricPositiveMatrix<float,float>::inverse() const {
  int n=size(0);
  SymmetricPositiveMatrix<float,float> *Ainv=
    OPERATOR_NEW SymmetricPositiveMatrix<float,float>(*this);
  int info;
  F77NAME(spotrf)('L',n,Ainv->addr(),n,info);
  CHECK_SAME(info,0);

  F77NAME(spotri)('L',n,Ainv->addr(),n,info);
  CHECK_SAME(info,0)
  return Ainv;
}
*/

template<> void SymmetricPositiveMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char) const {
  int n=size(0);
  SymmetricPositiveMatrix<float,float> AF(*this);
  int info;
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  x.copy(b);
  F77NAME(spotrf)('L',n,AF.addr(),n,info);
  CHECK_SAME(info,0)

  F77NAME(spotrs)('L',n,1,AF.addr(),n,x.addr(),n,info);
  CHECK_SAME(info,0)
}

template<> void SymmetricPositiveMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,char)
const {
  int n=size(0);
  SymmetricPositiveMatrix<float,float> AF(*this);
  int info;
  F77NAME(spotrf)('L',n,AF.addr(),n,info);
  CHECK_SAME(info,0)

  bool left_side=(side=='L' || side=='l');
  if (left_side) {
    int k=B.size(1);
    CHECK_SAME(k,X.size(1))
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    X.copy(B);
    F77NAME(spotrs)('L',n,1,AF.addr(),n,X.addr(),n,info);
    CHECK_SAME(info,0)
  } else {
    int k=B.size(0);
    CHECK_SAME(k,X.size(0))
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    float *t=OPERATOR_NEW_BRACKET(float,n);
    for (int i=0;i<k;i++) {
      F77NAME(scopy)(n,B.addr(i,0),k,t,1);
      F77NAME(spotrs)('L',n,1,AF.addr(),n,t,n,info);
      F77NAME(scopy)(n,t,1,X.addr(i,0),k);
      CHECK_SAME(info,0)
    }
    delete [] t; t=0;
  }
}

template class SymmetricPositiveMatrix<float,float>;
template void testSymmetricPositiveMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "BandMatrix.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float TridiagonalMatrix<float,float>::safety_ =
  float_zero_;
template<> const float TridiagonalMatrix<float,float>::outofbounds_ =
  float_zero_;
template<> const float TridiagonalMatrix<float,float>::undefined_ =
  HUGE_VAL;

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,float> *M=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  F77NAME(scopy)(n,D->addr(),1,M->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,M->addr(1,0),n+1);
  F77NAME(scopy)(n-1,U->addr(),1,M->addr(0,1),n+1);
  return M;
}

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::operator+(
const SymmetricMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int k=0;k<n;k++) {
    F77NAME(scopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) F77NAME(scopy)(n-k-1,M.addr(k+1,k),1,S->addr(k,k+1),n);
  }
  F77NAME(saxpy)(n,float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(scopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)=float_one_;
      if (k<n-1) F77NAME(scopy)(n-k-1,M.addr(k+1,k),1,S->addr(k+1,k),1);
    }
  }
  F77NAME(saxpy)(n,float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<float,float>*
TridiagonalMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<float,float> *S=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(scopy)(k+1,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(scopy)(k,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)=float_one_;
    }
  }
  F77NAME(saxpy)(n,float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  S->copy(M);
  F77NAME(saxpy)(n,float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,U->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::operator-(
const SymmetricMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(scopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    F77NAME(saxpy)(n-k,float_mone_,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) {
      F77NAME(saxpy)(n-k-1,float_mone_,M.addr(k+1,k),1,S->addr(k,k+1),n);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const SymmetricMatrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) { 
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int k=0;k<n;k++) {
    F77NAME(scopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    if (k<n-1) F77NAME(scopy)(n-k-1,M.addr(k+1,k),1,S->addr(k,k+1),n);
  }
  F77NAME(saxpy)(n,float_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0);
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(scopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(saxpy)(n-k,float_mone_,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)-=float_one_;
      if (k<n-1) {
        F77NAME(saxpy)(n-k-1,float_mone_,M.addr(k+1,k),1,
          S->addr(k+1,k),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const LowerTrapezoidalMatrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) { 
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(scopy)(n-k,M.addr(k,k),1,S->addr(k,k),1);
    } else {
      (*S)(k,k)=float_one_;
      if (k<n-1) {
        F77NAME(scopy)(n-k-1,M.addr(k+1,k),1,S->addr(k+1,k),1);
      }
    }
  }
  F77NAME(saxpy)(n,float_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<float,float>*
TridiagonalMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<float,float> *S=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&M)==0);
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(scopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(saxpy)(k+1,float_mone_,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(saxpy)(k,float_mone_,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)-=float_one_;
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<float,float>* operator-(
const UpperTrapezoidalMatrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  UpperHessenbergMatrix<float,float> *S=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&M)==0);
  for (int k=0;k<n;k++) {
    if (M_non_unit) {
      F77NAME(scopy)(k+1,M.addr(0,k),1,S->addr(0,k),1);
    } else {
      if (k>0) {
        F77NAME(scopy)(k,M.addr(0,k),1,S->addr(0,k),1);
      }
      (*S)(k,k)=float_one_;
    }
  }
  F77NAME(saxpy)(n,float_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(scopy)(n-1,U->addr(),1,S->addr(0,1),n+1);
  *S-=M;
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const Matrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  S->copy(M);
  F77NAME(saxpy)(n,float_mone_,T.addr(0,0),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(1,0),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(0,1),1,S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::operator*(
const SymmetricMatrix<float,float> &M) const {
// compute by bordering: note that
// [    tau     upsilon e_0^T ] [ sigma s^T ]
// [ e_0 lambda       T       ] [   s    S  ]
//   = [ tau sigma + upsilon e_0^T s , tau s^T + upsilon e_0^T S ]
//   = [ e_0 lambda sigma +      T s , e_0 lambda s^T + T S      ]
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  (*S)(n-2,n-2)=(*D)[n-2]*M(n-2,n-2)+(*U)[n-2]*M(n-1,n-2);
  (*S)(n-1,n-2)=(*L)[n-2]*M(n-2,n-2)+(*D)[n-1]*M(n-1,n-2);
  (*S)(n-2,n-1)=(*D)[n-2]*M(n-1,n-2)+(*U)[n-2]*M(n-1,n-1);
  (*S)(n-1,n-1)=(*L)[n-2]*M(n-1,n-2)+(*D)[n-1]*M(n-1,n-1);
  for (int k=n-3;k>=0;k--) {
    F77NAME(sgtmv)(n-k-1,float_one_,L->addr(k+1),D->addr(k+1),
      U->addr(k+1),M.addr(k+1,k),1,float_zero_,S->addr(k+1,k),1); // T s
    F77NAME(saxpy)(n-k-1,(*U)[k],M.addr(k+1,k+1),1,S->addr(k,k+1),n);
      // upsilon e_0^T S
    (*S)(k,k)=(*U)[k]*M(k+1,k); // upsilon e_0^T s
    F77NAME(saxpy)(n-k,(*D)[k],M.addr(k,k),1,S->addr(k,k),n);
      // tau [ sigma , s^T ]
    F77NAME(saxpy)(n-k,(*L)[k],M.addr(k,k),1,S->addr(k+1,k),n);
      // e_0 lambda [ sigma , s^T ]
  }
  return S;
}

template<> SquareMatrix<float,float>* operator*(
const SymmetricMatrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) {
// compute by bordering: note that
// [ sigma s^T ] [    tau     upsilon e_0^T ]
// [   s    S  ] [ e_0 lambda       T       ]
//   = [ sigma tau + s^T e_0 lambda , sigma upsilon e_0^T + s^T T ]
//   = [     s tau +   S e_0 lambda ,     s upsilon e_0^T +   S T ]
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  (*S)(n-2,n-2)=M(n-2,n-2)*T(n-2,n-2)+M(n-1,n-2)*T(n-1,n-2);
  (*S)(n-1,n-2)=M(n-1,n-2)*T(n-2,n-2)+M(n-1,n-1)*T(n-1,n-2);
  (*S)(n-2,n-1)=M(n-2,n-2)*T(n-2,n-1)+M(n-1,n-2)*T(n-1,n-1);
  (*S)(n-1,n-1)=M(n-1,n-2)*T(n-2,n-1)+M(n-1,n-1)*T(n-1,n-1);
  for (int k=n-3;k>=0;k--) {
    F77NAME(sgtmv)(n-k-1,float_one_,T.addr(k+1,k+2),T.addr(k+1,k+1),
      T.addr(k+2,k+1),M.addr(k+1,k),1,float_zero_,S->addr(k,k+1),n);
    F77NAME(saxpy)(n-k-1,T(k+1,k),M.addr(k+1,k+1),1,S->addr(k+1,k),1);
      // S e_0 lambda
    (*S)(k,k)=M(k+1,k)*T(k+1,k); // s^T e_0 lambda
    F77NAME(saxpy)(n-k,T(k,k),M.addr(k,k),1,S->addr(k,k),1);
      // [ sigma ] tau
      // [   s   ]
    F77NAME(saxpy)(n-k,T(k,k+1),M.addr(k,k),1,S->addr(k,k+1),1);
      // [ sigma ] upsilon
      // [   s   ]
  }
  return S;
}

template<> Matrix<float,float>*
TridiagonalMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &M) const {
// compute by columns: note that
// [      T_11       ,e_j upsilon e_0^T ][ 0 ] = [ e_j upsilon e_0^T ell ]
// [ e_0 lambda e_j^T,      T_22        ][ell] = [             T_22  ell ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0);
  if (M_non_unit) {
    F77NAME(sgtmv)(m,float_one_,L->addr(),D->addr(),U->addr(),
      M.addr(),1,float_zero_,S->addr(),1);
    for (int j=1;j<n;j++) {
      (*S)(j-1,j)=(*U)[j-1]*M(j,j); // e_j upsilon e_0^T ell
      if (j<m-1) { // T_22 ell
        F77NAME(sgtmv)(m-j,float_one_,L->addr(j),D->addr(j),
          U->addr(j),M.addr(j,j),1,float_one_,S->addr(j,j),1);
      } else (*S)(j,j)=(*D)[j]*M(j,j);
    }
  } else {
// note that
// [     tau    , upsilon e_0^T ] [ 1 ] = [ tau + upsilon e_0^T ell ]
// [ e_0 lambda ,        T      ] [ell] = [ e_0 lambda      + T ell ]
    (*S)(0,0)=(*D)[0]+(*U)[0]*M(1,0);
    (*S)(1,0)=(*L)[0];
    F77NAME(sgtmv)(m-1,float_one_,L->addr(1),D->addr(1),
      U->addr(1),M.addr(1,0),1,float_one_,S->addr(1,0),1);
    for (int j=1;j<n;j++) {
      (*S)(j-1,j)=(*U)[j-1]; // e_j upsilon e_0^T ell
      (*S)(j,j)=(*D)[j]; // tau
      if (j<m-1) {
        (*S)(j,j)+=(*U)[j]*M(j+1,j); // upsilon e_0^T ell
        (*S)(j+1,j)=(*L)[j]; // e_0 lambda
      }
      if (j<m-2) {
        F77NAME(sgtmv)(m-j-1,float_one_,L->addr(j+1),
          D->addr(j+1),U->addr(j+1),M.addr(j+1,j),1,float_one_,
          S->addr(j+1,j),1);
      } else if (j<m-1) (*S)(j+1,j)+=(*D)[j+1]*M(j+1,j);
    }
  }
  return S;
}

template<> Matrix<float,float>* operator*(
const LowerTrapezoidalMatrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) {
// compute by columns: note that
// L T e_j = L e_{j-1} T_{j-1,j} + L e_j T_{j,j} + L e_{j+1{ T){j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(m-j+1,T(j-1,j),M.addr(j-1,j-1),1,
          S->addr(j-1,j),1);
      }
      F77NAME(saxpy)(m-j,T(j,j),M.addr(j,j),1,S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(m-j-1,T(j+1,j),M.addr(j+1,j+1),1,S->addr(j+1,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        (*S)(j-1,j)=T(j-1,j);
        F77NAME(saxpy)(m-j,T(j-1,j),M.addr(j,j-1),1,S->addr(j,j),1);
      }
      (*S)(j,j)+=T(j,j);
      if (j<m-1) {
        F77NAME(saxpy)(m-j-1,T(j,j),M.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      if (j<n-1) {
        (*S)(j+1,j)+=T(j+1,j);
        if (j<m-2) {
          F77NAME(saxpy)(m-j-2,T(j+1,j),M.addr(j+2,j+1),1,
            S->addr(j+2,j),1);
        }
      }
    }
  }
  return S;
}

template<> Matrix<float,float>*
TridiagonalMatrix<float,float>::operator*(
const UpperTrapezoidalMatrix<float,float> &M) const {
// compute by columns: note that
// T [ U_1 , U_2 ] = [ T U_1 , T U_2 ]
// and that
// [      T_11       ,e_j upsilon e_0^T ][ u ] = [             T_11 u ]
// [ e_0 lambda e_j^T,      T_22        ][ 0 ] = [ e_0 lambda e_j^T u ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&M)==0);
  if (M_non_unit) {
    (*S)(0,0)=(*D)[0]*M(0,0);
    (*S)(1,0)=(*L)[0]*M(0,0);
    for (int j=1;j<m-1;j++) { // T U_1
      F77NAME(sgtmv)(j+1,float_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,float_zero_,S->addr(0,j),1); // T_11 u
      (*S)(j+1,j)=(*L)[j]*M(j,j); // e_0 lambda e_j^T u
    }
    for (int j=m-1;j<n;j++) { // T U_2
      F77NAME(sgtmv)(m,float_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,float_zero_,S->addr(0,j),1);
    }
  } else {
    (*S)(0,0)=(*D)[0];
    (*S)(1,0)=(*L)[0];
    for (int j=1;j<m;j++) {
      F77NAME(sgtmv)(j,float_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,float_zero_,S->addr(0,j),1);
      (*S)(j-1,j)+=(*U)[j-1];
      (*S)(j,j)=(*L)[j-1]*M(j-1,j)+(*D)[j];
      if (j<m-1) (*S)(j+1,j)=(*L)[j];
    }
    for (int j=m;j<n;j++) { // T U_2
      F77NAME(sgtmv)(m,float_one_,L->addr(),D->addr(),U->addr(),
        M.addr(0,j),1,float_zero_,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<float,float>* operator*(
const UpperTrapezoidalMatrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) {
// compute by columns: note that
// U T e_j = U e_{j-1} T_{j-1,j} + U e_j T_{j,j} + U e_{j+1{ T){j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j,m),T(j-1,j),M.addr(0,j-1),1,S->addr(0,j),1);
      }
      F77NAME(saxpy)(min(j+1,m),T(j,j),M.addr(0,j),1,S->addr(0,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(min(j+2,m),T(j+1,j),M.addr(0,j+1),1,
          S->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j-1,m),T(j-1,j),M.addr(0,j-1),1,
          S->addr(0,j),1);
        if (j<=m) (*S)(j-1,j)=T(j-1,j);
      }
      F77NAME(saxpy)(min(j,m),T(j,j),M.addr(0,j),1,S->addr(0,j),1);
      if (j<m) (*S)(j,j)+=T(j,j);
      if (j<n-1) {
        F77NAME(saxpy)(min(j+1,m),T(j+1,j),M.addr(0,j+1),1,
          S->addr(0,j),1);
        if (j<m-1) (*S)(j+1,j)+=T(j+1,j);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
TridiagonalMatrix<float,float>::operator*(
const SquareMatrix<float,float> &M) const {
  int m=M.size(0);
  CHECK_SAME(m,size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(m,float_zero_);
  for (int j=0;j<m;j++) {
    F77NAME(sgtmv)(m,float_one_,L->addr(),D->addr(),U->addr(),
      M.addr(0,j),1,float_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator*(
const SquareMatrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) {
  int m=M.size(0);
  CHECK_SAME(m,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(m,float_zero_);
  for (int j=0;j<m;j++) {
    for (int k=max(0,j-1);k<=min(m-1,j+1);k++) {
      F77NAME(saxpy)(m,T(k,j),M.addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<float,float>*
TridiagonalMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(sgtmv)(m,float_one_,L->addr(),D->addr(),U->addr(),
      M.addr(0,j),1,float_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,
const TridiagonalMatrix<float,float> &T) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
      F77NAME(saxpy)(m,T(k,j),M.addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Vector<float,float>*
TridiagonalMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<float,float> *w=
    OPERATOR_NEW Vector<float,float>(n,float_zero_);
  F77NAME(sgtmv)(n,float_one_,L->addr(),D->addr(),U->addr(),
    v.addr(),1,float_zero_,w->addr(),1);
  return w;
}

template<> float TridiagonalMatrix<float,float>::normFrobenius()
const {
  return F77NAME(slangt)('F',dim,L->addr(),D->addr(),U->addr());
}

template<> float TridiagonalMatrix<float,float>::normInfinity()
const {
  return F77NAME(slangt)('I',dim,L->addr(),D->addr(),U->addr());
}

template<> float TridiagonalMatrix<float,float>::normMaxEntry()
const {
  return F77NAME(slangt)('M',dim,L->addr(),D->addr(),U->addr());
}

template<> float TridiagonalMatrix<float,float>::normOne() const {
  return F77NAME(slangt)('O',dim,L->addr(),D->addr(),U->addr());
}

template<> float
TridiagonalMatrix<float,float>::reciprocalConditionNumber(char norm)
const {
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  Vector<float,float> *UF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *UF2=
    OPERATOR_NEW Vector<float,float>(max(0,dim-2));
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(sgttrf)(dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),ipiv,
    info);
  float rcond=undefined_;
  float *work=OPERATOR_NEW_BRACKET(float,2*dim);
  int *iwork=OPERATOR_NEW_BRACKET(int,dim);
  float anorm=F77NAME(slangt)(norm,dim,L->addr(),D->addr(),U->addr());
  F77NAME(sgtcon)(norm,dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
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

template<> void TridiagonalMatrix<float,float>::gtmv(float alpha,
const Vector<float,float> &x,float beta,Vector<float,float> &b,
char trans) const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  if (trans=='N' || trans=='n') {
    F77NAME(sgtmv)(n,alpha,L->addr(),D->addr(),U->addr(),
      x.addr(),1,beta,b.addr(),1);
  } else {
    F77NAME(sgtmv)(n,alpha,U->addr(),D->addr(),L->addr(),
      x.addr(),1,beta,b.addr(),1);
  }
}

template<> void TridiagonalMatrix<float,float>::gtmm(float alpha,
const Matrix<float,float> &X,float beta,Matrix<float,float> &B,
char side,char trans) const {
  if (side=='L' || side=='l') {
    int m=size(0),n=X.size(1);
    CHECK_SAME(n,B.size(1));
    CHECK_SAME(m,X.size(0));
    CHECK_SAME(m,B.size(0));
    if (trans=='N' || trans=='n') {
      for (int j=0;j<n;j++) {
        F77NAME(sgtmv)(m,alpha,L->addr(),D->addr(),U->addr(),
          X.addr(0,j),1,beta,B.addr(0,j),1);
      }
    } else {
      for (int j=0;j<n;j++) {
        F77NAME(sgtmv)(m,alpha,U->addr(),D->addr(),L->addr(),
          X.addr(0,j),1,beta,B.addr(0,j),1);
      }
    }
  } else {
    int m=X.size(0),n=X.size(1);
    CHECK_SAME(n,size(0));
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,B.size(1));
    if (abs(beta)==float_zero_) B=float_zero_;
    else if (beta!=float_one_) B*=beta;
    if (trans=='N' || trans=='n') {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(saxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<=min(n-1,j+1);k++) {
          F77NAME(saxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

template<> void TridiagonalMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  Vector<float,float> *UF=OPERATOR_NEW Vector<float,float>(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  x.copy(b);
  int info;
  if (trans!='N' && trans!='n') {
    F77NAME(sgtsv)(dim,1,UF->addr(),DF->addr(),LF->addr(),x.addr(),dim,
      info);
  } else {
    F77NAME(sgtsv)(dim,1,LF->addr(),DF->addr(),UF->addr(),x.addr(),dim,
      info);
  }
  CHECK_SAME(info,0);
  delete UF; UF=0;
  delete DF; DF=0;
  delete LF; LF=0;
}

template<> void TridiagonalMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  Vector<float,float> *UF=OPERATOR_NEW Vector<float,float>(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*U);
  int info;
  X.copy(B);
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    if (trans!='N' && trans!='n') {
      F77NAME(sgtsv)(m,n,UF->addr(),DF->addr(),LF->addr(),X.addr(),m,
        info);
    } else {
      F77NAME(sgtsv)(m,n,LF->addr(),DF->addr(),UF->addr(),X.addr(),m,
        info);
    }
  } else {
    CHECK_SAME(dim,n);
    Vector<float,float> *UF2=OPERATOR_NEW Vector<float,float>(dim-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
    F77NAME(sgttrf)(n,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
      ipiv,info);
    CHECK_TEST(info==0);
    float *t=OPERATOR_NEW_BRACKET(float,n);
    char rtrans=(trans!='N' && trans!='n' ? 'N' : 'T');
    for (int i=0;i<m;i++) {
      F77NAME(scopy)(n,B.addr(i,0),m,t,1);
      F77NAME(sgttrs)(rtrans,n,1,LF->addr(),DF->addr(),UF->addr(),
        UF2->addr(),ipiv,t,n,info);
      F77NAME(scopy)(n,t,1,X.addr(i,0),m);
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

template class TridiagonalMatrix<float,float>;
template void testTridiagonalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> float
  SymmetricTridiagonalMatrix<float,float>::safety_=float_zero_;
template<> const float
  SymmetricTridiagonalMatrix<float,float>::outofbounds_=
  float_zero_;
template<> const float
  SymmetricTridiagonalMatrix<float,float>::undefined_=
  HUGE_VAL;

template<> SymmetricTridiagonalMatrix<float,float>::
SymmetricTridiagonalMatrix(int n,float d) : dim(n) {
  L=OPERATOR_NEW Vector<float,float>(n-1,d);
  D=OPERATOR_NEW Vector<float,float>(n,d);
}

template<> SquareMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,float> *M=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  F77NAME(scopy)(n,D->addr(),1,M->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,M->addr(1,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,M->addr(0,1),n+1);
  return M;
}

template<> void SymmetricTridiagonalMatrix<float,float>::fillWith(
float scalar) {
  *L=scalar; *D=scalar;
}

template<> float
SymmetricTridiagonalMatrix<float,float>::upperDiagonalValue(int i)
const {
  return (*L)[i];
}

template<> SymmetricTridiagonalMatrix<float,float>&
SymmetricTridiagonalMatrix<float,float>::operator=(float scalar) {
  *L=scalar; *D=scalar; return *this;
}

template<> SymmetricTridiagonalMatrix<float,float>&
SymmetricTridiagonalMatrix<float,float>::operator*=(float scalar) {
  int n=this->size(0);
  (*D)*=scalar;
  (*L)*=scalar;
  return *this;
}

template<> SymmetricTridiagonalMatrix<float,float>&
SymmetricTridiagonalMatrix<float,float>::operator/=(float scalar) {
  int n=this->size(0);
  CHECK_NONZERO(scalar)
  (*D)/=scalar;
  (*L)/=scalar;
  return *this;
}

template<> TridiagonalMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator*(float d) const {
  int n=size(0);
  TridiagonalMatrix<float,float> *S=
    OPERATOR_NEW TridiagonalMatrix<float,float>(n);
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(0,1),1);
  *S*=d;
  return S;
}

template<> TridiagonalMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator/(float d) const {
  int n=size(0);
  TridiagonalMatrix<float,float> *S=
    OPERATOR_NEW TridiagonalMatrix<float,float>(n);
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(0,1),1);
  *S/=d;
  return S;
}

template<> TridiagonalMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator+(
const TridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,float> *S=
    OPERATOR_NEW TridiagonalMatrix<float,float>(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_one_,D->addr(),1,S->addr(0,0),1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(1,0),1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(0,1),1);
  return S;
}

template<> SymmetricMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator+(
const SymmetricMatrix<float,float> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricMatrix<float,float> *T=
    OPERATOR_NEW SymmetricMatrix<float,float>(n);
  T->copy(S);
  F77NAME(saxpy)(n,float_one_,D->addr(),1,T->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,T->addr(1,0),n+1);
  return T;
}

template<> SquareMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0); 
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=float_one_;
      if (j<n-1) {
        F77NAME(scopy)(n-j-1,M.addr(j+1,j),1,S->addr(j+1,j),1);
      }
    }
  }
  F77NAME(saxpy)(n,float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0); 
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(j+1,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(j,U.addr(0,j),1,H->addr(0,j),1);
      }
      (*H)(j,j)=float_one_;
    }
  }
  F77NAME(saxpy)(n,float_one_,D->addr(),1,H->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,H->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  S->copy(M);
  F77NAME(saxpy)(n,float_one_,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,S->addr(0,1),n+1);
  return S;
}

template<> SymmetricTridiagonalMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator-(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<float,float> *S=
    OPERATOR_NEW SymmetricTridiagonalMatrix<float,float>(n);
  S->copy(*this);
  F77NAME(saxpy)(n,float_mone_,T.D->addr(),1,S->D->addr(),1);
  F77NAME(saxpy)(n-1,float_mone_,T.L->addr(),1,S->L->addr(),1);
  return S;
}

template<> TridiagonalMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator-(
const TridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,float> *S=
    OPERATOR_NEW TridiagonalMatrix<float,float>(n);
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(0,1),1);
  F77NAME(saxpy)(n,float_mone_,T.addr(0,0),1,S->addr(0,0),1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(1,0),1,S->addr(1,0),1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(0,1),1,S->addr(0,1),1);
  return S;
}

template<> TridiagonalMatrix<float,float>* operator-(
const TridiagonalMatrix<float,float> &T,
const SymmetricTridiagonalMatrix<float,float> &St) {
  int n=St.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,float> *S=
    OPERATOR_NEW TridiagonalMatrix<float,float>(n);
  F77NAME(scopy)(n,T.addr(0,0),1,S->addr(0,0),1);
  F77NAME(scopy)(n-1,T.addr(1,0),1,S->addr(1,0),1);
  F77NAME(scopy)(n-1,T.addr(0,1),1,S->addr(0,1),1);
  F77NAME(saxpy)(n,float_mone_,St.diagonalAddr(0),1,S->addr(0,0),1);
  F77NAME(saxpy)(n-1,float_mone_,St.lowerDiagonalAddr(0),1,
    S->addr(1,0),1);
  F77NAME(saxpy)(n-1,float_mone_,St.lowerDiagonalAddr(0),1,
    S->addr(0,1),1);
  return S;
}

template<> SymmetricMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator-(
const SymmetricMatrix<float,float> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricMatrix<float,float> *T=
    OPERATOR_NEW SymmetricMatrix<float,float>(n,float_zero_);
  F77NAME(scopy)(n,D->addr(),1,T->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,T->addr(1,0),n+1);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,S.addr(j,j),1,T->addr(j,j),1);
  }
  return T;
}

template<> SymmetricMatrix<float,float>* operator-(
const SymmetricMatrix<float,float> &S,
const SymmetricTridiagonalMatrix<float,float> &T) {
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,float> *M=
    OPERATOR_NEW SymmetricMatrix<float,float>(n);
  M->copy(S);
  F77NAME(saxpy)(n,float_mone_,T.diagonalAddr(0),1,M->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    M->addr(1,0),n+1);
  return M;
}

template<> SquareMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0); 
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(0,1),n+1);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(n-j,float_mone_,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)-=float_one_;
      if (j<n-1) {
        F77NAME(saxpy)(n-j-1,float_mone_,M.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const LowerTrapezoidalMatrix<float,float> &M,
const SymmetricTridiagonalMatrix<float,float> &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0); 
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=float_one_;
      if (j<n-1) {
        F77NAME(scopy)(n-j-1,M.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  F77NAME(saxpy)(n,float_mone_,T.diagonalAddr(0),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(0,1),n+1);
  return S;
}

template<> UpperHessenbergMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0); 
  F77NAME(scopy)(n,D->addr(),1,H->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,H->addr(1,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,H->addr(0,1),n+1);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(j+1,float_mone_,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(j,float_mone_,U.addr(0,j),1,H->addr(0,j),1);
      }
      (*H)(j,j)-=float_one_;
    }
  }
  return H;
}

template<> UpperHessenbergMatrix<float,float>* operator-(
const UpperTrapezoidalMatrix<float,float> &U,
const SymmetricTridiagonalMatrix<float,float> &T) {
  int n=T.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0); 
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(j+1,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(scopy)(j,U.addr(0,j),1,H->addr(0,j),1);
      }
      (*H)(j,j)=float_one_;
    }
  }
  F77NAME(saxpy)(n,float_mone_,T.diagonalAddr(0),1,H->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  F77NAME(scopy)(n,D->addr(),1,S->addr(0,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(1,0),n+1);
  F77NAME(scopy)(n-1,L->addr(),1,S->addr(0,1),n+1);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n,float_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const Matrix<float,float> &M,
const SymmetricTridiagonalMatrix<float,float> &T) {
  int n=T.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  S->copy(M);
  F77NAME(saxpy)(n,float_mone_,T.diagonalAddr(0),1,S->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(0,1),n+1);
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator*(
const SymmetricMatrix<float,float> &M) const {
// compute by bordering: note that
// [     tau    , lambda e_0^T ] [ sigma S^T ]
// [ e_0 lambda ,      T       ] [   S    S  ]
//   = [ tau sigma + lambda e_0^T s , tau s^T        + lambda e_0^T S ]
//   = [ e_0 lambda sigma +     T s , e_0 lambda s^T +            T S ]
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  (*S)(n-2,n-2)=(*D)[n-2]*M(n-2,n-2)+(*L)[n-2]*M(n-1,n-2);
  (*S)(n-1,n-2)=(*L)[n-2]*M(n-2,n-2)+(*D)[n-1]*M(n-1,n-2);
  (*S)(n-2,n-1)=(*D)[n-2]*M(n-1,n-2)+(*L)[n-2]*M(n-1,n-1);
  (*S)(n-1,n-1)=(*L)[n-2]*M(n-1,n-2)+(*D)[n-1]*M(n-1,n-1);
  for (int k=n-3;k>=0;k--) {
    F77NAME(sstmv)(n-k-1,float_one_,L->addr(k+1),D->addr(k+1),
      M.addr(k+1,k),1,float_zero_,S->addr(k+1,k),1); // T s
    F77NAME(saxpy)(n-k-1,(*L)[k],M.addr(k+1,k+1),1,S->addr(k,k+1),n);
      // lambda e_0^T S = lambda ( S e_0 )^T
    (*S)(k,k)=(*L)[k]*M(k+1,k); // lambda e_0^T s
    F77NAME(saxpy)(n-k,(*D)[k],M.addr(k,k),1,S->addr(k,k),n);
      // tau [ sigma , s^T ]
    F77NAME(saxpy)(n-k,(*L)[k],M.addr(k,k),1,S->addr(k+1,k),n);
      // e_0 lambda [ sigma , s^T ]
  }
  return S;
}

template<> SquareMatrix<float,float>* operator*(
const SymmetricMatrix<float,float> &M,
const SymmetricTridiagonalMatrix<float,float> &T) {
// compute by bordering: note that
// [ sigma s^T ] [     tau    , lambda e_0^T ]
// [   s    S  ] [ e_0 lambda ,           T  ]
//   = [ sigma tau + s^T e_0 lambda , sigma lambda e_0^T + s^T T ]
//   = [   s   tau +  S  e_o lambda ,   s   lambda e_0^T +  S  T ]
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  (*S)(n-2,n-2)=M(n-2,n-2)*T.diagonalValue(n-2)
               +M(n-1,n-2)*T.lowerDiagonalValue(n-2);
  (*S)(n-1,n-2)=M(n-1,n-2)*T.diagonalValue(n-2)
               +M(n-1,n-1)*T.lowerDiagonalValue(n-2);
  (*S)(n-2,n-1)=M(n-2,n-2)*T.upperDiagonalValue(n-2)
               +M(n-1,n-2)*T.diagonalValue(n-1);
  (*S)(n-1,n-1)=M(n-1,n-2)*T.upperDiagonalValue(n-2)
               +M(n-1,n-1)*T.diagonalValue(n-1);
  float *t=OPERATOR_NEW_BRACKET(float,n);
  for (int k=n-3;k>=0;k--) {
    F77NAME(sstmv)(n-k-1,float_one_,T.lowerDiagonalAddr(k+1),
      T.diagonalAddr(k+1),M.addr(k+1,k),1,float_zero_,t,1);
    F77NAME(scopy)(n-k-1,t,1,S->addr(k,k+1),n); // s^T T = ( T s )^T
    F77NAME(saxpy)(n-k-1,T.lowerDiagonalValue(k),M.addr(k+1,k+1),1,
      S->addr(k+1,k),1); // S e_0 lambda
    (*S)(k,k)=M(k+1,k)*T.lowerDiagonalValue(k); // s^T e_0 lambda
    F77NAME(saxpy)(n-k,T.diagonalValue(k),M.addr(k,k),1,S->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(saxpy)(n-k,T.upperDiagonalValue(k),M.addr(k,k),1,
      S->addr(k,k+1),1);
      // [sigma ] bar(lambda)
      // [  s   ]
  }
  delete [] t; t=0;
  return S;
}

template<> Matrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &M) const {
// compute by columns: note that
// [     S_11        , e_j sigma e_0^T ] [ 0 ] = [ e_j sigma e_0^T ell ]
// [ e_0 sigma e_j^T ,    S_22         ] [ell] = [            S_22 ell ]
  int m=M.size(0),n=M.size(1);;
  CHECK_SAME(m,size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) (*S)(j-1,j)=upperDiagonalValue(j-1)*M(j,j);
        // sigma e_0^T ell
      if (j+1<m) { // S_22 ell
        F77NAME(sstmv)(m-j,float_one_,L->addr(j),D->addr(j),
          M.addr(j,j),1,float_zero_,S->addr(j,j),1);
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
          F77NAME(sstmv)(m-j-1,float_one_,L->addr(j+1),
            D->addr(j+1),M.addr(j+1,j),1,float_one_,
            S->addr(j+1,j),1);
        } else (*S)(j+1,j)+=diagonalValue(j+1)*M(j+1,j);
      }
    }
  }
  return S;
}

template<> Matrix<float,float>* operator*(
const LowerTrapezoidalMatrix<float,float> &M,
const SymmetricTridiagonalMatrix<float,float> &T) {
// compute by columns: note that
// L T e_j = L e_{j-1} T_{j-1,j} + L e_j T_{j,j} + L e_{j+1} T_{j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(m-j+1,T.upperDiagonalValue(j-1),M.addr(j-1,j-1),1,
          S->addr(j-1,j),1);
      }
      F77NAME(saxpy)(m-j,T.diagonalValue(j),M.addr(j,j),1,S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(m-j-1,T.lowerDiagonalValue(j),M.addr(j+1,j+1),1,
          S->addr(j+1,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        (*S)(j-1,j)=T.upperDiagonalValue(j-1);
        F77NAME(saxpy)(m-j,T.upperDiagonalValue(j-1),M.addr(j,j-1),1,
          S->addr(j,j),1);
      }
      (*S)(j,j)+=T.diagonalValue(j);
      if (j<m-1) {
        F77NAME(saxpy)(m-j-1,T.diagonalValue(j),M.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
      if (j<n-1) {
        (*S)(j+1,j)+=T.lowerDiagonalValue(j);
        if (j<m-2) {
          F77NAME(saxpy)(m-j-2,T.lowerDiagonalValue(j),M.addr(j+2,j+1),1,
            S->addr(j+2,j),1);
        }
      }
    }
  }
  return S;
}

template<> Matrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator*(
const UpperTrapezoidalMatrix<float,float> &M) const {
// compute by columns: note that
// [       S_11      , e_k sigma e_0^T ] [ u ]= [ S_11 u ]
// [ e_0 sigma e_k^T ,    S_22         ] [   ]= [ e_0 sigma e_k^T u ]
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<m-1;j++) {
      F77NAME(sstmv)(j+1,float_one_,L->addr(),D->addr(),
        M.addr(0,j),1,float_zero_,S->addr(0,j),1);
      (*S)(j+1,j)=lowerDiagonalValue(j)*M(j,j);
    }
    for (int j=m-1;j<n;j++) {
      F77NAME(sstmv)(m,float_one_,L->addr(),D->addr(),
        M.addr(0,j),1,float_zero_,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<m;j++) {
      if (j>0) {
        F77NAME(sstmv)(j,float_one_,L->addr(),D->addr(),
          M.addr(0,j),1,float_zero_,S->addr(0,j),1);
        (*S)(j-1,j)+=upperDiagonalValue(j-1);
        (*S)(j,j)=lowerDiagonalValue(j-1)*M(j-1,j);
      }
      (*S)(j,j)+=diagonalValue(j);
      if (j<m-1) (*S)(j+1,j)=lowerDiagonalValue(j);
    }
    for (int j=m;j<n;j++) {
      F77NAME(sstmv)(m,float_one_,L->addr(),D->addr(),
        M.addr(0,j),1,float_zero_,S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<float,float>* operator*(
const UpperTrapezoidalMatrix<float,float> &M,
const SymmetricTridiagonalMatrix<float,float> &T) {
// compute by columns: note that
// U T e_j = U e_{j-1} T_{j-1,j} + U e_j T_{j,j} + U e_{j+1} T_{j+1,j}
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool M_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&M)==0);
  if (M_non_unit) {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j,m),T.upperDiagonalValue(j-1),
          M.addr(0,j-1),1,S->addr(0,j),1);
      }
      F77NAME(saxpy)(min(j+1,m),T.diagonalValue(j),M.addr(0,j),1,
        S->addr(0,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(min(j+2,m),T.lowerDiagonalValue(j),
          M.addr(0,j+1),1,S->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(min(j-1,m),T.upperDiagonalValue(j-1),
          M.addr(0,j-1),1,S->addr(0,j),1);
        if (j<=m) (*S)(j-1,j)=T.upperDiagonalValue(j-1);
      }
      F77NAME(saxpy)(min(j,m),T.diagonalValue(j),M.addr(0,j),1,
        S->addr(0,j),1);
      if (j<m) (*S)(j,j)+=T.diagonalValue(j);
      if (j<n-1) {
        F77NAME(saxpy)(min(j+1,m),T.lowerDiagonalValue(j),
          M.addr(0,j+1),1,S->addr(0,j),1);
        if (j<m-1) (*S)(j+1,j)+=T.lowerDiagonalValue(j);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator*(
const SquareMatrix<float,float> &M) const {
  int n=M.size(0);
  CHECK_SAME(n,size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(sstmv)(n,float_one_,L->addr(),D->addr(),
      M.addr(0,j),1,float_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator*(
const SquareMatrix<float,float> &M,
const SymmetricTridiagonalMatrix<float,float> &T) {
  int n=M.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(saxpy)(n,T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(saxpy)(n,T.diagonalValue(j),M.addr(0,j),1,S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(n,T.lowerDiagonalValue(j),M.addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> Matrix<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(sstmv)(m,float_one_,L->addr(),D->addr(),
      M.addr(0,j),1,float_zero_,S->addr(0,j),1);
  }
  return S;
}

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,
const SymmetricTridiagonalMatrix<float,float> &T) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,T.size(0));
  Matrix<float,float> *S=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(saxpy)(m,T.upperDiagonalValue(j-1),M.addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(saxpy)(m,T.diagonalValue(j),M.addr(0,j),1,S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(m,T.lowerDiagonalValue(j),M.addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> Vector<float,float>*
SymmetricTridiagonalMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<float,float> *w=
    OPERATOR_NEW Vector<float,float>(n,float_zero_);
  F77NAME(sstmv)(n,float_one_,L->addr(),D->addr(),
    v.addr(),1,float_zero_,w->addr(),1);
  return w;
}

template<> float
SymmetricTridiagonalMatrix<float,float>::normFrobenius()
const {
  return F77NAME(slanst)('F',dim,D->addr(),L->addr());
}

template<> float
SymmetricTridiagonalMatrix<float,float>::normInfinity()
const {
  return F77NAME(slanst)('I',dim,D->addr(),L->addr());
}

template<> float
SymmetricTridiagonalMatrix<float,float>::normMaxEntry()
const {
  return F77NAME(slanst)('M',dim,D->addr(),L->addr());
}

template<> float
SymmetricTridiagonalMatrix<float,float>::normOne() const {
  return F77NAME(slanst)('O',dim,D->addr(),L->addr());
}

template<> float
SymmetricTridiagonalMatrix<float,float>::reciprocalConditionNumber(
char norm) const {
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  Vector<float,float> *UF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *UF2=
    OPERATOR_NEW Vector<float,float>(max(0,dim-2));
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*L);
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(sgttrf)(dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),ipiv,
    info);
  float rcond=undefined_;
  float *work=OPERATOR_NEW_BRACKET(float,2*dim);
  int *iwork=OPERATOR_NEW_BRACKET(int,dim);
  float anorm=F77NAME(slangt)(norm,dim,L->addr(),D->addr(),L->addr());
  F77NAME(sgtcon)(norm,dim,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
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

template<> void SymmetricTridiagonalMatrix<float,float>::stmv(
float alpha,const Vector<float,float> &x,float beta,
Vector<float,float> &b) const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  F77NAME(sstmv)(n,alpha,L->addr(),D->addr(),x.addr(),1,beta,b.addr(),1);
}

template<> void SymmetricTridiagonalMatrix<float,float>::stmm(
float alpha,const Matrix<float,float> &X,float beta,
Matrix<float,float> &B,char side) const {
  if (side=='L' || side=='l') {
    int m=size(0),n=X.size(1);
    CHECK_SAME(n,B.size(1));
    CHECK_SAME(m,X.size(0));
    CHECK_SAME(m,B.size(0));
    for (int j=0;j<n;j++) {
      F77NAME(sstmv)(m,alpha,L->addr(),D->addr(),X.addr(0,j),1,beta,
        B.addr(0,j),1);
    }
  } else {
    int m=X.size(0),n=X.size(1);
    CHECK_SAME(n,size(0));
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,B.size(1));
    if (abs(beta)==float_zero_) B=float_zero_;
    else if (beta!=float_one_) B*=beta;
    for (int j=0;j<n;j++) {
      if (j>0) {
        F77NAME(saxpy)(m,upperDiagonalValue(j-1)*alpha,X.addr(0,j-1),1,
          B.addr(0,j),1);
      }
      F77NAME(saxpy)(m,diagonalValue(j)*alpha,X.addr(0,j),1,
        B.addr(0,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(m,lowerDiagonalValue(j)*alpha,X.addr(0,j+1),1,
          B.addr(0,j),1);
      }
    }
  }
}

template<> Vector<float,float>* eigenvalues(
const SymmetricTridiagonalMatrix<float,float> &T,
OrthogonalMatrix<float,float> *&Q) {
  int dim=T.size(0);
  if (Q!=0) CHECK_SAME(dim,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'I');
  Vector<float,float> *lambda=OPERATOR_NEW Vector<float,float>(dim);
  lambda->copy(*T.diagonal());
  Vector<float,float> *L_copy=OPERATOR_NEW Vector<float,float>(dim-1);
  L_copy->copy(*T.lowerDiagonal());
  float *work=OPERATOR_NEW_BRACKET(float,2*dim-2);
  int info;
  float *qa=( Q==0 ? 0 : Q->addr() );
  F77NAME(ssteqr)(jobz,dim,lambda->addr(),L_copy->addr(),qa,dim,
    work,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete L_copy; L_copy=0;
  return lambda;
}

template<> void SymmetricTridiagonalMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  Vector<float,float> *UF=OPERATOR_NEW Vector<float,float>(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*L);
  x.copy(b);
  int info;
  F77NAME(sgtsv)(dim,1,LF->addr(),DF->addr(),UF->addr(),x.addr(),dim,
    info);
  CHECK_SAME(info,0);
  delete DF; DF=0;
  delete LF; LF=0;
  delete UF; UF=0;
}

template<> void SymmetricTridiagonalMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  Vector<float,float> *UF=OPERATOR_NEW Vector<float,float>(dim-1);
  LF->copy(*L);
  DF->copy(*D);
  UF->copy(*L);
  int info;
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    X.copy(B);
    F77NAME(sgtsv)(m,n,LF->addr(),DF->addr(),UF->addr(),X.addr(),m,
      info);
    CHECK_SAME(info,0);
  } else {
    CHECK_SAME(dim,n);
    float *t=OPERATOR_NEW_BRACKET(float,n);
    Vector<float,float> *UF2=OPERATOR_NEW Vector<float,float>(dim-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
    F77NAME(sgttrf)(n,LF->addr(),DF->addr(),UF->addr(),UF2->addr(),
      ipiv,info);
    CHECK_TEST(info==0);
    for (int i=0;i<m;i++) {
      F77NAME(scopy)(n,B.addr(i,0),m,t,1);
      F77NAME(sgttrs)('T',n,1,LF->addr(),DF->addr(),UF->addr(),
        UF2->addr(),ipiv,t,n,info);
      CHECK_TEST(info==0);
      F77NAME(scopy)(n,t,1,X.addr(i,0),m);
    }
    delete [] ipiv; ipiv=0;
    delete UF2; UF2=0;
    delete [] t; t=0;
  }
  delete LF; LF=0;
  delete DF; DF=0;
  delete UF; UF=0;
}
template class SymmetricTridiagonalMatrix<float,float>;
template void testSymmetricTridiagonalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveMatrix<float,float>*
SymmetricPositiveTridiagonalMatrix<float,float>::operator+(
const SymmetricPositiveMatrix<float,float> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SymmetricPositiveMatrix<float,float> *T=
    OPERATOR_NEW SymmetricPositiveMatrix<float,float>(n);
  T->copy(S);
  F77NAME(saxpy)(n,float_one_,D->addr(),1,T->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,L->addr(),1,T->addr(1,0),n+1);
  return T;
}

template<> float SymmetricPositiveTridiagonalMatrix<float,float>
::reciprocalConditionNumber(char norm) const {
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  LF->copy(*L);
  DF->copy(*D);
  int info;
  F77NAME(spttrf)(dim,DF->addr(),LF->addr(),info);
  CHECK_TEST(info==0);
  float rcond=undefined_;
  float *work=OPERATOR_NEW_BRACKET(float,2*dim);
  float anorm=normOne();
  F77NAME(sptcon)(dim,DF->addr(),LF->addr(),anorm,rcond,work,info);
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
  F77NAME(sstev)(jobz,dim,lambda->addr(),L_copy->addr(),qa,dim,work,info);
  delete [] work; work=0;
  delete L_copy; L_copy=0;
  return lambda;
}
*/

template<> void SymmetricPositiveTridiagonalMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  LF->copy(*L);
  DF->copy(*D);
  x.copy(b);
  int info;
  F77NAME(sptsv)(dim,1,DF->addr(),LF->addr(),x.addr(),dim,info);
  CHECK_SAME(info,0);
  delete DF; DF=0;
  delete LF; LF=0;
}

template<> void SymmetricPositiveTridiagonalMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side) const {
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(n,X.size(1));
  Vector<float,float> *LF=OPERATOR_NEW Vector<float,float>(dim-1);
  Vector<float,float> *DF=OPERATOR_NEW Vector<float,float>(dim);
  LF->copy(*L);
  DF->copy(*D);
  int info;
  if (side=='L' || side=='l') {
    CHECK_SAME(dim,m);
    X.copy(B);
    F77NAME(sptsv)(m,n,DF->addr(),LF->addr(),X.addr(),m,info);
    CHECK_SAME(info,0);
  } else {
    CHECK_SAME(dim,n);
    float *t=OPERATOR_NEW_BRACKET(float,n);
    F77NAME(spttrf)(n,DF->addr(),LF->addr(),info);
    CHECK_SAME(info,0);
    for (int i=0;i<m;i++) {
      F77NAME(scopy)(n,B.addr(i,0),m,t,1);
      F77NAME(spttrs)(n,1,DF->addr(),LF->addr(),t,n,info);
      CHECK_SAME(info,0);
      F77NAME(scopy)(n,t,1,X.addr(i,0),m);
    }
    delete [] t; t=0;
  }
  delete DF; DF=0;
  delete LF; LF=0;
}

template class SymmetricPositiveTridiagonalMatrix<float,float>;
template void testSymmetricPositiveTridiagonalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SquareMatrix<float,float>*
DiagonalMatrix<float,float>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,float> *M=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  F77NAME(scopy)(n,D->addr(),1,M->addr(0,0),n+1);
  return M;
}

template<> SymmetricTridiagonalMatrix<float,float>* operator+(
const DiagonalMatrix<float,float> &A,
const SymmetricTridiagonalMatrix<float,float> &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<float,float> *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<float,float>(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_one_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<float,float>*
DiagonalMatrix<float,float>::operator+(
const TridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,float> *S=OPERATOR_NEW
    TridiagonalMatrix<float,float>(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_one_,addr(),1,S->addr(0,0),1);
  return S;
}

template<> SymmetricMatrix<float,float>* operator+(
const DiagonalMatrix<float,float> &A,
const SymmetricMatrix<float,float> &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,float> *S=OPERATOR_NEW
    SymmetricMatrix<float,float>(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_one_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> UpperTriangularMatrix<float,float>*
DiagonalMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<float,float> *S=OPERATOR_NEW
    UpperTriangularMatrix<float,float>(n);
  S->copy(U);
  if (dynamic_cast<const UnitUpperTrapezoidalMatrix<float,float>*>(&U)
  ==0) {
    F77NAME(saxpy)(n,float_one_,addr(),1,S->addr(),n+1);
  } else {
    const float *di=D->addr();
    float *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di+float_one_;
  }
  return S;
}

template<> LowerTriangularMatrix<float,float>*
DiagonalMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<float,float> *S=OPERATOR_NEW
    LowerTriangularMatrix<float,float>(n);
  S->copy(L);
  if (dynamic_cast<const UnitLowerTrapezoidalMatrix<float,float>*>(&L)
  ==0) {
    F77NAME(saxpy)(n,float_one_,addr(),1,S->addr(),n+1);
  } else {
    const float *di=D->addr();
    float *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di+float_one_;
  }
  return S;
}

template<> SquareMatrix<float,float>*
DiagonalMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n);
  S->copy(M);
  F77NAME(saxpy)(n,float_one_,addr(),1,S->addr(),n+1);
  return S;
}

template<> SymmetricTridiagonalMatrix<float,float>* operator-(
const DiagonalMatrix<float,float> &A,
const SymmetricTridiagonalMatrix<float,float> &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<float,float> *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<float,float>(n);
  S->copy(T);
  *S*=float_mone_;
  F77NAME(saxpy)(n,float_one_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> SymmetricTridiagonalMatrix<float,float>* operator-(
const SymmetricTridiagonalMatrix<float,float> &T,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricTridiagonalMatrix<float,float> *S=OPERATOR_NEW
    SymmetricTridiagonalMatrix<float,float>(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_mone_,A.addr(),1,S->diagonalAddr(0),1);
  return S;
}

template<> TridiagonalMatrix<float,float>*
DiagonalMatrix<float,float>::operator-(
const TridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,float> *S=OPERATOR_NEW
    TridiagonalMatrix<float,float>(n);
  S->copy(T);
  *S*=float_mone_;
  F77NAME(saxpy)(n,float_one_,addr(),1,S->addr(0,0),1);
  return S;
}

template<> TridiagonalMatrix<float,float>* operator-(
const TridiagonalMatrix<float,float> &T,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,float> *S=OPERATOR_NEW
    TridiagonalMatrix<float,float>(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_mone_,A.addr(),1,S->addr(0,0),1);
  return S;
}

template<> SymmetricMatrix<float,float>* operator-(
const DiagonalMatrix<float,float> &A,
const SymmetricMatrix<float,float> &T) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,float> *S=OPERATOR_NEW
    SymmetricMatrix<float,float>(n);
  S->copy(T);
  *S*=float_mone_;
  F77NAME(saxpy)(n,float_one_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> SymmetricMatrix<float,float>* operator-(
const SymmetricMatrix<float,float> &T,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,float> *S=OPERATOR_NEW
    SymmetricMatrix<float,float>(n);
  S->copy(T);
  F77NAME(saxpy)(n,float_mone_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> UpperTriangularMatrix<float,float>*
DiagonalMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<float,float> *S=OPERATOR_NEW
    UpperTriangularMatrix<float,float>(n);
  S->copy(U);
  *S*=float_mone_;
  if (dynamic_cast<const UnitUpperTrapezoidalMatrix<float,float>*>(&U)
  ==0) {
    F77NAME(saxpy)(n,float_one_,addr(),1,S->addr(),n+1);
  } else {
    const float *di=D->addr();
    float *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di-float_one_;
  }
  return S;
}

template<> UpperTriangularMatrix<float,float>* operator-(
const UpperTrapezoidalMatrix<float,float> &U,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperTriangularMatrix<float,float> *S=OPERATOR_NEW
    UpperTriangularMatrix<float,float>(n);
  S->copy(U);
  if (dynamic_cast<const UnitUpperTrapezoidalMatrix<float,float>*>(&U)
  ==0) {
    F77NAME(saxpy)(n,float_mone_,A.addr(),1,S->addr(),n+1);
  } else {
    const float *di=A.addr();
    float *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=float_one_-*di;
  }
  return S;
}

template<> LowerTriangularMatrix<float,float>*
DiagonalMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<float,float> *S=OPERATOR_NEW
    LowerTriangularMatrix<float,float>(n);
  S->copy(L);
  *S*=float_mone_;
  if (dynamic_cast<const UnitLowerTrapezoidalMatrix<float,float>*>(&L)
  ==0) {
    F77NAME(saxpy)(n,float_one_,addr(),1,S->addr(),n+1);
  } else {
    const float *di=D->addr();
    float *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=*di-float_one_;
  }
  return S;
}

template<> LowerTriangularMatrix<float,float>* operator-(
const LowerTrapezoidalMatrix<float,float> &L,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  LowerTriangularMatrix<float,float> *S=OPERATOR_NEW
    LowerTriangularMatrix<float,float>(n);
  S->copy(L);
  if (dynamic_cast<const UnitLowerTrapezoidalMatrix<float,float>*>(&L)
  ==0) {
    F77NAME(saxpy)(n,float_mone_,A.addr(),1,S->addr(),n+1);
  } else {
    const float *di=A.addr();
    float *Sii=S->addr();
    for (int i=0;i<n;i++,di++,Sii+=n+1) *Sii=float_one_-*di;
  }
  return S;
}

template<> SquareMatrix<float,float>*
DiagonalMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n);
  S->copy(M);
  *S*=float_mone_;
  F77NAME(saxpy)(n,float_one_,addr(),1,S->addr(),n+1);
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const Matrix<float,float> &M,const DiagonalMatrix<float,float> &A)
{
  int n=A.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n);
  S->copy(M);
  F77NAME(saxpy)(n,float_mone_,A.addr(),1,S->addr(),n+1);
  return S;
}

template<> TridiagonalMatrix<float,float>* 
DiagonalMatrix<float,float>::operator*(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,float> *S=
    OPERATOR_NEW TridiagonalMatrix<float,float>(n);
  for (int i=0;i<n;i++) {
    float di=(*D)[i];
    if (i>0) (*S)(i,i-1)=di*T.lowerDiagonalValue(i-1);
    (*S)(i,i)=di*T.diagonalValue(i);
    if (i<n-1) (*S)(i,i+1)=di*T.lowerDiagonalValue(i);
  }
  return S;
}

template<> TridiagonalMatrix<float,float>* operator*(
const SymmetricTridiagonalMatrix<float,float> &T,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<float,float> *S=
    OPERATOR_NEW TridiagonalMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    float dj=A[j];
    if (j>0) (*S)(j-1,j)=T.lowerDiagonalValue(j-1)*dj;
    (*S)(j,j)=T.diagonalValue(j)*dj;
    if (j<n-1) (*S)(j+1,j)=T.lowerDiagonalValue(j)*dj;
  }
  return S;
}

template<> SquareMatrix<float,float>*
DiagonalMatrix<float,float>::operator*(
const SymmetricMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(dim);
  for (int i=0;i<n;i++) {
    float di=D->operator[](i);
    if (i>0) {
      const float *Tij=T.addr(i,0);
      float *Sij=S->addr(i,0);
      for (int j=0;j<i;j++,Tij+=n,Sij+=n) *Sij=di*(*Tij);
    }
    const float *Tji=T.addr(i,i);
    float *Sij=S->addr(i,i);
    for (int j=i;j<n;j++,Tji++,Sij+=n) *Sij=di*(*Tji);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator*(
const SymmetricMatrix<float,float> &T,
const DiagonalMatrix<float,float> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    float dj=A[j];
    if (j>0) {
      const float *Tji=T.addr(j,0);
      float *Sij=S->addr(0,j);
      for (int i=0;i<j;i++,Tji+=n,Sij++) *Sij=(*Tji)*dj;
    }
    const float *Tij=T.addr(j,j);
    float *Sij=S->addr(j,j);
    for (int i=j;i<n;i++,Tij++,Sij++) *Sij=(*Tij)*dj;
  }
  return S;
}

template class DiagonalMatrix<float,float>;
template void testDiagonalMatrix(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const float
  UpperHessenbergMatrix<float,float>::outofbounds_ = float_zero_;
template<> const float
  UpperHessenbergMatrix<float,float>::undefined_ = HUGE_VAL;
template<> float UpperHessenbergMatrix<float,float>::safety_ =
  float_zero_;

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::makeMatrix() const {
  int n=size(0);
  SquareMatrix<float,float> *M=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(n,j+2),addr(0,j),1,M->addr(0,j),1);
  }
  return M;
}

template<> UpperHessenbergMatrix<float,float>&
UpperHessenbergMatrix<float,float>::operator+=(
const UpperHessenbergMatrix<float,float> &H) {
  int n=size(0);
  CHECK_SAME(n,H.size(0));
  float *colj=addr();
  const float *Hcolj=H.addr();
  for (int j=0;j<n;j++,colj+=n,Hcolj+=n) {
    F77NAME(saxpy)(min(n,j+2),float_one_,Hcolj,1,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<float,float>&
UpperHessenbergMatrix<float,float>::operator-=(
const UpperHessenbergMatrix<float,float> &H) {
  int n=size(0);
  CHECK_SAME(n,H.size(0));
  float *colj=addr();
  const float *Hcolj=H.addr();
  for (int j=0;j<n;j++,colj+=n,Hcolj+=n) {
    F77NAME(saxpy)(min(n,j+2),float_mone_,Hcolj,1,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<float,float>&
UpperHessenbergMatrix<float,float>::operator*=(float scalar) {
  int n=size(0);
  float *colj=addr();
  for (int j=0;j<n;j++,colj+=n) {
    F77NAME(sscal)(min(n,j+2),scalar,colj,1);
  }
  return *this;
}

template<> UpperHessenbergMatrix<float,float>&
UpperHessenbergMatrix<float,float>::operator/=(float scalar) {
  int n=size(0);
  float *colj=addr();
  float s=float_one_/scalar;
  for (int j=0;j<n;j++,colj+=n) {
    F77NAME(sscal)(min(n,j+2),s,colj,1);
  }
  return *this;
}

template<> void UpperHessenbergMatrix<float,float>::copy(
const Matrix<float,float> &M) {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(j+2,n),M.addr(0,j),1,addr(0,j),1);
  }
}

/*
template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::transpose() const {
  int n=size(0);
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(j+2,n),addr(0,j),1,S->addr(j,0),n);
  }
  return S;
}
*/

template<> void UpperHessenbergMatrix<float,float>::copyFrom(
int m,const SquareMatrix<float,float> &S) {
  m=min(m,size(0));
  m=min(m,S.size(0));
  for (int j=0;j<m;j++) {
    F77NAME(scopy)(min(j+2,m),S.addr(0,j),1,addr(0,j),1);
  }
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator+(
const DiagonalMatrix<float,float> &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(*this);
  F77NAME(saxpy)(n,float_one_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator+(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(*this);
  F77NAME(saxpy)(n,float_one_,T.diagonalAddr(0),1,H->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,T.lowerDiagonalAddr(0),1,
    H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator+(
const TridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(*this);
  F77NAME(saxpy)(n,float_one_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator+(
const SymmetricMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    if (j<n-1) F77NAME(scopy)(n-j-1,M.addr(j+1,j),1,S->addr(j,j+1),n);
    F77NAME(saxpy)(min(j+2,n),float_one_,addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      F77NAME(saxpy)(n-j,float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)+=float_one_;
      if (j<n-1) {
        F77NAME(saxpy)(n-j-1,float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(*this);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(j+1,float_one_,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)+=float_one_;
      if (j>0) {
        F77NAME(saxpy)(j,float_one_,U.addr(0,j),1,H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(j+2,n),float_one_,addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator-(
const DiagonalMatrix<float,float> &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(*this);
  F77NAME(saxpy)(n,float_mone_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,float>* operator-(
const DiagonalMatrix<float,float> &D,
const UpperHessenbergMatrix<float,float> &U) {
  int n=U.size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(U);
  *H*=float_mone_;
  F77NAME(saxpy)(n,float_one_,D.addr(),1,H->addr(0,0),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator-(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(*this);
  F77NAME(saxpy)(n,float_mone_,T.diagonalAddr(0),1,H->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.lowerDiagonalAddr(0),1,
    H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,float>* operator-(
const SymmetricTridiagonalMatrix<float,float> &T,
const UpperHessenbergMatrix<float,float> &U) {
  int n=U.size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(U);
  *H*=float_mone_;
  F77NAME(saxpy)(n,float_one_,T.diagonalAddr(0),1,H->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,T.lowerDiagonalAddr(0),1,
    H->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,T.lowerDiagonalAddr(0),1,
    H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator-(
const TridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(*this);
  F77NAME(saxpy)(n,float_mone_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_mone_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> UpperHessenbergMatrix<float,float>* operator-(
const TridiagonalMatrix<float,float> &T,
const UpperHessenbergMatrix<float,float> &U) {
  int n=U.size(0);
  CHECK_SAME(n,T.size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(U);
  *H*=float_mone_;
  F77NAME(saxpy)(n,float_one_,T.addr(0,0),1,H->addr(0,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,T.addr(1,0),1,H->addr(1,0),n+1);
  F77NAME(saxpy)(n-1,float_one_,T.addr(0,1),1,H->addr(0,1),n+1);
  return H;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator-(
const SymmetricMatrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(n-j,float_mone_,M.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      F77NAME(saxpy)(n-j-1,float_mone_,M.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const SymmetricMatrix<float,float> &M,
const UpperHessenbergMatrix<float,float> &H) {
  int n=H.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,M.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      F77NAME(scopy)(n-j-1,M.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(j+2,n),float_mone_,H.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int n=size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      F77NAME(saxpy)(n-j,float_mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)-=float_one_;
      if (j<n-1) {
        F77NAME(saxpy)(n-j-1,float_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const LowerTrapezoidalMatrix<float,float> &L,
const UpperHessenbergMatrix<float,float> &H) {
  int n=H.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      F77NAME(saxpy)(min(j+2,n),float_mone_,H.addr(0,j),1,
        S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=float_one_;
      if (j<n-1) {
        F77NAME(scopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      F77NAME(saxpy)(min(j+2,n),float_mone_,H.addr(0,j),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &U) const {
  int n=size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(*this);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(j+1,float_mone_,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)-=float_one_;
      if (j>0) {
        F77NAME(saxpy)(j,float_mone_,U.addr(0,j),1,H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> UpperHessenbergMatrix<float,float>* operator-(
const UpperTrapezoidalMatrix<float,float> &U,
const UpperHessenbergMatrix<float,float> &M) {
  int n=M.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n,float_zero_);
  H->copy(M);
  *H*=float_mone_;
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  if (U_non_unit) {
    for (int j=0;j<n;j++) {
      F77NAME(saxpy)(j+1,float_one_,U.addr(0,j),1,H->addr(0,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*H)(j,j)+=float_one_;
      if (j>0) {
        F77NAME(saxpy)(j,float_one_,U.addr(0,j),1,H->addr(0,j),1);
      }
    }
  }
  return H;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(j+2,n),addr(0,j),1,S->addr(0,j),1);
    F77NAME(saxpy)(n,float_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const Matrix<float,float> &M,
const UpperHessenbergMatrix<float,float> &H) {
  int n=H.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n,M.addr(0,j),1,S->addr(0,j),1);
    F77NAME(saxpy)(min(j+2,n),float_mone_,H.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const UpperHessenbergMatrix<float,float> &H2) const {
  int n=size(0);
  CHECK_SAME(n,H2.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<=min(j+1,n-1);k++) {
      F77NAME(saxpy)(min(k+2,n),H2(k,j),addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> UpperHessenbergMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const DiagonalMatrix<float,float> &D) const {
  int n=size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,float> *S=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n);
  S->copy(*this);
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(min(j+2,n),D[j],S->addr(0,j),1);
  }
  return S;
}

template<> UpperHessenbergMatrix<float,float>* operator*(
const DiagonalMatrix<float,float> &D,
const UpperHessenbergMatrix<float,float> &H) {
  int n=H.size(0);
  CHECK_SAME(n,D.size(0));
  UpperHessenbergMatrix<float,float> *S=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(n);
  S->copy(H);
  for (int i=0;i<n;i++) {
    int j=max(0,i-1);
    F77NAME(sscal)(n-j,D[i],S->addr(i,j),n);
  }
  return S;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const SymmetricTridiagonalMatrix<float,float> &T) const {
//compute by columns: note that
// H T e_j = H e_{j-1} T_{j-1,j} + H e_j T_{j,j} + H e_{j+1} T_{j+1,j}
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    if (j>0) {
      F77NAME(saxpy)(min(j+1,n),T.upperDiagonalValue(j-1),addr(0,j-1),1,
        S->addr(0,j),1);
    }
    F77NAME(saxpy)(min(j+2,n),T.diagonalValue(j),addr(0,j),1,
      S->addr(0,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(min(j+3,n),T.lowerDiagonalValue(j),addr(0,j+1),1,
        S->addr(0,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator*(
const SymmetricTridiagonalMatrix<float,float> &T,
const UpperHessenbergMatrix<float,float> &H) {
// compute by rows: note that
// e_i^T T H
//   = T_{i,i-1} e_{i-1}^T H + T_{i,i} e_i^T H + T_{i,i+1} e_{i+1}^T H
  int n=H.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int i=0;i<n;i++) {
    if (i>0) {
      int j=max(0,i-2);
      F77NAME(saxpy)(n-j,T.lowerDiagonalValue(i-1),H.addr(i-1,j),n,
        S->addr(i,j),n);
    }
    int j=max(0,i-1);
    F77NAME(saxpy)(n-j,T.diagonalValue(i),H.addr(i,j),n,S->addr(i,j),n);
    if (i<n-1) {
      F77NAME(saxpy)(n-i,T.upperDiagonalValue(i),H.addr(i+1,i),n,
        S->addr(i,i),n);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const TridiagonalMatrix<float,float> &T) const {
//compute by columns: note that
// H T e_j = H e_{j-1} T_{j-1,j} + H e_j T_{j,j} + H e_{j+1} T_{j+1,j}
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(j-1,0);k<=min(j+1,n-1);k++) {
      F77NAME(saxpy)(min(k+2,n),T(k,j),addr(0,k),1,S->addr(0,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator*(
const TridiagonalMatrix<float,float> &T,
const UpperHessenbergMatrix<float,float> &H) {
// compute by rows: note that
// e_i^T T H
//   = T_{i,i-1} e_{i-1}^T H + T_{i,i} e_i^T H + T_{i,i+1} e_{i+1}^T H
  int n=H.size(0);
  CHECK_SAME(n,T.size(0));
  SquareMatrix<float,float> *S=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int i=0;i<n;i++) {
    for (int k=max(i-1,0);k<=min(i+1,n-1);k++) {
      int j=max(0,k-1);
      F77NAME(saxpy)(min(n-k+1,n),T(i,k),H.addr(k,j),n,S->addr(i,j),n);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const SymmetricMatrix<float,float> &S) const {
// compute by bordering: note that
// [     eta_11 h^T ] [ sigma s^T ]
// [ e_0 eta_21  H  ] [   s    S  ]
//   = [     eta_11 sigma + h^T s ,     eta_11 s^T + h^T S ]
//   = [ e_0 eta_21 sigma +  H  s , e_0 eta_21 s^T +  H  S ]
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *M=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  (*M)(n-1,n-1)=(*this)(n-1,n-1)*S(n-1,n-1);
  for (int k=n-2;k>=0;k--) {
    for (int j=k+1;j<n;j++) { // H s
      F77NAME(saxpy)(min(n-k-1,j-k+1),S(j,k),addr(k+1,j),1,
        M->addr(k+1,k),1);
    }
    F77NAME(ssymv)('L',n-k-1,float_one_,S.addr(k+1,k+1),n,
      addr(k,k+1),n,float_zero_,M->addr(k,k+1),n); // h^T S = ( S h )^T
    (*M)(k,k)=F77NAME(sdot)(n-k-1,addr(k,k+1),n,S.addr(k+1,k),1); // h^T s
    F77NAME(sger)(2,n-k,float_one_,addr(k,k),1,S.addr(k,k),1,
      M->addr(k,k),n);
  }
  return M;
}

template<> SquareMatrix<float,float>* operator*(
const SymmetricMatrix<float,float> &S,
const UpperHessenbergMatrix<float,float> &H) {
// compute by bordering: note that
// [ sigma s^T ] [     eta_11 h^T ]
// [   s    S  ] [ e_0 eta_21  H  ]
//   = [ sigma eta_11 + s^T e_0 eta_21 , sigma h^T + s^T H ]
//   = [   s   eta_11 +  S  e_0 eta_21 ,   s   h^T +  S  H ]
  int n=H.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *M=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  (*M)(n-1,n-1)=S(n-1,n-1)*H(n-1,n-1);
  for (int k=n-2;k>=0;k--) { // s^T H
    F77NAME(saxpy)(n-k-1,H(k+1,k),S.addr(k+1,k+1),1,M->addr(k+1,k),1);
      // S e_0 eta_21
    for (int j=k+1;j<n;j++) { // s^T H
      (*M)(k,j)=
        F77NAME(sdot)(min(n-k-1,j-k+1),S.addr(k+1,k),1,H.addr(k+1,j),1);
    }
    (*M)(k,k)=S(k+1,k)*H(k+1,k); // s^T e_0 eta_21
    F77NAME(sger)(n-k,n-k,float_one_,S.addr(k,k),1,H.addr(k,k),n,
      M->addr(k,k),n);
  }
  return M;
}

template<> Matrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &L) const {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0);
  if (L_non_unit) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<m;k++) {
        F77NAME(saxpy)(min(k+2,m),L(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(min(j+2,m),addr(0,j),1,M->addr(0,j),1);
      for (int k=j+1;k<m;k++) {
        F77NAME(saxpy)(min(k+2,m),L(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<float,float>* operator*(
const LowerTrapezoidalMatrix<float,float> &L,
const UpperHessenbergMatrix<float,float> &H) {
// compute by columns: note that
//   [ L_11      ] [ h ] = [ L_11 h ]
//   [ L_21 L_22 ] [   ] = [ L_21 h ]
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,H.size(0));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  char diag=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0 ? 'N' : 'U');
  for (int j=0;j<n;j++) {
    int k=min(j+2,n);
    F77NAME(scopy)(k,H.addr(0,j),1,M->addr(0,j),1);
    F77NAME(strmv)('L','N',diag,k,L.addr(),m,M->addr(0,j),1); // L_11 h
    if (k<m) {
      F77NAME(sgemv)('N',m-k,k,float_one_,L.addr(k,0),m,H.addr(0,j),1,
        float_zero_,M->addr(k,j),1); // L_21 h
    }
  }
  return M;
}

template<> Matrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const UpperTrapezoidalMatrix<float,float> &U) const {
// compute by columns: note that
//   H [ U_1 , U_2 ] = [ H U_1 , H U_2 ]
// and that
//   [      H_11     H_12 ] [ u ] = [     H_11      u ]
//   [ e_0 eta e_k^T H_22 ] [   ] = [ e_0 eta e_k^T u ]
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0);
  if (U_non_unit) { // H U_1
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,m-1);k++) {
        F77NAME(saxpy)(min(k+2,m),U(k,j),addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,m);k++) {
        F77NAME(saxpy)(min(k+2,m),U(k,j),addr(0,k),1,M->addr(0,j),1);
      }
      if (j<m) {
        F77NAME(saxpy)(min(j+2,m),float_one_,addr(0,j),1,M->addr(0,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<float,float>* operator*(
const UpperTrapezoidalMatrix<float,float> &U,
const UpperHessenbergMatrix<float,float> &H) {
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
    dynamic_cast<const UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0
    ? 'N' : 'U');
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<m;j++) { // U_1 H_11
    int k=min(j+2,m);
    F77NAME(scopy)(k,H.addr(0,j),1,M->addr(0,j),1);
    F77NAME(strmv)('U','N',diag,k,U.addr(),m,M->addr(0,j),1); // U_11 h
  }
  if (m<n) {
    F77NAME(saxpy)(m,H(m,m-1),U.addr(0,m),1,M->addr(0,m-1),1);
      // U_2 e_0 eta e_{m-1}^T
    F77NAME(slacpy)('A',m,n-m,H.addr(0,m),m,M->addr(0,m),m);
    F77NAME(strmm)('L','U','N',diag,m,n-m,float_one_,U.addr(),m,
      M->addr(0,m),m); // U_1 H_12
    F77NAME(sgemm)('N','N',m,n-m,n-m,float_one_,U.addr(0,m),m,
      H.addr(m,m),n,float_one_,M->addr(0,m),m); // U_2 H_22
  }
  return M;
}

template<> SquareMatrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const SquareMatrix<float,float> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *M=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=0;k<n;k++) {
      F77NAME(saxpy)(min(k+2,n),S(k,j),addr(0,k),1,M->addr(0,j),1);
    }
  }
  return M;
}

template<> SquareMatrix<float,float>* operator*(
const SquareMatrix<float,float> &S,
const UpperHessenbergMatrix<float,float> &H) {
  int n=H.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *M=
    OPERATOR_NEW SquareMatrix<float,float>(n);
  for (int j=0;j<n;j++) {
    F77NAME(sgemv)('N',n,min(j+2,n),float_one_,S.addr(),n,H.addr(0,j),1,
      float_zero_,M->addr(0,j),1);
  }
  return M;
}

template<> Matrix<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,size(0));
  Matrix<float,float> *R=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(j-1,0);k<m;k++) {
      F77NAME(saxpy)(min(k+2,m),M(k,j),addr(0,k),1,R->addr(0,j),1);
    }
  }
  return R;
}

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,
const UpperHessenbergMatrix<float,float> &H) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,H.size(0));
  Matrix<float,float> *R=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(sgemv)('N',m,min(j+2,n),float_one_,M.addr(),m,H.addr(0,j),1,
      float_zero_,R->addr(0,j),1);
  }
  return R;
}

template<> Vector<float,float>*
UpperHessenbergMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  int n=size(0);
  CHECK_SAME(n,v.size());
  Vector<float,float> *r=
    OPERATOR_NEW Vector<float,float>(n,float_zero_);
  for (int k=0;k<n;k++) {
    F77NAME(saxpy)(min(k+2,n),v[k],addr(0,k),1,r->addr(),1);
  }
  return r;
}

template<> void UpperHessenbergMatrix<float,float>::uhmv(float alpha,
const Vector<float,float> &x,float beta,Vector<float,float> &b,
char trans) const {
  int n=size(0);
  CHECK_SAME(n,x.size());
  CHECK_SAME(n,b.size());
  if (abs(beta)==float_zero_) b=float_zero_;
  else b*=beta;
  if (trans=='N' || trans=='n') { // b=H*x*alpha+b*beta
    float *xj=x.addr();
    for (int j=0;j<n;j++,xj++) {
      F77NAME(saxpy)(min(n,j+2),(*xj)*alpha,addr(0,j),1,b.addr(),1); 
    }
  } else { // b=H^T*x*alpha+b*beta
    float *bj=b.addr();
    for (int j=0;j<n;j++,bj++) {
      *bj+=F77NAME(sdot)(min(n,j+2),addr(0,j),1,x.addr(),1)*alpha;
    }
  }
}

template<> void UpperHessenbergMatrix<float,float>::uhmm(float alpha,
const Matrix<float,float> &X,float beta,Matrix<float,float> &B,
char side,char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (abs(beta)==float_zero_) B=float_zero_;
  else B*=beta;
  if (trans=='N' || trans=='n') {
    if (side=='L' || side=='l') { // B=H*X*alpha+B*beta
      CHECK_SAME(m,size(0));
      for (int j=0;j<n;j++) {
        for (int k=0;k<m;k++) {
          F77NAME(saxpy)(min(m,k+2),X(k,j)*alpha,addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else { // B=X*H*alpha+B*beta
      CHECK_SAME(n,size(0));
      for (int j=0;j<n;j++) {
        for (int k=0;k<=min(n-1,j+1);k++) {
          F77NAME(saxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  } else {
    if (side=='L' || side=='l') { // B=H^T*X*alpha+B*beta
      CHECK_SAME(m,size(0));
      for (int j=0;j<n;j++) {
        for (int i=0;i<m;i++) {
          B(i,j)+=F77NAME(sdot)(min(m,i+2),addr(0,i),1,X.addr(0,j),1)
                *alpha;
        }
      }
    } else { // B=X*H^T*alpha+B*beta
      CHECK_SAME(n,size(0));
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-1);k<n;k++) {
          F77NAME(saxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

template<> float UpperHessenbergMatrix<float,float>::normFrobenius()
const {
  int n=size(0);
  float *work=0;
  return F77NAME(slanhs)('F',n,addr(),n,work);
}

template<> float UpperHessenbergMatrix<float,float>::normInfinity()
const {
  int n=size(0);
  float *work=OPERATOR_NEW_BRACKET(float,n);
  float val=F77NAME(slanhs)('I',n,addr(),n,work);
  delete [] work; work=0;
  return val;
}

template<> float UpperHessenbergMatrix<float,float>::normMaxEntry()
const {
  int n=size(0);
  float *work=0;
  return F77NAME(slanhs)('M',n,addr(),n,work);
}

template<> float UpperHessenbergMatrix<float,float>::normOne()
const {
  int n=size(0);
  float *work=0;
  return F77NAME(slanhs)('O',n,addr(),n,work);
}

template<> Vector<float,complex<float> >*
UpperHessenbergMatrix<float,float>::eigenvalues(
SquareMatrix<float,complex<float> > *&V,
SquareMatrix<float,complex<float> > *&U) const {
  int n=size(0);
  if (V!=0) CHECK_SAME(n,V->size(0));
  if (U!=0) CHECK_SAME(n,U->size(0));
  UpperHessenbergMatrix<float,float> *H=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(*this);
  char job=(V==0 && U==0 ? 'E' : 'S');
  char compz=(V==0 && U==0 ? 'N' : 'I');
  float *wr=OPERATOR_NEW_BRACKET(float,n);
  float *wi=OPERATOR_NEW_BRACKET(float,n);
  OrthogonalMatrix<float,float> *Z=(V==0 && U==0 ? 0 :
    OPERATOR_NEW OrthogonalMatrix<float,float>(n,n));
  float *za=(Z==0 ? 0 : Z->addr());
  float w=float_undefined_;
  int lwork=-1;
  int info;
  F77NAME(shseqr)(job,compz,n,1,n,H->addr(),n,wr,wi,za,n,&w,lwork,info);
  CHECK_TEST(info==0);

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(shseqr)(job,compz,n,1,n,H->addr(),n,wr,wi,za,n,work,lwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;

  Vector<float,complex<float> > *lambda =
    OPERATOR_NEW Vector<float,complex<float> >(n);
  for (int j=0;j<n;j++) {
    complex<float> &ev=(*lambda)[j];
    ev.real()=wr[j];
    ev.imag()=wi[j];
  }

  if (Z!=0) {
    char side=(V==0 ? 'R' : (U==0 ? 'L' : 'B') );
    bool *select=0;
    float *vla=0;
    SquareMatrix<float,float> *Vl=0;
    if (V!=0) {
      Vl=OPERATOR_NEW SquareMatrix<float,float>(n);
      Vl->copy(*Z);
      vla=Vl->addr();
    }
    float *vra=0;
    SquareMatrix<float,float> *Vr=0;
    if (U!=0) {
      Vr=OPERATOR_NEW SquareMatrix<float,float>(n);
      Vr->copy(*Z);
      vra=Vr->addr();
    }
    work=OPERATOR_NEW_BRACKET(float,3*n);
    int mout=-1;
    F77NAME(strevc)(side,'B',select,n,H->addr(),n,vla,n,vra,n,n,mout,work,
      info);
    CHECK_TEST(info==0);
    for (int j=0;j<n;) {
      int jn=j*n;
      if (wi[j]>zero_) {
        int jp1n=jn+n;
        if (V!=0) {
          for (int i=0;i<n;i++) {
            complex<float> &z=V->operator()(i,j);
            z.real()=vla[i+jn];
            z.imag()=vla[i+jp1n];
            complex<float> &zz=V->operator()(i,j+1);
            zz.real()=vla[i+jn];
            zz.imag()=-vla[i+jp1n];
          }
        }
        if (U!=0) {
          for (int i=0;i<n;i++) {
            complex<float> &z=U->operator()(i,j);
            z.real()=vra[i+jn];
            z.imag()=vra[i+jp1n];
            complex<float> &zz=U->operator()(i,j+1);
            zz.real()=vra[i+jn];
            zz.imag()=-vra[i+jp1n];
          }
        }
        j+=2;
      } else {
        if (V!=0) {
          for (int i=0;i<n;i++) {
            complex<float> &z=V->operator()(i,j);
            z.real()=vla[i+jn];
            z.imag()=zero_;
          }
        }
        if (U!=0) {
          for (int i=0;i<n;i++) {
            complex<float> &z=U->operator()(i,j);
            z.real()=vra[i+jn];
            z.imag()=zero_;
          }
        }
        j++;
      }
    }
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

template<> int UpperHessenbergMatrix<float,float>::factor(int *ipiv) {
  int n=size(0);
  for (int k=0;k<n;k++) ipiv[k]=k;
  int info=0;
  for (int k=0;k<n-1;k++) {
    if (abs((*this)(k,k))<abs((*this)(k+1,k)) ) {
      ipiv[k]=k+1;
      F77NAME(sswap)(n-k,addr(k,k),n,addr(k+1,k),n);
    }
    if (abs((*this)(k,k))>zero_) {
      (*this)(k+1,k)/=(*this)(k,k);
      F77NAME(saxpy)(n-k-1,-(*this)(k+1,k),addr(k,k+1),n,
        addr(k+1,k+1),n);
    } else {
      info=k+1;
      break;
    }
  }
  return info;
}

template<> void UpperHessenbergMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  int n=size(0);
  CHECK_SAME(n,b.size());
  CHECK_SAME(n,x.size());
  x.copy(b);
  UpperHessenbergMatrix<float,float> *HF=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(*this);
  int *ipiv=OPERATOR_NEW_BRACKET(int,n);
  int info=HF->factor(ipiv);
  CHECK_TEST(info==0);
  if (trans!='N' && trans!='n') {
    x[0]/=(*HF)(0,0);
    for (int j=1;j<n;j++) {
      x[j]=(x[j]-F77NAME(sdot)(j,HF->addr(0,j),1,x.addr(),1))
          /(*HF)(j,j);
    }
    for (int i=n-2;i>=0;i--) {
      int ip=ipiv[i];
      float temp=x[i]-(*HF)(i+1,i)*x[i+1];
      x[i]=x[ip];
      x[ip]=temp;
    }
  } else {
    for (int i=0;i<n-1;i++) {
      int ip=ipiv[i];
      float temp=x[2*i-ip+1]-(*HF)(i+1,i)*x[ip];
      x[i]=x[ip];
      x[i+1]=temp;
    }
    for (int j=n-1;j>0;j--) {
      x[j]/=(*HF)(j,j);
      F77NAME(saxpy)(j,-x[j],HF->addr(0,j),1,x.addr(),1);
    }
    x[0]/=(*HF)(0,0);
  }
  delete ipiv; ipiv=0;
  delete HF; HF=0;
}

template<> void UpperHessenbergMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  int n=size(0);
  UpperHessenbergMatrix<float,float> *HF=
    OPERATOR_NEW UpperHessenbergMatrix<float,float>(*this);
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
          X(j,k)=(X(j,k)-F77NAME(sdot)(j,HF->addr(0,j),1,X.addr(0,k),1))
              /(*HF)(j,j);
        }
        for (int i=n-2;i>=0;i--) {
          int ip=ipiv[i];
          float temp=X(i,k)-(*HF)(i+1,i)*X(i+1,k);
          X(i,k)=X(ip,k);
          X(ip,k)=temp;
        }
      }
    } else {
      for (int k=0;k<nrhs;k++) {
        for (int i=0;i<n-1;i++) {
          int ip=ipiv[i];
          float temp=X(2*i-ip+1,k)-(*HF)(i+1,i)*X(ip,k);
          X(i,k)=X(ip,k);
          X(i+1,k)=temp;
        }
        for (int j=n-1;j>0;j--) {
          X(j,k)/=(*HF)(j,j);
          F77NAME(saxpy)(j,-X(j,k),HF->addr(0,j),1,X.addr(0,k),1);
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
          float temp=X(k,2*i-ip+1)-(*HF)(i+1,i)*X(k,ip);
          X(k,i)=X(k,ip);
          X(k,i+1)=temp;
        }
        for (int j=n-1;j>0;j--) {
          X(k,j)/=(*HF)(j,j);
          F77NAME(saxpy)(j,-X(k,j),HF->addr(0,j),1,X.addr(k,0),nrhs);
        }
        X(k,0)/=(*HF)(0,0);
      }
    } else {
      for (int k=0;k<nrhs;k++) {
        X(k,0)/=(*HF)(0,0);
        for (int j=1;j<n;j++) {
          X(k,j)=(X(k,j)
            -F77NAME(sdot)(j,HF->addr(0,j),1,X.addr(k,0),nrhs))
            /(*HF)(j,j);
        }
        for (int i=n-2;i>=0;i--) {
          int ip=ipiv[i];
          float temp=X(k,i)-(*HF)(i+1,i)*X(k,i+1);
          X(k,i)=X(k,ip);
          X(k,ip)=temp;
        }
      }
    }
  }
  delete ipiv; ipiv=0;
  delete HF; HF=0;
}

template class UpperHessenbergMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const float BandMatrix<float,float>::outofbounds_ =
  float_zero_;
template<> float BandMatrix<float,float>::safety_ = float_zero_;

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::makeMatrix() const {
  SquareMatrix<float,float> *M=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      M->addr(ibeg,j),1);
  }
  return M;
}

template<> BandMatrix<float,float>&
BandMatrix<float,float>::operator+=(const BandMatrix<float,float> &B)
{
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  CHECK_SAME(nsup,B.nsup);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,float_one_,
      B.addr(ibeg,j),1,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<float,float>&
BandMatrix<float,float>::operator-=(const BandMatrix<float,float> &B)
{
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  CHECK_SAME(nsup,B.nsup);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<float,float>&
BandMatrix<float,float>::operator*=(float d) {
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(sscal)(min(dim-1,j+nsub)-ibeg+1,d,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<float,float>&
BandMatrix<float,float>::operator/=(float d) {
  CHECK_TEST(abs(d)>float_zero_);
  float dinv=float_one_/d;
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(sscal)(min(dim-1,j+nsub)-ibeg+1,dinv,addr(ibeg,j),1);
  }
  return *this;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator+(const BandMatrix<float,float> &B)
const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,B.nsub),max(nsup,B.nsup),
    float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    ibeg=max(0,j-B.nsup);
    F77NAME(saxpy)(min(dim-1,j+B.nsub)-ibeg+1,float_one_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator+(
const UpperHessenbergMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    F77NAME(saxpy)(min(dim,j+2),float_one_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator+(
const DiagonalMatrix<float,float> &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<float,float> *S=
    OPERATOR_NEW BandMatrix<float,float>(*this);
  F77NAME(saxpy)(dim,float_one_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator+(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,1),max(nsup,1),float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(saxpy)(dim-1,float_one_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  F77NAME(saxpy)(dim,float_one_,T.diagonalAddr(0),1,
    S->addr(0,0),S->nt);
  F77NAME(saxpy)(dim-1,float_one_,T.lowerDiagonalAddr(0),1,
    S->addr(0,1),S->nt);
  return S;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator+(
const TridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,1),max(nsup,1),float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(saxpy)(dim,float_one_,T.diagonalAddr(),1,
    S->addr(0,0),S->nt);
  F77NAME(saxpy)(dim-1,float_one_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S->nt);
  F77NAME(saxpy)(dim-1,float_one_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S->nt);
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator+(
const SymmetricMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<dim;j++) {
    F77NAME(saxpy)(dim-j,float_one_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(saxpy)(dim-j-1,float_one_,H.addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(saxpy)(j+1,float_one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      if (j>0) {
        F77NAME(saxpy)(j,float_one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(saxpy)(dim-j,float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      (*S)(j,j)+=float_one_;
      if (j+1<dim) {
        F77NAME(saxpy)(dim-j-1,float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  S->copy(M);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,float_one_,
      addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator-(const BandMatrix<float,float> &B)
const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,B.nsub),max(nsup,B.nsup),
    float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    ibeg=max(0,j-B.nsup);
    F77NAME(saxpy)(min(B.dim-1,j+B.nsub)-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator-(
const UpperHessenbergMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
    F77NAME(saxpy)(min(dim,j+2),float_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const UpperHessenbergMatrix<float,float> &H,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(n,j+2),H.addr(0,j),1,S->addr(0,j),1);
    int ibeg=max(0,j-B.supDiags());
    F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator-(
const DiagonalMatrix<float,float> &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<float,float> *S=
    OPERATOR_NEW BandMatrix<float,float>(*this);
  F77NAME(saxpy)(dim,float_mone_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> BandMatrix<float,float>* operator-(
const DiagonalMatrix<float,float> &D,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(n,B.subDiags(),B.supDiags(),float_zero_);
  F77NAME(scopy)(n,D.addr(),1,S->addr(0,0),S->bands());
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator-(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,1),max(nsup,1),float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(saxpy)(dim-1,float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  F77NAME(saxpy)(dim,float_mone_,T.diagonalAddr(0),1,
    S->addr(0,0),S->nt);
  F77NAME(saxpy)(dim-1,float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(0,1),S->nt);
  return S;
}

template<> BandMatrix<float,float>* operator-(
const SymmetricTridiagonalMatrix<float,float> &T,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(n,max(B.subDiags(),1),max(B.supDiags(),1),
    float_zero_);
  int nb=S->bands();
  F77NAME(scopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(1,0),nb);
  F77NAME(scopy)(n,T.diagonalAddr(0),1,S->addr(0,0),nb);
  F77NAME(scopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(0,1),nb);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator-(
const TridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,1),max(nsup,1),float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  F77NAME(saxpy)(dim,float_mone_,T.diagonalAddr(),1,
    S->addr(0,0),nt);
  F77NAME(saxpy)(dim-1,float_mone_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),nt);
  F77NAME(saxpy)(dim-1,float_mone_,T.upperDiagonalAddr(),1,
    S->addr(0,1),nt);
  return S;
}

template<> BandMatrix<float,float>* operator-(
const TridiagonalMatrix<float,float> &T,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(n,max(B.subDiags(),1),max(B.supDiags(),1),
    float_zero_);
  F77NAME(scopy)(n,T.diagonalAddr(),1,S->addr(0,0),S->bands());
  F77NAME(scopy)(n-1,T.lowerDiagonalAddr(),1,S->addr(1,0),S->bands());
  F77NAME(scopy)(n-1,T.upperDiagonalAddr(),1,S->addr(0,1),S->bands());
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator-(
const SymmetricMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<dim;j++) {
    F77NAME(saxpy)(dim-j,float_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(saxpy)(dim-j-1,float_mone_,H.addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const SymmetricMatrix<float,float> &H,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,H.addr(j,j),1,S->addr(j,j),1);
    if (j+1<n) {
      F77NAME(scopy)(n-j-1,H.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(saxpy)(j+1,float_mone_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      if (j>0) {
        F77NAME(saxpy)(j,float_mone_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)-=float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const UpperTrapezoidalMatrix<float,float> &U,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(j+1,U.addr(0,j),1,S->addr(0,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-B.supDiags());
      if (j>0) F77NAME(scopy)(j,U.addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)=float_one_;
      F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      F77NAME(saxpy)(dim-j,float_mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
        S->addr(ibeg,j),1);
      (*S)(j,j)-=float_one_;
      if (j+1<dim) {
        F77NAME(saxpy)(dim-j-1,float_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const LowerTrapezoidalMatrix<float,float> &L,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=float_one_;
      if (j+1<n) F77NAME(scopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      int ibeg=max(0,j-B.supDiags());
      F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
        B.addr(ibeg,j),1,S->addr(ibeg,j),1);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-ibeg+1,addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  S->operator-=(M);
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const Matrix<float,float> &M,const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  S->copy(M);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(saxpy)(min(n-1,j+B.subDiags())-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator*(const BandMatrix<float,float> &B)
const {
  CHECK_SAME(dim,B.dim);
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(dim,nsub+B.nsub,nsup+B.nsup,float_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=max(0,j-B.nsup);k<=min(dim-1,j+B.nsub);k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,B(k,j),
        addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator*(
const UpperHessenbergMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=0;k<=min(dim-1,j+1);k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,H(k,j),
        addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<float,float>* operator*(
const UpperHessenbergMatrix<float,float> &H,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(saxpy)(min(n,k+2),B(k,j),H.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator*(
const DiagonalMatrix<float,float> &D) const {
  CHECK_SAME(dim,D.size(0));
  BandMatrix<float,float> *P=
    OPERATOR_NEW BandMatrix<float,float>(*this);
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsup);
    F77NAME(sscal)(min(dim-1,j+nsub)-ibeg+1,D[j],
      P->addr(ibeg,j),1);
  }
  return P;
}

template<> BandMatrix<float,float>* operator*(
const DiagonalMatrix<float,float> &D,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  BandMatrix<float,float> *P=
    OPERATOR_NEW BandMatrix<float,float>(n,B.subDiags(),B.supDiags());
  P->copy(B);
  int stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    int jbeg=max(0,i-B.subDiags());
    F77NAME(sscal)(min(n-1,i+B.supDiags())-jbeg+1,D[i],
      P->addr(i,jbeg),stride);
  }
  return P;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator*(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(dim,nsub+1,nsup+1,float_zero_);
  for (int j=0;j<dim;j++) {
    if (j>0) {
      int ibeg=max(0,j-1-nsup);
      F77NAME(saxpy)(min(dim-1,nsub+j-1)-ibeg+1,
        T.upperDiagonalValue(j-1),addr(ibeg,j-1),1,P->addr(ibeg,j),1);
    }
    int ibeg=max(0,j-nsup);
    F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,T.diagonalValue(j),
      addr(ibeg,j),1,P->addr(ibeg,j),1);
    if (j<dim-1) {
      int ibeg=max(0,j+1-nsup);
      F77NAME(saxpy)(min(dim-1,j+1+nsub)-ibeg+1,
        T.lowerDiagonalValue(j),addr(ibeg,j+1),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> BandMatrix<float,float>* operator*(
const SymmetricTridiagonalMatrix<float,float> &T,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(n,B.subDiags()+1,B.supDiags()+1,
    float_zero_);
  int B_stride=B.bands()-1;
  int P_stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    if (i>0) {
      int jbeg=max(0,i-1-B.subDiags());
      F77NAME(saxpy)(min(n-1,i-1+B.supDiags())-jbeg+1,
        T.lowerDiagonalValue(i-1),B.addr(i-1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
    int jbeg=max(0,i-B.subDiags());
    F77NAME(saxpy)(min(n-1,i+B.supDiags())-jbeg+1,T.diagonalValue(i),
      B.addr(i,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    if (i<n-1) {
      int jbeg=max(0,i+1-B.subDiags());
      F77NAME(saxpy)(min(n-1,i+1+B.supDiags())-jbeg+1,
        T.upperDiagonalValue(i),B.addr(i+1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
  }
  return P;
}

template<> BandMatrix<float,float>*
BandMatrix<float,float>::operator*(
const TridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(dim,nsub+1,nsup+1,float_zero_);
  for (int j=0;j<dim;j++) {
    if (j>0) {
      int ibeg=max(0,j-1-nsup);
      F77NAME(saxpy)(min(dim-1,nsub+j-1)-ibeg+1,T(j-1,j),
        addr(ibeg,j-1),1,P->addr(ibeg,j),1);
    }
    int ibeg=max(0,j-nsup);
    F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,T(j,j),addr(ibeg,j),1,
      P->addr(ibeg,j),1);
    if (j<dim-1) {
      int ibeg=max(0,j+1-nsup);
      F77NAME(saxpy)(min(dim-1,j+1+nsub)-ibeg+1,T(j+1,j),
        addr(ibeg,j+1),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> BandMatrix<float,float>* operator*(
const TridiagonalMatrix<float,float> &T,
const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(n,B.subDiags()+1,B.supDiags()+1,
    float_zero_);
  int B_stride=B.bands()-1;
  int P_stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    if (i>0) {
      int jbeg=max(0,i-1-B.subDiags());
      F77NAME(saxpy)(min(n-1,i-1+B.supDiags())-jbeg+1,T(i,i-1),
        B.addr(i-1,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    }
    int jbeg=max(0,i-B.subDiags());
    F77NAME(saxpy)(min(n-1,i+B.supDiags())-jbeg+1,T(i,i),
      B.addr(i,jbeg),B_stride,P->addr(i,jbeg),P_stride);
    if (i<n-1) {
      int jbeg=max(0,i+1-B.subDiags());
      F77NAME(saxpy)(min(n-1,i+1+B.supDiags())-jbeg+1,
        T(i,i+1),B.addr(i+1,jbeg),B_stride,
        P->addr(i,jbeg),P_stride);
    }
  }
  return P;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator*(
const SymmetricMatrix<float,float> &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<float,float> *P=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    for (int k=0;k<dim;k++) {
      int ibeg=max(0,k-nsup);
      F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,S(k,j),addr(ibeg,k),1,
        P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<float,float>* operator*(
const SymmetricMatrix<float,float> &S,
const BandMatrix<float,float> &B) {
//compute by bordering: note that
//  [ sigma s^T ] [ beta c^T ]
//  [   s    S  ] [   b   B  ]
//  = [ sigma beta + s^T b , sigma c^T + s^T B ]
//  = [   s   beta +  S  b ,   s  c^T +   S  B ]
  int n=B.size(0),nsub=B.subDiags(),nsup=B.supDiags();
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *P=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  for (int k=n-1;k>=0;k--) {
    int iend=min(n-1,k+nsub);
    if (k<n-1) {
      // s^T B = ( B^T s )^T
      F77NAME(sgbmv)('T',n-k-1,n-k-1,nsub,nsup,
        float_one_,B.addr(k+1-nsup,k+1),B.bands(),S.addr(k+1,k),1,
        float_zero_,P->addr(k,k+1),n);
      // S b: note that
      // [ S_11 , S_21^T ] [ b ] = [ S_11 b ]
      // [ S_21 ,  S_22  ] [ 0 ] = [ S_21 b ]
      if (nsub>0) { // S_11 b
        F77NAME(ssymv)('L',min(n-k-1,nsub),float_one_,S.addr(k+1,k+1),n,
          B.addr(k+1,k),1,float_zero_,P->addr(k+1,k),1);
        if (nsub<n-k-1) { // S_21 b
          F77NAME(sgemv)('N',n-k-1-nsub,nsub,float_one_,
            S.addr(k+1+nsub,k+1),n,B.addr(k+1,k),1,float_zero_,
            P->addr(k+1+nsub,k),1);
        }
      }
    }
    if (iend>k) { // s^T b
      (*P)(k,k)=F77NAME(sdot)(iend-k,S.addr(k+1,k),1,B.addr(k+1,k),1);
    }
    int jend=min(n-1,k+nsup);
    F77NAME(sger)(n-k,jend-k+1,float_one_,S.addr(k,k),1,
      B.addr(k,k),B.bands()-1,P->addr(k,k),n);
  }
  return P;
}

template<> Matrix<float,float>*
BandMatrix<float,float>::operator*(
const UpperTrapezoidalMatrix<float,float> &U) const {
  CHECK_SAME(dim,U.size(0));
  int n=U.size(1);
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(dim,U.size(1),float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,dim-1);k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,dim);k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
      if (j<dim) {
        int ibeg=max(0,j-nsup);
        F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,float_one_,
          addr(ibeg,j),1,M->addr(ibeg,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<float,float>* operator*(
const UpperTrapezoidalMatrix<float,float> &U,
const BandMatrix<float,float> &B) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(saxpy)(min(k+1,m),B(k,j),U.addr(0,k),1,M->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(saxpy)(min(k,m),B(k,j),U.addr(0,k),1,M->addr(0,j),1);
        if (k<m) (*M)(k,j)+=B(k,j);
      }
    }
  }
  return M;
}

template<> Matrix<float,float>*
BandMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &L) const {
  CHECK_SAME(dim,L.size(0));
  int n=L.size(1);
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(dim,L.size(1),float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<dim;k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-nsup);
      F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,float_one_,addr(ibeg,j),1,
        M->addr(ibeg,j),1);
      for (int k=j+1;k<dim;k++) {
        int ibeg=max(0,k-nsup);
        F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          M->addr(ibeg,j),1);
      }
    }
  }
  return M;
}

template<> Matrix<float,float>* operator*(
const LowerTrapezoidalMatrix<float,float> &L,
const BandMatrix<float,float> &B) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,float> *M=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(saxpy)(m-k,B(k,j),L.addr(k,k),1,M->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
        (*M)(k,j)+=B(k,j);
        if (k+1<m) {
          F77NAME(saxpy)(m-k-1,B(k,j),L.addr(k+1,k),1,M->addr(k+1,j),1);
        }
      }
    }
  }
  return M;
}

template<> SquareMatrix<float,float>*
BandMatrix<float,float>::operator*(
const SquareMatrix<float,float> &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(sgbmv)('N',dim,dim,nsub,nsup,float_one_,addr(),nt,
      S.addr(0,j),1,float_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> SquareMatrix<float,float>* operator*(
const SquareMatrix<float,float> &S,const BandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(saxpy)(n,B(k,j),S.addr(0,j),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Matrix<float,float>* BandMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  CHECK_SAME(dim,M.size(0));
  int n=M.size(1);
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(dim,n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(sgbmv)('N',dim,dim,nsub,nsup,float_one_,addr(),nt,
      M.addr(0,j),1,float_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,const BandMatrix<float,float> &B) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.supDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(saxpy)(m,B(k,j),M.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Vector<float,float>* BandMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  CHECK_SAME(dim,v.size());
  Vector<float,float> *w=
    OPERATOR_NEW Vector<float,float>(dim,float_zero_);
  F77NAME(sgbmv)('N',dim,dim,nsub,nsup,float_one_,addr(),nt,
    v.addr(),1,float_zero_,w->addr(),1);
  return w;
}

template<> void BandMatrix<float,float>::gbmv(float alpha,
const Vector<float,float> &x,float beta,Vector<float,float> &b,
char trans) const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  F77NAME(sgbmv)(trans,dim,dim,nsub,nsup,alpha,addr(),nt,x.addr(),1,
    beta,b.addr(),1);
}

template<> void BandMatrix<float,float>::gbmm(float alpha,
const Matrix<float,float> &X,float beta,Matrix<float,float> &B,
char side,char trans) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      F77NAME(sgbmv)(trans,dim,dim,nsub,nsup,alpha,addr(),nt,
        X.addr(0,j),1,beta,B.addr(0,j),1);
    }
  } else {
    CHECK_SAME(n,dim);
    if (abs(beta)==float_zero_) B=float_zero_;
    else B*=beta;
    if (trans=='N' || trans=='n') {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsup);k<=min(dim-1,j+nsub);k++) {
          F77NAME(saxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    } else {
      for (int j=0;j<n;j++) {
        for (int k=max(0,j-nsub);k<=min(dim-1,j+nsup);k++) {
          F77NAME(saxpy)(m,(*this)(j,k)*alpha,X.addr(0,k),1,
            B.addr(0,j),1);
        }
      }
    }
  }
}

/*
template<> BandMatrix<float,float>*
BandMatrix<float,float>::transpose() const {
  BandMatrix<float,float> *X=
    OPERATOR_NEW BandMatrix<float,float>(dim,nsup,nsub);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsub);
    F77NAME(scopy)(min(dim-1,j+nsup)-i+1,addr(j,i),nsub+nsup,
      X->addr(i,j),1);
  }
  return X;
}
*/

template<> float BandMatrix<float,float>::equilibrate(
Vector<float,float> &r,Vector<float,float> &c,float &rowcnd,
float &colcnd) const {
  CHECK_SAME(dim,r.size());
  CHECK_SAME(dim,c.size());
  float amax;
  int info;
  F77NAME(sgbequ)(dim,dim,nsub,nsup,addr(),nt,r.addr(),c.addr(),
    rowcnd,colcnd,amax,info);
  CHECK_TEST(info==0);
  return amax;
}

template<> float BandMatrix<float,float>::normFrobenius() const {
  float *work=0;
  return F77NAME(slangb)('F',dim,nsub,nsup,addr(),nt,work);
}

template<> float BandMatrix<float,float>::normInfinity() const {
  float *work=OPERATOR_NEW_BRACKET(float,dim);
  float val=F77NAME(slangb)('I',dim,nsub,nsup,addr(),nt,work);
  delete work;
  return val;
}

template<> float BandMatrix<float,float>::normMaxEntry() const {
  float *work=0;
  return F77NAME(slangb)('M',dim,nsub,nsup,addr(),nt,work);
}

template<> float BandMatrix<float,float>::normOne() const {
  float *work=0;
  return F77NAME(slangb)('O',dim,nsub,nsup,addr(),nt,work);
}

template<> float BandMatrix<float,float>::reciprocalConditionNumber(
char norm) const {
  BandMatrix<float,float> *BF=
    OPERATOR_NEW BandMatrix<float,float>(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(sgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);

  float anorm=(norm=='I' || norm=='i' ? normInfinity() : normOne());
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,3*dim);
  int *iwork=OPERATOR_NEW_BRACKET(int,dim);
  F77NAME(sgbcon)(norm,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,anorm,
    rcond,work,iwork,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete [] iwork; iwork=0;
  delete [] ipiv; ipiv=0;
  delete BF; BF=0;
  return rcond;
}

template<> void BandMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char trans)
const {
  CHECK_SAME(dim,b.size())
  CHECK_SAME(dim,x.size())
  BandMatrix<float,float> *BF=
    OPERATOR_NEW BandMatrix<float,float>(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(sgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);

  x.copy(b);
  F77NAME(sgbtrs)(trans,dim,nsub,nsup,1,BF->addr(),BF->nt,ipiv,
    x.addr(),dim,info);
  CHECK_TEST(info==0);
  delete [] ipiv;
  delete BF; BF=0;
}

template<> void BandMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,
char trans) const {
  BandMatrix<float,float> *BF=
    OPERATOR_NEW BandMatrix<float,float>(dim,nsub,nsub+nsup);
  for (int j=0;j<dim;j++) {
    int i=max(0,j-nsup);
    F77NAME(scopy)(min(dim-1,j+nsub)-i+1,addr(i,j),1,BF->addr(i,j),1);
  }
  int *ipiv=OPERATOR_NEW_BRACKET(int,dim);
  int info;
  F77NAME(sgbtrf)(dim,dim,nsub,nsup,BF->addr(),BF->nt,ipiv,info);
  CHECK_TEST(info==0);
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(dim,B.size(0))
    CHECK_SAME(dim,X.size(0))
    CHECK_SAME(nrhs,X.size(1))
    X.copy(B);
    F77NAME(sgbtrs)(trans,dim,nsub,nsup,nrhs,BF->addr(),BF->nt,ipiv,
      X.addr(),dim,info);
    CHECK_TEST(info==0);
  } else { // sgbtrs
    int nrhs=B.size(0);
    CHECK_SAME(dim,B.size(1))
    CHECK_SAME(dim,X.size(1))
    CHECK_SAME(nrhs,X.size(0))
    X.copy(B);
    int kd=nsub+nsup+1;
    if (trans!='N' && trans!='n') {
      if (nsub>0) { // solve L Y^T = B^T: sgbtrs
        for (int j=0;j<dim-1;j++) {
          int lm=min(nsub,dim-j-1);
          int jp=ipiv[j]-1;
          if (jp!=j) F77NAME(sswap)(nrhs,X.addr(0,jp),1,X.addr(0,j),1);
          F77NAME(sger)(nrhs,lm,float_mone_,BF->addr(j+1,j),1,
            X.addr(0,j),1,X.addr(0,j+1),1);
        }
      }
      for (int i=0;i<nrhs;i++) { // solve R X^T = Y^T: stbsv loop 20
        for (int j=dim-1;j>=0;j--) {
          if (abs(X(i,j))>float_zero_) {
            X(i,j)/=(*BF)(j,j);
            int ii=max(0,j-nsub-nsup);
            F77NAME(saxpy)(j-ii,-X(i,j),BF->addr(ii,j),1,
              X.addr(i,ii),nrhs);
          }
        }
      }
    } else {
      for (int i=0;i<nrhs;i++) { // solve R^T Y^T=B^T: stbsv loop 100
        for (int j=0;j<dim;j++) {
          int ii=max(0,j-nsub-nsup);
          X(i,j)=(X(i,j)
            -F77NAME(sdot)(j-i,BF->addr(i,j),1,X.addr(i,ii),1))
            /(*BF)(j,j);
        }
      }
      if (nsub>0) { // solve L^T X^T = Y^T: sgbtrs loop 40
        for (int j=dim-2;j>=0;j--) {
          int lm=min(nsub,dim-1-j);
          F77NAME(sgemv)('T',nrhs,lm,float_mone_,X.addr(0,j+1),nrhs,
            BF->addr(j+1,j),1,float_one_,X.addr(0,j),nrhs);
          int jp=ipiv[j]-1;
          if (jp!=j) {
            F77NAME(sswap)(nrhs,X.addr(0,jp),1,X.addr(0,j),1);
          }
        }
      }
    }
  }

  delete [] ipiv;
  delete BF; BF=0;
}

template class BandMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const float SymmetricBandMatrix<float,float>::outofbounds_
  = float_zero_;
template<> float SymmetricBandMatrix<float,float>::safety_ =
  float_zero_;

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::makeMatrix() const {
  SquareMatrix<float,float> *M=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(dim-j,nsub+1),addr(j,j),1,M->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(scopy)(min(dim-j-1,nsub),addr(j+1,j),1,M->addr(j,j+1),dim);
    }
  }
  return M;
}

template<> void SymmetricBandMatrix<float,float>::fillWith(float d) {
  float *colj=addr();
  for (int j=0;j<dim;j++,colj+=nt) {
    for (int i=0;i<min(dim-j,nt);i++) colj[i]=d;
  }
}

template<> float SymmetricBandMatrix<float,float>::operator()(int i,
int j) const {
  if (i-j>=0 && i-j<=nsub) return *AB->addr(i-j,j);
  if (j-i>0 && j-i<=nsub) return *AB->addr(j-i,i);
  return outofbounds_;
}

template<> SymmetricBandMatrix<float,float>::SymmetricBandMatrix(
const SymmetricTridiagonalMatrix<float,float> &T) : nsub(1), nt(2) {
  dim=T.size(0);
  AB=OPERATOR_NEW Matrix<float,float>(2,dim);
  F77NAME(scopy)(dim,T.diagonalAddr(0),1,addr(0,0),2);
  F77NAME(scopy)(dim-1,T.lowerDiagonalAddr(0),1,addr(1,0),2);
  (*AB)(1,dim-1)=HUGE_VAL;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::makeBandMatrix() const {
  BandMatrix<float,float> *B=
    OPERATOR_NEW BandMatrix<float,float>(dim,nsub,nsub);
  int stride=B->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,B->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,B->addr(j,j+1),stride);
    }
  }
  return B;
}

template<> SymmetricMatrix<float,float>*
SymmetricBandMatrix<float,float>::makeSymmetricMatrix() const {
  SymmetricMatrix<float,float> *S=
    OPERATOR_NEW SymmetricMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<float,float>&
SymmetricBandMatrix<float,float>::operator+=(
const SymmetricBandMatrix<float,float> &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  for (int j=0;j<dim;j++) {
    F77NAME(saxpy)(min(nt,dim-j),float_one_,B.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<float,float>&
SymmetricBandMatrix<float,float>::operator-=(
const SymmetricBandMatrix<float,float> &B) {
  CHECK_SAME(dim,B.dim);
  CHECK_SAME(nsub,B.nsub);
  for (int j=0;j<dim;j++) {
    F77NAME(saxpy)(min(nt,dim-j),float_mone_,B.addr(j,j),1,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<float,float>&
SymmetricBandMatrix<float,float>::operator*=(float d) {
  for (int j=0;j<dim;j++) {
    F77NAME(sscal)(min(nt,dim-j),d,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<float,float>&
SymmetricBandMatrix<float,float>::operator/=(float d) {
  CHECK_TEST(abs(d)>float_zero_);
  float dinv=float_one_/d;
  for (int j=0;j<dim;j++) {
    F77NAME(sscal)(min(nt,dim-j),dinv,addr(j,j),1);
  }
  return *this;
}

template<> SymmetricBandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const SymmetricBandMatrix<float,float> &B) const {
  CHECK_SAME(dim,B.dim);
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(dim,max(nsub,B.nsub),float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(saxpy)(min(B.nt,dim-j),float_one_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const BandMatrix<float,float> &B) const {
  CHECK_SAME(dim,B.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,B.subDiags()),
    max(nsub,B.supDiags()),float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
        S->addr(j,j+1),stride);
    }
    int ibeg=max(0,j-B.supDiags());
    F77NAME(saxpy)(min(dim-1,j+B.subDiags())-ibeg+1,float_one_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const UpperHessenbergMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j+1<dim) {
      F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
    F77NAME(saxpy)(min(dim,j+2),float_one_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const DiagonalMatrix<float,float> &D) const {
  CHECK_SAME(dim,D.size(0));
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(*this);
  F77NAME(saxpy)(dim,float_one_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> SymmetricBandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(dim,max(nsub,1),float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  F77NAME(saxpy)(dim-1,float_one_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  F77NAME(saxpy)(dim,float_one_,T.diagonalAddr(0),1,S->addr(0,0),S->nt);
  return S;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const TridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,1),max(nsub,1),float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,S->addr(j,j+1),stride);
    }
  }
  F77NAME(saxpy)(dim,float_one_,T.diagonalAddr(),1,
    S->addr(0,0),S->bands());
  F77NAME(saxpy)(dim-1,float_one_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S->bands());
  F77NAME(saxpy)(dim-1,float_one_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S->bands());
  return S;
}

template<> SymmetricMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const SymmetricMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SymmetricMatrix<float,float> *S=makeSymmetricMatrix();
  for (int j=0;j<dim;j++) {
    F77NAME(saxpy)(dim-j,float_one_,H.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const UpperTrapezoidalMatrix<float,float> &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      F77NAME(saxpy)(j+1,float_one_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      if (j>0) {
        F77NAME(saxpy)(j,float_one_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)+=float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const LowerTrapezoidalMatrix<float,float> &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      F77NAME(saxpy)(dim-j,float_one_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      (*S)(j,j)+=float_one_;
      if (j+1<dim) {
        F77NAME(saxpy)(dim-j-1,float_one_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator+(
const Matrix<float,float> &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  S->copy(M);
  for (int j=0;j<dim;j++) {
    F77NAME(saxpy)(min(nt,dim-j),float_one_,addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(saxpy)(min(nt,dim-j)-1,float_one_,addr(j+1,j),1,
        S->addr(j,j+1),dim);
    }
  }
  return S;
}

template<> SymmetricBandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const SymmetricBandMatrix<float,float> &B) const {
  CHECK_SAME(dim,B.dim);
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(dim,max(nsub,B.nsub),float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(saxpy)(min(B.nt,dim-j),float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const BandMatrix<float,float> &B) const {
  CHECK_SAME(dim,B.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,B.subDiags()),
      max(nsub,B.supDiags()),float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
        S->addr(j,j+1),stride);
    }
    int ibeg=max(0,j-B.supDiags());
    F77NAME(saxpy)(min(dim-1,j+B.subDiags())-ibeg+1,float_mone_,
      B.addr(ibeg,j),1,S->addr(ibeg,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>* operator-(
const BandMatrix<float,float> &B,
const SymmetricBandMatrix<float,float> &H) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(n,max(H.subDiags(),B.subDiags()),
    max(H.subDiags(),B.supDiags()),float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-B.supDiags());
    F77NAME(scopy)(min(n-1,j+B.subDiags())-ibeg+1,B.addr(ibeg,j),1,
      S->addr(ibeg,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(H.bands(),n-j),float_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(min(H.bands(),n-j)-1,float_mone_,
        H.addr(j+1,j),1,S->addr(j,j+1),stride);
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const UpperHessenbergMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,S->addr(j,j+1),dim);
    }
    F77NAME(saxpy)(min(dim,j+2),float_mone_,H.addr(0,j),1,
      S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const UpperHessenbergMatrix<float,float> &H,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(n,j+2),H.addr(0,j),1,S->addr(0,j),1);
  }
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(min(B.bands(),n-j)-1,float_mone_,
        B.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  return S;
}

template<> SymmetricBandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const DiagonalMatrix<float,float> &D) const {
  CHECK_SAME(dim,D.size(0));
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(*this);
  F77NAME(saxpy)(dim,float_mone_,D.addr(),1,S->addr(0,0),nt);
  return S;
}

template<> SymmetricBandMatrix<float,float>* operator-(
const DiagonalMatrix<float,float> &D,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,D.size(0));
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(n,B.subDiags(),float_zero_);
  int nt=B.bands();
  F77NAME(scopy)(n,D.addr(),1,S->addr(0,0),nt);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(nt,n-j),float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricBandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(dim,max(nsub,1),float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
  }
  F77NAME(saxpy)(dim-1,float_mone_,T.lowerDiagonalAddr(0),1,
    S->addr(1,0),S->nt);
  F77NAME(saxpy)(dim,float_mone_,T.diagonalAddr(0),1,S->addr(0,0),S->nt);
  return S;
}

template<> SymmetricBandMatrix<float,float>* operator-(
const SymmetricTridiagonalMatrix<float,float> &T,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(n,max(B.subDiags(),1),
    float_zero_);
  F77NAME(scopy)(n-1,T.lowerDiagonalAddr(0),1,S->addr(1,0),S->bands());
  F77NAME(scopy)(n,T.diagonalAddr(0),1,S->addr(0,0),S->bands());
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const TridiagonalMatrix<float,float> &T) const {
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(dim,max(nsub,1),max(nsub,1),float_zero_);
  int stride=S->bands()-1;
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,S->addr(j,j+1),stride);
    }
  }
  int S_nt=S->bands();
  F77NAME(saxpy)(dim,float_mone_,T.diagonalAddr(),1,
    S->addr(0,0),S_nt);
  F77NAME(saxpy)(dim-1,float_mone_,T.lowerDiagonalAddr(),1,
    S->addr(1,0),S_nt);
  F77NAME(saxpy)(dim-1,float_mone_,T.upperDiagonalAddr(),1,
    S->addr(0,1),S_nt);
  return S;
}

template<> BandMatrix<float,float>* operator-(
const TridiagonalMatrix<float,float> &T,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,T.size(0));
  BandMatrix<float,float> *S=OPERATOR_NEW
    BandMatrix<float,float>(n,max(B.subDiags(),1),max(B.subDiags(),1),
    float_zero_);
  int S_nt=S->bands();
  F77NAME(scopy)(n,T.diagonalAddr(),1,S->addr(0,0),S_nt);
  F77NAME(scopy)(n-1,T.lowerDiagonalAddr(),1,S->addr(1,0),S_nt);
  F77NAME(scopy)(n-1,T.upperDiagonalAddr(),1,S->addr(0,1),S_nt);
  int stride=S->bands()-1;
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(min(B.bands(),n-j)-1,float_mone_,B.addr(j+1,j),1,
        S->addr(j,j+1),stride);
    }
  }
  return S;
}

template<> SymmetricMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const SymmetricMatrix<float,float> &H) const {
  CHECK_SAME(dim,H.size(0));
  SymmetricMatrix<float,float> *S=makeSymmetricMatrix();
  for (int j=0;j<dim;j++) {
    F77NAME(saxpy)(dim-j,float_mone_,H.addr(j,j),1,
      S->addr(j,j),1);
//  if (j+1<dim) {
//    F77NAME(saxpy)(dim-j-1,float_mone_,H.addr(j+1,j),1,
//      S->addr(j,j+1),dim);
//  }
  }
  return S;
}

template<> SymmetricMatrix<float,float>* operator-(
const SymmetricMatrix<float,float> &H,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,H.size(0));
  SymmetricMatrix<float,float> *S=OPERATOR_NEW
    SymmetricMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(n-j,H.addr(j,j),1,S->addr(j,j),1);
    F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const UpperTrapezoidalMatrix<float,float> &U) const {
  CHECK_SAME(dim,U.size(0));
  CHECK_SAME(dim,U.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      F77NAME(saxpy)(j+1,float_mone_,U.addr(0,j),1,S->addr(0,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      if (j>0) {
        F77NAME(saxpy)(j,float_mone_,U.addr(0,j),1,S->addr(0,j),1);
      }
      (*S)(j,j)-=float_one_;
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const UpperTrapezoidalMatrix<float,float> &U,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,U.size(0));
  CHECK_SAME(n,U.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(j+1,U.addr(0,j),1,S->addr(0,j),1);
      F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
        S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(min(B.bands(),n-j)-1,float_mone_,
          B.addr(j+1,j),1,S->addr(j,j+1),n);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      if (j>0) F77NAME(scopy)(j,U.addr(0,j),1,S->addr(0,j),1);
      (*S)(j,j)=float_one_;
      F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
        S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(min(B.bands(),n-j)-1,float_mone_,
          B.addr(j+1,j),1,S->addr(j,j+1),n);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const LowerTrapezoidalMatrix<float,float> &L) const {
  CHECK_SAME(dim,L.size(0));
  CHECK_SAME(dim,L.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<dim;j++) {
      F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      F77NAME(saxpy)(dim-j,float_mone_,L.addr(j,j),1,S->addr(j,j),1);
    }
  } else {
    for (int j=0;j<dim;j++) {
      F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
      if (j<dim-1) {
        F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,
          S->addr(j,j+1),dim);
      }
      (*S)(j,j)-=float_one_;
      if (j+1<dim) {
        F77NAME(saxpy)(dim-j-1,float_mone_,L.addr(j+1,j),1,
          S->addr(j+1,j),1);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const LowerTrapezoidalMatrix<float,float> &L,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,L.size(0));
  CHECK_SAME(n,L.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      F77NAME(scopy)(n-j,L.addr(j,j),1,S->addr(j,j),1);
      F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
        S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(min(B.bands(),n-j)-1,float_mone_,
          B.addr(j+1,j),1,S->addr(j,j+1),n);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      (*S)(j,j)=float_one_;
      if (j+1<n) {
        F77NAME(scopy)(n-j-1,L.addr(j+1,j),1,S->addr(j+1,j),1);
      }
      F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
        S->addr(j,j),1);
      if (j<n-1) {
        F77NAME(saxpy)(min(B.bands(),n-j)-1,float_mone_,
          B.addr(j+1,j),1,S->addr(j,j+1),n);
      }
    }
  }
  return S;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator-(
const Matrix<float,float> &M) const {
  CHECK_SAME(dim,M.size(0));
  CHECK_SAME(dim,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(dim,float_zero_);
  for (int j=0;j<dim;j++) {
    F77NAME(scopy)(min(nt,dim-j),addr(j,j),1,S->addr(j,j),1);
    if (j<dim-1) {
      F77NAME(scopy)(min(nt,dim-j)-1,addr(j+1,j),1,S->addr(j,j+1),dim);
    }
    F77NAME(saxpy)(dim,float_mone_,M.addr(0,j),1,S->addr(0,j),1);
  }
  return S;
}

template<> SquareMatrix<float,float>* operator-(
const Matrix<float,float> &M,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,M.size(0));
  CHECK_SAME(n,M.size(1));
  SquareMatrix<float,float> *S=OPERATOR_NEW
    SquareMatrix<float,float>(n,float_zero_);
  S->copy(M);
  for (int j=0;j<n;j++) {
    F77NAME(saxpy)(min(B.bands(),n-j),float_mone_,B.addr(j,j),1,
      S->addr(j,j),1);
    if (j<n-1) {
      F77NAME(saxpy)(min(B.bands(),n-j)-1,float_mone_,
        B.addr(j+1,j),1,S->addr(j,j+1),n);
    }
  }
  return S;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const SymmetricBandMatrix<float,float> &B) const {
// compute by bordering: note that
// [ sigma s^T ] [ tau t^T ] = [ sigma tau + s^T t , sigma t^T + s^T T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^T +  S  T ]
  CHECK_SAME(dim,B.dim);
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(dim,nsub+B.nsub,nsub+B.nsub,float_zero_);
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
          F77NAME(saxpy)(kk-ibeg,B(kk,k),addr(kk,ibeg),nsub,
            P->addr(ibeg,k),1);
        }
//      for (int i=max(kk,ibeg);i<=min(kk+nsub,dim-1);i++) {
//        (*P)(i,k)+=(*this)(i,kk)*B(kk,k);
//      }
        ibeg=max(kk,ibeg);
        F77NAME(saxpy)(min(kk+nsub,dim-1)-ibeg+1,B(kk,k),
          addr(ibeg,kk),1,P->addr(ibeg,k),1);
      }
      for (int j=k+1;j<=min(k+nsub+B.nsub,dim-1);j++) { // s^T T
//      for (int kk=max(k+1,j-B.nsub);kk<min(j,k+nsub);kk++) {
//        (*P)(k,j)+=(*this)(kk,k)*B(j,kk);
//      }
        int kbeg=max(k+1,j-B.nsub);
        int kend=min(j-1,k+nsub);
        if (kbeg<=kend) {
          (*P)(k,j)=F77NAME(sdot)(kend-kbeg+1,addr(kbeg,k),1,
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
            F77NAME(sdot)(kend-kbeg+1,addr(kbeg,k),1,B.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^T t
        F77NAME(sdot)(min(nb,dim-k-1),addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(sger)(min(dim-k,nt),min(dim-k,B.nt),float_one_,
      addr(k,k),1,B.addr(k,k),1,P->addr(k,k),P->bands()-1);
  }
  return P;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const BandMatrix<float,float> &B) const {
// compute by bordering: note that
// [ sigma s^T ] [ beta c^T ] = [ sigma beta + s^T b , sigma c^T + s^T B ]
// [   s    S  ] [  b    B  ] = [   s   beta +  S  b ,   s   c^T +  S  B ]
  CHECK_SAME(dim,B.size(0));
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(dim,nsub+B.subDiags(),nsub+B.supDiags(),
    float_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*B(dim-1,dim-1);
  int nb=min(nsub,B.subDiags());
  for (int k=dim-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(dim-1,k+B.subDiags());kk++) { // S b
        int ibeg=max(k+1,kk-B.supDiags());
        if (kk>ibeg) {
          F77NAME(saxpy)(kk-ibeg,B(kk,k),addr(kk,ibeg),nsub,
            P->addr(ibeg,k),1);
        }
        ibeg=max(kk,ibeg);
        F77NAME(saxpy)(min(kk+nsub,dim-1)-ibeg+1,B(kk,k),
          addr(ibeg,kk),1,P->addr(ibeg,k),1);
      }
      for (int j=k+1;j<=min(k+P->supDiags(),dim-1);j++) { // s^T B
        int kbeg=max(k+1,j-B.supDiags());
        int kend=min(dim-1,min(k+nsub,j+B.subDiags()));
        if (kbeg<=kend) {
          (*P)(k,j)=F77NAME(sdot)(kend-kbeg+1,addr(kbeg,k),1,
            B.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^T b
        F77NAME(sdot)(min(nb,dim-k-1),addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(sger)(min(dim-k,nt),min(dim-k,B.supDiags()+1),float_one_,
      addr(k,k),1,B.addr(k,k),1,P->addr(k,k),P->bands()-1);
  }
  return P;
}

template<> BandMatrix<float,float>* operator*(
const BandMatrix<float,float> &B,
const SymmetricBandMatrix<float,float> &S) {
  int n=S.size(0);
  CHECK_SAME(n,B.size(0));
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(n,S.subDiags()+B.subDiags(),
      S.subDiags()+B.supDiags(),float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-S.subDiags());k<=min(n-1,j+S.subDiags());k++) {
      int ibeg=max(0,k-B.supDiags());
      F77NAME(saxpy)(min(n-1,k+B.subDiags())-ibeg+1,S(k,j),
        B.addr(ibeg,k),1,P->addr(ibeg,j),1);
    }
  }
  return P;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const UpperHessenbergMatrix<float,float> &H) const {
// compute by bordering: note that
// [ sigma s^T ] [     eta_11 h^T ]
// [   s    S  ] [ e_0 eta_21  H  ]
// = [ sigma eta_11 + s^T e_0 eta_21 , sigma h^T + s^T H ]
// = [   s   eta_11 +  S  e_0 eta_21 ,   s   h^T +  S  H ]
  CHECK_SAME(dim,H.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(dim,float_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*H(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nsub>0) {
      F77NAME(saxpy)(min(dim-k-1,nsub+1),H(k+1,k),addr(k+1,k+1),1,
        P->addr(k+1,k),1); // S  e_0 eta_21
      for (int j=k+1;j<dim;j++) { // s^T H
        (*P)(k,j)=
          F77NAME(sdot)(min(dim-k-1,min(j-k+1,nsub)),addr(k+1,k),1,
            H.addr(k+1,j),1);
      }
      (*P)(k,k)=(*this)(k+1,k)*H(k+1,k); // s^T e_0 eta_21
    }
    F77NAME(sger)(min(dim-k,nsub+1),dim-k,float_one_,addr(k,k),1,
      H.addr(k,k),dim,P->addr(k,k),dim);
  }
  return P;
}

template<> SquareMatrix<float,float>* operator*(
const UpperHessenbergMatrix<float,float> &H,
const SymmetricBandMatrix<float,float> &S) {
  int n=S.size(0);
  CHECK_SAME(n,H.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-S.subDiags());k<=min(n-1,j+S.subDiags());k++) {
      F77NAME(saxpy)(min(n,k+2),S(k,j),H.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const DiagonalMatrix<float,float> &D) const {
  BandMatrix<float,float> *P=makeBandMatrix();
  for (int j=0;j<dim;j++) {
    int ibeg=max(0,j-nsub);
    F77NAME(sscal)(min(dim-1,j+nsub)-ibeg+1,D[j],P->addr(ibeg,j),1);
  }
  return P;
}

template<> BandMatrix<float,float>* operator*(
const DiagonalMatrix<float,float> &D,
const SymmetricBandMatrix<float,float> &S) {
  int n=S.size(0);
  BandMatrix<float,float> *P=S.makeBandMatrix();
  int stride=P->bands()-1;
  for (int i=0;i<n;i++) {
    int jbeg=max(0,i-S.subDiags());
    F77NAME(sscal)(min(n-1,i+S.subDiags())-jbeg+1,D[i],
      P->addr(i,jbeg),stride);
  }
  return P;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const SymmetricTridiagonalMatrix<float,float> &T) const {
// compute by bordering: note that
// [ sigma s^T ] [     tau    , lambda e_0^T ]
// [   s    S  ] [ e_0 lambda ,           T  ]
//   = [ sigma tau + s^T e_0 lambda , sigma lambda e_0^T + s^T T ]
//   = [   s   tau +  S  e_o lambda ,   s   lambda e_0^T +  S  T ]
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(dim,nsub+1,nsub+1,float_zero_);
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
    F77NAME(saxpy)(min(dim-k-1,nsub+1),T.lowerDiagonalValue(k),
      addr(k+1,k+1),1,P->addr(k+1,k),1); // S e_0 lambda
    F77NAME(saxpy)(min(dim-k,nsub+1),T.diagonalValue(k),addr(k,k),1,
      P->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(saxpy)(min(dim-k,nsub+1),T.upperDiagonalValue(k),addr(k,k),1,
      P->addr(k,k+1),1);
     // [ sigma ] lambda
     // [   s   ]
  }
  return P;
}

template<> BandMatrix<float,float>* operator*(
const SymmetricTridiagonalMatrix<float,float> &T,
const SymmetricBandMatrix<float,float> &S) {
// compute by bordering: note that
// [     tau    , lambda e_0^T ] [ sigma S^T ]
// [ e_0 lambda ,      T       ] [   s    S  ]
//   = [ tau sigma + lambda e_0^T s , tau s^T        + lambda e_0^T S ]
//   = [ e_0 lambda sigma +     T s , e_0 lambda s^T +            T S ]
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=S.subDiags();
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(n,nb+1,nb+1,float_zero_);
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
    F77NAME(saxpy)(min(n-k-1,nb+1),T.lowerDiagonalValue(k),
      S.addr(k+1,k+1),1,P->addr(k,k+1),stride);
      // lambda e_0^T S = lambda ( S e_0 )^T
    F77NAME(saxpy)(min(n-k,nb+1),T.diagonalValue(k),
      S.addr(k,k),1,P->addr(k,k),stride);
      // tau [ sigma , s^T ]
    F77NAME(saxpy)(min(n-k,nb+1),T.lowerDiagonalValue(k),S.addr(k,k),1,
      P->addr(k+1,k),stride);
      // e_0 lambda [ sigma , s^T ]
  }
  return P;
}

template<> BandMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const TridiagonalMatrix<float,float> &T) const {
// compute by bordering: note that
// [ sigma s^T ] [    tau     upsilon e_0^T ]
// [   s    S  ] [ e_0 lambda       T       ]
//   = [ sigma tau + s^T e_0 lambda , sigma upsilon e_0^T + s^T T ]
//   = [     s tau +   S e_0 lambda ,     s upsilon e_0^T +   S T ]
  CHECK_SAME(dim,T.size(0));
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(dim,nsub+1,nsub+1,float_zero_);
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
    F77NAME(saxpy)(min(dim-k-1,nsub+1),T(k+1,k),addr(k+1,k+1),1,
      P->addr(k+1,k),1); // S e_0 lambda
    F77NAME(saxpy)(min(dim-k,nsub+1),T(k,k),addr(k,k),1,P->addr(k,k),1);
      // [sigma ] tau
      // [  s   ]
    F77NAME(saxpy)(min(dim-k,nsub+1),T(k,k+1),addr(k,k),1,
      P->addr(k,k+1),1);
      // [ sigma ] lambda
      // [   s   ]
  }
  return P;
}

template<> BandMatrix<float,float>* operator*(
const TridiagonalMatrix<float,float> &T,
const SymmetricBandMatrix<float,float> &S) {
// compute by bordering: note that
// [     tau    , upsilon e_0^T ] [ sigma S^T ]
// [ e_0 lambda ,       T       ] [   s    S  ]
//   = [ tau sigma + upsilon e_0^T s , tau s^T        + upsilon e_0^T S ]
//   = [ e_0 lambda sigma +      T s , e_0 lambda s^T +             T S ]
  int n=S.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=S.subDiags();
  BandMatrix<float,float> *P=OPERATOR_NEW
    BandMatrix<float,float>(n,nb+1,nb+1,float_zero_);
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
    F77NAME(saxpy)(min(n-k-1,nb+1),T(k,k+1),S.addr(k+1,k+1),1,
      P->addr(k,k+1),stride);
      // upsilon e_0^T S = upsilon ( S e_0 )^T
    F77NAME(saxpy)(min(n-k,nb+1),T(k,k),S.addr(k,k),1,
      P->addr(k,k),stride);
      // tau [ sigma , s^T ]
    F77NAME(saxpy)(min(n-k,nb+1),T(k+1,k),S.addr(k,k),1,
      P->addr(k+1,k),stride);
      // e_0 lambda [ sigma , s^T ]
  }
  return P;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const SymmetricMatrix<float,float> &S) const {
// compute by bordering: note that
// [ sigma s^T ] [ tau t^T ] = [ sigma tau + s^T t , sigma t^T + s^T T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^T +  S  T ]
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(dim,float_zero_);
  (*P)(dim-1,dim-1)=(*this)(dim-1,dim-1)*S(dim-1,dim-1);
  for (int k=dim-2;k>=0;k--) {
    if (nsub>0) {
      F77NAME(ssbmv)('L',dim-k-1,nsub,float_one_,addr(k+1,k+1),nt,
        S.addr(k+1,k),1,float_zero_,P->addr(k+1,k),1); // S t
      for (int j=k+1;j<dim;j++) { // s^T T
        int kend=min(j-1,k+nsub);
        if (k<kend) {
          (*P)(k,j)=F77NAME(sdot)(kend-k,addr(k+1,k),1,S.addr(j,k+1),dim);
        }
        int kbeg=max(k+1,j);
        kend=min(dim-1,k+nsub);
        if (kbeg<=kend) {
          (*P)(k,j)+=
            F77NAME(sdot)(kend-kbeg+1,addr(kbeg,k),1,S.addr(kbeg,j),1);
        }
      }
      (*P)(k,k)= // s^T t
        F77NAME(sdot)(min(nsub,dim-k-1),addr(k+1,k),1,S.addr(k+1,k),1);
    }
    F77NAME(sger)(min(dim-k,nt),dim-k,float_one_,addr(k,k),1,
      S.addr(k,k),1,P->addr(k,k),dim);
  }
  return P;
}

template<> SquareMatrix<float,float>* operator*(
const SymmetricMatrix<float,float> &S,
const SymmetricBandMatrix<float,float> &B) {
// compute by bordering: note that
// [ sigma s^T ] [ tau t^T ] = [ sigma tau + s^T t , sigma t^T + s^T T ]
// [   s    S  ] [  t   T  ] = [   s   tau +  S  t ,   s   t^T +  S  T ]
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  (*P)(n-1,n-1)=S(n-1,n-1)*B(n-1,n-1);
  int nb=B.subDiags();
  for (int k=n-2;k>=0;k--) {
    if (nb>0) {
      for (int kk=k+1;kk<=min(n-1,k+nb);kk++) { // S t 
        int ibeg=k+1;
        if (kk>ibeg) {
          F77NAME(saxpy)(kk-ibeg,B(kk,k),S.addr(kk,ibeg),n,
            P->addr(ibeg,k),1);
        }
        F77NAME(saxpy)(n-kk,B(kk,k),S.addr(kk,kk),1,P->addr(kk,k),1);
      }
      for (int j=k+1;j<n;j++) { // s^T T
        int kbeg=max(k+1,j-nb);
        if (kbeg<j) {
          (*P)(k,j)=F77NAME(sdot)(j-kbeg,S.addr(kbeg,k),1,
            B.addr(j,kbeg),B.bands()-1);
        }
        kbeg=max(kbeg,j);
        (*P)(k,j)+=F77NAME(sdot)(min(n-1,j+nb)-kbeg+1,
          S.addr(kbeg,k),1,B.addr(kbeg,j),1);
      }
      (*P)(k,k)= // s^T t
        F77NAME(sdot)(min(nb,n-k-1),S.addr(k+1,k),1,B.addr(k+1,k),1);
    }
    F77NAME(sger)(n-k,min(n-k,B.bands()),float_one_,
      S.addr(k,k),1,B.addr(k,k),1,P->addr(k,k),n);
  }
  return P;
}

template<> Matrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const UpperTrapezoidalMatrix<float,float> &U) const {
  CHECK_SAME(dim,U.size(0));
  int n=U.size(1);
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(dim,n,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=0;k<=min(j,dim-1);k++) {
        int ibeg=max(0,k-nsub);
        if (ibeg<k) {
          F77NAME(saxpy)(k-ibeg,U(k,j),addr(k,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,k);
        F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=0;k<min(j,dim);k++) {
        int ibeg=max(0,k-nsub);
        if (ibeg<k) {
          F77NAME(saxpy)(k-ibeg,U(k,j),addr(k,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,k);
        F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,U(k,j),addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
      if (j<dim) {
        int ibeg=max(0,j-nsub);
        if (ibeg<j) {
          F77NAME(saxpy)(j-ibeg,float_one_,addr(j,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,j);
        F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,float_one_,
          addr(ibeg,j),1,P->addr(ibeg,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<float,float>* operator*(
const UpperTrapezoidalMatrix<float,float> &U,
const SymmetricBandMatrix<float,float> &B) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  if (dynamic_cast<const
  UnitUpperTrapezoidalMatrix<float,float>*>(&U)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(saxpy)(min(k+1,m),B(k,j),U.addr(0,k),1,P->addr(0,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(saxpy)(min(k,m),B(k,j),U.addr(0,k),1,P->addr(0,j),1);
        if (k<m) (*P)(k,j)+=B(k,j);
      }
    }
  }
  return P;
}

template<> Matrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const LowerTrapezoidalMatrix<float,float> &L) const {
  CHECK_SAME(dim,L.size(0));
  int n=L.size(1);
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(dim,n,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=j;k<dim;k++) {
        int ibeg=max(0,k-nsub);
        if (ibeg<k) {
          F77NAME(saxpy)(k-ibeg,L(k,j),addr(k,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,k);
        F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      int ibeg=max(0,j-nsub);
      if (ibeg<j) {
        F77NAME(saxpy)(j-ibeg,float_one_,addr(j,ibeg),nt-1,
          P->addr(ibeg,j),1);
      }
      ibeg=max(ibeg,j);
      F77NAME(saxpy)(min(dim-1,j+nsub)-ibeg+1,float_one_,addr(ibeg,j),1,
        P->addr(ibeg,j),1);
      for (int k=j+1;k<dim;k++) {
        int ibeg=max(0,k-nsub);
        if (ibeg<k) {
          F77NAME(saxpy)(k-ibeg,L(k,j),addr(k,ibeg),nt-1,
            P->addr(ibeg,j),1);
        }
        ibeg=max(ibeg,k);
        F77NAME(saxpy)(min(dim-1,k+nsub)-ibeg+1,L(k,j),addr(ibeg,k),1,
          P->addr(ibeg,j),1);
      }
    }
  }
  return P;
}

template<> Matrix<float,float>* operator*(
const LowerTrapezoidalMatrix<float,float> &L,
const SymmetricBandMatrix<float,float> &B) {
int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  if (dynamic_cast<const
  UnitLowerTrapezoidalMatrix<float,float>*>(&L)==0) {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        F77NAME(saxpy)(m-k,B(k,j),L.addr(k,k),1,P->addr(k,j),1);
      }
    }
  } else {
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
        (*P)(k,j)+=B(k,j);
        if (k+1<m) {
          F77NAME(saxpy)(m-k-1,B(k,j),L.addr(k+1,k),1,P->addr(k+1,j),1);
        }
      }
    }
  }
  return P;
}

template<> SquareMatrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const SquareMatrix<float,float> &S) const {
  CHECK_SAME(dim,S.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(dim);
  for (int j=0;j<dim;j++) {
    F77NAME(ssbmv)('L',dim,nsub,float_one_,addr(),nt,
      S.addr(0,j),1,float_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> SquareMatrix<float,float>* operator*(
const SquareMatrix<float,float> &S,
const SymmetricBandMatrix<float,float> &B) {
  int n=B.size(0);
  CHECK_SAME(n,S.size(0));
  SquareMatrix<float,float> *P=
    OPERATOR_NEW SquareMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(saxpy)(n,B(k,j),S.addr(0,j),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Matrix<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const Matrix<float,float> &M) const {
  CHECK_SAME(dim,M.size(0));
  int n=M.size(1);
  Matrix<float,float> *P=OPERATOR_NEW Matrix<float,float>(dim,n);
  for (int j=0;j<dim;j++) {
    F77NAME(ssbmv)('L',dim,nsub,float_one_,addr(),nt,
      M.addr(0,j),1,float_zero_,P->addr(0,j),1);
  }
  return P;
}

template<> Matrix<float,float>* operator*(
const Matrix<float,float> &M,
const SymmetricBandMatrix<float,float> &B) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,B.size(0));
  Matrix<float,float> *P=
    OPERATOR_NEW Matrix<float,float>(m,n,float_zero_);
  for (int j=0;j<n;j++) {
    for (int k=max(0,j-B.subDiags());k<=min(n-1,j+B.subDiags());k++) {
      F77NAME(saxpy)(m,B(k,j),M.addr(0,k),1,P->addr(0,j),1);
    }
  }
  return P;
}

template<> Vector<float,float>*
SymmetricBandMatrix<float,float>::operator*(
const Vector<float,float> &v) const {
  CHECK_SAME(dim,v.size());
  Vector<float,float> *p=OPERATOR_NEW Vector<float,float>(dim);
  F77NAME(ssbmv)('L',dim,nsub,float_one_,addr(),nt,v.addr(),1,
    float_zero_,p->addr(),1);
  return p;
}

template<> void SymmetricBandMatrix<float,float>::sbmv(float alpha,
const Vector<float,float> &x,float beta,Vector<float,float> &b)
const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  F77NAME(ssbmv)('L',dim,nsub,alpha,addr(),nt,x.addr(),1,beta,
    b.addr(),1);
}

template<> void SymmetricBandMatrix<float,float>::sbmm(float alpha,
const Matrix<float,float> &X,float beta,Matrix<float,float> &B,
char side) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      F77NAME(ssbmv)('L',dim,nsub,alpha,addr(),nt,X.addr(0,j),1,beta,
        B.addr(0,j),1);
    }
  } else {
    CHECK_SAME(n,dim);
    if (abs(beta)==float_zero_) B=float_zero_;
    else B*=beta;
    for (int j=0;j<n;j++) {
      for (int k=max(0,j-nsub);k<=min(dim-1,j+nsub);k++) {
        F77NAME(saxpy)(m,(*this)(k,j)*alpha,X.addr(0,k),1,
          B.addr(0,j),1);
      }
    }
  }
}

template<> float SymmetricBandMatrix<float,float>::normFrobenius()
const {
  float *work=0;
  return F77NAME(slansb)('F','L',dim,nsub,addr(),nt,work);
}

template<> float SymmetricBandMatrix<float,float>::normInfinity()
const {
  float *work=OPERATOR_NEW_BRACKET(float,dim);;
  float val=F77NAME(slansb)('I','L',dim,nsub,addr(),nt,work);
  delete [] work; work=0;
  return val;
}

template<> float SymmetricBandMatrix<float,float>::normMaxEntry()
const {
  float *work=0;
  return F77NAME(slansb)('M','L',dim,nsub,addr(),nt,work);
}

template<> float SymmetricBandMatrix<float,float>::normOne() const {
  float *work=OPERATOR_NEW_BRACKET(float,dim);;
  float val=F77NAME(slansb)('O','L',dim,nsub,addr(),nt,work);
  delete [] work; work=0;
  return val;
}

template<> Vector<float,float>* 
SymmetricBandMatrix<float,float>::eigenvalues(
OrthogonalMatrix<float,float> *&Q) const {
  if (Q!=0) CHECK_SAME(dim,Q->size(0));
  char jobz=(Q==0 ? 'N' : 'V');
  Vector<float,float> *lambda=OPERATOR_NEW Vector<float,float>(dim);
  SymmetricBandMatrix<float,float> *copy=
    OPERATOR_NEW SymmetricBandMatrix<float,float>(*this);
  float *work=OPERATOR_NEW_BRACKET(float,3*dim-2);
  int info;
  float *qa=( Q==0 ? 0 : Q->addr() );
  F77NAME(ssbev)(jobz,'L',dim,nsub,copy->addr(),copy->bands(),
    lambda->addr(),qa,dim,work,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete copy; copy=0;
  return lambda;
}

template class SymmetricBandMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> SymmetricPositiveBandMatrix<float,float>*
SymmetricPositiveBandMatrix<float,float>::operator+(
const SymmetricPositiveBandMatrix<float,float> &B) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,B.size(0));
  SymmetricPositiveBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricPositiveBandMatrix<float,float>(n,max(ns,B.subDiags()),
    float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(saxpy)(min(B.bands(),n-j),float_one_,B.addr(j,j),1,
      S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricPositiveBandMatrix<float,float>*
SymmetricPositiveBandMatrix<float,float>::operator+(
const SymmetricPositiveTridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricPositiveBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricPositiveBandMatrix<float,float>(n,max(ns,1),
    float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=T(j,j);
    if (j+1<n) (*S)(j+1,j)+=T(j+1,j);
  }
  return S;
}

template<> SymmetricBandMatrix<float,float>*
SymmetricPositiveBandMatrix<float,float>::operator+(
const SymmetricTridiagonalMatrix<float,float> &T) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricBandMatrix<float,float> *S=OPERATOR_NEW
    SymmetricBandMatrix<float,float>(n,max(ns,1),float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    (*S)(j,j)+=T(j,j);
    if (j+1<n) (*S)(j+1,j)+=T(j+1,j);
  }
  return S;
}

template<> SymmetricPositiveMatrix<float,float>*
SymmetricPositiveBandMatrix<float,float>::operator+(
const SymmetricPositiveMatrix<float,float> &M) const {
  int n=size(0);
  int nb=bands();
  CHECK_SAME(n,M.size(0));
  SymmetricPositiveMatrix<float,float> *S=OPERATOR_NEW
    SymmetricPositiveMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(saxpy)(n-j,float_one_,M.addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> SymmetricMatrix<float,float>*
SymmetricPositiveBandMatrix<float,float>::operator+(
const SymmetricMatrix<float,float> &T) const {
  int n=size(0);
  int nb=bands();
  CHECK_SAME(n,T.size(0));
  SymmetricMatrix<float,float> *S=OPERATOR_NEW
    SymmetricMatrix<float,float>(n,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(nb,n-j),addr(j,j),1,S->addr(j,j),1);
    F77NAME(saxpy)(n-j,float_one_,T.addr(j,j),1,S->addr(j,j),1);
  }
  return S;
}

template<> float
SymmetricPositiveBandMatrix<float,float>::equilibrate(
Vector<float,float> &s,float &scond) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,s.size());
  float amax;
  int info;
  F77NAME(spbequ)('L',n,ns,addr(),nb,s.addr(),scond,amax,info);
  CHECK_TEST(info==0);
  return amax;
}

template<> float
SymmetricPositiveBandMatrix<float,float>::reciprocalConditionNumber()
const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  SymmetricPositiveBandMatrix<float,float> *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<float,float>(*this);
  int info;
  F77NAME(spbtrf)('L',n,ns,BF->addr(),nb,info);
  CHECK_TEST(info==0);

  float anorm=normOne();
  float rcond=HUGE_VAL;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(spbcon)('L',n,ns,BF->addr(),nb,anorm,rcond,work,iwork,info);
  delete [] work; work=0;
  delete [] iwork; iwork=0;
  delete BF; BF=0;
  return rcond;
}

template<> void SymmetricPositiveBandMatrix<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,char) const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  CHECK_SAME(n,b.size())
  CHECK_SAME(n,x.size())
  SymmetricPositiveBandMatrix<float,float> *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<float,float>(*this);
  x.copy(b);
  int info;
  F77NAME(spbsv)('L',n,ns,1,BF->addr(),nb,x.addr(),n,info);
  CHECK_TEST(info==0);
  delete BF; BF=0;
}

template<> void SymmetricPositiveBandMatrix<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,char side,char)
const {
  int n=size(0);
  int ns=subDiags();
  int nb=bands();
  SymmetricPositiveBandMatrix<float,float> *BF=
    OPERATOR_NEW SymmetricPositiveBandMatrix<float,float>(*this);
  int info;
  if (side=='L' || side=='l') {
    int nrhs=B.size(1);
    CHECK_SAME(n,B.size(0))
    CHECK_SAME(n,X.size(0))
    CHECK_SAME(nrhs,X.size(1))
    X.copy(B);
    F77NAME(spbsv)('L',n,ns,1,BF->addr(),nb,X.addr(),n,info);
    CHECK_TEST(info==0);
  } else {
    int nrhs=B.size(0);
    CHECK_SAME(n,B.size(1))
    CHECK_SAME(n,X.size(1))
    CHECK_SAME(nrhs,X.size(0))
    X.copy(B);
    F77NAME(spbtrf)('L',n,ns,BF->addr(),nb,info);
    CHECK_TEST(info==0);
    for (int j=0;j<nrhs;j++) { // spbtrs loop 20:
      F77NAME(stbsv)('L','N','N',n,ns,BF->addr(),nb,X.addr(0,j),nrhs);
      F77NAME(stbsv)('L','T','N',n,ns,BF->addr(),nb,X.addr(0,j),nrhs);
    }
  }
  delete BF; BF=0;
}

template class SymmetricPositiveBandMatrix<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "GaussianFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//template<> char GaussianFactorization<float,float>::norm='1';
template<>
GaussianFactorization<float,float,SquareMatrix<float,float> >::
GaussianFactorization(const SquareMatrix<float,float>& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),piv_op(po),equed('N'),
colcnd(HUGE_VAL),rowcnd(HUGE_VAL) {
  TRACER_CALL(t,"GaussianFactorization::GaussianFactorization");
  int n=A.size(0);
  LU=OPERATOR_NEW SquareMatrix<float,float>(n);
  LU->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    r=OPERATOR_NEW Vector<float,float>(n);
    c=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    F77NAME(sgeequ)(n,n,LU->addr(),n,r->addr(),c->addr(),rowcnd,colcnd,
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
      F77NAME(slaqge)(n,n,LU->addr(),n,r->addr(),c->addr(),
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
	float alpha=float_one_/(*LU)(k,k);
	int nrows=n-k-1;
	if (nrows>0) {
	  F77NAME(sscal)(nrows,alpha,LU->addr(k+1,k),1);
	  for (int j=k+1;j<n;j++) {
	    alpha=-(*LU)(k,j);
	    F77NAME(saxpy)(nrows,alpha,LU->addr(k+1,k),1,
			   LU->addr(k+1,j),1);
	  }
	}
      }
      break;
    }
    case Factorization::PIVOT_ROWS_AND_COLUMNS: {
      TRACER_CALL(t,"GaussianFactorization::GaussianFactorization sgetc2");
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      jpiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(sgetc2)(n,LU->addr(),n,ipiv,jpiv,info);
      break;
    }
    default: {
      TRACER_CALL(t,"GaussianFactorization::GaussianFactorization sgetrf");
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(sgetrf)(n,n,LU->addr(),n,ipiv,info);
      break;
    }
  }
  float *work=0;
  if (info>0) { // see sgesvx
    rpvgrw=F77NAME(slantr)('M','U','N',info,info,LU->addr(),n,work);
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
    rpvgrw=F77NAME(slantr)('M','U','N',n,n,LU->addr(),n,work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
  }
}

template<> void
GaussianFactorization<float,float,SquareMatrix<float,float> >::solve(
const Vector<float,float> &b,Vector<float,float> &x,
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
      const float *ri=r->addr();
      float *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
    if (ipiv!=0) F77NAME(slaswp)(1,x.addr(),n,1,n-1,ipiv,1);
    F77NAME(strsm)('L','L','N','U',n,1,float_one_,LU->addr(),n,
      x.addr(),n);
    F77NAME(strsm)('L','U','N','N',n,1,float_one_,LU->addr(),n,
      x.addr(),n);
    if (jpiv!=0) F77NAME(slaswp)(1,x.addr(),n,1,n-1,jpiv,-1);
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ci=c->addr();
      float *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
  } else { // A^T X = B ==> U^T L^T Q^T R^{-1} X = P^T C^T B
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ci=c->addr();
      float *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
    if (jpiv!=0) F77NAME(slaswp)(1,x.addr(),n,1,n-1,jpiv,1);
    F77NAME(strsm)('L','U','T','N',n,1,float_one_,LU->addr(),n,
      x.addr(),n);
    F77NAME(strsm)('L','L','T','U',n,1,float_one_,LU->addr(),n,
      x.addr(),n);
    if (ipiv!=0) F77NAME(slaswp)(1,x.addr(),n,1,n-1,ipiv,-1);
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ri=r->addr();
      float *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
  }
}

template<> void
GaussianFactorization<float,float,SquareMatrix<float,float> >::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,
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
          const float *ri=r->addr();
          float *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
      if (ipiv!=0) F77NAME(slaswp)(n,X.addr(),m,1,m-1,ipiv,1);
      F77NAME(strsm)('L','L','N','U',m,n,float_one_,
                     LU->addr(),m,X.addr(),m);
      F77NAME(strsm)('L','U','N','N',m,n,float_one_,
                     LU->addr(),m,X.addr(),m);
      if (jpiv!=0) F77NAME(slaswp)(n,X.addr(),m,1,m-1,jpiv,-1);
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ci=c->addr();
          float *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
    } else { // A^T X = B ==> U^T L^T Q^T R^{-1} X = P^T C^T B
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ci=c->addr();
          float *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
      if (jpiv!=0) F77NAME(slaswp)(n,X.addr(),m,1,m-1,jpiv,1);
      F77NAME(strsm)('L','U','T','N',m,n,float_one_,
                     LU->addr(),m,X.addr(),m);
      F77NAME(strsm)('L','L','T','U',m,n,float_one_,
                     LU->addr(),m,X.addr(),m);
      if (ipiv!=0) F77NAME(slaswp)(n,X.addr(),m,1,m-1,ipiv,-1);
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ri=r->addr();
          float *Xij=X.addr(0,j);
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
        for (int j=0;j<n;j++) F77NAME(sscal)(m,(*c)[j],X.addr(0,j),1);
      }
      if (jpiv!=0) {
        for (int j=n-2;j>=0;j--) {
          cout << "jpiv[" << j << "] = " << jpiv[j] << endl;
          if (j!=jpiv[j]-1) {
            F77NAME(sswap)(m,X.addr(0,j),1,X.addr(0,jpiv[j]-1),1);
          }
        }
      }
      F77NAME(strsm)('R','U','N','N',m,n,float_one_,
                     LU->addr(),n,X.addr(),m);
      F77NAME(strsm)('R','L','N','U',m,n,float_one_,
                     LU->addr(),n,X.addr(),m);
      if (ipiv!=0) {
        for (int i=0;i<n-1;i++) {
          if (i!=ipiv[i]-1) {
            F77NAME(sswap)(m,X.addr(0,i),1,X.addr(0,ipiv[i]-1),1);
          }
        }
      }
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(sscal)(m,(*r)[j],X.addr(0,j),1);
      }
    } else { // X A^T = B ==> X C^{-1} P U^T L^T = B R Q
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(sscal)(m,(*r)[j],X.addr(0,j),1);
      }
      if (ipiv!=0) {
        for (int i=n-2;i>=0;i--) {
          if (i!=ipiv[i]-1) {
            F77NAME(sswap)(m,X.addr(0,i),1,X.addr(0,ipiv[i]-1),1);
          }
        }
      }
      F77NAME(strsm)('R','L','T','U',m,n,float_one_,
                     LU->addr(),n,X.addr(),m);
      F77NAME(strsm)('R','U','T','N',m,n,float_one_,
                     LU->addr(),n,X.addr(),m);
      if (jpiv!=0) {
        for (int j=0;j<n-1;j++) {
          if (j!=jpiv[j]-1) {
            F77NAME(sswap)(m,X.addr(0,j),1,X.addr(0,jpiv[j]-1),1);
          }
        }
      }
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(sscal)(m,(*c)[j],X.addr(0,j),1);
      }
    }
  }
}

template<> float
GaussianFactorization<float,float,SquareMatrix<float,float> >::
reciprocalConditionNumber(Factorization::CONDITION_NUMBER_NORM cnn) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  float anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,4*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  F77NAME(sgecon)(norm,n,LU->addr(),n,anorm,rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void
GaussianFactorization<float,float,SquareMatrix<float,float> >::
improve(const Vector<float,float> &b,Vector<float,float> &x,
float &berr,float &ferr,Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(n+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,float> residual(n);
  Vector<float,float> work(n);
  Vector<float,float> v(n);
  char trans=(to==Factorization::TRANSPOSE ? 'T' : 'N');
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(sgemv)(trans,n,n,float_mone_,A_original->addr(),n,
      x.addr(),1,float_one_,residual.addr(),1);
    work.copy(b);
    F77_NAME(sla_geamv)(F77NAME(ilatrans)(trans),n,n,float_one_,
      A_original->addr(),n,x.addr(),1,float_one_,work.addr(),1);

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
    F77NAME(saxpy)(n,float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  float *residuali=residual.addr();
  float *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(slacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual,to);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual,to);
    }
  }
  int i=F77NAME(isamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void
GaussianFactorization<float,float,SquareMatrix<float,float> >::
improve(const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &berr,Vector<float,float> &ferr,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
//compare to sgesvx.f
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
  Vector<float,float> x(k);
  Vector<float,float> rhs(k);
  Vector<float,float> residual(k);
  Vector<float,float> work(k);
  Vector<float,float> v(k);
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
      F77NAME(scopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(scopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(scopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(scopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(sgemv)(trans,k,k,float_mone_,A_original->addr(),k,
        x.addr(),1,float_one_,residual.addr(),1);
      work.copy(rhs);
      F77_NAME(sla_geamv)(F77NAME(ilatrans)(trans),k,k,float_one_,
        A_original->addr(),k,x.addr(),1,float_one_,work.addr(),1);

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
      F77NAME(saxpy)(m,float_one_,residual.addr(),1,x.addr(),1);
      lstres=s;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(scopy)(k,x.addr(),1,X.addr(0,j),1);
      } else F77NAME(scopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(scopy)(k,x.addr(),1,X.addr(j,0),m);
      } else F77NAME(scopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    float *residuali=residual.addr();
    float *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(slacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual,lto);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual,lto);
      }
    }
    int i=F77NAME(isamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr[j]/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<>
GaussianFactorization<float,float,TridiagonalMatrix<float,float> >::
GaussianFactorization(const TridiagonalMatrix<float,float>& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),
piv_op(Factorization::PIVOT_ROWS),equ_op(Factorization::NO_EQUILIBRATION),
equed('N'),colcnd(HUGE_VAL),rowcnd(HUGE_VAL) {
  int n=A.size(0);
//LU=OPERATOR_NEW SquareMatrix<float,float>(n);
//LU->copy(A);

  anormi=A.normInfinity();
  anormm=A.normMaxEntry();
  anormo=A.normOne();
  rpvgrw=float_one_; // not computed
}

template<> void
GaussianFactorization<float,float,TridiagonalMatrix<float,float> >::
solve(const Vector<float,float> &b,Vector<float,float> &x,
Factorization::TRANSPOSE_OPTION to) {
//constructor did not factor
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  LU=OPERATOR_NEW TridiagonalMatrix<float,float>(n);
  LU->copy(*A_original);
  if (&x!=&b) x.copy(b);
  int info;
  if (to==Factorization::NO_TRANSPOSE) {
    if (piv_op==Factorization::PIVOT_ROWS) {
      F77NAME(sgtsv)(n,1,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),x.addr(),n,info);
    } else {
      F77NAME(sgtsvnp)(n,1,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),x.addr(),n,info);
    }
  } else {
    if (piv_op==Factorization::PIVOT_ROWS) {
      F77NAME(sgtsv)(n,1,LU->upperDiagonalAddr(),LU->diagonalAddr(),
        LU->lowerDiagonalAddr(),x.addr(),n,info);
    } else {
      F77NAME(sgtsvnp)(n,1,LU->upperDiagonalAddr(),LU->diagonalAddr(),
        LU->lowerDiagonalAddr(),x.addr(),n,info);
    }
  }
  CHECK_TEST(info==0);
  delete LU; LU=0;
}

template<> void
GaussianFactorization<float,float,TridiagonalMatrix<float,float> >::
solve(const Matrix<float,float> &B,Matrix<float,float> &X,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
//constructor did not factor
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  LU=OPERATOR_NEW TridiagonalMatrix<float,float>(A_original->size(0));
  LU->copy(*A_original);
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) { // A X = B
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(sgtsv)(m,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(sgtsvnp)(m,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      }
    } else { // A^T X = B
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(sgtsv)(m,n,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(sgtsvnp)(m,n,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      }
    }
  } else {
    CHECK_SAME(n,LU->size(0));
    if (to==Factorization::NO_TRANSPOSE) { // X A = B ==> A^T X^T = B^T
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(sgtsvr)(n,m,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(sgtsvrnp)(n,m,LU->upperDiagonalAddr(),LU->diagonalAddr(),
          LU->lowerDiagonalAddr(),X.addr(),m,info);
      }
    } else { // X A^T = B ==> A X^T = B^T
      if (piv_op==Factorization::PIVOT_ROWS) {
        F77NAME(sgtsvr)(n,m,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      } else {
        F77NAME(sgtsvrnp)(n,m,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
          LU->upperDiagonalAddr(),X.addr(),m,info);
      }
    }
  }
  CHECK_TEST(info==0);
  delete LU; LU=0;
}

template<> float
GaussianFactorization<float,float,TridiagonalMatrix<float,float> >::
reciprocalConditionNumber(Factorization::CONDITION_NUMBER_NORM cnn)
{
  LU=OPERATOR_NEW TridiagonalMatrix<float,float>(A_original->size(0));
  LU->copy(*A_original);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  float anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<float,float> *u2=OPERATOR_NEW Vector<float,float>(n-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(sgttrf)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    F77NAME(sgtcon)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,iwork,
      info);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(sgttrfnp)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    F77NAME(sgtconnp)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,iwork,info);
  }
  CHECK_SAME(info,0)
  delete LU; LU=0;
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void
GaussianFactorization<float,float,TridiagonalMatrix<float,float> >::
improve(const Vector<float,float> &b,Vector<float,float> &x,
float &berr,float &ferr,Factorization::TRANSPOSE_OPTION to) {
  LU=OPERATOR_NEW TridiagonalMatrix<float,float>(A_original->size(0));
  LU->copy(*A_original);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  char norm=(to==Factorization::NO_TRANSPOSE ? 'O' : 'I');
  float anorm=(to==Factorization::NO_TRANSPOSE ? anormo : anormi);
  char trans=(to==Factorization::NO_TRANSPOSE ? 'N' : 'T');
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info;
  float rcond;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<float,float> *u2=OPERATOR_NEW Vector<float,float>(n-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,n);
    F77NAME(sgttrf)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    CHECK_TEST(info==0);
    F77NAME(sgtcon)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,iwork,
      info);
    CHECK_TEST(info==0);
    F77NAME(sgtrfs)(trans,n,1,A_original->lowerDiagonalAddr(),
      A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
      LU->lowerDiagonalAddr(),LU->diagonalAddr(),LU->upperDiagonalAddr(),
      u2->addr(),ipiv,b.addr(),n,x.addr(),n,&ferr,&berr,work,iwork,info);
    CHECK_TEST(info==0);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(sgttrfnp)(n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    CHECK_TEST(info==0);
    F77NAME(sgtconnp)(norm,n,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,iwork,info);
    CHECK_TEST(info==0);
    F77NAME(sgtrfsnp)(trans,n,1,A_original->lowerDiagonalAddr(),
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
GaussianFactorization<float,float,TridiagonalMatrix<float,float> >::
improve(const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &berr,Vector<float,float> &ferr,
Factorization::TRANSPOSE_OPTION to,Factorization::SIDE_OPTION so) {
  LU=OPERATOR_NEW TridiagonalMatrix<float,float>(A_original->size(0));
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
  char trans=(to==Factorization::NO_TRANSPOSE ? 'N' : 'T');
  float *work=OPERATOR_NEW_BRACKET(float,3*k);
  int *iwork=OPERATOR_NEW_BRACKET(int,k);
  int info;
  float rcond;
  if (piv_op==Factorization::PIVOT_ROWS) {
    Vector<float,float> *u2=OPERATOR_NEW Vector<float,float>(k-2);
    int *ipiv=OPERATOR_NEW_BRACKET(int,k);
    F77NAME(sgttrf)(k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,info);
    CHECK_TEST(info==0);
    F77NAME(sgtcon)(norm,k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),u2->addr(),ipiv,anorm,rcond,work,iwork,
      info);
    CHECK_TEST(info==0);
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(sgtrfs)(trans,k,n,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),u2->addr(),ipiv,B.addr(),m,X.addr(),m,
        ferr.addr(),berr.addr(),work,iwork,info);
    } else {
      F77NAME(sgtrfsr)(trans,k,m,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),u2->addr(),ipiv,B.addr(),m,X.addr(),m,
        ferr.addr(),berr.addr(),work,iwork,info);
    }
    CHECK_TEST(info==0);
    delete u2; u2=0;
    delete [] ipiv;ipiv=0;
  } else {
    F77NAME(sgttrfnp)(k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),info);
    CHECK_TEST(info==0);
    F77NAME(sgtconnp)(norm,k,LU->lowerDiagonalAddr(),LU->diagonalAddr(),
      LU->upperDiagonalAddr(),anorm,rcond,work,iwork,info);
    CHECK_TEST(info==0);
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(sgtrfsnp)(trans,k,n,A_original->lowerDiagonalAddr(),
        A_original->diagonalAddr(),A_original->upperDiagonalAddr(),
        LU->lowerDiagonalAddr(),LU->diagonalAddr(),
        LU->upperDiagonalAddr(),B.addr(),m,X.addr(),m,ferr.addr(),
        berr.addr(),work,iwork,info);
    } else {
      F77NAME(sgtrfsrnp)(trans,k,m,A_original->lowerDiagonalAddr(),
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
GaussianFactorization<float,float,BandMatrix<float,float> >::
GaussianFactorization(const BandMatrix<float,float>& A,
Factorization::PIVOT_OPTION po,Factorization::EQUILIBRATE_OPTION eo) :
A_original(&A),LU(0),r(0),c(0),ipiv(0),jpiv(0),piv_op(po),equ_op(eo),
equed('N'),colcnd(HUGE_VAL),rowcnd(HUGE_VAL) {
  int n=A.size(0),nsub=A.subDiags(),nsup=A.supDiags();
  LU=OPERATOR_NEW BandMatrix<float,float>(n,nsub,
    (po==Factorization::NO_PIVOTING ? nsup : nsub+nsup),float_zero_);
  for (int j=0;j<n;j++) {
    int ibeg=max(0,j-nsup);
    int iend=min(n-1,j+nsub);
    F77NAME(scopy)(iend-ibeg+1,A.addr(ibeg,j),1,LU->addr(ibeg,j),1);
  }

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    r=OPERATOR_NEW Vector<float,float>(n);
    c=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    F77NAME(sgbequ)(n,n,LU->subDiags(),LU->supDiags(),LU->addr(),
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
      F77NAME(slaqgb)(n,n,nsub,nsub+nsup,LU->addr(),LU->bands(),r->addr(),
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
      F77NAME(sgbtf2np)(n,n,nsub,nsup,LU->addr(),LU->bands(),info);
      break;
    }
    default: {
      ipiv=OPERATOR_NEW_BRACKET(int,n);
      F77NAME(sgbtrf)(n,n,nsub,nsup,LU->addr(),LU->bands(),ipiv,info);
    }
  }
  float *work=0;
  if (info>0) { // see dgbsvx
    rpvgrw=F77NAME(slantb)('M','U','N',info,min(info-1,LU->supDiags()),
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
    rpvgrw=F77NAME(slantb)('M','U','N',n,LU->supDiags(),LU->addr(),
      LU->bands(),work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
  }
}

template<> void
GaussianFactorization<float,float,BandMatrix<float,float> >::
solve(const Vector<float,float> &b,Vector<float,float> &x,
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
      const float *ri=r->addr();
      float *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
    if (ipiv!=0) {
      F77NAME(sgbtrs)('N',n,A_original->subDiags(),A_original->supDiags(),
        1,LU->addr(),LU->bands(),ipiv,x.addr(),n,info);
    } else {
      F77NAME(sgbtrsnp)('N',n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),x.addr(),n,info);
    }
    CHECK_TEST(info==0);
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ci=c->addr();
      float *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
  } else { // A^T X = B ==> U^T L^T Q^T R^{-1} X = C^T B
    if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ci=c->addr();
      float *xi=x.addr();
      for (int i=0;i<n;i++,ci++,xi++) *xi*=*ci;
    }
    if (ipiv!=0) {
      F77NAME(sgbtrs)('T',n,A_original->subDiags(),A_original->supDiags(),
        1,LU->addr(),LU->bands(),ipiv,x.addr(),n,info);
    } else {
      F77NAME(sgbtrsnp)('T',n,A_original->subDiags(),
        A_original->supDiags(),1,LU->addr(),LU->bands(),x.addr(),n,info);
    }
    CHECK_TEST(info==0);
    if (equ_op==Factorization::EQUILIBRATE_ROWS ||
    equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      const float *ri=r->addr();
      float *xi=x.addr();
      for (int i=0;i<n;i++,ri++,xi++) *xi*=*ri;
    }
  }
}

template<> void
GaussianFactorization<float,float,BandMatrix<float,float> >::
solve(const Matrix<float,float> &B,Matrix<float,float> &X,
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
          const float *ri=r->addr();
          float *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ri++,Xij++) *Xij*=*ri;
        }
      }
      if (ipiv!=0) {
        F77NAME(sgbtrs)('N',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),ipiv,
          X.addr(),m,info);
      } else {
        F77NAME(sgbtrsnp)('N',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),X.addr(),m,
          info);
      }
      CHECK_TEST(info==0);
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ci=c->addr();
          float *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
    } else { // A^T X = B ==> U^T L^T Q^T R^{-1} X = C^T B
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ci=c->addr();
          float *Xij=X.addr(0,j);
          for (int i=0;i<m;i++,ci++,Xij++) *Xij*=*ci;
        }
      }
      if (ipiv!=0) {
        F77NAME(sgbtrs)('T',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),ipiv,
          X.addr(),m,info);
      } else {
        F77NAME(sgbtrsnp)('T',m,A_original->subDiags(),
          A_original->supDiags(),n,LU->addr(),LU->bands(),X.addr(),m,
          info);
      }
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) {
          const float *ri=r->addr();
          float *Xij=X.addr(0,j);
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
        for (int j=0;j<n;j++) F77NAME(sscal)(m,(*c)[j],X.addr(0,j),1);
      }
      Vector<float,float> *x=OPERATOR_NEW Vector<float,float>(n);
      for (int i=0;i<m;i++) {
        F77NAME(scopy)(n,X.addr(i,0),m,x->addr(),1);
        if (ipiv!=0) {
          F77NAME(sgbtrs)('T',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
            x->addr(),n,info);
        } else {
          F77NAME(sgbtrsnp)('T',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),x->addr(),n,
            info);
        }
        F77NAME(scopy)(n,x->addr(),1,X.addr(i,0),m);
      }
      delete x; x=0;
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(sscal)(m,(*r)[j],X.addr(0,j),1);
      }
    } else {
      // X A^T = B ==> X C^{-1} U^T L^T = B R Q
      //           ==> L U C^{-1} X^T = Q^T R B^T
      if (equ_op==Factorization::EQUILIBRATE_ROWS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(sscal)(m,(*r)[j],X.addr(0,j),1);
      }
      Vector<float,float> *x=OPERATOR_NEW Vector<float,float>(n);
      for (int i=0;i<m;i++) {
        F77NAME(scopy)(n,X.addr(i,0),m,x->addr(),1);
        if (ipiv!=0) {
          F77NAME(sgbtrs)('N',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),ipiv,
            x->addr(),n,info);
        } else {
          F77NAME(sgbtrsnp)('N',n,A_original->subDiags(),
            A_original->supDiags(),1,LU->addr(),LU->bands(),x->addr(),n,
            info);
        }
        F77NAME(scopy)(n,x->addr(),1,X.addr(i,0),m);
      }
      delete x; x=0;
      if (equ_op==Factorization::EQUILIBRATE_COLUMNS ||
      equ_op==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
        for (int j=0;j<n;j++) F77NAME(sscal)(m,(*c)[j],X.addr(0,j),1);
      }
    }
  }
}

template<> float
GaussianFactorization<float,float,BandMatrix<float,float> >::
reciprocalConditionNumber(Factorization::CONDITION_NUMBER_NORM cnn) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  char norm=(cnn==Factorization::INFINITY_NORM ? 'I' : 'O');
  float anorm=(cnn==Factorization::INFINITY_NORM ? anormi : anormo);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  if (piv_op==Factorization::PIVOT_ROWS) {
    F77NAME(sgbcon)(norm,n,A_original->subDiags(),A_original->supDiags(),
      LU->addr(),LU->bands(),ipiv,anorm,rcond,work,iwork,info);
  } else {
    F77NAME(sgbconnp)(norm,n,A_original->subDiags(),
      A_original->supDiags(),LU->addr(),LU->bands(),anorm,rcond,work,
      iwork,info);
  }
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void
GaussianFactorization<float,float,BandMatrix<float,float> >::improve(
const Vector<float,float> &b,Vector<float,float> &x,float &berr,
float &ferr,Factorization::TRANSPOSE_OPTION to) {
  CHECK_POINTER(LU);
  int n=LU->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(A_original->bands()+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,float> residual(n);
  Vector<float,float> work(n);
  char trans=(to==Factorization::TRANSPOSE ? 'T' : 'N');
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(sgbmv)(trans,n,n,A_original->subDiags(),
      A_original->supDiags(),float_mone_,A_original->addr(),
      A_original->bands(),x.addr(),1,float_one_,residual.addr(),1);
    work.copy(b);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(sgbamv)(trans,n,n,A_original->subDiags(),
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
    F77NAME(saxpy)(n,float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }
  float *residuali=residual.addr();
  float *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  Vector<float,float> v(n);
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(slacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual,to);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual,to);
    }
  }
  int i=F77NAME(isamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void
GaussianFactorization<float,float,BandMatrix<float,float> >::
improve(const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &berr,Vector<float,float> &ferr,
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
  float nz=static_cast<float>(A_original->bands()+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,float> x(k);
  Vector<float,float> rhs(k);
  Vector<float,float> residual(k);
  Vector<float,float> work(k);
  Vector<float,float> v(k);
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
      F77NAME(scopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(scopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(scopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(scopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(sgbmv)(trans,k,k,A_original->subDiags(),
        A_original->supDiags(),float_mone_,
        A_original->addr(),A_original->bands(),x.addr(),1,float_one_,
        residual.addr(),1);
      for (int i=0;i<k;i++) work[i]=abs(rhs[i]);
      F77NAME(sgbamv)(trans,k,k,A_original->subDiags(),
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
      solve(residual,residual,to);
      F77NAME(saxpy)(k,float_one_,residual.addr(),1,x.addr(),1);
      lstres=s;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(scopy)(k,x.addr(),1,X.addr(0,j),1);
      } else F77NAME(scopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      if (to==Factorization::NO_TRANSPOSE) {
        F77NAME(scopy)(k,x.addr(),1,X.addr(j,0),m);
      } else F77NAME(scopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    float *residuali=residual.addr();
    float *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(slacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],
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
    int i=F77NAME(isamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template class
  GaussianFactorization<float,float,SquareMatrix<float,float> >;
template void testGaussianFactorization(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "CholeskyFactorization.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> CholeskyFactorization<float,float,
SymmetricPositiveMatrix<float,float> >::CholeskyFactorization(
const SymmetricPositiveMatrix<float,float>& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equed('N'),scond(HUGE_VAL) {
  int n=A.size(0);
  L=OPERATOR_NEW SymmetricPositiveMatrix<float,float>(n);
  L->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    F77NAME(spoequb)(n,L->addr(),n,s->addr(),scond,amax,info);
    equed='N';
    F77NAME(slaqsy)('L',n,L->addr(),n,s->addr(),scond,amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=L->normInfinity();
  anormm=L->normMaxEntry();
  anormo=L->normOne();

  int info=0;
  F77NAME(spotrf)('L',n,L->addr(),n,info);
  float *work=OPERATOR_NEW_BRACKET(float,4*n);
  if (info!=0) { // see sposvxx
    rpvgrw=F77_NAME(sla_porpvgrw)('L',info,A_original->addr(),n,
      L->addr(),n, work);
    delete L; L=0;
    if (s!=0) delete s; s=0;
  } else {
    rpvgrw=F77_NAME(sla_porpvgrw)('L',n,A_original->addr(),n,
      L->addr(),n, work);
  }
  delete [] work; work=0;
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveMatrix<float,float> >::solve(
const Vector<float,float> &b,Vector<float,float> &x) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^T
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,L->size(0));
  if (&x!=&b) x.copy(b);
//A x = b ==> L L^T S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    float *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  int info;
  F77NAME(spotrs)('L',n,1,L->addr(),n,x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    float *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveMatrix<float,float> >::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,
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
        const float *si=s->addr();
        float *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    int info;
    F77NAME(spotrs)('L',m,n,L->addr(),m,X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const float *si=s->addr();
        float *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,L->size(0));
//  X A = B ==> X S^{-1} L L^T = B S
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(sscal)(m,(*s)[j],X.addr(0,j),1);
    }
    F77NAME(strsm)('R','L','T','N',m,n,float_one_,L->addr(),n,
      X.addr(),m);
    F77NAME(strsm)('R','L','N','N',m,n,float_one_,L->addr(),n,
      X.addr(),m);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(sscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> float CholeskyFactorization<float,float,
SymmetricPositiveMatrix<float,float> >::reciprocalConditionNumber() {
  CHECK_POINTER(L);
  int n=L->size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  F77NAME(spocon)('L',n,L->addr(),n,anormo,rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveMatrix<float,float> >::improve(const Vector<float,
float> &b,Vector<float,float> &x,float &berr,float &ferr) {
  CHECK_POINTER(L);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(n+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,float> residual(n);
  Vector<float,float> work(n);
  Vector<float,float> v(n);
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(ssymv)('L',n,float_mone_,A_original->addr(),n,
      x.addr(),1,float_one_,residual.addr(),1);
    work.copy(b);
    F77_NAME(sla_syamv)(F77NAME(ilauplo)('L'),n,float_one_,
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
    F77NAME(saxpy)(n,float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  float *residuali=residual.addr();
  float *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(slacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(isamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveMatrix<float,float> >::improve(
const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &berr,Vector<float,float> &ferr,
Factorization::SIDE_OPTION so) {
//compare to sposvx.f
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
  Vector<float,float> x(k);
  Vector<float,float> rhs(k);
  Vector<float,float> residual(k);
  Vector<float,float> work(k);
  Vector<float,float> v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(scopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(scopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(scopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(scopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(ssymv)('L',k,float_mone_,A_original->addr(),k,
        x.addr(),1,float_one_,residual.addr(),1);
      work.copy(rhs);
      F77_NAME(sla_syamv)(F77NAME(ilauplo)('L'),k,float_one_,
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
      F77NAME(saxpy)(k,float_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(scopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      F77NAME(scopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    float *residuali=residual.addr();
    float *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(slacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(isamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr[j]/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> CholeskyFactorization<float,float,
SymmetricPositiveTridiagonalMatrix<float,float> >::
CholeskyFactorization(
const SymmetricPositiveTridiagonalMatrix<float,float>& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equ_op(Factorization::NO_EQUILIBRATION),equed('N'),scond(HUGE_VAL) {
  int n=A.size(0);
  anormi=A.normInfinity();
  anormm=A.normMaxEntry();
  anormo=A.normOne();
  rpvgrw=float_one_; // not computed
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveTridiagonalMatrix<float,float> >::solve(
const Vector<float,float> &b,Vector<float,float> &x) {
//constructor did not factor
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<float,float>(n);
  L->copy(*A_original);
  if (&x!=&b) x.copy(b);
  int info;
  F77NAME(sptsv)(n,1,L->diagonalAddr(0),L->lowerDiagonalAddr(0),
    x.addr(),n,info);
  CHECK_TEST(info==0);
  delete L; L=0;
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveTridiagonalMatrix<float,float> >::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,
Factorization::SIDE_OPTION so) {
//constructor did not factor
  int m=B.size(0),n=B.size(1);
  CHECK_SAME(m,X.size(0))
  CHECK_SAME(n,X.size(1))
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<float,float>(
    A_original->size(0));
  L->copy(*A_original);
  if (&X!=&B) X.copy(B);
  int info;
  if (so==Factorization::LEFT_SIDE) {
    CHECK_SAME(m,L->size(0));
    F77NAME(sptsv)(m,n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),
      X.addr(),m,info);
  } else {
    CHECK_SAME(n,L->size(0));
    F77NAME(spttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
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

template<> float CholeskyFactorization<float,float,
SymmetricPositiveTridiagonalMatrix<float,float> >::
reciprocalConditionNumber() {
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<float,float>(
    A_original->size(0));
  L->copy(*A_original);
  int n=L->size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,n);
  int info=0;
  F77NAME(spttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  F77NAME(sptcon)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,work,info);
  CHECK_SAME(info,0)
  delete L; L=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveTridiagonalMatrix<float,float> >::improve(
const Vector<float,float> &b,Vector<float,float> &x,
float &berr,float &ferr) {
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<float,float>(
    A_original->size(0));
  L->copy(*A_original);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float *work=OPERATOR_NEW_BRACKET(float,2*n);
  int info;
  float rcond;
  F77NAME(spttrf)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  CHECK_TEST(info==0);
  F77NAME(sptcon)(n,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,work,info);
  CHECK_TEST(info==0);
  F77NAME(sptrfs)(n,1,A_original->diagonalAddr(0),
    A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
    L->lowerDiagonalAddr(0),b.addr(),n,x.addr(),n,&ferr,&berr,work,info);
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete L; L=0;
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveTridiagonalMatrix<float,float> >::improve(
const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &berr,Vector<float,float> &ferr,
Factorization::SIDE_OPTION so) {
  L=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<float,float>(
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
  float *work=OPERATOR_NEW_BRACKET(float,2*k);
  int info;
  float rcond;
  F77NAME(spttrf)(k,L->diagonalAddr(0),L->lowerDiagonalAddr(0),info);
  CHECK_TEST(info==0);
  F77NAME(sptcon)(k,L->diagonalAddr(0),L->lowerDiagonalAddr(0),anormo,
    rcond,work,info);
  CHECK_TEST(info==0);
  if (so==Factorization::LEFT_SIDE) {
    F77NAME(sptrfs)(k,n,A_original->diagonalAddr(0),
      A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
      L->lowerDiagonalAddr(0),B.addr(),m,X.addr(),m,ferr.addr(),
      berr.addr(),work,info);
  } else {
    F77NAME(sptrfsr)(k,m,A_original->diagonalAddr(0),
      A_original->lowerDiagonalAddr(0),L->diagonalAddr(0),
      L->lowerDiagonalAddr(0),B.addr(),m,X.addr(),m,ferr.addr(),
      berr.addr(),work,info);
  }
  CHECK_TEST(info==0);
  delete [] work; work=0;
  delete L; L=0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template<> CholeskyFactorization<float,float,
SymmetricPositiveBandMatrix<float,float> >::CholeskyFactorization(
const SymmetricPositiveBandMatrix<float,float>& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),L(0),s(0),
equ_op(eo),equed('N'),scond(HUGE_VAL) {
  int n=A.size(0),nsub=A.subDiags();
  L=OPERATOR_NEW
    SymmetricPositiveBandMatrix<float,float>(n,nsub,float_zero_);
  for (int j=0;j<n;j++) {
    F77NAME(scopy)(min(n-j,nsub+1),A.addr(j,j),1,L->addr(j,j),1);
  }

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    F77NAME(spbequ)('L',n,L->subDiags(),L->addr(),L->bands(),s->addr(),
      scond,amax,info);
    equed='N';
    F77NAME(slaqsb)('L',n,nsub,L->addr(),L->bands(),s->addr(),scond,
      amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
  }
  anormi=L->normInfinity();
  anormm=L->normMaxEntry();
  anormo=L->normOne();

  int info=0;
  F77NAME(spbtrf)('L',n,nsub,L->addr(),L->bands(),info);
  float *work=0;
  if (info>0) { // see dgbsvx
    rpvgrw=F77NAME(slansb)('M','L',info,min(info-1,L->subDiags()),
      L->addr(),L->bands(),work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
    cerr << "zero pivot in CholeskyFactorization::CholeskyFactorization"
       << "\n zero pivot number,reciprocal pivot growth = " << info
       << " " << rpvgrw << endl;
    if (L) delete L; L=0;
    if (s) delete s; s=0;
    A_original=0;
  } else {
    rpvgrw=F77NAME(slansb)('M','L',n,L->subDiags(),L->addr(),
      L->bands(),work);
    rpvgrw=(rpvgrw==float_zero_ ? float_one_ : anormm / rpvgrw);
  }
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveBandMatrix<float,float> >::solve(
const Vector<float,float> &b,Vector<float,float> &x) {
  CHECK_POINTER(L);
//constructor factored S A S = L L^T
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,A_original->size(0));
  if (&x!=&b) x.copy(b);
  int info;
//A x = b ==> L l^T S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    float *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  F77NAME(spbtrs)('L',n,L->subDiags(),1,L->addr(),L->bands(),
    x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    float *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveBandMatrix<float,float> >::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,
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
        const float *si=s->addr();
        float *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    F77NAME(spbtrs)('L',m,L->subDiags(),n,L->addr(),L->bands(),
      X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const float *si=s->addr();
        float *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,L->size(0));
//  X A = B ==> X R^{-1} Q L U = B C
//          ==> U^T L^T Q^T R^{-1} X^T = C^T B^T
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(sscal)(m,(*s)[j],X.addr(0,j),1);
    }
    Vector<float,float> *x=OPERATOR_NEW Vector<float,float>(n);
    for (int i=0;i<m;i++) {
      F77NAME(scopy)(n,X.addr(i,0),m,x->addr(),1);
      F77NAME(spbtrs)('L',n,L->subDiags(),1,L->addr(),L->bands(),
        x->addr(),n,info);
      F77NAME(scopy)(n,x->addr(),1,X.addr(i,0),m);
    }
    delete x; x=0;
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(sscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> float CholeskyFactorization<float,float,
SymmetricPositiveBandMatrix<float,float> >::reciprocalConditionNumber()
{
  CHECK_POINTER(L);
  int n=L->size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  F77NAME(spbcon)('L',n,L->subDiags(),L->addr(),L->bands(),anormo,
    rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveBandMatrix<float,float> >::improve(
const Vector<float,float> &b,Vector<float,float> &x,float &berr,
float &ferr) {
  CHECK_POINTER(L);
  int n=L->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(A_original->bands()+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,float> residual(n);
  Vector<float,float> work(n);
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(ssbmv)('L',n,A_original->subDiags(),float_mone_,
      A_original->addr(),A_original->bands(),x.addr(),1,
      float_one_,residual.addr(),1);
    for (int i=0;i<n;i++) work[i]=abs(b[i]);
    F77NAME(ssbamv)('L',n,A_original->subDiags(),float_one_,
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
    F77NAME(saxpy)(n,float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }
  float *residuali=residual.addr();
  float *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  Vector<float,float> v(n);
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(slacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(isamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void CholeskyFactorization<float,float,
SymmetricPositiveBandMatrix<float,float> >::improve(
const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &berr,Vector<float,float> &ferr,
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
  float nz=static_cast<float>(A_original->bands()+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,float> x(k);
  Vector<float,float> rhs(k);
  Vector<float,float> residual(k);
  Vector<float,float> work(k);
  Vector<float,float> v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(scopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(scopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(scopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(scopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(ssbmv)('L',k,A_original->subDiags(),float_mone_,
        A_original->addr(),A_original->bands(),x.addr(),1,float_one_,
        residual.addr(),1);
      for (int i=0;i<n;i++) work[i]=abs(rhs[i]);
      F77NAME(ssbamv)('L',k,A_original->subDiags(),float_one_,
        A_original->addr(),A_original->bands(),x.addr(),1,float_one_,
        work.addr(),1);
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
      F77NAME(saxpy)(k,float_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(scopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      F77NAME(scopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    float *residuali=residual.addr();
    float *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(slacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],
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
    int i=F77NAME(isamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template class CholeskyFactorization<float,float,
  SymmetricPositiveMatrix<float,float> >;
template void testCholeskyFactorization(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "MDMtFactorization.C"
template<> MDMtFactorization<float,float>::MDMtFactorization(
const SymmetricMatrix<float,float>& A,
Factorization::EQUILIBRATE_OPTION eo) : A_original(&A),MD(0),ipiv(0),s(0),
equed('N'),scond(HUGE_VAL) {
  int n=A.size(0);
  MD=OPERATOR_NEW SymmetricMatrix<float,float>(n);
  MD->copy(A);

  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    s=OPERATOR_NEW Vector<float,float>(n);
    float amax;
    int info;
    float *work=OPERATOR_NEW_BRACKET(float,3*n);
    F77NAME(ssyequb)('L',n,MD->addr(),n,s->addr(),scond,amax,work,info);
    equed='N';
    F77NAME(slaqsy)('L',n,MD->addr(),n,s->addr(),scond,amax,equed);
    equ_op=(equed=='Y' ? Factorization::EQUILIBRATE_ROWS_AND_COLUMNS :
      Factorization::NO_EQUILIBRATION);
    delete [] work; work=0;
  }
  anormi=MD->normInfinity();
  anormm=MD->normMaxEntry();
  anormo=MD->normOne();

  int info=0;
  ipiv=OPERATOR_NEW_BRACKET(int,n);
  float workq=float_zero_;
  int lwork=-1;
  F77NAME(ssytrf)('L',n,MD->addr(),n,ipiv,&workq,lwork,info);
  lwork=max(2*n,static_cast<int>(workq));
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(ssytrf)('L',n,MD->addr(),n,ipiv,work,lwork,info);
  if (info!=0) { // see ssysvxx
    rpvgrw=F77_NAME(sla_syrpvgrw)('L',n,info,A_original->addr(),n,
      MD->addr(),n,ipiv,work);
    delete MD; MD=0;
    if (s!=0) delete s; s=0;
    if (ipiv!=0) delete ipiv; ipiv=0;
  } else {
    rpvgrw=F77_NAME(sla_syrpvgrw)('L',n,n,A_original->addr(),n,
      MD->addr(),n,ipiv,work);
  }
  delete [] work; work=0;
}

template<> void MDMtFactorization<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x) {
  CHECK_POINTER(MD);
  int n=b.size();
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,MD->size(0));
  if (&x!=&b) x.copy(b);
//A x = b ==> M D M^T S^{-1} x = S b
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    float *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
  int info;
  F77NAME(ssytrs)('L',n,1,MD->addr(),n,ipiv,x.addr(),n,info);
  CHECK_TEST(info==0);
  if (equ_op!=Factorization::NO_EQUILIBRATION) {
    const float *si=s->addr();
    float *xi=x.addr();
    for (int i=0;i<n;i++,si++,xi++) *xi*=*si;
  }
}

template<> void MDMtFactorization<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,
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
        const float *si=s->addr();
        float *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
    int info;
    F77NAME(ssytrs)('L',m,n,MD->addr(),m,ipiv,X.addr(),m,info);
    CHECK_TEST(info==0);
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) {
        const float *si=s->addr();
        float *Xij=X.addr(0,j);
        for (int i=0;i<m;i++,si++,Xij++) *Xij*=*si;
      }
    }
  } else {
    CHECK_SAME(n,MD->size(0));
//  X A = B ==> X S^{-1} P M D M^T = B S P
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(sscal)(m,(*s)[j],X.addr(0,j),1);
    }
    Vector<float,float> x(n);
    for (int i=0;i<m;i++) {
      float *Xij=X.addr(i,0); 
      float *xj=x.addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *xj=*Xij;
      int info;
      F77NAME(ssytrs)('L',n,1,MD->addr(),n,ipiv,x.addr(),n,info);
      CHECK_TEST(info==0);
      Xij=X.addr(i,0); 
      xj=x.addr();
      for (int j=0;j<n;j++,Xij+=m,xj++) *Xij=*xj;
    }
    if (equ_op!=Factorization::NO_EQUILIBRATION) {
      for (int j=0;j<n;j++) F77NAME(sscal)(m,(*s)[j],X.addr(0,j),1);
    }
  }
}

template<> float
MDMtFactorization<float,float>::reciprocalConditionNumber() {
  CHECK_POINTER(MD);
  int n=MD->size(0);
  float rcond;
  float *work=OPERATOR_NEW_BRACKET(float,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  F77NAME(ssycon)('L',n,MD->addr(),n,ipiv,anormo,rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork; iwork=0;
  delete [] work; work=0;
  return rcond;
}

template<> void MDMtFactorization<float,float>::improve(
const Vector<float,float> &b,Vector<float,float> &x,float &berr,
float &ferr) {
  CHECK_POINTER(MD);
  int n=MD->size(0);
  CHECK_SAME(n,x.size())
  CHECK_SAME(n,b.size())
  float nz=static_cast<float>(n+1);
  float eps=F77NAME(slamch)('E');
  float safe1=nz*F77NAME(slamch)('S');
  float safe2=safe1/eps;
  int itmax=5;
  Vector<float,float> residual(n);
  Vector<float,float> work(n);
  Vector<float,float> v(n);
  float lstres=3.;
  int count=1;
  while (true) {
    residual.copy(b);
    F77NAME(ssymv)('L',n,float_mone_,A_original->addr(),n,
      x.addr(),1,float_one_,residual.addr(),1);
    work.copy(b);
    F77_NAME(sla_syamv)(F77NAME(ilauplo)('L'),n,float_one_,
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
    F77NAME(saxpy)(n,float_one_,residual.addr(),1,x.addr(),1);
    lstres=berr;
    count++;
  }

  float *residuali=residual.addr();
  float *worki=work.addr();
  for (int i=0;i<n;i++) {
    work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
      : abs(residual[i])+work[i]*nz*eps+safe1);
  }
  int kase=0;
  int *isgn=OPERATOR_NEW_BRACKET(int,n);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  int info;
  while (true) {
    F77NAME(slacn2)(n,v.addr(),residual.addr(),isgn,ferr,kase,isave);
    if (kase==0) break;
    if (kase==1) {
      solve(residual,residual);
      for (int i=0;i<n;i++) residual[i]*=work[i];
    } else {
      for (int i=0;i<n;i++) residual[i]*=work[i];
      solve(residual,residual);
    }
  }
  int i=F77NAME(isamax)(n,x.addr(),1)-1;
  lstres=abs(x[i]);
  if (lstres!=float_zero_) ferr/=lstres;
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template<> void MDMtFactorization<float,float>::improve(
const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &berr,Vector<float,float> &ferr,
Factorization::SIDE_OPTION so) {
//compare to ssysvx.f
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
  Vector<float,float> x(k);
  Vector<float,float> rhs(k);
  Vector<float,float> residual(k);
  Vector<float,float> work(k);
  Vector<float,float> v(k);
  int *isgn=OPERATOR_NEW_BRACKET(int,k);
  int *isave=OPERATOR_NEW_BRACKET(int,3);
  for (int j=0;j<berr.size();j++) {
    if (so==Factorization::LEFT_SIDE) { // note that k=m
      F77NAME(scopy)(k,X.addr(0,j),1,x.addr(),1);
      F77NAME(scopy)(k,B.addr(0,j),1,rhs.addr(),1);
    } else { // note that k=n
      F77NAME(scopy)(k,X.addr(j,0),m,x.addr(),1);
      F77NAME(scopy)(k,B.addr(j,0),m,rhs.addr(),1);
    }
    float lstres=3.;
    int count=1;
    while (true) {
      residual.copy(rhs);
      F77NAME(ssymv)('L',k,float_mone_,A_original->addr(),k,
        x.addr(),1,float_one_,residual.addr(),1);
      work.copy(rhs);
      F77_NAME(sla_syamv)(F77NAME(ilauplo)('L'),k,float_one_,
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
      F77NAME(saxpy)(k,float_one_,residual.addr(),1,x.addr(),1);
      lstres=berrj;
      count++;
    }
    if (so==Factorization::LEFT_SIDE) {
      F77NAME(scopy)(k,x.addr(),1,X.addr(0,j),1);
    } else {
      F77NAME(scopy)(k,x.addr(),1,X.addr(j,0),m);
    }

    float *residuali=residual.addr();
    float *worki=work.addr();
    for (int i=0;i<k;i++) {
      work[i]=(work[i]>safe2 ? abs(residual[i])+work[i]*nz*eps
        : abs(residual[i])+work[i]*nz*eps+safe1);
    }
    int kase=0;
    int info;
    while (true) {
      F77NAME(slacn2)(k,v.addr(),residual.addr(),isgn,ferr[j],kase,isave);
      if (kase==0) break;
      if (kase==1) {
        solve(residual,residual);
        for (int i=0;i<k;i++) residual[i]*=work[i];
      } else {
        for (int i=0;i<k;i++) residual[i]*=work[i];
        solve(residual,residual);
      }
    }
    int i=F77NAME(isamax)(k,x.addr(),1)-1;
    lstres=abs(x[i]);
    if (lstres!=float_zero_) ferr[j]/=lstres;
  }
  delete [] isgn; isgn=0;
  delete [] isave; isave=0;
}

template class MDMtFactorization<float,float>;
template void testMDMtFactorization(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "HouseholderQRFactorization.C"
template<>
HouseholderQRFactorization<float,float>::HouseholderQRFactorization(
const Matrix<float,float> &A) : QR(0),tau(0),A_original(&A),iascl(0),
ascl(float_one_),anrm(HUGE_VAL) {
//TRACER_CALL(t,"HouseholderQRFactorization::HouseholderQRFactorization");
  int m=A.size(0),n=A.size(1);
  QR=OPERATOR_NEW Matrix<float,float>(m,n);
  QR->copy(A);
  float w=HUGE_VAL;
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);
  anrm=F77NAME(slange)('M',m,n,A.addr(),m,&w);
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
    F77NAME(slascl)('G',0,0,anrm,ascl,m,n,QR->addr(),m,info);
    CHECK_SAME(info,0);
  }

  int lwork=-1;
  if (m>=n) {
    tau=OPERATOR_NEW Vector<float,float>(n);
    F77NAME(sgeqrf)(m,n,QR->addr(),m,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    float *work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sgeqrf)(m,n,QR->addr(),m,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work;
  } else {
    tau=OPERATOR_NEW Vector<float,float>(m);
    F77NAME(sgelqf)(m,n,QR->addr(),m,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    float *work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sgelqf)(m,n,QR->addr(),m,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work;
  }
}

template<> OrthogonalMatrix<float,float>* 
HouseholderQRFactorization<float,float>::orthogonalPart() const {
//TRACER_CALL(t,"HouseholderQRFactorization::orthogonalPart");
  int m=QR->size(0),n=QR->size(1);
  if (m>=n) {
    OrthogonalMatrix<float,float> *Q=
      OPERATOR_NEW OrthogonalMatrix<float,float>(m,m);
    Q->copyFrom('A',m,n,*QR);

    float w=HUGE_VAL;
    int lwork=-1;
    int info;
    F77NAME(sorgqr)(m,m,n,Q->addr(),m,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    float *work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sorgqr)(m,m,n,Q->addr(),m,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    return Q;
  } else {
    OrthogonalMatrix<float,float> *Q=
      OPERATOR_NEW OrthogonalMatrix<float,float>(n,n);
    Q->copyFrom('A',m,n,*QR);

    float w=HUGE_VAL;
    int lwork=-1;
    int info;
    F77NAME(sorglq)(n,n,m,Q->addr(),n,tau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    float *work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sorglq)(n,n,m,Q->addr(),n,tau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
    return Q;
  }
}

template<> float HouseholderQRFactorization<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::solve");
  int m=QR->size(0),n=QR->size(1);
  float w=HUGE_VAL;
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(slange)('M',brow,1,b.addr(),b.size(),&w);
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
  float residual_norm=HUGE_VAL;
  if (m>=n) {
    if (tr==Factorization::NO_TRANSPOSE) { // min || b - A x ||
//    if A=Q[R], r=b-Ax and 0=A'r then
//          [0]
//      [v]=Q'r=Q'b-[R]x=[y]-[Rx] and 0=[R',0]Q'r=R'v ==> v=0, w=z
//      [w]         [0]  [z] [0]
//    so compute [y]=Q'b, solve Rx=y
//               [z]
      CHECK_SAME(m,b.size())
      CHECK_SAME(n,x.size())
      Vector<float,float> xtmp(m);
      xtmp.copy(b);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    Q' b = [ y \\ z ]:
      F77NAME(sormqr)('L','T',m,1,n,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,&w,lwork,info);
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormqr)('L','T',m,1,n,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,work,lwork,info);
      CHECK_SAME(info,0)

//    solve R x = y:
      F77NAME(strsv)('U','N','N',n,QR->addr(),m,xtmp.addr(),1);
      delete [] work; work=0;
      x.copyFrom(n,xtmp);
      residual_norm=F77NAME(snrm2)(m-n,xtmp.addr(n),1);
      scllen=n;
    } else { // min || x || s.t. A' x = b
//    if A=Q[R], b=A'x and x=Ax then
//          [0]
//      b=[R',0]Q'x=[R',0][v]=R'v and [v]=Q'x=[R]s ==> w=0
//                        [w]         [w]     [0]
//    so solve R'v=b, compute x=Q[v]
//                               [0]
      CHECK_SAME(n,b.size())
      CHECK_SAME(m,x.size())

      x=float_zero_;
      x.copyFrom(m,b);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,x.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    solve R' v = b
      F77NAME(strsv)('U','T','N',n,QR->addr(),m,x.addr(),1);
//    Q [v\\0]
      F77NAME(sormqr)('L','N',m,1,n,QR->addr(),m,tau->addr(),
                      x.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormqr)('L','N',m,1,n,QR->addr(),m,tau->addr(),
                      x.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      residual_norm=float_zero_;
      scllen=m;
    }
  } else {
    if (tr==Factorization::NO_TRANSPOSE) {// min || x || s.t. A x = b
//    if A=[L,0] Q, b=Ax and x=A'x then
//      b=[L,0]Qx=[L,0][v]=Lv and [v]=Qx=[L']s ==> w=0
//                     [w]        [w]    [0 ]
//    so solve Lv=b, compute x=Q'[v]
//                               [0]
      CHECK_SAME(m,b.size());
      CHECK_SAME(n,x.size());
      x=float_zero_;
      x.copyFrom(m,b);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,x.addr(),m,info);
        CHECK_SAME(info,0)
      }
      F77NAME(strsv)('L','N','N',m,QR->addr(),m,x.addr(),1);
      F77NAME(sormlq)('L','T',n,1,m,QR->addr(),m,tau->addr(),
                      x.addr(),n,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormlq)('L','T',n,1,m,QR->addr(),m,tau->addr(),
                      x.addr(),n,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      residual_norm=float_zero_;
      scllen=n;
    } else { // min || b - A' x ||
//    if A=[L,0]Q, r=b-A'x and 0=Ar then
//      [v]=Qr=Qb-[L']x=[y]-[L'x] and 0=[L,0]Qr=Lv ==> v=0, w=z
//      [w]       [0 ]  [z] [ 0 ]
//    so compute [y]=Qb, solve L'x=y
//               [z]
      CHECK_SAME(n,b.size());
      CHECK_SAME(m,x.size());
      Vector<float,float> xtmp(n);
      xtmp.copy(b);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
        CHECK_SAME(info,0)
      }
      F77NAME(sormlq)('L','N',m,1,m,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormlq)('L','N',m,1,m,QR->addr(),m,tau->addr(),
                      xtmp.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      F77NAME(strsv)('L','T','N',m,QR->addr(),m,xtmp.addr(),1);
      x.copyFrom(m,xtmp);
      residual_norm=F77NAME(snrm2)(n-m,xtmp.addr(m),1);
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(slascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(slascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
  return residual_norm;
}

template<> void HouseholderQRFactorization<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &residual_norm,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::solve");
  int m=QR->size(0),n=QR->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  CHECK_SAME(k,residual_norm.size())
  float w=HUGE_VAL;
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(slange)('M',brow,k,B.addr(),B.size(0),&w);
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
  if (m>=n) {
    if (tr==Factorization::NO_TRANSPOSE) { // min || B - A x ||
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<float,float> Xtmp(m,k);
      Xtmp.copy(B);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    Q' B = [ Y \\ Z ]:
      F77NAME(sormqr)('L','T',m,k,n,QR->addr(),m,tau->addr(),
                      Xtmp.addr(),m,&w,lwork,info);
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormqr)('L','T',m,k,n,QR->addr(),m,tau->addr(),
                      Xtmp.addr(),m,work,lwork,info);
      CHECK_SAME(info,0)

//    solve R X = Y:
      F77NAME(strsm)('L','U','N','N',n,k,float_one_,QR->addr(),m,
        Xtmp.addr(),m);
      delete [] work; work=0;
      X.copyFrom('A',n,k,Xtmp);
      for (int j=0;j<k;j++) {
        residual_norm[j]=F77NAME(snrm2)(m-n,Xtmp.addr(n,j),1);
      }
      scllen=n;
    } else { // min || x || s.t. A^t x = B
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))

      X=float_zero_;
      X.copyFrom('A',m,k,B);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,k,X.addr(),m,info);
        CHECK_SAME(info,0)
      }
      F77NAME(strsm)('L','U','T','N',n,k,float_one_,QR->addr(),m,
        X.addr(),m);
      F77NAME(sormqr)('L','N',m,k,n,QR->addr(),m,tau->addr(),
                      X.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormqr)('L','N',m,k,n,QR->addr(),m,tau->addr(),
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
      X=float_zero_;
      X.copyFrom('A',m,k,B);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,k,X.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    solve L X = B:
      F77NAME(strsm)('L','L','N','N',m,k,float_one_,QR->addr(),m,
        X.addr(),n);
//    Q' V
      F77NAME(sormlq)('L','T',n,k,m,QR->addr(),m,tau->addr(),
                      X.addr(),n,&w,lwork,info);
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormlq)('L','T',n,k,m,QR->addr(),m,tau->addr(),
                      X.addr(),n,work,lwork,info);
      CHECK_SAME(info,0)
      delete [] work; work=0;
      residual_norm=float_zero_;
      scllen=n;
    } else { // // min || B - A' X ||
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))

      Matrix<float,float> Xtmp(n,k);
      Xtmp.copy(B);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
        CHECK_SAME(info,0)
      }
      F77NAME(sormlq)('L','N',m,k,m,QR->addr(),m,tau->addr(),
                      Xtmp.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormlq)('L','N',m,k,m,QR->addr(),m,tau->addr(),
                      Xtmp.addr(),m,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;

      F77NAME(strsm)('L','L','T','N',m,k,float_one_,QR->addr(),m,
        Xtmp.addr(),n);
      X.copyFrom('A',m,k,Xtmp);
      residual_norm=float_zero_;
      for (int j=0;j<k;j++) {
        residual_norm[j]=F77NAME(snrm2)(n-m,Xtmp.addr(m,j),1);
      }
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(slascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(slascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
}

template<> void HouseholderQRFactorization<float,float>::improve(
const Vector<float,float> &b,Vector<float,float> &x,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::improve");
  int m=QR->size(0),n=QR->size(1);
  float w=HUGE_VAL;
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(slange)('M',brow,1,b.addr(),b.size(),&w);
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
    F77NAME(slascl)('G',0,0,bnrm,bscl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (iascl!=0) {
    F77NAME(slascl)('G',0,0,ascl,anrm,scllen,1,x.addr(),x.size(),info);
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
      Vector<float,float> r(m);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    r=b-Ax
      F77NAME(sgemv)('N',m,n,float_mone_,A_original->addr(),m,x.addr(),1,
        float_one_,r.addr(),1);
      Vector<float,float> c(n,float_zero_);
//    -dc=A'r
      F77NAME(sgemv)('T',m,n,float_one_,A_original->addr(),m,r.addr(),1,
        float_zero_,c.addr(),1);
//    solve R'(-dv) =-dc:
      F77NAME(strsv)('U','T','N',n,QR->addr(),m,c.addr(),1);
//    solve Rdx=-v:
      F77NAME(strsv)('U','N','N',n,QR->addr(),m,c.addr(),1);
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
      Vector<float,float> r(n);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),m,info);
        CHECK_SAME(info,0)
      }
//    r=b-A'x
      F77NAME(sgemv)('T',m,n,float_mone_,A_original->addr(),m,x.addr(),1,
        float_one_,r.addr(),1);
//    solve R'dv=r:
      Vector<float,float> c(m,float_zero_);
      c.copyFrom(n,r);
      F77NAME(strsv)('U','T','N',n,QR->addr(),m,c.addr(),1);
//    dx=Q[dv\\0]
      F77NAME(sormqr)('L','T',m,1,n,QR->addr(),m,tau->addr(),
                      c.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormqr)('L','T',m,1,n,QR->addr(),m,tau->addr(),
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
      Vector<float,float> v(n,float_zero_);
      v.copyFrom(m,b);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,v.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    r=b-Ax
      F77NAME(sgemv)('N',m,n,float_mone_,A_original->addr(),m,x.addr(),1,
        float_one_,v.addr(),1);
//    solve Ldv=r:
      F77NAME(strsv)('L','N','N',m,QR->addr(),m,v.addr(),1);
//    dx=A'[dv\\dg]:
      F77NAME(sormlq)('L','T',n,1,m,QR->addr(),m,tau->addr(),
                      v.addr(),n,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormlq)('L','T',n,1,m,QR->addr(),m,tau->addr(),
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
      Vector<float,float> r(n);
      r.copy(b);
      if (ibscl!=0) {
        F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,r.addr(),n,info);
        CHECK_SAME(info,0)
      }
//    r=b-A'x
      F77NAME(sgemv)('T',m,n,float_mone_,A_original->addr(),m,x.addr(),1,
        float_one_,r.addr(),1);
      Vector<float,float> c(m,float_zero_);
//    -dc=Ar
      F77NAME(sgemv)('N',m,n,float_one_,A_original->addr(),m,r.addr(),1,
        float_zero_,c.addr(),1);
//    solve L(-dv)=-dc:
      F77NAME(strsv)('L','N','N',m,QR->addr(),m,c.addr(),1);
//    solve L'dx=-dv:
      F77NAME(strsv)('L','T','N',m,QR->addr(),m,c.addr(),1);
      x+=c;
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(slascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(slascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
}

template<> void HouseholderQRFactorization<float,float>::improve(
const Matrix<float,float> &B,Matrix<float,float> &X,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"HouseholderQRFactorization::improve");
  int m=QR->size(0),n=QR->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  float w=HUGE_VAL;
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);
  
  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(slange)('M',brow,k,B.addr(),B.size(0),&w);
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
    F77NAME(slascl)('G',0,0,bnrm,bscl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (iascl!=0) {
    F77NAME(slascl)('G',0,0,ascl,anrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (m>=n) {
    if (tr==Factorization::NO_TRANSPOSE) { // min || B - A x ||
      CHECK_SAME(m,B.size(0))
      CHECK_SAME(n,X.size(0))
      Matrix<float,float> R(m,k);
      R.copy(B);
//    R=B-AX
      F77NAME(sgemm)('N','N',m,k,n,float_mone_,A_original->addr(),m,
        X.addr(),n,float_one_,R.addr(),m);
      Matrix<float,float> C(n,k,float_zero_);
//    -dC=A'R
      F77NAME(sgemm)('T','N',n,k,m,float_one_,A_original->addr(),m,
        R.addr(),m,float_zero_,C.addr(),n);
//    solve R'(-dV)=-dC:
      F77NAME(strsm)('L','U','T','N',n,k,float_one_,QR->addr(),m,
        C.addr(),n);
//    solve RdX=-V:
      F77NAME(strsm)('L','U','N','N',n,k,float_one_,QR->addr(),m,
        C.addr(),n);
      X+=C;
      scllen=n;
    } else { // min || X || s.t. A^t X = B
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<float,float> R(n,k);
      R.copy(B);
//    R=B-A'X
      F77NAME(sgemm)('T','N',n,k,m,float_mone_,A_original->addr(),m,
        X.addr(),m,float_one_,R.addr(),n);
//    solve R'dV=R:
      Matrix<float,float> C(m,k,float_zero_);
      C.copyFrom('A',n,k,R);
      F77NAME(strsm)('L','U','T','N',n,k,float_one_,QR->addr(),m,
        C.addr(),n);
//    dX=Q[dV\\0]
      F77NAME(sormqr)('L','T',m,k,n,QR->addr(),m,tau->addr(),
                      C.addr(),m,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormqr)('L','T',m,k,n,QR->addr(),m,tau->addr(),
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
      Matrix<float,float> V(n,k,float_zero_);
      V.copyFrom('A',m,k,B);
//    R=B-AX
      F77NAME(sgemm)('N','N',m,k,n,float_mone_,A_original->addr(),m,
        X.addr(),n,float_one_,V.addr(),m);
//    solve LdV=R:
      F77NAME(strsm)('L','L','N','N',m,k,float_one_,QR->addr(),m,
        V.addr(),n);
//    dX=Q'[dV\\0]:
      F77NAME(sormlq)('L','T',n,k,m,QR->addr(),m,tau->addr(),
                      V.addr(),n,&w,lwork,info );
      CHECK_SAME(info,0)

      lwork=static_cast<int>(w);
      float *work=OPERATOR_NEW_BRACKET(float,lwork);
      F77NAME(sormlq)('L','T',n,k,m,QR->addr(),m,tau->addr(),
                      V.addr(),n,work,lwork,info );
      CHECK_SAME(info,0)
      delete [] work; work=0;
      X+=V;
      scllen=n;
    } else { // min || B - A' X ||
      CHECK_SAME(n,B.size(0))
      CHECK_SAME(m,X.size(0))
      Matrix<float,float> R(n,k);
      R.copy(B);
//    R=B-A'X
      F77NAME(sgemm)('T','N',n,k,m,float_mone_,A_original->addr(),m,
        X.addr(),m,float_one_,R.addr(),n);
      Matrix<float,float> C(m,k,float_zero_);
//    (-dC)=A R
      F77NAME(sgemm)('N','N',m,k,n,float_one_,A_original->addr(),m,
        R.addr(),n,float_zero_,C.addr(),m);
//    L(-dV)=-dC
      F77NAME(strsm)('L','L','N','N',m,k,float_one_,QR->addr(),m,
        C.addr(),m);
//    L'dX=-dV
      F77NAME(strsm)('L','L','T','N',m,k,float_one_,QR->addr(),m,
        C.addr(),m);
      X+=C;
      scllen=m;
    }
  }
  if (iascl!=0) {
    F77NAME(slascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(slascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
}

template class HouseholderQRFactorization<float,float>;
//template void testHouseholderQRFactorization(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "CompleteOrthogonalDecomposition.C"

template<> CompleteOrthogonalDecomposition<float,float>::
CompleteOrthogonalDecomposition(const Matrix<float,float> &A,
float rcond) : rank(0),URV(0),utau(0),vtau(0),jpvt(0),iascl(0),
ascl(float_one_),anrm(HUGE_VAL),A_original(&A) {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::CompleteOrthogonalDecomposition");
  int m=A.size(0),n=A.size(1);
  URV=OPERATOR_NEW Matrix<float,float>(m,n);
  URV->copy(A);
  float w=HUGE_VAL;
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);
  anrm=F77NAME(slange)('M',m,n,A.addr(),m,&w);
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
    F77NAME(slascl)('G',0,0,anrm,ascl,m,n,URV->addr(),m,info);
    CHECK_SAME(info,0);
  }

  int mn=min(m,n);
  int nb1=F77NAME(ilaenv)(1,"DGEQRF"," ",m,n,-1,1);
  int nb2=F77NAME(ilaenv)(1,"DGERQF"," ",m,n,-1,1);
  int nb=max(nb1,nb2);
  int lwkmin=mn+max(2*mn,n+1);
  int lwork=max(lwkmin,mn+2*n+nb*(n+1));
  float *work=OPERATOR_NEW_BRACKET(float,lwork);

  utau=OPERATOR_NEW Vector<float,float>(mn);
  jpvt=OPERATOR_NEW_BRACKET(int,n);
  F77NAME(sgeqp3)(m,n,URV->addr(),m,jpvt,utau->addr(),work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  float smax=abs((*URV)(0,0));
  if (smax<=float_zero_) return;
  float smin=smax;

  rank=1;
  Vector<float,float> xmin(mn);
  xmin[0]=float_one_;
  Vector<float,float> xmax(mn);
  xmax[0]=float_one_;
  while (rank<mn) {
    int i=rank+1;
    float sminpr=HUGE_VAL,s1=HUGE_VAL,c1=HUGE_VAL;
//                 (job,j,x,sest,w,gamma,sestpr,s,c)
    F77NAME(slaic1)(2,rank,xmin.addr(),smin,URV->addr(0,i-1),
      (*URV)(i-1,i-1),sminpr,s1,c1);
    float smaxpr=HUGE_VAL,s2=HUGE_VAL,c2=HUGE_VAL;
    F77NAME(slaic1)(1,rank,xmax.addr(),smax,URV->addr(0,i-1),
      (*URV)(i-1,i-1),smaxpr,s2,c2);
    if (smaxpr*rcond>sminpr) break;
    for (int i=0;i<rank;i++) {
      xmin[i]=s1*xmin[i];
      xmax[i]=s2*xmax[i];
    }
    xmin[rank-1]=c1;
    xmax[rank-1]=c2;
    smin=sminpr;
    smax=smaxpr;
    rank++;
  }
  if (rank<n) {
    vtau=OPERATOR_NEW Vector<float,float>(rank);
    w=HUGE_VAL;
    int lwork=-1;
    F77NAME(stzrzf)(rank,n,URV->addr(),m,vtau->addr(),&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(stzrzf)(rank,n,URV->addr(),m,vtau->addr(),work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
  }
}

template<> OrthogonalMatrix<float,float>* 
CompleteOrthogonalDecomposition<float,float>::leftOrthogonalPart()
const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::leftOrthogonalPart");
  int m=URV->size(0),n=URV->size(1);
  OrthogonalMatrix<float,float> *Q=
    OPERATOR_NEW OrthogonalMatrix<float,float>(m,m);
  Q->copyFrom('A',m,rank,*URV);

  float w=HUGE_VAL;
  int lwork=-1;
  int info;
  F77NAME(sorgqr)(m,m,rank,Q->addr(),m,utau->addr(),&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(sorgqr)(m,m,rank,Q->addr(),m,utau->addr(),work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work; work=0;
  return Q;
}

template<> float CompleteOrthogonalDecomposition<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve");
  int m=URV->size(0),n=URV->size(1);
  float w=HUGE_VAL;
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(slange)('M',brow,1,b.addr(),b.size(),&w);
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
  float residual_norm=HUGE_VAL;
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
    Vector<float,float> xtmp(m);
    xtmp.copy(b);
    if (ibscl!=0) {
      F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),m,info);
      CHECK_SAME(info,0)
    }
//  U' b = [ y \\ z ]:
    F77NAME(sormqr)('L','T',m,1,mn,URV->addr(),m,utau->addr(),
                    xtmp.addr(),m,&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    float *work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sormqr)('L','T',m,1,mn,URV->addr(),m,utau->addr(),
                    xtmp.addr(),m,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R v = y:
    F77NAME(strsv)('U','N','N',rank,URV->addr(),m,xtmp.addr(),1);
    Vector<float,float> px(n,float_zero_);
    px.copyFrom(rank,xtmp);
    float w=HUGE_VAL;
    lwork=-1;
//  P'x=V[v\\0]:
    F77NAME(sormrz)('L','T',n,1,rank,n-rank,URV->addr(),m,vtau->addr(),
      px.addr(),n,&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sormrz)('L','T',n,1,rank,n-rank,URV->addr(),m,vtau->addr(),
      px.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  permute
    for (int i=0;i<n;i++) x[jpvt[i]-1]=px[i];
    residual_norm=F77NAME(snrm2)(m-rank,xtmp.addr(rank),1);
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
    Vector<float,float> xtmp(n);
//  P'b:
    for (int i=0;i<n;i++) xtmp[i]=b[jpvt[i]-1];
    if (ibscl!=0) {
      F77NAME(slascl)('G',0,0,bnrm,bscl,brow,1,xtmp.addr(),n,info);
      CHECK_SAME(info,0)
    }
    float w=HUGE_VAL;
    lwork=-1;
//  [y\\z]=V'(P'b):
    F77NAME(sormrz)('L','N',n,1,rank,n-rank,URV->addr(),m,vtau->addr(),
      xtmp.addr(),n,&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    float *work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sormrz)('L','N',n,1,rank,n-rank,URV->addr(),m,vtau->addr(),
      xtmp.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R' v = b
    F77NAME(strsv)('U','T','N',rank,URV->addr(),m,xtmp.addr(),1);

    x=float_zero_;
    x.copyFrom(rank,xtmp);
//  U [v\\0]
    w=HUGE_VAL;
    lwork=-1;
    F77NAME(sormqr)('L','N',m,1,mn,URV->addr(),m,utau->addr(),
                    x.addr(),m,&w,lwork,info );
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sormqr)('L','N',m,1,mn,URV->addr(),m,utau->addr(),
                    x.addr(),m,work,lwork,info );
    CHECK_SAME(info,0)
    delete [] work; work=0;
    residual_norm=F77NAME(snrm2)(n-rank,xtmp.addr(rank),1);
  }
  if (iascl!=0) {
    F77NAME(slascl)('G',0,0,anrm,ascl,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(slascl)('G',0,0,bscl,bnrm,scllen,1,x.addr(),x.size(),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
  return residual_norm;
}

template<> void CompleteOrthogonalDecomposition<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,
Vector<float,float> &residual_norm,
Factorization::TRANSPOSE_OPTION tr) const {
//TRACER_CALL(t,"CompleteOrthogonalDecomposition::solve");
  int m=URV->size(0),n=URV->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1))
  CHECK_SAME(k,residual_norm.size())
  float w=HUGE_VAL;
  float smlnum=F77NAME(slamch)('S')/F77NAME(slamch)('P');
  float bignum=float_one_/smlnum;
  F77NAME(slabad)(smlnum,bignum);

  int brow=(tr==Factorization::NO_TRANSPOSE ? m : n);
  float bnrm=F77NAME(slange)('M',brow,k,B.addr(),B.size(0),&w);
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
  residual_norm=HUGE_VAL;
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
    Matrix<float,float> Xtmp(m,k);
    Xtmp.copy(B);
    if (ibscl!=0) {
      F77NAME(slascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),m,info);
      CHECK_SAME(info,0)
    }
//  U' B = [ Y \\ Z ]:
    F77NAME(sormqr)('L','T',m,k,mn,URV->addr(),m,utau->addr(),
                    Xtmp.addr(),m,&w,lwork,info);
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    float *work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sormqr)('L','T',m,k,mn,URV->addr(),m,utau->addr(),
                    Xtmp.addr(),m,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R V = Y:
    F77NAME(strsv)('U','N','N',rank,URV->addr(),m,Xtmp.addr(),1);
    Matrix<float,float> PX(n,k,float_zero_);
    PX.copyFrom('A',rank,k,Xtmp);
    float w=HUGE_VAL;
    lwork=-1;
//  P'X=V[V\\0]:
    F77NAME(sormrz)('L','T',n,k,rank,n-rank,URV->addr(),m,vtau->addr(),
      PX.addr(),n,&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sormrz)('L','T',n,k,rank,n-rank,URV->addr(),m,vtau->addr(),
      PX.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  permute
    for (int j=0;j<k;j++) {
      for (int i=0;i<n;i++) X(jpvt[i]-1,j)=PX(i,j);
      residual_norm[j]=F77NAME(snrm2)(m-rank,Xtmp.addr(rank,j),1);
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
    Matrix<float,float> Xtmp(n,k);
//  P'B:
    for (int j=0;j<k;j++) {
      for (int i=0;i<n;i++) Xtmp(i,j)=B(jpvt[i]-1,j);
    }
    if (ibscl!=0) {
      F77NAME(slascl)('G',0,0,bnrm,bscl,brow,k,Xtmp.addr(),n,info);
      CHECK_SAME(info,0)
    }
    float w=HUGE_VAL;
    lwork=-1;
//  [Y\\Z]=V'(P'B):
    F77NAME(sormrz)('L','N',n,k,rank,n-rank,URV->addr(),m,vtau->addr(),
      Xtmp.addr(),n,&w,lwork,info);
    CHECK_SAME(info,0)
    lwork=static_cast<int>(w);
    float *work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sormrz)('L','N',n,k,rank,n-rank,URV->addr(),m,vtau->addr(),
      Xtmp.addr(),n,work,lwork,info);
    CHECK_SAME(info,0)
    delete [] work; work=0;
//  solve R' V = B
    F77NAME(strsv)('U','T','N',rank,URV->addr(),m,Xtmp.addr(),1);

    X=float_zero_;
    X.copyFrom('A',rank,k,Xtmp);
//  U [V\\0]
    w=HUGE_VAL;
    lwork=-1;
    F77NAME(sormqr)('L','N',m,k,mn,URV->addr(),m,utau->addr(),
                    X.addr(),m,&w,lwork,info );
    CHECK_SAME(info,0)

    lwork=static_cast<int>(w);
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sormqr)('L','N',m,k,mn,URV->addr(),m,utau->addr(),
                    X.addr(),m,work,lwork,info );
    CHECK_SAME(info,0)
    delete [] work; work=0;
    for (int j=0;j<k;j++) {
      residual_norm[j]=F77NAME(snrm2)(n-rank,Xtmp.addr(rank,j),1);
    }
  }
  if (iascl!=0) {
    F77NAME(slascl)('G',0,0,anrm,ascl,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
  }
  if (ibscl!=0) {
    F77NAME(slascl)('G',0,0,bscl,bnrm,scllen,k,X.addr(),X.size(0),info);
    CHECK_SAME(info,0)
    residual_norm*=bnrm/bscl;
  }
}

template class CompleteOrthogonalDecomposition<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "GramSchmidtQRFactorization.C"
template<> float GramSchmidtQRFactorization<float,float>::tol=
  float_one_/sqrt(2.);

template<> void 
GramSchmidtQRFactorization<float,float>::reorthogonalize(int j) {
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,R->size(0))
  CHECK_BOUNDS(j,0,n)
  for (int k=0;k<j;k++) {
    float dot=F77NAME(sdot)(m,Q->addr(0,k),1,Q->addr(0,j),1)/(*R)(k,k);
    (*R)(k,j)+=dot;
    F77NAME(saxpy)(m,-dot,Q->addr(0,k),1,Q->addr(0,j),1);
  }
  (*R)(j,j)=F77NAME(sdot)(m,Q->addr(0,j),1,Q->addr(0,j),1);
}

template<> 
GramSchmidtQRFactorization<float,float>::GramSchmidtQRFactorization(
const Matrix<float,float> &A) : Q(0),R(0),A_original(&A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::GramSchmidtQRFactorization");
  int m=A.size(0),n=A.size(1);
  assert(m>=n);
  Q=OPERATOR_NEW OrthogonalMatrix<float,float>(m,n);
  R=OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(n,n);
  Q->copy(A);

  float *norm=OPERATOR_NEW_BRACKET(float,n);
  for (int j=0;j<n;j++) {
    norm[j]=F77NAME(sdot)(m,Q->addr(0,j),1,Q->addr(0,j),1);
  }
  float tol2=tol*tol;
  for (int k=0;k<n;k++) {
    (*R)(k,k)=F77NAME(sdot)(m,Q->addr(0,k),1,Q->addr(0,k),1);
    if ((*R)(k,k)<tol2*norm[k]) reorthogonalize(k);
    if ((*R)(k,k)>float_zero_ && k<n-1) {
      F77NAME(sgemv)('T',m,n-k-1,float_one_/(*R)(k,k),Q->addr(0,k+1),m,
        Q->addr(0,k),1,float_zero_,R->addr(k,k+1),n);
      F77NAME(sger)(m,n-k-1,float_mone_,Q->addr(0,k),1,R->addr(k,k+1),n,
        Q->addr(0,k+1),m);
    }
  }
  delete [] norm;
}

template<> void
GramSchmidtQRFactorization<float,float>::solveOverdetermined(
const Vector<float,float> &b,Vector<float,float> &x,
Vector<float,float> &residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,b.size());
  CHECK_SAME(n,x.size());
  CHECK_SAME(m,residual.size());
  residual.copy(b);
  for (int j=0;j<n;j++) {
    float dot=-F77NAME(sdot)(m,Q->addr(0,j),1,residual.addr(),1)
             /(*R)(j,j);
    x[j]=-dot;
    F77NAME(saxpy)(m,dot,Q->addr(0,j),1,residual.addr(),1);
  }
  F77NAME(strsv)('U','N','U',n,R->addr(),n,x.addr(),1);
}

template<> void
GramSchmidtQRFactorization<float,float>::solveOverdetermined(
const Matrix<float,float> &B,Matrix<float,float> &X,
Matrix<float,float> &Residual) const { 
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,X.size(0));
  CHECK_SAME(k,X.size(1));
  CHECK_SAME(m,Residual.size(0));
  CHECK_SAME(k,Residual.size(1));
  Residual.copy(B);
  for (int j=0;j<n;j++) {
    F77NAME(sgemv)('T',m,k,float_one_/(*R)(j,j),Residual.addr(),m,
      Q->addr(0,j),1,float_zero_,X.addr(j,0),n);
    F77NAME(sger)(m,k,float_mone_,Q->addr(0,j),1,X.addr(j,0),n,
      Residual.addr(),m);
  }
  F77NAME(strsm)('L','U','N','U',n,k,float_one_,R->addr(),n,X.addr(),n);
}

template<> void
GramSchmidtQRFactorization<float,float>::solveUnderdetermined(
const Vector<float,float> &b,Vector<float,float> &x) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveUnderdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,b.size());
  CHECK_SAME(m,x.size());

  Vector<float,float> *v=OPERATOR_NEW Vector<float,float>(n);
  v->copy(b);

  F77NAME(strsv)('U','T','U',n,R->addr(),n,v->addr(),1);
  x=float_zero_;
  for (int j=0;j<n;j++) {
    float omega=((*v)[j]-F77NAME(sdot)(m,Q->addr(0,j),1,x.addr(),1))
               /(*R)(j,j);
    F77NAME(saxpy)(m,omega,Q->addr(0,j),1,x.addr(),1);
  }
  delete v; v=0;
}

template<>
void GramSchmidtQRFactorization<float,float>::solveUnderdetermined(
const Matrix<float,float> &B,Matrix<float,float> &X) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveUnderdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(n,B.size(0));
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(k,X.size(1));

  Matrix<float,float> *V=OPERATOR_NEW Matrix<float,float>(n,k);
  V->copy(B);
  F77NAME(strsm)('L','U','T','U',n,k,float_one_,R->addr(),n,
    V->addr(),n);
  X=float_zero_;
  float *work=OPERATOR_NEW_BRACKET(float,k);
  for (int j=0;j<n;j++) {
    F77NAME(sgemv)('T',m,k,float_one_,X.addr(),m,Q->addr(0,j),1,
      float_zero_,work,1);
    for (int kk=0;kk<k;kk++) work[kk]=((*V)(j,kk)-work[kk])/(*R)(j,j);
    F77NAME(sger)(m,k,float_one_,Q->addr(0,j),1,work,1,X.addr(),m);
  }
  delete [] work; work=0;
  delete V; V=0;
}

template<> void
GramSchmidtQRFactorization<float,float>::improveOverdetermined(
const Vector<float,float> &b,Vector<float,float> &x,
Vector<float,float> &residual) const { 
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveOverdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,b.size());
  CHECK_SAME(n,x.size());
  CHECK_SAME(m,residual.size());
  Vector<float,float> *delta_b=OPERATOR_NEW Vector<float,float>(m);
  Vector<float,float> *delta_z=OPERATOR_NEW Vector<float,float>(n,0.);
  Vector<float,float> *delta_y=OPERATOR_NEW Vector<float,float>(n,0.);
  F77NAME(sgemv)('N',m,n,float_mone_,A_original->addr(),m,x.addr(),1,
    float_zero_,delta_b->addr(),1);
  F77NAME(sgemv)('T',m,n,float_mone_,A_original->addr(),m,
    residual.addr(),1,float_zero_,delta_z->addr(),1);
  for (int i=0;i<m;i++) (*delta_b)[i]+=b[i]-residual[i];
  F77NAME(strsv)('U','T','U',n,R->addr(),n,delta_z->addr(),1);
  for (int j=0;j<n;j++) {
    (*delta_y)[j]=(F77NAME(sdot)(m,Q->addr(0,j),1,delta_b->addr(),1)
                 -(*delta_z)[j])/(*R)(j,j);
    F77NAME(saxpy)(m,-(*delta_y)[j],Q->addr(0,j),1,delta_b->addr(),1);
  }
  F77NAME(strsv)('U','N','U',n,R->addr(),n,delta_y->addr(),1);
  x+=(*delta_y);
  residual+=(*delta_b);
  delete delta_y; delta_y=0;
  delete delta_z; delta_z=0;
  delete delta_b; delta_b=0;
}

template<> void
GramSchmidtQRFactorization<float,float>::improveOverdetermined(
const Matrix<float,float> &B,Matrix<float,float> &X,
Matrix<float,float> &Residual) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::solveOverdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,X.size(0));
  CHECK_SAME(k,X.size(1));
  CHECK_SAME(m,Residual.size(0));
  CHECK_SAME(k,Residual.size(1));
  Matrix<float,float> *delta_B=OPERATOR_NEW Matrix<float,float>(m,k);
  Matrix<float,float> *delta_Z=
    OPERATOR_NEW Matrix<float,float>(n,k,0.);
  Matrix<float,float> *delta_Y=
    OPERATOR_NEW Matrix<float,float>(n,k,0.);
  F77NAME(sgemm)('N','N',m,k,n,float_mone_,A_original->addr(),m,
    X.addr(),n,float_zero_,delta_B->addr(),m);
  F77NAME(sgemm)('T','N',n,k,m,float_mone_,A_original->addr(),m,
    Residual.addr(),m,float_zero_,delta_Z->addr(),n);
  for (int j=0;j<k;j++) {
    for (int i=0;i<m;i++) (*delta_B)(i,j)+=B(i,j)-Residual(i,j);
  }
  F77NAME(strsm)('L','U','T','U',n,k,float_one_,R->addr(),n,
    delta_Z->addr(),n);
  float *work=OPERATOR_NEW_BRACKET(float,k);
  for (int j=0;j<n;j++) {
    F77NAME(sgemv)('T',m,k,float_one_,delta_B->addr(),m,Q->addr(0,j),1,
      float_zero_,work,1);
    for (int kk=0;kk<k;kk++) {
      (*delta_Y)(j,kk)=(work[kk]-(*delta_Z)(j,kk))/(*R)(j,j);
    }
    F77NAME(sger)(m,k,float_mone_,Q->addr(0,j),1,delta_Y->addr(j,0),n,
      delta_B->addr(),m);
  }
  delete [] work; work=0;
  F77NAME(strsm)('L','U','N','U',n,k,float_one_,R->addr(),n,
    delta_Y->addr(),n);
  X+=(*delta_Y);
  Residual+=(*delta_B);
  delete delta_Y; delta_Y=0;
  delete delta_Z; delta_Z=0;
  delete delta_B; delta_B=0;
}

template<>
void GramSchmidtQRFactorization<float,float>::improveUnderdetermined(
const Vector<float,float> &b,Vector<float,float> &x) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveUnderdetermined");
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,b.size());
  CHECK_SAME(m,x.size());

  Vector<float,float> *y=OPERATOR_NEW Vector<float,float>(n);
  Vector<float,float> *delta_x=OPERATOR_NEW Vector<float,float>(m);
  Vector<float,float> *delta_by=OPERATOR_NEW Vector<float,float>(n);
  y->copy(b);
  delta_x->copy(x);
  delta_by->copy(b);
  F77NAME(strsv)('U','T','U',n,R->addr(),n,y->addr(),1);
  for (int j=0;j<n;j++) (*y)[j]/=(*R)(j,j);
  F77NAME(sgemv)('N',m,n,float_one_,Q->addr(),m,y->addr(),1,
    float_mone_,delta_x->addr(),1);
  F77NAME(sgemv)('T',m,n,float_mone_,A_original->addr(),m,x.addr(),1,
    float_one_,delta_by->addr(),1);
  F77NAME(strsv)('U','T','U',n,R->addr(),n,delta_by->addr(),1);
  for (int j=0;j<n;j++) {
    float omega=(F77NAME(sdot)(m,Q->addr(0,j),1,delta_x->addr(),1)
      -(*delta_by)[j])/(*R)(j,j);
    F77NAME(saxpy)(m,-omega,Q->addr(0,j),1,delta_x->addr(),1);
  }
  x+=(*delta_x);
  delete delta_by; delta_by=0;
  delete delta_x; delta_x=0;
  delete y; y=0;
}

template<>
void GramSchmidtQRFactorization<float,float>::improveUnderdetermined(
const Matrix<float,float> &B,Matrix<float,float> &X) const {
//TRACER_CALL(t,"GramSchmidtQRFactorization::improveUnderdetermined");
  int m=Q->size(0),n=Q->size(1),k=B.size(1);
  CHECK_SAME(n,B.size(0));
  CHECK_SAME(m,X.size(0));
  CHECK_SAME(k,X.size(1));

  Matrix<float,float> *Y=OPERATOR_NEW Matrix<float,float>(n,k);
  Matrix<float,float> *delta_X=OPERATOR_NEW Matrix<float,float>(m,k);
  Matrix<float,float> *delta_BY=OPERATOR_NEW Matrix<float,float>(n,k);
  Y->copy(B);
  delta_X->copy(X);
  delta_BY->copy(B);

  F77NAME(strsm)('L','U','T','U',n,k,float_one_,R->addr(),n,Y->addr(),n);
  for (int j=0;j<n;j++) {
    F77NAME(sscal)(k,float_one_/(*R)(j,j),Y->addr(j,0),n);
  }
  F77NAME(sgemm)('N','N',m,k,n,float_one_,Q->addr(),m,Y->addr(),n,
    float_mone_,delta_X->addr(),m);
  F77NAME(sgemm)('T','N',n,k,m,float_mone_,A_original->addr(),m,
    X.addr(),m,float_one_,delta_BY->addr(),n);
  F77NAME(strsm)('L','U','T','U',n,k,float_one_,R->addr(),n,
    delta_BY->addr(),n);
  float *work=OPERATOR_NEW_BRACKET(float,k);
  for (int j=0;j<n;j++) {
    F77NAME(sgemv)('T',m,k,float_one_,delta_X->addr(),m,Q->addr(0,j),1,
      float_zero_,work,1); // w^T = delta_X^T q_j
    for (int kk=0;kk<k;kk++) {
      work[kk]=(work[kk]-(*delta_BY)(j,kk))/(*R)(j,j);
    }
    F77NAME(sger)(m,k,float_mone_,Q->addr(0,j),1,work,1,
      delta_X->addr(),m);
  }
  X+=(*delta_X);
  delete [] work; work=0;
  delete delta_BY; delta_BY=0;
  delete delta_X; delta_X=0;
  delete Y; Y=0;
}

template<> void GramSchmidtQRFactorization<float,float>::exchangeColumn(
int j,int jAin,int jAout,const Matrix<float,float> &A) {
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
  Vector<float,float> dif(m);
  for (int i=0;i<m;i++) dif[i]=A(i,jAin)-A(i,jAout);
#ifdef DEBUG
//cout << "\tdif = " << endl;
//dif.printOn(cout);
#endif
  Vector<float,float> QTdif(n,float_zero_);
  Q->gemv(float_one_,dif,float_zero_,QTdif,'T');
#ifdef DEBUG
//cout << "\tQTdif = " << endl;
//QTdif.printOn(cout);
#endif
  Vector<float,float> subdiag(n-1,float_zero_);
  for (int k=n-2;k>j;k--) {
    float c,s;
    F77NAME(srotg)(QTdif[k],QTdif[k+1],c,s);
#ifdef DEBUG
//  cout << "\tk,c,s = " << k << " " << c << " " << s << endl;
#endif
    int ncols=n-k-1;
    F77NAME(srot)(1,R->addr(k,k),n,subdiag.addr(k),1,c,s);
    F77NAME(srot)(n-k-1,R->addr(k,k+1),n,R->addr(k+1,k+1),n,c,s);
    F77NAME(srot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
#ifdef DEBUG
//  printOn(cout);
#endif
  }
  F77NAME(saxpy)(j+1,float_one_,QTdif.addr(),1,R->addr(0,j),1);
  if (j+1<n) {
    subdiag[j]=QTdif[j+1];
    for (int k=j;k<n-1;k++) {
      float c,s;
      F77NAME(srotg)((*R)(k,k),subdiag[k],c,s);
#ifdef DEBUG
//    cout << "\tk,c,s = " << k << " " << c << " " << s << endl;
#endif
      F77NAME(srot)(n-k-1,R->addr(k,k+1),n,R->addr(k+1,k+1),n,c,s);
      F77NAME(srot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
#ifdef DEBUG
//    printOn(cout);
#endif
    }
  }
#ifdef DEBUG
//printOn(cout);
#endif
}

template<> void GramSchmidtQRFactorization<float,float>::exchangeRow(
int i,int iAin,int iAout,const Matrix<float,float> &A) {
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
  Vector<float,float> dif(n);
  for (int j=0;j<n;j++) dif[j]=A(iAin,j)-A(iAout,j);
#ifdef DEBUG
//cout << "\tdif = " << endl;
//dif.printOn(cout);
#endif
  Vector<float,float> axis(m,float_zero_);
  axis[i]=float_one_;
  Vector<float,float> QTaxis(n,float_zero_);
  Q->gemv(float_one_,axis,float_zero_,QTaxis,'T');
#ifdef DEBUG
//cout << "\tQTaxis = " << endl;
//QTaxis.printOn(cout);
#endif
  Vector<float,float> subdiag(n-1,float_zero_);
  for (int k=n-2;k>=0;k--) {
    float c,s;
    F77NAME(srotg)(QTaxis[k],QTaxis[k+1],c,s);
#ifdef DEBUG
//  cout << "\tk,c,s = " << k << " " << c << " " << s << endl;
#endif
    int ncols=n-k-1;
    F77NAME(srot)(1,R->addr(k,k),n,subdiag.addr(k),1,c,s);
    F77NAME(srot)(n-k-1,R->addr(k,k+1),n,R->addr(k+1,k+1),n,c,s);
    F77NAME(srot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
#ifdef DEBUG
//  printOn(cout);
#endif
  }
  F77NAME(saxpy)(n,QTaxis[0],dif.addr(),1,R->addr(),n);
  for (int k=0;k<n-1;k++) {
    float c,s;
    F77NAME(srotg)((*R)(k,k),subdiag[k],c,s);
#ifdef DEBUG
//    cout << "\tk,c,s = " << k << " " << c << " " << s << endl;
#endif
      F77NAME(srot)(n-k-1,R->addr(k,k+1),n,R->addr(k+1,k+1),n,c,s);
      F77NAME(srot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
#ifdef DEBUG
//    printOn(cout);
#endif
    }
#ifdef DEBUG
//printOn(cout);
#endif
}

template<> void GramSchmidtQRFactorization<float,float>::dropColumn(
int j) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropColumn");
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(j,0,n)
  for (int k=j+1;k<n;k++) {
    float c,s;
    F77NAME(srotg)(R->operator()(k,k),R->operator()(j,k),c,s);
    int ncols=n-k-1;
    if (ncols>0) {
      F77NAME(srot)(ncols,R->addr(k,k+1),n,R->addr(j,k+1),n,c,s);
    }
    F77NAME(srot)(m,Q->addr(0,k),1,Q->addr(0,j),1,c,s);
  }

  for (int k=j+1;k<n;k++) {
    for (int i=0;i<j;i++) (*R)(i,k-1)=(*R)(i,k);
    for (int i=j;i<k;i++) (*R)(i,k-1)=(*R)(i+1,k);
    F77NAME(scopy)(m,Q->addr(0,k),1,Q->addr(0,k-1),1);
  }

  OrthogonalMatrix<float,float> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,float>(m,n-1);
  new_Q->copyFrom('A',m,n-1,*Q);
  delete Q;
  Q=new_Q;
  UpperTrapezoidalMatrix<float,float> *new_R
    =OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(n-1,n-1);
  new_R->copyFrom(n-1,n-1,*R);
  delete R;
  R=new_R;
}

template<> void GramSchmidtQRFactorization<float,float>::addColumn(
int j,const Matrix<float,float> &A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::addColumn");
#ifdef DEBUG
//printOn(cout);
//cout << "\tj = " << j << endl;
//cout << "\tA = " << endl;
//A.printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(m,A.size(0))

  OrthogonalMatrix<float,float> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,float>(m,n+1);
  new_Q->copyFrom('A',m,n,*Q);
  delete Q;
  Q=new_Q;
  F77NAME(scopy)(m,A.addr(0,j),1,Q->addr(0,n),1);

  UpperTrapezoidalMatrix<float,float> *new_R
    =OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(n+1,n+1);
  new_R->copyFrom(n,n,*R);
  delete R;
  R=new_R;

  float norm=F77NAME(snrm2)(m,Q->addr(0,n),1);
  for (j=0;j<n;j++) {
    float dot=-F77NAME(sdot)(m,Q->addr(0,j),1,Q->addr(0,n),1);
    (*R)(j,n)=-dot;
    F77NAME(saxpy)(m,dot,Q->addr(0,j),1,Q->addr(0,n),1);
  }
  (*R)(n,n)=F77NAME(snrm2)(m,Q->addr(0,n),1);
  if ((*R)(n,n)<tol*norm) reorthogonalize(n);
  float a=float_one_/(*R)(n,n);
  F77NAME(sscal)(m,a,Q->addr(0,n),1);
#ifdef DEBUG
//printOn(cout);
#endif
}

template<> void GramSchmidtQRFactorization<float,float>::dropRow(int i)
{
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropRow");
#ifdef DEBUG
//printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(i,0,m)
//permute row i to last
  for (int k=i+1;k<m;k++) {
    F77NAME(sswap)(n,Q->addr(k-1,0),m,Q->addr(k,0),m);
  }
#ifdef DEBUG
//Q->printOn(cout);
#endif
//the following from Daniel, Gragg, Kaufman and Stewart fails if
//  the last axis vector was already a column of Q
  Vector<float,float> work_column(m);
  int Mmin1=m-1;
  F77NAME(sgemv)('N',Mmin1,n,-float_one_,Q->addr(),m,Q->addr(m-1,0),m,
		 float_zero_,work_column.addr(),1);
  work_column[m-1]=sqrt(float_one_
	  -F77NAME(sdot)(n,Q->addr(m-1,0),n,Q->addr(m-1,0),m));
  float alpha=float_one_/work_column[m-1];
  F77NAME(sscal)(Mmin1,alpha,work_column.addr(),1);

  Vector<float,float> work_row(n,float_zero_);
  for (int j=n-1;j>=0;j--) {
    float c,s;
    F77NAME(srotg)(work_column[Mmin1],Q->operator()(Mmin1,j),c,s);
    F77NAME(srot)(Mmin1,work_column.addr(),1,Q->addr(0,j),1,c,s);
    int ncols=n-j;
    F77NAME(srot)(ncols,work_row.addr(j),1,R->addr(j,j),n,c,s);
  }

  OrthogonalMatrix<float,float> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,float>(Mmin1,n);
  new_Q->copyFrom('A',Mmin1,n,*Q);
  delete Q;
  Q=new_Q;
#ifdef DEBUG
//printOn(cout);
#endif
}

template<> void
GramSchmidtQRFactorization<float,float>::dropRowAndColumn(int i,int j)
{
//TRACER_CALL(t,"GramSchmidtQRFactorization::dropRowAndColumn");
#ifdef DEBUG
//cout << "\ti,j = " << i << " " << j << endl;
//printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_BOUNDS(i,0,m)
  CHECK_BOUNDS(j,0,n)
  Vector<float,float> rwork(n-1);
  UpperTrapezoidalMatrix<float,float> *new_R
    =OPERATOR_NEW UpperTrapezoidalMatrix<float,float>(n-1,n-1);
  for (int k=0;k<j;k++) {
    memcpy(new_R->addr(0,k),R->addr(0,k),(k+1)*sizeof(float));
  }
  for (int k=j+1;k<n;k++) {
    memcpy(new_R->addr(0,k-1),R->addr(0,k),k*sizeof(float));
    rwork[k-1]=(*R)(k,k);
  }
  for (int k=j;k<n-1;k++) {
    float c,s;
    F77NAME(srotg)((*new_R)(k,k),rwork[k],c,s);
    int ncols=n-k-2;
    if (ncols>0) {
      F77NAME(srot)(ncols,new_R->addr(k,k+1),n-1,
        new_R->addr(k+1,k+1),n-1,c,s);
    }
    F77NAME(srot)(m,Q->addr(0,k),1,Q->addr(0,k+1),1,c,s);
  }
  OrthogonalMatrix<float,float> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,float>(m-1,n-1);
  for (int k=0;k<i;k++) {
    F77NAME(scopy)(n-1,Q->addr(k,0),m,new_Q->addr(k,0),m-1);
  }
  for (int k=i+1;k<m;k++) {
    F77NAME(scopy)(n-1,Q->addr(k-1,0),m,new_Q->addr(k-1,0),m-1);
  }
  Vector<float,float> qwork(m-1); // last column of Q
  memcpy(qwork.addr(),Q->addr(0,n-1),(m-1)*sizeof(float));
  rwork=float_zero_;
  for (int k=n-2;k>=0;k--) {
    float c,s;
    F77NAME(srotg)((*Q)(n-1,n-1),(*Q)(n-1,k),c,s);
    F77NAME(srot)(m-1,qwork.addr(),1,new_Q->addr(0,k),1,c,s);
    F77NAME(srot)(n-1-k,rwork.addr(k),1,new_R->addr(k,k),n-1,c,s);
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

template<> void GramSchmidtQRFactorization<float,float>::addRow(int i,
const Matrix<float,float> &A) {
//TRACER_CALL(t,"GramSchmidtQRFactorization::addRow");
#ifdef DEBUG
//printOn(cout);
//cout << "\ti = " << i << endl;
//cout << "\tA = " << endl;
//A.printOn(cout);
#endif
  int m=Q->size(0),n=Q->size(1);
  CHECK_SAME(n,A.size(1))

  OrthogonalMatrix<float,float> *new_Q=
    OPERATOR_NEW OrthogonalMatrix<float,float>(m+1,n);
  new_Q->copyFrom('A',m,n,*Q);
  delete Q;
  Q=new_Q;
  for (int j=0;j<n;j++) (*Q)(m,j)=float_zero_;

  Vector<float,float> work_row(n);
  F77NAME(scopy)(n,A.addr(i,0),A.size(0),work_row.addr(),1);

  Vector<float,float> work_column(m+1,float_zero_);
  work_column[m]=float_one_;

  int Mp1=m+1;
  for (int j=0;j<n;j++) {
    float c,s;
    F77NAME(srotg)((*R)(j,j),work_row[j],c,s);
    int ncols=n-j-1;
    if (ncols>0) {
      F77NAME(srot)(ncols,R->addr(j,j+1),n,work_row.addr(j+1),1,c,s);
    }
    F77NAME(srot)(Mp1,Q->addr(0,j),1,work_column.addr(),1,c,s);
  }
#ifdef DEBUG
//printOn(cout);
#endif
}

template class GramSchmidtQRFactorization<float,float>;
template void testGramSchmidtQRFactorization(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "SingularValueDecomposition.C"

template<> SingularValueDecomposition<float,float>::
SingularValueDecomposition(const Matrix<float,float> &A) : U(0),
Vtranspose(0),s(0),A_original(&A) {
//TRACER_CALL(t,"SingularValueDecomposition::SingularValueDecomposition");
  int m=A.size(0),n=A.size(1);
  int mn=min(m,n);
  U=OPERATOR_NEW OrthogonalMatrix<float,float>(m,mn);
  Vtranspose=OPERATOR_NEW OrthogonalMatrix<float,float>(mn,n);
  s=OPERATOR_NEW Vector<float,float>(mn);
  Matrix<float,float> Acopy(m,n);
  Acopy.copy(A);

  int info=0;
  float w=HUGE_VAL;
  int lwork=-1;
  F77NAME(sgesvd)('S','S',m,n,Acopy.addr(),m,s->addr(),U->addr(),m,
    Vtranspose->addr(),mn,&w,lwork,info);
  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(sgesvd)('S','S',m,n,Acopy.addr(),m,s->addr(),U->addr(),m,
    Vtranspose->addr(),mn,work,lwork,info);
  CHECK_SAME(info,0);
  delete [] work; work=0;
}

template<> void SingularValueDecomposition<float,float>::solve(
const Vector<float,float> &b,Vector<float,float> &x,
float rcond,Factorization::TRANSPOSE_OPTION to) const {
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
  float *y=OPERATOR_NEW_BRACKET(float,rank);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(sgemv)('C',m,rank,float_one_,U->addr(),m,b.addr(),1,
      float_zero_,y,1);
    for (int i=0;i<rank;i++) y[i]/=(*s)[i];
    F77NAME(sgemv)('C',rank,n,float_one_,Vtranspose->addr(),mn,y,1,
      float_zero_,x.addr(),1);
  } else {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(sgemv)('N',rank,n,float_one_,Vtranspose->addr(),mn,
      b.addr(),1,float_zero_,y,1);
    for (int i=0;i<rank;i++) y[i]/=(*s)[i];
    F77NAME(sgemv)('N',m,rank,float_one_,U->addr(),m,y,1,float_zero_,
      x.addr(),1);
  }
  delete [] y; y=0;
}

template<> void SingularValueDecomposition<float,float>::regularize(
const Vector<float,float> &b,Vector<float,float> &x,
float ridge,Factorization::TRANSPOSE_OPTION to) const {
  CHECK_TEST(ridge>=float_zero_);
  int m=A_original->size(0),n=A_original->size(1);
  int mn=min(m,n);
  float sr=sqrt(ridge);
  float *y=OPERATOR_NEW_BRACKET(float,mn);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,b.size());
    CHECK_SAME(n,x.size());
    F77NAME(sgemv)('C',m,mn,float_one_,U->addr(),m,b.addr(),1,
      float_zero_,y,1);
    for (int i=0;i<mn;i++) {
      float si=(*s)[i];
      float scale=max(si,sr);
      si/=scale;
      float ri=sr/scale;
      y[i]*=si/(scale*(si*si+ri*ri));
    }
    F77NAME(sgemv)('C',mn,n,float_one_,Vtranspose->addr(),mn,y,1,
      float_zero_,x.addr(),1);
  } else {
    CHECK_SAME(n,b.size());
    CHECK_SAME(m,x.size());
    F77NAME(sgemv)('N',mn,n,float_one_,Vtranspose->addr(),mn,b.addr(),1,
      float_zero_,y,1);
    for (int i=0;i<mn;i++) {
      float si=(*s)[i];
      float scale=max(si,sr);
      si/=scale;
      float ri=sr/scale;
      y[i]*=si/(scale*(si*si+ri*ri));
    }
    F77NAME(sgemv)('N',m,mn,float_one_,U->addr(),m,y,1,float_zero_,
      x.addr(),1);
  }
  delete [] y; y=0;
}

template<> void SingularValueDecomposition<float,float>::solve(
const Matrix<float,float> &B,Matrix<float,float> &X,
float rcond,Factorization::TRANSPOSE_OPTION to) const {
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
  float *Y=OPERATOR_NEW_BRACKET(float,rank*k);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,X.size(0));
    F77NAME(sgemm)('C','N',rank,k,m,float_one_,U->addr(),m,B.addr(),m,
      float_zero_,Y,rank);
    for (int i=0;i<rank;i++) F77NAME(srscl)(k,(*s)[i],Y+i,rank);
    F77NAME(sgemm)('C','N',n,k,rank,float_one_,Vtranspose->addr(),mn,
      Y,rank,float_zero_,X.addr(),n);
  } else {
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(m,X.size(0));
    F77NAME(sgemm)('N','N',rank,k,n,float_one_,Vtranspose->addr(),mn,
      B.addr(),n,float_zero_,Y,rank);
    for (int i=0;i<rank;i++) F77NAME(srscl)(k,(*s)[i],Y+i,rank);
    F77NAME(sgemm)('N','N',m,k,rank,float_one_,U->addr(),m,Y,rank,
      float_zero_,X.addr(),m);
  }
  delete [] Y; Y=0;
}

template<> void SingularValueDecomposition<float,float>::regularize(
const Matrix<float,float> &B,Matrix<float,float> &X,
float ridge,Factorization::TRANSPOSE_OPTION to) const {
  CHECK_TEST(ridge>=float_zero_);
  int m=A_original->size(0),n=A_original->size(1),k=B.size(1);
  CHECK_SAME(k,X.size(1));
  int mn=min(m,n);
  float sr=sqrt(ridge);
  float *Y=OPERATOR_NEW_BRACKET(float,mn*k);
  if (to==Factorization::NO_TRANSPOSE) {
    CHECK_SAME(m,B.size(0));
    CHECK_SAME(n,X.size(0));
    F77NAME(sgemm)('C','N',mn,k,m,float_one_,U->addr(),m,B.addr(),m,
      float_zero_,Y,mn);
    for (int i=0;i<mn;i++) {
      float si=(*s)[i];
      float scale=max(abs(si),sr);
      si/=scale;
      float ri=sr/scale;
      F77NAME(sscal)(k,si/(scale*(si*si+ri*ri)),Y+i,mn);
    }
    F77NAME(sgemm)('C','N',n,k,mn,float_one_,Vtranspose->addr(),mn,Y,mn,
      float_zero_,X.addr(),n);
  } else {
    CHECK_SAME(n,B.size(0));
    CHECK_SAME(m,X.size(0));
    F77NAME(sgemm)('N','N',mn,k,n,float_one_,Vtranspose->addr(),mn,
      B.addr(),n,float_zero_,Y,mn);
    for (int i=0;i<mn;i++) {
      float si=(*s)[i];
      float scale=max(abs(si),sr);
      si/=scale;
      float ri=sr/scale;
      F77NAME(sscal)(k,si/(scale*(si*si+ri*ri)),Y+i,mn);
    }
    F77NAME(sgemm)('N','N',m,k,mn,float_one_,U->addr(),m,Y,mn,
      float_zero_,X.addr(),m);
  }
  delete [] Y; Y=0;
}

template class SingularValueDecomposition<float,float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> CholeskyFactorization<float,float>::CholeskyFactorization(
const SymmetricPositiveMatrix<float,float> &A) : anorm(float_zero_) {
//TRACER_CALL(t,"CholeskyFactorization::CF");
  int n=A.size(0);
  int inc=1;
  for (int j=0;j<n;j++) {
    int k=n-j;
    float col_sum=F77NAME(sasum)(k,A.addr(j,j),inc);
    col_sum+=F77NAME(sasum)(j,A.addr(j,0),n);
    if (col_sum>anorm) anorm=col_sum;
  }
  SymmetricPositiveMatrix<float,float> *Acopy=
    OPERATOR_NEW SymmetricPositiveMatrix<float,float>(n);
  Acopy->copy(A);
  int info;
  char uplo='L';
  F77NAME(spotrf)(uplo,n,Acopy->addr(),n,info);
  CHECK_SAME(info,0)
  L=OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(n,n);
  for (int j=0;j<n;j++) {
    for (int i=j;i<n;i++) (*L)(i,j)=(*Acopy)(i,j);
  }
  delete Acopy;
}

template<> Matrix<float,float>*
CholeskyFactorization<float,float>::solveAXeqB(
const Matrix<float,float> &B) const {
//TRACER_CALL(t,"CholeskyFactorization::solveAXeqB");
  int n=L->size(0); 
  CHECK_SAME(n,B.size(0))
  int k=B.size(1);

  Matrix<float,float> *X=OPERATOR_NEW Matrix<float,float>(n,k);
  X->copy(B);
  char uplo='L';
  int info;
  F77NAME(spotrs)(uplo,n,k,L->addr(),n,X->addr(),n,info);
  CHECK_SAME(info,0)
  return X;
}

template<> float
CholeskyFactorization<float,float>::conditionNumber() const {
  int n=L->size(0);
  float rcond=float_zero_;
  float *work=OPERATOR_NEW_BRACKET(float,3*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  char uplo='L';
  F77NAME(spocon)(uplo,n,L->addr(),n,anorm,rcond,work,iwork,
                  info);
  CHECK_SAME(info,0)
  delete [] iwork;
  delete [] work;
  if (rcond>float_zero_) return float_one_/rcond;
  return HUGE_VAL;
}

template class CholeskyFactorization<float,float>;
//template ostream& operator<<(ostream&,
//  const CholeskyFactorization<float,float>&);
template void testCholeskyFactorization(float,float);
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> MDMtFactorization<float,float>::MDMtFactorization(
const SymmetricMatrix<float,float> &A) : anorm(float_zero_) {
//TRACER_CALL(t,"MDMtFactorization::MDMtF");
  int n=A.size(0);
  ipiv=OPERATOR_NEW_BRACKET(int,n);
  int inc=1;
  for (int j=0;j<n;j++) {
    int k=n-j;
    float col_sum=F77NAME(sasum)(k,A.addr(j,j),inc);
    col_sum+=F77NAME(sasum)(j,A.addr(j,0),n);
    if (col_sum>anorm) anorm=col_sum;
  }
  SymmetricMatrix<float,float> *Acopy=
    OPERATOR_NEW SymmetricMatrix<float,float>(n);
  Acopy->copy(A);
  int info;
  char uplo='L';
  float w=HUGE_VAL;
  int lwork=-1;
  F77NAME(ssytrf)(uplo,n,Acopy->addr(),n,ipiv,&w,lwork,info);
  CHECK_SAME(info,0)

  lwork=static_cast<int>(w);
  float *work=OPERATOR_NEW_BRACKET(float,lwork);
  F77NAME(ssytrf)(uplo,n,Acopy->addr(),n,ipiv,work,lwork,info);
  CHECK_SAME(info,0)
  delete [] work;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<float,float>(n,n);
  for (int j=0;j<n;j++) {
    for (int i=j;i<n;i++) (*L)(i,j)=(*Acopy)(i,j);
  }
  delete Acopy;
}

template<> Matrix<float,float>*
MDMtFactorization<float,float>::solveAXeqB(
const Matrix<float,float> &B) const {
//TRACER_CALL(t,"MDMtFactorization::solveAXeqB");
  int n=L->size(0);
  CHECK_SAME(n,B.size(0))
  int k=B.size(1);

  Matrix<float,float> *X=OPERATOR_NEW Matrix<float,float>(n,k);
  X->copy(B);
  char uplo='L';
  int info;
  F77NAME(ssytrs)(uplo,n,k,L->addr(),n,ipiv,X->addr(),n,info);
  CHECK_SAME(info,0)
  return X;
}

template<> float MDMtFactorization<float,float>::conditionNumber()
const {
  int n=L->size(0);
  float rcond=float_zero_;
  float *work=OPERATOR_NEW_BRACKET(float,2*n);
  int *iwork=OPERATOR_NEW_BRACKET(int,n);
  int info=0;
  char uplo='L';
  F77NAME(ssycon)(uplo,n,L->addr(),n,ipiv,anorm,rcond,work,iwork,info);
  CHECK_SAME(info,0)
  delete [] iwork;
  delete [] work;
  if (rcond>float_zero_) return float_one_/rcond;
  return HUGE_VAL;
}

template class MDMtFactorization<float,float>;
//template ostream& operator<<(ostream&,
//  const MDMtFactorization<float,float>&);
template void testMDMtFactorization(float,float);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<>
SingularValueDecomposition<float,float>::SingularValueDecomposition(
const Matrix<float,float> &A) {
//TRACER_CALL(t,"SingularValueDecomposition::SVD");
  int m=A.size(0),n=A.size(1),minmn=min(m,n);
  Sigma=OPERATOR_NEW Vector<float,float>(minmn);
  int info=0;
  float *work=0;
  int lwork=5*max(m,n);
  if (m>=n) {
    U=OPERATOR_NEW Matrix<float,float>(m,n);
    U->copy(A);
    V_transpose=OPERATOR_NEW Matrix<float,float>(n,n);
    char jobu='O';
    char jobvt='S';
//  int lwork=max(1,max( 3*n + m, 5*n - 4));
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sgesvd)(jobu,jobvt,m,n,U->addr(),m,Sigma->addr(),
      0,m,V_transpose->addr(),n,work,lwork,info);
  } else {
    U=OPERATOR_NEW Matrix<float,float>(m,m);
    V_transpose=OPERATOR_NEW Matrix<float,float>(m,n);
    V_transpose->copy(A);
    char jobu='S';
    char jobvt='O';
//  int lwork=max(1,max( 3*m + n, 5*m - 4));
    work=OPERATOR_NEW_BRACKET(float,lwork);
    F77NAME(sgesvd)(jobu,jobvt,m,n,V_transpose->addr(),m,
      Sigma->addr(),U->addr(),m,0,m,work,lwork,info);
  }
  CHECK_SAME(info,0)
  if (work) delete work;
}

template<> float SingularValueDecomposition<float,float>::conditionNumber(
float cutoff) const {
  CHECK_POSITIVE(cutoff)
  int r=rank(cutoff);
  if (r<1) return HUGE_VAL;
  if (r==1) return float_one_;
  return (*Sigma)[0]/(*Sigma)[r-1];
}

template<> Matrix<float,float>* SingularValueDecomposition<float,float>::solveAXeqB(
const Matrix<float,float> &B,float cutoff,transpose_option option) const {
//TRACER_CALL(t,"SingularValueDecomposition::solveAXeqB");
  int m=U->size(0),n=V_transpose->size(1),minmn=Sigma->size();
  int K=B.size(1);
  int r=rank(cutoff);
  float alpha=float_one_;
  float beta=float_zero_;
  char trans='N';
  switch (option) {
    case no_matrix_transpose: {
      CHECK_SAME(m,B.size(0))
      Matrix<float,float> C(r,K);
      char transa='T';
      F77NAME(sgemm)(transa,trans,r,K,m,alpha,U->addr(),m,
        B.addr(),m,beta,C.addr(),r);
      for (int i=0;i<r;i++) {
        float temp=1./(*Sigma)[i];
        F77NAME(sscal)(K,temp,C.addr(i,0),r);
      }
      Matrix<float,float> *X=OPERATOR_NEW Matrix<float,float>(n,K);
      F77NAME(sgemm)(transa,trans,n,K,r,alpha,
        V_transpose->addr(),minmn,C.addr(),r,beta,X->addr(),n);
      return X;
      break;
    }
    case matrix_transpose: {
      CHECK_SAME(n,B.size(0))
      Matrix<float,float> C(r,K);
      F77NAME(sgemm)(trans,trans,r,K,n,alpha,V_transpose->addr(),
        n,B.addr(),n,beta,C.addr(),r);
      for (int i=0;i<r;i++) {
        float temp=1./(*Sigma)[i];
        F77NAME(sscal)(K,temp,C.addr(i,0),r);
      }
      Matrix<float,float> *X=OPERATOR_NEW Matrix<float,float>(m,K);
      F77NAME(sgemm)(trans,trans,m,K,r,alpha,U->addr(),m,
        C.addr(),r,beta,X->addr(),m);
      return X;
      break;
    }
    default:
      return 0;
  }
}

template class SingularValueDecomposition<float,float>;
//template ostream& operator<<(ostream&,
//  const SingularValueDecomposition<float,float>&);
template void testSingularValueDecomposition(float,float);
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//#include "LinearProgram.C"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//template<> float VirtualLinearProgram<float>::zero_ = float_zero_;
//template<> float VirtualLinearProgram<float>::one_ = float_one_;
//template<> float VirtualLinearProgram<float>::huge_ = HUGE_VAL;

//template class VirtualLinearProgram<float>;
//template ostream& operator<<(ostream&,
//  const VirtualLinearProgram<float>&);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> int SFLinearProgram<float>::smallestDualSlack() {
//TRACER_CALL(t,"SFLinearProgram::smallestDualSlack");
  int M=A->size(0),N=A->size(1);
  int jmin=0;
  float rmin=float_zero_;
  for (int j=0;j<N-M;j++) {
    int colj=column[j+M];
    float rj=(*c)[colj]-F77NAME(sdot)(M,y->addr(),1,A->addr(0,colj),1);
    if (rj < rmin) { jmin=j; rmin=rj; }
    (*r)[j]=rj;
  }
#ifdef DEBUG
//cout << "\tr = " << endl;
//r->printOn(cout);
#endif
  return F77NAME(ismin)(N-M,r->addr(),1)-1;
}
*/

/*
//if perturbation in single component c_i between lower[i] and upper[i],
//then optimal x unchanged
template<> void SFLinearProgram<float>::costBounds(
Vector<float,float> &lower,Vector<float,float> &upper) const {
  int M=A->size(0),N=A->size(1);
  CHECK_SAME(N,lower.size())
  CHECK_SAME(N,upper.size())

//perturbation in basic variable
  for (int i=0;i<M;i++) {
    float low=-huge_;
    float high=huge_;
    Vector<float,float> axisi(M,this->zero_);
    axisi[i]=this->one_;
    Vector<float,float> hi(M);
    QR->solveUnderdetermined(axisi,hi);
    for (int j=M;j<N;j++) {
      int colj=column[j];
      float denom=F77NAME(sdot)(M,hi.addr(),1,A->addr(0,colj),1);
      float rj=(*r)[j-M];
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
template<> void SFLinearProgram<float>::arrayBounds(int i,int j,
float &lower,float &upper) const {
//TRACER_CALL(t,"SFLinearProgram::arrayBounds");
  int m=A->size(0),n=A->size(1);
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
    if (column[jj]==j) {
      float ainv_ji=(*basic_inverse)(jj,i);
      if (fabs(ainv_ji)>sfmin) {
	float bound=-float_one_/ainv_ji;
	if (ainv_ji>float_zero_) lower=bound;
	else upper=bound;
      }
      for (int k=0;k<m;k++) {
        if (k!=jj) {
	  float xk=(*x)(column[k],0);
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
	float term=F77NAME(sdot)(m,basic_inverse->addr(jj,0),m,
			   A->addr(0,column[kk]),incA);
	float denom=yi*term-ainv_ji*rk;
	if (fabs(denom)>rk*sfmin) {
	  float bound=rk/denom;
	  if (denom>float_zero_ && bound<upper) upper=bound;
	  if (denom<float_zero_ && bound>lower) lower=bound;
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
      float rj=(*r_trans)(0,jj-m);
      if (fabs(yi)>rj*sfmin) {
	float bound=(*A)(i,j)+rj/yi;
	if (yi>float_zero_ && bound<upper) upper=bound;
	if (yi<float_zero_ && bound>lower) bound=lower;
      }
    }
  } }
}
*/
//template class SFLinearProgram<float>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
template<> LinearProgram<float>::LinearProgram(
const Matrix<float,float> &Ain,const Vector<float,float> &bin,
const Vector<float,float> &cin) : 
VirtualLinearProgram<float>(Ain,bin,cin),column(0),row(0),QR(0),
bbasic(0),cbasic(0),xbasic(0),ybasic(0),r(0),s(0),basic_number(0) {
//TRACER_CALL(t,"LinearProgram::LinearProgram");
  int M=A->size(0),N=A->size(1);
#ifdef DEBUG
//cout << "\tA = " << endl;
//Ain.printOn(cout);
//cout << "\tb = " << endl;
//bin.printOn(cout);
//cout << "\tissumn(b) = " << F77NAME(issumn)(M,bin.addr(),1) << endl;
//cout << "\tc = " << endl;
//cin.printOn(cout);
//cout << "\tissump(c) = " << F77NAME(issump)(N,cin.addr(),1) << endl;
#endif

  row=OPERATOR_NEW_BRACKET(int,M);
  for (int i=0;i<M;i++) row[i]=i;
  column=OPERATOR_NEW_BRACKET(int,N);
  for (int j=0;j<N;j++) column[j]=j;
  if (F77NAME(issumn)(M,bin.addr(),1)>=M) {
    current_status=PRIMAL_FEASIBLE; // x=0 is feasible, all nonbasic
  }
  if (F77NAME(issump)(N,cin.addr(),1)>=N) {
    current_status=DUAL_FEASIBLE; // y=0 is feasible, all nonbasic
  }
#ifdef DEBUG
//cout << "\tcurrent_status = " << current_status << endl;
#endif
  if (current_status!=UNKNOWN) computeSolution();
}
*/

/*
template<> void LinearProgram<float>::specifyBasicVariables(
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
  SquareMatrix<float,float> Abasic(basic_no);
  bbasic=OPERATOR_NEW Vector<float,float>(basic_no);
  cbasic=OPERATOR_NEW Vector<float,float>(basic_no);
  for (int j=0;j<basic_no;j++) {
    (*bbasic)[j]=(*this->b)[basic_row_index[j]];
    (*cbasic)[j]=(*this->c)[basic_column_index[j]];
    for (int i=0;i<basic_no;i++) {
      Abasic(i,j)=(*this->A)(basic_row_index[i],basic_column_index[j]);
    }
  }
  QR=OPERATOR_NEW GramSchmidtQRFactorization<float,float>(Abasic);

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
  if (F77NAME(issump)(basic_no,xbasic->addr(),1)>=basic_no) {
    if (basic_no>=M) {
      this->current_status=PRIMAL_FEASIBLE;
      return;
    }
    if (F77NAME(issump)(M-basic_no,s->addr(),1)>=M-basic_no) {
      this->current_status=PRIMAL_FEASIBLE;
      return;
    }
  }
  if (F77NAME(issump)(basic_no,ybasic->addr(),1)>=basic_no) {
    if (basic_no>=N) {
      this->current_status=DUAL_FEASIBLE;
      return;
    }
    if (F77NAME(issump)(N-basic_no,r->addr(),1)>=N-basic_no) {
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
template<> STATUS_OPTION LinearProgram<float>::primalSimplexStep() {
//TRACER_CALL(t,"LinearProgram::primalSimplexStep");
  if (this->current_status==OPTIMAL) return OPTIMAL;
  CHECK_TEST(this->current_status==PRIMAL_FEASIBLE);
#ifdef DEBUG
//printOn(cout);
#endif

//find largest feasible simplex step
  int M=this->A->size(0),N=this->A->size(1);
  int i_basic=-1,j_non_basic=-1;
  float rmin=this->huge_,ymin=this->huge_;
  if (basic_number>0) {
    i_basic=F77NAME(ismin)(basic_number,ybasic->addr(),1)-1;
    ymin=(*ybasic)[i_basic];
  }
  if (basic_number<N) {
    j_non_basic=F77NAME(ismin)(N-basic_number,r->addr(),1)-1;
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
    Vector<float,float> axis(basic_number,float_zero_);
    axis[i_basic]=-this->one_;
    Vector<float,float> gb(basic_number,this->zero_);
    Vector<float,float> resid(basic_number);
    QR->solveOverdetermined(axis,gb,resid);
    float gb_norm=gb.nrm2()*static_cast<float>(basic_number);
    float basic_ratio=float_undefined_;
    int j_basic=basicPrimalPivot(gb,basic_ratio);

    float non_basic_ratio=float_undefined_;
    int i_non_basic=-1;
    if (basic_number<M) {
      Vector<float,float> gn(M-basic_number);
      for (int i=basic_number;i<M;i++) {
        float gni=this->zero_;
        int rowi=row[i];
        float Ai_norm=float_zero_;
        for (int j=0;j<basic_number;j++) {
          float Aij=(*A)(rowi,column[j]);
          gni+=Aij*gb[j];
          Ai_norm+=Aij*Aij;
        }
        if (abs(gni)<gb_norm*sqrt(Ai_norm)*DBL_EPSILON) gni=float_zero_;
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
    Vector<float,float> *hb=0;
    int j_basic=-1,i_non_basic=-1;
    float basic_ratio=float_undefined_;
    float non_basic_ratio=float_undefined_;
    float hb_norm=float_zero_;
    if (basic_number>0) {
      Vector<float,float> rhs(basic_number);
      for (int i=0;i<basic_number;i++) {
        rhs[i]=(*A)(row[i],coljn);
      }
      hb=OPERATOR_NEW Vector<float,float>(basic_number);
      Vector<float,float> resid(basic_number);
      QR->solveOverdetermined(rhs,*hb,resid);
      j_basic=basicPrimalPivot(*hb,basic_ratio);
      hb_norm=hb->nrm2()*static_cast<float>(basic_number);
    }
    if (basic_number<M) {
      Vector<float,float> hn(M-basic_number);
      for (int i=basic_number;i<M;i++) {
        int rowi=row[i];
        float hni=-(*this->A)(rowi,coljn);
        float Ai_norm=float_zero_;
        for (int j=0;j<basic_number;j++) {
          float Aij=(*A)(rowi,column[j]);
          hni+=Aij*(*hb)[j];
          Ai_norm+=Aij*Aij;
        }
        if (abs(hni)<(abs((*A)(rowi,coljn))+hb_norm*Ai_norm)*DBL_EPSILON){
          hni=float_zero_;
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
template<> STATUS_OPTION LinearProgram<float>::dualSimplexStep() {
//TRACER_CALL(t,"LinearProgram::dualSimplexStep");
  if (current_status==OPTIMAL) return OPTIMAL;
  CHECK_TEST(current_status==DUAL_FEASIBLE);
#ifdef DEBUG
//printOn(cout);
#endif
  int M=A->size(0),N=A->size(1);
//find largest feasible simplex step
  float smin=float_undefined_,xmin=float_undefined_;
  int i_non_basic=-1,j_basic=-1;
  if (basic_number<M) {
    i_non_basic=F77NAME(ismin)(M-basic_number,s->addr(),1)-1;
    smin=(*s)[i_non_basic];
    i_non_basic+=basic_number;
  }
  if (basic_number>0) {
    j_basic=F77NAME(ismin)(basic_number,xbasic->addr(),1)-1;
    xmin=(*xbasic)[j_basic];
  }
#ifdef DEBUG
//cout << "\tj_basic,xmin = " << j_basic << " " << xmin << endl;
//cout << "\ti_non_basic,smin = " << i_non_basic << " " << smin << endl;
#endif
  if (min(xmin,smin)>=float_zero_) return current_status=OPTIMAL;

  if (xmin<smin) {
//  TRACER_CALL(t,"LinearProgram::dualSimplexStep xmin<smin");
    Vector<float,float> axis(basic_number,float_zero_);
    axis[j_basic]=float_one_;
    Vector<float,float> gb(basic_number);
    QR->solveUnderdetermined(axis,gb);
    float gb_norm=gb.nrm2()*static_cast<float>(basic_number);
    float basic_ratio=float_undefined_;
    int i_basic=basicDualPivot(gb,basic_ratio);

    float non_basic_ratio=float_undefined_;
    int j_non_basic=-1;
    if (basic_number<N) {
      Vector<float,float> gn(N-basic_number);
      for (int j=basic_number;j<N;j++) {
        float gnj=float_zero_;
        int colj=column[j];
        float Aj_norm=float_zero_;
        for (int i=0;i<basic_number;i++) {
          float Aij=(*A)(row[i],colj);
          gnj+=(gb)[i]*Aij;
          Aj_norm+=Aij*Aij;
        }
        if (abs(gnj)<gb_norm*sqrt(Aj_norm)*DBL_EPSILON) gnj=float_zero_;
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
    float basic_ratio=float_undefined_;
    float non_basic_ratio=float_undefined_;
    int i_basic=-1,j_non_basic=-1;
    int rowin=row[i_non_basic];
    Vector<float,float> *hb=0;
    float hb_norm=float_zero_;
    if (basic_number>0) {
      Vector<float,float> rhs(basic_number);
      for (int j=0;j<basic_number;j++) {
        rhs[j]=(*A)(rowin,column[j]);
      }
      hb=OPERATOR_NEW Vector<float,float>(basic_number);
      QR->solveUnderdetermined(rhs,*hb);
      i_basic=basicDualPivot(*hb,basic_ratio);
      hb_norm=hb->nrm2()*static_cast<float>(basic_number);
    }
    if (basic_number<N) {
      Vector<float,float> hn(N-basic_number);
      for (int j=basic_number;j<N;j++) {
        int colj=column[j];
        float hnj=(*A)(rowin,colj);
        float Aj_norm=float_zero_;
        for (int i=0;i<basic_number;i++) {
          float Aij=(*A)(row[i],colj);
          hnj-=(*hb)[i]*Aij;
          Aj_norm+=Aij*Aij;
        }
        if (abs(hnj)<(abs((*A)(rowin,colj))+hb_norm*Aj_norm)*DBL_EPSILON){
          hnj=float_zero_;
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
LinearProgram<float>::findPrimalBasicFeasibleGuess() {
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
  int ib=F77NAME(issumn)(m,b->addr(),1);
  Matrix<float,float> *temp_A=
    OPERATOR_NEW Matrix<float,float>(m,n+m-ib,float_zero_);
  temp_A->copyFrom('A',m,n,*A);
  Vector<float,float> *temp_c=
    OPERATOR_NEW Vector<float,float>(n+m-ib,float_zero_);
  int pos_b=n;
  for (int i=0;i<m;i++) {
    if ((*b)[i]>float_zero_) {
      (*temp_A)(i,pos_b)=float_one_;
      (*temp_c)[pos_b]=float_one_;
      pos_b++;
    }
  }
  LinearProgram LP(*temp_A,*b,*temp_c);
//int *basic_row_index=OPERATOR_NEW_BRACKET(int,m-ib);
//int ip=0;
//for (int i=0;i<m;i++) {
//  if ((*b)[i]>float_zero_) basic_row_index[ip++]=i;
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
  QR=OPERATOR_NEW GramSchmidtQRFactorization<float,float>(*LP.QR);

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
LinearProgram<float>::findDualBasicFeasibleGuess() {
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
  int ic=F77NAME(issump)(n,c->addr(),1);
  Matrix<float,float> *temp_A=
    OPERATOR_NEW Matrix<float,float>(m+n-ic,n,float_zero_);
  temp_A->copyFrom('A',m,n,*A);
  Vector<float,float> *temp_b=
    OPERATOR_NEW Vector<float,float>(m+n-ic,float_zero_);
  int neg_c=m;
  for (int j=0;j<n;j++) {
    if ((*c)[j]<float_zero_) {
      (*temp_A)(neg_c,j)=-float_one_;
      (*temp_b)[neg_c]=-float_one_;
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
//  if ((*c)[j]<float_zero_) basic_column_index[jn++]=j;
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
  QR=OPERATOR_NEW GramSchmidtQRFactorization<float,float>(*LP.QR);

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

//template class LinearProgram<float>;
//template LinearProgram<float>::LinearProgram(const Matrix<float,float>&,const Matrix<float,float>&,const Matrix<float,float>&);
//template float LinearProgram<float>::currentValue() const;
//template status_option LinearProgram<float>::simplexStep();
//template void LinearProgram<float>::costBounds(Matrix<float,float>&,
//   Matrix<float,float>&) const;
//template void LinearProgram<float>::constraintBounds(Matrix<float,float>&,
//   Matrix<float,float>&) const;
//template void LinearProgram<float>::printOn(ostream&) const;

/*
template<> STATUS_OPTION LinearProgram<float>::findBasicFeasibleGuess() {
//TRACER_CALL(t,"LinearProgram::findBasicFeasibleGuess");
//if we have already dinked with this program, get out now
  CHECK_TEST(current_status!=INFEASIBLE && current_status!=UNBOUNDED);
  if (current_status!=UNKNOWN) return current_status;
  int m=A->size(0),n=A->size(1);
  int ib=F77NAME(issumn)(m,b->addr(),1);
  int ic=F77NAME(issump)(n,c->addr(),1);
#ifdef DEBUG
//cout << "ib,ic = " << ib << " " << ic << endl;
#endif
  if (m-ib<=n-ic) return findPrimalBasicFeasibleGuess();
  else return findDualBasicFeasibleGuess();
}
*/
