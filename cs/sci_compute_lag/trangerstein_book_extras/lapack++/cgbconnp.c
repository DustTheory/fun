/* cgbconnp.f -- translated by f2c (version 20090411).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* modified from cgbcon to avoid pivoting */
/* Subroutine */ int cgbconnp_(char *norm, integer *n, integer *kl, integer *
	ku, complex *ab, integer *ldab, real *anorm, real *rcond, complex *
	work, real *rwork, integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    real r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double r_imag(complex *);

    /* Local variables */
    static integer j;
    static complex t;
    static integer kd, lm, ix, kase, kase1;
    static real scale;
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static logical lnoti;
    extern /* Subroutine */ int clacn2_(integer *, complex *, complex *, real 
	    *, integer *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int clatbs_(char *, char *, char *, char *, 
	    integer *, integer *, complex *, integer *, complex *, real *, 
	    real *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), xerbla_(char *
	    , integer *, ftnlen);
    static real ainvnm;
    extern integer izamax_(integer *, complex *, integer *);
    extern /* Subroutine */ int csrscl_(integer *, real *, complex *, integer 
	    *);
    static logical onenrm;
    static char normin[1];
    static real smlnum;

    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --work;
    --rwork;

    /* Function Body */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*ldab < *kl + *ku + 1) {
	*info = -6;
    } else if (*anorm < 0.f) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("cgbconnp", &i__1, (ftnlen)8);
	return 0;
    }
    *rcond = 0.f;
    if (*n == 0) {
	*rcond = 1.f;
	return 0;
    } else if (*anorm == 0.f) {
	return 0;
    }
    smlnum = slamch_("Safe minimum", (ftnlen)12);
    ainvnm = 0.f;
    *(unsigned char *)normin = 'n';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kd = *ku + 1;
    lnoti = *kl > 0;
    kase = 0;
L10:
    clacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == kase1) {
	    if (lnoti) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		    i__2 = *kl, i__3 = *n - j;
		    lm = min(i__2,i__3);
		    i__2 = j;
		    t.r = work[i__2].r, t.i = work[i__2].i;
		    q__1.r = -t.r, q__1.i = -t.i;
		    caxpy_(&lm, &q__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		}
	    }
	    clatbs_("Upper", "No transpose", "Non-unit", normin, n, ku, &ab[
		    ab_offset], ldab, &work[1], &scale, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
	} else {
	    clatbs_("Upper", "Conjugate transpose", "Non-unit", normin, n, ku,
		     &ab[ab_offset], ldab, &work[1], &scale, &rwork[1], info, 
		    (ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);
	    if (lnoti) {
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		    i__1 = *kl, i__2 = *n - j;
		    lm = min(i__1,i__2);
		    i__1 = j;
		    i__2 = j;
		    cdotc_(&q__2, &lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		    q__1.r = work[i__2].r - q__2.r, q__1.i = work[i__2].i - 
			    q__2.i;
		    work[i__1].r = q__1.r, work[i__1].i = q__1.i;
		}
	    }
	}
	*(unsigned char *)normin = 'Y';
	if (scale != 1.f) {
	    ix = izamax_(n, &work[1], &c__1);
	    i__1 = ix;
	    if (scale < ((r__1 = work[i__1].r, dabs(r__1)) + (r__2 = r_imag(&
		    work[ix]), dabs(r__2))) * smlnum || scale == 0.f) {
		goto L40;
	    }
	    csrscl_(n, &scale, &work[1], &c__1);
	}
	goto L10;
    }
    if (ainvnm != 0.f) {
	*rcond = 1.f / ainvnm / *anorm;
    }
L40:
    return 0;
} /* cgbconnp_ */

