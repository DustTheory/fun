/* dgbconnp.f -- translated by f2c (version 20090411).
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

/* modified from dgbcon to avoid pivoting */
/* Subroutine */ int dgbconnp_(char *norm, integer *n, integer *kl, integer *
	ku, doublereal *ab, integer *ldab, doublereal *anorm, doublereal *
	rcond, doublereal *work, integer *iwork, integer *info, ftnlen 
	norm_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t;
    static integer kd, lm, ix, kase;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lnoti;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dlacn2_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    static logical onenrm;
    static char normin[1];
    static doublereal smlnum;

    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    --work;
    --iwork;

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
    } else if (*anorm < 0.) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("dgbconnp", &i__1, (ftnlen)8);
	return 0;
    }
    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return 0;
    } else if (*anorm == 0.) {
	return 0;
    }
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
    ainvnm = 0.;
    *(unsigned char *)normin = 'N';
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kd = *ku + 1;
    lnoti = *kl > 0;
    kase = 0;
L10:
    dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == kase1) {
	    if (lnoti) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		    i__2 = *kl, i__3 = *n - j;
		    lm = min(i__2,i__3);
		    t = work[j];
		    d__1 = -t;
		    daxpy_(&lm, &d__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		}
	    }
/* lnoti */
	    dlatbs_("Upper", "No transpose", "Non-unit", normin, n, ku, &ab[
		    ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
	} else {
/* kase.ne.kase1 */
	    dlatbs_("Upper", "Transpose", "non-unit", normin, n, ku, &ab[
		    ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);
	    if (lnoti) {
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		    i__1 = *kl, i__2 = *n - j;
		    lm = min(i__1,i__2);
		    work[j] -= ddot_(&lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		}
	    }
/* lnoti */
	}
	*(unsigned char *)normin = 'Y';
	if (scale != 1.) {
	    ix = idamax_(n, &work[1], &c__1);
	    if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
		    {
		goto L40;
	    }
	    drscl_(n, &scale, &work[1], &c__1);
	}
	goto L10;
    }
    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }
L40:
    return 0;
} /* dgbconnp_ */

