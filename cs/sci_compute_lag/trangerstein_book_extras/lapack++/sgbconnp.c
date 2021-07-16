/* sgbconnp.f -- translated by f2c (version 20090411).
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

/* modified from sgbcon to avoid pivoting */
/* Subroutine */ int sgbconnp_(char *norm, integer *n, integer *kl, integer *
	ku, real *ab, integer *ldab, real *anorm, real *rcond, real *work, 
	integer *iwork, integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer j;
    static real t;
    static integer kd, lm, ix, kase;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static integer kase1;
    static real scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    static logical lnoti;
    extern /* Subroutine */ int srscl_(integer *, real *, real *, integer *), 
	    saxpy_(integer *, real *, real *, integer *, real *, integer *), 
	    slacn2_(integer *, real *, real *, integer *, real *, integer *, 
	    integer *);
    extern doublereal slamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, real *, integer *);
    static real ainvnm;
    extern /* Subroutine */ int slatbs_(char *, char *, char *, char *, 
	    integer *, integer *, real *, integer *, real *, real *, real *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical onenrm;
    static char normin[1];
    static real smlnum;

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
    } else if (*anorm < 0.f) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("sgbconnp", &i__1, (ftnlen)8);
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
    slacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == kase1) {
	    if (lnoti) {
		i__1 = *n - 1;
		for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		    i__2 = *kl, i__3 = *n - j;
		    lm = min(i__2,i__3);
		    t = work[j];
		    r__1 = -t;
		    saxpy_(&lm, &r__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		}
	    }
/* lnoti */
	    slatbs_("Upper", "No transpose", "Non-unit", normin, n, ku, &ab[
		    ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
	} else {
/* kase.ne.kase1 */
	    slatbs_("Upper", "Transpose", "non-unit", normin, n, ku, &ab[
		    ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);
	    if (lnoti) {
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		    i__1 = *kl, i__2 = *n - j;
		    lm = min(i__1,i__2);
		    work[j] -= sdot_(&lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		}
	    }
/* lnoti */
	}
	*(unsigned char *)normin = 'Y';
	if (scale != 1.f) {
	    ix = isamax_(n, &work[1], &c__1);
	    if (scale < (r__1 = work[ix], dabs(r__1)) * smlnum || scale == 
		    0.f) {
		goto L40;
	    }
	    srscl_(n, &scale, &work[1], &c__1);
	}
	goto L10;
    }
    if (ainvnm != 0.f) {
	*rcond = 1.f / ainvnm / *anorm;
    }
L40:
    return 0;
} /* sgbconnp_ */

