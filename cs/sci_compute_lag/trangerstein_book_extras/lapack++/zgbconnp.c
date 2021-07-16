/* zgbconnp.f -- translated by f2c (version 20090411).
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

/* modified from zgbcon to avoid pivoting */
/* Subroutine */ int zgbconnp_(char *norm, integer *n, integer *kl, integer *
	ku, doublecomplex *ab, integer *ldab, doublereal *anorm, doublereal *
	rcond, doublecomplex *work, doublereal *rwork, integer *info, ftnlen 
	norm_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer j;
    static doublecomplex t;
    static integer kd, lm, ix, kase, kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical lnoti;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zlacn2_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal ainvnm;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    static logical onenrm;
    extern /* Subroutine */ int zlatbs_(char *, char *, char *, char *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen), zdrscl_(integer *, doublereal *, doublecomplex *, 
	    integer *);
    static char normin[1];
    static doublereal smlnum;

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
    } else if (*anorm < 0.) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("zgbconnp", &i__1, (ftnlen)8);
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
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
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
		    z__1.r = -t.r, z__1.i = -t.i;
		    zaxpy_(&lm, &z__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		}
	    }
	    zlatbs_("Upper", "No transpose", "Non-unit", normin, n, ku, &ab[
		    ab_offset], ldab, &work[1], &scale, &rwork[1], info, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
	} else {
	    zlatbs_("Upper", "Conjugate transpose", "Non-unit", normin, n, ku,
		     &ab[ab_offset], ldab, &work[1], &scale, &rwork[1], info, 
		    (ftnlen)5, (ftnlen)19, (ftnlen)8, (ftnlen)1);
	    if (lnoti) {
		for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		    i__1 = *kl, i__2 = *n - j;
		    lm = min(i__1,i__2);
		    i__1 = j;
		    i__2 = j;
		    zdotc_(&z__2, &lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
			    work[j + 1], &c__1);
		    z__1.r = work[i__2].r - z__2.r, z__1.i = work[i__2].i - 
			    z__2.i;
		    work[i__1].r = z__1.r, work[i__1].i = z__1.i;
		}
	    }
	}
	*(unsigned char *)normin = 'Y';
	if (scale != 1.) {
	    ix = izamax_(n, &work[1], &c__1);
	    i__1 = ix;
	    if (scale < ((d__1 = work[i__1].r, abs(d__1)) + (d__2 = d_imag(&
		    work[ix]), abs(d__2))) * smlnum || scale == 0.) {
		goto L40;
	    }
	    zdrscl_(n, &scale, &work[1], &c__1);
	}
	goto L10;
    }
    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }
L40:
    return 0;
} /* zgbconnp_ */

