/* zgtconnp.f -- translated by f2c (version 20090411).
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

/* Subroutine */ int zgtconnp_(char *norm, integer *n, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublereal *anorm, doublereal *
	rcond, doublecomplex *work, integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int zgttrsnp_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, integer *, ftnlen);
    static integer kase, kase1;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int zlacn2_(integer *, doublecomplex *, 
	    doublecomplex *, doublereal *, integer *, integer *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal ainvnm;
    static logical onenrm;

    /* Parameter adjustments */
    --work;
    --du;
    --d__;
    --dl;

    /* Function Body */
    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*anorm < 0.) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("zgtconnp", &i__1, (ftnlen)8);
	return 0;
    }
    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return 0;
    } else if (*anorm == 0.) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	if (d__[i__2].r == 0. && d__[i__2].i == 0.) {
	    return 0;
	}
    }
    ainvnm = 0.;
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kase = 0;
L20:
    zlacn2_(n, &work[*n + 1], &work[1], &ainvnm, &kase, isave);
    if (kase != 0) {
	if (kase == kase1) {
	    zgttrsnp_("No transpose", n, &c__1, &dl[1], &d__[1], &du[1], &
		    work[1], n, info, (ftnlen)12);
	} else {
	    zgttrsnp_("Conjugate transpose", n, &c__1, &dl[1], &d__[1], &du[1]
		    , &work[1], n, info, (ftnlen)19);
	}
	goto L20;
    }
    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }
    return 0;
} /* zgtconnp_ */

