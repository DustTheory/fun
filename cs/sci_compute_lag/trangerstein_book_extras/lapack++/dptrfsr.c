/* dptrfsr.f -- translated by f2c (version 20090411).
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
static doublereal c_b7 = 1.;

/* modified from dptrfs to improve X A=B <==> A X^T=B^T */
/* Subroutine */ int dptrfsr_(integer *n, integer *nrhs, doublereal *d__, 
	doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer 
	*ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	 doublereal *work, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j;
    static doublereal s, bi, cx, dx, ex;
    static integer ix, nz;
    static doublereal eps, safe1, safe2;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer count;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int dpttrs_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *);

    /* Parameter adjustments */
    --d__;
    --e;
    --df;
    --ef;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*ldb < max(1,*nrhs)) {
	*info = -8;
    } else if (*ldx < max(1,*nrhs)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("dptrfsr", &i__1, (ftnlen)7);
	return 0;
    }
    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    ferr[j] = 0.;
	    berr[j] = 0.;
	}
	return 0;
    }
    nz = 4;
    eps = dlamch_("Epsilon", (ftnlen)7);
    safmin = dlamch_("Safe minimum", (ftnlen)12);
    safe1 = nz * safmin;
    safe2 = safe1 / eps;
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
	count = 1;
	lstres = 3.;
L20:
	if (*n == 1) {
	    bi = b[j + b_dim1];
	    dx = d__[1] * x[j + x_dim1];
	    work[*n + 1] = bi - dx;
	    work[1] = abs(bi) + abs(dx);
	} else {
	    bi = b[j + b_dim1];
	    dx = d__[1] * x[j + x_dim1];
	    ex = e[1] * x[j + (x_dim1 << 1)];
	    work[*n + 1] = bi - dx - ex;
	    work[1] = abs(bi) + abs(dx) + abs(ex);
	    i__2 = *n - 1;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		bi = b[j + i__ * b_dim1];
		cx = e[i__ - 1] * x[j + (i__ - 1) * x_dim1];
		dx = d__[i__] * x[j + i__ * x_dim1];
		ex = e[i__] * x[j + (i__ + 1) * x_dim1];
		work[*n + i__] = bi - cx - dx - ex;
		work[i__] = abs(bi) + abs(cx) + abs(dx) + abs(ex);
	    }
	    bi = b[j + *n * b_dim1];
	    cx = e[*n - 1] * x[j + (*n - 1) * x_dim1];
	    dx = d__[*n] * x[j + *n * x_dim1];
	    work[*n + *n] = bi - cx - dx;
	    work[*n] = abs(bi) + abs(cx) + abs(dx);
	}
	s = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (work[i__] > safe2) {
/* Computing MAX */
		d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
			i__];
		s = max(d__2,d__3);
	    } else {
/* Computing MAX */
		d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
			/ (work[i__] + safe1);
		s = max(d__2,d__3);
	    }
	}
	berr[j] = s;
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {
	    dpttrs_(n, &c__1, &df[1], &ef[1], &work[*n + 1], n, info);
	    daxpy_(n, &c_b7, &work[*n + 1], &c__1, &x[j + x_dim1], ldx);
	    lstres = berr[j];
	    ++count;
	    goto L20;
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (work[i__] > safe2) {
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__];
	    } else {
		work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
			work[i__] + safe1;
	    }
	}
	ix = idamax_(n, &work[1], &c__1);
	ferr[j] = work[ix];
	work[1] = 1.;
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    work[i__] = work[i__ - 1] * (d__1 = ef[i__ - 1], abs(d__1)) + 1.;
	}
	work[*n] /= df[*n];
	for (i__ = *n - 1; i__ >= 1; --i__) {
	    work[i__] = work[i__] / df[i__] + work[i__ + 1] * (d__1 = ef[i__],
		     abs(d__1));
	}
	ix = idamax_(n, &work[1], &c__1);
	ferr[j] *= (d__1 = work[ix], abs(d__1));
	lstres = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__2 = lstres, d__3 = (d__1 = x[j + i__ * x_dim1], abs(d__1));
	    lstres = max(d__2,d__3);
	}
	if (lstres != 0.) {
	    ferr[j] /= lstres;
	}
    }
    return 0;
} /* dptrfsr_ */

