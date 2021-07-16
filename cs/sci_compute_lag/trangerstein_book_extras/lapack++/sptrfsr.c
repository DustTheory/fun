/* sptrfsr.f -- translated by f2c (version 20090411).
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
static real c_b7 = 1.f;

/* modified from sptrfs to improve X A=B <==> A X^T=B^T */
/* Subroutine */ int sptrfsr_(integer *n, integer *nrhs, real *d__, real *e, 
	real *df, real *ef, real *b, integer *ldb, real *x, integer *ldx, 
	real *ferr, real *berr, real *work, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__, j;
    static real s, bi, cx, dx, ex;
    static integer ix, nz;
    static real eps, safe1, safe2;
    static integer count;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);
    extern doublereal slamch_(char *, ftnlen);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, real *, integer *);
    static real lstres;
    extern /* Subroutine */ int spttrs_(integer *, integer *, real *, real *, 
	    real *, integer *, integer *);

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
	xerbla_("sptrfsr", &i__1, (ftnlen)7);
	return 0;
    }
    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    ferr[j] = 0.f;
	    berr[j] = 0.f;
	}
	return 0;
    }
    nz = 4;
    eps = slamch_("Epsilon", (ftnlen)7);
    safmin = slamch_("Safe minimum", (ftnlen)12);
    safe1 = nz * safmin;
    safe2 = safe1 / eps;
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
	count = 1;
	lstres = 3.f;
L20:
	if (*n == 1) {
	    bi = b[j + b_dim1];
	    dx = d__[1] * x[j + x_dim1];
	    work[*n + 1] = bi - dx;
	    work[1] = dabs(bi) + dabs(dx);
	} else {
	    bi = b[j + b_dim1];
	    dx = d__[1] * x[j + x_dim1];
	    ex = e[1] * x[j + (x_dim1 << 1)];
	    work[*n + 1] = bi - dx - ex;
	    work[1] = dabs(bi) + dabs(dx) + dabs(ex);
	    i__2 = *n - 1;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		bi = b[j + i__ * b_dim1];
		cx = e[i__ - 1] * x[j + (i__ - 1) * x_dim1];
		dx = d__[i__] * x[j + i__ * x_dim1];
		ex = e[i__] * x[j + (i__ + 1) * x_dim1];
		work[*n + i__] = bi - cx - dx - ex;
		work[i__] = dabs(bi) + dabs(cx) + dabs(dx) + dabs(ex);
	    }
	    bi = b[j + *n * b_dim1];
	    cx = e[*n - 1] * x[j + (*n - 1) * x_dim1];
	    dx = d__[*n] * x[j + *n * x_dim1];
	    work[*n + *n] = bi - cx - dx;
	    work[*n] = dabs(bi) + dabs(cx) + dabs(dx);
	}
	s = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (work[i__] > safe2) {
/* Computing MAX */
		r__2 = s, r__3 = (r__1 = work[*n + i__], dabs(r__1)) / work[
			i__];
		s = dmax(r__2,r__3);
	    } else {
/* Computing MAX */
		r__2 = s, r__3 = ((r__1 = work[*n + i__], dabs(r__1)) + safe1)
			 / (work[i__] + safe1);
		s = dmax(r__2,r__3);
	    }
	}
	berr[j] = s;
	if (berr[j] > eps && berr[j] * 2.f <= lstres && count <= 5) {
	    spttrs_(n, &c__1, &df[1], &ef[1], &work[*n + 1], n, info);
	    saxpy_(n, &c_b7, &work[*n + 1], &c__1, &x[j + x_dim1], ldx);
	    lstres = berr[j];
	    ++count;
	    goto L20;
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (work[i__] > safe2) {
		work[i__] = (r__1 = work[*n + i__], dabs(r__1)) + nz * eps * 
			work[i__];
	    } else {
		work[i__] = (r__1 = work[*n + i__], dabs(r__1)) + nz * eps * 
			work[i__] + safe1;
	    }
	}
	ix = isamax_(n, &work[1], &c__1);
	ferr[j] = work[ix];
	work[1] = 1.f;
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    work[i__] = work[i__ - 1] * (r__1 = ef[i__ - 1], dabs(r__1)) + 
		    1.f;
	}
	work[*n] /= df[*n];
	for (i__ = *n - 1; i__ >= 1; --i__) {
	    work[i__] = work[i__] / df[i__] + work[i__ + 1] * (r__1 = ef[i__],
		     dabs(r__1));
	}
	ix = isamax_(n, &work[1], &c__1);
	ferr[j] *= (r__1 = work[ix], dabs(r__1));
	lstres = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = lstres, r__3 = (r__1 = x[j + i__ * x_dim1], dabs(r__1));
	    lstres = dmax(r__2,r__3);
	}
	if (lstres != 0.f) {
	    ferr[j] /= lstres;
	}
    }
    return 0;
} /* sptrfsr_ */

