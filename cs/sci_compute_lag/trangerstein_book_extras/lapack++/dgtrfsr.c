/* dgtrfsr.f -- translated by f2c (version 20090411).
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
static doublereal c_b14 = 1.;

/*     modified from dgtrfs to multiply by matrix on the right */
/* Subroutine */ int dgtrfsr_(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf, 
	doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv, 
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer count;
    extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static char transn[1];
    extern /* Subroutine */ int dgttrs_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen);
    static char transt[1];
    static doublereal lstres;

    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    --dlf;
    --df;
    --duf;
    --du2;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;
    --iwork;

    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(1,*nrhs)) {
	*info = -13;
    } else if (*ldx < max(1,*nrhs)) {
	*info = -15;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("dgtrfsr", &i__1, (ftnlen)7);
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
    if (! notran) {
	*(unsigned char *)transn = 'N';
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transn = 'T';
	*(unsigned char *)transt = 'N';
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
	if (! notran) {
/* X A^T = B <==> A X^T = B^T */
	    if (*n == 1) {
		work[1] = (d__1 = b[j + b_dim1], abs(d__1)) + (d__2 = d__[1] *
			 x[j + x_dim1], abs(d__2));
		work[*n + 1] = b[j + b_dim1] - d__[1] * x[j + x_dim1];
	    } else {
		work[1] = (d__1 = b[j + b_dim1], abs(d__1)) + (d__2 = d__[1] *
			 x[j + x_dim1], abs(d__2)) + (d__3 = du[1] * x[j + (
			x_dim1 << 1)], abs(d__3));
		work[*n + 1] = b[j + b_dim1] - d__[1] * x[j + x_dim1] - du[1] 
			* x[j + (x_dim1 << 1)];
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    work[i__] = (d__1 = b[j + i__ * b_dim1], abs(d__1)) + (
			    d__2 = dl[i__ - 1] * x[j + (i__ - 1) * x_dim1], 
			    abs(d__2)) + (d__3 = d__[i__] * x[j + i__ * 
			    x_dim1], abs(d__3)) + (d__4 = du[i__] * x[j + (
			    i__ + 1) * x_dim1], abs(d__4));
		    work[*n + i__] = b[j + i__ * b_dim1] - dl[i__ - 1] * x[j 
			    + (i__ - 1) * x_dim1] - d__[i__] * x[j + i__ * 
			    x_dim1] - du[i__] * x[j + (i__ + 1) * x_dim1];
		}
		work[*n] = (d__1 = b[j + *n * b_dim1], abs(d__1)) + (d__2 = 
			dl[*n - 1] * x[j + (*n - 1) * x_dim1], abs(d__2)) + (
			d__3 = d__[*n] * x[j + *n * x_dim1], abs(d__3));
		work[*n + *n] = b[j + *n * b_dim1] - dl[*n - 1] * x[j + (*n - 
			1) * x_dim1] - d__[*n] * x[j + *n * x_dim1];
	    }
	} else {
/* X A = B <==> A^T X^T = B^T */
	    if (*n == 1) {
		work[1] = (d__1 = b[j + b_dim1], abs(d__1)) + (d__2 = d__[1] *
			 x[j + x_dim1], abs(d__2));
		work[*n + 1] = b[j + b_dim1] - d__[1] * x[j + x_dim1];
	    } else {
		work[1] = (d__1 = b[j + b_dim1], abs(d__1)) + (d__2 = d__[1] *
			 x[j + x_dim1], abs(d__2)) + (d__3 = dl[1] * x[j + (
			x_dim1 << 1)], abs(d__3));
		work[*n + 1] = b[j + b_dim1] - d__[1] * x[j + x_dim1] - dl[1] 
			* x[j + (x_dim1 << 1)];
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    work[i__] = (d__1 = b[j + i__ * b_dim1], abs(d__1)) + (
			    d__2 = du[i__ - 1] * x[j + (i__ - 1) * x_dim1], 
			    abs(d__2)) + (d__3 = d__[i__] * x[j + i__ * 
			    x_dim1], abs(d__3)) + (d__4 = dl[i__] * x[j + (
			    i__ + 1) * x_dim1], abs(d__4));
		    work[*n + i__] = b[j + i__ * b_dim1] - du[i__ - 1] * x[j 
			    + (i__ - 1) * x_dim1] - d__[i__] * x[j + i__ * 
			    x_dim1] - dl[i__] * x[j + (i__ + 1) * x_dim1];
		}
		work[*n] = (d__1 = b[j + *n * b_dim1], abs(d__1)) + (d__2 = 
			du[*n - 1] * x[j + (*n - 1) * x_dim1], abs(d__2)) + (
			d__3 = d__[*n] * x[j + *n * x_dim1], abs(d__3));
		work[*n + *n] = b[j + *n * b_dim1] - du[*n - 1] * x[j + (*n - 
			1) * x_dim1] - d__[*n] * x[j + *n * x_dim1];
	    }
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
	    dgttrs_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
		    ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
	    daxpy_(n, &c_b14, &work[*n + 1], &c__1, &x[j + x_dim1], ldx);
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
	kase = 0;
L70:
	dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
		kase, isave);
	if (kase != 0) {
	    if (kase == 1) {
		dgttrs_(transt, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
			ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    work[*n + i__] = work[i__] * work[*n + i__];
		}
	    } else {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    work[*n + i__] = work[i__] * work[*n + i__];
		}
		dgttrs_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &du2[1], &
			ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
	    }
	    goto L70;
	}
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
} /* dgtrfsr_ */

