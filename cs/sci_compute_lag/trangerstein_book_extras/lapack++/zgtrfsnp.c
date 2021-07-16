/* zgtrfsnp.f -- translated by f2c (version 20090411).
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
static doublereal c_b16 = -1.;
static doublereal c_b17 = 1.;
static doublecomplex c_b21 = {1.,0.};

/* modified from zgtrfs to avoid pivoting */
/* Subroutine */ int zgtrfsnp_(char *trans, integer *n, integer *nrhs, 
	doublecomplex *dl, doublecomplex *d__, doublecomplex *du, 
	doublecomplex *dlf, doublecomplex *df, doublecomplex *duf, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int zgttrsnp_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, doublecomplex *
	    , integer *, integer *, ftnlen);
    static doublereal s;
    static integer nz;
    static doublereal eps;
    static integer kase;
    static doublereal safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3], count;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zlacn2_(
	    integer *, doublecomplex *, doublecomplex *, doublereal *, 
	    integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), zlagtm_(
	    char *, integer *, integer *, doublereal *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, doublecomplex *, integer *, 
	    doublereal *, doublecomplex *, integer *, ftnlen);
    static logical notran;
    static char transn[1], transt[1];
    static doublereal lstres;

    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    --dlf;
    --df;
    --duf;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --ferr;
    --berr;
    --work;
    --rwork;

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
    } else if (*ldb < max(1,*n)) {
	*info = -11;
    } else if (*ldx < max(1,*n)) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("zgtrfsnp", &i__1, (ftnlen)8);
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
    if (notran) {
	*(unsigned char *)transn = 'N';
	*(unsigned char *)transt = 'C';
    } else {
	*(unsigned char *)transn = 'C';
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
	zcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
	zlagtm_(trans, n, &c__1, &c_b16, &dl[1], &d__[1], &du[1], &x[j * 
		x_dim1 + 1], ldx, &c_b17, &work[1], n, (ftnlen)1);
	if (notran) {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		rwork[1] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			j * b_dim1 + 1]), abs(d__2)) + ((d__3 = d__[1].r, abs(
			d__3)) + (d__4 = d_imag(&d__[1]), abs(d__4))) * ((
			d__5 = x[i__3].r, abs(d__5)) + (d__6 = d_imag(&x[j * 
			x_dim1 + 1]), abs(d__6)));
	    } else {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		i__4 = j * x_dim1 + 2;
		rwork[1] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			j * b_dim1 + 1]), abs(d__2)) + ((d__3 = d__[1].r, abs(
			d__3)) + (d__4 = d_imag(&d__[1]), abs(d__4))) * ((
			d__5 = x[i__3].r, abs(d__5)) + (d__6 = d_imag(&x[j * 
			x_dim1 + 1]), abs(d__6))) + ((d__7 = du[1].r, abs(
			d__7)) + (d__8 = d_imag(&du[1]), abs(d__8))) * ((d__9 
			= x[i__4].r, abs(d__9)) + (d__10 = d_imag(&x[j * 
			x_dim1 + 2]), abs(d__10)));
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    i__3 = i__ + j * b_dim1;
		    i__4 = i__ - 1;
		    i__5 = i__ - 1 + j * x_dim1;
		    i__6 = i__;
		    i__7 = i__ + j * x_dim1;
		    i__8 = i__;
		    i__9 = i__ + 1 + j * x_dim1;
		    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = 
			    d_imag(&b[i__ + j * b_dim1]), abs(d__2)) + ((d__3 
			    = dl[i__4].r, abs(d__3)) + (d__4 = d_imag(&dl[i__ 
			    - 1]), abs(d__4))) * ((d__5 = x[i__5].r, abs(d__5)
			    ) + (d__6 = d_imag(&x[i__ - 1 + j * x_dim1]), abs(
			    d__6))) + ((d__7 = d__[i__6].r, abs(d__7)) + (
			    d__8 = d_imag(&d__[i__]), abs(d__8))) * ((d__9 = 
			    x[i__7].r, abs(d__9)) + (d__10 = d_imag(&x[i__ + 
			    j * x_dim1]), abs(d__10))) + ((d__11 = du[i__8].r,
			     abs(d__11)) + (d__12 = d_imag(&du[i__]), abs(
			    d__12))) * ((d__13 = x[i__9].r, abs(d__13)) + (
			    d__14 = d_imag(&x[i__ + 1 + j * x_dim1]), abs(
			    d__14)));
		}
		i__2 = *n + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n - 1 + j * x_dim1;
		i__5 = *n;
		i__6 = *n + j * x_dim1;
		rwork[*n] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			*n + j * b_dim1]), abs(d__2)) + ((d__3 = dl[i__3].r, 
			abs(d__3)) + (d__4 = d_imag(&dl[*n - 1]), abs(d__4))) 
			* ((d__5 = x[i__4].r, abs(d__5)) + (d__6 = d_imag(&x[*
			n - 1 + j * x_dim1]), abs(d__6))) + ((d__7 = d__[i__5]
			.r, abs(d__7)) + (d__8 = d_imag(&d__[*n]), abs(d__8)))
			 * ((d__9 = x[i__6].r, abs(d__9)) + (d__10 = d_imag(&
			x[*n + j * x_dim1]), abs(d__10)));
	    }
	} else {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		rwork[1] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			j * b_dim1 + 1]), abs(d__2)) + ((d__3 = d__[1].r, abs(
			d__3)) + (d__4 = d_imag(&d__[1]), abs(d__4))) * ((
			d__5 = x[i__3].r, abs(d__5)) + (d__6 = d_imag(&x[j * 
			x_dim1 + 1]), abs(d__6)));
	    } else {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		i__4 = j * x_dim1 + 2;
		rwork[1] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			j * b_dim1 + 1]), abs(d__2)) + ((d__3 = d__[1].r, abs(
			d__3)) + (d__4 = d_imag(&d__[1]), abs(d__4))) * ((
			d__5 = x[i__3].r, abs(d__5)) + (d__6 = d_imag(&x[j * 
			x_dim1 + 1]), abs(d__6))) + ((d__7 = dl[1].r, abs(
			d__7)) + (d__8 = d_imag(&dl[1]), abs(d__8))) * ((d__9 
			= x[i__4].r, abs(d__9)) + (d__10 = d_imag(&x[j * 
			x_dim1 + 2]), abs(d__10)));
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    i__3 = i__ + j * b_dim1;
		    i__4 = i__ - 1;
		    i__5 = i__ - 1 + j * x_dim1;
		    i__6 = i__;
		    i__7 = i__ + j * x_dim1;
		    i__8 = i__;
		    i__9 = i__ + 1 + j * x_dim1;
		    rwork[i__] = (d__1 = b[i__3].r, abs(d__1)) + (d__2 = 
			    d_imag(&b[i__ + j * b_dim1]), abs(d__2)) + ((d__3 
			    = du[i__4].r, abs(d__3)) + (d__4 = d_imag(&du[i__ 
			    - 1]), abs(d__4))) * ((d__5 = x[i__5].r, abs(d__5)
			    ) + (d__6 = d_imag(&x[i__ - 1 + j * x_dim1]), abs(
			    d__6))) + ((d__7 = d__[i__6].r, abs(d__7)) + (
			    d__8 = d_imag(&d__[i__]), abs(d__8))) * ((d__9 = 
			    x[i__7].r, abs(d__9)) + (d__10 = d_imag(&x[i__ + 
			    j * x_dim1]), abs(d__10))) + ((d__11 = dl[i__8].r,
			     abs(d__11)) + (d__12 = d_imag(&dl[i__]), abs(
			    d__12))) * ((d__13 = x[i__9].r, abs(d__13)) + (
			    d__14 = d_imag(&x[i__ + 1 + j * x_dim1]), abs(
			    d__14)));
		}
		i__2 = *n + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n - 1 + j * x_dim1;
		i__5 = *n;
		i__6 = *n + j * x_dim1;
		rwork[*n] = (d__1 = b[i__2].r, abs(d__1)) + (d__2 = d_imag(&b[
			*n + j * b_dim1]), abs(d__2)) + ((d__3 = du[i__3].r, 
			abs(d__3)) + (d__4 = d_imag(&du[*n - 1]), abs(d__4))) 
			* ((d__5 = x[i__4].r, abs(d__5)) + (d__6 = d_imag(&x[*
			n - 1 + j * x_dim1]), abs(d__6))) + ((d__7 = d__[i__5]
			.r, abs(d__7)) + (d__8 = d_imag(&d__[*n]), abs(d__8)))
			 * ((d__9 = x[i__6].r, abs(d__9)) + (d__10 = d_imag(&
			x[*n + j * x_dim1]), abs(d__10)));
	    }
	}
	s = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (rwork[i__] > safe2) {
/* Computing MAX */
		i__3 = i__;
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2))) / rwork[i__];
		s = max(d__3,d__4);
	    } else {
/* Computing MAX */
		i__3 = i__;
		d__3 = s, d__4 = ((d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + safe1) / (rwork[i__] 
			+ safe1);
		s = max(d__3,d__4);
	    }
	}
	berr[j] = s;
	if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {
	    zgttrsnp_(trans, n, &c__1, &dlf[1], &df[1], &duf[1], &work[1], n, 
		    info, (ftnlen)1);
	    zaxpy_(n, &c_b21, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
	    lstres = berr[j];
	    ++count;
	    goto L20;
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (rwork[i__] > safe2) {
		i__3 = i__;
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			;
	    } else {
		i__3 = i__;
		rwork[i__] = (d__1 = work[i__3].r, abs(d__1)) + (d__2 = 
			d_imag(&work[i__]), abs(d__2)) + nz * eps * rwork[i__]
			 + safe1;
	    }
	}
	kase = 0;
L70:
	zlacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
	if (kase != 0) {
	    if (kase == 1) {
		zgttrsnp_(transt, n, &c__1, &dlf[1], &df[1], &duf[1], &work[1]
			, n, info, (ftnlen)1);
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__;
		    i__4 = i__;
		    i__5 = i__;
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
		}
	    } else {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__;
		    i__4 = i__;
		    i__5 = i__;
		    z__1.r = rwork[i__4] * work[i__5].r, z__1.i = rwork[i__4] 
			    * work[i__5].i;
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
		}
		zgttrsnp_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &work[1]
			, n, info, (ftnlen)1);
	    }
	    goto L70;
	}
	lstres = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = i__ + j * x_dim1;
	    d__3 = lstres, d__4 = (d__1 = x[i__3].r, abs(d__1)) + (d__2 = 
		    d_imag(&x[i__ + j * x_dim1]), abs(d__2));
	    lstres = max(d__3,d__4);
	}
	if (lstres != 0.) {
	    ferr[j] /= lstres;
	}
    }
    return 0;
} /* zgtrfsnp_ */

