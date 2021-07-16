/* cptrfsr.f -- translated by f2c (version 20090411).
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

/* modified from cptrfs for X A = B <==> A X^H = B^H */
/* Subroutine */ int cptrfsr_(char *uplo, integer *n, integer *nrhs, real *
	d__, complex *e, real *df, complex *ef, complex *b, integer *ldb, 
	complex *x, integer *ldx, real *ferr, real *berr, complex *work, real 
	*rwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10, r__11, 
	    r__12;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    double r_imag(complex *), c_abs(complex *);

    /* Local variables */
    static integer i__, j;
    static real s;
    static complex bi, cx, dx, ex;
    static integer ix, nz;
    static real eps, safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer count;
    static logical upper;
    extern doublereal slamch_(char *, ftnlen);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer isamax_(integer *, real *, integer *);
    static real lstres;
    extern /* Subroutine */ int cpttrs_(char *, integer *, integer *, real *, 
	    complex *, complex *, integer *, integer *, ftnlen);

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
    --rwork;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(1,*nrhs)) {
	*info = -9;
    } else if (*ldx < max(1,*nrhs)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("cptrfsr", &i__1, (ftnlen)7);
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
	if (upper) {
	    if (*n == 1) {
		r_cnjg(&q__1, &b[j + b_dim1]);
		bi.r = q__1.r, bi.i = q__1.i;
		r_cnjg(&q__2, &x[j + x_dim1]);
		q__1.r = d__[1] * q__2.r, q__1.i = d__[1] * q__2.i;
		dx.r = q__1.r, dx.i = q__1.i;
		q__1.r = bi.r - dx.r, q__1.i = bi.i - dx.i;
		work[1].r = q__1.r, work[1].i = q__1.i;
		rwork[1] = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = dx.r, dabs(r__3)) + (r__4 = 
			r_imag(&dx), dabs(r__4)));
	    } else {
		r_cnjg(&q__1, &b[j + b_dim1]);
		bi.r = q__1.r, bi.i = q__1.i;
		r_cnjg(&q__2, &x[j + x_dim1]);
		q__1.r = d__[1] * q__2.r, q__1.i = d__[1] * q__2.i;
		dx.r = q__1.r, dx.i = q__1.i;
		r_cnjg(&q__2, &x[j + (x_dim1 << 1)]);
		q__1.r = e[1].r * q__2.r - e[1].i * q__2.i, q__1.i = e[1].r * 
			q__2.i + e[1].i * q__2.r;
		ex.r = q__1.r, ex.i = q__1.i;
		q__2.r = bi.r - dx.r, q__2.i = bi.i - dx.i;
		q__1.r = q__2.r - ex.r, q__1.i = q__2.i - ex.i;
		work[1].r = q__1.r, work[1].i = q__1.i;
		i__2 = j + (x_dim1 << 1);
		rwork[1] = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = dx.r, dabs(r__3)) + (r__4 = 
			r_imag(&dx), dabs(r__4))) + ((r__5 = e[1].r, dabs(
			r__5)) + (r__6 = r_imag(&e[1]), dabs(r__6))) * ((r__7 
			= x[i__2].r, dabs(r__7)) + (r__8 = r_imag(&x[j + (
			x_dim1 << 1)]), dabs(r__8)));
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    r_cnjg(&q__1, &b[j + i__ * b_dim1]);
		    bi.r = q__1.r, bi.i = q__1.i;
		    r_cnjg(&q__2, &e[i__ - 1]);
		    r_cnjg(&q__3, &x[j + (i__ - 1) * x_dim1]);
		    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = 
			    q__2.r * q__3.i + q__2.i * q__3.r;
		    cx.r = q__1.r, cx.i = q__1.i;
		    i__3 = i__;
		    r_cnjg(&q__2, &x[j + i__ * x_dim1]);
		    q__1.r = d__[i__3] * q__2.r, q__1.i = d__[i__3] * q__2.i;
		    dx.r = q__1.r, dx.i = q__1.i;
		    i__3 = i__;
		    r_cnjg(&q__2, &x[j + (i__ + 1) * x_dim1]);
		    q__1.r = e[i__3].r * q__2.r - e[i__3].i * q__2.i, q__1.i =
			     e[i__3].r * q__2.i + e[i__3].i * q__2.r;
		    ex.r = q__1.r, ex.i = q__1.i;
		    i__3 = i__;
		    q__3.r = bi.r - cx.r, q__3.i = bi.i - cx.i;
		    q__2.r = q__3.r - dx.r, q__2.i = q__3.i - dx.i;
		    q__1.r = q__2.r - ex.r, q__1.i = q__2.i - ex.i;
		    work[i__3].r = q__1.r, work[i__3].i = q__1.i;
		    i__3 = i__ - 1;
		    i__4 = j + (i__ - 1) * x_dim1;
		    i__5 = i__;
		    i__6 = i__ + 1 + j * x_dim1;
		    rwork[i__] = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&
			    bi), dabs(r__2)) + ((r__3 = e[i__3].r, dabs(r__3))
			     + (r__4 = r_imag(&e[i__ - 1]), dabs(r__4))) * ((
			    r__5 = x[i__4].r, dabs(r__5)) + (r__6 = r_imag(&x[
			    j + (i__ - 1) * x_dim1]), dabs(r__6))) + ((r__7 = 
			    dx.r, dabs(r__7)) + (r__8 = r_imag(&dx), dabs(
			    r__8))) + ((r__9 = e[i__5].r, dabs(r__9)) + (
			    r__10 = r_imag(&e[i__]), dabs(r__10))) * ((r__11 =
			     x[i__6].r, dabs(r__11)) + (r__12 = r_imag(&x[i__ 
			    + 1 + j * x_dim1]), dabs(r__12)));
		}
		r_cnjg(&q__1, &b[j + *n * b_dim1]);
		bi.r = q__1.r, bi.i = q__1.i;
		r_cnjg(&q__2, &e[*n - 1]);
		r_cnjg(&q__3, &x[j + (*n - 1) * x_dim1]);
		q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * 
			q__3.i + q__2.i * q__3.r;
		cx.r = q__1.r, cx.i = q__1.i;
		i__2 = *n;
		r_cnjg(&q__2, &x[j + *n * x_dim1]);
		q__1.r = d__[i__2] * q__2.r, q__1.i = d__[i__2] * q__2.i;
		dx.r = q__1.r, dx.i = q__1.i;
		i__2 = *n;
		q__2.r = bi.r - cx.r, q__2.i = bi.i - cx.i;
		q__1.r = q__2.r - dx.r, q__1.i = q__2.i - dx.i;
		work[i__2].r = q__1.r, work[i__2].i = q__1.i;
		i__2 = *n - 1;
		i__3 = j + (*n - 1) * x_dim1;
		rwork[*n] = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = e[i__2].r, dabs(r__3)) + (r__4 
			= r_imag(&e[*n - 1]), dabs(r__4))) * ((r__5 = x[i__3]
			.r, dabs(r__5)) + (r__6 = r_imag(&x[j + (*n - 1) * 
			x_dim1]), dabs(r__6))) + ((r__7 = dx.r, dabs(r__7)) + 
			(r__8 = r_imag(&dx), dabs(r__8)));
	    }
	} else {
	    if (*n == 1) {
		r_cnjg(&q__1, &b[j + b_dim1]);
		bi.r = q__1.r, bi.i = q__1.i;
		r_cnjg(&q__2, &x[j + x_dim1]);
		q__1.r = d__[1] * q__2.r, q__1.i = d__[1] * q__2.i;
		dx.r = q__1.r, dx.i = q__1.i;
		q__1.r = bi.r - dx.r, q__1.i = bi.i - dx.i;
		work[1].r = q__1.r, work[1].i = q__1.i;
		rwork[1] = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = dx.r, dabs(r__3)) + (r__4 = 
			r_imag(&dx), dabs(r__4)));
	    } else {
		r_cnjg(&q__1, &b[j + b_dim1]);
		bi.r = q__1.r, bi.i = q__1.i;
		r_cnjg(&q__2, &x[j + x_dim1]);
		q__1.r = d__[1] * q__2.r, q__1.i = d__[1] * q__2.i;
		dx.r = q__1.r, dx.i = q__1.i;
		r_cnjg(&q__2, &e[1]);
		r_cnjg(&q__3, &x[j + (x_dim1 << 1)]);
		q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * 
			q__3.i + q__2.i * q__3.r;
		ex.r = q__1.r, ex.i = q__1.i;
		q__2.r = bi.r - dx.r, q__2.i = bi.i - dx.i;
		q__1.r = q__2.r - ex.r, q__1.i = q__2.i - ex.i;
		work[1].r = q__1.r, work[1].i = q__1.i;
		i__2 = j + (x_dim1 << 1);
		rwork[1] = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = dx.r, dabs(r__3)) + (r__4 = 
			r_imag(&dx), dabs(r__4))) + ((r__5 = e[1].r, dabs(
			r__5)) + (r__6 = r_imag(&e[1]), dabs(r__6))) * ((r__7 
			= x[i__2].r, dabs(r__7)) + (r__8 = r_imag(&x[j + (
			x_dim1 << 1)]), dabs(r__8)));
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    r_cnjg(&q__1, &b[j + i__ * b_dim1]);
		    bi.r = q__1.r, bi.i = q__1.i;
		    i__3 = i__ - 1;
		    r_cnjg(&q__2, &x[j + (i__ - 1) * x_dim1]);
		    q__1.r = e[i__3].r * q__2.r - e[i__3].i * q__2.i, q__1.i =
			     e[i__3].r * q__2.i + e[i__3].i * q__2.r;
		    cx.r = q__1.r, cx.i = q__1.i;
		    i__3 = i__;
		    r_cnjg(&q__2, &x[j + i__ * x_dim1]);
		    q__1.r = d__[i__3] * q__2.r, q__1.i = d__[i__3] * q__2.i;
		    dx.r = q__1.r, dx.i = q__1.i;
		    r_cnjg(&q__2, &e[i__]);
		    r_cnjg(&q__3, &x[j + (i__ + 1) * x_dim1]);
		    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = 
			    q__2.r * q__3.i + q__2.i * q__3.r;
		    ex.r = q__1.r, ex.i = q__1.i;
		    i__3 = i__;
		    q__3.r = bi.r - cx.r, q__3.i = bi.i - cx.i;
		    q__2.r = q__3.r - dx.r, q__2.i = q__3.i - dx.i;
		    q__1.r = q__2.r - ex.r, q__1.i = q__2.i - ex.i;
		    work[i__3].r = q__1.r, work[i__3].i = q__1.i;
		    i__3 = i__ - 1;
		    i__4 = j + (i__ - 1) * x_dim1;
		    i__5 = i__;
		    i__6 = j + (i__ + 1) * x_dim1;
		    rwork[i__] = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&
			    bi), dabs(r__2)) + ((r__3 = e[i__3].r, dabs(r__3))
			     + (r__4 = r_imag(&e[i__ - 1]), dabs(r__4))) * ((
			    r__5 = x[i__4].r, dabs(r__5)) + (r__6 = r_imag(&x[
			    j + (i__ - 1) * x_dim1]), dabs(r__6))) + ((r__7 = 
			    dx.r, dabs(r__7)) + (r__8 = r_imag(&dx), dabs(
			    r__8))) + ((r__9 = e[i__5].r, dabs(r__9)) + (
			    r__10 = r_imag(&e[i__]), dabs(r__10))) * ((r__11 =
			     x[i__6].r, dabs(r__11)) + (r__12 = r_imag(&x[j + 
			    (i__ + 1) * x_dim1]), dabs(r__12)));
		}
		r_cnjg(&q__1, &b[j + *n * b_dim1]);
		bi.r = q__1.r, bi.i = q__1.i;
		i__2 = *n - 1;
		r_cnjg(&q__2, &x[j + (*n - 1) * x_dim1]);
		q__1.r = e[i__2].r * q__2.r - e[i__2].i * q__2.i, q__1.i = e[
			i__2].r * q__2.i + e[i__2].i * q__2.r;
		cx.r = q__1.r, cx.i = q__1.i;
		i__2 = *n;
		r_cnjg(&q__2, &x[j + *n * x_dim1]);
		q__1.r = d__[i__2] * q__2.r, q__1.i = d__[i__2] * q__2.i;
		dx.r = q__1.r, dx.i = q__1.i;
		i__2 = *n;
		q__2.r = bi.r - cx.r, q__2.i = bi.i - cx.i;
		q__1.r = q__2.r - dx.r, q__1.i = q__2.i - dx.i;
		work[i__2].r = q__1.r, work[i__2].i = q__1.i;
		i__2 = *n - 1;
		i__3 = j + (*n - 1) * x_dim1;
		rwork[*n] = (r__1 = bi.r, dabs(r__1)) + (r__2 = r_imag(&bi), 
			dabs(r__2)) + ((r__3 = e[i__2].r, dabs(r__3)) + (r__4 
			= r_imag(&e[*n - 1]), dabs(r__4))) * ((r__5 = x[i__3]
			.r, dabs(r__5)) + (r__6 = r_imag(&x[j + (*n - 1) * 
			x_dim1]), dabs(r__6))) + ((r__7 = dx.r, dabs(r__7)) + 
			(r__8 = r_imag(&dx), dabs(r__8)));
	    }
	}
	s = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (rwork[i__] > safe2) {
/* Computing MAX */
		i__3 = i__;
		r__3 = s, r__4 = ((r__1 = work[i__3].r, dabs(r__1)) + (r__2 = 
			r_imag(&work[i__]), dabs(r__2))) / rwork[i__];
		s = dmax(r__3,r__4);
	    } else {
/* Computing MAX */
		i__3 = i__;
		r__3 = s, r__4 = ((r__1 = work[i__3].r, dabs(r__1)) + (r__2 = 
			r_imag(&work[i__]), dabs(r__2)) + safe1) / (rwork[i__]
			 + safe1);
		s = dmax(r__3,r__4);
	    }
	}
	berr[j] = s;
	if (berr[j] > eps && berr[j] * 2.f <= lstres && count <= 5) {
	    cpttrs_(uplo, n, &c__1, &df[1], &ef[1], &work[1], n, info, (
		    ftnlen)1);
/*         call caxpy(n,cmplx(one),work,1,x(1,j),1) */
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = j + i__ * x_dim1;
		i__4 = j + i__ * x_dim1;
		r_cnjg(&q__2, &work[i__]);
		q__1.r = x[i__4].r + q__2.r, q__1.i = x[i__4].i + q__2.i;
		x[i__3].r = q__1.r, x[i__3].i = q__1.i;
	    }
	    lstres = berr[j];
	    ++count;
	    goto L20;
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (rwork[i__] > safe2) {
		i__3 = i__;
		rwork[i__] = (r__1 = work[i__3].r, dabs(r__1)) + (r__2 = 
			r_imag(&work[i__]), dabs(r__2)) + nz * eps * rwork[
			i__];
	    } else {
		i__3 = i__;
		rwork[i__] = (r__1 = work[i__3].r, dabs(r__1)) + (r__2 = 
			r_imag(&work[i__]), dabs(r__2)) + nz * eps * rwork[
			i__] + safe1;
	    }
	}
	ix = isamax_(n, &rwork[1], &c__1);
	ferr[j] = rwork[ix];
	rwork[1] = 1.f;
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    rwork[i__] = rwork[i__ - 1] * c_abs(&ef[i__ - 1]) + 1.f;
	}
	rwork[*n] /= df[*n];
	for (i__ = *n - 1; i__ >= 1; --i__) {
	    rwork[i__] = rwork[i__] / df[i__] + rwork[i__ + 1] * c_abs(&ef[
		    i__]);
	}
	ix = isamax_(n, &rwork[1], &c__1);
	ferr[j] *= (r__1 = rwork[ix], dabs(r__1));
	lstres = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__1 = lstres, r__2 = c_abs(&x[j + i__ * x_dim1]);
	    lstres = dmax(r__1,r__2);
	}
	if (lstres != 0.f) {
	    ferr[j] /= lstres;
	}
    }
    return 0;
} /* cptrfsr_ */

