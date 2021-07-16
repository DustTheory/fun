/* zptrfsr.f -- translated by f2c (version 20090411).
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

/* modified from zptrfs for X A = B <==> A X^H = B^H */
/* Subroutine */ int zptrfsr_(char *uplo, integer *n, integer *nrhs, 
	doublereal *d__, doublecomplex *e, doublereal *df, doublecomplex *ef, 
	doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, 
	doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *
	rwork, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *), z_abs(doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublereal s;
    static doublecomplex bi, cx, dx, ex;
    static integer ix, nz;
    static doublereal eps, safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer count;
    static logical upper;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal lstres;
    extern /* Subroutine */ int zpttrs_(char *, integer *, integer *, 
	    doublereal *, doublecomplex *, doublecomplex *, integer *, 
	    integer *, ftnlen);

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
	xerbla_("zptrfsr", &i__1, (ftnlen)7);
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
	if (upper) {
	    if (*n == 1) {
		d_cnjg(&z__1, &b[j + b_dim1]);
		bi.r = z__1.r, bi.i = z__1.i;
		d_cnjg(&z__2, &x[j + x_dim1]);
		z__1.r = d__[1] * z__2.r, z__1.i = d__[1] * z__2.i;
		dx.r = z__1.r, dx.i = z__1.i;
		z__1.r = bi.r - dx.r, z__1.i = bi.i - dx.i;
		work[1].r = z__1.r, work[1].i = z__1.i;
		rwork[1] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = dx.r, abs(d__3)) + (d__4 = 
			d_imag(&dx), abs(d__4)));
	    } else {
		d_cnjg(&z__1, &b[j + b_dim1]);
		bi.r = z__1.r, bi.i = z__1.i;
		d_cnjg(&z__2, &x[j + x_dim1]);
		z__1.r = d__[1] * z__2.r, z__1.i = d__[1] * z__2.i;
		dx.r = z__1.r, dx.i = z__1.i;
		d_cnjg(&z__2, &x[j + (x_dim1 << 1)]);
		z__1.r = e[1].r * z__2.r - e[1].i * z__2.i, z__1.i = e[1].r * 
			z__2.i + e[1].i * z__2.r;
		ex.r = z__1.r, ex.i = z__1.i;
		z__2.r = bi.r - dx.r, z__2.i = bi.i - dx.i;
		z__1.r = z__2.r - ex.r, z__1.i = z__2.i - ex.i;
		work[1].r = z__1.r, work[1].i = z__1.i;
		i__2 = j + (x_dim1 << 1);
		rwork[1] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = dx.r, abs(d__3)) + (d__4 = 
			d_imag(&dx), abs(d__4))) + ((d__5 = e[1].r, abs(d__5))
			 + (d__6 = d_imag(&e[1]), abs(d__6))) * ((d__7 = x[
			i__2].r, abs(d__7)) + (d__8 = d_imag(&x[j + (x_dim1 <<
			 1)]), abs(d__8)));
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    d_cnjg(&z__1, &b[j + i__ * b_dim1]);
		    bi.r = z__1.r, bi.i = z__1.i;
		    d_cnjg(&z__2, &e[i__ - 1]);
		    d_cnjg(&z__3, &x[j + (i__ - 1) * x_dim1]);
		    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
			    z__2.r * z__3.i + z__2.i * z__3.r;
		    cx.r = z__1.r, cx.i = z__1.i;
		    i__3 = i__;
		    d_cnjg(&z__2, &x[j + i__ * x_dim1]);
		    z__1.r = d__[i__3] * z__2.r, z__1.i = d__[i__3] * z__2.i;
		    dx.r = z__1.r, dx.i = z__1.i;
		    i__3 = i__;
		    d_cnjg(&z__2, &x[j + (i__ + 1) * x_dim1]);
		    z__1.r = e[i__3].r * z__2.r - e[i__3].i * z__2.i, z__1.i =
			     e[i__3].r * z__2.i + e[i__3].i * z__2.r;
		    ex.r = z__1.r, ex.i = z__1.i;
		    i__3 = i__;
		    z__3.r = bi.r - cx.r, z__3.i = bi.i - cx.i;
		    z__2.r = z__3.r - dx.r, z__2.i = z__3.i - dx.i;
		    z__1.r = z__2.r - ex.r, z__1.i = z__2.i - ex.i;
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
		    i__3 = i__ - 1;
		    i__4 = j + (i__ - 1) * x_dim1;
		    i__5 = i__;
		    i__6 = i__ + 1 + j * x_dim1;
		    rwork[i__] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&
			    bi), abs(d__2)) + ((d__3 = e[i__3].r, abs(d__3)) 
			    + (d__4 = d_imag(&e[i__ - 1]), abs(d__4))) * ((
			    d__5 = x[i__4].r, abs(d__5)) + (d__6 = d_imag(&x[
			    j + (i__ - 1) * x_dim1]), abs(d__6))) + ((d__7 = 
			    dx.r, abs(d__7)) + (d__8 = d_imag(&dx), abs(d__8))
			    ) + ((d__9 = e[i__5].r, abs(d__9)) + (d__10 = 
			    d_imag(&e[i__]), abs(d__10))) * ((d__11 = x[i__6]
			    .r, abs(d__11)) + (d__12 = d_imag(&x[i__ + 1 + j *
			     x_dim1]), abs(d__12)));
		}
		d_cnjg(&z__1, &b[j + *n * b_dim1]);
		bi.r = z__1.r, bi.i = z__1.i;
		d_cnjg(&z__2, &e[*n - 1]);
		d_cnjg(&z__3, &x[j + (*n - 1) * x_dim1]);
		z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
			z__3.i + z__2.i * z__3.r;
		cx.r = z__1.r, cx.i = z__1.i;
		i__2 = *n;
		d_cnjg(&z__2, &x[j + *n * x_dim1]);
		z__1.r = d__[i__2] * z__2.r, z__1.i = d__[i__2] * z__2.i;
		dx.r = z__1.r, dx.i = z__1.i;
		i__2 = *n;
		z__2.r = bi.r - cx.r, z__2.i = bi.i - cx.i;
		z__1.r = z__2.r - dx.r, z__1.i = z__2.i - dx.i;
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
		i__2 = *n - 1;
		i__3 = j + (*n - 1) * x_dim1;
		rwork[*n] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = e[i__2].r, abs(d__3)) + (d__4 = 
			d_imag(&e[*n - 1]), abs(d__4))) * ((d__5 = x[i__3].r, 
			abs(d__5)) + (d__6 = d_imag(&x[j + (*n - 1) * x_dim1])
			, abs(d__6))) + ((d__7 = dx.r, abs(d__7)) + (d__8 = 
			d_imag(&dx), abs(d__8)));
	    }
	} else {
	    if (*n == 1) {
		d_cnjg(&z__1, &b[j + b_dim1]);
		bi.r = z__1.r, bi.i = z__1.i;
		d_cnjg(&z__2, &x[j + x_dim1]);
		z__1.r = d__[1] * z__2.r, z__1.i = d__[1] * z__2.i;
		dx.r = z__1.r, dx.i = z__1.i;
		z__1.r = bi.r - dx.r, z__1.i = bi.i - dx.i;
		work[1].r = z__1.r, work[1].i = z__1.i;
		rwork[1] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = dx.r, abs(d__3)) + (d__4 = 
			d_imag(&dx), abs(d__4)));
	    } else {
		d_cnjg(&z__1, &b[j + b_dim1]);
		bi.r = z__1.r, bi.i = z__1.i;
		d_cnjg(&z__2, &x[j + x_dim1]);
		z__1.r = d__[1] * z__2.r, z__1.i = d__[1] * z__2.i;
		dx.r = z__1.r, dx.i = z__1.i;
		d_cnjg(&z__2, &e[1]);
		d_cnjg(&z__3, &x[j + (x_dim1 << 1)]);
		z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = z__2.r * 
			z__3.i + z__2.i * z__3.r;
		ex.r = z__1.r, ex.i = z__1.i;
		z__2.r = bi.r - dx.r, z__2.i = bi.i - dx.i;
		z__1.r = z__2.r - ex.r, z__1.i = z__2.i - ex.i;
		work[1].r = z__1.r, work[1].i = z__1.i;
		i__2 = j + (x_dim1 << 1);
		rwork[1] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = dx.r, abs(d__3)) + (d__4 = 
			d_imag(&dx), abs(d__4))) + ((d__5 = e[1].r, abs(d__5))
			 + (d__6 = d_imag(&e[1]), abs(d__6))) * ((d__7 = x[
			i__2].r, abs(d__7)) + (d__8 = d_imag(&x[j + (x_dim1 <<
			 1)]), abs(d__8)));
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    d_cnjg(&z__1, &b[j + i__ * b_dim1]);
		    bi.r = z__1.r, bi.i = z__1.i;
		    i__3 = i__ - 1;
		    d_cnjg(&z__2, &x[j + (i__ - 1) * x_dim1]);
		    z__1.r = e[i__3].r * z__2.r - e[i__3].i * z__2.i, z__1.i =
			     e[i__3].r * z__2.i + e[i__3].i * z__2.r;
		    cx.r = z__1.r, cx.i = z__1.i;
		    i__3 = i__;
		    d_cnjg(&z__2, &x[j + i__ * x_dim1]);
		    z__1.r = d__[i__3] * z__2.r, z__1.i = d__[i__3] * z__2.i;
		    dx.r = z__1.r, dx.i = z__1.i;
		    d_cnjg(&z__2, &e[i__]);
		    d_cnjg(&z__3, &x[j + (i__ + 1) * x_dim1]);
		    z__1.r = z__2.r * z__3.r - z__2.i * z__3.i, z__1.i = 
			    z__2.r * z__3.i + z__2.i * z__3.r;
		    ex.r = z__1.r, ex.i = z__1.i;
		    i__3 = i__;
		    z__3.r = bi.r - cx.r, z__3.i = bi.i - cx.i;
		    z__2.r = z__3.r - dx.r, z__2.i = z__3.i - dx.i;
		    z__1.r = z__2.r - ex.r, z__1.i = z__2.i - ex.i;
		    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
		    i__3 = i__ - 1;
		    i__4 = j + (i__ - 1) * x_dim1;
		    i__5 = i__;
		    i__6 = j + (i__ + 1) * x_dim1;
		    rwork[i__] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&
			    bi), abs(d__2)) + ((d__3 = e[i__3].r, abs(d__3)) 
			    + (d__4 = d_imag(&e[i__ - 1]), abs(d__4))) * ((
			    d__5 = x[i__4].r, abs(d__5)) + (d__6 = d_imag(&x[
			    j + (i__ - 1) * x_dim1]), abs(d__6))) + ((d__7 = 
			    dx.r, abs(d__7)) + (d__8 = d_imag(&dx), abs(d__8))
			    ) + ((d__9 = e[i__5].r, abs(d__9)) + (d__10 = 
			    d_imag(&e[i__]), abs(d__10))) * ((d__11 = x[i__6]
			    .r, abs(d__11)) + (d__12 = d_imag(&x[j + (i__ + 1)
			     * x_dim1]), abs(d__12)));
		}
		d_cnjg(&z__1, &b[j + *n * b_dim1]);
		bi.r = z__1.r, bi.i = z__1.i;
		i__2 = *n - 1;
		d_cnjg(&z__2, &x[j + (*n - 1) * x_dim1]);
		z__1.r = e[i__2].r * z__2.r - e[i__2].i * z__2.i, z__1.i = e[
			i__2].r * z__2.i + e[i__2].i * z__2.r;
		cx.r = z__1.r, cx.i = z__1.i;
		i__2 = *n;
		d_cnjg(&z__2, &x[j + *n * x_dim1]);
		z__1.r = d__[i__2] * z__2.r, z__1.i = d__[i__2] * z__2.i;
		dx.r = z__1.r, dx.i = z__1.i;
		i__2 = *n;
		z__2.r = bi.r - cx.r, z__2.i = bi.i - cx.i;
		z__1.r = z__2.r - dx.r, z__1.i = z__2.i - dx.i;
		work[i__2].r = z__1.r, work[i__2].i = z__1.i;
		i__2 = *n - 1;
		i__3 = j + (*n - 1) * x_dim1;
		rwork[*n] = (d__1 = bi.r, abs(d__1)) + (d__2 = d_imag(&bi), 
			abs(d__2)) + ((d__3 = e[i__2].r, abs(d__3)) + (d__4 = 
			d_imag(&e[*n - 1]), abs(d__4))) * ((d__5 = x[i__3].r, 
			abs(d__5)) + (d__6 = d_imag(&x[j + (*n - 1) * x_dim1])
			, abs(d__6))) + ((d__7 = dx.r, abs(d__7)) + (d__8 = 
			d_imag(&dx), abs(d__8)));
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
	    zpttrs_(uplo, n, &c__1, &df[1], &ef[1], &work[1], n, info, (
		    ftnlen)1);
/*         call zaxpy(n,dcmplx(one),work,1,x(1,j),1) */
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = j + i__ * x_dim1;
		i__4 = j + i__ * x_dim1;
		d_cnjg(&z__2, &work[i__]);
		z__1.r = x[i__4].r + z__2.r, z__1.i = x[i__4].i + z__2.i;
		x[i__3].r = z__1.r, x[i__3].i = z__1.i;
	    }
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
	ix = idamax_(n, &rwork[1], &c__1);
	ferr[j] = rwork[ix];
	rwork[1] = 1.;
	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    rwork[i__] = rwork[i__ - 1] * z_abs(&ef[i__ - 1]) + 1.;
	}
	rwork[*n] /= df[*n];
	for (i__ = *n - 1; i__ >= 1; --i__) {
	    rwork[i__] = rwork[i__] / df[i__] + rwork[i__ + 1] * z_abs(&ef[
		    i__]);
	}
	ix = idamax_(n, &rwork[1], &c__1);
	ferr[j] *= (d__1 = rwork[ix], abs(d__1));
	lstres = 0.;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    d__1 = lstres, d__2 = z_abs(&x[j + i__ * x_dim1]);
	    lstres = max(d__1,d__2);
	}
	if (lstres != 0.) {
	    ferr[j] /= lstres;
	}
    }
    return 0;
} /* zptrfsr_ */

