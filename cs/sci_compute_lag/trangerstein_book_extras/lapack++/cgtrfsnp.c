/* cgtrfsnp.f -- translated by f2c (version 20090411).
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
static real c_b16 = -1.f;
static real c_b17 = 1.f;
static complex c_b21 = {1.f,0.f};

/* modified from cgtrfs to avoid pivoting */
/* Subroutine */ int cgtrfsnp_(char *trans, integer *n, integer *nrhs, 
	complex *dl, complex *d__, complex *du, complex *dlf, complex *df, 
	complex *duf, complex *b, integer *ldb, complex *x, integer *ldx, 
	real *ferr, real *berr, complex *work, real *rwork, integer *info, 
	ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10, r__11, 
	    r__12, r__13, r__14;
    complex q__1;

    /* Builtin functions */
    double r_imag(complex *);

    /* Local variables */
    extern /* Subroutine */ int cgttrsnp_(char *, integer *, integer *, 
	    complex *, complex *, complex *, complex *, integer *, integer *, 
	    ftnlen);
    static integer i__, j;
    static real s;
    static integer nz;
    static real eps;
    static integer kase;
    static real safe1, safe2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer isave[3];
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    static integer count;
    extern /* Subroutine */ int clacn2_(integer *, complex *, complex *, real 
	    *, integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int clagtm_(char *, integer *, integer *, real *, 
	    complex *, complex *, complex *, complex *, integer *, real *, 
	    complex *, integer *, ftnlen);
    static real safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical notran;
    static char transn[1], transt[1];
    static real lstres;

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
	xerbla_("cgtrfsnp", &i__1, (ftnlen)8);
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
	lstres = 3.f;
L20:
	ccopy_(n, &b[j * b_dim1 + 1], &c__1, &work[1], &c__1);
	clagtm_(trans, n, &c__1, &c_b16, &dl[1], &d__[1], &du[1], &x[j * 
		x_dim1 + 1], ldx, &c_b17, &work[1], n, (ftnlen)1);
	if (notran) {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		rwork[1] = (r__1 = b[i__2].r, dabs(r__1)) + (r__2 = r_imag(&b[
			j * b_dim1 + 1]), dabs(r__2)) + ((r__3 = d__[1].r, 
			dabs(r__3)) + (r__4 = r_imag(&d__[1]), dabs(r__4))) * 
			((r__5 = x[i__3].r, dabs(r__5)) + (r__6 = r_imag(&x[j 
			* x_dim1 + 1]), dabs(r__6)));
	    } else {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		i__4 = j * x_dim1 + 2;
		rwork[1] = (r__1 = b[i__2].r, dabs(r__1)) + (r__2 = r_imag(&b[
			j * b_dim1 + 1]), dabs(r__2)) + ((r__3 = d__[1].r, 
			dabs(r__3)) + (r__4 = r_imag(&d__[1]), dabs(r__4))) * 
			((r__5 = x[i__3].r, dabs(r__5)) + (r__6 = r_imag(&x[j 
			* x_dim1 + 1]), dabs(r__6))) + ((r__7 = du[1].r, dabs(
			r__7)) + (r__8 = r_imag(&du[1]), dabs(r__8))) * ((
			r__9 = x[i__4].r, dabs(r__9)) + (r__10 = r_imag(&x[j *
			 x_dim1 + 2]), dabs(r__10)));
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    i__3 = i__ + j * b_dim1;
		    i__4 = i__ - 1;
		    i__5 = i__ - 1 + j * x_dim1;
		    i__6 = i__;
		    i__7 = i__ + j * x_dim1;
		    i__8 = i__;
		    i__9 = i__ + 1 + j * x_dim1;
		    rwork[i__] = (r__1 = b[i__3].r, dabs(r__1)) + (r__2 = 
			    r_imag(&b[i__ + j * b_dim1]), dabs(r__2)) + ((
			    r__3 = dl[i__4].r, dabs(r__3)) + (r__4 = r_imag(&
			    dl[i__ - 1]), dabs(r__4))) * ((r__5 = x[i__5].r, 
			    dabs(r__5)) + (r__6 = r_imag(&x[i__ - 1 + j * 
			    x_dim1]), dabs(r__6))) + ((r__7 = d__[i__6].r, 
			    dabs(r__7)) + (r__8 = r_imag(&d__[i__]), dabs(
			    r__8))) * ((r__9 = x[i__7].r, dabs(r__9)) + (
			    r__10 = r_imag(&x[i__ + j * x_dim1]), dabs(r__10))
			    ) + ((r__11 = du[i__8].r, dabs(r__11)) + (r__12 = 
			    r_imag(&du[i__]), dabs(r__12))) * ((r__13 = x[
			    i__9].r, dabs(r__13)) + (r__14 = r_imag(&x[i__ + 
			    1 + j * x_dim1]), dabs(r__14)));
		}
		i__2 = *n + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n - 1 + j * x_dim1;
		i__5 = *n;
		i__6 = *n + j * x_dim1;
		rwork[*n] = (r__1 = b[i__2].r, dabs(r__1)) + (r__2 = r_imag(&
			b[*n + j * b_dim1]), dabs(r__2)) + ((r__3 = dl[i__3]
			.r, dabs(r__3)) + (r__4 = r_imag(&dl[*n - 1]), dabs(
			r__4))) * ((r__5 = x[i__4].r, dabs(r__5)) + (r__6 = 
			r_imag(&x[*n - 1 + j * x_dim1]), dabs(r__6))) + ((
			r__7 = d__[i__5].r, dabs(r__7)) + (r__8 = r_imag(&d__[
			*n]), dabs(r__8))) * ((r__9 = x[i__6].r, dabs(r__9)) 
			+ (r__10 = r_imag(&x[*n + j * x_dim1]), dabs(r__10)));
	    }
	} else {
	    if (*n == 1) {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		rwork[1] = (r__1 = b[i__2].r, dabs(r__1)) + (r__2 = r_imag(&b[
			j * b_dim1 + 1]), dabs(r__2)) + ((r__3 = d__[1].r, 
			dabs(r__3)) + (r__4 = r_imag(&d__[1]), dabs(r__4))) * 
			((r__5 = x[i__3].r, dabs(r__5)) + (r__6 = r_imag(&x[j 
			* x_dim1 + 1]), dabs(r__6)));
	    } else {
		i__2 = j * b_dim1 + 1;
		i__3 = j * x_dim1 + 1;
		i__4 = j * x_dim1 + 2;
		rwork[1] = (r__1 = b[i__2].r, dabs(r__1)) + (r__2 = r_imag(&b[
			j * b_dim1 + 1]), dabs(r__2)) + ((r__3 = d__[1].r, 
			dabs(r__3)) + (r__4 = r_imag(&d__[1]), dabs(r__4))) * 
			((r__5 = x[i__3].r, dabs(r__5)) + (r__6 = r_imag(&x[j 
			* x_dim1 + 1]), dabs(r__6))) + ((r__7 = dl[1].r, dabs(
			r__7)) + (r__8 = r_imag(&dl[1]), dabs(r__8))) * ((
			r__9 = x[i__4].r, dabs(r__9)) + (r__10 = r_imag(&x[j *
			 x_dim1 + 2]), dabs(r__10)));
		i__2 = *n - 1;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    i__3 = i__ + j * b_dim1;
		    i__4 = i__ - 1;
		    i__5 = i__ - 1 + j * x_dim1;
		    i__6 = i__;
		    i__7 = i__ + j * x_dim1;
		    i__8 = i__;
		    i__9 = i__ + 1 + j * x_dim1;
		    rwork[i__] = (r__1 = b[i__3].r, dabs(r__1)) + (r__2 = 
			    r_imag(&b[i__ + j * b_dim1]), dabs(r__2)) + ((
			    r__3 = du[i__4].r, dabs(r__3)) + (r__4 = r_imag(&
			    du[i__ - 1]), dabs(r__4))) * ((r__5 = x[i__5].r, 
			    dabs(r__5)) + (r__6 = r_imag(&x[i__ - 1 + j * 
			    x_dim1]), dabs(r__6))) + ((r__7 = d__[i__6].r, 
			    dabs(r__7)) + (r__8 = r_imag(&d__[i__]), dabs(
			    r__8))) * ((r__9 = x[i__7].r, dabs(r__9)) + (
			    r__10 = r_imag(&x[i__ + j * x_dim1]), dabs(r__10))
			    ) + ((r__11 = dl[i__8].r, dabs(r__11)) + (r__12 = 
			    r_imag(&dl[i__]), dabs(r__12))) * ((r__13 = x[
			    i__9].r, dabs(r__13)) + (r__14 = r_imag(&x[i__ + 
			    1 + j * x_dim1]), dabs(r__14)));
		}
		i__2 = *n + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n - 1 + j * x_dim1;
		i__5 = *n;
		i__6 = *n + j * x_dim1;
		rwork[*n] = (r__1 = b[i__2].r, dabs(r__1)) + (r__2 = r_imag(&
			b[*n + j * b_dim1]), dabs(r__2)) + ((r__3 = du[i__3]
			.r, dabs(r__3)) + (r__4 = r_imag(&du[*n - 1]), dabs(
			r__4))) * ((r__5 = x[i__4].r, dabs(r__5)) + (r__6 = 
			r_imag(&x[*n - 1 + j * x_dim1]), dabs(r__6))) + ((
			r__7 = d__[i__5].r, dabs(r__7)) + (r__8 = r_imag(&d__[
			*n]), dabs(r__8))) * ((r__9 = x[i__6].r, dabs(r__9)) 
			+ (r__10 = r_imag(&x[*n + j * x_dim1]), dabs(r__10)));
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
	    cgttrsnp_(trans, n, &c__1, &dlf[1], &df[1], &duf[1], &work[1], n, 
		    info, (ftnlen)1);
	    caxpy_(n, &c_b21, &work[1], &c__1, &x[j * x_dim1 + 1], &c__1);
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
	kase = 0;
L70:
	clacn2_(n, &work[*n + 1], &work[1], &ferr[j], &kase, isave);
	if (kase != 0) {
	    if (kase == 1) {
		cgttrsnp_(transt, n, &c__1, &dlf[1], &df[1], &duf[1], &work[1]
			, n, info, (ftnlen)1);
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__;
		    i__4 = i__;
		    i__5 = i__;
		    q__1.r = rwork[i__4] * work[i__5].r, q__1.i = rwork[i__4] 
			    * work[i__5].i;
		    work[i__3].r = q__1.r, work[i__3].i = q__1.i;
		}
	    } else {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    i__3 = i__;
		    i__4 = i__;
		    i__5 = i__;
		    q__1.r = rwork[i__4] * work[i__5].r, q__1.i = rwork[i__4] 
			    * work[i__5].i;
		    work[i__3].r = q__1.r, work[i__3].i = q__1.i;
		}
		cgttrsnp_(transn, n, &c__1, &dlf[1], &df[1], &duf[1], &work[1]
			, n, info, (ftnlen)1);
	    }
	    goto L70;
	}
	lstres = 0.f;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    i__3 = i__ + j * x_dim1;
	    r__3 = lstres, r__4 = (r__1 = x[i__3].r, dabs(r__1)) + (r__2 = 
		    r_imag(&x[i__ + j * x_dim1]), dabs(r__2));
	    lstres = dmax(r__3,r__4);
	}
	if (lstres != 0.f) {
	    ferr[j] /= lstres;
	}
    }
    return 0;
} /* cgtrfsnp_ */

