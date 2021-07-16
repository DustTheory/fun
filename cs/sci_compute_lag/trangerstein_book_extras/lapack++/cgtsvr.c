/* cgtsvr.f -- translated by f2c (version 20090411).
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

/* modified from cgtsv to solve for rows */
/* Subroutine */ int cgtsvr_(integer *n, integer *nrhs, complex *dl, complex *
	d__, complex *du, complex *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    static integer j, k;
    static complex temp, mult;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*ldb < max(1,*nrhs)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("cgtsvr", &i__1, (ftnlen)6);
	return 0;
    }
    if (*n == 0) {
	return 0;
    }
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if (dl[i__2].r == 0.f && dl[i__2].i == 0.f) {
	    i__2 = k;
	    if (d__[i__2].r == 0.f && d__[i__2].i == 0.f) {
		*info = k;
		return 0;
	    }
	} else /* if(complicated condition) */ {
	    i__2 = k;
	    i__3 = k;
	    if ((r__1 = d__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&d__[k]), 
		    dabs(r__2)) >= (r__3 = dl[i__3].r, dabs(r__3)) + (r__4 = 
		    r_imag(&dl[k]), dabs(r__4))) {
		c_div(&q__1, &dl[k], &d__[k]);
		mult.r = q__1.r, mult.i = q__1.i;
		i__2 = k + 1;
		i__3 = k + 1;
		i__4 = k;
		q__2.r = mult.r * du[i__4].r - mult.i * du[i__4].i, q__2.i = 
			mult.r * du[i__4].i + mult.i * du[i__4].r;
		q__1.r = d__[i__3].r - q__2.r, q__1.i = d__[i__3].i - q__2.i;
		d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
		i__2 = *nrhs;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = j + (k + 1) * b_dim1;
		    i__4 = j + (k + 1) * b_dim1;
		    i__5 = j + k * b_dim1;
		    q__2.r = mult.r * b[i__5].r - mult.i * b[i__5].i, q__2.i =
			     mult.r * b[i__5].i + mult.i * b[i__5].r;
		    q__1.r = b[i__4].r - q__2.r, q__1.i = b[i__4].i - q__2.i;
		    b[i__3].r = q__1.r, b[i__3].i = q__1.i;
		}
		if (k < *n - 1) {
		    i__2 = k;
		    dl[i__2].r = 0.f, dl[i__2].i = 0.f;
		}
	    } else {
		c_div(&q__1, &d__[k], &dl[k]);
		mult.r = q__1.r, mult.i = q__1.i;
		i__2 = k;
		i__3 = k;
		d__[i__2].r = dl[i__3].r, d__[i__2].i = dl[i__3].i;
		i__2 = k + 1;
		temp.r = d__[i__2].r, temp.i = d__[i__2].i;
		i__2 = k + 1;
		i__3 = k;
		q__2.r = mult.r * temp.r - mult.i * temp.i, q__2.i = mult.r * 
			temp.i + mult.i * temp.r;
		q__1.r = du[i__3].r - q__2.r, q__1.i = du[i__3].i - q__2.i;
		d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
		if (k < *n - 1) {
		    i__2 = k;
		    i__3 = k + 1;
		    dl[i__2].r = du[i__3].r, dl[i__2].i = du[i__3].i;
		    i__2 = k + 1;
		    q__2.r = -mult.r, q__2.i = -mult.i;
		    i__3 = k;
		    q__1.r = q__2.r * dl[i__3].r - q__2.i * dl[i__3].i, 
			    q__1.i = q__2.r * dl[i__3].i + q__2.i * dl[i__3]
			    .r;
		    du[i__2].r = q__1.r, du[i__2].i = q__1.i;
		}
		i__2 = k;
		du[i__2].r = temp.r, du[i__2].i = temp.i;
		i__2 = *nrhs;
		for (j = 1; j <= i__2; ++j) {
		    i__3 = j + k * b_dim1;
		    temp.r = b[i__3].r, temp.i = b[i__3].i;
		    i__3 = j + k * b_dim1;
		    i__4 = j + (k + 1) * b_dim1;
		    b[i__3].r = b[i__4].r, b[i__3].i = b[i__4].i;
		    i__3 = j + (k + 1) * b_dim1;
		    i__4 = j + (k + 1) * b_dim1;
		    q__2.r = mult.r * b[i__4].r - mult.i * b[i__4].i, q__2.i =
			     mult.r * b[i__4].i + mult.i * b[i__4].r;
		    q__1.r = temp.r - q__2.r, q__1.i = temp.i - q__2.i;
		    b[i__3].r = q__1.r, b[i__3].i = q__1.i;
		}
	    }
	}
    }
    i__1 = *n;
    if (d__[i__1].r == 0.f && d__[i__1].i == 0.f) {
	*info = *n;
	return 0;
    }
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + *n * b_dim1;
	c_div(&q__1, &b[j + *n * b_dim1], &d__[*n]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	if (*n > 1) {
	    i__2 = j + (*n - 1) * b_dim1;
	    i__3 = j + (*n - 1) * b_dim1;
	    i__4 = *n - 1;
	    i__5 = j + *n * b_dim1;
	    q__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, q__3.i =
		     du[i__4].r * b[i__5].i + du[i__4].i * b[i__5].r;
	    q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
	    c_div(&q__1, &q__2, &d__[*n - 1]);
	    b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	}
	for (k = *n - 2; k >= 1; --k) {
	    i__2 = j + k * b_dim1;
	    i__3 = j + k * b_dim1;
	    i__4 = k;
	    i__5 = j + (k + 1) * b_dim1;
	    q__4.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, q__4.i =
		     du[i__4].r * b[i__5].i + du[i__4].i * b[i__5].r;
	    q__3.r = b[i__3].r - q__4.r, q__3.i = b[i__3].i - q__4.i;
	    i__6 = k;
	    i__7 = j + (k + 2) * b_dim1;
	    q__5.r = dl[i__6].r * b[i__7].r - dl[i__6].i * b[i__7].i, q__5.i =
		     dl[i__6].r * b[i__7].i + dl[i__6].i * b[i__7].r;
	    q__2.r = q__3.r - q__5.r, q__2.i = q__3.i - q__5.i;
	    c_div(&q__1, &q__2, &d__[k]);
	    b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	}
    }
    return 0;
} /* cgtsvr_ */

