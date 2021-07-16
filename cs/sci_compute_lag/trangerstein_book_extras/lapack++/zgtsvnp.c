/* zgtsvnp.f -- translated by f2c (version 20090411).
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

/* modified from zgtsv to avoid pivoting */
/* Subroutine */ int zgtsvnp_(integer *n, integer *nrhs, doublecomplex *dl, 
	doublecomplex *d__, doublecomplex *du, doublecomplex *b, integer *ldb,
	 integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex mult;
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
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("zgtsvnp", &i__1, (ftnlen)7);
	return 0;
    }
    if (*n == 0) {
	return 0;
    }
    i__1 = *n - 1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if (dl[i__2].r == 0. && dl[i__2].i == 0.) {
	    i__2 = k;
	    if (d__[i__2].r == 0. && d__[i__2].i == 0.) {
		*info = k;
		return 0;
	    }
	}
	z_div(&z__1, &dl[k], &d__[k]);
	mult.r = z__1.r, mult.i = z__1.i;
	i__2 = k + 1;
	i__3 = k + 1;
	i__4 = k;
	z__2.r = mult.r * du[i__4].r - mult.i * du[i__4].i, z__2.i = mult.r * 
		du[i__4].i + mult.i * du[i__4].r;
	z__1.r = d__[i__3].r - z__2.r, z__1.i = d__[i__3].i - z__2.i;
	d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
	i__2 = *nrhs;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = k + 1 + j * b_dim1;
	    i__4 = k + 1 + j * b_dim1;
	    i__5 = k + j * b_dim1;
	    z__2.r = mult.r * b[i__5].r - mult.i * b[i__5].i, z__2.i = mult.r 
		    * b[i__5].i + mult.i * b[i__5].r;
	    z__1.r = b[i__4].r - z__2.r, z__1.i = b[i__4].i - z__2.i;
	    b[i__3].r = z__1.r, b[i__3].i = z__1.i;
	}
	if (k < *n - 1) {
	    i__2 = k;
	    dl[i__2].r = 0., dl[i__2].i = 0.;
	}
    }
    i__1 = *n;
    if (d__[i__1].r == 0. && d__[i__1].i == 0.) {
	*info = *n;
	return 0;
    }
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n + j * b_dim1;
	z_div(&z__1, &b[*n + j * b_dim1], &d__[*n]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	if (*n > 1) {
	    i__2 = *n - 1 + j * b_dim1;
	    i__3 = *n - 1 + j * b_dim1;
	    i__4 = *n - 1;
	    i__5 = *n + j * b_dim1;
	    z__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, z__3.i =
		     du[i__4].r * b[i__5].i + du[i__4].i * b[i__5].r;
	    z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
	    z_div(&z__1, &z__2, &d__[*n - 1]);
	    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	}
	for (k = *n - 2; k >= 1; --k) {
	    i__2 = k + j * b_dim1;
	    i__3 = k + j * b_dim1;
	    i__4 = k;
	    i__5 = k + 1 + j * b_dim1;
	    z__4.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, z__4.i =
		     du[i__4].r * b[i__5].i + du[i__4].i * b[i__5].r;
	    z__3.r = b[i__3].r - z__4.r, z__3.i = b[i__3].i - z__4.i;
	    i__6 = k;
	    i__7 = k + 2 + j * b_dim1;
	    z__5.r = dl[i__6].r * b[i__7].r - dl[i__6].i * b[i__7].i, z__5.i =
		     dl[i__6].r * b[i__7].i + dl[i__6].i * b[i__7].r;
	    z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	    z_div(&z__1, &z__2, &d__[k]);
	    b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	}
    }
    return 0;
} /* zgtsvnp_ */

