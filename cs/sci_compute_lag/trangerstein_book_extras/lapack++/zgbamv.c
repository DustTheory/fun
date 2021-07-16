/* zgbamv.f -- translated by f2c (version 20090411).
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

/* modifed from zgbmv to use absolute values */
/* Subroutine */ int zgbamv_(char *trans, integer *m, integer *n, integer *kl,
	 integer *ku, doublecomplex *alpha, doublecomplex *a, integer *lda, 
	doublecomplex *x, integer *incx, doublecomplex *beta, doublereal *y, 
	integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, ix, iy, jx, jy, kx, ky, kup1, info;
    static doublecomplex temp;
    static integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconj;

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
	    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
	    ) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*kl < 0) {
	info = 4;
    } else if (*ku < 0) {
	info = 5;
    } else if (*lda < *kl + *ku + 1) {
	info = 8;
    } else if (*incx == 0) {
	info = 10;
    } else if (*incy == 0) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("zgbamv ", &info, (ftnlen)7);
	return 0;
    }
    if (*m == 0 || *n == 0 || alpha->r == 0. && alpha->i == 0. && (beta->r == 
	    1. && beta->i == 0.)) {
	return 0;
    }
    noconj = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }
    if (beta->r != 1. || beta->i != 0.) {
	if (*incy == 1) {
	    if (beta->r == 0. && beta->i == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = i__;
		    y[i__2] = 0.;
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = i__;
		    i__3 = i__;
		    z__1.r = y[i__3] * beta->r, z__1.i = y[i__3] * beta->i;
		    y[i__2] = z__1.r;
		}
	    }
	} else {
	    iy = ky;
	    if (beta->r == 0. && beta->i == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = iy;
		    y[i__2] = 0.;
		    iy += *incy;
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = iy;
		    i__3 = iy;
		    z__1.r = y[i__3] * beta->r, z__1.i = y[i__3] * beta->i;
		    y[i__2] = z__1.r;
		    iy += *incy;
		}
	    }
	}
    }
    if (alpha->r == 0. && alpha->i == 0.) {
	return 0;
    }
    kup1 = *ku + 1;
    if (lsame_(trans, "n", (ftnlen)1, (ftnlen)1)) {
	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = jx;
		if (x[i__2].r != 0. || x[i__2].i != 0.) {
		    d__1 = z_abs(&x[jx]);
		    z__1.r = d__1 * alpha->r, z__1.i = d__1 * alpha->i;
		    temp.r = z__1.r, temp.i = z__1.i;
		    k = kup1 - j;
/* Computing MAX */
		    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__4 = min(i__5,i__6);
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
			i__2 = i__;
			i__3 = i__;
			d__1 = z_abs(&a[k + i__ + j * a_dim1]);
			z__2.r = d__1 * temp.r, z__2.i = d__1 * temp.i;
			z__1.r = y[i__3] + z__2.r, z__1.i = z__2.i;
			y[i__2] = z__1.r;
		    }
		}
		jx += *incx;
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__4 = jx;
		if (x[i__4].r != 0. || x[i__4].i != 0.) {
		    d__1 = z_abs(&x[jx]);
		    z__1.r = d__1 * alpha->r, z__1.i = d__1 * alpha->i;
		    temp.r = z__1.r, temp.i = z__1.i;
		    iy = ky;
		    k = kup1 - j;
/* Computing MAX */
		    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__3 = min(i__5,i__6);
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
			i__4 = iy;
			i__2 = iy;
			d__1 = z_abs(&a[k + i__ + j * a_dim1]);
			z__2.r = d__1 * temp.r, z__2.i = d__1 * temp.i;
			z__1.r = y[i__2] + z__2.r, z__1.i = z__2.i;
			y[i__4] = z__1.r;
			iy += *incy;
		    }
		}
		jx += *incx;
		if (j > *ku) {
		    ky += *incy;
		}
	    }
	}
    } else {
	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp.r = 0., temp.i = 0.;
		k = kup1 - j;
		if (noconj) {
/* Computing MAX */
		    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__2 = min(i__5,i__6);
		    for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
			i__3 = k + i__ + j * a_dim1;
			i__4 = i__;
			z__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, z__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
			d__1 = z_abs(&z__2);
			z__1.r = temp.r + d__1, z__1.i = temp.i;
			temp.r = z__1.r, temp.i = z__1.i;
		    }
		} else {
/* Computing MAX */
		    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__4 = min(i__5,i__6);
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
			d_cnjg(&z__3, &a[k + i__ + j * a_dim1]);
			i__2 = i__;
			z__2.r = z__3.r * x[i__2].r - z__3.i * x[i__2].i, 
				z__2.i = z__3.r * x[i__2].i + z__3.i * x[i__2]
				.r;
			d__1 = z_abs(&z__2);
			z__1.r = temp.r + d__1, z__1.i = temp.i;
			temp.r = z__1.r, temp.i = z__1.i;
		    }
		}
		i__4 = jy;
		i__2 = jy;
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		z__1.r = y[i__2] + z__2.r, z__1.i = z__2.i;
		y[i__4] = z__1.r;
		jy += *incy;
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp.r = 0., temp.i = 0.;
		ix = kx;
		k = kup1 - j;
		if (noconj) {
/* Computing MAX */
		    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__3 = min(i__5,i__6);
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
			i__4 = k + i__ + j * a_dim1;
			i__2 = ix;
			z__2.r = a[i__4].r * x[i__2].r - a[i__4].i * x[i__2]
				.i, z__2.i = a[i__4].r * x[i__2].i + a[i__4]
				.i * x[i__2].r;
			d__1 = z_abs(&z__2);
			z__1.r = temp.r + d__1, z__1.i = temp.i;
			temp.r = z__1.r, temp.i = z__1.i;
			ix += *incx;
		    }
		} else {
/* Computing MAX */
		    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__2 = min(i__5,i__6);
		    for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
			d_cnjg(&z__3, &a[k + i__ + j * a_dim1]);
			i__3 = ix;
			z__2.r = z__3.r * x[i__3].r - z__3.i * x[i__3].i, 
				z__2.i = z__3.r * x[i__3].i + z__3.i * x[i__3]
				.r;
			d__1 = z_abs(&z__2);
			z__1.r = temp.r + d__1, z__1.i = temp.i;
			temp.r = z__1.r, temp.i = z__1.i;
			ix += *incx;
		    }
		}
		i__2 = jy;
		i__3 = jy;
		z__2.r = alpha->r * temp.r - alpha->i * temp.i, z__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		z__1.r = y[i__3] + z__2.r, z__1.i = z__2.i;
		y[i__2] = z__1.r;
		jy += *incy;
		if (j > *ku) {
		    kx += *incx;
		}
	    }
	}
    }
    return 0;
} /* zgbamv_ */

