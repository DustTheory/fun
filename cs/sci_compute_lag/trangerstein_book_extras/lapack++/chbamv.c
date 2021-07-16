/* chbamv.f -- translated by f2c (version 20090411).
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

/* modifed from chbmv to work with absolute values */
/* Subroutine */ int chbamv_(char *uplo, integer *n, integer *k, real *alpha, 
	complex *a, integer *lda, complex *x, integer *incx, real *beta, real 
	*y, integer *incy, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double c_abs(complex *);

    /* Local variables */
    static integer i__, j, l, ix, iy, jx, jy, kx, ky, info;
    static real temp1, temp2;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer kplus1;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U", (ftnlen)1, (ftnlen)1) && ! lsame_(uplo, "L", (
	    ftnlen)1, (ftnlen)1)) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*k < 0) {
	info = 3;
    } else if (*lda < *k + 1) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("chbamv ", &info, (ftnlen)7);
	return 0;
    }
    if (*n == 0 || *alpha == 0.f && *beta == 1.f) {
	return 0;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (*n - 1) * *incy;
    }
    if (*beta != 1.f) {
	if (*incy == 1) {
	    if (*beta == 0.f) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = i__;
		    y[i__2] = 0.f;
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.f) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = iy;
		    y[i__2] = 0.f;
		    iy += *incy;
		}
	    } else {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
		}
	    }
	}
    }
    if (*alpha == 0.f) {
	return 0;
    }
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
	kplus1 = *k + 1;
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * c_abs(&x[j]);
		temp2 = 0.f;
		l = kplus1 - j;
/* Computing MAX */
		i__2 = 1, i__3 = j - *k;
		i__4 = j - 1;
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		    y[i__] += temp1 * c_abs(&a[l + i__ + j * a_dim1]);
		    temp2 += c_abs(&a[l + i__ + j * a_dim1]) * c_abs(&x[i__]);
		}
		y[j] = y[j] + temp1 * c_abs(&a[kplus1 + j * a_dim1]) + *alpha 
			* temp2;
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * c_abs(&x[jx]);
		temp2 = 0.f;
		ix = kx;
		iy = ky;
		l = kplus1 - j;
/* Computing MAX */
		i__4 = 1, i__2 = j - *k;
		i__3 = j - 1;
		for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
		    y[iy] += temp1 * c_abs(&a[l + i__ + j * a_dim1]);
		    temp2 += c_abs(&a[l + i__ + j * a_dim1]) * c_abs(&x[ix]);
		    ix += *incx;
		    iy += *incy;
		}
		y[jy] = y[jy] + temp1 * c_abs(&a[kplus1 + j * a_dim1]) + *
			alpha * temp2;
		jx += *incx;
		jy += *incy;
		if (j > *k) {
		    kx += *incx;
		    ky += *incy;
		}
	    }
	}
    } else {
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * c_abs(&x[j]);
		temp2 = 0.f;
		y[j] += temp1 * c_abs(&a[j * a_dim1 + 1]);
		l = 1 - j;
/* Computing MIN */
		i__4 = *n, i__2 = j + *k;
		i__3 = min(i__4,i__2);
		for (i__ = j + 1; i__ <= i__3; ++i__) {
		    y[i__] += temp1 * c_abs(&a[l + i__ + j * a_dim1]);
		    temp2 += c_abs(&a[l + i__ + j * a_dim1]) * c_abs(&x[i__]);
		}
		y[j] += *alpha * temp2;
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp1 = *alpha * c_abs(&x[jx]);
		temp2 = 0.f;
		y[jy] += temp1 * c_abs(&a[j * a_dim1 + 1]);
		l = 1 - j;
		ix = jx;
		iy = jy;
/* Computing MIN */
		i__4 = *n, i__2 = j + *k;
		i__3 = min(i__4,i__2);
		for (i__ = j + 1; i__ <= i__3; ++i__) {
		    ix += *incx;
		    iy += *incy;
		    y[iy] += temp1 * c_abs(&a[l + i__ + j * a_dim1]);
		    temp2 += c_abs(&a[l + i__ + j * a_dim1]) * c_abs(&x[ix]);
		}
		y[jy] += *alpha * temp2;
		jx += *incx;
		jy += *incy;
	    }
	}
    }
    return 0;
} /* chbamv_ */

