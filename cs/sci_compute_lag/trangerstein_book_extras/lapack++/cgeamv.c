/* cgeamv.f -- translated by f2c (version 20090411).
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

/* modified from cgemv to use absolute values */
/* Subroutine */ int cgeamv_(char *trans, integer *m, integer *n, real *alpha,
	 complex *a, integer *lda, complex *x, integer *incx, real *beta, 
	real *y, integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    complex q__1, q__2;

    /* Builtin functions */
    double c_abs(complex *);
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer i__, j, ix, iy, jx, jy, kx, ky, info;
    static real temp;
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
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("cgeamv ", &info, (ftnlen)7);
	return 0;
    }
    if (*m == 0 || *n == 0 || *alpha == 0.f && *beta == 1.f) {
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
    if (*beta != 1.f) {
	if (*incy == 1) {
	    if (*beta == 0.f) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = i__;
		    y[i__2] = 0.f;
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.f) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = iy;
		    y[i__2] = 0.f;
		    iy += *incy;
		}
	    } else {
		i__1 = leny;
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
    if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = jx;
		if (x[i__2].r != 0.f || x[i__2].i != 0.f) {
		    temp = *alpha * c_abs(&x[jx]);
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[i__] += temp * c_abs(&a[i__ + j * a_dim1]);
		    }
		}
		jx += *incx;
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = jx;
		if (x[i__2].r != 0.f || x[i__2].i != 0.f) {
		    temp = *alpha * c_abs(&x[jx]);
		    iy = ky;
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[iy] += temp * c_abs(&a[i__ + j * a_dim1]);
			iy += *incy;
		    }
		}
		jx += *incx;
	    }
	}
    } else {
	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.f;
		if (noconj) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * a_dim1;
			i__4 = i__;
			q__1.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, q__1.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
			temp += c_abs(&q__1);
		    }
		} else {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			r_cnjg(&q__2, &a[i__ + j * a_dim1]);
			i__3 = i__;
			q__1.r = q__2.r * x[i__3].r - q__2.i * x[i__3].i, 
				q__1.i = q__2.r * x[i__3].i + q__2.i * x[i__3]
				.r;
			temp += c_abs(&q__1);
		    }
		}
		y[jy] += *alpha * temp;
		jy += *incy;
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.f;
		ix = kx;
		if (noconj) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * a_dim1;
			i__4 = ix;
			q__1.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, q__1.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
			temp += c_abs(&q__1);
			ix += *incx;
		    }
		} else {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			r_cnjg(&q__2, &a[i__ + j * a_dim1]);
			i__3 = ix;
			q__1.r = q__2.r * x[i__3].r - q__2.i * x[i__3].i, 
				q__1.i = q__2.r * x[i__3].i + q__2.i * x[i__3]
				.r;
			temp += c_abs(&q__1);
			ix += *incx;
		    }
		}
		y[jy] += *alpha * temp;
		jy += *incy;
	    }
	}
    }
    return 0;
} /* cgeamv_ */

