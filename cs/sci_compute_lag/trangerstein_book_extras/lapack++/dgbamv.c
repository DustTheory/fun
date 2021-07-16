/* dgbamv.f -- translated by f2c (version 20090411).
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

/* modified from dgbmv to use abs(a) and abs(y) */
/* Subroutine */ int dgbamv_(char *trans, integer *m, integer *n, integer *kl,
	 integer *ku, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *x, integer *incx, doublereal *beta, doublereal *y, 
	integer *incy, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, ix, iy, jx, jy, kx, ky, kup1, info;
    static doublereal temp;
    static integer lenx, leny;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

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
	xerbla_("dgbamv ", &info, (ftnlen)7);
	return 0;
    }
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }
    if (lsame_(trans, "n", (ftnlen)1, (ftnlen)1)) {
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
    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
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
    if (*alpha == 0.) {
	return 0;
    }
    kup1 = *ku + 1;
    if (lsame_(trans, "n", (ftnlen)1, (ftnlen)1)) {
	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * (d__1 = x[jx], abs(d__1));
		    k = kup1 - j;
/* Computing MAX */
		    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__4 = min(i__5,i__6);
		    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
			y[i__] += temp * (d__1 = a[k + i__ + j * a_dim1], abs(
				d__1));
		    }
		}
		jx += *incx;
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * (d__1 = x[jx], abs(d__1));
		    iy = ky;
		    k = kup1 - j;
/* Computing MAX */
		    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__3 = min(i__5,i__6);
		    for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
			y[iy] += temp * (d__1 = a[k + i__ + j * a_dim1], abs(
				d__1));
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
		temp = 0.;
		k = kup1 - j;
/* Computing MAX */
		i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
		i__5 = *m, i__6 = j + *kl;
		i__2 = min(i__5,i__6);
		for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
		    temp += (d__1 = a[k + i__ + j * a_dim1] * x[i__], abs(
			    d__1));
		}
		y[jy] += *alpha * temp;
		jy += *incy;
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		ix = kx;
		k = kup1 - j;
/* Computing MAX */
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		i__5 = *m, i__6 = j + *kl;
		i__4 = min(i__5,i__6);
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		    temp += (d__1 = a[k + i__ + j * a_dim1] * x[ix], abs(d__1)
			    );
		    ix += *incx;
		}
		y[jy] += *alpha * temp;
		jy += *incy;
		if (j > *ku) {
		    kx += *incx;
		}
	    }
	}
    }
    return 0;
} /* dgbamv_ */

