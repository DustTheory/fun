/* dla_geamv.f -- translated by f2c (version 20090411).
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

/* > \brief \b DLA_GEAMV */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DLA_GEAMV + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dla_gea
mv.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dla_gea
mv.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dla_gea
mv.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DLA_GEAMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, */
/*                              Y, INCY ) */

/*       .. Scalar Arguments .. */
/*       DOUBLE PRECISION   ALPHA, BETA */
/*       INTEGER            INCX, INCY, LDA, M, N, TRANS */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DLA_GEAMV  performs one of the matrix-vector operations */
/* > */
/* >         y := alpha*abs(A)*abs(x) + beta*abs(y), */
/* >    or   y := alpha*abs(A)**T*abs(x) + beta*abs(y), */
/* > */
/* > where alpha and beta are scalars, x and y are vectors and A is an */
/* > m by n matrix. */
/* > */
/* > This function is primarily used in calculating error bounds. */
/* > To protect against underflow during evaluation, components in */
/* > the resulting vector are perturbed away from zero by (N+1) */
/* > times the underflow threshold.  To prevent unnecessarily large */
/* > errors for block-structure embedded in general matrices, */
/* > "symbolically" zero components are not perturbed.  A zero */
/* > entry is considered "symbolic" if all multiplications involved */
/* > in computing that entry have at least one zero multiplicand. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] TRANS */
/* > \verbatim */
/* >          TRANS is INTEGER */
/* >           On entry, TRANS specifies the operation to be performed as */
/* >           follows: */
/* > */
/* >             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y) */
/* >             BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y) */
/* >             BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y) */
/* > */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >           On entry, M specifies the number of rows of the matrix A. */
/* >           M must be at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >           On entry, N specifies the number of columns of the matrix A. */
/* >           N must be at least zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] ALPHA */
/* > \verbatim */
/* >          ALPHA is DOUBLE PRECISION */
/* >           On entry, ALPHA specifies the scalar alpha. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ) */
/* >           Before entry, the leading m by n part of the array A must */
/* >           contain the matrix of coefficients. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >           On entry, LDA specifies the first dimension of A as declared */
/* >           in the calling (sub) program. LDA must be at least */
/* >           max( 1, m ). */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] X */
/* > \verbatim */
/* >          X is DOUBLE PRECISION array, dimension */
/* >           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/* >           Before entry, the incremented array X must contain the */
/* >           vector x. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] INCX */
/* > \verbatim */
/* >          INCX is INTEGER */
/* >           On entry, INCX specifies the increment for the elements of */
/* >           X. INCX must not be zero. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in] BETA */
/* > \verbatim */
/* >          BETA is DOUBLE PRECISION */
/* >           On entry, BETA specifies the scalar beta. When BETA is */
/* >           supplied as zero then Y need not be set on input. */
/* >           Unchanged on exit. */
/* > \endverbatim */
/* > */
/* > \param[in,out] Y */
/* > \verbatim */
/* >          Y is DOUBLE PRECISION */
/* >           Array of DIMENSION at least */
/* >           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/* >           and at least */
/* >           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/* >           Before entry with BETA non-zero, the incremented array Y */
/* >           must contain the vector y. On exit, Y is overwritten by the */
/* >           updated vector y. */
/* > \endverbatim */
/* > */
/* > \param[in] INCY */
/* > \verbatim */
/* >          INCY is INTEGER */
/* >           On entry, INCY specifies the increment for the elements of */
/* >           Y. INCY must not be zero. */
/* >           Unchanged on exit. */
/* > */
/* >  Level 2 Blas routine. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup doubleGEcomputational */

/*  ===================================================================== */
/* Subroutine */ int dla_geamv__(integer *trans, integer *m, integer *n, 
	doublereal *alpha, doublereal *a, integer *lda, doublereal *x, 
	integer *incx, doublereal *beta, doublereal *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    extern integer ilatrans_(char *, ftnlen);
    static integer i__, j;
    static logical symb_zero__;
    static integer iy, jx, kx, ky, info;
    static doublereal temp;
    static integer lenx, leny;
    static doublereal safe1;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     November 2011 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! (*trans == ilatrans_("N", (ftnlen)1) || *trans == ilatrans_("T", (
	    ftnlen)1) || *trans == ilatrans_("C", (ftnlen)1))) {
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
	xerbla_("DLA_GEAMV ", &info, (ftnlen)10);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    if (*trans == ilatrans_("N", (ftnlen)1)) {
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

/*     Set SAFE1 essentially to be the underflow threshold times the */
/*     number of additions in each row. */

    safe1 = dlamch_("Safe minimum", (ftnlen)12);
    safe1 = (*n + 1) * safe1;

/*     Form  y := alpha*abs(A)*abs(x) + beta*abs(y). */

/*     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to */
/*     the inexact flag.  Still doesn't help change the iteration order */
/*     to per-column. */

    iy = ky;
    if (*incx == 1) {
	if (*trans == ilatrans_("N", (ftnlen)1)) {
	    i__1 = leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*beta == 0.) {
		    symb_zero__ = TRUE_;
		    y[iy] = 0.;
		} else if (y[iy] == 0.) {
		    symb_zero__ = TRUE_;
		} else {
		    symb_zero__ = FALSE_;
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
		}
		if (*alpha != 0.) {
		    i__2 = lenx;
		    for (j = 1; j <= i__2; ++j) {
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
		    }
		}
		if (! symb_zero__) {
		    y[iy] += d_sign(&safe1, &y[iy]);
		}
		iy += *incy;
	    }
	} else {
	    i__1 = leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*beta == 0.) {
		    symb_zero__ = TRUE_;
		    y[iy] = 0.;
		} else if (y[iy] == 0.) {
		    symb_zero__ = TRUE_;
		} else {
		    symb_zero__ = FALSE_;
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
		}
		if (*alpha != 0.) {
		    i__2 = lenx;
		    for (j = 1; j <= i__2; ++j) {
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[j] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[j], abs(d__1)) * temp;
		    }
		}
		if (! symb_zero__) {
		    y[iy] += d_sign(&safe1, &y[iy]);
		}
		iy += *incy;
	    }
	}
    } else {
	if (*trans == ilatrans_("N", (ftnlen)1)) {
	    i__1 = leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*beta == 0.) {
		    symb_zero__ = TRUE_;
		    y[iy] = 0.;
		} else if (y[iy] == 0.) {
		    symb_zero__ = TRUE_;
		} else {
		    symb_zero__ = FALSE_;
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
		}
		if (*alpha != 0.) {
		    jx = kx;
		    i__2 = lenx;
		    for (j = 1; j <= i__2; ++j) {
			temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[jx] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
			jx += *incx;
		    }
		}
		if (! symb_zero__) {
		    y[iy] += d_sign(&safe1, &y[iy]);
		}
		iy += *incy;
	    }
	} else {
	    i__1 = leny;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		if (*beta == 0.) {
		    symb_zero__ = TRUE_;
		    y[iy] = 0.;
		} else if (y[iy] == 0.) {
		    symb_zero__ = TRUE_;
		} else {
		    symb_zero__ = FALSE_;
		    y[iy] = *beta * (d__1 = y[iy], abs(d__1));
		}
		if (*alpha != 0.) {
		    jx = kx;
		    i__2 = lenx;
		    for (j = 1; j <= i__2; ++j) {
			temp = (d__1 = a[j + i__ * a_dim1], abs(d__1));
			symb_zero__ = symb_zero__ && (x[jx] == 0. || temp == 
				0.);
			y[iy] += *alpha * (d__1 = x[jx], abs(d__1)) * temp;
			jx += *incx;
		    }
		}
		if (! symb_zero__) {
		    y[iy] += d_sign(&safe1, &y[iy]);
		}
		iy += *incy;
	    }
	}
    }

    return 0;

/*     End of DLA_GEAMV */

} /* dla_geamv__ */
