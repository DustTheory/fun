/* zgttrfnp.f -- translated by f2c (version 20090411).
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

/* modified from zgttrf to avoid pivoting */
/* Subroutine */ int zgttrfnp_(integer *n, doublecomplex *dl, doublecomplex *
	d__, doublecomplex *du, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__;
    static doublecomplex fact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

    /* Parameter adjustments */
    --du;
    --d__;
    --dl;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("zgttrfnp", &i__1, (ftnlen)8);
	return 0;
    }
    if (*n == 0) {
	return 0;
    }
    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) != 0.) {
	    z_div(&z__1, &dl[i__], &d__[i__]);
	    fact.r = z__1.r, fact.i = z__1.i;
	    i__2 = i__;
	    dl[i__2].r = fact.r, dl[i__2].i = fact.i;
	    i__2 = i__ + 1;
	    i__3 = i__ + 1;
	    i__4 = i__;
	    z__2.r = fact.r * du[i__4].r - fact.i * du[i__4].i, z__2.i = 
		    fact.r * du[i__4].i + fact.i * du[i__4].r;
	    z__1.r = d__[i__3].r - z__2.r, z__1.i = d__[i__3].i - z__2.i;
	    d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
	}
    }
    if (*n > 1) {
	i__ = *n - 1;
	i__1 = i__;
	if ((d__1 = d__[i__1].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) != 0.) {
	    z_div(&z__1, &dl[i__], &d__[i__]);
	    fact.r = z__1.r, fact.i = z__1.i;
	    i__1 = i__;
	    dl[i__1].r = fact.r, dl[i__1].i = fact.i;
	    i__1 = i__ + 1;
	    i__2 = i__ + 1;
	    i__3 = i__;
	    z__2.r = fact.r * du[i__3].r - fact.i * du[i__3].i, z__2.i = 
		    fact.r * du[i__3].i + fact.i * du[i__3].r;
	    z__1.r = d__[i__2].r - z__2.r, z__1.i = d__[i__2].i - z__2.i;
	    d__[i__1].r = z__1.r, d__[i__1].i = z__1.i;
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	if ((d__1 = d__[i__2].r, abs(d__1)) + (d__2 = d_imag(&d__[i__]), abs(
		d__2)) == 0.) {
	    *info = i__;
	    goto L50;
	}
    }
L50:
    return 0;
} /* zgttrfnp_ */

