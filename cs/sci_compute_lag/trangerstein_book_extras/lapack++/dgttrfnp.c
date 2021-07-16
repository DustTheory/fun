/* dgttrfnp.f -- translated by f2c (version 20090411).
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

/*     modified from dgttrf to avoid pivoting */
/* Subroutine */ int dgttrfnp_(integer *n, doublereal *dl, doublereal *d__, 
	doublereal *du, integer *info)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static doublereal fact;
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
	xerbla_("dgttrfnp", &i__1, (ftnlen)8);
	return 0;
    }
    if (*n == 0) {
	return 0;
    }
    i__1 = *n - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] != 0.) {
	    fact = dl[i__] / d__[i__];
	    dl[i__] = fact;
	    d__[i__ + 1] -= fact * du[i__];
	}
    }
    if (*n > 1) {
	i__ = *n - 1;
	if (d__[i__] != 0.) {
	    fact = dl[i__] / d__[i__];
	    dl[i__] = fact;
	    d__[i__ + 1] -= fact * du[i__];
	}
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (d__[i__] == 0.) {
	    *info = i__;
	    goto L50;
	}
    }
L50:
    return 0;
} /* dgttrfnp_ */

