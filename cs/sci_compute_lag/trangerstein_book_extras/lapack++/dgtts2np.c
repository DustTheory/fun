/* dgtts2np.f -- translated by f2c (version 20090411).
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

/*     modified from dgtts2 to avoid pivoting */
/* Subroutine */ int dgtts2np_(integer *itrans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *b, 
	integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal temp;

    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    if (*n == 0 || *nrhs == 0) {
	return 0;
    }
    if (*itrans == 0) {
	if (*nrhs <= 1) {
	    j = 1;
L10:
	    i__1 = *n - 1;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		temp = b[i__ + 1 + j * b_dim1] - dl[i__] * b[i__ + j * b_dim1]
			;
		b[i__ + j * b_dim1] = b[i__ + j * b_dim1];
		b[i__ + 1 + j * b_dim1] = temp;
	    }
	    b[*n + j * b_dim1] /= d__[*n];
	    if (*n > 1) {
		b[*n - 1 + j * b_dim1] = (b[*n - 1 + j * b_dim1] - du[*n - 1] 
			* b[*n + j * b_dim1]) / d__[*n - 1];
	    }
	    for (i__ = *n - 2; i__ >= 1; --i__) {
		b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[i__ 
			+ 1 + j * b_dim1]) / d__[i__];
	    }
	    if (j < *nrhs) {
		++j;
		goto L10;
	    }
	} else {
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n - 1;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    b[i__ + 1 + j * b_dim1] -= dl[i__] * b[i__ + j * b_dim1];
		}
		b[*n + j * b_dim1] /= d__[*n];
		if (*n > 1) {
		    b[*n - 1 + j * b_dim1] = (b[*n - 1 + j * b_dim1] - du[*n 
			    - 1] * b[*n + j * b_dim1]) / d__[*n - 1];
		}
		for (i__ = *n - 2; i__ >= 1; --i__) {
		    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__] * b[
			    i__ + 1 + j * b_dim1]) / d__[i__];
		}
	    }
	}
    } else {
	if (*nrhs <= 1) {
	    j = 1;
L70:
	    b[j * b_dim1 + 1] /= d__[1];
	    if (*n > 1) {
		b[j * b_dim1 + 2] = (b[j * b_dim1 + 2] - du[1] * b[j * b_dim1 
			+ 1]) / d__[2];
	    }
	    i__1 = *n;
	    for (i__ = 3; i__ <= i__1; ++i__) {
		b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__ - 1] * b[
			i__ - 1 + j * b_dim1]) / d__[i__];
	    }
	    for (i__ = *n - 1; i__ >= 1; --i__) {
		temp = b[i__ + j * b_dim1] - dl[i__] * b[i__ + 1 + j * b_dim1]
			;
		b[i__ + j * b_dim1] = temp;
	    }
	    if (j < *nrhs) {
		++j;
		goto L70;
	    }
	} else {
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		b[j * b_dim1 + 1] /= d__[1];
		if (*n > 1) {
		    b[j * b_dim1 + 2] = (b[j * b_dim1 + 2] - du[1] * b[j * 
			    b_dim1 + 1]) / d__[2];
		}
		i__2 = *n;
		for (i__ = 3; i__ <= i__2; ++i__) {
		    b[i__ + j * b_dim1] = (b[i__ + j * b_dim1] - du[i__ - 1] *
			     b[i__ - 1 + j * b_dim1]) / d__[i__];
		}
		for (i__ = *n - 1; i__ >= 1; --i__) {
		    b[i__ + j * b_dim1] -= dl[i__] * b[i__ + 1 + j * b_dim1];
		}
	    }
	}
    }
    return 0;
} /* dgtts2np_ */

