/* dgtsvrnp.f -- translated by f2c (version 20090411).
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

/* modified from dgtsv to solve for rows of b, rather than columns */
/* Subroutine */ int dgtsvrnp_(integer *n, integer *nrhs, doublereal *dl, 
	doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer 
	*info)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal fact;
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
	xerbla_("dgtsvr", &i__1, (ftnlen)6);
	return 0;
    }
    if (*n == 0) {
	return 0;
    }
    if (*nrhs == 1) {
	i__1 = *n - 2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (d__[i__] != 0.) {
		fact = dl[i__] / d__[i__];
		d__[i__ + 1] -= fact * du[i__];
		b[(i__ + 1) * b_dim1 + 1] -= fact * b[i__ * b_dim1 + 1];
	    } else {
		*info = i__;
		return 0;
	    }
	    dl[i__] = 0.;
	}
	if (*n > 1) {
	    i__ = *n - 1;
	    if (d__[i__] != 0.) {
		fact = dl[i__] / d__[i__];
		d__[i__ + 1] -= fact * du[i__];
		b[(i__ + 1) * b_dim1 + 1] -= fact * b[i__ * b_dim1 + 1];
	    } else {
		*info = i__;
		return 0;
	    }
	}
	if (d__[*n] == 0.) {
	    *info = *n;
	    return 0;
	}
    } else {
	i__1 = *n - 2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (d__[i__] != 0.) {
		fact = dl[i__] / d__[i__];
		d__[i__ + 1] -= fact * du[i__];
		i__2 = *nrhs;
		for (j = 1; j <= i__2; ++j) {
		    b[j + (i__ + 1) * b_dim1] -= fact * b[j + i__ * b_dim1];
		}
	    } else {
		*info = i__;
		return 0;
	    }
	    dl[i__] = 0.;
	}
	if (*n > 1) {
	    i__ = *n - 1;
	    if (d__[i__] != 0.) {
		fact = dl[i__] / d__[i__];
		d__[i__ + 1] -= fact * du[i__];
		i__1 = *nrhs;
		for (j = 1; j <= i__1; ++j) {
		    b[j + (i__ + 1) * b_dim1] -= fact * b[j + i__ * b_dim1];
		}
	    } else {
		*info = i__;
		return 0;
	    }
	}
	if (d__[*n] == 0.) {
	    *info = *n;
	    return 0;
	}
    }
    if (*nrhs <= 2) {
	j = 1;
L70:
	b[j + *n * b_dim1] /= d__[*n];
	if (*n > 1) {
	    b[j + (*n - 1) * b_dim1] = (b[j + (*n - 1) * b_dim1] - du[*n - 1] 
		    * b[j + *n * b_dim1]) / d__[*n - 1];
	}
	for (i__ = *n - 2; i__ >= 1; --i__) {
	    b[j + i__ * b_dim1] = (b[j + i__ * b_dim1] - du[i__] * b[j + (i__ 
		    + 1) * b_dim1] - dl[i__] * b[j + (i__ + 2) * b_dim1]) / 
		    d__[i__];
	}
	if (j < *nrhs) {
	    ++j;
	    goto L70;
	}
    } else {
	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j) {
	    b[j + *n * b_dim1] /= d__[*n];
	    if (*n > 1) {
		b[j + (*n - 1) * b_dim1] = (b[j + (*n - 1) * b_dim1] - du[*n 
			- 1] * b[j + *n * b_dim1]) / d__[*n - 1];
	    }
	    for (i__ = *n - 2; i__ >= 1; --i__) {
		b[j + i__ * b_dim1] = (b[j + i__ * b_dim1] - du[i__] * b[j + (
			i__ + 1) * b_dim1] - dl[i__] * b[j + (i__ + 2) * 
			b_dim1]) / d__[i__];
	    }
	}
    }
    return 0;
} /* dgtsvrnp_ */

