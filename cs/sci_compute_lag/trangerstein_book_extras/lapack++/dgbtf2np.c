/* dgbtf2np.f -- translated by f2c (version 20090411).
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

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b4 = -1.;

/* modified from dgbtf2 to avoid pivoting */
/* Subroutine */ int dgbtf2np_(integer *m, integer *n, integer *kl, integer *
	ku, doublereal *ab, integer *ldab, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;

    /* Local variables */
    static integer j, km, ju;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *), xerbla_(char *, integer *, ftnlen);

    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*ldab < *kl + *ku + 1) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("dgbtf2np", &i__1, (ftnlen)8);
	return 0;
    }
    if (*m == 0 || *n == 0) {
	return 0;
    }
    ju = 1;
    i__1 = min(*m,*n);
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	i__2 = *kl, i__3 = *m - j;
	km = min(i__2,i__3);
	if (ab[*ku + j + j * ab_dim1] != 0.) {
/* Computing MAX */
/* Computing MIN */
	    i__4 = j + *ku + j - 1;
	    i__2 = ju, i__3 = min(i__4,*n);
	    ju = max(i__2,i__3);
	    if (km > 0) {
		d__1 = 1. / ab[*ku + 1 + j * ab_dim1];
		dscal_(&km, &d__1, &ab[*ku + 2 + j * ab_dim1], &c__1);
		if (ju > j) {
		    i__2 = ju - j;
		    i__3 = *ldab - 1;
		    i__4 = *ldab - 1;
		    dger_(&km, &i__2, &c_b4, &ab[*ku + 2 + j * ab_dim1], &
			    c__1, &ab[*ku + (j + 1) * ab_dim1], &i__3, &ab[*
			    ku + 1 + (j + 1) * ab_dim1], &i__4);
		}
	    }
	} else {
	    if (*info == 0) {
		*info = j;
	    }
	}
    }
    return 0;
} /* dgbtf2np_ */

