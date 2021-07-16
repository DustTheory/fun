/* zgbtf2np.f -- translated by f2c (version 20090411).
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

static doublecomplex c_b1 = {1.,0.};
static integer c__1 = 1;

/* modified from zgbtf2 to avoid pivoting */
/* Subroutine */ int zgbtf2np_(integer *m, integer *n, integer *kl, integer *
	ku, doublecomplex *ab, integer *ldab, integer *info)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, km, ju;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zgeru_(integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), xerbla_(char *, integer *,
	     ftnlen);

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
	xerbla_("zgbtf2np", &i__1, (ftnlen)8);
	return 0;
    }
    if (*m == 0 || *n == 0) {
	return 0;
    }
/*     do j=ku+2,min(ku,n) */
/*       do i=ku-j+2,kl */
/*         ab(i,j)=zero */
/*       enddo */
/*     enddo */
    ju = 1;
    i__1 = min(*m,*n);
    for (j = 1; j <= i__1; ++j) {
/*       if (j+ku.le.n) then */
/*         do i=1,kl */
/*           ab(i,j+ku)=zero */
/*         enddo */
/*       endif */
/* Computing MIN */
	i__2 = *kl, i__3 = *m - j;
	km = min(i__2,i__3);
	i__2 = *ku + 1 + j * ab_dim1;
	if (ab[i__2].r != 0. || ab[i__2].i != 0.) {
/* Computing MAX */
/* Computing MIN */
	    i__4 = j + *ku;
	    i__2 = ju, i__3 = min(i__4,*n);
	    ju = max(i__2,i__3);
	    if (km > 0) {
		z_div(&z__1, &c_b1, &ab[*ku + 1 + j * ab_dim1]);
		zscal_(&km, &z__1, &ab[*ku + 2 + j * ab_dim1], &c__1);
		if (ju > j) {
		    i__2 = ju - j;
		    z__1.r = -1., z__1.i = -0.;
		    i__3 = *ldab - 1;
		    i__4 = *ldab - 1;
		    zgeru_(&km, &i__2, &z__1, &ab[*ku + 2 + j * ab_dim1], &
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
} /* zgbtf2np_ */

