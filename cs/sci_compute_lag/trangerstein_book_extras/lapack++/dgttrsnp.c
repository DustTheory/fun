/* dgttrsnp.f -- translated by f2c (version 20090411).
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
static integer c_n1 = -1;

/*     modified from dgttrs to avoid pivoting */
/* Subroutine */ int dgttrsnp_(char *trans, integer *n, integer *nrhs, 
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *b, 
	integer *ldb, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, jb, nb;
    extern /* Subroutine */ int dgtts2_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer itrans;
    static logical notran;

    /* Parameter adjustments */
    --dl;
    --d__;
    --du;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    notran = *(unsigned char *)trans == 'N' || *(unsigned char *)trans == 'n';
    if (! notran && ! (*(unsigned char *)trans == 'T' || *(unsigned char *)
	    trans == 't') && ! (*(unsigned char *)trans == 'C' || *(unsigned 
	    char *)trans == 'c')) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < max(*n,1)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("dgttrsnp", &i__1, (ftnlen)8);
	return 0;
    }
    if (*n == 0 || *nrhs == 0) {
	return 0;
    }
    if (notran) {
	itrans = 0;
    } else {
	itrans = 1;
    }
    if (*nrhs == 1) {
	nb = 1;
    } else {
/* Computing MAX */
	i__1 = 1, i__2 = ilaenv_(&c__1, "dgttrs", trans, n, nrhs, &c_n1, &
		c_n1, (ftnlen)6, (ftnlen)1);
	nb = max(i__1,i__2);
    }
    if (nb >= *nrhs) {
	dgtts2_(&itrans, n, nrhs, &dl[1], &d__[1], &du[1], &b[b_offset], ldb);
    } else {
	i__1 = *nrhs;
	i__2 = nb;
	for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
/* Computing MIN */
	    i__3 = *nrhs - j + 1;
	    jb = min(i__3,nb);
	    dgtts2_(&itrans, n, &jb, &dl[1], &d__[1], &du[1], &b[j * b_dim1 + 
		    1], ldb);
	}
    }
    return 0;
} /* dgttrsnp_ */

