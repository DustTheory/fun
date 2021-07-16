/* sgbtrsnp.f -- translated by f2c (version 20090411).
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

static real c_b6 = -1.f;
static integer c__1 = 1;
static real c_b19 = 1.f;

/* modified from sgbtrs to avoid pivoting */
/* Subroutine */ int sgbtrsnp_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, real *ab, integer *ldab, real *b, integer *ldb, 
	integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, kd, lm;
    extern /* Subroutine */ int sger_(integer *, integer *, real *, real *, 
	    integer *, real *, integer *, real *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sgemv_(char *, integer *, integer *, real *, 
	    real *, integer *, real *, integer *, real *, real *, integer *, 
	    ftnlen);
    static logical lnoti;
    extern /* Subroutine */ int stbsv_(char *, char *, char *, integer *, 
	    integer *, real *, integer *, real *, integer *, ftnlen, ftnlen, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical notran;

    /* Parameter adjustments */
    ab_dim1 = *ldab;
    ab_offset = 1 + ab_dim1;
    ab -= ab_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*nrhs < 0) {
	*info = -5;
    } else if (*ldab < *kl + *ku + 1) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("sgbtrsnp", &i__1, (ftnlen)8);
	return 0;
    }
    if (*n == 0 || *nrhs == 0) {
	return 0;
    }
    kd = *ku + 1;
    lnoti = *kl > 0;
    if (notran) {
	if (lnoti) {
	    i__1 = *n - 1;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
		i__2 = *kl, i__3 = *n - j;
		lm = min(i__2,i__3);
		sger_(&lm, nrhs, &c_b6, &ab[kd + 1 + j * ab_dim1], &c__1, &b[
			j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
	    }
	}
	i__1 = *nrhs;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    stbsv_("Upper", "No transpose", "Non-unit", n, ku, &ab[ab_offset],
		     ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12,
		     (ftnlen)8);
	}
    } else {
	i__1 = *nrhs;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    stbsv_("Upper", "Transpose", "Non-unit", n, ku, &ab[ab_offset], 
		    ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)9, (
		    ftnlen)8);
	}
	if (lnoti) {
	    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		i__1 = *kl, i__2 = *n - j;
		lm = min(i__1,i__2);
		sgemv_("Transpose", &lm, nrhs, &c_b6, &b[j + 1 + b_dim1], ldb,
			 &ab[kd + 1 + j * ab_dim1], &c__1, &c_b19, &b[j + 
			b_dim1], ldb, (ftnlen)9);
	    }
	}
    }
    return 0;
} /* sgbtrsnp_ */

