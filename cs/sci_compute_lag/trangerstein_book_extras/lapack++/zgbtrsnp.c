/* zgbtrsnp.f -- translated by f2c (version 20090411).
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

/* modified from zgbtrs to avoid pivoting */
/* Subroutine */ int zgbtrsnp_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *b,
	 integer *ldb, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, kd, lm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lnoti;
    extern /* Subroutine */ int zgemv_(char *, integer *, integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, ftnlen), 
	    zgeru_(integer *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *)
	    , ztbsv_(char *, char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen), zlacgv_(
	    integer *, doublecomplex *, integer *);
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
    if (! notran && ! lsame_(trans, "t", (ftnlen)1, (ftnlen)1) && ! lsame_(
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
	xerbla_("zgbtrsnp", &i__1, (ftnlen)8);
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
		z__1.r = -1., z__1.i = -0.;
		zgeru_(&lm, nrhs, &z__1, &ab[kd + 1 + j * ab_dim1], &c__1, &b[
			j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
	    }
	}
	i__1 = *nrhs;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ztbsv_("Upper", "No transpose", "Non-unit", n, ku, &ab[ab_offset],
		     ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)12,
		     (ftnlen)8);
	}
    } else if (lsame_(trans, "T", (ftnlen)1, (ftnlen)1)) {
	i__1 = *nrhs;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ztbsv_("Upper", "Transpose", "Non-unit", n, ku, &ab[ab_offset], 
		    ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)9, (
		    ftnlen)8);
	}
	if (lnoti) {
	    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		i__1 = *kl, i__2 = *n - j;
		lm = min(i__1,i__2);
		z__1.r = -1., z__1.i = -0.;
		zgemv_("Transpose", &lm, nrhs, &z__1, &b[j + 1 + b_dim1], ldb,
			 &ab[kd + 1 + j * ab_dim1], &c__1, &c_b1, &b[j + 
			b_dim1], ldb, (ftnlen)9);
	    }
	}
    } else {
	i__1 = *nrhs;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ztbsv_("Upper", "Conjugate transpose", "Non-unit", n, ku, &ab[
		    ab_offset], ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, 
		    (ftnlen)19, (ftnlen)8);
	}
	if (lnoti) {
	    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		i__1 = *kl, i__2 = *n - j;
		lm = min(i__1,i__2);
		zlacgv_(nrhs, &b[j + b_dim1], ldb);
		z__1.r = -1., z__1.i = -0.;
		zgemv_("Conjugate transpose", &lm, nrhs, &z__1, &b[j + 1 + 
			b_dim1], ldb, &ab[kd + 1 + j * ab_dim1], &c__1, &c_b1,
			 &b[j + b_dim1], ldb, (ftnlen)19);
		zlacgv_(nrhs, &b[j + b_dim1], ldb);
	    }
	}
    }
    return 0;
} /* zgbtrsnp_ */

