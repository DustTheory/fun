/* cgtts2np.f -- translated by f2c (version 20090411).
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

/* modifed from cgtts2 to avoid pivoting */
/* Subroutine */ int cgtts2np_(integer *itrans, integer *n, integer *nrhs, 
	complex *dl, complex *d__, complex *du, complex *b, integer *ldb)
{
    /* System generated locals */
    integer b_dim1, b_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    static integer i__, j;

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
		i__2 = i__ + 1 + j * b_dim1;
		i__3 = i__ + 1 + j * b_dim1;
		i__4 = i__;
		i__5 = i__ + j * b_dim1;
		q__2.r = dl[i__4].r * b[i__5].r - dl[i__4].i * b[i__5].i, 
			q__2.i = dl[i__4].r * b[i__5].i + dl[i__4].i * b[i__5]
			.r;
		q__1.r = b[i__3].r - q__2.r, q__1.i = b[i__3].i - q__2.i;
		b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	    }
	    i__1 = *n + j * b_dim1;
	    c_div(&q__1, &b[*n + j * b_dim1], &d__[*n]);
	    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
	    if (*n > 1) {
		i__1 = *n - 1 + j * b_dim1;
		i__2 = *n - 1 + j * b_dim1;
		i__3 = *n - 1;
		i__4 = *n + j * b_dim1;
		q__3.r = du[i__3].r * b[i__4].r - du[i__3].i * b[i__4].i, 
			q__3.i = du[i__3].r * b[i__4].i + du[i__3].i * b[i__4]
			.r;
		q__2.r = b[i__2].r - q__3.r, q__2.i = b[i__2].i - q__3.i;
		c_div(&q__1, &q__2, &d__[*n - 1]);
		b[i__1].r = q__1.r, b[i__1].i = q__1.i;
	    }
	    for (i__ = *n - 2; i__ >= 1; --i__) {
		i__1 = i__ + j * b_dim1;
		i__2 = i__ + j * b_dim1;
		i__3 = i__;
		i__4 = i__ + 1 + j * b_dim1;
		q__3.r = du[i__3].r * b[i__4].r - du[i__3].i * b[i__4].i, 
			q__3.i = du[i__3].r * b[i__4].i + du[i__3].i * b[i__4]
			.r;
		q__2.r = b[i__2].r - q__3.r, q__2.i = b[i__2].i - q__3.i;
		c_div(&q__1, &q__2, &d__[i__]);
		b[i__1].r = q__1.r, b[i__1].i = q__1.i;
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
		    i__3 = i__ + 1 + j * b_dim1;
		    i__4 = i__ + 1 + j * b_dim1;
		    i__5 = i__;
		    i__6 = i__ + j * b_dim1;
		    q__2.r = dl[i__5].r * b[i__6].r - dl[i__5].i * b[i__6].i, 
			    q__2.i = dl[i__5].r * b[i__6].i + dl[i__5].i * b[
			    i__6].r;
		    q__1.r = b[i__4].r - q__2.r, q__1.i = b[i__4].i - q__2.i;
		    b[i__3].r = q__1.r, b[i__3].i = q__1.i;
		}
		i__2 = *n + j * b_dim1;
		c_div(&q__1, &b[*n + j * b_dim1], &d__[*n]);
		b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		if (*n > 1) {
		    i__2 = *n - 1 + j * b_dim1;
		    i__3 = *n - 1 + j * b_dim1;
		    i__4 = *n - 1;
		    i__5 = *n + j * b_dim1;
		    q__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, 
			    q__3.i = du[i__4].r * b[i__5].i + du[i__4].i * b[
			    i__5].r;
		    q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
		    c_div(&q__1, &q__2, &d__[*n - 1]);
		    b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		}
		for (i__ = *n - 2; i__ >= 1; --i__) {
		    i__2 = i__ + j * b_dim1;
		    i__3 = i__ + j * b_dim1;
		    i__4 = i__;
		    i__5 = i__ + 1 + j * b_dim1;
		    q__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, 
			    q__3.i = du[i__4].r * b[i__5].i + du[i__4].i * b[
			    i__5].r;
		    q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
		    c_div(&q__1, &q__2, &d__[i__]);
		    b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		}
	    }
	}
    } else if (*itrans == 1) {
	if (*nrhs <= 1) {
	    j = 1;
L70:
	    i__1 = j * b_dim1 + 1;
	    c_div(&q__1, &b[j * b_dim1 + 1], &d__[1]);
	    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
	    if (*n > 1) {
		i__1 = j * b_dim1 + 2;
		i__2 = j * b_dim1 + 2;
		i__3 = j * b_dim1 + 1;
		q__3.r = du[1].r * b[i__3].r - du[1].i * b[i__3].i, q__3.i = 
			du[1].r * b[i__3].i + du[1].i * b[i__3].r;
		q__2.r = b[i__2].r - q__3.r, q__2.i = b[i__2].i - q__3.i;
		c_div(&q__1, &q__2, &d__[2]);
		b[i__1].r = q__1.r, b[i__1].i = q__1.i;
	    }
	    i__1 = *n;
	    for (i__ = 3; i__ <= i__1; ++i__) {
		i__2 = i__ + j * b_dim1;
		i__3 = i__ + j * b_dim1;
		i__4 = i__ - 1;
		i__5 = i__ - 1 + j * b_dim1;
		q__3.r = du[i__4].r * b[i__5].r - du[i__4].i * b[i__5].i, 
			q__3.i = du[i__4].r * b[i__5].i + du[i__4].i * b[i__5]
			.r;
		q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
		c_div(&q__1, &q__2, &d__[i__]);
		b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	    }
	    for (i__ = *n - 1; i__ >= 1; --i__) {
		i__1 = i__ + j * b_dim1;
		i__2 = i__ + j * b_dim1;
		i__3 = i__;
		i__4 = i__ + 1 + j * b_dim1;
		q__2.r = dl[i__3].r * b[i__4].r - dl[i__3].i * b[i__4].i, 
			q__2.i = dl[i__3].r * b[i__4].i + dl[i__3].i * b[i__4]
			.r;
		q__1.r = b[i__2].r - q__2.r, q__1.i = b[i__2].i - q__2.i;
		b[i__1].r = q__1.r, b[i__1].i = q__1.i;
	    }
	    if (j < *nrhs) {
		++j;
		goto L70;
	    }
	} else {
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j * b_dim1 + 1;
		c_div(&q__1, &b[j * b_dim1 + 1], &d__[1]);
		b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		if (*n > 1) {
		    i__2 = j * b_dim1 + 2;
		    i__3 = j * b_dim1 + 2;
		    i__4 = j * b_dim1 + 1;
		    q__3.r = du[1].r * b[i__4].r - du[1].i * b[i__4].i, 
			    q__3.i = du[1].r * b[i__4].i + du[1].i * b[i__4]
			    .r;
		    q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
		    c_div(&q__1, &q__2, &d__[2]);
		    b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		}
		i__2 = *n;
		for (i__ = 3; i__ <= i__2; ++i__) {
		    i__3 = i__ + j * b_dim1;
		    i__4 = i__ + j * b_dim1;
		    i__5 = i__ - 1;
		    i__6 = i__ - 1 + j * b_dim1;
		    q__3.r = du[i__5].r * b[i__6].r - du[i__5].i * b[i__6].i, 
			    q__3.i = du[i__5].r * b[i__6].i + du[i__5].i * b[
			    i__6].r;
		    q__2.r = b[i__4].r - q__3.r, q__2.i = b[i__4].i - q__3.i;
		    c_div(&q__1, &q__2, &d__[i__]);
		    b[i__3].r = q__1.r, b[i__3].i = q__1.i;
		}
		for (i__ = *n - 1; i__ >= 1; --i__) {
		    i__2 = i__ + j * b_dim1;
		    i__3 = i__ + j * b_dim1;
		    i__4 = i__;
		    i__5 = i__ + 1 + j * b_dim1;
		    q__2.r = dl[i__4].r * b[i__5].r - dl[i__4].i * b[i__5].i, 
			    q__2.i = dl[i__4].r * b[i__5].i + dl[i__4].i * b[
			    i__5].r;
		    q__1.r = b[i__3].r - q__2.r, q__1.i = b[i__3].i - q__2.i;
		    b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		}
	    }
	}
    } else {
	if (*nrhs <= 1) {
	    j = 1;
L130:
	    i__1 = j * b_dim1 + 1;
	    r_cnjg(&q__2, &d__[1]);
	    c_div(&q__1, &b[j * b_dim1 + 1], &q__2);
	    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
	    if (*n > 1) {
		i__1 = j * b_dim1 + 2;
		i__2 = j * b_dim1 + 2;
		r_cnjg(&q__4, &du[1]);
		i__3 = j * b_dim1 + 1;
		q__3.r = q__4.r * b[i__3].r - q__4.i * b[i__3].i, q__3.i = 
			q__4.r * b[i__3].i + q__4.i * b[i__3].r;
		q__2.r = b[i__2].r - q__3.r, q__2.i = b[i__2].i - q__3.i;
		r_cnjg(&q__5, &d__[2]);
		c_div(&q__1, &q__2, &q__5);
		b[i__1].r = q__1.r, b[i__1].i = q__1.i;
	    }
	    i__1 = *n;
	    for (i__ = 3; i__ <= i__1; ++i__) {
		i__2 = i__ + j * b_dim1;
		i__3 = i__ + j * b_dim1;
		r_cnjg(&q__4, &du[i__ - 1]);
		i__4 = i__ - 1 + j * b_dim1;
		q__3.r = q__4.r * b[i__4].r - q__4.i * b[i__4].i, q__3.i = 
			q__4.r * b[i__4].i + q__4.i * b[i__4].r;
		q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
		r_cnjg(&q__5, &d__[i__]);
		c_div(&q__1, &q__2, &q__5);
		b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	    }
	    for (i__ = *n - 1; i__ >= 1; --i__) {
		i__1 = i__ + j * b_dim1;
		i__2 = i__ + j * b_dim1;
		r_cnjg(&q__3, &dl[i__]);
		i__3 = i__ + 1 + j * b_dim1;
		q__2.r = q__3.r * b[i__3].r - q__3.i * b[i__3].i, q__2.i = 
			q__3.r * b[i__3].i + q__3.i * b[i__3].r;
		q__1.r = b[i__2].r - q__2.r, q__1.i = b[i__2].i - q__2.i;
		b[i__1].r = q__1.r, b[i__1].i = q__1.i;
	    }
	    if (j < *nrhs) {
		++j;
		goto L130;
	    }
	} else {
	    i__1 = *nrhs;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = j * b_dim1 + 1;
		r_cnjg(&q__2, &d__[1]);
		c_div(&q__1, &b[j * b_dim1 + 1], &q__2);
		b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		if (*n > 1) {
		    i__2 = j * b_dim1 + 2;
		    i__3 = j * b_dim1 + 2;
		    r_cnjg(&q__4, &du[1]);
		    i__4 = j * b_dim1 + 1;
		    q__3.r = q__4.r * b[i__4].r - q__4.i * b[i__4].i, q__3.i =
			     q__4.r * b[i__4].i + q__4.i * b[i__4].r;
		    q__2.r = b[i__3].r - q__3.r, q__2.i = b[i__3].i - q__3.i;
		    r_cnjg(&q__5, &d__[2]);
		    c_div(&q__1, &q__2, &q__5);
		    b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		}
		i__2 = *n;
		for (i__ = 3; i__ <= i__2; ++i__) {
		    i__3 = i__ + j * b_dim1;
		    i__4 = i__ + j * b_dim1;
		    r_cnjg(&q__4, &du[i__ - 1]);
		    i__5 = i__ - 1 + j * b_dim1;
		    q__3.r = q__4.r * b[i__5].r - q__4.i * b[i__5].i, q__3.i =
			     q__4.r * b[i__5].i + q__4.i * b[i__5].r;
		    q__2.r = b[i__4].r - q__3.r, q__2.i = b[i__4].i - q__3.i;
		    r_cnjg(&q__5, &d__[i__]);
		    c_div(&q__1, &q__2, &q__5);
		    b[i__3].r = q__1.r, b[i__3].i = q__1.i;
		}
		for (i__ = *n - 1; i__ >= 1; --i__) {
		    i__2 = i__ + j * b_dim1;
		    i__3 = i__ + j * b_dim1;
		    r_cnjg(&q__3, &dl[i__]);
		    i__4 = i__ + 1 + j * b_dim1;
		    q__2.r = q__3.r * b[i__4].r - q__3.i * b[i__4].i, q__2.i =
			     q__3.r * b[i__4].i + q__3.i * b[i__4].r;
		    q__1.r = b[i__3].r - q__2.r, q__1.i = b[i__3].i - q__2.i;
		    b[i__2].r = q__1.r, b[i__2].i = q__1.i;
		}
	    }
	}
    }
    return 0;
} /* cgtts2np_ */
