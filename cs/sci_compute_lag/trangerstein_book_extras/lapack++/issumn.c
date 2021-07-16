/* issumn.f -- translated by f2c (version 20090411).
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

integer issumn_(integer *n, real *a, integer *inca)
{
    /* System generated locals */
    integer a_dim1, a_offset, ret_val, i__1;

    /* Local variables */
    static integer i__;

    /* Parameter adjustments */
    a_dim1 = *inca;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    ret_val = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (a[i__ * a_dim1 + 1] <= 0.) {
	    ++ret_val;
	}
    }
    return ret_val;
} /* issumn_ */

