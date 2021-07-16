#line 1 "./finite_difference.f"
/*  -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

#line 1 "./finite_difference.f"
/* Common Block Declarations */

struct {
    doublereal roundoff, small, huge__, undefind;
} machine_;

#define machine_1 machine_

/* Subroutine */ int initialize_(fi, la, ncells, matrix, rhs, solution)
integer *fi, *la, *ncells;
doublereal *matrix, *rhs, *solution;
{
    /* System generated locals */
    integer matrix_dim1, matrix_offset, rhs_offset, solution_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sin();

    /* Local variables */
    static integer i__;
    static doublereal x, pi, dx, left;
    extern /* Subroutine */ int abort_();
    static doublereal right, dxinv;

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     bvp: */
/*       - d/dx( p du/dx ) + r u = f, 0 < x < 1 */
/*       u(0)=left, u(1)=right */
#line 23 "./finite_difference.f"
    /* Parameter adjustments */
#line 23 "./finite_difference.f"
    solution_offset = *fi;
#line 23 "./finite_difference.f"
    solution -= solution_offset;
#line 23 "./finite_difference.f"
    rhs_offset = *fi;
#line 23 "./finite_difference.f"
    rhs -= rhs_offset;
#line 23 "./finite_difference.f"
    matrix_dim1 = *la - *fi + 1;
#line 23 "./finite_difference.f"
    matrix_offset = *fi + matrix_dim1 * -1;
#line 23 "./finite_difference.f"
    matrix -= matrix_offset;
#line 23 "./finite_difference.f"

#line 23 "./finite_difference.f"
    /* Function Body */
#line 23 "./finite_difference.f"
    left = (float)0.;
#line 24 "./finite_difference.f"
    right = (float)0.;
#line 25 "./finite_difference.f"
    pi = 3.14159265358979323846;
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 27 "./finite_difference.f"
    if (*fi > 1 || *la < *ncells - 1) {
#line 27 "./finite_difference.f"
	abort_();
#line 27 "./finite_difference.f"
    }
#line 29 "./finite_difference.f"
    dxinv = (doublereal) (*ncells);
#line 30 "./finite_difference.f"
    dx = 1. / dxinv;
#line 32 "./finite_difference.f"
    x = dx * .5;
#line 33 "./finite_difference.f"
    matrix[-matrix_dim1 + 1] = -(1. * dxinv);
#line 34 "./finite_difference.f"
    i__1 = *ncells - 1;
#line 34 "./finite_difference.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 35 "./finite_difference.f"
	x += dx;
#line 36 "./finite_difference.f"
	matrix[i__ + (-matrix_dim1)] = -(1. * dxinv);
#line 37 "./finite_difference.f"
	matrix[i__ - 1 + matrix_dim1] = matrix[i__ + (-matrix_dim1)];
#line 38 "./finite_difference.f"
    }
#line 39 "./finite_difference.f"
    x += dx;
#line 40 "./finite_difference.f"
    matrix[*ncells - 1 + matrix_dim1] = -(1. * dxinv);
#line 42 "./finite_difference.f"
    x = 0.;
#line 43 "./finite_difference.f"
    i__1 = *ncells - 1;
#line 43 "./finite_difference.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 44 "./finite_difference.f"
	x += dx;
#line 45 "./finite_difference.f"
	matrix[i__] = -matrix[i__ + (-matrix_dim1)] - matrix[i__ + 
		matrix_dim1] + 0. * dx;
/* Computing 2nd power */
#line 46 "./finite_difference.f"
	d__1 = pi;
#line 46 "./finite_difference.f"
	rhs[i__] = sin(pi * x) * (d__1 * d__1) * dx;
#line 47 "./finite_difference.f"
	solution[i__] = sin(pi * x);
#line 48 "./finite_difference.f"
    }
#line 50 "./finite_difference.f"
    rhs[1] -= matrix[-matrix_dim1 + 1] * left;
#line 51 "./finite_difference.f"
    rhs[*ncells - 1] -= matrix[*ncells - 1 + matrix_dim1] * right;
/*     call print_loc(rhs(0)) */
/*     do i=1,ncells-1 */
/*       print *, "matrix[",i,"] = ",matrix(i,-1),matrix(i,0),matrix(i,1) */
/*     enddo */
/*     do i=1,ncells-1 */
/*       print *, "rhs[",i,"] = ",rhs(i) */
/*     enddo */
/*     do i=1,ncells-1 */
/*       print *, "solution[",i,"] = ",solution(i) */
/*     enddo */
/*     print *, "leaving initialize" */
/*     call flush(6) */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 66 "./finite_difference.f"
    return 0;
} /* initialize_ */

/* *********************************************************************** */
/* Subroutine */ int solve_(fi, la, ncells, matrix, rhs_in_soln_out__)
integer *fi, *la, *ncells;
doublereal *matrix, *rhs_in_soln_out__;
{
    /* System generated locals */
    integer matrix_dim1, matrix_offset, rhs_in_soln_out_offset, i__1;

    /* Local variables */
    static integer i__;
    extern /* Subroutine */ int abort_();

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 79 "./finite_difference.f"
    /* Parameter adjustments */
#line 79 "./finite_difference.f"
    rhs_in_soln_out_offset = *fi;
#line 79 "./finite_difference.f"
    rhs_in_soln_out__ -= rhs_in_soln_out_offset;
#line 79 "./finite_difference.f"
    matrix_dim1 = *la - *fi + 1;
#line 79 "./finite_difference.f"
    matrix_offset = *fi + matrix_dim1 * -1;
#line 79 "./finite_difference.f"
    matrix -= matrix_offset;
#line 79 "./finite_difference.f"

#line 79 "./finite_difference.f"
    /* Function Body */
#line 79 "./finite_difference.f"
    if (*fi > 1 || *la < *ncells - 1) {
#line 79 "./finite_difference.f"
	abort_();
#line 79 "./finite_difference.f"
    }
/*     print *, "entering solve" */
/*     do i=1,ncells-1 */
/*       print *, "matrix[",i,"] = ",matrix(i,-1),matrix(i,0),matrix(i,1) */
/*     enddo */
/*     do i=1,ncells-1 */
/*       print *, "rhs[",i,"] = ",rhs_in_soln_out(i) */
/*     enddo */
/*     call flush(6) */
#line 89 "./finite_difference.f"
    i__1 = *ncells - 1;
#line 89 "./finite_difference.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 90 "./finite_difference.f"
	matrix[i__ + (-matrix_dim1)] /= matrix[i__ - 1];
#line 91 "./finite_difference.f"
	matrix[i__] -= matrix[i__ + (-matrix_dim1)] * matrix[i__ - 1 + 
		matrix_dim1];
#line 92 "./finite_difference.f"
	rhs_in_soln_out__[i__] -= matrix[i__ + (-matrix_dim1)] * 
		rhs_in_soln_out__[i__ - 1];
#line 94 "./finite_difference.f"
    }
#line 95 "./finite_difference.f"
    rhs_in_soln_out__[*ncells - 1] /= matrix[*ncells - 1];
#line 97 "./finite_difference.f"
    for (i__ = *ncells - 2; i__ >= 1; --i__) {
#line 98 "./finite_difference.f"
	rhs_in_soln_out__[i__] = (rhs_in_soln_out__[i__] - matrix[i__ + 
		matrix_dim1] * rhs_in_soln_out__[i__ + 1]) / matrix[i__];
#line 100 "./finite_difference.f"
    }
/*     print *, "leaving solve" */
/*     do i=1,ncells-1 */
/*       print *, "soln[",i,"] = ",rhs_in_soln_out(i) */
/*     enddo */
/*     call flush(6) */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 108 "./finite_difference.f"
    return 0;
} /* solve_ */

