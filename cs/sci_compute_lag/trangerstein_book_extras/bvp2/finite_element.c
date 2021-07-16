#line 1 "./finite_element.f"
/*  -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

#line 1 "./finite_element.f"
/* Common Block Declarations */

struct {
    doublereal roundoff, small, huge__, undefind, pi;
} machine_;

#define machine_1 machine_

/* Table of constant values */

static integer c__1 = 1;
static real c_b24 = (float)10.;

/* Subroutine */ int canonical_(continuity, order, basis, basis_deriv__, 
	gauss, weight)
integer *continuity, *order;
doublereal *basis, *basis_deriv__, *gauss, *weight;
{
    /* System generated locals */
    integer basis_dim1, basis_offset, basis_deriv_dim1, basis_deriv_offset, 
	    i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer j;
    static doublereal g2, g4, omg, omg2, om2g, om3g, tm3g;

/*     integer i */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 13 "./finite_element.f"
    /* Parameter adjustments */
#line 13 "./finite_element.f"
    --weight;
#line 13 "./finite_element.f"
    --gauss;
#line 13 "./finite_element.f"
    basis_deriv_dim1 = *order + 1;
#line 13 "./finite_element.f"
    basis_deriv_offset = 1 + basis_deriv_dim1 * 1;
#line 13 "./finite_element.f"
    basis_deriv__ -= basis_deriv_offset;
#line 13 "./finite_element.f"
    basis_dim1 = *order + 1;
#line 13 "./finite_element.f"
    basis_offset = 1 + basis_dim1 * 1;
#line 13 "./finite_element.f"
    basis -= basis_offset;
#line 13 "./finite_element.f"

#line 13 "./finite_element.f"
    /* Function Body */
#line 13 "./finite_element.f"
    if (*order == 1) {
/*       continuous piecewise linear */
#line 15 "./finite_element.f"
	gauss[1] = .5;
#line 16 "./finite_element.f"
	weight[1] = 1.;
#line 18 "./finite_element.f"
	basis[basis_dim1 + 1] = 1. - gauss[1];
#line 19 "./finite_element.f"
	basis[basis_dim1 + 2] = gauss[1];
#line 21 "./finite_element.f"
	basis_deriv__[basis_deriv_dim1 + 1] = -1.;
#line 22 "./finite_element.f"
	basis_deriv__[basis_deriv_dim1 + 2] = 1.;
#line 23 "./finite_element.f"
    } else if (*order == 2) {
/*       continuous piecewise quadratic */
#line 25 "./finite_element.f"
	gauss[1] = (1. - sqrt(.33333333333333331)) * .5;
#line 26 "./finite_element.f"
	gauss[2] = (sqrt(.33333333333333331) + 1.) * .5;
#line 27 "./finite_element.f"
	weight[1] = .5;
#line 28 "./finite_element.f"
	weight[2] = .5;
#line 30 "./finite_element.f"
	i__1 = *order;
#line 30 "./finite_element.f"
	for (j = 1; j <= i__1; ++j) {
#line 31 "./finite_element.f"
	    omg = 1. - gauss[j];
#line 32 "./finite_element.f"
	    om2g = 1. - gauss[j] * 2.;
#line 33 "./finite_element.f"
	    g4 = gauss[j] * 4.;
#line 35 "./finite_element.f"
	    basis[j * basis_dim1 + 1] = om2g * omg;
#line 36 "./finite_element.f"
	    basis[j * basis_dim1 + 2] = g4 * omg;
#line 37 "./finite_element.f"
	    basis[j * basis_dim1 + 3] = -(gauss[j] * om2g);
#line 39 "./finite_element.f"
	    basis_deriv__[j * basis_deriv_dim1 + 1] = g4 - 3.;
#line 40 "./finite_element.f"
	    basis_deriv__[j * basis_deriv_dim1 + 2] = om2g * 4.;
#line 41 "./finite_element.f"
	    basis_deriv__[j * basis_deriv_dim1 + 3] = g4 - 1.;
#line 42 "./finite_element.f"
	}
#line 43 "./finite_element.f"
    } else if (*order == 3) {
#line 44 "./finite_element.f"
	gauss[1] = (1. - sqrt(.6)) * .5;
#line 45 "./finite_element.f"
	gauss[2] = .5;
#line 46 "./finite_element.f"
	gauss[3] = (sqrt(.6) + 1.) * .5;
#line 47 "./finite_element.f"
	weight[1] = .27777777777777779;
#line 48 "./finite_element.f"
	weight[2] = .44444444444444442;
#line 49 "./finite_element.f"
	weight[3] = weight[1];
#line 50 "./finite_element.f"
	if (*continuity == 0) {
/*         continuous piecewise cubic */
#line 52 "./finite_element.f"
	    i__1 = *order;
#line 52 "./finite_element.f"
	    for (j = 1; j <= i__1; ++j) {
#line 53 "./finite_element.f"
		omg = 1. - gauss[j];
#line 54 "./finite_element.f"
		om3g = 1. - gauss[j] * 3.;
#line 55 "./finite_element.f"
		tm3g = om3g + 1.;
#line 57 "./finite_element.f"
		basis[j * basis_dim1 + 1] = om3g * .5 * tm3g * omg;
#line 58 "./finite_element.f"
		basis[j * basis_dim1 + 2] = gauss[j] * 4.5 * tm3g * omg;
#line 59 "./finite_element.f"
		basis[j * basis_dim1 + 3] = -(gauss[j] * 4.5 * om3g * omg);
#line 60 "./finite_element.f"
		basis[j * basis_dim1 + 4] = gauss[j] * .5 * om3g * tm3g;
#line 62 "./finite_element.f"
		basis_deriv__[j * basis_deriv_dim1 + 1] = gauss[j] * (18. - 
			gauss[j] * 13.5) - 5.5;
#line 63 "./finite_element.f"
		basis_deriv__[j * basis_deriv_dim1 + 2] = gauss[j] * (gauss[j]
			 * 40.5 - 45.) + 9.;
#line 64 "./finite_element.f"
		basis_deriv__[j * basis_deriv_dim1 + 3] = gauss[j] * (36. - 
			gauss[j] * 40.5) - 4.5;
#line 65 "./finite_element.f"
		basis_deriv__[j * basis_deriv_dim1 + 4] = gauss[j] * (gauss[j]
			 * 13.5 - 9.) + 1.;
#line 66 "./finite_element.f"
	    }
#line 67 "./finite_element.f"
	} else if (*continuity == 1) {
/*         Hermite piecewise cubic */
#line 69 "./finite_element.f"
	    i__1 = *order;
#line 69 "./finite_element.f"
	    for (j = 1; j <= i__1; ++j) {
#line 70 "./finite_element.f"
		omg = 1. - gauss[j];
/* Computing 2nd power */
#line 71 "./finite_element.f"
		d__1 = omg;
#line 71 "./finite_element.f"
		omg2 = d__1 * d__1;
/* Computing 2nd power */
#line 72 "./finite_element.f"
		d__1 = gauss[j];
#line 72 "./finite_element.f"
		g2 = d__1 * d__1;
/*           interpolate value at 0, slope at 0, value at 1, slope at 1 */
#line 75 "./finite_element.f"
		basis[j * basis_dim1 + 1] = (gauss[j] * 2. + 1.) * omg2;
#line 76 "./finite_element.f"
		basis[j * basis_dim1 + 2] = gauss[j] * omg2;
#line 77 "./finite_element.f"
		basis[j * basis_dim1 + 3] = (3. - gauss[j] * 2.) * g2;
#line 78 "./finite_element.f"
		basis[j * basis_dim1 + 4] = -(omg * g2);
#line 80 "./finite_element.f"
		basis_deriv__[j * basis_deriv_dim1 + 1] = -(gauss[j] * 6. * 
			omg);
#line 81 "./finite_element.f"
		basis_deriv__[j * basis_deriv_dim1 + 2] = omg * (1. - gauss[j]
			 * 3.);
#line 82 "./finite_element.f"
		basis_deriv__[j * basis_deriv_dim1 + 3] = gauss[j] * 6. * omg;
#line 83 "./finite_element.f"
		basis_deriv__[j * basis_deriv_dim1 + 4] = -(gauss[j] * (2. - 
			gauss[j] * 3.));
#line 84 "./finite_element.f"
	    }
#line 85 "./finite_element.f"
	}
/*     else if (order.eq.4) then */
/*       g1=sqrt(525.d0-70.d0*sqrt(30.d0))/35.d0 */
/*       g2=sqrt(525.d0+70.d0*sqrt(30.d0))/35.d0 */
/*       gauss(1)=0.5d0*(1.d0-g2) */
/*       gauss(1)=0.5d0*(1.d0-g1) */
/*       gauss(1)=0.5d0*(1.d0+g1) */
/*       gauss(1)=0.5d0*(1.d0+g2) */
/*       weight(1)=(18.d0-sqrt(30.d0))/72.d0 */
/*       weight(2)=(18.d0+sqrt(30.d0))/72.d0 */
/*       weight(3)=weight(2) */
/*       weight(4)=weight(1) */
/*       if (continuity.eq.0) then */
/*       else */
/*         do j=1,order */
/*           basis(1,j)=(1.d0-2.d0*gauss(j))*(1.d0+4.d0*gauss(j)) */
/*    &                *(1.d0-gauss(j))**2 */
/*           basis(2,j)=gauss(j)*(2.d0*gauss(j)-1.d0)*(1.d0-gauss(j))**2 */
/*           basis(3,j)=16.d0*(gauss(j)*(1.d0-gauss(j)))**2 */
/*           basis(4,j)=(-1.d0+2.d0*gauss(j))*(5.d0-4.d0*gauss(j)) */
/*    &                *gauss(j)**2 */
/*           basis(5,j)=(1.d0-2.d0*gauss(j))*(1.d0-gauss(j))*gauss(j)**2 */
/*           basis_deriv(1,j)=2.d0*gauss(j)*(1.d0-gauss(j)) */
/*    &                      *(11.d0-16.d0*gauss(j)) */
/*           basis_deriv(2,j)=(1.d0-gauss(j)) */
/*    &                      *(1.d0+gauss(j)*(-7.d0+gauss(j)*8.d0)) */
/*           basis_deriv(3,j)=32.d0*gauss(j)*(1.d0-gauss(j)) */
/*    &                      *(1.d0-2.d0*gauss(j)) */
/*           basis_deriv(4,j)=2.d0*gauss(j)*(1.d0-gauss(j)) */
/*    &                      *(16.d0*gauss(j)-5.d0) */
/*           basis_deriv(5,j)=gauss(j) */
/*    &                      *(2.d0+gauss(j)*(-9.d0+gauss(j)*8.d0)) */
/*         enddo */
/*       endif */
#line 119 "./finite_element.f"
    }
/*     do i=1,order */
/*       print *, "gauss[",i,"] = ",gauss(i) */
/*     enddo */
/*     do i=1,order */
/*       print *, "weight[",i,"] = ",weight(i) */
/*     enddo */
/*     do i=1,order+1 */
/*       print *, "basis[",i,"] = ",(basis(i,j),j=1,order) */
/*     enddo */
/*     do i=1,order+1 */
/*       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=1,order) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 133 "./finite_element.f"
    return 0;
} /* canonical_ */

/* *********************************************************************** */
doublereal approximation_(continuity, nelements, nnodes, order, x, 
	element_to_node__, mesh, soln)
integer *continuity, *nelements, *nnodes, *order;
doublereal *x;
integer *element_to_node__;
doublereal *mesh, *soln;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static doublereal dx, xi, omg, om3g, tm3g;
    static integer element;

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 150 "./finite_element.f"
    /* Parameter adjustments */
#line 150 "./finite_element.f"
    --mesh;
#line 150 "./finite_element.f"
    --soln;
#line 150 "./finite_element.f"
    element_to_node_dim1 = *nelements;
#line 150 "./finite_element.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 1;
#line 150 "./finite_element.f"
    element_to_node__ -= element_to_node_offset;
#line 150 "./finite_element.f"

#line 150 "./finite_element.f"
    /* Function Body */
#line 150 "./finite_element.f"
    dx = (float)1. / (doublereal) (*nelements);
#line 151 "./finite_element.f"
    element = (integer) (*x / dx);
#line 152 "./finite_element.f"
    if (element < *nelements) {
#line 152 "./finite_element.f"
	++element;
#line 152 "./finite_element.f"
    }
#line 153 "./finite_element.f"
    xi = (*x - mesh[element]) / dx;
#line 154 "./finite_element.f"
    if (*order == 1) {
/*       print *, "x = ",x,mesh(element),dx */
/*       print *, "element = ",element */
/*       print *, "xi = ",xi */
/*       print *, "nodes = ",element_to_node(element,1), */
/*    &    element_to_node(element,2) */
#line 160 "./finite_element.f"
	ret_val = (1. - xi) * soln[element_to_node__[element + 
		element_to_node_dim1]] + xi * soln[element_to_node__[element 
		+ (element_to_node_dim1 << 1)]];
#line 162 "./finite_element.f"
    } else if (*order == 2) {
#line 163 "./finite_element.f"
	omg = 1. - xi;
#line 164 "./finite_element.f"
	ret_val = (1. - xi * 2.) * (omg * soln[element_to_node__[element + 
		element_to_node_dim1]] - xi * soln[element_to_node__[element 
		+ element_to_node_dim1 * 3]]) + xi * 4. * omg * soln[
		element_to_node__[element + (element_to_node_dim1 << 1)]];
#line 168 "./finite_element.f"
    } else if (*order == 3) {
#line 169 "./finite_element.f"
	if (*continuity == 0) {
#line 170 "./finite_element.f"
	    omg = 1. - xi;
#line 171 "./finite_element.f"
	    om3g = 1. - xi * 3.;
#line 172 "./finite_element.f"
	    tm3g = om3g + 1.;
#line 173 "./finite_element.f"
	    ret_val = om3g * .5 * tm3g * (omg * soln[element_to_node__[
		    element + element_to_node_dim1]] + xi * soln[
		    element_to_node__[element + (element_to_node_dim1 << 2)]])
		     + xi * 4.5 * omg * (tm3g * soln[element_to_node__[
		    element + (element_to_node_dim1 << 1)]] - om3g * soln[
		    element_to_node__[element + element_to_node_dim1 * 3]]);
#line 178 "./finite_element.f"
	} else {
#line 179 "./finite_element.f"
	    omg = 1. - xi;
/* Computing 2nd power */
#line 180 "./finite_element.f"
	    d__1 = omg;
/* Computing 2nd power */
#line 180 "./finite_element.f"
	    d__2 = xi;
#line 180 "./finite_element.f"
	    ret_val = ((xi * 2. + 1.) * soln[element_to_node__[element + 
		    element_to_node_dim1]] + xi * soln[element_to_node__[
		    element + (element_to_node_dim1 << 1)]]) * (d__1 * d__1) 
		    + ((3. - xi * 2.) * soln[element_to_node__[element + 
		    element_to_node_dim1 * 3]] - omg * soln[element_to_node__[
		    element + (element_to_node_dim1 << 2)]]) * (d__2 * d__2);
#line 185 "./finite_element.f"
	}
#line 186 "./finite_element.f"
    }
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 188 "./finite_element.f"
    return ret_val;
} /* approximation_ */

/* *********************************************************************** */
/* Subroutine */ int grid_(continuity, nelements, nnodes, order, 
	element_to_node__, mesh)
integer *continuity, *nelements, *nnodes, *order, *element_to_node__;
doublereal *mesh;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, i__1;

    /* Local variables */
    static integer i__, e2, e3;
    static doublereal dx;
    extern /* Subroutine */ int abort_();
    static integer element;

/*     integer j */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 203 "./finite_element.f"
    /* Parameter adjustments */
#line 203 "./finite_element.f"
    --mesh;
#line 203 "./finite_element.f"
    element_to_node_dim1 = *nelements;
#line 203 "./finite_element.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 1;
#line 203 "./finite_element.f"
    element_to_node__ -= element_to_node_offset;
#line 203 "./finite_element.f"

#line 203 "./finite_element.f"
    /* Function Body */
#line 203 "./finite_element.f"
    if (*nnodes != (*order - *continuity) * *nelements + *continuity + 1) {
#line 204 "./finite_element.f"
	abort_();
#line 205 "./finite_element.f"
    }
#line 207 "./finite_element.f"
    dx = 1. / (doublereal) (*nelements);
#line 208 "./finite_element.f"
    i__1 = *nelements + 1;
#line 208 "./finite_element.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "./finite_element.f"
	mesh[i__] = (doublereal) (i__ - 1) * dx;
#line 210 "./finite_element.f"
    }
#line 212 "./finite_element.f"
    if (*continuity >= *order) {
#line 212 "./finite_element.f"
	abort_();
#line 212 "./finite_element.f"
    }
#line 213 "./finite_element.f"
    if (*order == 1) {
#line 214 "./finite_element.f"
	if (*continuity != 0) {
#line 214 "./finite_element.f"
	    abort_();
#line 214 "./finite_element.f"
	}
#line 215 "./finite_element.f"
	i__1 = *nelements;
#line 215 "./finite_element.f"
	for (element = 1; element <= i__1; ++element) {
#line 216 "./finite_element.f"
	    element_to_node__[element + element_to_node_dim1] = element;
#line 217 "./finite_element.f"
	    element_to_node__[element + (element_to_node_dim1 << 1)] = 
		    element + 1;
#line 218 "./finite_element.f"
	}
#line 219 "./finite_element.f"
    } else if (*order == 2) {
#line 220 "./finite_element.f"
	if (*continuity != 0) {
#line 220 "./finite_element.f"
	    abort_();
#line 220 "./finite_element.f"
	}
#line 221 "./finite_element.f"
	i__1 = *nelements;
#line 221 "./finite_element.f"
	for (element = 1; element <= i__1; ++element) {
#line 222 "./finite_element.f"
	    e2 = element << 1;
#line 223 "./finite_element.f"
	    element_to_node__[element + element_to_node_dim1] = e2 - 1;
#line 224 "./finite_element.f"
	    element_to_node__[element + (element_to_node_dim1 << 1)] = e2;
#line 225 "./finite_element.f"
	    element_to_node__[element + element_to_node_dim1 * 3] = e2 + 1;
#line 226 "./finite_element.f"
	}
#line 227 "./finite_element.f"
	e2 = *nelements << 1;
#line 228 "./finite_element.f"
    } else if (*order == 3) {
#line 229 "./finite_element.f"
	if (*continuity == 0) {
#line 230 "./finite_element.f"
	    i__1 = *nelements;
#line 230 "./finite_element.f"
	    for (element = 1; element <= i__1; ++element) {
#line 231 "./finite_element.f"
		e3 = element * 3;
#line 232 "./finite_element.f"
		element_to_node__[element + element_to_node_dim1] = e3 - 2;
#line 233 "./finite_element.f"
		element_to_node__[element + (element_to_node_dim1 << 1)] = e3 
			- 1;
#line 234 "./finite_element.f"
		element_to_node__[element + element_to_node_dim1 * 3] = e3;
#line 235 "./finite_element.f"
		element_to_node__[element + (element_to_node_dim1 << 2)] = e3 
			+ 1;
#line 236 "./finite_element.f"
	    }
#line 237 "./finite_element.f"
	    e3 = *nelements * 3;
#line 238 "./finite_element.f"
	} else {
#line 239 "./finite_element.f"
	    if (*continuity != 1) {
#line 239 "./finite_element.f"
		abort_();
#line 239 "./finite_element.f"
	    }
#line 240 "./finite_element.f"
	    i__1 = *nelements;
#line 240 "./finite_element.f"
	    for (element = 1; element <= i__1; ++element) {
#line 241 "./finite_element.f"
		e2 = element << 1;
#line 242 "./finite_element.f"
		element_to_node__[element + element_to_node_dim1] = e2 - 1;
#line 243 "./finite_element.f"
		element_to_node__[element + (element_to_node_dim1 << 1)] = e2;
#line 244 "./finite_element.f"
		element_to_node__[element + element_to_node_dim1 * 3] = e2 + 
			1;
#line 245 "./finite_element.f"
		element_to_node__[element + (element_to_node_dim1 << 2)] = e2 
			+ 2;
#line 246 "./finite_element.f"
	    }
#line 247 "./finite_element.f"
	    e2 = *nelements << 1;
#line 248 "./finite_element.f"
	}
#line 249 "./finite_element.f"
    }
/*     do i=1,nelements+1 */
/*       print *, "mesh[",i,"] = ",mesh(i) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "element_to_node[",i,"] = ", */
/*    &    (element_to_node(i,j),j=1,order+1) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 259 "./finite_element.f"
    return 0;
} /* grid_ */

/* *********************************************************************** */
doublereal solution_(x)
doublereal *x;
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sin();

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 271 "./finite_element.f"
    ret_val = sin(machine_1.pi * *x);
/*     print *, "x = ",x */
/*     print *, "pi = ",pi */
/*     print *, "solution = ",solution */
#line 275 "./finite_element.f"
    return ret_val;
} /* solution_ */

/* *********************************************************************** */
/* Subroutine */ int mult_(nelements, nnodes, order, basis, basis_deriv__, 
	dirichlet, element_to_node__, pgauss, rgauss, x, ax)
integer *nelements, *nnodes, *order;
doublereal *basis, *basis_deriv__;
logical1 *dirichlet;
integer *element_to_node__;
doublereal *pgauss, *rgauss, *x, *ax;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, basis_dim1, 
	    basis_offset, basis_deriv_dim1, basis_deriv_offset, pgauss_dim1, 
	    pgauss_offset, rgauss_dim1, rgauss_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j;
    static doublereal bp, br;
    static integer node_i__, node_j__, ngauss, element;

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     do i=1,order+1 */
/*       print *, "basis[",i,"] = ",(basis(i,j),j=1,order) */
/*     enddo */
/*     do i=1,order+1 */
/*       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=1,order) */
/*     enddo */
/*     do i=1,nnodes */
/*       print *, "dirichlet[",i,"] = ",dirichlet(i) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "element_to_node[",i,"] = ", */
/*    &    (element_to_node(i,j),j=1,order+1) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "pgauss[",i,"] = ",(pgauss(i,j),j=1,order) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "rgauss[",i,"] = ",(rgauss(i,j),j=1,order) */
/*     enddo */
/*     do i=1,nnodes */
/*       print *, "x[",i,"] = ",x(i) */
/*     enddo */
#line 318 "./finite_element.f"
    /* Parameter adjustments */
#line 318 "./finite_element.f"
    --ax;
#line 318 "./finite_element.f"
    --x;
#line 318 "./finite_element.f"
    --dirichlet;
#line 318 "./finite_element.f"
    rgauss_dim1 = *nelements;
#line 318 "./finite_element.f"
    rgauss_offset = 1 + rgauss_dim1 * 1;
#line 318 "./finite_element.f"
    rgauss -= rgauss_offset;
#line 318 "./finite_element.f"
    pgauss_dim1 = *nelements;
#line 318 "./finite_element.f"
    pgauss_offset = 1 + pgauss_dim1 * 1;
#line 318 "./finite_element.f"
    pgauss -= pgauss_offset;
#line 318 "./finite_element.f"
    element_to_node_dim1 = *nelements;
#line 318 "./finite_element.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 1;
#line 318 "./finite_element.f"
    element_to_node__ -= element_to_node_offset;
#line 318 "./finite_element.f"
    basis_deriv_dim1 = *order + 1;
#line 318 "./finite_element.f"
    basis_deriv_offset = 1 + basis_deriv_dim1 * 1;
#line 318 "./finite_element.f"
    basis_deriv__ -= basis_deriv_offset;
#line 318 "./finite_element.f"
    basis_dim1 = *order + 1;
#line 318 "./finite_element.f"
    basis_offset = 1 + basis_dim1 * 1;
#line 318 "./finite_element.f"
    basis -= basis_offset;
#line 318 "./finite_element.f"

#line 318 "./finite_element.f"
    /* Function Body */
#line 318 "./finite_element.f"
    i__1 = *nnodes;
#line 318 "./finite_element.f"
    for (node_i__ = 1; node_i__ <= i__1; ++node_i__) {
#line 319 "./finite_element.f"
	ax[node_i__] = 0.;
#line 320 "./finite_element.f"
    }
#line 321 "./finite_element.f"
    i__1 = *nelements;
#line 321 "./finite_element.f"
    for (element = 1; element <= i__1; ++element) {
#line 322 "./finite_element.f"
	i__2 = *order;
#line 322 "./finite_element.f"
	for (ngauss = 1; ngauss <= i__2; ++ngauss) {
#line 323 "./finite_element.f"
	    i__3 = *order + 1;
#line 323 "./finite_element.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 324 "./finite_element.f"
		node_i__ = element_to_node__[element + i__ * 
			element_to_node_dim1];
#line 325 "./finite_element.f"
		bp = basis_deriv__[i__ + ngauss * basis_deriv_dim1] * pgauss[
			element + ngauss * pgauss_dim1];
#line 326 "./finite_element.f"
		br = basis[i__ + ngauss * basis_dim1] * rgauss[element + 
			ngauss * rgauss_dim1];
#line 327 "./finite_element.f"
		i__4 = *order + 1;
#line 327 "./finite_element.f"
		for (j = 1; j <= i__4; ++j) {
#line 328 "./finite_element.f"
		    node_j__ = element_to_node__[element + j * 
			    element_to_node_dim1];
#line 329 "./finite_element.f"
		    ax[node_i__] += (bp * basis_deriv__[j + ngauss * 
			    basis_deriv_dim1] + br * basis[j + ngauss * 
			    basis_dim1]) * x[node_j__];
#line 331 "./finite_element.f"
		}
#line 332 "./finite_element.f"
	    }
#line 333 "./finite_element.f"
	}
#line 334 "./finite_element.f"
    }
#line 335 "./finite_element.f"
    i__1 = *nnodes;
#line 335 "./finite_element.f"
    for (node_i__ = 1; node_i__ <= i__1; ++node_i__) {
#line 336 "./finite_element.f"
	if (dirichlet[node_i__]) {
#line 336 "./finite_element.f"
	    ax[node_i__] = 0.;
#line 336 "./finite_element.f"
	}
#line 337 "./finite_element.f"
    }
/*     do i=1,nnodes */
/*       print *,"x,Ax[",i,"] = ",x(i),Ax(i) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 342 "./finite_element.f"
    return 0;
} /* mult_ */

/* *********************************************************************** */
/* Subroutine */ int initialize_(continuity, nelements, nnodes, order, basis, 
	basis_deriv__, element_to_node__, gauss, mesh, weight, dirichlet, 
	pgauss, rgauss, residual, soln)
integer *continuity, *nelements, *nnodes, *order;
doublereal *basis, *basis_deriv__;
integer *element_to_node__;
doublereal *gauss, *mesh, *weight;
logical1 *dirichlet;
doublereal *pgauss, *rgauss, *residual, *soln;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, basis_dim1, 
	    basis_offset, basis_deriv_dim1, basis_deriv_offset, pgauss_dim1, 
	    pgauss_offset, rgauss_dim1, rgauss_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double sin();

    /* Local variables */
    static integer i__, j;
    static doublereal fg, dx, xg, left;
    extern /* Subroutine */ int mult_(), abort_();
    static doublereal right;
    static integer node_i__, ngauss, element;

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     bvp: */
/*       - d/dx( p du/dx ) + r u = f, 0 < x < 1 */
/*       u(0)=left, u(1)=right */
#line 378 "./finite_element.f"
    /* Parameter adjustments */
#line 378 "./finite_element.f"
    --mesh;
#line 378 "./finite_element.f"
    --soln;
#line 378 "./finite_element.f"
    --residual;
#line 378 "./finite_element.f"
    --dirichlet;
#line 378 "./finite_element.f"
    rgauss_dim1 = *nelements;
#line 378 "./finite_element.f"
    rgauss_offset = 1 + rgauss_dim1 * 1;
#line 378 "./finite_element.f"
    rgauss -= rgauss_offset;
#line 378 "./finite_element.f"
    pgauss_dim1 = *nelements;
#line 378 "./finite_element.f"
    pgauss_offset = 1 + pgauss_dim1 * 1;
#line 378 "./finite_element.f"
    pgauss -= pgauss_offset;
#line 378 "./finite_element.f"
    --weight;
#line 378 "./finite_element.f"
    --gauss;
#line 378 "./finite_element.f"
    element_to_node_dim1 = *nelements;
#line 378 "./finite_element.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 1;
#line 378 "./finite_element.f"
    element_to_node__ -= element_to_node_offset;
#line 378 "./finite_element.f"
    basis_deriv_dim1 = *order + 1;
#line 378 "./finite_element.f"
    basis_deriv_offset = 1 + basis_deriv_dim1 * 1;
#line 378 "./finite_element.f"
    basis_deriv__ -= basis_deriv_offset;
#line 378 "./finite_element.f"
    basis_dim1 = *order + 1;
#line 378 "./finite_element.f"
    basis_offset = 1 + basis_dim1 * 1;
#line 378 "./finite_element.f"
    basis -= basis_offset;
#line 378 "./finite_element.f"

#line 378 "./finite_element.f"
    /* Function Body */
#line 378 "./finite_element.f"
    left = (float)0.;
#line 379 "./finite_element.f"
    right = (float)0.;
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 381 "./finite_element.f"
    if (*nnodes != (*order - *continuity) * *nelements + *continuity + 1) {
#line 382 "./finite_element.f"
	abort_();
#line 383 "./finite_element.f"
    }
#line 385 "./finite_element.f"
    i__1 = *nnodes;
#line 385 "./finite_element.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 386 "./finite_element.f"
	soln[i__] = 0.;
/*       let subroutine mult put inhomogeneities into initial residual */
#line 388 "./finite_element.f"
	dirichlet[i__] = FALSE_;
#line 389 "./finite_element.f"
    }
/*     dirichlet boundary conditions */
#line 392 "./finite_element.f"
    j = 1;
#line 393 "./finite_element.f"
    soln[j] = left;
#line 394 "./finite_element.f"
    j = *nnodes - *continuity;
#line 395 "./finite_element.f"
    soln[j] = right;
/*     do i=1,nnodes */
/*       print *, "soln[",i,"] = ",soln(i) */
/*     enddo */
#line 401 "./finite_element.f"
    i__1 = *nelements;
#line 401 "./finite_element.f"
    for (element = 1; element <= i__1; ++element) {
#line 402 "./finite_element.f"
	dx = mesh[element + 1] - mesh[element];
#line 403 "./finite_element.f"
	i__2 = *order;
#line 403 "./finite_element.f"
	for (ngauss = 1; ngauss <= i__2; ++ngauss) {
#line 404 "./finite_element.f"
	    xg = mesh[element] + gauss[ngauss] * dx;
#line 405 "./finite_element.f"
	    pgauss[element + ngauss * pgauss_dim1] = 1. * weight[ngauss] / dx;
#line 406 "./finite_element.f"
	    rgauss[element + ngauss * rgauss_dim1] = 0. * weight[ngauss] * dx;
#line 407 "./finite_element.f"
	}
#line 408 "./finite_element.f"
    }
#line 410 "./finite_element.f"
    mult_(nelements, nnodes, order, &basis[basis_offset], &basis_deriv__[
	    basis_deriv_offset], &dirichlet[1], &element_to_node__[
	    element_to_node_offset], &pgauss[pgauss_offset], &rgauss[
	    rgauss_offset], &soln[1], &residual[1]);
#line 413 "./finite_element.f"
    i__1 = *nnodes;
#line 413 "./finite_element.f"
    for (node_i__ = 1; node_i__ <= i__1; ++node_i__) {
#line 414 "./finite_element.f"
	residual[node_i__] = -residual[node_i__];
#line 415 "./finite_element.f"
    }
#line 417 "./finite_element.f"
    i__1 = *nelements;
#line 417 "./finite_element.f"
    for (element = 1; element <= i__1; ++element) {
#line 418 "./finite_element.f"
	dx = mesh[element + 1] - mesh[element];
#line 419 "./finite_element.f"
	i__2 = *order;
#line 419 "./finite_element.f"
	for (ngauss = 1; ngauss <= i__2; ++ngauss) {
#line 420 "./finite_element.f"
	    xg = mesh[element] + gauss[ngauss] * dx;
/* Computing 2nd power */
#line 421 "./finite_element.f"
	    d__1 = machine_1.pi;
#line 421 "./finite_element.f"
	    fg = sin(machine_1.pi * xg) * (d__1 * d__1) * weight[ngauss] * dx;
#line 422 "./finite_element.f"
	    i__3 = *order + 1;
#line 422 "./finite_element.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 423 "./finite_element.f"
		node_i__ = element_to_node__[element + i__ * 
			element_to_node_dim1];
#line 424 "./finite_element.f"
		residual[node_i__] += basis[i__ + ngauss * basis_dim1] * fg;
#line 425 "./finite_element.f"
	    }
#line 426 "./finite_element.f"
	}
#line 427 "./finite_element.f"
    }
/*     essential boundary conditions: */
#line 430 "./finite_element.f"
    j = 1;
#line 431 "./finite_element.f"
    dirichlet[j] = TRUE_;
#line 432 "./finite_element.f"
    if (dirichlet[j]) {
#line 432 "./finite_element.f"
	residual[j] = 0.;
#line 432 "./finite_element.f"
    }
#line 434 "./finite_element.f"
    j = *nnodes - *continuity;
#line 435 "./finite_element.f"
    dirichlet[j] = TRUE_;
#line 436 "./finite_element.f"
    if (dirichlet[j]) {
#line 436 "./finite_element.f"
	residual[j] = 0.;
#line 436 "./finite_element.f"
    }
/*     do i=1,nnodes */
/*       print *, "residual[",i,"] = ",residual(i) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 442 "./finite_element.f"
    return 0;
} /* initialize_ */

/* *********************************************************************** */
/* Subroutine */ int pre_(nnodes, x, px)
integer *nnodes;
doublereal *x, *px;
{
    extern /* Subroutine */ int dcopy_();

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 453 "./finite_element.f"
    /* Parameter adjustments */
#line 453 "./finite_element.f"
    --px;
#line 453 "./finite_element.f"
    --x;
#line 453 "./finite_element.f"

#line 453 "./finite_element.f"
    /* Function Body */
#line 453 "./finite_element.f"
    dcopy_(nnodes, &x[1], &c__1, &px[1], &c__1);
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 455 "./finite_element.f"
    return 0;
} /* pre_ */

/* *********************************************************************** */
/* Subroutine */ int precg_(limit, ndigit, nelements, nnodes, order, basis, 
	basis_deriv__, dirichlet, element_to_node__, pgauss, rgauss, b, x, w)
integer *limit, *ndigit, *nelements, *nnodes, *order;
doublereal *basis, *basis_deriv__;
logical1 *dirichlet;
integer *element_to_node__;
doublereal *pgauss, *rgauss, *b, *x, *w;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, basis_dim1, 
	    basis_offset, basis_deriv_dim1, basis_deriv_offset, pgauss_dim1, 
	    pgauss_offset, rgauss_dim1, rgauss_offset, w_dim1, w_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal r__, s, t;
    static integer it;
    static doublereal dif;
    extern /* Subroutine */ int pre_();
    extern doublereal ddot_();
    static doublereal size;
    extern /* Subroutine */ int mult_();
    extern doublereal dasum_();
    extern /* Subroutine */ int daxpy_(), stopit_();

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     initialize variables */
/*     print *, "limit,ndigit = ",limit,ndigit */
/*     print *, "nelements = ",nelements */
/*     do i=1,nnodes */
/*       print *,"x,b[",i,"] = ",x(i),b(i) */
/*     enddo */
/*     call mult(nnodes,order, */
/*    &  dirichlet,matrix,x, */
/*    &  w) */
/*     do i=1,nnodes */
/*       print *,"Ax[",i,"] = ",w(i,1) */
/*     enddo */
/*     assume that b contains the initial residual */
#line 494 "./finite_element.f"
    /* Parameter adjustments */
#line 494 "./finite_element.f"
    w_dim1 = *nnodes;
#line 494 "./finite_element.f"
    w_offset = 1 + w_dim1 * 1;
#line 494 "./finite_element.f"
    w -= w_offset;
#line 494 "./finite_element.f"
    --x;
#line 494 "./finite_element.f"
    --b;
#line 494 "./finite_element.f"
    --dirichlet;
#line 494 "./finite_element.f"
    rgauss_dim1 = *nelements;
#line 494 "./finite_element.f"
    rgauss_offset = 1 + rgauss_dim1 * 1;
#line 494 "./finite_element.f"
    rgauss -= rgauss_offset;
#line 494 "./finite_element.f"
    pgauss_dim1 = *nelements;
#line 494 "./finite_element.f"
    pgauss_offset = 1 + pgauss_dim1 * 1;
#line 494 "./finite_element.f"
    pgauss -= pgauss_offset;
#line 494 "./finite_element.f"
    element_to_node_dim1 = *nelements;
#line 494 "./finite_element.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 1;
#line 494 "./finite_element.f"
    element_to_node__ -= element_to_node_offset;
#line 494 "./finite_element.f"
    basis_deriv_dim1 = *order + 1;
#line 494 "./finite_element.f"
    basis_deriv_offset = 1 + basis_deriv_dim1 * 1;
#line 494 "./finite_element.f"
    basis_deriv__ -= basis_deriv_offset;
#line 494 "./finite_element.f"
    basis_dim1 = *order + 1;
#line 494 "./finite_element.f"
    basis_offset = 1 + basis_dim1 * 1;
#line 494 "./finite_element.f"
    basis -= basis_offset;
#line 494 "./finite_element.f"

#line 494 "./finite_element.f"
    /* Function Body */
#line 494 "./finite_element.f"
    i__1 = *nnodes;
#line 494 "./finite_element.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       w(i,1) = b(i) - w(i,1) */
#line 496 "./finite_element.f"
	w[i__ + w_dim1] = b[i__];
#line 497 "./finite_element.f"
    }
/*     do i=1,nnodes */
/*       print *,"b-Ax[",i,"] = ",w(i,1) */
/*     enddo */
#line 502 "./finite_element.f"
    pre_(nnodes, &w[w_offset], &w[(w_dim1 << 1) + 1]);
/*     do i=1,nnodes */
/*       print *,"P*(b-Ax)[",i,"] = ",w(i,2) */
/*     enddo */
#line 507 "./finite_element.f"
    s = ddot_(nnodes, &w[w_dim1 + 1], &c__1, &w[(w_dim1 << 1) + 1], &c__1);
/*     print *, "r . P r = ",s */
#line 510 "./finite_element.f"
    if (s <= 0.) {
#line 510 "./finite_element.f"
	return 0;
#line 510 "./finite_element.f"
    }
/*     start iterations */
#line 512 "./finite_element.f"
    dif = 1.;
#line 513 "./finite_element.f"
    it = 0;
#line 514 "./finite_element.f"
    while(dif > 0.) {
#line 515 "./finite_element.f"
	mult_(nelements, nnodes, order, &basis[basis_offset], &basis_deriv__[
		basis_deriv_offset], &dirichlet[1], &element_to_node__[
		element_to_node_offset], &pgauss[pgauss_offset], &rgauss[
		rgauss_offset], &w[(w_dim1 << 1) + 1], &b[1]);
/*       print *, " " */
/*       print *, "it = ",it */
/*       do i=1,nnodes */
/*         print *,"Ap[",i,"] = ",b(i) */
/*       enddo */
#line 525 "./finite_element.f"
	t = ddot_(nnodes, &w[(w_dim1 << 1) + 1], &c__1, &b[1], &c__1);
/*       print *, "p . A p = ",t */
#line 528 "./finite_element.f"
	if (t <= 0.) {
#line 528 "./finite_element.f"
	    return 0;
#line 528 "./finite_element.f"
	}
#line 529 "./finite_element.f"
	r__ = s / t;
#line 530 "./finite_element.f"
	dif = 0.;
#line 531 "./finite_element.f"
	size = 0.;
#line 532 "./finite_element.f"
	daxpy_(nnodes, &r__, &w[(w_dim1 << 1) + 1], &c__1, &x[1], &c__1);
/*       do i=1,nnodes */
/*         print *,"x[",i,"] = ",x(i) */
/*       enddo */
#line 537 "./finite_element.f"
	d__1 = -r__;
#line 537 "./finite_element.f"
	daxpy_(nnodes, &d__1, &b[1], &c__1, &w[w_dim1 + 1], &c__1);
/*       do i=1,nnodes */
/*         print *,"r[",i,"] = ",w(i,1) */
/*       enddo */
#line 542 "./finite_element.f"
	dif = r__ * dasum_(nnodes, &w[(w_dim1 << 1) + 1], &c__1);
#line 543 "./finite_element.f"
	size = dasum_(nnodes, &x[1], &c__1);
/*       print *, "dif = ",dif */
/*       print *, "size = ",size */
#line 547 "./finite_element.f"
	pre_(nnodes, &w[w_offset], &b[1]);
/*       do i=1,nnodes */
/*         print *,"P r[",i,"] = ",b(i) */
/*       enddo */
#line 552 "./finite_element.f"
	t = ddot_(nnodes, &b[1], &c__1, &w[w_dim1 + 1], &c__1);
/*       print *, "r . P r = ",t */
#line 555 "./finite_element.f"
	r__ = t / s;
#line 556 "./finite_element.f"
	i__1 = *nnodes;
#line 556 "./finite_element.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 557 "./finite_element.f"
	    w[i__ + (w_dim1 << 1)] = b[i__] + r__ * w[i__ + (w_dim1 << 1)];
#line 558 "./finite_element.f"
	}
/*       do i=1,nnodes */
/*         print *,"p[",i,"] = ",w(i,2) */
/*       enddo */
#line 563 "./finite_element.f"
	s = t;
#line 564 "./finite_element.f"
	stopit_(&dif, &size, ndigit, limit);
/*       print *, "dif = ",dif */
#line 566 "./finite_element.f"
	++it;
#line 567 "./finite_element.f"
    }
/*     do i=1,nnodes */
/*       print *,"x[",i,"] = ",x(i) */
/*     enddo */
/*     call mult(nelements,nnodes,order, */
/*    &  basis,basis_deriv,dirichlet,element_to_node,pgauss,rgauss, */
/*    &    x, */
/*    &  b) */
/*     do i=1,nnodes */
/*       print *,"Ax[",i,"] = ",b(i) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 579 "./finite_element.f"
    return 0;
} /* precg_ */

/* *********************************************************************** */
/* Subroutine */ int stopit_(dif, size, ndigit, limit)
doublereal *dif, *size;
integer *ndigit, *limit;
{
    /* Initialized data */

    static integer i__ = 0;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double pow_ri();

    /* Local variables */
    static doublereal e, t;

#line 589 "./finite_element.f"
    *dif = abs(*dif);
#line 590 "./finite_element.f"
    *size = abs(*size);
/*     initialization during first iteration */
#line 592 "./finite_element.f"
    if (i__ <= 0) {
#line 592 "./finite_element.f"
	i__1 = -(*ndigit);
#line 592 "./finite_element.f"
	t = pow_ri(&c_b24, &i__1);
#line 592 "./finite_element.f"
    }
#line 593 "./finite_element.f"
    ++i__;
#line 594 "./finite_element.f"
    e = (doublereal) i__ * 3.;
/*     print *, "i,t,e = ",i,t,e */
#line 597 "./finite_element.f"
    if (*dif <= t * *size) {
/*       stopping criterion I */
#line 599 "./finite_element.f"
	e += 1.;
#line 600 "./finite_element.f"
	*dif = -(*dif);
#line 601 "./finite_element.f"
	i__ = 0;
#line 602 "./finite_element.f"
    } else if (i__ >= *limit) {
/*       stopping criterion II */
#line 604 "./finite_element.f"
	e += 2.;
#line 605 "./finite_element.f"
	*dif = -(*dif);
#line 606 "./finite_element.f"
	i__ = 0;
#line 607 "./finite_element.f"
    }
#line 608 "./finite_element.f"
    *size = e;
#line 609 "./finite_element.f"
    return 0;
} /* stopit_ */

