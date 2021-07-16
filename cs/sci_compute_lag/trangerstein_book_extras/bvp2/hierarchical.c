#line 1 "./hierarchical.f"
/*  -- translated by f2c (version 20000817).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

#line 1 "./hierarchical.f"
/* Common Block Declarations */

struct {
    doublereal roundoff, small, huge__, undefind, pi;
} machine_;

#define machine_1 machine_

/* Table of constant values */

static integer c__1 = 1;
static real c_b22 = (float)10.;

/* Subroutine */ int canonical_(order, basis, basis_deriv__, lobatto, weight)
integer *order;
doublereal *basis, *basis_deriv__, *lobatto, *weight;
{
    /* System generated locals */
    integer basis_dim1, basis_offset, basis_deriv_dim1, basis_deriv_offset, 
	    i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal legendre[11];
    static integer i__, j;
    static doublereal w, jl, j2m1;
    extern /* Subroutine */ int abort_();
    static doublereal rt2j2m1, legendre_deriv__[11];

/*     integer i */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 14 "./hierarchical.f"
    /* Parameter adjustments */
#line 14 "./hierarchical.f"
    basis_deriv_dim1 = *order - 0 + 1;
#line 14 "./hierarchical.f"
    basis_deriv_offset = 0 + basis_deriv_dim1 * 0;
#line 14 "./hierarchical.f"
    basis_deriv__ -= basis_deriv_offset;
#line 14 "./hierarchical.f"
    basis_dim1 = *order - 0 + 1;
#line 14 "./hierarchical.f"
    basis_offset = 0 + basis_dim1 * 0;
#line 14 "./hierarchical.f"
    basis -= basis_offset;
#line 14 "./hierarchical.f"

#line 14 "./hierarchical.f"
    /* Function Body */
#line 14 "./hierarchical.f"
    if (*order < 1 || *order > 10) {
#line 14 "./hierarchical.f"
	abort_();
#line 14 "./hierarchical.f"
    }
/*     lobatto holds zeros of derivative of legendre polynomial */
/*       of degree = order */
#line 18 "./hierarchical.f"
    lobatto[0] = -1.;
#line 19 "./hierarchical.f"
    lobatto[*order] = 1.;
#line 20 "./hierarchical.f"
    if (*order == 2) {
#line 21 "./hierarchical.f"
	lobatto[1] = 0.;
#line 22 "./hierarchical.f"
    } else if (*order == 3) {
#line 23 "./hierarchical.f"
	lobatto[2] = sqrt(.2);
#line 24 "./hierarchical.f"
	lobatto[1] = -sqrt(.2);
#line 25 "./hierarchical.f"
    } else if (*order == 4) {
#line 26 "./hierarchical.f"
	lobatto[3] = sqrt(21.) / 7.;
#line 27 "./hierarchical.f"
	lobatto[2] = 0.;
#line 28 "./hierarchical.f"
	lobatto[1] = -lobatto[3];
#line 29 "./hierarchical.f"
    } else if (*order == 5) {
#line 30 "./hierarchical.f"
	lobatto[4] = sqrt(sqrt(7.) * 42. + 147.) / 21.;
#line 31 "./hierarchical.f"
	lobatto[3] = sqrt(147. - sqrt(7.) * 42.) / 21.;
#line 32 "./hierarchical.f"
	lobatto[2] = -lobatto[3];
#line 33 "./hierarchical.f"
	lobatto[1] = -lobatto[4];
#line 34 "./hierarchical.f"
    } else if (*order == 6) {
#line 35 "./hierarchical.f"
	lobatto[5] = sqrt(sqrt(15.) * 66. + 495.) / 33.;
#line 36 "./hierarchical.f"
	lobatto[4] = sqrt(495. - sqrt(15.) * 66.) / 33.;
#line 37 "./hierarchical.f"
	lobatto[3] = 0.;
#line 38 "./hierarchical.f"
	lobatto[2] = -lobatto[4];
#line 39 "./hierarchical.f"
	lobatto[1] = -lobatto[5];
#line 40 "./hierarchical.f"
    } else if (*order == 7) {
#line 41 "./hierarchical.f"
	lobatto[6] = .87174014850960661533;
#line 42 "./hierarchical.f"
	lobatto[5] = .59170018143314230214;
#line 43 "./hierarchical.f"
	lobatto[4] = .20929921790247886873;
#line 44 "./hierarchical.f"
	lobatto[3] = -lobatto[4];
#line 45 "./hierarchical.f"
	lobatto[2] = -lobatto[5];
#line 46 "./hierarchical.f"
	lobatto[1] = -lobatto[6];
#line 47 "./hierarchical.f"
    } else if (*order == 8) {
#line 48 "./hierarchical.f"
	lobatto[7] = .89975799541146015735;
#line 49 "./hierarchical.f"
	lobatto[6] = .67718627951073775343;
#line 50 "./hierarchical.f"
	lobatto[5] = .36311746382617815869;
#line 51 "./hierarchical.f"
	lobatto[4] = 0.;
#line 52 "./hierarchical.f"
	lobatto[3] = -lobatto[5];
#line 53 "./hierarchical.f"
	lobatto[2] = -lobatto[6];
#line 54 "./hierarchical.f"
	lobatto[1] = -lobatto[7];
#line 55 "./hierarchical.f"
    } else if (*order == 9) {
#line 56 "./hierarchical.f"
	lobatto[8] = .91953390816645881383;
#line 57 "./hierarchical.f"
	lobatto[7] = .73877386510550507501;
#line 58 "./hierarchical.f"
	lobatto[6] = .47792494981044449566;
#line 59 "./hierarchical.f"
	lobatto[5] = .16527895766638702463;
#line 60 "./hierarchical.f"
	lobatto[4] = -lobatto[5];
#line 61 "./hierarchical.f"
	lobatto[3] = -lobatto[6];
#line 62 "./hierarchical.f"
	lobatto[2] = -lobatto[7];
#line 63 "./hierarchical.f"
	lobatto[1] = -lobatto[8];
#line 64 "./hierarchical.f"
    } else if (*order == 10) {
#line 65 "./hierarchical.f"
	lobatto[9] = .93400143040805913433;
#line 66 "./hierarchical.f"
	lobatto[8] = .78448347366314441862;
#line 67 "./hierarchical.f"
	lobatto[7] = .56523532699620500647;
#line 68 "./hierarchical.f"
	lobatto[6] = .29575813558693939143;
#line 69 "./hierarchical.f"
	lobatto[5] = 0.;
#line 70 "./hierarchical.f"
	lobatto[4] = -lobatto[6];
#line 71 "./hierarchical.f"
	lobatto[3] = -lobatto[7];
#line 72 "./hierarchical.f"
	lobatto[2] = -lobatto[8];
#line 73 "./hierarchical.f"
	lobatto[1] = -lobatto[9];
#line 74 "./hierarchical.f"
    }
#line 76 "./hierarchical.f"
    legendre[0] = 1.;
#line 77 "./hierarchical.f"
    legendre_deriv__[0] = 0.;
#line 78 "./hierarchical.f"
    legendre_deriv__[1] = 1.;
#line 79 "./hierarchical.f"
    w = 2. / (doublereal) (*order * (*order + 1));
#line 80 "./hierarchical.f"
    i__1 = *order;
#line 80 "./hierarchical.f"
    for (j = 0; j <= i__1; ++j) {
#line 81 "./hierarchical.f"
	legendre[1] = lobatto[j];
#line 82 "./hierarchical.f"
	basis[j * basis_dim1] = (1. - lobatto[j]) * .5;
#line 83 "./hierarchical.f"
	basis[j * basis_dim1 + 1] = (lobatto[j] + 1.) * .5;
#line 84 "./hierarchical.f"
	basis_deriv__[j * basis_deriv_dim1] = -.5;
#line 85 "./hierarchical.f"
	basis_deriv__[j * basis_deriv_dim1 + 1] = .5;
#line 86 "./hierarchical.f"
	i__2 = *order;
#line 86 "./hierarchical.f"
	for (i__ = 2; i__ <= i__2; ++i__) {
#line 87 "./hierarchical.f"
	    j2m1 = (doublereal) i__ * 2. - 1;
#line 88 "./hierarchical.f"
	    jl = j2m1 * legendre[i__ - 1];
#line 89 "./hierarchical.f"
	    rt2j2m1 = 1. / sqrt(j2m1 * 2.);
#line 90 "./hierarchical.f"
	    legendre[i__] = (jl * lobatto[j] - (doublereal) (i__ - 1) * 
		    legendre[i__ - 2]) / (doublereal) i__;
#line 91 "./hierarchical.f"
	    legendre_deriv__[i__] = jl + legendre_deriv__[i__ - 2];
#line 92 "./hierarchical.f"
	    basis[i__ + j * basis_dim1] = (legendre[i__] - legendre[i__ - 2]) 
		    * rt2j2m1;
#line 93 "./hierarchical.f"
	    basis_deriv__[i__ + j * basis_deriv_dim1] = (legendre_deriv__[i__]
		     - legendre_deriv__[i__ - 2]) * rt2j2m1;
#line 95 "./hierarchical.f"
	}
/* Computing 2nd power */
#line 96 "./hierarchical.f"
	d__1 = legendre[*order];
#line 96 "./hierarchical.f"
	weight[j] = w / (d__1 * d__1);
#line 97 "./hierarchical.f"
    }
/*     do i=0,order */
/*       print *, "lobatto[",i,"] = ",lobatto(i) */
/*     enddo */
/*     do i=0,order */
/*       print *, "weight[",i,"] = ",weight(i) */
/*     enddo */
/*     do i=0,order */
/*       print *, "basis[",i,"] = ",(basis(i,j),j=0,order) */
/*     enddo */
/*     do i=0,order */
/*       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=0,order) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 111 "./hierarchical.f"
    return 0;
} /* canonical_ */

/* *********************************************************************** */
/* Subroutine */ int grid_(nelements, nnodes, order, element_to_node__, mesh)
integer *nelements, *nnodes, *order, *element_to_node__;
doublereal *mesh;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal dx;
    static integer node;
    extern /* Subroutine */ int abort_();
    static integer element;

/*     integer j */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 126 "./hierarchical.f"
    /* Parameter adjustments */
#line 126 "./hierarchical.f"
    element_to_node_dim1 = *nelements;
#line 126 "./hierarchical.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 0;
#line 126 "./hierarchical.f"
    element_to_node__ -= element_to_node_offset;
#line 126 "./hierarchical.f"

#line 126 "./hierarchical.f"
    /* Function Body */
#line 126 "./hierarchical.f"
    if (*nnodes != *order * *nelements + 1) {
#line 126 "./hierarchical.f"
	abort_();
#line 126 "./hierarchical.f"
    }
#line 128 "./hierarchical.f"
    dx = 1. / (doublereal) (*nelements);
#line 129 "./hierarchical.f"
    i__1 = *nelements;
#line 129 "./hierarchical.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 130 "./hierarchical.f"
	mesh[i__] = (doublereal) i__ * dx;
#line 131 "./hierarchical.f"
    }
#line 133 "./hierarchical.f"
    node = 1;
#line 134 "./hierarchical.f"
    i__1 = *nelements;
#line 134 "./hierarchical.f"
    for (element = 1; element <= i__1; ++element) {
#line 135 "./hierarchical.f"
	element_to_node__[element] = node;
#line 136 "./hierarchical.f"
	element_to_node__[element + element_to_node_dim1] = node + *order;
#line 137 "./hierarchical.f"
	i__2 = *order - 1;
#line 137 "./hierarchical.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 138 "./hierarchical.f"
	    element_to_node__[element + (i__ + 1) * element_to_node_dim1] = 
		    node + i__;
#line 139 "./hierarchical.f"
	}
#line 140 "./hierarchical.f"
	node += *order;
#line 141 "./hierarchical.f"
    }
/*     do i=0,nelements */
/*       print *, "mesh[",i,"] = ",mesh(i) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "element_to_node[",i,"] = ", */
/*    &    (element_to_node(i,j),j=0,order) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 151 "./hierarchical.f"
    return 0;
} /* grid_ */

/* *********************************************************************** */
/* Subroutine */ int mult_(nelements, nnodes, order, basis, basis_deriv__, 
	dirichlet, element_to_node__, plobatto, rlobatto, x, ax)
integer *nelements, *nnodes, *order;
doublereal *basis, *basis_deriv__;
logical1 *dirichlet;
integer *element_to_node__;
doublereal *plobatto, *rlobatto, *x, *ax;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, basis_dim1, 
	    basis_offset, basis_deriv_dim1, basis_deriv_offset, plobatto_dim1,
	     plobatto_offset, rlobatto_dim1, rlobatto_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer nlobatto, i__, j;
    static doublereal bp, br;
    static integer node_i__, node_j__, element;

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     do i=0,order */
/*       print *, "basis[",i,"] = ",(basis(i,j),j=0,order) */
/*     enddo */
/*     do i=0,order */
/*       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=0,order) */
/*     enddo */
/*     do i=1,nnodes */
/*       print *, "dirichlet[",i,"] = ",dirichlet(i) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "element_to_node[",i,"] = ", */
/*    &    (element_to_node(i,j),j=0,order) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "plobatto[",i,"] = ",(plobatto(i,j),j=0,order) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "rlobatto[",i,"] = ",(rlobatto(i,j),j=0,order) */
/*     enddo */
/*     do i=1,nnodes */
/*       print *, "x[",i,"] = ",x(i) */
/*     enddo */
#line 198 "./hierarchical.f"
    /* Parameter adjustments */
#line 198 "./hierarchical.f"
    --ax;
#line 198 "./hierarchical.f"
    --x;
#line 198 "./hierarchical.f"
    --dirichlet;
#line 198 "./hierarchical.f"
    rlobatto_dim1 = *nelements;
#line 198 "./hierarchical.f"
    rlobatto_offset = 1 + rlobatto_dim1 * 0;
#line 198 "./hierarchical.f"
    rlobatto -= rlobatto_offset;
#line 198 "./hierarchical.f"
    plobatto_dim1 = *nelements;
#line 198 "./hierarchical.f"
    plobatto_offset = 1 + plobatto_dim1 * 0;
#line 198 "./hierarchical.f"
    plobatto -= plobatto_offset;
#line 198 "./hierarchical.f"
    element_to_node_dim1 = *nelements;
#line 198 "./hierarchical.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 0;
#line 198 "./hierarchical.f"
    element_to_node__ -= element_to_node_offset;
#line 198 "./hierarchical.f"
    basis_deriv_dim1 = *order - 0 + 1;
#line 198 "./hierarchical.f"
    basis_deriv_offset = 0 + basis_deriv_dim1 * 0;
#line 198 "./hierarchical.f"
    basis_deriv__ -= basis_deriv_offset;
#line 198 "./hierarchical.f"
    basis_dim1 = *order - 0 + 1;
#line 198 "./hierarchical.f"
    basis_offset = 0 + basis_dim1 * 0;
#line 198 "./hierarchical.f"
    basis -= basis_offset;
#line 198 "./hierarchical.f"

#line 198 "./hierarchical.f"
    /* Function Body */
#line 198 "./hierarchical.f"
    i__1 = *nnodes;
#line 198 "./hierarchical.f"
    for (node_i__ = 1; node_i__ <= i__1; ++node_i__) {
#line 199 "./hierarchical.f"
	ax[node_i__] = 0.;
#line 200 "./hierarchical.f"
    }
/*     do i=1,nnodes */
/*       print *,"x,Ax[",i,"] = ",x(i),Ax(i) */
/*     enddo */
#line 204 "./hierarchical.f"
    i__1 = *nelements;
#line 204 "./hierarchical.f"
    for (element = 1; element <= i__1; ++element) {
/*       print *, " " */
/*       print *, "element = ",element */
#line 207 "./hierarchical.f"
	i__2 = *order;
#line 207 "./hierarchical.f"
	for (nlobatto = 0; nlobatto <= i__2; ++nlobatto) {
/*         print *, "nlobatto = ",nlobatto */
#line 209 "./hierarchical.f"
	    i__3 = *order;
#line 209 "./hierarchical.f"
	    for (i__ = 0; i__ <= i__3; ++i__) {
#line 210 "./hierarchical.f"
		node_i__ = element_to_node__[element + i__ * 
			element_to_node_dim1];
#line 211 "./hierarchical.f"
		bp = basis_deriv__[i__ + nlobatto * basis_deriv_dim1] * 
			plobatto[element + nlobatto * plobatto_dim1];
#line 212 "./hierarchical.f"
		br = basis[i__ + nlobatto * basis_dim1] * rlobatto[element + 
			nlobatto * rlobatto_dim1];
/*           print *, "node_i = ",node_i */
/*           print *, "bp,br = ",bp,br */
#line 215 "./hierarchical.f"
		i__4 = *order;
#line 215 "./hierarchical.f"
		for (j = 0; j <= i__4; ++j) {
#line 216 "./hierarchical.f"
		    node_j__ = element_to_node__[element + j * 
			    element_to_node_dim1];
/*             print *, "x[",node_j,"] = ",x(node_j) */
#line 218 "./hierarchical.f"
		    ax[node_i__] += (bp * basis_deriv__[j + nlobatto * 
			    basis_deriv_dim1] + br * basis[j + nlobatto * 
			    basis_dim1]) * x[node_j__];
/*             print *, "Ax[",node_i,"] = ",Ax(node_i) */
#line 222 "./hierarchical.f"
		}
#line 223 "./hierarchical.f"
	    }
#line 224 "./hierarchical.f"
	}
#line 225 "./hierarchical.f"
    }
#line 226 "./hierarchical.f"
    i__1 = *nnodes;
#line 226 "./hierarchical.f"
    for (node_i__ = 1; node_i__ <= i__1; ++node_i__) {
#line 227 "./hierarchical.f"
	if (dirichlet[node_i__]) {
#line 227 "./hierarchical.f"
	    ax[node_i__] = 0.;
#line 227 "./hierarchical.f"
	}
#line 228 "./hierarchical.f"
    }
/*     do i=1,nnodes */
/*       print *,"x,Ax[",i,"] = ",x(i),Ax(i) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 233 "./hierarchical.f"
    return 0;
} /* mult_ */

/* *********************************************************************** */
/* Subroutine */ int initialize_(nelements, nnodes, order, basis, 
	basis_deriv__, element_to_node__, lobatto, mesh, weight, dirichlet, 
	plobatto, rlobatto, residual, soln)
integer *nelements, *nnodes, *order;
doublereal *basis, *basis_deriv__;
integer *element_to_node__;
doublereal *lobatto, *mesh, *weight;
logical1 *dirichlet;
doublereal *plobatto, *rlobatto, *residual, *soln;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, basis_dim1, 
	    basis_offset, basis_deriv_dim1, basis_deriv_offset, plobatto_dim1,
	     plobatto_offset, rlobatto_dim1, rlobatto_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Builtin functions */
    double sin();

    /* Local variables */
    static integer nlobatto, i__;
    static doublereal fg, dx, xg, left;
    extern /* Subroutine */ int mult_(), abort_();
    static doublereal right;
    static integer node_i__, element;

/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     integer j */
/*     bvp: */
/*       - d/dx( p du/dx ) + r u = f, 0 < x < 1 */
/*       u(0)=left, u(1)=right */
#line 274 "./hierarchical.f"
    /* Parameter adjustments */
#line 274 "./hierarchical.f"
    --soln;
#line 274 "./hierarchical.f"
    --residual;
#line 274 "./hierarchical.f"
    --dirichlet;
#line 274 "./hierarchical.f"
    rlobatto_dim1 = *nelements;
#line 274 "./hierarchical.f"
    rlobatto_offset = 1 + rlobatto_dim1 * 0;
#line 274 "./hierarchical.f"
    rlobatto -= rlobatto_offset;
#line 274 "./hierarchical.f"
    plobatto_dim1 = *nelements;
#line 274 "./hierarchical.f"
    plobatto_offset = 1 + plobatto_dim1 * 0;
#line 274 "./hierarchical.f"
    plobatto -= plobatto_offset;
#line 274 "./hierarchical.f"
    element_to_node_dim1 = *nelements;
#line 274 "./hierarchical.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 0;
#line 274 "./hierarchical.f"
    element_to_node__ -= element_to_node_offset;
#line 274 "./hierarchical.f"
    basis_deriv_dim1 = *order - 0 + 1;
#line 274 "./hierarchical.f"
    basis_deriv_offset = 0 + basis_deriv_dim1 * 0;
#line 274 "./hierarchical.f"
    basis_deriv__ -= basis_deriv_offset;
#line 274 "./hierarchical.f"
    basis_dim1 = *order - 0 + 1;
#line 274 "./hierarchical.f"
    basis_offset = 0 + basis_dim1 * 0;
#line 274 "./hierarchical.f"
    basis -= basis_offset;
#line 274 "./hierarchical.f"

#line 274 "./hierarchical.f"
    /* Function Body */
#line 274 "./hierarchical.f"
    left = (float)0.;
#line 275 "./hierarchical.f"
    right = (float)0.;
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     do i=0,order */
/*       print *, "basis[",i,"] = ",(basis(i,j),j=0,order) */
/*     enddo */
/*     do i=0,order */
/*       print *, "basis_deriv[",i,"] = ",(basis_deriv(i,j),j=0,order) */
/*     enddo */
/*     do i=1,nelements */
/*       print *, "element_to_node[",i,"] = ", */
/*    &    (element_to_node(i,j),j=0,order) */
/*     enddo */
/*     do i=0,order */
/*       print *, "lobatto[",i,"] = ",lobatto(i) */
/*     enddo */
/*     do i=0,nelements */
/*       print *, "mesh[",i,"] = ",mesh(i) */
/*     enddo */
/*     do i=0,order */
/*       print *, "weight[",i,"] = ",weight(i) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 297 "./hierarchical.f"
    if (*nnodes != *order * *nelements + 1) {
#line 298 "./hierarchical.f"
	abort_();
#line 299 "./hierarchical.f"
    }
#line 301 "./hierarchical.f"
    i__1 = *nnodes;
#line 301 "./hierarchical.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 302 "./hierarchical.f"
	soln[i__] = 0.;
/*       let subroutine mult put inhomogeneities into initial residual */
#line 304 "./hierarchical.f"
	dirichlet[i__] = FALSE_;
#line 305 "./hierarchical.f"
    }
/*     dirichlet boundary conditions */
#line 308 "./hierarchical.f"
    soln[1] = left;
#line 309 "./hierarchical.f"
    soln[*nnodes] = right;
/*     do i=1,nnodes */
/*       print *, "soln[",i,"] = ",soln(i) */
/*     enddo */
#line 315 "./hierarchical.f"
    i__1 = *nelements;
#line 315 "./hierarchical.f"
    for (element = 1; element <= i__1; ++element) {
#line 316 "./hierarchical.f"
	dx = (mesh[element] - mesh[element - 1]) * .5;
#line 317 "./hierarchical.f"
	i__2 = *order;
#line 317 "./hierarchical.f"
	for (nlobatto = 0; nlobatto <= i__2; ++nlobatto) {
#line 318 "./hierarchical.f"
	    xg = mesh[element - 1] + (lobatto[nlobatto] + 1.) * dx;
#line 319 "./hierarchical.f"
	    plobatto[element + nlobatto * plobatto_dim1] = 1. * weight[
		    nlobatto] / dx;
#line 320 "./hierarchical.f"
	    rlobatto[element + nlobatto * rlobatto_dim1] = 0. * weight[
		    nlobatto] * dx;
#line 321 "./hierarchical.f"
	}
#line 322 "./hierarchical.f"
    }
#line 324 "./hierarchical.f"
    mult_(nelements, nnodes, order, &basis[basis_offset], &basis_deriv__[
	    basis_deriv_offset], &dirichlet[1], &element_to_node__[
	    element_to_node_offset], &plobatto[plobatto_offset], &rlobatto[
	    rlobatto_offset], &soln[1], &residual[1]);
#line 328 "./hierarchical.f"
    i__1 = *nnodes;
#line 328 "./hierarchical.f"
    for (node_i__ = 1; node_i__ <= i__1; ++node_i__) {
#line 329 "./hierarchical.f"
	residual[node_i__] = -residual[node_i__];
#line 330 "./hierarchical.f"
    }
#line 332 "./hierarchical.f"
    i__1 = *nelements;
#line 332 "./hierarchical.f"
    for (element = 1; element <= i__1; ++element) {
#line 333 "./hierarchical.f"
	dx = (mesh[element] - mesh[element - 1]) * .5;
/*       print *, " " */
/*       print *, "element = ",element */
#line 336 "./hierarchical.f"
	i__2 = *order;
#line 336 "./hierarchical.f"
	for (nlobatto = 0; nlobatto <= i__2; ++nlobatto) {
#line 337 "./hierarchical.f"
	    xg = mesh[element - 1] + (lobatto[nlobatto] + 1.) * dx;
/* Computing 2nd power */
#line 338 "./hierarchical.f"
	    d__1 = machine_1.pi;
#line 338 "./hierarchical.f"
	    fg = sin(machine_1.pi * xg) * (d__1 * d__1) * weight[nlobatto] * 
		    dx;
/*         print *, "fg(",xg,") = ",fg */
#line 340 "./hierarchical.f"
	    i__3 = *order;
#line 340 "./hierarchical.f"
	    for (i__ = 0; i__ <= i__3; ++i__) {
#line 341 "./hierarchical.f"
		node_i__ = element_to_node__[element + i__ * 
			element_to_node_dim1];
#line 342 "./hierarchical.f"
		residual[node_i__] += basis[i__ + nlobatto * basis_dim1] * fg;
/*           print *, "basis(",i,",",nlobatto,") = ",basis(i,nlobatto) */
/*           print *, "residual(",node_i,") = ",residual(node_i) */
#line 345 "./hierarchical.f"
	    }
#line 346 "./hierarchical.f"
	}
#line 347 "./hierarchical.f"
    }
/*     essential boundary conditions: */
#line 350 "./hierarchical.f"
    dirichlet[1] = TRUE_;
#line 351 "./hierarchical.f"
    if (dirichlet[1]) {
#line 351 "./hierarchical.f"
	residual[1] = 0.;
#line 351 "./hierarchical.f"
    }
#line 353 "./hierarchical.f"
    dirichlet[*nnodes] = TRUE_;
#line 354 "./hierarchical.f"
    if (dirichlet[*nnodes]) {
#line 354 "./hierarchical.f"
	residual[*nnodes] = 0.;
#line 354 "./hierarchical.f"
    }
/*     do i=1,nnodes */
/*       print *, "residual[",i,"] = ",residual(i) */
/*     enddo */
/*     do j=1,nnodes */
/*       do i=1,nnodes */
/*         soln(i)=0.d0 */
/*         dirichlet(i)=.false. */
/*       enddo */
/*       dirichlet(1)=.true. */
/*       dirichlet(nnodes)=.true. */
/*       soln(j)=1.d0 */
/*       call mult(nelements,nnodes,order, */
/*    &    basis,basis_deriv,dirichlet,element_to_node,plobatto,rlobatto, */
/*    &      soln, */
/*    &    residual) */
/*       print *, " " */
/*       print *, "j = ",j */
/*       do i=1,nnodes */
/*         print *, "residual[",i,"] = ",residual(i) */
/*       enddo */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 379 "./hierarchical.f"
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
#line 390 "./hierarchical.f"
    /* Parameter adjustments */
#line 390 "./hierarchical.f"
    --px;
#line 390 "./hierarchical.f"
    --x;
#line 390 "./hierarchical.f"

#line 390 "./hierarchical.f"
    /* Function Body */
#line 390 "./hierarchical.f"
    dcopy_(nnodes, &x[1], &c__1, &px[1], &c__1);
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 392 "./hierarchical.f"
    return 0;
} /* pre_ */

/* *********************************************************************** */
/* Subroutine */ int precg_(limit, ndigit, nelements, nnodes, order, basis, 
	basis_deriv__, dirichlet, element_to_node__, plobatto, rlobatto, b, x,
	 w)
integer *limit, *ndigit, *nelements, *nnodes, *order;
doublereal *basis, *basis_deriv__;
logical1 *dirichlet;
integer *element_to_node__;
doublereal *plobatto, *rlobatto, *b, *x, *w;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, basis_dim1, 
	    basis_offset, basis_deriv_dim1, basis_deriv_offset, plobatto_dim1,
	     plobatto_offset, rlobatto_dim1, rlobatto_offset, w_dim1, 
	    w_offset, i__1;
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
/*     print *, "limit,ndigit = ",limit,ndigit */
/*     print *, "nelements = ",nelements */
/*     do i=1,nnodes */
/*       print *,"x,b[",i,"] = ",x(i),b(i) */
/*     enddo */
/*     initialize variables */
/*     assume that b contains the initial residual */
#line 427 "./hierarchical.f"
    /* Parameter adjustments */
#line 427 "./hierarchical.f"
    w_dim1 = *nnodes;
#line 427 "./hierarchical.f"
    w_offset = 1 + w_dim1 * 1;
#line 427 "./hierarchical.f"
    w -= w_offset;
#line 427 "./hierarchical.f"
    --x;
#line 427 "./hierarchical.f"
    --b;
#line 427 "./hierarchical.f"
    --dirichlet;
#line 427 "./hierarchical.f"
    rlobatto_dim1 = *nelements;
#line 427 "./hierarchical.f"
    rlobatto_offset = 1 + rlobatto_dim1 * 0;
#line 427 "./hierarchical.f"
    rlobatto -= rlobatto_offset;
#line 427 "./hierarchical.f"
    plobatto_dim1 = *nelements;
#line 427 "./hierarchical.f"
    plobatto_offset = 1 + plobatto_dim1 * 0;
#line 427 "./hierarchical.f"
    plobatto -= plobatto_offset;
#line 427 "./hierarchical.f"
    element_to_node_dim1 = *nelements;
#line 427 "./hierarchical.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 0;
#line 427 "./hierarchical.f"
    element_to_node__ -= element_to_node_offset;
#line 427 "./hierarchical.f"
    basis_deriv_dim1 = *order - 0 + 1;
#line 427 "./hierarchical.f"
    basis_deriv_offset = 0 + basis_deriv_dim1 * 0;
#line 427 "./hierarchical.f"
    basis_deriv__ -= basis_deriv_offset;
#line 427 "./hierarchical.f"
    basis_dim1 = *order - 0 + 1;
#line 427 "./hierarchical.f"
    basis_offset = 0 + basis_dim1 * 0;
#line 427 "./hierarchical.f"
    basis -= basis_offset;
#line 427 "./hierarchical.f"

#line 427 "./hierarchical.f"
    /* Function Body */
#line 427 "./hierarchical.f"
    i__1 = *nnodes;
#line 427 "./hierarchical.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       w(i,1) = b(i) - w(i,1) */
#line 429 "./hierarchical.f"
	w[i__ + w_dim1] = b[i__];
#line 430 "./hierarchical.f"
    }
/*     do i=1,nnodes */
/*       print *,"b-Ax[",i,"] = ",w(i,1) */
/*     enddo */
#line 435 "./hierarchical.f"
    pre_(nnodes, &w[w_offset], &w[(w_dim1 << 1) + 1]);
/*     do i=1,nnodes */
/*       print *,"P*(b-Ax)[",i,"] = ",w(i,2) */
/*     enddo */
#line 440 "./hierarchical.f"
    s = ddot_(nnodes, &w[w_dim1 + 1], &c__1, &w[(w_dim1 << 1) + 1], &c__1);
/*     print *, "r . P r = ",s */
#line 443 "./hierarchical.f"
    if (s <= 0.) {
#line 443 "./hierarchical.f"
	return 0;
#line 443 "./hierarchical.f"
    }
/*     start iterations */
#line 445 "./hierarchical.f"
    dif = 1.;
#line 446 "./hierarchical.f"
    it = 0;
#line 447 "./hierarchical.f"
    while(dif > 0.) {
#line 448 "./hierarchical.f"
	mult_(nelements, nnodes, order, &basis[basis_offset], &basis_deriv__[
		basis_deriv_offset], &dirichlet[1], &element_to_node__[
		element_to_node_offset], &plobatto[plobatto_offset], &
		rlobatto[rlobatto_offset], &w[(w_dim1 << 1) + 1], &b[1]);
/*       print *, " " */
/*       print *, "it = ",it */
/*       do i=1,nnodes */
/*         print *,"Ap[",i,"] = ",b(i) */
/*       enddo */
#line 458 "./hierarchical.f"
	t = ddot_(nnodes, &w[(w_dim1 << 1) + 1], &c__1, &b[1], &c__1);
/*       print *, "p . A p = ",t */
#line 461 "./hierarchical.f"
	if (t <= 0.) {
#line 461 "./hierarchical.f"
	    return 0;
#line 461 "./hierarchical.f"
	}
#line 462 "./hierarchical.f"
	r__ = s / t;
#line 463 "./hierarchical.f"
	dif = 0.;
#line 464 "./hierarchical.f"
	size = 0.;
#line 465 "./hierarchical.f"
	daxpy_(nnodes, &r__, &w[(w_dim1 << 1) + 1], &c__1, &x[1], &c__1);
/*       do i=1,nnodes */
/*         print *,"x[",i,"] = ",x(i) */
/*       enddo */
#line 470 "./hierarchical.f"
	d__1 = -r__;
#line 470 "./hierarchical.f"
	daxpy_(nnodes, &d__1, &b[1], &c__1, &w[w_dim1 + 1], &c__1);
/*       do i=1,nnodes */
/*         print *,"r[",i,"] = ",w(i,1) */
/*       enddo */
#line 475 "./hierarchical.f"
	dif = r__ * dasum_(nnodes, &w[(w_dim1 << 1) + 1], &c__1);
#line 476 "./hierarchical.f"
	size = dasum_(nnodes, &x[1], &c__1);
/*       print *, "dif = ",dif */
/*       print *, "size = ",size */
#line 480 "./hierarchical.f"
	pre_(nnodes, &w[w_offset], &b[1]);
/*       do i=1,nnodes */
/*         print *,"P r[",i,"] = ",b(i) */
/*       enddo */
#line 485 "./hierarchical.f"
	t = ddot_(nnodes, &b[1], &c__1, &w[w_dim1 + 1], &c__1);
/*       print *, "r . P r = ",t */
#line 488 "./hierarchical.f"
	r__ = t / s;
#line 489 "./hierarchical.f"
	i__1 = *nnodes;
#line 489 "./hierarchical.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 490 "./hierarchical.f"
	    w[i__ + (w_dim1 << 1)] = b[i__] + r__ * w[i__ + (w_dim1 << 1)];
#line 491 "./hierarchical.f"
	}
/*       do i=1,nnodes */
/*         print *,"p[",i,"] = ",w(i,2) */
/*       enddo */
#line 496 "./hierarchical.f"
	s = t;
#line 497 "./hierarchical.f"
	stopit_(&dif, &size, ndigit, limit);
/*       print *, "dif = ",dif */
#line 499 "./hierarchical.f"
	++it;
#line 500 "./hierarchical.f"
    }
/*     do i=1,nnodes */
/*       print *,"x[",i,"] = ",x(i) */
/*     enddo */
/*     call mult(nelements,nnodes,order, */
/*    &  basis,basis_deriv,dirichlet,element_to_node,plobatto,rlobatto, */
/*    &    x, */
/*    &  b) */
/*     do i=1,nnodes */
/*       print *,"Ax[",i,"] = ",b(i) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 512 "./hierarchical.f"
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

#line 522 "./hierarchical.f"
    *dif = abs(*dif);
#line 523 "./hierarchical.f"
    *size = abs(*size);
/*     initialization during first iteration */
#line 525 "./hierarchical.f"
    if (i__ <= 0) {
#line 525 "./hierarchical.f"
	i__1 = -(*ndigit);
#line 525 "./hierarchical.f"
	t = pow_ri(&c_b22, &i__1);
#line 525 "./hierarchical.f"
    }
#line 526 "./hierarchical.f"
    ++i__;
#line 527 "./hierarchical.f"
    e = (doublereal) i__ * 3.;
/*     print *, "i,t,e = ",i,t,e */
#line 530 "./hierarchical.f"
    if (*dif <= t * *size) {
/*       stopping criterion I */
#line 532 "./hierarchical.f"
	e += 1.;
#line 533 "./hierarchical.f"
	*dif = -(*dif);
#line 534 "./hierarchical.f"
	i__ = 0;
#line 535 "./hierarchical.f"
    } else if (i__ >= *limit) {
/*       stopping criterion II */
#line 537 "./hierarchical.f"
	e += 2.;
#line 538 "./hierarchical.f"
	*dif = -(*dif);
#line 539 "./hierarchical.f"
	i__ = 0;
#line 540 "./hierarchical.f"
    }
#line 541 "./hierarchical.f"
    *size = e;
#line 542 "./hierarchical.f"
    return 0;
} /* stopit_ */

/* *********************************************************************** */
doublereal approximation_(nelements, nnodes, order, x, element_to_node__, 
	mesh, soln)
integer *nelements, *nnodes, *order;
doublereal *x;
integer *element_to_node__;
doublereal *mesh, *soln;
{
    /* System generated locals */
    integer element_to_node_dim1, element_to_node_offset, i__1;
    doublereal ret_val;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal legendre[11];
    static integer i__;
    static doublereal dx, xi, j2m1, basis;
    extern /* Subroutine */ int abort_();
    static integer element;

/*     integer j */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/*     do i=1,nelements */
/*       print *, "element_to_node[",i,"] = ", */
/*    &    (element_to_node(i,j),j=0,order) */
/*     enddo */
/*     do i=0,nelements */
/*       print *, "mesh[",i,"] = ",mesh(i) */
/*     enddo */
/*     do i=1,nnodes */
/*       print *,"soln[",i,"] = ",soln(i) */
/*     enddo */
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 572 "./hierarchical.f"
    /* Parameter adjustments */
#line 572 "./hierarchical.f"
    --soln;
#line 572 "./hierarchical.f"
    element_to_node_dim1 = *nelements;
#line 572 "./hierarchical.f"
    element_to_node_offset = 1 + element_to_node_dim1 * 0;
#line 572 "./hierarchical.f"
    element_to_node__ -= element_to_node_offset;
#line 572 "./hierarchical.f"

#line 572 "./hierarchical.f"
    /* Function Body */
#line 572 "./hierarchical.f"
    if (*order < 1 || *order > 10) {
#line 572 "./hierarchical.f"
	abort_();
#line 572 "./hierarchical.f"
    }
#line 574 "./hierarchical.f"
    dx = (float)1. / (doublereal) (*nelements);
#line 575 "./hierarchical.f"
    element = (integer) (*x / dx);
#line 576 "./hierarchical.f"
    if (element < *nelements) {
#line 576 "./hierarchical.f"
	++element;
#line 576 "./hierarchical.f"
    }
#line 577 "./hierarchical.f"
    xi = (*x - mesh[element - 1]) * 2. / dx - 1.;
#line 579 "./hierarchical.f"
    ret_val = ((1. - xi) * soln[element_to_node__[element]] + (xi + 1.) * 
	    soln[element_to_node__[element + element_to_node_dim1]]) * .5;
#line 581 "./hierarchical.f"
    if (*order > 1) {
#line 582 "./hierarchical.f"
	legendre[0] = 1.;
#line 583 "./hierarchical.f"
	legendre[1] = xi;
#line 584 "./hierarchical.f"
	i__1 = *order;
#line 584 "./hierarchical.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 585 "./hierarchical.f"
	    j2m1 = (doublereal) i__ * 2. - 1;
#line 586 "./hierarchical.f"
	    legendre[i__] = (j2m1 * legendre[i__ - 1] * xi - (doublereal) (
		    i__ - 1) * legendre[i__ - 2]) / (doublereal) i__;
#line 588 "./hierarchical.f"
	    basis = (legendre[i__] - legendre[i__ - 2]) / sqrt(j2m1 * 2.);
#line 589 "./hierarchical.f"
	    ret_val += basis * soln[element_to_node__[element + i__ * 
		    element_to_node_dim1]];
#line 591 "./hierarchical.f"
	}
#line 592 "./hierarchical.f"
    }
/*     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#line 594 "./hierarchical.f"
    return ret_val;
} /* approximation_ */

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
#line 606 "./hierarchical.f"
    ret_val = sin(machine_1.pi * *x);
/*     print *, "x = ",x */
/*     print *, "pi = ",pi */
/*     print *, "solution = ",solution */
#line 610 "./hierarchical.f"
    return ret_val;
} /* solution_ */

