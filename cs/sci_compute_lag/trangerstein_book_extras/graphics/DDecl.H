//**********************************************************************
// Copyright 2006 John A. Trangenstein
//
// This software is made available for research and instructional use 
// only. 
// You may copy and use this software without charge for these 
// non-commercial purposes, provided that the copyright notice and 
// associated text is reproduced on all copies.  
// For all other uses (including distribution of modified versions), 
// please contact the author at
//   John A. Trangenstein
//   Department of Mathematics
//   Duke University
//   Durham, NC 27708-0320
//   USA
// or
//   johnt@math.duke.edu
// 
// This software is made available "as is" without any assurance that it
// is completely correct, or that it will work for your purposes.  
// Use the software at your own risk.
//**********************************************************************
#ifndef DDECL_H
#define DDECL_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/DDecl.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"

//examples:
//  for surly fortran compilers,
//    fort_(D_EXPR(fi[0],fi[1],fi[2]),...)
//  D_EXPR(vect[0] *= s, vect[1] *= s, vect[2] *= s);
//  IntVect(D_DECL(p[0] + s, p[1] + s, p[2] + s));
//  return D_TERM(len[0], *len[1], *len[2]);
#if (SPACEDIM==1)
# define D_EXPR(a,b,c) ((void)((a),0))
# define D_DECL(a,b,c) a
# define D_TERM(a,b,c) a
# define D_ASSIGN(a,b,c) a;
# define D_CHOOSE(a,b,c) a
# define DECL_ARG(d) const int&
# define DECL_DARG(d) const double*
#endif
#if (SPACEDIM==2)
# define D_EXPR(a,b,c) ((void)((a),(b),0))
# define D_DECL(a,b,c) a,b
# define D_TERM(a,b,c) a b
# define D_ASSIGN(a,b,c) a;b;
# define D_CHOOSE(a,b,c) b
# define DECL_ARG(d) const int&,const int&
# define DECL_DARG(d) const double*,const double*
#endif
#if (SPACEDIM==3)
# define D_EXPR(a,b,c) ((void)((a),(b),(c),0))
# define D_DECL(a,b,c) a,b,c
# define D_TERM(a,b,c) a b c
# define D_ASSIGN(a,b,c) a;b;c;
# define D_CHOOSE(a,b,c) c
# define DECL_ARG(d) const int&,const int&,const int&
# define DECL_DARG(d) const double*,const double*,const double*
#endif
# define DIM_ARG(d) D_DECL(d[0],d[1],d[2])

#endif
