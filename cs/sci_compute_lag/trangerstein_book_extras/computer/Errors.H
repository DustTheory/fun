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
#ifndef ERRORS_H
#define ERRORS_H "$Header:$"
#include <cstdlib>
#include <iostream>

using namespace std;

//assert no longer dumps core, so we define our own macro
# define ASSERT(s) { if (!(s)) { cerr << "assertion " << __STRING(s) << " failed in file " << __FILE__ << " at line " << __LINE__ << endl; abort(); } }
# define OBSOLETE(s) { cerr << s << " in file " << __FILE__ << " at line " << __LINE__ << endl; abort(); }
# ifdef DEBUG
//  if we throw but don't catch, then all destructors are called,
//    and the debugger traces back to the Memory Debugger
//    THANK YOU, BJARNE
//  don't use throw for copy constructors: Stroustrup p 373 
#   define CHECK_TEST(s) { ASSERT(s) }
#   define CHECK_POINTER(p) { ASSERT(p!=0) }
#   define CHECK_DIRECTION(d) { ASSERT(d>=0 && d<SPACEDIM) } 
#   define CHECK_SAME_DIRECTIONS(d1,d2) { ASSERT(d1==d2) }
#   define CHECK_DIFFERENT_DIRECTIONS(d1,d2) { ASSERT(d1!=d2) }
#   define CHECK_NONNEGATIVE(r) { ASSERT(r>=0) }
#   define CHECK_POSITIVE(r) { ASSERT(r>0) }
#   define CHECK_BOUNDS(v,l,u) { ASSERT(v>=l && v<u) }
#   define CHECK_NONZERO(q) { ASSERT(q!=0) }
#   define CHECK_SAME(v1,v2) { ASSERT((v1)==(v2)) }
#   define CHECK_CONTAINS(T,b,v) { ASSERT(b.contains(v)) }
#   define CHECK_STRICTLY_CONTAINS(T,b,v) { ASSERT((b).strictlyContains(v)) }
# else
#   define CHECK_TEST(s) ;
#   define CHECK_POINTER(p) ;
#   define CHECK_DIRECTION(d) ;
#   define CHECK_SAME_DIRECTIONS(d1,d2) ;
#   define CHECK_DIFFERENT_DIRECTIONS(d1,d2) ;
#   define CHECK_NONNEGATIVE(r) ;
#   define CHECK_POSITIVE(r) ;
#   define CHECK_BOUNDS(v,l,u) ;
#   define CHECK_NONZERO(q) ;
#   define CHECK_SAME(v1,v2) ;
#   define CHECK_CONTAINS(T,b,v) ;
#   define CHECK_STRICTLY_CONTAINS(T,b,v) ;
# endif
#endif
