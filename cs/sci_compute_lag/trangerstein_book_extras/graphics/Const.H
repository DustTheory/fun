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
#ifndef CONST_H
#define CONST_H "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/Const.H,v 1.1 2009/08/20 17:33:32 johnt Exp $"
#include <float.h>
#include <limits.h>

#define SMALL DBL_MIN
#define LARGE DBL_MAX
#define ROUNDOFF DBL_EPSILON
#define UNDEFINED HUGE_VAL

#define ZERO 0.e0
#define FOURTH 0.25e0
#define THIRD 0.33333333333333e0
#define HALF 0.5e0
#define ONE 1.e0
#define FOUR_THIRDS 1.33333333333333e0
#define THREE_HALVES 1.5e0
#define TWO 2.e0
#define THREE 3.e0
#define FOUR 4.e0

#define MAX_GRIDS 1
#define LENGTH_SUFFIX 6
#define LENGTH_OUTPUT_FILE_NAME_SUFFIX 2
#define LENGTH_HOST_SUFFIX 5
#define LENGTH_NAME 100
#define LENGTH_COMMENT 120
#define MAX_REPLICATION_FACTOR 12
#define COLOR_SIZE 256
#define MAX_COLORS 30
#define MAX_VECTORS 100

#define ARROW_SIZE 0.3e0
#endif
