/***********************************************************************
*  Copyright 2006 John A. Trangenstein
*
*  This software is made available for research and instructional use 
*  only. 
*  You may copy and use this software without charge for these 
*  non-commercial purposes, provided that the copyright notice and 
*  associated text is reproduced on all copies.  
*  For all other uses (including distribution of modified versions), 
*  please contact the author at
*    John A. Trangenstein
*    Department of Mathematics
*    Duke University
*    Durham, NC 27708-0320
*    USA
*  or
*    johnt@math.duke.edu
*  
*  This software is made available "as is" without any assurance that it
*  is completely correct, or that it will work for your purposes.  
*  Use the software at your own risk.
***********************************************************************/
/* "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/second.c,v 1.1 2009/08/20 17:33:34 johnt Exp $" */
#include <sys/times.h>
#include <unistd.h>

void second(t) double *t;
{
  struct tms buffer;

  times(&buffer);
  *t = (buffer.tms_utime + buffer.tms_stime)
     /(1.0*sysconf(_SC_CLK_TCK));
}

/* for Fortran calls */
void second_(t) double *t;
{
  second(t);
}
