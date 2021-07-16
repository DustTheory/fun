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
#ifndef _ARCH_H_
#define _ARCH_H_ "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/Arch.H,v 1.1 2009/08/20 17:33:32 johnt Exp $" 

//using namespace std;

//#if defined (__digital__) && defined (__unix__)
//#define F77NAME(x) x ## _
//#define F77_NAME(x) x ## _
//#else
#define F77NAME(x) x ## _
#define F77_NAME(x) x ## __
//#endif
#endif
