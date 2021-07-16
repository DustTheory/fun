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
#ifndef TRACER_H
#define TRACER_H "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/Tracer.H,v 1.1 2009/08/20 17:33:33 johnt Exp $"

#include "Errors.H"

//using namespace std;

class Tracer {
  private:
    char *str;
  public:
    Tracer(const char*);
    Tracer(const char*,const char*,int);
    Tracer(const char*,int,const char*,int);
    ~Tracer();
};

#define TRACER_CALL(n,s) Tracer n(s,__FILE__,__LINE__);
#define TRACER_CALLV(n,s,v) Tracer n(s,v,__FILE__,__LINE__);
#endif
