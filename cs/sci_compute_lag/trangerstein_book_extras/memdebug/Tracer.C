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
// "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/Tracer.C,v 1.1 2009/08/20 17:33:33 johnt Exp $"
#include "Tracer.H"
#include <stdio.h>
#include <string.h>
#include "Debug.H"
#include "MemoryDebugger.H"

#define INTEGER_STRING_LENGTH 5
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Tracer::Tracer(const char *s) {
  int length=strlen(s)+1;
  str=OPERATOR_NEW_BRACKET(char,length);
  strcpy(str,s);
  cout << "entering " << str << endl;
  DEBUG_ONLY(mem_check();)
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Tracer::Tracer(const char *s,const char *f,int l) {
  int length=strlen(s)+strlen(f)+18+INTEGER_STRING_LENGTH+1;
  str=OPERATOR_NEW_BRACKET(char,length);
  snprintf(str,length,"%s in file %s at line %d",s,f,l);
  cout << "entering " << str << endl; 
  DEBUG_ONLY(mem_check();)
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Tracer::Tracer(const char *s,int v,const char *f,int l) {
  int length=strlen(s)+strlen(f)+18+2*INTEGER_STRING_LENGTH+1;
  str=OPERATOR_NEW_BRACKET(char,length);
  snprintf(str,length,"%s%d in file %s at line %d",s,v,f,l);
  cout << "entering " << str << endl;
  DEBUG_ONLY(mem_check();)
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Tracer::~Tracer()  {
  cout << "leaving " << str << endl;
  delete [] str;
  DEBUG_ONLY(mem_check();)
}
