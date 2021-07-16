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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/InputParameter.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "InputParameter.H"
#include <stdio.h>
#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const char* InputParameter<bool>::getValue() const {
  snprintf(vstring,LENGTH_NAME,"%s",(*ptr_data?"true":"false"));
  return vstring;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const char* InputParameter<int>::getValue() const {
  snprintf(vstring,LENGTH_NAME,"%d",*ptr_data);
  return vstring;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> const char* InputParameter<double>::getValue() const {
  snprintf(vstring,LENGTH_NAME,"%g",*ptr_data);
  return vstring;
}
