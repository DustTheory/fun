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
// "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/TimedObject.C,v 1.1 2009/08/20 17:33:33 johnt Exp $"
#include "TimedObject.H"
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <string.h>
#include "MemoryDebugger.H"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
TimedObject::TimedObject(const char *na) : number_calls(0),is_on(false),
total_run_time(0.),start(0) {
  int length_name=strlen(na);
  name=OPERATOR_NEW_BRACKET(char,length_name+1);
  strcpy(name,na);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
TimedObject::TimedObject(const char *na,TimedObjectList &tol) :
number_calls(0),is_on(false),total_run_time(0.),start(0) {
  int length_name=strlen(na);
  name=OPERATOR_NEW_BRACKET(char,length_name+1);
  strcpy(name,na);
  tol.append(*this);
}
*/
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void TimedObject::turnOn() {
  ASSERT(!is_on);
  is_on=true;
  number_calls++;
  times(&usage);
  start=usage.tms_utime;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double TimedObject::turnOff() {
  if (!is_on) return 0.;
  times(&usage);
  double elapsed=static_cast<double>(usage.tms_utime-start)
                /static_cast<double>(sysconf(_SC_CLK_TCK));
  total_run_time += elapsed;
  is_on=false;
  return elapsed;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void TimedObject::printOn(ostream &os) const {
  os << "TimedObject: name = " << name << endl;
  if (number_calls>0) {
    os << "\ttotal_run_time = " << total_run_time
       << " number_calls = " << number_calls
       << "\n\t\ttime/call = "
       << total_run_time/static_cast<double>(number_calls) 
       << endl;
  } else {
    os << "\ttotal_run_time = " << total_run_time
       << " number_calls = " << number_calls << endl;
  }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
/*
TimedObjectList::~TimedObjectList() {
  while (notEmpty()) delete delAfter(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void TimedObjectList::reset() const {
  for (TimedObject *p=first();p;p=next(p)) p->reset();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void TimedObjectList::printOn(ostream &os) const {
  for (TimedObject *p=first();p;p=next(p)) {
    if (p->numberCalls()>0) p->printOn(os);
  }
}
*/
