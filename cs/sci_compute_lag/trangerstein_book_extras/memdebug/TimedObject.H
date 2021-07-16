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
#ifndef TIMEDOBJECT_H
#define TIMEDOBJECT_H "$Header:$"
#include <stdlib.h>
#include <sys/times.h>
#include <sys/param.h>
#include "Errors.H"

class TimedObject {
  private:
  protected:
    char *name;
    int number_calls;
    bool is_on;
    double total_run_time;
    clock_t start;
    struct tms usage;

    TimedObject( const TimedObject& ) { 
      OBSOLETE("copy constructor disallowed"); 
    }
  public:
    TimedObject(const char *);
    virtual ~TimedObject() { if (name) delete name; name=0; }

    void turnOn();
    double turnOff();
    bool isOn() const { return is_on; }
    int numberCalls() const { return number_calls; }
    double totalRunTime() const { return total_run_time; }
    virtual void reset() {
      ASSERT(!is_on);
      number_calls=0; is_on=false; total_run_time=0.; start=0;
    }
    virtual void printOn(ostream&) const;
};

class Timer {
  private:
  protected:
    TimedObject *timed_object;
  public:
    Timer(TimedObject *t) : timed_object(t) { timed_object->turnOn(); }
    virtual ~Timer() { timed_object->turnOff(); }
    virtual void printOn(ostream &os) const {
      timed_object->printOn(os);
    }
};
#endif
