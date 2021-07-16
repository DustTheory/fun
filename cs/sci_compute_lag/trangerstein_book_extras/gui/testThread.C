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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/testThread.C,v 1.1 2009/08/20 17:32:37 johnt Exp $"
#include "Thread.H"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "MemoryDebugger.H"
#include "Tracer.H"
#include "InputParameter.C"

template class InputParameter<bool>;
template class InputParameter<int>;
template class InputParameter<double>;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void a() {
  sched_yield();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void b() {  
#ifndef NO_THREAD
  Thread t2(a,"b calling a");
  t2.start();
  sched_yield();
//suspend b to give t2 a chance to run
  sleep(1);
#endif
}

int main(int argc,char **argv) {
  cout << boolalpha;
  {
#ifdef MEM_DEBUG
    MemoryDebugger md(1);
#endif

#ifndef NO_THREAD
    Thread *t0=OPERATOR_NEW Thread(a,"main calling a");
//  start is separate from the constructor so that we can control
//    when the processes start; useful for AMR
    t0->start();

    Thread *t1=OPERATOR_NEW Thread(b,"main calling b");
    t1->start();

//  halt the execution of the thread until resume is called
    t0->suspend();
//  mark this process as low priority
    t0->yield();

//  suspend main to give t1 a chance to run
    sleep(1);

    t0->resume();
  //t0->join();

//  suspend main to give t0 a chance to run
    sleep(1);

    t0->stop();
    delete t1;
    delete t0;
#endif
  }
  return EXIT_SUCCESS;
}
