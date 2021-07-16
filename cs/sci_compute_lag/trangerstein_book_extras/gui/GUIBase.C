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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIBase.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "GUIBase.H"
#include "ClassThread.H"
#include "MemoryDebugger.H"
#include "Tracer.H"

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GUIBase::GUIBase(char* dn,void (*m)(bool),void (*cm)(),void (*cu)(),
void (*sd)(),bool ut) : display_name(dn),
runMain(m),checkMainInput(cm),cleanup(cu),shutdown(sd),
#ifndef NO_THREAD
run_main_thread(0),
#endif
run_main_called(false),quit_called(false),use_thread(ut) {
#ifdef NO_THREAD
  use_thread=false;
#else
  if (use_thread) {
    run_main_thread=OPERATOR_NEW ClassThread<GUIBase>(this,
      &GUIBase::runMainAndPause,const_cast<char*>("runMainAndPause"));
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GUIBase::~GUIBase() {
#ifndef NO_THREAD
  if (run_main_thread) delete run_main_thread; 
  run_main_thread=0;
#endif
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUIBase::printOn(ostream &os) const {
  os << "GUIBase: display_name = " << display_name << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUIBase::suspendRun() {
#ifndef NO_THREAD
  CHECK_POINTER(run_main_thread)
  run_main_thread->suspend();
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUIBase::resumeRun() {
#ifndef NO_THREAD
  CHECK_POINTER(run_main_thread)
  run_main_thread->resume();
#endif
}

#ifndef NO_THREAD
template class ClassThread<GUIBase>;
#endif
