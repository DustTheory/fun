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
#ifndef GUI_BASE_H
#define GUI_BASE_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIBase.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
using namespace std;
#include <iostream>
#include "ClassThread.H"

class GUIInputParameterList;
class GUIBase {
  protected :
    char* display_name;
    void (*runMain)(bool);
    void (*checkMainInput)();
    void (*cleanup)();
    void (*shutdown)();
#ifndef NO_THREAD
    ClassThread<GUIBase> *run_main_thread;
#endif
    bool run_main_called,quit_called;
    bool use_thread;

    void runMainProc() { 
#ifndef NO_THREAD
      if (use_thread && run_main_thread!=0) run_main_thread->start(); 
      else 
#endif
      runMain(run_main_called);
      run_main_called=true; 
    }
    void runMainAndPause() const { 
      runMain(run_main_called); 
    }
  protected:
    void suspendRun();
    void resumeRun();
    bool mainRunning() const { 
#ifndef NO_THREAD
      if (use_thread && run_main_thread!=0) {
        return !run_main_thread->finished();
      }
#endif
      return false;
    }
    bool mainFinished() const { 
#ifndef NO_THREAD
      if (use_thread && run_main_thread!=0) {
        return run_main_thread->finished();
      }
#endif
      return true;
    }
   
    virtual void hideBusyDialog() const=0;
    virtual void showWarningDialog(const char* str)=0;
  public:
    GUIBase(char*,void (*)(bool),void (*)()=0,void (*)()=0,
      void (*)()=0,bool=false);
    virtual ~GUIBase();
    virtual void createMainWindow(int,char*[])=0;
    virtual void createFileMenu()=0;
    virtual void eventLoop()=0;
    void createHelpMenu() {;}
    virtual void printOn(ostream &os = cerr) const;
};
#endif
