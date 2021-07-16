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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/Thread.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#ifndef NO_THREAD
#include <signal.h>
#include <string.h>
#include "Thread.H"
#include "MemoryDebugger.H"
#include "Tracer.H"

#include <stdlib.h>

#ifdef __linux__
#define PTHREAD_KILL    21
#define PTHREAD_SUSPEND 22
#else
#define PTHREAD_KILL    SIGUSR1
#define PTHREAD_SUSPEND SIGUSR2
#endif

pthread_mutex_t VirtualThread::mutex=PTHREAD_MUTEX_INITIALIZER;
pthread_key_t VirtualThread::key;
pthread_attr_t VirtualThread::attr;
int VirtualThread::nthreads=0;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void pthread_signal_handler(int signum) {
  switch(signum) {
    case PTHREAD_SUSPEND:
      VirtualThread::keySuspend();
      break;
    case PTHREAD_KILL:
      VirtualThread::keyKill();
      break;
    case SIGTERM:
    case SIGINT:
    case SIGQUIT:
    case SIGABRT:
#ifdef DEBUG
      cout << "Threads received abort signal. Bye!"<<endl;
#endif    
      pthread_exit(0);
      break;
    default: 
#ifdef DEBUG
      cout << "Threads received unknown signal. Bye!"<<endl;
#endif    
      pthread_exit(0);
      break;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualThread::VirtualThread(char *tname) : name(0),thread(0),
suspend_flag(FALSE),run_status(NOT_STARTED),joinable(FALSE) {
  pthread_mutex_init(&suspend_mutex,0);
  pthread_cond_init(&suspend_cond,0);

  if (tname) {
    name=OPERATOR_NEW_BRACKET(char,strlen(tname)+1);
    strcpy(name,tname);
  }

  gblock(); nthreads++; gunblock();

  if (nthreads == 1) { // could do pthread_once() for this
    gblock();
    if (nthreads == 1) {
      struct sigaction new_action;
      
      new_action.sa_handler = pthread_signal_handler;
      new_action.sa_flags = 0;
#ifdef __linux__
      sigdelset(&(new_action.sa_mask),SIGTERM);
      sigdelset(&(new_action.sa_mask),SIGINT);
      sigdelset(&(new_action.sa_mask),SIGABRT);
      sigdelset(&(new_action.sa_mask),SIGQUIT);
#endif
      sigaction(PTHREAD_SUSPEND, &new_action,0);
      sigaction(PTHREAD_KILL, &new_action,0);

      sigaction(SIGTERM, &new_action,0);
      sigaction(SIGINT, &new_action,0);
      sigaction(SIGABRT, &new_action,0);
      sigaction(SIGQUIT, &new_action,0);

      pthread_key_create(&key,0);
      
      pthread_attr_init(&attr);
      pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
    }
    gunblock();
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualThread::~VirtualThread() {
  if (thread) {
    stop();
    pthread_detach(thread);
    pthread_cancel(thread);
    pthread_mutex_destroy(&suspend_mutex);
    pthread_cond_destroy(&suspend_cond);
    thread=0;
  }
  gblock(); nthreads--; gunblock();
  if (nthreads == 0) {
    gblock();
    if (nthreads == 0) {
      pthread_mutex_destroy(&mutex);
      pthread_key_delete(key);
      pthread_attr_destroy(&attr);
      nthreads=0;
    }
    gunblock();
  }
  if (name) delete [] name; name=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualThread::start(bool j){
  joinable=j; 
  pthread_attr_t *__attr= (joinable ? 0 : &attr);
  void* __arg=reinterpret_cast<void*>(this);
  ASSERT(!pthread_create(&thread,__attr,staticRun,__arg));
  while (!isAlive() && !finished()) {;}
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualThread::stop() { 
  if (isSelf()) {
    pthread_signal_handler(PTHREAD_KILL);
  } else if (thread) {
    pthread_detach(thread);
    pthread_kill(thread,PTHREAD_KILL); 
  }
  gblock(); run_status=FINISHED; gunblock();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualThread::suspend() {
  gblock();
  if (run_status!=RUNNING){
    block(); suspend_flag=TRUE; unblock();
    gunblock();
    return;
  }
  gunblock();
  
  block();
  if (suspend_flag) {
    unblock();
  } else {
    if (isSelf()) { // self-suspend
      unblock();
      pthread_signal_handler(PTHREAD_SUSPEND);
    } else {
      pthread_kill(thread, PTHREAD_SUSPEND);
      cvWait(); 
      yield();
      unblock();
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualThread::resume() {
  block();
  if (!suspend_flag) {
    unblock();
    return;
  }
  suspend_flag=FALSE;
  signal();
  unblock();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualThread::printOn(ostream &os) const {
  os << "VirtualThread: name = " << name 
     << ", this = " << this << ", tid = " << tid() << endl;
  os << "\tthread = " << thread << endl;
//os << "\tsuspend_mutex = " << &suspend_mutex << endl;
//os << "\tsuspend_cond = " << &suspend_cond << endl;
  os << "\tsuspend_flag = " << suspend_flag << endl;
//os << "\trun_status = " << run_status << endl;
//os << "\n\tcond = " << &cond << endl;
//os << "\tmutex = " << &mutex << endl;
//os << "\tonce = " << once << endl;
  os << "\tkey = " << key << endl;
  os << "\tkeyThread = " << keyThread() << endl;
//os << "\tsignal_thread = " << signal_thread << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Thread::printOn(ostream &os) const {
  os << "Thread:" << endl;
  VirtualThread::printOn(os);
}
#endif
