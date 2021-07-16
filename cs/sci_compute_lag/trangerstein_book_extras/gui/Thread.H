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
#ifndef THREAD_H
#define THREAD_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/Thread.H,v 1.1 2009/08/20 17:32:36 johnt Exp $"

#ifndef NO_THREAD
#include <pthread.h>
#include <sched.h>
#include <iostream>
#include "Errors.H"

class VirtualThread {
  friend void pthread_signal_handler(int);
  private:
    enum RUN_STATUS { NOT_STARTED, RUNNING, FINISHED };

    static pthread_mutex_t mutex; //shared among all threads
    static pthread_key_t key;
    static pthread_attr_t attr;
    static int nthreads;

    char *name;
    pthread_t thread;
    pthread_mutex_t suspend_mutex;
    pthread_cond_t  suspend_cond;
    volatile bool suspend_flag;
    RUN_STATUS run_status;
    bool joinable;

    static void* staticRun(void *thrd) {
      VirtualThread* t=reinterpret_cast<VirtualThread*>(thrd);
      CHECK_POINTER(t);
      t->baseRun();
      return static_cast<void*>(0);
    }
    static VirtualThread* keyThread() {
      return reinterpret_cast<VirtualThread*>(pthread_getspecific(key));
    }

    bool isSelf() const {
      return pthread_equal(thread,pthread_self());
    }
    void block() { pthread_mutex_lock(&suspend_mutex); }
    void unblock() { pthread_mutex_unlock(&suspend_mutex); }
    void gblock() { pthread_mutex_lock(&mutex); }
    void gunblock() { pthread_mutex_unlock(&mutex); }
    pthread_t tid() const { return thread; }
    void baseRun() { 
      pthread_setspecific(VirtualThread::key,this);
      gblock(); run_status=RUNNING; gunblock();
      block();
      if (suspend_flag) {
        unblock();
        suspend();
      } else {
        unblock();
      }
      run(); 
      gblock(); 
      run_status=FINISHED; thread=0; gunblock();
    }

    static void keyKill() { // kill calling thread
      VirtualThread *k=0;
      k=keyThread();
      if (k) {
        k->block(); k->suspend_flag=FALSE; k->unblock();
        k->gblock(); k->run_status=FINISHED; k->gunblock();
        pthread_detach(k->tid());
      }
      pthread_exit(0);
    }

    static void keySuspend() {
      VirtualThread* k=keyThread();
      k->block(); k->suspend_flag=TRUE; k->unblock();
      k->signal();
      k->block();
      while (k->suspend_flag) k->cvWait();
      k->unblock();
    }
  public:
    VirtualThread(char* =0);
    virtual ~VirtualThread();

    bool isAlive() const { return run_status==RUNNING; }// gblock()?
    bool finished() const { return run_status==FINISHED; }
    bool isSuspended() const { return suspend_flag; }
    void resetStatus() { if (finished()) run_status=NOT_STARTED; }
    char* getName() const { return name; }
    static void yield() { 
      sched_yield(); 
    }
    void signal() { 
      pthread_cond_signal(&suspend_cond); 
    }
    void cvWait() { 
      pthread_cond_wait(&suspend_cond, &suspend_mutex);
    }

    void start(bool=FALSE); 
    void suspend();
    void resume(); // valid only after a suspend() call

//  wait for the specified joinable thread to end
    void join() { 
      pthread_join(thread,0); 
    }
    void stop(); //kill the thread
    virtual void run() const = 0;

    virtual void printOn(ostream&) const;
};

class Thread : public VirtualThread {
  private:
    void (*exec)();

    virtual void run() const { if (exec) exec(); }
  public:
    Thread(void (*f)(),char *tname=0) : VirtualThread(tname),exec(f) {;}
    virtual ~Thread() {;}
    virtual void printOn(ostream&) const;
};
#endif
#endif
