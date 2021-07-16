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
#ifndef CLASS_THREAD_H
#define CLASS_THREAD_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/ClassThread.H,v 1.1 2009/08/20 17:32:34 johnt Exp $"
#ifndef NO_THREAD
#include "Thread.H"

template<class T> class Thread1 : public VirtualThread {
  private:
    void (*exec)(T);
    T data;

    virtual void run() const { if (exec) (*exec)(data); }
  public:
    Thread1(void (*f)(T),T t,char *tname=0) : VirtualThread(tname),
      exec(f),data(t) {;}
    virtual ~Thread1() {;}
    virtual void printOn(ostream &os) const {
      os << "\tThread1:" << endl;
      VirtualThread::printOn(os);
    }
};

template<class C> class ClassThread : public VirtualThread {
  private:
    void (C::*exec)() const;
    C *c_this;

    virtual void run() const { if (exec) (c_this->*exec)(); }
  public:
    ClassThread(C *ct,void (C::*f)() const,char *n=0) : 
      VirtualThread(n), c_this(ct),exec(f) {;}
    virtual ~ClassThread() {;}
    virtual void printOn(ostream &os) const {
      os << "\tClassThread:" << endl;
      VirtualThread::printOn(os);
    }
};

template<class C,class T> class ClassThread1 : public VirtualThread {
  private:
    void (C::*exec)(T) const;
    C *c_this;
    T data;

    virtual void run() const { if (exec) (c_this->*exec)(data); }
  public:
    ClassThread1(C *ct,void (C::*f)(T) const,T t,char *n=0) : 
      VirtualThread(n),c_this(ct),exec(f),data(t) {;}
    virtual ~ClassThread1() {;}
    virtual void printOn(ostream &os) const {
      os << "\tClassThread1:" << endl;
      VirtualThread::printOn(os);
    }
};
#endif
#endif
