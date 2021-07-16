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
// "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/ISLList.C,v 1.1 2009/08/20 17:33:32 johnt Exp $"
#include "ISLList.H"
#include <limits.h>
#include <stdlib.h>
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> int ISLList<T>::length() const {
  int l = 0;
  ISLListNode *t = lastnode;
  if (t) do { ++l; t = t->after(); } while (t != lastnode);
  return l;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> bool ISLList<T>::owns(const T *p) const {
  ISLListNode* t = lastnode;
  if (t != 0 && p != 0) {
    do {
      if (t == static_cast<const ISLListNode* >(p)) return true;
      t = t->after();
    } while (t != lastnode);
  }
  return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> bool ISLList<T>::ok() const {
  int v = true;
  if (lastnode) {
    ISLListNode* t = lastnode;
    long count = INT_MAX;
    do {
      count--;
      t = t->after();
    } while (count > 0 && t != lastnode);
    v &= count > 0;
  }
  CHECK_POINTER(v);
  return v;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void ISLList<T>::printOn(ostream& s) const {
  s << "ISLList: List with " << length() << " elements\n";
  if (ok()) {
    int l=0;
    for(T* p = first(); p; p=next(p), l++ ) {
      s << l << " : ";  p->printOn(s); s << "\n";
    }
  } else {
    s << "list is not ok" << endl;
  }
  s << "\n";
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> ostream& operator<<(ostream& s,const ISLList<T> &sl) {
  s << "List with " << sl.length() << " elements\n";
  T *p=sl.first();
  if (p) {
    for (int l=0; p; p=sl.next(p), l++ ) {
      s << l << " : ";  p->printOn(s); s << endl;
    }
  } else s << "list is empty" << endl;
  return s;
}
