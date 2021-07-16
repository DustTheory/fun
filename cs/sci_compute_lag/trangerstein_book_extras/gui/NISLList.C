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
// "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/NISLList.C,v 1.1 2009/08/20 17:33:32 johnt Exp $"
#include "NISLList.H"
#include "Errors.H"
#include "MemoryDebugger.H"
#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> ostream& operator<<(ostream& strm,
const NISLListNode<T> &node) {
  node.printOn(strm); return strm;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> ofstream& operator<<(ofstream& strm,
const NISLListNode<T> &node) {
  strm << *node.getData();
  return strm;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NISLList<T>::prepend(T *item) {
  if (&item != 0) {
    NISLListNode<T> *t = OPERATOR_NEW NISLListNode<T>(item);
    if (lastnode == 0) lastnode = t->selfConnect();
    else t->placeAfter(lastnode);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NISLList<T>::append(T *item) {
  if (&item != 0) {
    NISLListNode<T> *t=OPERATOR_NEW NISLListNode<T>(item);
    if (lastnode == 0) lastnode = t->selfConnect();
    else lastnode=t->placeAfter(lastnode);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NISLList<T>::insAfter(NISLListNode<T>* p,T *item)
{
  if (&item != 0) {
    NISLListNode<T> *t=OPERATOR_NEW NISLListNode<T>(item);
    if (lastnode == 0) lastnode = t->selfConnect();
    else if (p == 0) t->placeAfter(lastnode);
    else {
      t->placeAfter(p);
      if (p == lastnode) lastnode = t;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> T* NISLList<T>::delAfter(NISLListNode<T> *p) {
//test for empty list
  CHECK_POINTER(lastnode);
//make sure p not at end of list
  CHECK_TEST( p != lastnode );

//delAfter(0) means delete first
  if (!p) p = lastnode;      
  NISLListNode<T>* t = p->next();
  if (p == t) lastnode = 0;
  else if (lastnode == t) lastnode = p;
  p->disconnectNext();
  T *d = t->getData();
  delete t;
  return d;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NISLList<T>::join(NISLList<T> &b){
  NISLListNode<T>* t = b.lastnode;
  if (lastnode == 0) lastnode = t;
  else if (t != 0) {
    lastnode->switchConnectionWith(t);
    lastnode = t;
  }
  b.lastnode = 0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> int NISLList<T>::length() const {
  int l = 0;
  NISLListNode<T>* t = lastnode;
  if (t != 0) do { ++l; t = t->next(); } while (t != lastnode);
  return l;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> bool NISLList<T>::owns(NISLListNode<T>* p) const {
  NISLListNode<T>* t = lastnode;
  if (t != 0 && p != 0) {
    do {
      if (t == p) return true;
      t = t->next();
    } while (t != lastnode);
  }
  return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> bool NISLList<T>::ok() const {
  bool v = true;
  if (lastnode != 0) {
    NISLListNode<T>* t = lastnode;
    long count = 10000;      // Lots of chances to find last!
    do {
      count--;
      t = t->next();
    } while (count > 0 && t != lastnode);
    v &= count > 0;
  }
  CHECK_TEST(v);
  return v;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> void NISLList<T>::printOn( ostream& strm ) const {
  if (ok()) strm << *this;
  else strm << "list is not ok" << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> ostream& operator<<(ostream& strm,
const NISLList<T> &list) {
  NISLListNode<T>* p = list.first();
  if (p) {
    for( int l=0; p; p=list.next(p), l++ ){
      strm << l << " : " << *p << endl;
    }
  } else strm << "list is empty" << endl;
  return strm;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<class T> ofstream& operator<<(ofstream& strm,
const NISLList<T> &list) {
  strm << list.length();
  for( NISLListNode<T>* p=list.first(); p; p=list.next(p) ) strm << *p;
  return strm;
}
