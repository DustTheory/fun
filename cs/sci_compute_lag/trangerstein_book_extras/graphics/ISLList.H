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
#ifndef ISLLIST_H
#define ISLLIST_H "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/ISLList.H,v 1.1 2009/08/20 17:33:32 johnt Exp $"
//#include "Errors.H"
#include <stdlib.h>
using namespace std;
# ifndef ASSERT
# define ASSERT(s) { if (!(s)) { cerr << "assertion " << __STRING(s) << " failed in file " << __FILE__ << " at line " << __LINE__ << endl; abort(); } }
# endif
# ifndef OBSOLETE
# define OBSOLETE(s) { cerr << s << " in file " << __FILE__ << " at line " << __LINE__ << endl; abort(); }
# endif
# ifdef DEBUG
#   define CHECK_TEST(s) { ASSERT(s) }
#   define CHECK_POINTER(p) { ASSERT(p!=0) }
#   define CHECK_NONNEGATIVE(r) { ASSERT(r>=0) }
#   define CHECK_POSITIVE(r) { ASSERT(r>0) }
#   define CHECK_BOUNDS(v,l,u) { ASSERT(v>=l && v<u) }
#   define CHECK_NONZERO(q) { ASSERT(q!=0) }
#   define CHECK_SAME(v1,v2) { ASSERT((v1)==(v2)) }
# else
#   define CHECK_TEST(s) ;
#   define CHECK_POINTER(p) ;
#   define CHECK_NONNEGATIVE(r) ;
#   define CHECK_POSITIVE(r) ;
#   define CHECK_BOUNDS(v,l,u) ;
#   define CHECK_NONZERO(q) ;
#   define CHECK_SAME(v1,v2) ;
# endif

#include "ISLListNode.H"

//this will not work unless class T is publically derived from 
//  ISLListNode
template<class T> class ISLList {
  private:
    ISLListNode *lastnode;
  public:
    ISLList() : lastnode(0) {;}
    bool empty() const { return lastnode==0; }
    bool notEmpty() const { return lastnode != 0; }

    T* first() const {
      return static_cast<T*>((!lastnode) ? 0 : lastnode->after());
    }
    T* next(T *p) const {
      ISLListNode *pn=p;
      return ( !p ? first() : 
                    static_cast<T*>( pn==lastnode ? 0 : pn->after() ) );
    }
    T* last() const { return static_cast<T*>(lastnode); }

    void prepend(T &item) {
      if (&item != 0) {
        if (!lastnode) lastnode=item.selfConnect();
        else item.placeAfter(lastnode);
      }
    }
    void append(T &item) {
      if (&item != 0) {
        if (!lastnode) item.selfConnect();
        else item.placeAfter(lastnode);
        lastnode = &item;
      }
    }
    void join(ISLList &b) {
      ISLListNode *t=b.lastnode;
      if (!lastnode) lastnode=t;
      else if (t) {
        lastnode->switchConnectionsWith(t);
        lastnode=t;
      }
      b.lastnode = 0;
    }

    bool owns(const T*) const;
//  insAfter(0) means place first
    void insAfter(T *p,T &item) {
      if (!lastnode) lastnode=item.selfConnect();
      else if (!p) item.placeAfter(lastnode);
      else { item.placeAfter(p); if (p == lastnode) lastnode = &item; }
    }
//  delAfter(0) means delete first
    T* delAfter(T *p) {
      CHECK_POINTER(lastnode);
      CHECK_TEST(p != lastnode);
      if (!p) p = static_cast<T*>(lastnode);
      T *t = static_cast<T*>( p->after() );
      if (p == t) lastnode = 0;
      else if (lastnode == t) lastnode = p;
      return dynamic_cast<T*>(p->disconnectNext());
    }
    int length() const;
    bool ok() const;

//  a procedure to delete all members of the list has not been
//    provided, because it may be necessary to delete the SLListNode
//    data as well

    void printOn(ostream& os = cerr) const;
};
template<class T> ostream& operator<<(ostream&,const ISLList<T> &);
#endif
