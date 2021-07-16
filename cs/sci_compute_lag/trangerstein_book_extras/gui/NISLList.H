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
#ifndef NISLLIST_H
#define NISLLIST_H "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/NISLList.H,v 1.1 2009/08/20 17:33:32 johnt Exp $"
#include <fstream>
#include <iostream>
#include "Types.H"

template<class T> class NISLList;
template<class T> class NISLListNode {
  private:
    NISLListNode *nextnode;
    T *data;
  public:
    inline NISLListNode(T *item) : nextnode(0), data(item) {;}
    inline virtual ~NISLListNode() {;}
    inline T* getData() { return data; }
    inline const T* getData() const { return data; }
    inline NISLListNode* next() const { return nextnode; }
    inline NISLListNode* selfConnect() { return nextnode=this; }
    inline NISLListNode* placeAfter(NISLListNode *n) {
      nextnode=n->nextnode; n->nextnode=this; return this;
    }
    inline NISLListNode* disconnectNext() { 
      NISLListNode *t=nextnode; nextnode=t->nextnode; t->nextnode=0;
      return t;
    }
    inline void switchConnectionWith(NISLListNode* n) {
      NISLListNode *t=nextnode; nextnode=n->nextnode; n->nextnode=t;
    }

    virtual void printOn(ostream &os= cerr) const {
      os << "NISLListNode: this = " << this
         << " nextnode = " << nextnode 
         << " data = " << reinterpret_cast<void*>(data);
    }
};
template<class T> ostream& operator<<(ostream&,const NISLListNode<T>&);
template<class T> ofstream& operator<<(ofstream&,const NISLListNode<T>&);

template<class T> class NISLList {
  private:
    NISLListNode<T> *lastnode;
  public:
    inline NISLList() : lastnode(0) {;}
    inline virtual ~NISLList() {;}
    inline bool empty() const { return lastnode==0; }
    inline bool notEmpty() const { return lastnode != 0; }

    inline NISLListNode<T>* first() const {
      return (lastnode==0) ? static_cast<NISLListNode<T>* >(0) 
                              : lastnode->next();
    }
    inline NISLListNode<T>* next(const NISLListNode<T> *p) const {
      return (p==0 || p==lastnode) ? static_cast<NISLListNode<T>* >(0) 
                                      : p->next();
    }
    inline NISLListNode<T>* last() const { return lastnode; }

    void prepend(T*);
    void append(T*);
    void join(NISLList&); 

    bool owns(NISLListNode<T>*) const;
    void insAfter(NISLListNode<T>*,T*);
    T* delAfter(NISLListNode<T>*);
    int length() const;
    bool ok() const;
    virtual void printOn(ostream& = cerr) const ;

//  a procedure to delete all members of the list has not been
//    provided, because it may be necessary to delete the T
//    data as well

};
template<class T> ofstream& operator<<(ofstream&,const NISLList<T>&);
template<class T> ostream& operator<<(ostream&,const NISLList<T>&);
#endif
