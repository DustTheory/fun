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
#ifndef ISLLIST_NODE_H
#define ISLLIST_NODE_H "$Header: /home/faculty/johnt/cvs/deal_new/memdebug/ISLListNode.H,v 1.1 2009/08/20 17:33:32 johnt Exp $"
#include <iostream>

using namespace std;

class ISLListNode {
  private:
    ISLListNode *nextnode;
  public:
    inline ISLListNode() : nextnode(0) {;}
    inline virtual ~ISLListNode() {;}
    inline ISLListNode* after() const { return nextnode; }
    inline ISLListNode* selfConnect() { 
      return nextnode=this;
    }
    inline void placeAfter(ISLListNode *n) {
      nextnode=n->nextnode;
      n->nextnode=this;
    }
    inline void disconnect() { nextnode=0; }
    inline ISLListNode* disconnectNext() {
      ISLListNode *n=nextnode;
      nextnode=n->nextnode;
      n->disconnect();
      return n;
    }
//if this and n are in different lists, joins the lists
//if this and n are in same list, separates the lists
    inline void switchConnectionsWith(ISLListNode *n) {
      ISLListNode *t=nextnode;
      nextnode=n->nextnode;
      n->nextnode=t;
    }
    virtual void printOn(ostream &os = cerr) const {
      os << "ISLListNode: this = " << this 
         << " nextnode = " << nextnode;
    }
};
inline ostream& operator<<(ostream &os,const ISLListNode &node) {
  node.printOn(os); return os;
}
#endif
