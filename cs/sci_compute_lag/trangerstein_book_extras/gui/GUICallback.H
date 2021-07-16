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
#ifndef GUI_CALLBACK_H
#define GUI_CALLBACK_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUICallback.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "ISLList.H"

#ifdef USE_GTK
#define GUI_TYPE GTKGUI
#else
#define GUI_TYPE GUI
#endif

class GUI_TYPE;
template<class T> class GUICallbackStructure : public ISLListNode {
  private:
    T *data;
    GUI_TYPE *gui;
  public:
    GUICallbackStructure(T *d,GUI_TYPE *g) : data(d),gui(g) {;}
    GUI_TYPE *getGUI() const { return gui; }
    T *getData() const { return data; }
    virtual void printOn(ostream &os) const {
      os << "GUICallbackStructure:" << endl;
    }
};

template<class T> class GUICallbackStructureList : 
public ISLList< GUICallbackStructure<T> > {
  public:
    GUICallbackStructureList() {;}
    virtual ~GUICallbackStructureList() {
      while (notEmpty()) delete delAfter(0);
    }
    bool notEmpty() const { 
      return ISLList< GUICallbackStructure<T> >::notEmpty(); 
    }
    GUICallbackStructure<T>* delAfter(GUICallbackStructure<T> *p) {
      return ISLList< GUICallbackStructure<T> >::delAfter(p);
    }
    virtual void printOn(ostream &os) const {
      os << "GUICallbackStructureList:" << endl;
      ISLList< GUICallbackStructure<T> >::printOn(os);
    }
};
#endif
