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
#ifndef GUI_PLOT_ITEM_H
#define GUI_PLOT_ITEM_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIPlotItem.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "ISLList.H"
#include "VirtualInput.H"

class RealFlowVariable;
class GUIPlotItem : public ISLListNode {
  private:
    GUI_WIDGET_POINTER w;
    const RealFlowVariable* v;
  public:
    GUIPlotItem() : ISLListNode(),w(0),v(0) {;}
    GUIPlotItem(GUI_WIDGET_POINTER wp,const RealFlowVariable* vp) : 
      ISLListNode(),w(wp),v(vp) {}
    GUI_WIDGET_POINTER getWidget() { return w; }
    const RealFlowVariable* getVar() { return v; }
    virtual void printOn(ostream &os) const {
      os << "GUIPlotItem:" << endl;
      ISLListNode::printOn(os);
    }
};
#endif
