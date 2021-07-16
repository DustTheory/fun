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
#ifndef TYPED_PLOT_OBJ_H
#define TYPED_PLOT_OBJ_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/TypedPlotObj.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#if (SPACEDIM>1)

#include <iostream>
#include <X11/Intrinsic.h>
#include <GL/glx.h>
#include "VirtualGLWindow.H"

class VirtualVolGraphTool;

template<class T> class TypedPlotObj : public VirtualGLWindow::PlotObj {
  private:
    T *dataset;
    VirtualVolGraphTool *graph_tool;
  public:
//  t must be constructed via operator new:
    TypedPlotObj(T *t) : dataset(t),graph_tool(t->getGraphTool()) { 
      CHECK_POINTER(dataset)
    }
    ~TypedPlotObj() { /* if (dataset!=0) delete dataset; dataset=0; */ }
    virtual bool drawSurface() const {return dataset->drawSurface();}
    virtual void plot(bool do_best) const { dataset->plot(do_best); }
# if (SPACEDIM==3)
    virtual void plot(bool do_best,const AxisClipPlane *cp) const {
      dataset->plot(do_best,cp);
    }
# endif
    virtual void printOn(ostream &os=cout) const {
      os << "TypedPlotObj: " << this << endl;
      os << "\tdataset = " << dataset << endl;
      dataset->printOn(os);
    }
};
#endif
#endif
