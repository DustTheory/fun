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
#ifndef XY_GRAPH_TOOL_H
#define XY_GRAPH_TOOL_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/XYGraphTool.H,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#include <complex>
#include <stdio.h>
#ifdef USE_GTK
#include "GTKWindow.H"
#else
#include "XTWindow.H"
#endif
#include "VGT.H"

class VirtualColormap;
class XYGraphTool : public VirtualGraphTool12 {
  public:
#ifdef USE_GTK
    typedef GTKWindow WINDOW_TYPE;
#else
    typedef XTWindow WINDOW_TYPE;
#endif
    class Buffer;
    friend class Buffer;
  protected:
    WINDOW_TYPE *win;

    void drawVector(complex<double> base,complex<double> head,complex<double> arrow) {
      OBSOLETE("drawVector meaningless in XYGraphTool");
    }
    void colorQuad(double,double*,double*) {
      OBSOLETE("colorQuad meaningless in XYGraphTool");
    }
  public:
    XYGraphTool(const char*,const char*,const char*,double,double,double,double,
      VirtualColormap*,char* =0,double=HALF,bool =FALSE);
    ~XYGraphTool();

    WINDOW_TYPE* getWindow() { return win; }
    const WINDOW_TYPE* getWindow() const { return win; } 

    void copyFrom(const XYGraphTool::Buffer&);
    void writeXPM(const char *filename) const { win->writeXPM(filename); }
    void expose() const { win->expose(); }
    void eventLoop() { win->eventLoop(); }

    virtual void printOn(ostream& os = cerr) const;
};

class XYGraphTool::Buffer : public VirtualGraphTool12 {
  friend class XYGraphTool;
  private:
    WINDOW_TYPE::Buffer* win;

#if (SPACEDIM>1)
    void drawVector(complex<double>,complex<double>,complex<double>) {;}
    void colorQuad(double,double*,double*) {;}
#endif
  public:
    Buffer(XYGraphTool*);
    ~Buffer();
    WINDOW_TYPE::Buffer* getWindow() { return win; } 
    const WINDOW_TYPE::Buffer* getWindow() const { return win; } 
    virtual void printOn(ostream& os = cerr) const;
};
#endif
