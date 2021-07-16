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
#ifndef AREA_GRAPH_TOOL_H
#define AREA_GRAPH_TOOL_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/AreaGraphTool.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#if (SPACEDIM>1)
#include <iostream>
#include <complex>
#include <stdio.h>
#include "VGT.H"
#ifdef USE_GTK
#include "GTKWindow.H"
#else
#include "XTWindow.H"
#endif
class VirtualColormap;

class AreaGraphTool : public VirtualGraphTool12 {
  public:
#ifdef USE_GTK
    typedef GTKWindow WINDOW_TYPE;
#else
    typedef XTWindow WINDOW_TYPE;
#endif
    class Buffer;
    friend class Buffer;
  private:
    WINDOW_TYPE *win;  
  protected:
  public:
    AreaGraphTool(const char*,const char*,const char*,const double*,
      const double*,VirtualColormap* =0,char* =0,double=HALF,bool=FALSE);
    AreaGraphTool(const char*,const char*,const char*,
      VirtualColormap* =0,char* =0,double=HALF,double=HALF,
      bool mono=FALSE);
    ~AreaGraphTool() { if (win) delete win; win=0; }
    WINDOW_TYPE* getWindow() { return win; }
    const WINDOW_TYPE* getWindow() const { return win; }
    void writeXPM(const char *filename) const { win->writeXPM(filename); }
    virtual void printOn(ostream& = cerr) const;
};

class AreaGraphTool::Buffer : public VirtualGraphTool12 {
  friend class AreaGraphTool;
  private:
    const AreaGraphTool *parent;
    WINDOW_TYPE::Buffer* win;

#if (SPACEDIM>1)
    void drawVector(complex<double>,complex<double>,complex<double>) { OBSOLETE(0);}
#endif
  public: 
    Buffer(AreaGraphTool*);
    ~Buffer() { if (win) delete win; win=0; }
    WINDOW_TYPE::Buffer* getWindow() { return win; } 
    const WINDOW_TYPE::Buffer* getWindow() const { return win; } 
    void colorQuad(double,double*,double*);
    void printOn(ostream& = cerr) const;
};
#endif
#endif
