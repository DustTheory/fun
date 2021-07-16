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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/testVolGraphTool3.C,v 1.1 2009/08/20 17:31:48 johnt Exp $"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include "GTKWindow.H"
//#include "MemoryDebugger.H"
#include "Palette.H"
//#include "Tracer.H"
#include "TypedPlotObj.H"
#include "Vector3.H"
#include "VolGraphTool.H"
#include <unistd.h>

class SphereDataSet {
  private:
    VolGraphTool *graph_tool;
    const Vector3 &low;
    const Vector3 &high;
    int npts;
    double radius;
  public:
    SphereDataSet(VolGraphTool *gt,const Vector3 &l,const Vector3 &h,
    int n) : graph_tool(gt),low(l),high(h),npts(n) {
      radius=sumSquares((high-low)*0.5);
    }
    virtual ~SphereDataSet() { }
    VolGraphTool* getGraphTool() const { return graph_tool; }
    bool drawSurface() const { return false; }
    void plot(bool,const AxisClipPlane*);
    void plot(bool) {}
    virtual void printOn(ostream&) const {}
};


void SphereDataSet::plot(bool,const AxisClipPlane *cp) {
//GLWindow::DrawDelimiter dd(graph_tool->getWindow(),cp);
  double location=graph_tool->getLocation(cp);
  unsigned int dir=cp->getDirection();
  double fudge=graph_tool->getLength()[dir]*0.0009765625;
  if (cp->getHand()==LEFT_HAND) location += fudge;
  else location -= fudge;
  int d1=(dir+1)%3,d2=(dir+2)%3;
  
  Vector3 dif=high-low;
  Vector3 l,h;
  l[dir]=h[dir]=location;

  double dv1=dif[d1]/double(2*npts);
  double dv2=dif[d2]/double(2*npts);
  for (int i2=-npts;i2<=npts;i2++) {
    double v2l=low[d2]+dv2*double(i2+npts);
    double v2=v2l+0.5*dv2;
    l[d2]=v2l;
    h[d2]=v2l+dv2;
    for (int i1=-npts;i1<=npts;i1++) {
      double v1l=low[d1]+dv1*double(i1+npts);
      double v1=v1l+0.5*dv1;
      l[d1]=v1l; 
      h[d1]=v1l+dv1;
      double f=location*location+v1*v1+v2*v2;
      graph_tool->colorRect(f/radius,l,h,dir);
    }
  }
}

int main(int argc,char* argv[]) {
  cout << boolalpha;
  {
#ifdef MEM_DEBUG
//  MemoryDebugger md(1);
#endif
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

//  window size, as fraction of min(screen width, screen height)
    double winsize=0.5;

//  user coordinates
    Vector3 low(-1.,-1.,-1.);
    Vector3 high(1.,1.,1.);

//open graphics output window
    Palette *pal=OPERATOR_NEW Palette();
    bool use_lighting=false;
    VolGraphTool *g= OPERATOR_NEW VolGraphTool("VolGraphTool test",
      low,high,*pal,use_lighting,winsize,false);

//OpenGL needs to know which display to work on
    { VolGraphTool::DrawDelimiter dd(true,g);
      g->drawBox(low,high);
    }
    g->newPage();
    g->expose();

    SphereDataSet *ds=OPERATOR_NEW SphereDataSet(g,low,high,10);
    TypedPlotObj<SphereDataSet> *tpo=
      OPERATOR_NEW TypedPlotObj<SphereDataSet>(ds);
//  graph_tool->setIsoSurfaceFrac(ZERO);
    g->initialize(tpo);
    g->plotClipPlanes();
    g->eventLoop();
    g->terminatePlotObj();

    if (ds!=0) delete ds; ds=0;
    if (g!=0) delete g; g=0;
    if (pal!=0) delete pal; pal=0;
  }
  return EXIT_SUCCESS;

//procedures not tested:
}

template class TypedPlotObj<SphereDataSet>;
