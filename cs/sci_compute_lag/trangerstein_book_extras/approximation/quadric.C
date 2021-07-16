//**********************************************************************
// Copyright 2014 John A. Trangenstein
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
// "$Header:$"
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>
#include "GTKWindow.H"
//#include "MemoryDebugger.H"
#include "SurfaceDataSet.H"
#include "Surface.H"
#include "ToggleCallback.H"
#include "VolGraphTool.H"

class QuadricPlotObj : public VirtualGLWindow::PlotObj {
  private:
    VirtualVolGraphTool *graph_tool;
  public:
//  t must be constructed via operator new:
    QuadricPlotObj(VirtualVolGraphTool *gt) : graph_tool(gt) { 
      CHECK_POINTER(gt)
      GLuint startList=glGenLists(4);
      GLUquadricObj *qobj=gluNewQuadric();
      gluQuadricDrawStyle(qobj, GLU_FILL); /* smooth shaded */
      gluQuadricNormals(qobj, GLU_SMOOTH);
      glNewList(startList, GL_COMPILE);
        gluSphere(qobj, 0.75, 15, 10);
      glEndList();


    }
    ~QuadricPlotObj() { }
    virtual bool drawSurface() const { return true;}
    virtual void plot(bool do_best) const { 
      dataset->plot(do_best);
    }
# if (SPACEDIM==3)
    virtual void plot(bool do_best,const AxisClipPlane *cp) const {
    }
# endif
};

int main(int argc,char* argv[]) {
  cout << boolalpha;
#ifdef MEM_DEBUG
//MemoryDebugger md(1);
#endif
  {
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    int nptsf=15,nptsc=5;
    double vmin=1.0, vmax=5.0, iso_value=3.0;

//  open graphics output window
    Vector3 low(-5.0,-5.0,-5.0);
    Vector3 high(5.0,5.0,5.0);
    double winsize=0.5; // fraction of min(screen width, height)

    Palette *pal=OPERATOR_NEW Palette();
    bool use_lighting=true;
    VolGraphTool *g= OPERATOR_NEW VolGraphTool("VolGraphTool test", 
      low, high,*pal,use_lighting,winsize);
    g->newPage();

    {
      SurfaceDataSet dataset(g,low,high,nptsf,nptsc,surface_fun, 
        vmin,vmax);
      Surface<SurfaceDataSet> *object=
        OPERATOR_NEW Surface<SurfaceDataSet>(&dataset, g, 16);

      g->initialize(object);
      g->eventLoop();
      g->terminatePlotObj();
    }

    if (g) delete g; g=0;
    if (pal) delete pal; pal=0;
  }
}

template class Surface<SurfaceDataSet>;
template class ToggleButtonArrayCallbackStructure< Surface<SurfaceDataSet> >;
