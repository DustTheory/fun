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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/testVolGraphTool5.C,v 1.1 2009/08/20 17:31:49 johnt Exp $"
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>
#include "GTKWindow.H"
//#include "MemoryDebugger.H"
#include "SurfaceDataSet.H"
#include "Surface.H"
#include "ToggleCallback.H"
#include "VolGraphTool.H"

double surface_fun(const Vector3 &x) {
//return x[0]*x[0]+x[1]*x[1]-x[2]*x[2];
//return x[2]*x[2]-x[0]*x[0]-x[1]*x[1];
  return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); // sphere
}

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
