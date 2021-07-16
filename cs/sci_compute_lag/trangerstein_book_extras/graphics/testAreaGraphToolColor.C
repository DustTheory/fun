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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/testAreaGraphToolColor.C,v 1.1 2009/08/20 17:31:48 johnt Exp $"
#include <limits.h>
#include <iostream>
#include <stdlib.h>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
#include "AreaGraphTool.H"
//#include "MemoryDebugger.H"
# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif
#include "Palette.H"
#if (SPACEDIM!=2)
#include "barf"
#endif

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
    double rlo[2]; rlo[0]=0.; rlo[1]=0.;
    double rhi[2]; rhi[0]=1.; rhi[1]=1.;
    double zmin=HUGE_VAL,zmax=-HUGE_VAL;

//  define the surface via a uniform lattice
    int npts=50;
    double *z=OPERATOR_NEW_BRACKET(double,npts*npts);
    double dx=(rhi[0]-rlo[0])/double(npts);
    double dy=(rhi[1]-rlo[1])/double(npts);
    for (int j=0;j<npts;j++) {
      double y=rlo[1]+dy*double(j);
      for (int i=0;i<npts;i++) {
        double x=rlo[0]+dx*double(i);
        double zz=x*x+y*y;                      // definition of surface
        z[i+j*npts]=zz;
        if (zmin>zz) zmin=zz;
        if (zmax<zz) zmax=zz;
      }
    }
    CHECK_TEST(zmax>zmin);

//open graphics output window
    Palette pal;
    AreaGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    AreaGraphTool g("AreaGraphTool Color Test","x","y",rlo,rhi,&cmap,0,
      winsize);
    g.drawAxes();

//color the image
    double xpos[4],ypos[4];
    for (int j=0;j<npts;j++) {
      ypos[1]=ypos[0]=rlo[1]+dy*double(j);
      ypos[3]=ypos[2]=ypos[1]+dy;
      for (int i=0;i<npts;i++) {
        xpos[3]=xpos[0]=rlo[0]+dx*double(i);
        xpos[2]=xpos[1]=xpos[0]+dx;
        double ratio=(z[i+j*npts]-zmin)/(zmax-zmin);
        ratio=max(0.,min(1.,ratio));
        g.colorQuad(ratio,xpos,ypos);
      }
    }
    g.flush();
    AreaGraphTool::WINDOW_TYPE::QuitButton qb;

    if (z) delete [] z;
  }
  return EXIT_SUCCESS;
}
