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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/testVolGraphTool2.C,v 1.1 2009/08/20 17:31:48 johnt Exp $"
#include <limits.h>
#include <iostream>
#include <stdlib.h>
#include "GTKWindow.H"
//#include "MemoryDebugger.H"
#include "Palette.H"
//#include "Tracer.H"
#include "Vector3.H"
#include "VolGraphTool.H"

void plot(VolGraphTool *g,int npts,const Vector3 &low,
const Vector3 &high,double dx,double dy,double *z) {
  int n2p1=2*npts+1;
  double valdif=ONE/(high[2]-low[2]);
  for (int j=-npts;j<npts;j++) {
    int jn=j+npts;
    double y=low[1]+dy*double(jn);
    for (int i=-npts;i<npts;i++) {
      int in=i+npts;
      double x=low[0]+dx*double(in);
      int k=in+jn*n2p1;
      Vector3 v0(x,y,z[k]);
      Vector3 v1(x+dx,y,z[k+1]);
      Vector3 v2(x+dx,y+dy,z[k+1+n2p1]);
      Vector3 v3(x,y+dy,z[k+n2p1]);
      Vector3 vave=(v0+v1+v2+v3)*FOURTH;

      double frac=((v0[2]+v1[2]+vave[2])*THIRD-low[2])*valdif;
      frac=max(ZERO,min(ONE,frac));
      g->colorTriangle(frac,v0,v1,vave);

      frac=((v1[2]+v2[2]+vave[2])*THIRD-low[2])*valdif;
      frac=max(ZERO,min(ONE,frac));
      g->colorTriangle(frac,v1,v2,vave);

      frac=((v2[2]+v3[2]+vave[2])*THIRD-low[2])*valdif;
      frac=max(ZERO,min(ONE,frac));
      g->colorTriangle(frac,v2,v3,vave);

      frac=((v3[2]+v0[2]+vave[2])*THIRD-low[2])*valdif;
      frac=max(ZERO,min(ONE,frac));
      g->colorTriangle(frac,v3,v0,vave);
    }
  }
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

//  window size, as fraction of min(screen width, screen height)
    double winsize=0.5;

//  user coordinates
    Vector3 low(-1.,-1.,DBL_MAX);
    Vector3 high(1.,1.,-DBL_MAX);

//  define the surface via a uniform lattice
    int npts=5;
    int n2p1=2*npts+1;
    double *z=OPERATOR_NEW_BRACKET(double,n2p1*n2p1);
    double dx=(high[0]-low[0])/double(2*npts);
    double dy=(high[1]-low[1])/double(2*npts);
    for (int j=-npts;j<=npts;j++) {
      int jn=j+npts;
      double y=low[1]+dy*double(jn);
      for (int i=-npts;i<=npts;i++) {
        int in=i+npts;
        double x=low[0]+dx*double(in);
        double zz=x*x+y*y;          // definition of surface
        z[in+jn*n2p1]=zz;
        if (low[2]>zz) low[2]=zz;
        if (high[2]<zz) high[2]=zz;
      }
    }

//open graphics output window
    Palette *pal=OPERATOR_NEW Palette();
    bool use_lighting=false;
    VolGraphTool *g= OPERATOR_NEW VolGraphTool("testVolGraphTool2",
      low,high,*pal,use_lighting,winsize);

//OpenGL needs to know which DisplayList to work on
    { VolGraphTool::DrawDelimiter dd(true,g);
      g->drawBox(low,high);
    }
    g->newPage();
    g->expose();

    { bool do_best=false;
      VolGraphTool::DrawDelimiter dd(do_best,g);
      plot(g,npts,low,high,dx,dy,z);
    }
    { bool do_best=true;
      VolGraphTool::DrawDelimiter dd(do_best,g);
      plot(g,npts,low,high,dx,dy,z);
    }
    g->eventLoop(false);

    if (g) delete g; g=0;
    if (pal) delete pal; pal=0;
    if (z) delete [] z;
  }
  return EXIT_SUCCESS;
}
