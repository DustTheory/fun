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
#include <limits.h>
#include <iostream>
#include <stdlib.h>
//#define SPACEDIM 3
#include "GTKWindow.H"
//#include "MemoryDebugger.H"
#include "Palette.H"
//#include "Tracer.H"
#include "Vector3.H"
#include "VolGraphTool.H"

double c[3]={2.,2.,1.}; // assume components are nonincreasing

void plot(VolGraphTool *g) {
  double cmin=c[2];
  double cmax=c[0];
  double cdif=(cmax-cmin)/6.;
  double val=cmin+cdif;
  for (int i=0;i<5;i++,val+=cdif) {
    switch (i) {
      case 4 :
        g->setfgColor("blue");
        break;
      case 3 :
        g->setfgColor("cyan");
        break;
      case 2 :
        g->setfgColor("green");
        break;
      case 1 :
        g->setfgColor("yellow");
        break;
      case 0 :
        g->setfgColor("red");
        break;
    }
    if (val<c[1]) {
      double x1=(val-c[2])/(c[1]-c[2]);
      Vector3 start(0.,x1,1.-x1);
      double x0=(val-c[2])/(c[0]-c[2]);
      Vector3 finish(x0,0.,1.-x0);
      g->drawLine(start,finish);
    } else {
      double x0=(val-c[1])/(c[0]-c[1]);
      Vector3 start(x0,1.-x0,0.);
      x0=(val-c[2])/(c[0]-c[2]);
      Vector3 finish(x0,0.,1.-x0);
      g->drawLine(start,finish);
    }
  }
  g->setfgColor("black");
  Vector3 origin(0.,0.,0.);
  Vector3 x0(1.,0.,0.);
  Vector3 x1(0.,1.,0.);
  Vector3 x2(0.,0.,1.);
  g->drawLine(x0,x1);
  g->drawLine(x1,x2);
  g->drawLine(x2,x0);
  g->drawLine(origin,x0);
  g->drawLine(origin,x1);
  g->drawLine(origin,x2);
  x0[0]=0.9;
  x1[1]=0.9;
  x2[2]=0.9;
  g->putString(&x0,"x");
  g->putString(&x1,"y");
  g->putString(&x2,"z");
}
#define OPERATOR_NEW new

int main(int argc,char* argv[]) {
  cout << boolalpha;
#ifdef MEM_DEBUG
//MemoryDebugger md(1);
#endif
  { 
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    double winsize=0.5;
    Vector3 low(-0.05,-0.05,-0.05);
    Vector3 high(1.05,1.05,1.05);

    Palette *pal=OPERATOR_NEW Palette();
    bool use_lighting=false;
//  bool rotate_best=false;
//  bool show_bounding_box=true;
    VolGraphTool *g=
      OPERATOR_NEW VolGraphTool("LP standard form objective",
      low,high,*pal,use_lighting,winsize
      /*,rotate_best,show_bounding_box */);
    { VolGraphTool::DrawDelimiter dd(true,g);
      g->drawBox(low,high);
    }
    g->newPage();
    g->expose();

    { bool do_best=false;
      VolGraphTool::DrawDelimiter dd(do_best,g);
      plot(g);
    }
    { bool do_best=true;
      VolGraphTool::DrawDelimiter dd(do_best,g);
      plot(g);
    }
    g->eventLoop(false);

    if (g) delete g; g=0;
    if (pal) delete pal; pal=0;
  }
  return EXIT_SUCCESS;
}
