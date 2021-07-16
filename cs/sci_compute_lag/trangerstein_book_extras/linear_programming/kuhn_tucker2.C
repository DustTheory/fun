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
#include <float.h>
#include <fstream>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
#include <iostream>
#include <math.h>
#include <unistd.h>

#include "AreaGraphTool.H"
//#include "MemoryDebugger.H"
#include "Palette.H"
//#include "SetTraps.H"
//#include "TimedObject.H"
//#include "Tracer.H"
#include "VGT.H"

#define npts 100
#define nc 20

inline double f1(double x) { return pow(1.-x*x*x,1./3.); }

int main(int argc,char* argv[]) {
  cout << boolalpha;
#ifdef DEBUG
//setTraps();
#endif

  { 
//  MemoryDebugger md(1);
    {
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    Palette pal;
    AreaGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    double rlo[2]={0.,0.};
    double rhi[2]={1.,1.};
    AreaGraphTool gt("Kuhn Tucker","x","y",rlo,rhi,&cmap,0,winsize);
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

    double x=2./pow(16.,1./3.);
    double dval=x*sqrt(2.)/static_cast<double>(nc);

    complex<double> base(x,f1(x));
    complex<double> head(x,f1(x));
    complex<double> head1(3.*x*x,3.*f1(x)*f1(x));
    complex<double> arrow(-0.3*cos(0.3),0.3*sin(0.3));
    head*=0.25/abs(head);
    head1*=0.25/abs(head1);
//  see Palette::defaultColors
    gt.setfgColor(79,80);
    double dx=1./static_cast<double>(npts);
    x=0.;
    gt.movePen(x,f1(x));
    for (int i=0;i<=npts;i++,x+=dx) {
      gt.drawLine(x,f1(x));
    }
    gt.drawVector(base,head1,arrow);

    double val=dval;
    for (int i=1;i<=nc;i++,val+=dval) {
      dx=val/static_cast<double>(npts);
      x=0.;
      gt.setfgColor(nc-i+1,nc);
      double y=sqrt(val*val-x*x);
      bool inbounds=(y<=f1(x));
      if (inbounds) gt.movePen(x,sqrt(val*val-x*x));
      for (int j=0;j<=npts;j++,x+=dx) {
        y=sqrt(max(0.,val*val-x*x));
        if (inbounds) {
          inbounds=(y<=f1(x));
          if (inbounds) gt.drawLine(x,y);
        } else {
          inbounds=(y<=f1(x));
          if (inbounds) gt.movePen(x,y);
        }
      }
    }
    gt.setfgColor(2,80);
    gt.drawVector(base,head,arrow);
    gt.flush();
    AreaGraphTool::WINDOW_TYPE::QuitButton qb;
    } // pal,cmap,gt go out of scope here 
  }
  return EXIT_SUCCESS;
}
