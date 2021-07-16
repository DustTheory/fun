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

inline double f1(double x) { return sqrt(4.-2.*x*x); }
inline double f2(double x) { return sqrt((9.-x*x)/3.); }
//double A=-4.,B=4.;

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
    double rhi[2]={3.,2.};
    AreaGraphTool gt("Kuhn Tucker","x","y",rlo,rhi,&cmap,0,winsize);
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

    double x=sqrt(0.6);
    complex<double> base(x,f1(x));
    complex<double> head(x,f1(x));
    complex<double> head1(4.*x,2.*f1(x));
    complex<double> head2(0.05*2.*x,0.05*6.*f2(x));
    complex<double> arrow(-0.3*cos(0.3),0.3*sin(0.3));
    head*=0.25/abs(head);
    head1*=0.25/abs(head1);
    head2*=0.25/abs(head2);
    gt.setfgColor(79,80);
    double dx=sqrt(2.)/static_cast<double>(npts);
    x=0.;
    gt.movePen(x,f1(x));
    for (int i=0;i<=npts;i++,x+=dx) {
      gt.drawLine(x,f1(x));
    }
    gt.drawVector(base,head1,arrow);
    gt.setfgColor(39,80);
    dx=3./static_cast<double>(npts);
    x=0.;
    gt.movePen(x,f2(x));
    for (int i=0;i<=npts;i++,x+=dx) {
      gt.drawLine(x,f2(x));
    }
    gt.drawVector(base,head2,arrow);

    double dval=sqrt(3.2)/static_cast<double>(nc);
    double val=dval;
    for (int i=1;i<=nc;i++,val+=dval) {
      dx=val/static_cast<double>(npts);
      x=0.;
      gt.setfgColor(nc-i+1,nc);
      double y=sqrt(val*val-x*x);
      bool inbounds=(y<=f1(x) && y<=f2(x));
      if (inbounds) gt.movePen(x,sqrt(val*val-x*x));
      for (int j=0;j<=npts;j++,x+=dx) {
        y=sqrt(max(0.,val*val-x*x));
        if (inbounds) {
          inbounds=(y<=f1(x) && y<=f2(x));
          if (inbounds) gt.drawLine(x,y);
        } else {
          inbounds=(y<=f1(x) && y<=f2(x));
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
