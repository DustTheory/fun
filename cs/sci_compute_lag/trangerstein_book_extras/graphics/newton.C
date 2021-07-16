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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/newton.C,v 1.1 2009/08/20 17:31:48 johnt Exp $"
#include <float.h>
#include <fstream>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
#include <iostream>
#include <math.h>
#include <unistd.h>

//#include "MemoryDebugger.H"
#include "Palette.H"
//#include "SetTraps.H"
//#include "TimedObject.H"
//#include "Tracer.H"
#include "VGT.H"
#include "XYGraphTool.H"

#define npts 100

inline double f(double x) { return x*x-4.*sin(x); }
inline double fp(double x) { return 2.*x-4.*cos(x); }
double A=-4.,B=4.;

//inline double f(double x) { return atan(x); }
//inline double fp(double x) { return 1./(1.+x*x); }
//double A=-4.,B=4.;

//inline double f(double x) { return exp(x)*cos(x); }
//inline double fp(double x) { return exp(x)*(cos(x)-sin(x)); }
//double A=-4.,B=4.;

void computePlotPoints(double a,double b,double &fmin,double &fmax,
double *xarray,double *farray) {
//TRACER_CALL(t,"computePlotPoints");
  double dx=(b-a)*0.0625;
  a-=dx; b+=dx;

  dx=(b-a)/double(npts);
  xarray[0]=a; farray[0]=f(a);
  fmax=fmin=f(a);
  for (int i=1;i<=npts;i++) {
    xarray[i]=xarray[i-1]+dx; farray[i]=f(xarray[i]);
    fmin=min(fmin,farray[i]);
    fmax=max(fmax,farray[i]);
  }
}

void plotPoints(const double *xarray,const double *farray,
XYGraphTool &gt) {
//TRACER_CALL(t,"plotPoints");
  gt.movePen(xarray[0],farray[0]);
  for (int i=1;i<=npts;i++) {
    gt.drawLine(xarray[i],farray[i]);
  }
}

void plotTangentLine(double x,double fx,double fpx,XYGraphTool &gt) {
//TRACER_CALL(t,"plotTangentLine")
//cout << "\tx = " << x << endl;
//cout << "\tfx = " << fx << endl;
//cout << "\tfpx = " << fpx << endl;
  double xlo=gt.getLowX();
  double xhi=gt.getHighX();
  double ylo=gt.getLowY();
  double yhi=gt.getHighY();
  gt.movePen(x,fx);
  double y=fx+fpx*(xlo-x);
  if (y<ylo) {
    gt.drawLine((ylo-fx)/fpx+x,ylo);
  } else if (y>yhi) {
    gt.drawLine((yhi-fx)/fpx+x,yhi);
  } else {
    gt.drawLine(xlo,y);
  }
  gt.movePen(x,fx);
  y=fx+fpx*(xhi-x);
  if (y<ylo) {
    gt.drawLine((ylo-fx)/fpx+x,ylo);
  } else if (y>yhi) {
    gt.drawLine((yhi-fx)/fpx+x,yhi);
  } else {
    gt.drawLine(xhi,y);
  }
}

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

    int nsteps;
    cout << "enter number steps" << endl;
//  the following line generates a memory leak
    cin >> nsteps;
    if (nsteps<0) nsteps=0;
    cout << "number steps = " << nsteps << endl;

    double x;
//  cout << "enter A" << endl;
//  cin >> A;
//  cout << "enter B" << endl;
//  cin >> B;
    if (A>B) { x=A; A=B; B=x; }
    cout << "a,b = " << A << " " << B << endl;

    double dx=(B-A)*0.0078125;
    double xarray[npts+1],farray[npts+1];
    double fmin,fmax;
    computePlotPoints(A,B,fmin,fmax,xarray,farray);
//  cout << "a,b = " << A << " " << B << endl;

//  setup interactive graphics
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("newton","x","f",A,B,fmin,fmax,&cmap,0,winsize);

//  initialize graphics display
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();
//  draw function
    gt.setfgColor("blue");
    plotPoints(xarray,farray,gt);
//  force X to perform the requests
    gt.flush();

    double y;
    cout << "click mouse in window to choose starting guess" <<endl;
    int button=gt.getMouse(x,y);
//  cout << "enter x" << endl;
//  cin >> x;
//  x=max(A,min(B,x));

//  loop over iterations
    for (int i=0;i<nsteps;i++) {
      double fx=f(x), fpx=fp(x);
      cout << "f(" << x << ") = " << fx << ", f'(x) = " << fpx << endl;

      gt.newPage();
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      plotPoints(xarray,farray,gt);

//    mark points
      gt.setfgColor("red");
      plotTangentLine(x,fx,fpx,gt);
      gt.setfgColor("green");
      gt.drawCross(x,fx,dx);
      gt.flush();
      gt.writeXPM("newton");
      XYGraphTool::WINDOW_TYPE::QuitButton qb;

      x=x-fx/fpx;
    }
  } // pal,cmap,gt go out of scope here 
  }
  return EXIT_SUCCESS;
}
