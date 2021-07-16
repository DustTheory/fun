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
#include <float.h>
#include <fstream>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
#include <iostream>
#include <limits>
#include <math.h>
#include <unistd.h>

#include "AreaGraphTool.H"
//#include "MemoryDebugger.H"
#include "Palette.H"
//#include "Tracer.H"

int ncells=100;
int ncorn=101;

inline double f1(double x) { return sqrt(max(0.,4.-2.*x*x)); }
inline double f2(double x) { return sqrt(max(0.,(9.-x*x)/3.)); }

double rho=10000.;
inline double a1(double x,double y) { return 2.*x*x+y*y-4.; }
inline double a2(double x,double y) { return x*x+3.*y*y-9.; }
inline double phi(double x,double y) { return -x*x-y*y; }
inline double penalty(double x,double y) {
  return phi(x,y)+rho*(pow(max(0.,a1(x,y)),2)+pow(max(0.,a2(x,y)),2)
    + pow(max(0.,-x),2)+pow(max(0.,-y),2));
}

int main(int argc,char* argv[]) {
  cout << boolalpha;
  { 
//  MemoryDebugger md(1);
    {
    double rlo[2]={-0.1,-0.1};
    double rhi[2]={3.1,2.1};

    double zmin=numeric_limits<double>::max();
    double zmax=-zmin;
    double z[ncorn*ncorn];
    double dx=(rhi[0]-rlo[0])/static_cast<double>(ncells);
    double dy=(rhi[1]-rlo[1])/static_cast<double>(ncells);
    for (int j=0;j<ncorn;j++) {
      double y=rlo[1]+dy*double(j);
      for (int i=0;i<ncorn;i++) {
        double x=rlo[0]+dx*double(i);
        double zz=penalty(x,y);
        z[i+j*ncorn]=zz;
        if (zmin>zz) zmin=zz;
        if (zmax<zz) zmax=zz;
      }
    }
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    Palette pal;
    AreaGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    AreaGraphTool gt("Penalty Function Contours","x","y",rlo,rhi,&cmap,0,
      winsize);
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

    double ddx=sqrt(2.)/static_cast<double>(ncells);
    double x=0.;
    gt.movePen(x,f1(x));
    for (int i=0;i<=ncells;i++,x+=ddx) {
      double y=f1(x);
      if (y<0.) break;
      gt.drawLine(x,y);
    }
    ddx=3./static_cast<double>(ncells);
    x=0.;
    gt.movePen(x,f2(x));
    for (int i=0;i<=ncells;i++,x+=ddx) {
      double y=f2(x);
      if (y<0.) break;
      gt.drawLine(x,y);
    }

    double dif[5],xpos[5],ypos[5],xy[4][2];
    int ncontours=30;
    double ratio=exp(log(zmax-zmin)/static_cast<double>(ncontours));
    for (int ic=0;ic<ncontours;ic++) {
      double value=zmin+pow(ratio,ic);
      gt.setfgColor(ic,ncontours);
      for (int j=0;j<ncells;j++) {
        ypos[4]=ypos[1]=ypos[0]=rlo[1]+dy*double(j);
        ypos[3]=ypos[2]=ypos[1]+dy;
        for (int i=0;i<ncells;i++) {
          xpos[4]=xpos[3]=xpos[0]=rlo[0]+dx*double(i);
          xpos[2]=xpos[1]=xpos[0]+dx;
          int ij=i+j*ncorn;
          dif[0]=z[ij]-value;
          dif[1]=z[ij+1]-value;
          dif[2]=z[ij+ncorn+1]-value;
          dif[3]=z[ij+ncorn]-value;
          dif[4]=dif[0];
          int count=0,first,last;
          for (int k=0;k<4;k++) {
            if (dif[k]*dif[k+1]<=0. && abs(dif[k+1])>0.) {
              if (count==0) first=k;
              else last=k;
              count++;
              double pos=dif[k]/(dif[k]-dif[k+1]);
              xy[k][0]=xpos[k]+pos*(xpos[k+1]-xpos[k]);
              xy[k][1]=ypos[k]+pos*(ypos[k+1]-ypos[k]);
            }
          }
          if (count==4) {
            gt.movePen(xy[0][0],xy[0][1]);
            gt.drawLine(xy[2][0],xy[2][1]);
            gt.movePen(xy[1][0],xy[1][1]);
            gt.drawLine(xy[3][0],xy[3][1]);
          } else if (count==2) {
            gt.movePen(xy[first][0],xy[first][1]);
            gt.drawLine(xy[last][0],xy[last][1]);
          }
        }
      }   
    }   
    gt.flush();
    AreaGraphTool::WINDOW_TYPE::QuitButton qb;
    } // pal,cmap,gt go out of scope here 
  }
  return EXIT_SUCCESS;
}
