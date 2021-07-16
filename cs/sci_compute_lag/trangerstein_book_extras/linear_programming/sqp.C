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

int ncells=512;
int ncorn=513;

inline double f1(double x) { return sqrt(max(0.,4.-2.*x*x)); }
inline double f2(double x) { return sqrt(max(0.,(9.-x*x)/3.)); }

double safety=1.e-3;
double tau=10.;
inline double a1(const double *x) { return 2.*x[0]*x[0]+x[1]*x[1]-4.; }
inline void da1dx(const double *x,double *g) { g[0]=4.*x[0]; g[1]=2.*x[1]; }
inline double a2(const double *x) { return x[0]*x[0]+3.*x[1]*x[1]-9.; }
inline void da2dx(const double *x,double *g) { g[0]=2.*x[0]; g[1]=6.*x[1]; }
inline double phi(const double *x) { return -x[0]*x[0]-x[1]*x[1]; }
inline void dphidx(const double *x,double *g) {
  g[0]=-2.*x[0]; g[1]=-2.*x[1];
}
inline void d2phidx2(const double *x,double *H) {
  H[0]=-2.; H[1]=0.; H[2]=0.; H[3]=-2.;
}
inline double quadratic(const double *x,const double *g,const double *H) {
  double Hx0=H[0]*x[0]+H[2]*x[1];
  double Hx1=H[1]*x[0]+H[3]*x[1];
  return g[0]*x[0]+g[1]*x[1]+0.5*(x[0]*Hx0+x[1]*Hx1);
}
inline double linear(const double *x,double v,const double *g) {
  return v+g[0]*x[0]+g[1]*x[1];
}

int main(int argc,char* argv[]) {
  cout << boolalpha;
  { 
//  MemoryDebugger md(1);
    {
//  double xref[2]={sqrt(0.6),sqrt(2.8)};
    double xref[2]={2.0,0.5};
//  double xref[2]={0.5,1.7};
    double gref[2]; dphidx(xref,gref);
    double Href[4]; d2phidx2(xref,Href);
    double a1ref=a1(xref);
    double a2ref=a2(xref);
    double da1refdx[2]; da1dx(xref,da1refdx);
    double da2refdx[2]; da2dx(xref,da2refdx);

    double rlo[2]={0.,0.};
//  double rhi[2]={1.42,1.74};
    double rhi[2]={xref[0]-max((a1ref-da1refdx[1]*xref[1])/da1refdx[0],
                               (a2ref-da2refdx[1]*xref[1])/da2refdx[0]),
                   xref[1]-max((a1ref-da1refdx[0]*xref[0])/da1refdx[1],
                               (a2ref-da2refdx[0]*xref[0])/da2refdx[1])};

    double zmin=numeric_limits<double>::max();
    double zmax=-zmin;
    double z[ncorn*ncorn];
    double dx=(rhi[0]-rlo[0])/static_cast<double>(ncells);
    double dy=(rhi[1]-rlo[1])/static_cast<double>(ncells);
    double x[2];
    for (int j=0;j<ncorn;j++) {
      x[1]=rlo[1]+dy*double(j);
      for (int i=0;i<ncorn;i++) {
        x[0]=rlo[0]+dx*double(i);
        double dif[2]={x[0]-xref[0],x[1]-xref[1]};
        double zz=quadratic(dif,gref,Href);
        z[i+j*ncorn]=zz;
        double c1=linear(dif,a1ref,da1refdx);
        double c2=linear(dif,a2ref,da2refdx);
        if (c1<=0. && c2<=0.) {
          zmin=min(zmin,zz);
          zmax=max(zmax,zz);
        }
//      cout << "\ti,j,x = " << i << " " << j << " " << x[0] << " "
//        << x[1] << endl;
//      cout << "\tz,c1,c2 = " << zz << " " << c1 << " " << c2 << endl;
      }
    }
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    Palette pal;
    AreaGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    AreaGraphTool gt("Quadratic Function Contours","x","y",rlo,rhi,&cmap,0,
      winsize);
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

    double yy=xref[1]-(a1ref+da1refdx[0]*(rlo[0]-xref[0]))/da1refdx[1];

    if (yy<rlo[1]) {
      gt.movePen(xref[0]-(a1ref+da1refdx[1]*(rlo[1]-xref[1]))/da1refdx[0],
        rlo[1]);
    } else if (yy>rhi[1]) {
      gt.movePen(xref[0]-(a1ref+da1refdx[1]*(rhi[1]-xref[1]))/da1refdx[0],
        rhi[1]);
    } else gt.movePen(rlo[0],yy);
    yy=xref[1]-(a1ref+da1refdx[0]*(rhi[0]-xref[0]))/da1refdx[1];
    if (yy<rlo[1]) {
      gt.drawLine(xref[0]-(a1ref+da1refdx[1]*(rlo[1]-xref[1]))/da1refdx[0],
        rlo[1]);
    } else if (yy > rhi[1]) {
      gt.drawLine(xref[0]-(a1ref+da1refdx[1]*(rhi[1]-xref[1]))/da1refdx[0],
        rhi[1]);
    } else gt.drawLine(rhi[0],yy);

    yy=xref[1]-(a2ref+da2refdx[0]*(rlo[0]-xref[0]))/da2refdx[1];
    if (yy<rlo[1]) {
      gt.movePen(xref[0]-(a2ref+da2refdx[1]*(rlo[1]-xref[1]))/da2refdx[0],
        rlo[1]);
    } else if (yy > rhi[1]) {
      gt.movePen(xref[0]-(a2ref+da2refdx[1]*(rhi[1]-xref[1]))/da2refdx[0],
        rhi[1]);
    } else gt.movePen(rlo[0],yy);
    yy=xref[1]-(a2ref+da2refdx[0]*(rhi[0]-xref[0]))/da2refdx[1];
    if (yy<rlo[1]) {
      gt.drawLine(xref[0]-(a2ref+da2refdx[1]*(rlo[1]-xref[1]))/da2refdx[0],
        rlo[1]);
    } else if (yy > rhi[1]) {
      gt.drawLine(xref[0]-(a2ref+da2refdx[1]*(rhi[1]-xref[1]))/da2refdx[0],
        rhi[1]);
    } else gt.drawLine(rhi[0],yy);

    double dif[5],xpos[5],ypos[5],xy[4][2];
    int ncontours=30;
    double ratio=(zmax-zmin)/static_cast<double>(ncontours+1);
    double value=zmin+ratio;
    for (int ic=1;ic<=ncontours;ic++,value+=ratio) {
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
//        cout << "\ti,j,count = " << i << " " << j << " " << count
//          << endl;
          double dd[2]={xpos[2]-xref[0],ypos[2]-xref[1]};
          double c1=linear(dd,a1ref,da1refdx);
          double c2=linear(dd,a2ref,da2refdx);
          if (c1<=0. && c2<=0.) {
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
    }   
    gt.flush();
    AreaGraphTool::WINDOW_TYPE::QuitButton qb;
    } // pal,cmap,gt go out of scope here 
  }
  return EXIT_SUCCESS;
}
