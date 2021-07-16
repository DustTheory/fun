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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/testAreaGraphToolContour.C,v 1.1 2009/08/20 17:31:48 johnt Exp $"
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
    int ncells=50;
    int ncorn=ncells+1;
    double *z=OPERATOR_NEW_BRACKET(double,ncorn*ncorn);
    double dx=(rhi[0]-rlo[0])/double(ncells);
    double dy=(rhi[1]-rlo[1])/double(ncells);
    for (int j=0;j<ncorn;j++) {
      double y=rlo[1]+dy*double(j);
      for (int i=0;i<ncorn;i++) {
        double x=rlo[0]+dx*double(i);
        double zz=x*x+y*y;                      // definition of surface
        z[i+j*ncorn]=zz;
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

    g.setfgColor("black");
    g.movePen(rlo[0],rlo[1]);
    g.drawLine(rhi[0],rlo[1]);
    g.drawLine(rhi[0],rhi[1]);
    g.drawLine(rlo[0],rhi[1]);
    g.drawLine(rlo[0],rlo[1]);

//contour the image
    double dif[5],xpos[5],ypos[5],xy[4][2];
    int ncontours=30;
    for (int ic=0;ic<ncontours;ic++) {
      double value=zmin+double(ic+1)*(zmax-zmin)/double(ncontours+1);
      g.setfgColor(ic,ncontours);
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
            g.movePen(xy[0][0],xy[0][1]);
            g.drawLine(xy[2][0],xy[2][1]);
            g.movePen(xy[1][0],xy[1][1]);
            g.drawLine(xy[3][0],xy[3][1]);
          } else if (count==2) {
            g.movePen(xy[first][0],xy[first][1]);
            g.drawLine(xy[last][0],xy[last][1]);
          }
        }
      }
    }
    g.flush();
    AreaGraphTool::WINDOW_TYPE::QuitButton qb;

    delete [] z;
  }
  return EXIT_SUCCESS;
}
