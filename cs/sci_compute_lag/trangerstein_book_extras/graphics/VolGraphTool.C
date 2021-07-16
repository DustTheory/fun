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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/VolGraphTool.C,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#include "VolGraphTool.H"
#if (SPACEDIM>1)
#include <iostream>
#include <iomanip>
#include <math.h>
#define LABEL_LENGTH 80
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "Debug.H"
#include "LookupTable.H"
//#include "Tracer.H"
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int cubeIndex(const double value[2][2][2]) {
  int m=1, cube_index=0;
  for (int i0=0;i0<2;i0++) {
    for (int i1=0;i1<2;i1++) {
      for (int i2=0;i2<2;i2++) {
        if (value[i0][i1][i2]>ZERO) cube_index+=m;
        m*=2;
      }
    }
  }
  return cube_index;
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Check if the value at some side center has the same sign as 
// the corresponding corner.
bool sideCenterSameSign(const int &side,const double value[2][2][2]) {
  int hand=side%2;
  int dir=side/2;
  int id1=(dir+1)%3;
  int id2=(dir+2)%3;

  int h0[3]; h0[dir]=hand; h0[id1]=0; h0[id2]=0;
  int h1[3]; h1[dir]=hand; h1[id1]=1; h1[id2]=0;
  int h2[3]; h2[dir]=hand; h2[id1]=0; h2[id2]=1;

  return value[h0[0]][h0[1]][h0[2]]
       *(value[h1[0]][h1[1]][h1[2]]+value[h2[0]][h2[1]][h2[2]]) > ZERO;
}
#endif
#if (SPACEDIM==3)
// :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int subTableIndex(int index,const double value[2][2][2]) {
  int indx=0;
  int indx0=marchingCubesTable[index][2];
  int m=-marchingCubesTable[index][1];
  const double *val=&value[0][0][0];
  switch(m) {
    case 1: {
      int k0=marchingCubesTable[index][5];
      int k1=marchingCubesTable[index][3];
      int k2=marchingCubesTable[index][4];
      int k3=(k1+k2)-k0;
      if((val[k0]*val[k3])>(val[k1]*val[k2])) {
        indx=2*indx0;
      } else {
        indx=2*indx0+1;
      }
      break;
    }
    case 2: {
      int s1=0, s2=0;
      int k1=marchingCubesTable[index][3];
      int k2=marchingCubesTable[index][4];
      int k0=marchingCubesTable[index][5];
      int k3=(k1+k2)-k0;
      if((val[k0]*val[k3])>(val[k1]*val[k2])) {
        s1=0;
      } else {
        s1=1;
      }
      k1=marchingCubesTable[index][6];
      k2=marchingCubesTable[index][7];
      k0=marchingCubesTable[index][8];
      k3=(k1+k2)-k0;
      if((val[k0]*val[k3])>(val[k1]*val[k2])) {
        s2=0;
      } else {
        s2=1;
      }
      indx=2*indx0+s1*2+s2;
      break;
    }
    case 3: {
      int n=3, s[3];
      for(int i=0;i<3;i++) {
        int k1=marchingCubesTable[index][n++];
        int k2=marchingCubesTable[index][n++];
        int k0=marchingCubesTable[index][n++];
        int k3=(k1+k2)-k0;
        if((val[k0]*val[k3])>(val[k1]*val[k2])) {
          s[i]=0;
        } else {
          s[i]=1;
        }
      }
      indx=2*indx0+s[0]*4+s[1]*2+s[2];
      break;
    }
    case 6: {
      int n=3, s[6];
      for(int i=0;i<6;i++) {
        int k1=marchingCubesTable[index][n++];
        int k2=marchingCubesTable[index][n++];
        int k0=marchingCubesTable[index][n++];
        int k3=(k1+k2)-k0;
        if((val[k0]*val[k3])>(val[k1]*val[k2])) {
          s[i]=0;
        } else {
          s[i]=1;
        }
      }
      indx=2*indx0+s[0]*32+s[1]*16+s[2]*8+s[3]*4+s[4]*2+s[5];
      break;
    }
    default: {
      cerr<<"Error in SubTableIndex()!"<<endl;
      exit(1);
    }
  }
  return indx;
}
#endif

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualVolGraphTool::DrawDelimiter::DrawDelimiter(bool do_best,
VirtualVolGraphTool *gt) : VirtualGLWindow::DrawDelimiter(do_best,
dynamic_cast<WINDOW_TYPE*>(gt->getWindow())) {
}
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualVolGraphTool::DrawDelimiter::DrawDelimiter(
VirtualVolGraphTool *gt,const VirtualGLWindow::ClipPlane *cp) : 
VirtualGLWindow::DrawDelimiter(dynamic_cast<WINDOW_TYPE*>(gt->getWindow()),
cp) {
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualVolGraphTool::drawPolygon(const WINDOW_TYPE::COLOR_TYPE &xc,
const int &edge_number,const Vector3* const vertex) {
  Vector3 centroid;
  for(int m0=0; m0<edge_number; m0++) {
    centroid+=vertex[m0];
  }
  centroid/=double(edge_number);

  for (int m=0;m<edge_number-1;m++) {
    colorTriangle(xc,centroid,vertex[m],vertex[m+1]);
  }
  colorTriangle(xc,centroid,vertex[edge_number-1],vertex[0]);
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualVolGraphTool::drawSurfaceWithinCube(double frac,
const Vector3 intersected[3][2][2],const double value[2][2][2]) {
  int entry[18];
  int index=cubeIndex(value);
  int n=marchingCubesTable[index][1];
  if (n>0) {
    for (int i=0;i<18;i++)
      entry[i]=marchingCubesTable[index][i];
  } else if (n<0) {
    int index0=subTableIndex(index, value);
    for (int i=0;i<18;i++) {
      entry[i]=marchingCubesSubTable[index0][i];
    }
  }

  if (n) {
    int current=6, k=1;
    int i=current;
    int size=entry[k];
    WINDOW_TYPE::COLOR_TYPE xc=color_array.getColor(frac);
    const Vector3 *inters=&intersected[0][0][0];
    while (size>0) {
      Vector3 vertex[8];
      for (int m=0; m<size; m++) {
        int n=entry[i++];
        vertex[m]=inters[n];
      }
      drawPolygon(xc,size,vertex);

      current+=size;
      size=entry[++k];
    }
  }
}
#endif

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool::VolGraphTool(const char *name,const Vector3 &l,
const Vector3 &h,const Palette &pal,bool use_lighting,
double maxwinsize,bool rotate_best,bool show_bounding_box) :
VirtualVolGraphTool(name,l,h,pal),win(0),raster_file(0) {
  Vector3 wc=windowCoords(high)*1.0625;
  win=OPERATOR_NEW WINDOW_TYPE(name,wc,use_lighting,
    maxwinsize,maxwinsize,rotate_best,false,true,show_bounding_box);
  rescale(l,h);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VolGraphTool::~VolGraphTool() {
  if (win!=0) delete win; win=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VolGraphTool::eventLoop(bool loop_only_if_escaped) {
  win->newPage();
  win->expose();
  win->eventLoop(loop_only_if_escaped);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VolGraphTool::createRasterFile(bool compress) {
  CHECK_TEST(raster_file==0);
  char label[LABEL_LENGTH];
  snprintf(label,LABEL_LENGTH,"%s.raster",getWindow()->getName());
  if (compress) {
    char label2[LABEL_LENGTH];
    snprintf(label2,LABEL_LENGTH,"gzip -c > %s.gz",label);
    raster_file=popen(label2,"w");
  } else {
    raster_file=fopen(label,"w");
  }
  getWindow()->createRaster(raster_file);
}

#include <numeric>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VolGraphTool::printOn(ostream &os) const {
  os << "VolGraphTool: maxlen = " << maxlen
     << "\n\twin = " << (void*) win
//   << "\n\tlow = " << low
//   << "\n\thigh = " << high
//   << "\n\tcenter = " << center 
     << endl;
  if (win!=0) win->printOn(os);
//color_array.printOn(os);
  VirtualGraphTool::printOn(os);
}
#endif
