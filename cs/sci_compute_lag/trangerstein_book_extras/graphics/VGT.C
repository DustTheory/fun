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
 // "$Header: /home/faculty/johnt/cvs/deal_new/graphics/VGT.C,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#include "VGT.H"
#include <algorithm>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
//#include <values.h>
//#include "Tracer.H"
#define LABEL_LENGTH 10
#define TICK_LENGTH .01

# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualGraphTool::VirtualGraphTool(const char *str) : name(0) {
  name = OPERATOR_NEW_BRACKET(char,strlen(str) + 1);
  strcpy(name,str);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool::printOn(ostream &os) const {
  os << "VirtualGraphTool: " << this 
     << "\n\tname = " << name << endl;
  ISLListNode::printOn(os);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//x and y are in user coordinates
//window functions are called in terms of fractions of window width
//  and height
VirtualGraphTool12::VirtualGraphTool12(const char *str,const char *xl,
const char *yl) : VirtualGraphTool(str),raster_file(0),
black_and_white(false),xlabel(0),ylabel(0) {
  xlabel=OPERATOR_NEW_BRACKET(char,strlen(xl)+1);
  ylabel=OPERATOR_NEW_BRACKET(char,strlen(yl)+1);
  strcpy(xlabel,xl);
  strcpy(ylabel,yl);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
VirtualGraphTool12::~VirtualGraphTool12() {
  if (xlabel) { delete [] xlabel; xlabel=0; }
  if (ylabel) { delete [] ylabel; ylabel=0; }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::delayedConstructor(double xxlo,double xxhi,
double yylo,double yyhi,bool keep_aspect_ratio,bool mono) {
  rescale(xxlo,xxhi,yylo,yyhi);
  initGT(mono);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::delayedConstructor(double ww,double wh,
bool mono) {
  rescale(ZERO,ONE,ZERO,ONE);
  initGT(mono);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::initGT(bool mono) {
  black_and_white=mono || getWindow()->monochrome();
  getWindow()->newPage();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::rescale(double xxlo,double xxhi,double yylo,
double yyhi){
  CHECK_TEST(xxhi>xxlo);
  CHECK_TEST(xxlo>-LARGE);
  CHECK_TEST(xxhi<LARGE);
  CHECK_TEST(yylo>-LARGE);
  CHECK_TEST(yyhi<LARGE);

  double dx=(xxhi-xxlo)*0.05;
  xlo=xxlo-dx;
  xhi=xxhi+dx;

  if (yyhi<=yylo) {
    yyhi=yylo+HALF;
    yylo-=HALF;
  }
  double dy=(yyhi-yylo)*0.05;
  ylo=yylo-dy;
  yhi=yyhi+dy;
  xlen=xhi-xlo;
  ylen=yhi-ylo;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool VirtualGraphTool12::buttonIsPressed(double &x, double &y) {
  double xfrac,yfrac;
  bool pressed = getWindow()->buttonIsPressed(xfrac,yfrac);
  if (pressed) {
    x = xlo + (xhi-xlo)*xfrac;
    y = ylo + (yhi-ylo)*yfrac;
  }
  return pressed;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int VirtualGraphTool12::getMouse(double &x, double &y) {
  int button=-1;
  double xfrac,yfrac;
  button = getWindow()->getMouse(xfrac,yfrac);
  x = xlo + (xhi-xlo)*xfrac;
  y = ylo + (yhi-ylo)*yfrac;
  return button;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::createRasterFile(bool compress) {
  CHECK_TEST(raster_file==0);
  char rn[LENGTH_NAME];
  snprintf(rn,LENGTH_NAME,"%s.raster",getWindow()->getName());
  if (compress) {
    char cmd[LENGTH_NAME];
    snprintf(cmd,LENGTH_NAME,"gzip -c > %s.gz",rn);
    raster_file=popen(cmd,"w");
  } else raster_file=fopen(rn,"w");
//these correspond to GlobalMain::runMovie
  fwrite(reinterpret_cast<char*>(&xlo), sizeof(double), 1, raster_file);
  fwrite(reinterpret_cast<char*>(&ylo), sizeof(double), 1, raster_file);
  fwrite(reinterpret_cast<char*>(&xhi), sizeof(double), 1, raster_file);
  fwrite(reinterpret_cast<char*>(&yhi), sizeof(double), 1, raster_file);

  getWindow()->writePaletteName(raster_file);

//this corresponds to code in VGT12::openRaster
  getWindow()->createRaster(raster_file);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::openRasterFile(FILE *rf,bool compress) {
  CHECK_TEST(raster_file==0);
  raster_file=rf;
  getWindow()->openRaster(raster_file);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::writeRasterFile() {
  CHECK_POINTER(raster_file);
  getWindow()->writeRaster(raster_file);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::readRasterFile() {
  CHECK_POINTER(raster_file);
  getWindow()->readRaster(raster_file);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::closeRasterFile(bool compress) {
  CHECK_POINTER(raster_file);
  getWindow()->closeRaster();
  if (compress) pclose(raster_file);
  else fclose(raster_file);
  raster_file=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double ticSpacing(double low,double high) { 
  double sig_digit=pow(10.,static_cast<int>(floor(log10(high-low)))); 
  double digit_count=floor(high/sig_digit)-ceil(low/sig_digit);
  if (digit_count>10.) return 2.*sig_digit;
  else if (digit_count>=5.) return sig_digit;
  else if (digit_count>=2.) return 0.5*sig_digit;
  else return 0.25*sig_digit;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::drawAxes() {
  double fudge=0.05*(xhi-xlo)/1.1;
  double xxlo=xlo+fudge;
  double xxhi=xhi-fudge;
  double xtic=ticSpacing(xxlo,xxhi);
  int xfirst_tic=static_cast<int>(ceil(xlo/xtic));
  int xlast_tic=static_cast<int>(floor(xhi/xtic));
  double xaxis=max(xtic*static_cast<double>(xfirst_tic),
                 min(ZERO,xtic*static_cast<double>(xlast_tic)));

  fudge=0.05*(yhi-ylo)/1.1;
  double yylo=ylo+fudge;
  double yyhi=yhi-fudge;
  double ytic=ticSpacing(yylo,yyhi);
  ytic=max(ytic,-TWO*ylo/static_cast<double>(INT_MAX));
  ytic=max(ytic,TWO*yhi/static_cast<double>(INT_MAX));
  int yfirst_tic=static_cast<int>(ceil(ylo/ytic));
  int ylast_tic=static_cast<int>(floor(yhi/ytic));
  double yaxis=max(ytic*static_cast<double>(yfirst_tic),
                 min(ZERO,ytic*static_cast<double>(ylast_tic)));

  movePen(xlo,yaxis);
  drawLine(xhi,yaxis);

  char label[LABEL_LENGTH];
  for (int ix=xfirst_tic;ix<=xlast_tic;ix++) {
//  TRACER_CALL(t,"VirtualGraphTool12::drawAxes x axis");
    double x=static_cast<double>(ix)*xtic;
    movePen(x,yaxis-TICK_LENGTH*ylen);
    drawLine(x,yaxis+TICK_LENGTH*ylen);
    snprintf(label,LABEL_LENGTH,"%g",x);
    putString(x,yaxis-0.03e0*ylen,label,0.);
  }

  movePen(xaxis,ylo);
  drawLine(xaxis,yhi);

  for (int iy=yfirst_tic;iy<=ylast_tic;iy++) {
    double y=static_cast<double>(iy)*ytic;
    movePen(xaxis-TICK_LENGTH*xlen,y);
    drawLine(xaxis+TICK_LENGTH*xlen,y);
    snprintf(label,LABEL_LENGTH,"%g",y);
    putString(xaxis-0.03e0*xlen,y,label,90.);
  }
  getWindow()->writeXBorder(xlabel);
  getWindow()->writeYBorder(ylabel);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::colorPolygon(int npts,double *xpos,double *ypos) {
  double *xfrac=OPERATOR_NEW_BRACKET(double,npts);
  double *yfrac=OPERATOR_NEW_BRACKET(double,npts);
  for (int i=0;i<npts;i++) {
    xfrac[i]=(xpos[i]-xlo)/(xhi-xlo);
    yfrac[i]=(ypos[i]-ylo)/(yhi-ylo);
  }
  getWindow()->colorPolygon(&npts,xfrac,yfrac);
  delete [] yfrac;
  delete [] xfrac;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::colorPolygon(double ratio,int npts,double *xpos,
double *ypos) {
  double *xfrac=OPERATOR_NEW_BRACKET(double,npts);
  double *yfrac=OPERATOR_NEW_BRACKET(double,npts);
  for (int i=0;i<npts;i++) {
    xfrac[i]=(xpos[i]-xlo)/(xhi-xlo);
    yfrac[i]=(ypos[i]-ylo)/(yhi-ylo);
  }
  getWindow()->colorPolygon(&ratio,&npts,xfrac,yfrac);
  delete [] yfrac;
  delete [] xfrac;
}
#if (SPACEDIM>1)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::colorQuad(double ratio,double *xpos,double *ypos) {
  double xfrac[4],yfrac[4];
  for (int i=0;i<4;i++) {
    xfrac[i]=(xpos[i]-xlo)/(xhi-xlo);
    yfrac[i]=(ypos[i]-ylo)/(yhi-ylo);
  }
  int four=4;
  getWindow()->colorPolygon(&ratio,&four,xfrac,yfrac);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void VirtualGraphTool12::printOn(ostream &os) const {
  os << "VirtualGraphTool12: " << this << "\n"
     << "xlo = " << xlo << "\n"
     << "ylo = " << ylo << "\n"
     << "xhi = " << xhi << "\n"
     << "yhi = " << yhi << "\n"
     << "xlen = " << xlen << "\n"
     << "ylen = " << ylen << "\n" 
     << "black_and_white = " << black_and_white 
     << endl;
  VirtualGraphTool::printOn(os);
}
