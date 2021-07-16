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
#ifndef VIRTUAL_GRAPH_TOOL_H
#define VIRTUAL_GRAPH_TOOL_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/VGT.H,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#include <complex>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "Const.H"
#include "ISLList.H"
#include "VirtualWindow.H"

class VirtualColormap;
class VirtualGraphTool : public ISLListNode {
  private:
    char* name;
  protected:
    VirtualGraphTool(const VirtualGraphTool&) { 
      OBSOLETE("copy constructor disallowed"); 
    }
  public:
    VirtualGraphTool(const char*);
    virtual ~VirtualGraphTool() { 
      if (name) delete [] name; name=0; 
    }
    virtual const char* getName() const { return name; }
    virtual VirtualWindow* getWindow() { return 0; }
    virtual const VirtualWindow* getWindow() const { return 0; }
    virtual void newPage() { 
      VirtualWindow *w=getWindow();
      CHECK_POINTER(w)
      w->newPage(); 
    }
    virtual void setLineWidth(int lw) { 
      VirtualWindow *w=getWindow();
      CHECK_POINTER(w)
      w->setLineWidth(&lw); 
    }
    virtual void setfgColor(const char *color) {
      VirtualWindow *w=getWindow();
      CHECK_POINTER(w)
      w->setfgColor(color);
    }
    virtual void setbgColor(const char *color) {
      VirtualWindow *w=getWindow();
      CHECK_POINTER(w)
      w->setbgColor(color);
    }
    virtual void flush(bool use_buffer=true) { 
      VirtualWindow *w=getWindow();
      CHECK_POINTER(w)
      w->flush(use_buffer); 
    }
    virtual void createRasterFile(bool) {;}
    virtual void writeRasterFile() {;}
    virtual void closeRasterFile(bool) {;}
    virtual void openRasterFile(FILE*,bool) {;}
    virtual void readRasterFile() {;}
    virtual void printOn(ostream& os = cerr) const;
};

class VirtualGraphTool12 : public VirtualGraphTool {
  private:
  protected:
    FILE *raster_file;
    double xlo,ylo,xhi,yhi; //user coordinates for window corners
    double xlen,ylen;       //width and height in user coordinates
    bool black_and_white;
    char *xlabel,*ylabel;

    virtual void initGT(bool);
    VirtualGraphTool12(const VirtualGraphTool12 &g) : 
      VirtualGraphTool(g) { OBSOLETE("copy constructor disallowed"); }
  public:
    VirtualGraphTool12(const char *str,const char *xl,const char *yl);
    void delayedConstructor(double,double,double,double,bool =true,
      bool =false);
    void delayedConstructor(double =HALF,double =HALF,bool =false);
    virtual ~VirtualGraphTool12();
    virtual VirtualWindow12* getWindow() { return 0; }
    virtual const VirtualWindow12* getWindow() const { return 0; }
    void setFont(char *fn) { getWindow()->setFont(fn); }
    void putString(double x,double y,const char *str,double angle=0.) {
      double xx=(x-xlo)/xlen;
      double yy=(y-ylo)/ylen;
      getWindow()->putString(&xx,&yy,str,&angle);
    }
    const char* getXLabel() const { return xlabel; }
    const char* getYLabel() const { return ylabel; }
    double getLowX() const { return xlo; }
    double getLowY() const { return ylo; }
    double getHighX() const { return xhi; }
    double getHighY() const { return yhi; }
    double getPixelWidth() const { 
      return (xhi-xlo)/(1.1*static_cast<double>(getWindow()->getWidth()));
    }
    double getPixelHeight() const { 
     return (xhi-xlo)/(1.1*static_cast<double>(getWindow()->getHeight()));
    }
    bool monochrome() const { 
      return getWindow()->monochrome(); 
    }
    void setfgColor(const char *color) { 
      VirtualGraphTool::setfgColor(color); 
    }
    void setfgColor(int ic,int im) {
      getWindow()->setfgColor(&ic,&im);
    }

    virtual void rescale(double,double,double,double);
    virtual bool buttonIsPressed(double&,double&);
    virtual int getMouse(double&,double&);
    void movePen(double x,double y) {
      double xfrac=(x-xlo)/xlen; double yfrac=(y-ylo)/ylen;
      getWindow()->movePen(&xfrac,&yfrac);
    }
    void drawLine(double x,double y,bool use_buffer=true) {
      double xfrac=(x-xlo)/xlen; double yfrac=(y-ylo)/ylen;
      getWindow()->drawLine(&xfrac,&yfrac,&use_buffer);
    }
    void drawPlusGivenCorners(double xxlo,double yylo,double xxhi,double yyhi) {
      double xmid=HALF*(xxlo+xxhi),ymid=HALF*(yylo+yyhi);
      movePen(xxlo,ymid); drawLine(xxhi,ymid);
      movePen(xmid,yylo); drawLine(xmid,yyhi);
    }
    void drawPlusGivenCenter(double x,double y,double dx,double dy) {
      movePen(x-dx,y); drawLine(x+dx,y);
      movePen(x,y-dy); drawLine(x,y+dy);
    }
    void drawPlus(double x,double y,double dx=HALF) {
      drawPlusGivenCenter(x,y,dx,dx*(yhi-ylo)/(xhi-xlo));
    }
    void drawCrossGivenCorners(double xxlo,double yylo,double xxhi,double yyhi){
      movePen(xxlo,yylo); drawLine(xxhi,yyhi);
      movePen(xxhi,yylo); drawLine(xxlo,yyhi);
    }
    void drawCrossGivenCenter(double x,double y,double dx,double dy) {
      drawCrossGivenCorners(x-dx,y-dy,x+dx,y+dy);
    }
    void drawCross(double x,double y,double dx=HALF) {
      drawCrossGivenCenter(x,y,dx,dx*(yhi-ylo)/(xhi-xlo));
    }
    void drawDiamondGivenCorners(double xxlo,double yylo,double xxhi,
    double yyhi) {
      double xmid=HALF*(xxlo+xxhi),ymid=HALF*(yylo+yyhi);
      movePen(xxhi,ymid);
      drawLine(xmid,yyhi); drawLine(xxlo,ymid);
      drawLine(xmid,yylo); drawLine(xxhi,ymid);
    }
    void drawDiamondGivenCenter(double x,double y,double dx=HALF) {
      double xxhi=x+dx;
      double dy=dx*(yhi-ylo)/(xhi-xlo);
      movePen(xxhi,y);
      drawLine(x,y+dy); drawLine(x-dx,y);
      drawLine(x,y-dy); drawLine(xxhi,y);
    }
    void drawBoxGivenCorners(double xxlo,double yylo,double xxhi,double yyhi) {
      movePen(xxlo,yylo);
      drawLine(xxhi,yylo); drawLine(xxhi,yyhi);
      drawLine(xxlo,yyhi); drawLine(xxlo,yylo);
    }
    void drawBoxGivenCenter(double x,double y,double dx=HALF) {
      double dy=dx*(yhi-ylo)/(xhi-xlo);
      drawBoxGivenCorners(x-dx,y-dy,x+dx,y+dy);
    }
    void colorPolygon(int,double*,double*);
    void colorPolygon(double,int,double*,double*);
    void newPage() { getWindow()->newPage(); }
    void createRasterFile(bool);
    void writeRasterFile();
    void closeRasterFile(bool);
    void openRasterFile(FILE*,bool);
    void readRasterFile();
#if (SPACEDIM>1)
    virtual void drawVector(complex<double> base,complex<double> head,complex<double> arrow) {
      complex<double> position=base+head; movePen(real(base),imag(base));
      drawLine(real(position),imag(position));
      position=position+head*arrow;
      drawLine(real(position),imag(position));
    }
    virtual void colorQuad(double,double*,double*);
#endif
    void drawAxes();
    virtual void printOn(ostream& os = cerr) const;
};

double ticSpacing(double,double);
#endif
