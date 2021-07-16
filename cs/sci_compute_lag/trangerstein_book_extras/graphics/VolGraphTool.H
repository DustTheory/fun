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
#ifndef VOL_GRAPHTOOL_H
#define VOL_GRAPHTOOL_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/VolGraphTool.H,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#if (SPACEDIM>1)

#include <iostream>
#include <X11/Intrinsic.h>
#include <GL/glx.h>

#include "VGT.H"
#include "Vector3.H"

#ifdef USE_GTK
#include "GDKColorArray.H"
#include "GTKGLWindow.H"
#else
#include "XColorArray.H"
#include "GLWindow.H"
#endif

class VirtualVolGraphTool : public VirtualGraphTool {
  public:
#ifdef USE_GTK
    typedef GDKColorArray COLOR_ARRAY_TYPE;
    typedef GTKGLWindow WINDOW_TYPE;
#else
    typedef XColorArray COLOR_ARRAY_TYPE;
    typedef GLWindow WINDOW_TYPE;
#endif
    class DrawDelimiter;
    friend class DrawDelimiter;
  protected:
    COLOR_ARRAY_TYPE color_array;
    Vector3 low,high,center;
    double maxlen;
#if (SPACEDIM==3)
    double iso_surface_frac;
#endif
    Vector3 windowCoords(const Vector3 &user_coord) const {
      Vector3 c=(user_coord-center)/maxlen;
#if (SPACEDIM==2)
      c[2]*=maxlen/(high[2]-low[2]);
#endif
      return c;
    }
    Vector3 userCoords(const Vector3 &window_coord) const {
      Vector3 u=center+window_coord*maxlen;
#if (SPACEDIM==2)
      u[2]=center[2]+window_coord[2]*(high[2]-low[2]);
#endif
      return u;
    }
  public:
    VirtualVolGraphTool(const char *name,const Vector3 &l,
    const Vector3 &h,const Palette &pal) : 
    VirtualGraphTool(name),color_array(pal),maxlen(ZERO)
#if (SPACEDIM==3)
    ,iso_surface_frac(UNDEFINED)
#endif
    {
      rescale(l,h);
    }
    virtual void eventLoop(bool=FALSE)=0;
    virtual void rescale(const Vector3 &l,const Vector3 &h) {
      CHECK_TEST(l<h);
      low=l;
      high=h;
      center=(low+high)*HALF;
#if (SPACEDIM==2)
      maxlen=max(high[0]-low[0],high[1]-low[1]);
#endif
#if (SPACEDIM==3)
      maxlen=(high-low).max();
#endif
    }
    Vector3 getLength() const { return high-low; }
    const Vector3& getLow() const { return low; }
    const Vector3& getHigh() const { return high; }
    virtual double getLocation(const AxisClipPlane *cp) const=0;
    virtual bool monochrome() const { return FALSE; }
    virtual VirtualWindow* getWindow()=0;
    virtual const VirtualWindow* getWindow() const=0;
    double getPixelSize() const { 
      return maxlen/static_cast<double>(getWindow()->getWidth());
    }
    virtual const COLOR_ARRAY_TYPE* getColorArray() const=0;
    virtual void setfgColor(const char*)=0;
    virtual void newPage()=0;
    virtual void expose()=0;
#if (SPACEDIM==3)
    virtual WINDOW_TYPE::SlideBar* slideBar() const=0;
    virtual bool hasSlideBar() const=0;
    virtual void deleteSlideBar() const=0;
    virtual WINDOW_TYPE::ToggleButtonArray* toggleButtonArray(int ns=2) 
      const=0;
    virtual bool hasToggleButtonArray() const=0;
    virtual void deleteToggleButtonArray() const=0;
    virtual void setSurfaceNormalDirection(
      WINDOW_TYPE::SURFACE_NORMAL_DIRECTION)=0;
#endif
    virtual int numberSurfaces() const=0;
    virtual void colorRect(const WINDOW_TYPE::COLOR_TYPE&,
      const Vector3&,const Vector3&,int)=0;
    virtual void colorRect(const char*,const Vector3&,const Vector3&,
      int)=0;
    virtual void colorRect(double,const Vector3&,const Vector3&,int)=0;
#if (SPACEDIM==2)
    virtual void colorTriangle(double,const Vector3&,const Vector3&,
      const Vector3&)=0;
#endif
#if (SPACEDIM==3)
    virtual double isoSurfaceFrac() const {
      return iso_surface_frac;
    }
    void setIsoSurfaceFrac(double frac) { 
      CHECK_TEST(frac!=UNDEFINED)
      iso_surface_frac=max(ZERO,min(ONE,frac)); 
    }
    void drawPolygon(const WINDOW_TYPE::COLOR_TYPE&,const int&,
      const Vector3* const);
    virtual void colorTriangle(const WINDOW_TYPE::COLOR_TYPE&,
      const Vector3&,const Vector3&,const Vector3&)=0;
    virtual void drawVector(const Vector3&,const Vector3&)=0;
    virtual void drawVector(const Vector3&,const Vector3&,
      const AxisClipPlane*)=0;
#endif
    virtual void drawBox(const Vector3&,const Vector3&) const=0;
    virtual void drawLine(const Vector3&,const Vector3&) const=0;
    virtual bool isDrawing()=0;
    virtual void initialize(const VirtualGLWindow::PlotObj*)=0;
    virtual void terminatePlotObj()=0;
    virtual void plotClipPlanes() const=0;
    virtual void debugClipPlanes() const=0;
#if (SPACEDIM==3)
    virtual void drawSurfaceWithinCube(double,const Vector3[3][2][2],
      const double[2][2][2]);
#endif
    virtual int getMouse(double &x,double &y,double &z) const=0;
    virtual void printOn(ostream &os=cerr) const {
      os << "VirtualVolGraphTool:" << endl;
      os << "\tlow = " << low << endl;
      os << "\thigh = " << high << endl;
      os << "\tcenter = " << center << endl;
      os << "\tmaxlen = " << maxlen << endl;
      VirtualGraphTool::printOn(os);
    }
};

class VirtualVolGraphTool::DrawDelimiter : 
public VirtualGLWindow::DrawDelimiter {
  public:
    DrawDelimiter(bool,VirtualVolGraphTool*);
#if (SPACEDIM==3)
    DrawDelimiter(VirtualVolGraphTool*,
      const VirtualGLWindow::ClipPlane*);
#endif
    ~DrawDelimiter() {;}
};

class VolGraphTool : public VirtualVolGraphTool {
  private:
    WINDOW_TYPE *win;
    FILE *raster_file;

    VolGraphTool();
  public:
    VolGraphTool(const char*,const Vector3&,const Vector3&,
      const Palette&,bool,double=HALF,bool=FALSE,bool=FALSE);
    virtual ~VolGraphTool();
    virtual void rescale(const Vector3 &l,const Vector3 &h) {
      Vector3 d=(h-l)*.015625;
      VirtualVolGraphTool::rescale(l-d,h+d);
      win->rescale(l-d,h+d);
    }
    void eventLoop(bool=FALSE);
    double getLocation(const AxisClipPlane *cp) const {
      double frac=cp->getLocation();
      unsigned int dir=cp->getDirection();
      return center[dir]+frac*maxlen;
    }
    WINDOW_TYPE* getWindow() { return win; } 
    const WINDOW_TYPE* getWindow() const { return win; } 
    const COLOR_ARRAY_TYPE* getColorArray() const { 
      return &color_array; 
    }
    virtual void putString(const Vector3 *user_coords,const char *str) 
    const {
      Vector3 window_coords=windowCoords(*user_coords);
      win->putString(&window_coords,str);
    }
    void setfgColor(const char *str) {
      win->setfgColor(str);
    }
    void newPage() { win->clearDisplayLists(); win->newPage(); }
    void expose() { win->expose(); }
    int numberSurfaces() const { return win->numberSurfaces(); }
#if (SPACEDIM==3)
    WINDOW_TYPE::SlideBar* slideBar() const { return win->slideBar(); }
    bool hasSlideBar() const { return win->hasSlideBar(); }
    void deleteSlideBar() const { win->deleteSlideBar(); }
    WINDOW_TYPE::ToggleButtonArray* toggleButtonArray(int ns=2) const {
      return win->toggleButtonArray(ns);
    }
    bool hasToggleButtonArray() const {
      return win->hasToggleButtonArray();
    }
    void deleteToggleButtonArray() const {
      win->deleteToggleButtonArray();
    }
    virtual void setSurfaceNormalDirection(
    WINDOW_TYPE::SURFACE_NORMAL_DIRECTION snd) {
      win->setSurfaceNormalDirection(snd);
    }
#endif

//  these must be inside scope of a DrawDelimiter:
    void colorRect(const WINDOW_TYPE::COLOR_TYPE &xc,const Vector3 &l,
    const Vector3 &h,int dir) {
      Vector3 winlo=windowCoords(l);
      Vector3 winhi=windowCoords(h);
      win->colorRect(&xc,&winlo,&winhi,&dir);
    }
    void colorRect(const char *str,const Vector3 &l,const Vector3 &h,
    int dir) {
      win->colorRect(str,windowCoords(l),windowCoords(h),dir);
    }
    void colorRect(double frac,const Vector3 &l,const Vector3 &h,int dir){
      WINDOW_TYPE::COLOR_TYPE xc=color_array.getColor(frac);
      Vector3 winlo=windowCoords(l);
      Vector3 winhi=windowCoords(h);
      win->colorRect(&xc,&winlo,&winhi,&dir);
    }
    void colorTriangle(double frac,const Vector3 &u0,
    const Vector3 &u1,const Vector3 &u2) {
      frac=max(ZERO,min(ONE,frac));
      WINDOW_TYPE::COLOR_TYPE xc=color_array.getColor(frac);
      Vector3 w0=windowCoords(u0);
      Vector3 w1=windowCoords(u1);
      Vector3 w2=windowCoords(u2);
      win->colorTriangle(&xc,&w0,&w1,&w2);
    }
    void drawVector(const Vector3 &base,const Vector3 &head) {
      Vector3 wb=windowCoords(base);
      Vector3 wh=windowCoords(head);
      win->drawVector(&wb,&wh);
    }
#if (SPACEDIM==3)
    void colorTriangle(const WINDOW_TYPE::COLOR_TYPE &xc,
    const Vector3 &u0,const Vector3 &u1,const Vector3 &u2) {
      Vector3 w0=windowCoords(u0);
      Vector3 w1=windowCoords(u1);
      Vector3 w2=windowCoords(u2);
      win->colorTriangle(&xc,&w0,&w1,&w2);
    }
    void drawVector(const Vector3 &base,const Vector3 &head,
    const AxisClipPlane *cp) {
      Vector3 wb=windowCoords(base);
      Vector3 wh=windowCoords(head);
      win->drawVector(&wb,&wh,cp);
    }
#endif
    void drawBox(const Vector3 &l,const Vector3 &h) const {
      Vector3 winlo=windowCoords(l);
      Vector3 winhi=windowCoords(h);
      GLenum mode=GL_RENDER;
      GLenum type=GL_QUADS;
      win->drawBox(&winlo,&winhi,&mode,&type);
    }
    void drawLine(const Vector3 &start,const Vector3 &end) const {
      Vector3 wins=windowCoords(start);
      Vector3 wine=windowCoords(end);
      win->drawLine(&wins,&wine);
    }
    bool isDrawing() { return win->isDrawing(); }
    void initialize(const VirtualGLWindow::PlotObj *po) {
      win->initialize(po);
    }
    void terminatePlotObj() { win->terminatePlotObj(); }
    void plotClipPlanes() const {
      win->plotClipPlanes();
    }
    void debugClipPlanes() const {
      win->debugClipPlanes();
    }
    bool buttonIsPressed(Vector3 &user_coord) const {
      bool button_pressed=FALSE;
      Vector3 window_coord;
      if (win) button_pressed=win->buttonIsPressed(window_coord);
      user_coord=userCoords(window_coord);
      return button_pressed;
    }
    bool buttonIsPressed(double &x,double &y,double &z) const {
      Vector3 user_coord;
      bool button_pressed=buttonIsPressed(user_coord);
      x=user_coord[0]; y=user_coord[1]; z=user_coord[2];
      return button_pressed;
    }
    int getMouse(Vector3 &user_coord) const {
      Vector3 window_coord=win->getMouse();
      user_coord=userCoords(window_coord);
      return 2;
    }
    int getMouse(double &x,double &y,double &z) const {
      Vector3 user_coord;
      int button=getMouse(user_coord);
      x=user_coord[0]; y=user_coord[1]; z=user_coord[2];
      return button;
    }
    void createRasterFile(bool);
    void writeRasterFile() {
      CHECK_POINTER(raster_file);
      getWindow()->writeRaster(raster_file);
    }
    void closeRasterFile(bool compress) {
      CHECK_POINTER(raster_file);
      getWindow()->closeRaster();
      if (compress) pclose(raster_file);
      else fclose(raster_file);
      raster_file=0;
    }
    void openRasterFile(FILE *rf,bool compress) {
      CHECK_TEST(raster_file==0);
      raster_file=rf;
      getWindow()->openRaster(raster_file);
    }
    void readRasterFile() {
      CHECK_POINTER(raster_file);
      getWindow()->readRaster(raster_file);
    }
    void printOn(ostream& = cerr) const;
};
#endif
#endif
