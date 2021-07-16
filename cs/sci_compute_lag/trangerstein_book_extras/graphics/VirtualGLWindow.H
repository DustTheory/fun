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
#ifndef VIRTUAL_GL_WINDOW_H
#define VIRTUAL_GL_WINDOW_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/VirtualGLWindow.H,v 1.1 2009/08/20 17:31:47 johnt Exp $"

#if (SPACEDIM>1)
# include <GL/gl.h>
# if defined (__digital__) && defined (__unix__)
#endif
//#include "MemoryDebugger.H"
//#include "Tracer.H"
#include "Types.H"
#include "Vector3.H"
#include "VirtualWindow.H"

class VirtualGLWindow;
class MyDisplayList {
  private:
    GLuint base;
    GLsizei number;
  public:
    MyDisplayList();
    virtual ~MyDisplayList();
    GLuint getBase() const { return base; }
    void call() const;
    virtual void printOn(ostream &os) const {
      os << "\tMyDisplayList:" << endl;
    }
};

class AxisClipPlane;
class VirtualGLWindow : public VirtualWindow {
  private:
    VirtualWindow* getWindow() { CHECK_POINTER(0) return 0; }
    const VirtualWindow* getWindow() const {CHECK_POINTER(0) return 0;}
  public:
    class ClipPlane;
    class DrawDelimiter;
    class PlotObj;
    friend class VirtualGLWindow::DrawDelimiter;
  protected:
    char *name;
    const PlotObj *plot_obj;
    bool rotate_best;
    AxisClipPlane *clip_plane[3][2];
#if (SPACEDIM==3)
    int number_surfaces;
#endif
    Vector3 bound;// low corner=-bound, high corner = bound
    Vector3 user_low,user_high,user_center;
    double user_maxlen;
    bool inside_draw_delimiter;

    VirtualGLWindow(const VirtualGLWindow&) { 
      OBSOLETE("copy constructor disallowed"); 
    }
//  using (x,H-y) was inconsistent with GL's gluProject
    Vector3 userCoords(const Vector3 &window_coord) const {
      Vector3 u=user_center+window_coord*user_maxlen;
#if (SPACEDIM==2)
      u[2]=user_center[2]+window_coord[2]*(user_high[2]-user_low[2]);
#endif
      return u;
    }
    GLdouble detMatrix(GLdouble*) const;
    void clearClipPlanes();
  public:
    VirtualGLWindow(const char*,const Vector3&,bool=false);
    virtual ~VirtualGLWindow();
    virtual int numberSurfaces() const { 
      return D_CHOOSE(1,1,number_surfaces); 
    }
    virtual void startDrawing(bool) { inside_draw_delimiter=true; }
    virtual bool isDrawing() const { return inside_draw_delimiter; }
    virtual void finishDrawing() { inside_draw_delimiter=false; }
    void disableClipPlanes() const;
    void enableClipPlanes() const;
    double clipPlaneLocation(unsigned int d,HAND h) const;
    virtual const char* getName() const { return name; }
  
    virtual void ginit();
    virtual void rescale(const Vector3&,const Vector3&);
    virtual void plotObject(bool)=0;
    virtual void initialize(const VirtualGLWindow::PlotObj*)=0;
    virtual void plotClipPlanes() const=0;
    virtual void terminatePlotObj()=0;
    virtual void setLineWidth(const int*) const=0;
    virtual void setfgColor(const char*) const=0;
    virtual void setbgColor(const char*) const=0;
    virtual void newPage()=0;
    virtual void flush(bool=true)=0;
    virtual void colorRect(const COLOR_TYPE*,const Vector3*,
      const Vector3*,const int*) const=0;
    virtual void colorRect(const char *str,const Vector3 &l,
      const Vector3 &h,int d) const=0;
    virtual void drawLine(const Vector3*,const Vector3*) const=0;
    virtual void putString(const Vector3*,const char*)=0;
    virtual void drawBox(const Vector3*,const Vector3*,const GLenum*,
      const GLenum*) const=0;
#if (SPACEDIM==3)
    virtual void drawVector(const Vector3*,const Vector3*) const=0;
    virtual void drawVector(const Vector3*,const Vector3*,
      const AxisClipPlane*) const=0;
#endif
    virtual Vector3 getMouse()=0;

    virtual void printOn(ostream &os=cerr) const;
    void printMatrix(GLdouble*) const;
    void debugClipPlanes() const;
}; 

class VirtualGLWindow::ClipPlane {
  private:
    GLenum plane_number;
    GLdouble equation_coefs[4];
#if (SPACEDIM==3)
    MyDisplayList *display_list;
#endif
  protected:
    void fixCoefs() { 
      glClipPlane(plane_number,equation_coefs);
    }
  public:
//  p. 107: equation_coefs' object_coords >=0
    ClipPlane(GLenum pn,GLdouble x,GLdouble y,GLdouble z,GLdouble w);
    virtual ~ClipPlane() {
#if (SPACEDIM==3)
      if (display_list!=0) delete display_list; display_list=0;
#endif
    }
    GLdouble* equationCoefs() { return equation_coefs; } 
    const GLdouble* equationCoefs() const { return equation_coefs; } 
    bool isEnabled() const { return glIsEnabled(plane_number); }
    void enable() { fixCoefs(); glEnable(plane_number); }
    void disable() { glDisable(plane_number); }
#if (SPACEDIM==3)
    GLuint getDisplayList() const { return display_list->getBase(); }
    void callPlotList() const { display_list->call(); }
#endif
    void printOn(ostream &os) const {
      os << "ClipPlane: plane_number = " << plane_number
         << " enabled = " << isEnabled() << endl;
      GLdouble coefs[4];
      glGetClipPlane(plane_number,coefs);
      os << "\tequation_coefs = <" << coefs[0] << "," << coefs[1] << ","
         << coefs[2] << "," << coefs[3] << ">" << endl;
    }
};

class VirtualGLWindow::DrawDelimiter {
  private:
    bool is_enabled,clipping;
    VirtualGLWindow *win;
  public:
    DrawDelimiter(bool,VirtualGLWindow*);
#if (SPACEDIM==3)
    DrawDelimiter(VirtualGLWindow*,const ClipPlane*);
#endif
    ~DrawDelimiter();
};

class AxisClipPlane : public VirtualGLWindow::ClipPlane {
  private:
    unsigned int dir;
    HAND hand;
  public:
    AxisClipPlane(unsigned int d,HAND h,GLdouble val) : dir(d),hand(h),
      ClipPlane(GLenum(GL_CLIP_PLANE0+h+2*d),
                (d==0) ? (h==LEFT_HAND?ONE:-ONE) : ZERO,
                (d==1) ? (h==LEFT_HAND?ONE:-ONE) : ZERO,
                (d==2) ? (h==LEFT_HAND?ONE:-ONE) : ZERO,
                         (h==LEFT_HAND?-val:val)       ) {
    }
    ~AxisClipPlane() {}
    unsigned int getDirection() const { return dir; }
    HAND getHand() const { return hand; }
    GLdouble getLocation() const { 
      return -equationCoefs()[3]/equationCoefs()[dir]; 
    }
    void setLocation(GLdouble val) { 
      equationCoefs()[3]=-val*equationCoefs()[dir]; 
      fixCoefs();
    }
    void printOn(ostream &os) const {
      os << "AxisClipPlane: dir = " << dir << ", hand = " << hand
         << "\tlocation = " << getLocation() << endl; 
      VirtualGLWindow::ClipPlane::printOn(os);
    }
};

class VirtualGLWindow::PlotObj {
  public:
    virtual ~PlotObj() {}
    virtual bool drawSurface() const = 0;
    virtual void plot(bool) const = 0;
#if (SPACEDIM==3)
    virtual void plot(bool,const AxisClipPlane*) const = 0;
#endif
    virtual void printOn(ostream &os=cout) const {
      os << "GLWindow::PlotObj" << endl;
    }
};

#endif
#endif
