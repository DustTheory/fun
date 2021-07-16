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
#ifndef GTKGLWINDOW_H
#define GTKGLWINDOW_H "$Header: /home/faculty/johnt/cvs/graphics/GTKGLWindow.H,v 1.2 2005/08/04 18:20:12 johnt Exp $"
#ifdef USE_GTK
#if (SPACEDIM>1)
#include <iostream>
#include <gdk/gdkkeysyms.h>
#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <GL/gl.h>
#if defined (__digital__) && defined (__unix__)
#else
#include <GL/glut.h>
#endif
#include <GL/glx.h>

#include "Const.H"
#include "GTKColormap.H"
#include "Line.H"
//#include "MemoryDebugger.H"
#include "NumPtr.H"
#include "Pair.H"
#include "Palette.H"
#include "Quaternion.H"
#include "ToggleCallback.H"
//#include "Tracer.H"
#include "Vector3.H"
#include "VGT.H"
#include "VirtualGLWindow.H"

// middle button: clip
// left button: rotate
// right button: cursor position
// extra widgets: level surface selection, colormap editor

class GTKGLWindow : public VirtualGLWindow {
  private:
  public:
    class AttributeDelimiter;
    class CullDelimiter;
    class DisablerDelimiter;
    class FeedbackRecord;
    class FeedbackDelimiter;
    class HitRecord;
    class LineDelimiter;
    class ListDelimiter;
    class MatrixDelimiter;
    class PrimitiveDelimiter;
    class QuadDelimiter;
    class SelectDelimiter;
    class StackCount;
#if (SPACEDIM==3)
    class SlideBar;
    class ToggleButtonArray;
    enum SURFACE_NORMAL_DIRECTION{INCREASING_FUNCTION,
      DECREASING_FUNCTION,INCREASING_OR_DECREASING_FUNCTION};
#endif
    class TriangleDelimiter;
    friend gboolean onRealizeGTKGLWindow(GtkWidget*,gpointer);
    friend gboolean onButtonPressEventGTKGLWindow(GtkWidget*,
      GdkEventButton*,gpointer);
    friend gboolean onButtonReleaseEventGTKGLWindow(GtkWidget*,
      GdkEventButton*,gpointer);
    friend gboolean onMotionNotifyEventGTKGLWindow(GtkWidget*,
      GdkEventMotion*,gpointer);
  private:
    GtkWidget *window,*vbox,*drawing_area;
    GdkGLDrawable *gldrawable;
    GdkGLContext *glcontext;

    GdkDisplay *display;
    GdkScreen *screen;
    GdkGC *gc;
    GdkPixbuf *image;
    int raster_number;

#if (SPACEDIM==3)
    SlideBar *slide_bar;
    ToggleButtonArray *toggle_button_array;
    SURFACE_NORMAL_DIRECTION surface_normal_direction;
#endif
    
    unsigned int W,H;
    double rwidth,rheight;

    GLfloat radius; // radius of the boundary box
    
    unsigned int active_clip_direction;
    HAND active_clip_hand;
    Line closest_edge;

    Pair<GLdouble> wuvec[3];

    bool use_lighting,perspective,three_dimensional;
    bool button_pressed;
    bool break_flag;
    bool finding_mouse;
    bool best_list_exposed;
    bool show_bounding_box;
    bool found_escape;
    bool double_buffer;

    Vector3 cursor;

    Quaternion<double> curquat;
    Pair<GLdouble> old_screen_coord;
    bool new_model,mat_valid,moved,cursor_selected;
    GLdouble ModelMat[16], ProjMat[16];
    GLint viewport[4];
    
    GLuint fontbase;

    MyDisplayList *bounding_box_list,*mouse_list,*user_list; 
    MyDisplayList *best_list;
    Pair<gdouble> xyToScreen(gdouble x,gdouble y) const {
      return Pair<gdouble>(x,static_cast<gdouble>(H)-y);
    }
    Pair<GLdouble> screenCoords(const Vector3&);

    GLdouble screenToX(Pair<GLdouble> screen) const {
      return screen[0];
    }
//  using H-screen[1] was inconsistent with GL's gluProject
    GLdouble screenToY(Pair<GLdouble> screen) const {
      return rheight-screen[1];
    }
    Line objectCoords(Pair<GLdouble>);

    Vector3 clipCoord(const Vector3&);

//  see Open GL Ref Manual, p. 427
    Pair<GLdouble> screenToUnitCircle(Pair<GLdouble> screen) const {
      return Pair<GLdouble>(TWO*screen[0]/rwidth-ONE,
                            TWO*screen[1]/rheight-ONE);
    }
    void clearErrors() const;
    void markFaceLocation(const Pair<GLdouble>&);
    int textWidth(char *str) const {
//    int dir,ascent,descent;
      ASSERT(0);
      return 0;
    }
    void resize();
    void makeBoundingBox();
    void recalcModelView();

    void makeCurrent();
    
    void selectClipPlane(const Pair<GLdouble>&);
    void moveClipPlane(const Pair<GLdouble>&);
    
    void makeRasterFont();

    Vector3 faceLocation(const Pair<GLdouble>&);
    void crossHairs(const Vector3&);
    void rotateImage(const Pair<GLdouble>&);
    void moveCursor(const Pair<GLdouble>&);
    void updateMat();

    void printMatrix(GLdouble*) const;

  public:
    GTKGLWindow(const char*,const Vector3&,bool,double=HALF,double=HALF,
      bool=false,bool=false,bool=true,bool=false);
    ~GTKGLWindow();
    void ginit();
    void eventLoop(bool=false);
    GdkColor getColor(const char *str) const {
      GdkColor color;
      gint success=gdk_color_parse(str,&color);
      return color;
    }
    void rescale(const Vector3&,const Vector3&);
    int getWidth() const { return static_cast<int>(W); }
    int getHeight() const { return static_cast<int>(H); }

    void setfgColor(const GdkColor &gdkc) const {
      glColor3us(gdkc.red,gdkc.green,gdkc.blue);
    }
    void setfgColor(const char *str) const {
      GTKGLWindow *copy=const_cast<GTKGLWindow*>(this);
      copy->setfgColor(getColor(str));
    }
    void setfgColor(GLclampf r,GLclampf g,GLclampf b) const {
      glColor3f(r,g,b);
    }

    void setbgColor(GLclampf r,GLclampf g,GLclampf b,GLclampf a=0.)
    const {
      glClearColor(r,g,b,a);
    }
    void setbgColor(const GdkColor &gdkc) {
      setbgColor(gdkc.red,gdkc.green,gdkc.blue);
    }
    void setbgColor(const char *str) const {
      GTKGLWindow *copy=const_cast<GTKGLWindow*>(this);
      copy->setbgColor(getColor(str));
    }

    void setLineWidth(const int *lw) {
//    see glGetFloatv(GL_LINE_WIDTH_RANGE)
      glLineWidth(static_cast<GLfloat>(*lw));
    }
    void setLineWidth(const int *lw) const {
      glLineWidth(static_cast<GLfloat>(*lw));
    }
    void setPointSize(int ps) const {
//    see glGetFloatv(GL_POINT_SIZE_RANGE)
      glPointSize(static_cast<GLfloat>(ps));
    }
    void setPolygonMode(GLenum face=GL_FRONT_AND_BACK,
    GLenum draw_mode=GL_FILL,GLenum face_mode=GL_CCW) const {
      CHECK_TEST(face==GL_FRONT_AND_BACK || face==GL_FRONT ||face==GL_BACK);
      CHECK_TEST(draw_mode==GL_POINT || draw_mode==GL_LINE || draw_mode==GL_FILL);
      CHECK_TEST(face_mode==GL_CCW || face_mode==GL_CW);
      glPolygonMode(face,draw_mode);
      glFrontFace(face_mode);
    }
    void newPage() { 
      glClear(GL_COLOR_BUFFER_BIT); 
    }
    void flush(bool) { 
      glFlush(); 
    }
    void putString(const Vector3 *window_coords,const char *str) {
      glRasterPos3f((*window_coords)[0],(*window_coords)[1],
        (*window_coords)[2]);
      printString(str);
    }
    void drawBox(const Vector3*,const Vector3*,const GLenum*,
      const GLenum*) const;
    void expose();

    void drawLine(const Vector3*,const Vector3*) const;
    void colorRect(const GdkColor*,const Vector3*,const Vector3*,
      const int*) const;
    void colorRect(const char *str,const Vector3 &l,const Vector3 &h,
    int d) const {
      GdkColor gdkc;
      gint success=gdk_color_parse(str,&gdkc);
      colorRect(&gdkc,&l,&h,&d);
    }
    void colorTriangle(const GdkColor*,const Vector3*,const Vector3*,
      const Vector3*) const;
    void drawVector(const Vector3*,const Vector3*) const;
    void drawVector(const Vector3*,const Vector3*,const AxisClipPlane*) 
      const;

    void printString(const char*);

    GLuint getDisplayList(bool do_best) { 
      if (do_best || rotate_best) {
        if (!best_list) best_list=OPERATOR_NEW MyDisplayList();
        return best_list->getBase(); 
      }
      if (!user_list) user_list=OPERATOR_NEW MyDisplayList();
      return user_list->getBase(); 
    }
    void freeDisplayList(bool do_best) { 
      if (do_best || rotate_best) {
        if (best_list) delete best_list; best_list=0;
      } else {
        if (user_list) delete user_list; user_list=0;
      }
    }
    void clearDisplayLists() {
      freeDisplayList(false);
      freeDisplayList(true);
    }
    
    const Vector3& getCursor() const { return cursor; } 

    void initialize(const VirtualGLWindow::PlotObj*);
    void terminatePlotObj();
    void plotClipPlanes() const;

    void plotObject(bool);

    int numberSurfaces() const;
#if (SPACEDIM==3)
    GTKGLWindow::SlideBar* slideBar();
    bool hasSlideBar() const { return slide_bar!=0; }
    void deleteSlideBar();
    GTKGLWindow::ToggleButtonArray* toggleButtonArray(int=2);
    bool hasToggleButtonArray() const { 
      return toggle_button_array!=0; 
    }
    void deleteToggleButtonArray();
    void setSurfaceNormalDirection(SURFACE_NORMAL_DIRECTION snd) {
      surface_normal_direction=snd;
    }
#endif

    void startDrawing(bool do_best) {
      CHECK_TEST(!isDrawing())
      VirtualGLWindow::startDrawing(do_best);
      glNewList(getDisplayList(do_best || rotate_best),GL_COMPILE);
      if (use_lighting) {
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
      }
    }
    void finishDrawing() {
      CHECK_TEST(isDrawing())
      if (use_lighting) {
        glDisable(GL_LIGHT0);
        glDisable(GL_LIGHTING);
      }
      VirtualGLWindow::finishDrawing();
      glEndList();
    }

    Vector3 getMouse();
    bool buttonIsPressed(Vector3 &where) const {
      where=cursor; return button_pressed;
    }
    void createRaster(FILE*);
    void writeRaster(FILE*);
    void closeRaster();
    void openRaster(FILE*);
    void readRaster(FILE*);

    void printOn(ostream &os=cerr) const;
};

gboolean onConfigureEventGTKGLWindow(GtkWidget*,GdkEventConfigure*,
  gpointer);
gboolean onExposeEventGTKGLWindow(GtkWidget*,GdkEventExpose*,gpointer);
gboolean onButtonPressEventGTKGLWindow(GtkWidget*,GdkEventButton*,
  gpointer);
gboolean onButtonReleaseEventGTKGLWindow(GtkWidget*,GdkEventButton*,
  gpointer);
gboolean onMotionNotifyEventGTKGLWindow(GtkWidget*,GdkEventMotion*,
  gpointer);
gboolean onMapEventGTKGLWindow(GtkWidget*,GdkEventAny*,gpointer);
gboolean onUnmapEventGTKGLWindow(GtkWidget*,GdkEventAny*,gpointer);
gboolean onVisibilityNotifyEventGTKGLWindow(GtkWidget*,
  GdkEventVisibility*,gpointer);
gboolean onKeyPressEventGTKGLWindow(GtkWidget*,GdkEventKey*,gpointer);
void onClickBreakButtonGTKGLWindow(GtkWidget*,gpointer);
void onClickQuitButtonGTKGLWindow(GtkWidget*,gpointer);
void onClickedNoGTKGLWindow(GtkWidget*,gpointer);
void onClickedYesGTKGLWindow(GtkWidget*,gpointer);

// p. 360
class GTKGLWindow::HitRecord {
  private:
    GLuint number_names;
    const GLuint *names;
  public:
    inline HitRecord(const GLuint *pos) : number_names(pos[0]),
      names(pos+3) {;}
    inline virtual ~HitRecord() {;}
    inline GLuint bottomName() const { return names[0]; }
    inline GLuint topName() const { return names[number_names-1]; }
    inline HitRecord next() const { 
      return HitRecord(names+number_names); 
    }
    virtual void printOn(ostream &os) const {
      os << "\tGLWindow::HitRecord:" << endl;
    }
};

// p. 378
class GTKGLWindow::FeedbackRecord {
  private:
    GLfloat token;
    const GLfloat *vertices;
  public:
    static GLenum feedback_type;

    inline FeedbackRecord(const GLfloat *pos):token(pos[0]),
      vertices(pos+1){;}
    inline virtual ~FeedbackRecord() {;}
    inline const GLfloat* position() const { return &token; }
    inline GLfloat getToken() const { return token; }
    FeedbackRecord next() const { 
      int nvertices=0;
      int extraneous=0;
      if (token==GL_POINT_TOKEN || token==GL_BITMAP_TOKEN ||
      token==GL_DRAW_PIXEL_TOKEN || token==GL_COPY_PIXEL_TOKEN) {
        nvertices=1;
      } else if (token==GL_LINE_TOKEN || token==GL_LINE_RESET_TOKEN) {
        nvertices=2;
      } else if (token==GL_POLYGON_TOKEN) {
        OBSOLETE("what were we trying to do?")
      } else if (token==GL_PASS_THROUGH_TOKEN) extraneous=1;
      switch (feedback_type) {
        case GL_2D:
          return FeedbackRecord(vertices+2*nvertices+extraneous);
        case GL_3D:
          return FeedbackRecord(vertices+3*nvertices+extraneous);
        case GL_3D_COLOR:
        case GL_3D_COLOR_TEXTURE:
        case GL_4D_COLOR_TEXTURE:
        default:
          OBSOLETE("have to determine whether in RGBA or color-index");
          return FeedbackRecord(&token);
      }
    }
    virtual void printOn(ostream &os) const {
      os << "\tGLWindow::BeedbackRecord:" << endl;
    }
};

class GTKGLWindow::StackCount {
  private:
    static unsigned int count;
  public:
    StackCount() { count++; }
    virtual ~StackCount() { 
      CHECK_POSITIVE(count)
      count--;
    }
    virtual void printOn(ostream &os) const {
      os << "\tStackCount:" << endl;
    }
};

class GTKGLWindow::AttributeDelimiter : public GTKGLWindow::StackCount {
  public:
    AttributeDelimiter(GLbitfield mask) { glPushAttrib(mask); }
    virtual ~AttributeDelimiter() { glPopAttrib(); }
};

class GTKGLWindow::MatrixDelimiter : public GTKGLWindow::StackCount {
  public:
    MatrixDelimiter() { glPushMatrix(); }
    virtual ~MatrixDelimiter() { glPopMatrix(); }
};

// GL_COMPILE              = 4864 = 0x1300
// GL_COMPILE_AND_EXECUTE  = 4865 = 0x1301
// GL_RENDER   = 7168 = 0x1C00
// GL_FEEDBACK = 7169 = 0x1C01
// GL_SELECT   = 7170 = 0x1C02
class GTKGLWindow::PrimitiveDelimiter {
  private:
    GLint list_mode,render_mode;
  public:
    PrimitiveDelimiter(GLenum mode) {
      glGetIntegerv(GL_RENDER_MODE,&render_mode);
      glGetIntegerv(GL_LIST_MODE,&list_mode);
      if (render_mode==GL_RENDER && list_mode==GL_COMPILE) 
        glBegin(mode); 
    }
    ~PrimitiveDelimiter() { 
      if (render_mode==GL_RENDER && list_mode==GL_COMPILE) 
        glEnd(); 
    }
};

//eg, for cap=GL_DITHER
class GTKGLWindow::DisablerDelimiter {
  private:
    GLenum capability;
  public:
    DisablerDelimiter(GLenum c) : capability(c) {glDisable(capability);}
    ~DisablerDelimiter() { glEnable(capability); }
};

class GTKGLWindow::LineDelimiter : 
public GTKGLWindow::PrimitiveDelimiter {
  public:
    LineDelimiter() : PrimitiveDelimiter(GL_LINES) {;}
};

// GL_QUADS = 7
class GTKGLWindow::QuadDelimiter : 
public GTKGLWindow::PrimitiveDelimiter {
  public:
    QuadDelimiter() : PrimitiveDelimiter(GL_QUADS) {;}
};

class GTKGLWindow::TriangleDelimiter : 
public GTKGLWindow::PrimitiveDelimiter
{
  public:
    TriangleDelimiter() : PrimitiveDelimiter(GL_TRIANGLES) {;}
};

class GTKGLWindow::FeedbackDelimiter {
  public:
    FeedbackDelimiter() { glRenderMode(GL_FEEDBACK); }
    ~FeedbackDelimiter() { }
};

class GTKGLWindow::SelectDelimiter {
  public:
    SelectDelimiter() {
      glRenderMode(GL_SELECT);
      glInitNames();
      glPushName(UINT_MAX);
    }
    ~SelectDelimiter() { }
};

class GTKGLWindow::ListDelimiter {
  private:
    GLuint active_list;
  public:
    ListDelimiter(GLuint list,GLenum mode) : active_list(list) {
      CHECK_TEST(mode==GL_COMPILE || mode==GL_COMPILE_AND_EXECUTE);
      CHECK_TEST(glIsList(list));
      glNewList(list,mode); 
    }
    ~ListDelimiter() { 
      glEndList(); 
    }
};

#if (SPACEDIM==3)
class GTKGLWindow::SlideBar {
  private:
    double frac,var_min,var_max;
    GtkWidget *slide_bar;
  public:
    SlideBar(GtkWidget*,gint);
    void rescale(double,double);
    void setValue(double);
    double getValue() const { return var_min+frac*(var_max-var_min); }
    void setCallback(void (*)(GtkWidget*,gpointer),gpointer) ;
    void removeCallback(void (*)(GtkWidget*,gpointer),gpointer) const {}
    void printOn(ostream&) const;
};

class GTKGLWindow::ToggleButtonArray {
  private:
    GtkWidget *hbutton_box;
    NumPtr<GtkWidget*> toggle_button_array;
  public:
    ToggleButtonArray(GtkWidget*,int);
    int numberSurfaces () const { 
      return toggle_button_array.getNumber(); 
      ASSERT(0);
    }
    gboolean drawSurface(int n) const {
      return gtk_toggle_button_get_active(
        GTK_TOGGLE_BUTTON(toggle_button_array[n]));
    }
    void setCallback(void (*)(GtkWidget*,gpointer),gpointer) ;
    void removeCallback(void (*)(GtkWidget*,gpointer),gpointer) const {}
    void printOn(ostream&) const;
};
#endif
#endif
#endif
#endif
