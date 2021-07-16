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
#ifndef GTK_WINDOW_H
#define GTK_WINDOW_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/GTKWindow.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#ifdef USE_GTK
#include <iostream>
#include <gdk/gdkdrawable.h>
#include <gdk/gdkpixmap.h>
#include <gtk/gtkbutton.h>
#include <gtk/gtkdrawingarea.h>
#include <gtk/gtklabel.h>
#include <gtk/gtkmain.h>
#include <gtk/gtkwindow.h>
#include "GTKColormap.H"
#include "VirtualWindow.H"

class GTKWindow : public VirtualWindow12 {
  public:
    typedef GTKColormap COLOR_MAP_TYPE;
    class Buffer;
    friend class Buffer;
    class QuitButton;
    friend void onRealizeGTKWindow(GtkObject*,gpointer);
    friend gint onExposeEventGTKWindow(GtkObject*,GdkEventExpose*,
      gpointer);
    friend gint onButtonPressGTKWindow(GtkObject*,GdkEventButton*,
      gpointer);
    friend gint onButtonReleaseGTKWindow(GtkObject*,GdkEventButton*,
      gpointer);
    friend gint onMouseMotionGTKWindow(GtkObject*,GdkEventMotion*,
      gpointer);
  private:
    char *name;

    GtkWidget *widget;
    GtkWidget *frame;
    GtkWidget *drawa;
    GdkColormap *cmap;
    GTKColormap *gtk_colormap;
    GdkDisplay *display;
    GdkScreen *screen;
    GdkGC *gc;
    GdkGC *buffer_gc;
    GdkPixmap *buffer_map;
    gdouble mouseX,mouseY;
    bool getting_mouse,monitor_mouse;
    guint which_button_pressed;

    GdkPixbuf *image;
    unsigned int raster_number,xpm_number;

    unsigned int xborder,yborder;
    unsigned int width,height;
    double rwidth, rheight;
    unsigned int x0, y0;
    bool grey_scale;
    gint depth;

    int myY(int y) const { return height - 1 - y; }
    unsigned int ixLoc(double x) const { 
      return xborder+static_cast<unsigned int>(min(ONE,max(ZERO,x))
                    *static_cast<double>(width-xborder)); 
    }
    unsigned int iyLoc(double y) const { 
      return myY(yborder+static_cast<unsigned int>(min(ONE,max(ZERO,y))
                        *static_cast<unsigned int>(height-yborder))); 
    }

    GdkPixmap *makePixmap() {
      return gdk_pixmap_new(widget->window,width,height,-1);
    }
    GdkColor getForeground(GdkGC *the_gc) const {
      GdkGCValues gcvalues;
      gdk_gc_get_values(the_gc,&gcvalues);
      return gcvalues.foreground;
    }
    void setfgColour(const GdkColor &color,GdkGC *the_gc) const {
      gdk_gc_set_foreground(the_gc,&color);
    }
    void setfgColor(const char *colorname,GdkGC *the_gc){
      GdkColor color=allocColor(colorname);
      setfgColour(color,the_gc);
    }
    void setfgColour(const GdkColor &color) const {
      setfgColour(color,buffer_gc);
    }
    void setfgColor(double frac,GdkGC *the_gc) const {
      frac=max(ZERO,min(ONE,frac));
      int color=pickColor(frac);
      setfgColour(gtk_colormap->getColor(color),the_gc);
    }
#if (SPACEDIM == 1)
    void setfgColor(int ic,int,GdkGC *the_gc) 
    const {
      ic=max(1,min(gtk_colormap->getNumberColors(),ic));
      setfgColour(gtk_colormap->getColor(ic),the_gc);
    }
#else
    void setfgColor(int ic,int im,GdkGC *the_gc) 
    const {
      setfgColor(double(ic+1)/double(im+1),the_gc);
    }
#endif
    void setbgColour(const GdkColor &color,GdkPixmap *drawable,
    GdkGC *the_gc) const {
      GdkColor foreground=getForeground(the_gc);
      setfgColour(color,the_gc);
      gboolean filled=true;
      gdk_draw_rectangle(drawable,the_gc,filled,0,0,width,height);
      gdk_gc_set_foreground(the_gc,&foreground);
    }
    void setbgColor(const char *colorname,GdkPixmap *drawable,
    GdkGC *the_gc) {
      GdkColor color=allocColor(colorname);
      setbgColour(color,drawable,the_gc);
    }
    void setbgColour(const GdkColor &color) const {
      setbgColour(color,buffer_map,buffer_gc);
    }
    void setbgColor(double frac,GdkPixmap *drawable,GdkGC *the_gc) const {
      frac=max(ZERO,min(ONE,frac));
      setbgColour(gtk_colormap->getColor(pickColor(frac)),drawable,
        the_gc);
    }
    void newPage(GdkPixmap*,GdkGC*) const;
    void colorRect(int color,int w,int h,int x,int y,
    GdkPixmap *drawable,GdkGC *the_gc) const {
      color=max(0,min(color,gtk_colormap->getNumberColors()));
      const GdkColor &c=gtk_colormap->getColor(color);
      gdk_gc_set_foreground(the_gc,&c);
      gboolean filled=true;
      gdk_draw_rectangle(drawable,the_gc,filled,0,0,w,h);
    }
    void colorRect(double ratio,double rw,double rh,double rx,double ry,
    GdkPixmap *drawable,GdkGC *the_gc) const {
      ratio=max(ZERO,min(ONE,ratio));
      int color=1+int(ratio*(gtk_colormap->getNumberColors()-2));
      rw=min(ONE,max(ZERO,rw));
      rh=min(ONE,max(ZERO,rh));
      rx=min(ONE,max(ZERO,rx));
      ry=min(ONE,max(ZERO,ry));
      colorRect(color,static_cast<int>(rw*rwidth),
        static_cast<int>(rh*rheight),static_cast<int>(rx*rwidth),
        static_cast<int>(ry*rheight),drawable,the_gc);
    }
    void colorRect(const char *colorname,double rw,double rh,double rx,
    double ry,GdkPixmap *drawable,GdkGC *the_gc) {
      GdkColor color=allocColor(colorname);
      gdk_gc_set_foreground(the_gc,&color);
      gboolean filled=true;
      gdk_draw_rectangle(drawable,the_gc,filled,
        static_cast<gint>(rx*rwidth),static_cast<gint>(ry*rheight),
        static_cast<gint>(rw*rwidth),static_cast<gint>(rh*rheight));
    }
    void colorPolygon(int,int,GdkPoint*,GdkPixmap*,GdkGC*) const;
    void colorPolygon(int,GdkPoint*,GdkPixmap*,GdkGC*) const;
    void colorPolygon(double,int,const double*,const double*,GdkPixmap*,
      GdkGC*) const;
    void colorPolygon(int,const double*,const double*,GdkPixmap*,GdkGC*) 
      const;

    void monitorMouse();
    guint getMouse();

  public:
    GTKWindow(const char*,GTKColormap*,double=HALF,double=HALF,
      bool=false);
    ~GTKWindow();
    virtual GTKWindow* getWindow() { return this; }
    virtual const GTKWindow* getWindow() const { return this; }
    const char* getName() const { return name; }
    const gchar* getDisplayName() const { 
      return gdk_display_get_name(display);
    }
    bool doubleBuffer() const { return true; }
    int getWidth() const { return width; }
    int getHeight() const { return height; }

    bool buttonIsPressed(double&,double&);
    void newPage();
    void movePen(const double *px,const double *py) const {
      GTKWindow *copy=const_cast<GTKWindow*>(this);
      copy->x0 = ixLoc(*px);
      copy->y0 = iyLoc(*py);
    }
    void drawLine(const double*,const double*,const bool* =0) const;
    bool monochrome() const { return grey_scale; }
    const GdkDisplay* getDisplay() const { return display; }
    int getMouse(double&,double&);
    void colorRect(const double*,const double*,const double*,const double*,
      const double*) const;
    void colorRect(const char *name, double rw,double rh,double rx,double ry);
    void colorPolygon(const double*,const int*,const double*,
      const double*) const;
    void colorPolygon(const int*,const double*,const double*) const;
    void setLineWidth(const int *plw) const {
      gdk_gc_set_line_attributes(buffer_gc,(*plw+1)/2,GDK_LINE_SOLID,
        GDK_CAP_ROUND,GDK_JOIN_ROUND);
    }
    void expose();
    void setFont(const char*) const;
    GdkColor allocColor(const char*);
    void putString(const double *px,const double *py,const char *str,
      const double *angle) const;
    void writeXBorder(const char *str) const;
    void writeYBorder(const char *str) const;
    void setfgColor(const char*) const;
    void setfgColor(const double*) const;
    void setfgColor(const int*,const int*) const;
    void setbgColor(const char*) const;
    void setbgColor(const double*) const;
    int pickColor(double) const ;
    void refresh(int,int,int,int) ;
    void refreshFrom(const Buffer&);
    void flush(bool=true);
    const GdkColormap* getColormap() {return cmap;}

    void writePaletteName(FILE*);
    void createRaster(FILE*);
    void writeRaster(FILE*);
    void closeRaster();
    void openRaster(FILE*);
    void readRaster(FILE*);
    void writeXPM(const char*);
    void printOn(ostream& os = cerr) const;
    static void printScreenInfo(ostream&,GdkDisplay *d,GdkScreen *s);
    static void gtkInit(int argc,char* argv[]) {
#ifndef NO_THREAD
      g_thread_init(0);
      gdk_threads_init();
      gdk_threads_enter();
#endif
      gtk_init(&argc,&argv);
    }
    static void gtkMain() {
      gtk_main();
#ifndef NO_THREAD
      gdk_threads_leave();
#endif
    }
    void eventLoop() { gtkMain(); }
};
void onRealizeGTKWindow(GtkObject*,gpointer);
gint onExposeEventGTKWindow(GtkObject*,GdkEventExpose*,gpointer);
gint onButtonPressGTKWindow(GtkObject*,GdkEventButton*,gpointer);
gint onButtonReleaseGTKWindow(GtkObject*,GdkEventButton*,gpointer);
gint onMouseMotionGTKWindow(GtkObject*,GdkEventMotion*,gpointer);

class GTKWindow::Buffer : public VirtualWindow12 {
  friend class GTKWindow;
  private:
    GTKWindow*const parent;
    GdkPixmap *buffer_map;
    GdkGC *buffer_gc;
    int x0, y0;

    GdkPixmap *getPixmap() const {return buffer_map; }
    void setFont(const char*) const {;}
    void flush(bool) {;}
    bool buttonIsPressed(double &x,double &y) { x=y=ZERO; return false; }
    int getMouse(double &x,double &y) { x=y=ZERO; return -1; }
    void writePaletteName(FILE*) {;}
    void createRaster(FILE*) {;}
    void writeRaster(FILE*) {;}
    void closeRaster() {;}
    void openRaster(FILE*) {;}
    void readRaster(FILE*) {;}
  protected:
    Buffer(const Buffer&) : parent(0) { 
      OBSOLETE("copy constructor disallowed"); 
    }
  public:
    Buffer(GTKWindow*);
    virtual ~Buffer();
    GTKWindow::Buffer* getWindow() { return this; }
    const GTKWindow::Buffer* getWindow() const { return this; }
    const char* getName() const { return parent->getName(); }
    const char* getDisplayName() const { 
      return parent->getDisplayName(); 
    }
    int getWidth() const { return parent->getWidth(); }
    int getHeight() const { return parent->getHeight(); }
  
    const GdkDisplay* getDisplay() const { 
      return parent->getDisplay(); 
    }
    bool monochrome() const { return parent->monochrome(); }
    void setLineWidth(const int*) const;
    void setfgColour(const GdkColor &color) { 
      parent->setfgColour(color,buffer_gc); 
    }
    void setfgColor(const char *colorname) const {
      parent->setfgColor(colorname,buffer_gc);
    }
    void setfgColor(const double *pfrac) const { 
      parent->setfgColor(*pfrac,buffer_gc); 
    }
    void setfgColor(const int *ic,const int* im) const {
      parent->setfgColor(*ic,*im,buffer_gc);
    }
    void setbgColour(const GdkColor &color) {
      parent->setbgColour(color,buffer_map,buffer_gc);
    }
    void setbgColor(const char *colorname) const {
      parent->setbgColor(colorname,buffer_map,buffer_gc);
    }
    void setbgColor(const double *pfrac) const {
      parent->setbgColor(*pfrac,buffer_map,buffer_gc); 
    }
    void putString(const double*,const double*,const char*,
      const double *angle) const;
    void writeXBorder(const char *str) const;
    void writeYBorder(const char *str) const;
    void movePen(const double *px,const double *py) const {
      GTKWindow::Buffer *copy=const_cast<GTKWindow::Buffer*>(this);
      copy->x0=parent->ixLoc(*px); copy->y0=parent->iyLoc(*py);
    }
    void drawLine(const double *px,const double *py,const bool* =false) 
    const {
      GTKWindow::Buffer *copy=const_cast<GTKWindow::Buffer*>(this);
      int x1=parent->ixLoc(*px); int y1=parent->iyLoc(*py);
      gdk_draw_line(buffer_map,buffer_gc,x0,y0,x1,y1);
      copy->x0 = x1; copy->y0 = y1;
    }

    void newPage() { parent->newPage(buffer_map,buffer_gc); }
    void colorPolygon(int color,int npts,GdkPoint *p) const;
    void colorPolygon(const double*,const int*,const double*,
      const double*) const;
    void colorPolygon(const int*,const double*,const double*) const;
    void colorRect(const double*,const double*,const double*,const double*,
      const double*) const;
    void colorRect(const char*,double,double,double,double);
    int pickColor(double frac) const { return parent->pickColor(frac); }
    void printOn(ostream& os = cerr) const;
};

class GTKWindow::QuitButton {
  private:
    GtkWidget *win;
    GtkWidget *button;
    GtkWidget *label;
  public:
    QuitButton();
    ~QuitButton();
};

#endif
#endif
