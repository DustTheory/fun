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
// /usr/include/gtk-2.0
// /usr/share/doc/gtk+-devel-1.2.10/examples
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/GTKWindow.C,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#ifdef USE_GTK
#include "GTKWindow.H"
#include <gtk/gtk.h>
#include <string.h>
//#include "Tracer.H"

#define NTSC_WIDTH 600
#define NTSC_HEIGHT 420
#define LABEL_LENGTH 80
#define SEQUENCE_LENGTH 4

# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKWindow::GTKWindow(const char *str,GTKColormap *cm,double rw,double rh,
bool mono) : name(0),widget(0),frame(0),drawa(0),cmap(0),
gtk_colormap(cm),display(0),screen(0),gc(0),buffer_gc(0),buffer_map(0),
mouseX(0.),mouseY(0.),getting_mouse(false),monitor_mouse(false),
which_button_pressed(0),image(0),raster_number(0),xpm_number(0),
width(0),height(0),rwidth(ZERO),rheight(ZERO),x0(0),y0(0),
grey_scale(mono),xborder(0),yborder(0) {
  if (str) {
    name=OPERATOR_NEW_BRACKET(char,strlen(str)+1);
    strcpy(name,str);
  } else {
    name=OPERATOR_NEW_BRACKET(char,1);
    name[0]='\0';
  }

  ASSERT(cm!=0);
  if (cm) {
    display=cm->getDisplay();
    screen=cm->getScreen();
    cmap=cm->getColormap();
  } else {
    GdkDisplayManager *display_manager=gdk_display_manager_get();
    display=gdk_display_manager_get_default_display(display_manager);
    screen=gdk_display_get_default_screen(display);
    cmap=gdk_screen_get_default_colormap(screen);
  }
  ASSERT(display);

  gint dw=gdk_screen_get_width(screen);
  gint dh=gdk_screen_get_height(screen);
  int dm=min(dw,dh);
  if (rw>ONE) rw=ONE;
  if (rh>ONE) rh=ONE;
  if (rw<=ZERO || rh<=ZERO) {
    width=NTSC_WIDTH;
    height=NTSC_HEIGHT;
  } else {
    width=int(rw*dm);
    height=int(rh*dm);
  }
  xborder=width/20;
  yborder=height/20;
  rwidth=double(width);
  rheight=double(height);

  widget=gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(widget),name);
  gtk_container_set_reallocate_redraws(GTK_CONTAINER(widget),true);
  g_signal_connect(G_OBJECT(widget),"delete_event",
    G_CALLBACK(gtk_main_quit),0);

  if (!gtk_colormap) depth=1;
  else depth=gdk_visual_get_best_depth();
  grey_scale= (depth==1);

  frame=gtk_frame_new(0);
  gtk_container_add(GTK_CONTAINER(widget),frame);
  gtk_widget_show(frame);

  drawa=gtk_drawing_area_new();
  gtk_widget_set_size_request(drawa,width,height);
  gtk_container_add(GTK_CONTAINER(frame),drawa);
  gtk_widget_show(drawa);

// signals in /usr/include/gtk-2.0/gtk/gtksignal.h
// event masks in /usr/include/gtk-2.0/gdk/gdkevents.h
  gtk_widget_set_events(drawa, GDK_EXPOSURE_MASK |
    GDK_BUTTON_MOTION_MASK | 
    GDK_BUTTON_PRESS_MASK | 
    GDK_BUTTON_RELEASE_MASK);
  gtk_signal_connect_after(GTK_OBJECT(drawa),"realize",
    G_CALLBACK(onRealizeGTKWindow),static_cast<gpointer>(this));
  gtk_signal_connect(GTK_OBJECT(drawa),"expose_event",
    G_CALLBACK(onExposeEventGTKWindow),static_cast<gpointer>(this));
  gtk_signal_connect(GTK_OBJECT(drawa),"button_press_event",
    G_CALLBACK(onButtonPressGTKWindow),static_cast<gpointer>(this));
  gtk_signal_connect(GTK_OBJECT(drawa),"button_release_event",
    G_CALLBACK(onButtonReleaseGTKWindow),static_cast<gpointer>(this));
  gtk_signal_connect(GTK_OBJECT(drawa),"motion_notify_event",
    G_CALLBACK(onMouseMotionGTKWindow),static_cast<gpointer>(this));

  gtk_widget_show(widget);

  x0=0;
  y0=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKWindow::~GTKWindow() {
  delete [] name; name=0;
  if (gc) g_object_unref(gc); gc=0;
  if (buffer_gc) g_object_unref(buffer_gc); buffer_gc=0;
  if (buffer_map) g_object_unref(buffer_map); buffer_map=0;
  if (widget) gtk_widget_destroy(widget); widget=0;
  frame=0;
  drawa=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void onRealizeGTKWindow(GtkObject *o,gpointer d) {
  GtkWidget *w=reinterpret_cast<GtkWidget*>(o);
  GTKWindow *gtkw=reinterpret_cast<GTKWindow*>(d);
  gtkw->gc=gdk_gc_new(w->window);
  ASSERT(gtkw->gc);

  gtkw->buffer_map=
    gdk_pixmap_new(w->window,gtkw->width,gtkw->height,-1);

  gtkw->buffer_gc=gdk_gc_new(gtkw->buffer_map);
  ASSERT(gtkw->buffer_gc);
  gdk_window_clear(w->window);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::refresh(int x,int y,int w,int h) {
  gdk_draw_drawable(drawa->window,buffer_gc,buffer_map,x,y,x,y,w,h);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::expose() {
  refresh(0,0,width,height);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gint onExposeEventGTKWindow(GtkObject *o,GdkEventExpose *event,
gpointer d) {
  GtkWidget *w=reinterpret_cast<GtkWidget*>(o);
  GTKWindow *gtkw=reinterpret_cast<GTKWindow*>(d);
  gdk_draw_drawable(w->window,gtkw->buffer_gc,gtkw->buffer_map,
    event->area.x,event->area.y,event->area.x,event->area.y,
    event->area.width,event->area.height);
  return true;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gint onButtonPressGTKWindow(GtkObject *o,GdkEventButton *event,
gpointer d) {
  GtkWidget *w=reinterpret_cast<GtkWidget*>(o);
  GTKWindow *gtkw=reinterpret_cast<GTKWindow*>(d);
  if (gtkw->getting_mouse || gtkw->monitor_mouse) {
    gtkw->mouseX=event->x;
    gtkw->mouseY=event->y;
    gtkw->which_button_pressed=event->button;
    gtk_main_quit();
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gint onButtonReleaseGTKWindow(GtkObject *o,GdkEventButton *event,
gpointer d) {
  GTKWindow *gtkw=reinterpret_cast<GTKWindow*>(d);
  if (gtkw->monitor_mouse) {
    gtkw->monitor_mouse=false;
    gtkw->mouseX=event->x;
    gtkw->mouseY=event->y;
    gtk_main_quit();
  }
}  
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gint onMouseMotionGTKWindow(GtkObject *o,GdkEventMotion *event,
gpointer d) {
  GTKWindow *gtkw=reinterpret_cast<GTKWindow*>(d);
  if (gtkw->monitor_mouse) {
    gtkw->mouseX=event->x;
    gtkw->mouseY=event->y;
    gtk_main_quit();
  }
}  
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::setfgColor(const double *frac) const {
  setfgColor(*frac,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::setbgColor(const double *frac) const {
  setbgColor(*frac,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#if (SPACEDIM == 1)
void GTKWindow::setfgColor(const int *ic,const int *im) const {
  setfgColor(*ic,*im,buffer_gc);
}
#else
void GTKWindow::setfgColor(const int *ic,const int *im) const {
  double val=double(*ic+1)/double(*im+1);
  setfgColor(val,buffer_gc);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void  GTKWindow::setFont(const char *font) const {
  ASSERT(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GdkColor GTKWindow::allocColor(const char *colorname) {
  CHECK_POINTER(colorname);
  return gtk_colormap->allocColor(colorname);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//see http://maemo.org/lxr/source/gtk+/demos/gtk-demo/rotated_text.c
void gdk_rot_draw_string(GdkDrawable *drawable,GtkWidget *drawa,
gfloat angle,GdkGC *gc,gint x,gint y,const gchar *string) {
  while (angle<0) angle+=360.;
  while(angle>=360) angle-=360.;
  if (angle==0.) {
    GdkFont *text_font=gdk_font_load("8x13");
    gdk_draw_string(drawable,text_font,gc,x,y,string);
  } else {
    PangoContext *context=gtk_widget_create_pango_context(drawa);
    PangoLayout *layout = pango_layout_new(context);
    pango_layout_set_text(layout,string,-1);
    PangoFontDescription *desc=
      pango_font_description_from_string("8x13");
    pango_layout_set_font_description(layout,desc);
    pango_font_description_free(desc);
    PangoMatrix matrix=PANGO_MATRIX_INIT;
    pango_matrix_rotate(&matrix,angle);
    pango_context_set_matrix(context,&matrix);
    pango_layout_context_changed(layout);
    gdk_draw_layout(drawable,gc,x,y,layout);
    g_object_unref(layout);
    g_object_unref(context);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::putString(const double *px,const double *py,const char *str,
const double *angle) const {
  gfloat ang=static_cast<gfloat>(*angle);
  gdk_rot_draw_string(buffer_map,drawa,ang,buffer_gc,ixLoc(*px),
    iyLoc(*py),str);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::writeXBorder(const char *str) const {
  GdkFont *text_font=gdk_font_load("8x13");
  gdk_draw_string(buffer_map,text_font,buffer_gc,width/2,myY(yborder/2),
    str);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::writeYBorder(const char *str) const {
  GdkFont *text_font=gdk_font_load("8x13");
  gdk_rot_draw_string(buffer_map,drawa,90.,buffer_gc,xborder/2,
    myY(height/2),str);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::setfgColor(const char *colorname) const {
  GTKWindow *copy=const_cast<GTKWindow*>(this);
  copy->setfgColour(copy->allocColor(colorname));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::setbgColor(const char *colorname) const {
  GTKWindow *copy=const_cast<GTKWindow*>(this);
  copy->setbgColour(copy->allocColor(colorname));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int GTKWindow::pickColor(double frac) const {
  if (gtk_colormap) {
    int ifrac=static_cast<int>(
      frac*static_cast<double>(gtk_colormap->getNumberColors())-HALF);
    return max(0,min(gtk_colormap->getNumberColors()-1,ifrac));
  } else return 1;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::newPage(GdkPixmap *drawable,GdkGC *the_gc) const {
  setbgColour(gtk_colormap->getColor(0),drawable,the_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::newPage() {
  setbgColour(gtk_colormap->getColor(0),buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::drawLine(const double *px,const double *py,
const bool *puse_buffer) const {
  GTKWindow *copy=const_cast<GTKWindow*>(this);
  bool use_buffer=(puse_buffer ? *puse_buffer : true);
  int x1=ixLoc(*px); 
  int y1=iyLoc(*py);
  if (use_buffer) {
    gdk_draw_line(buffer_map,buffer_gc,
      x0,y0,x1,y1);
  } else {
    gdk_draw_line(drawa->window,gc, x0, y0, x1, y1 );
  }
  copy->x0 = x1; copy->y0 = y1;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::refreshFrom(const Buffer &xb) {
  gdk_draw_drawable(drawa->window,gc,xb.getPixmap(),0,0,0,0,
    width,height);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::flush(bool use_buffer) { 
  if (use_buffer) refresh(0,0,width,height);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//this should be called before createRaster
void GTKWindow::writePaletteName(FILE *os) {
  const char* pal_name=gtk_colormap->getPalette()->getName();
  int length=strlen(pal_name)+1;
  fwrite(reinterpret_cast<char*>(&length), sizeof(int), 1, os);
  fwrite( pal_name, sizeof(char), length, os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::createRaster(FILE *os) {
  CHECK_TEST(!image);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::writeRaster(FILE *os) {
  CHECK_TEST(!image);

  ASSERT(gdk_pixbuf_get_from_drawable(image,buffer_map,
    gtk_colormap->getColormap(),0,0,0,0,width,height));
  ASSERT(image);

  char label[LABEL_LENGTH];
  for (int j=0;j<LABEL_LENGTH;j++) label[j]=0;

  char labelc[SEQUENCE_LENGTH+1];
  labelc[SEQUENCE_LENGTH]='\0';
  unsigned int n=raster_number;
  for (int i=SEQUENCE_LENGTH-1;i>=0;i--,n/=10) labelc[i]='0'+n %10;

  snprintf(label,LABEL_LENGTH,"%s.raster.%s",name,labelc);
  GError **error=0;
  ASSERT(gdk_pixbuf_save(image,label,"png",error,0));
  g_object_unref(image); image=0;
  raster_number++;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::closeRaster() {
  CHECK_TEST(!image);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::openRaster(FILE *is) {
  CHECK_TEST(!image);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::readRaster(FILE *is) {
  CHECK_TEST(!image);

  char label[LABEL_LENGTH];
  for (int j=0;j<LABEL_LENGTH;j++) label[j]=0;
  snprintf(label,LABEL_LENGTH,"%s.raster.%d",name,raster_number);

  GError **error=0;
  image=gdk_pixbuf_new_from_file(label,error);
  ASSERT(image);
  gdk_draw_pixbuf(drawa->window,gc,image,0,0,0,0,
    width,height,GDK_RGB_DITHER_NONE,0,0);
  g_object_unref(image); image=0;
  raster_number++;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::writeXPM(const char *filename) {
  CHECK_TEST(!image);

  image=gdk_pixbuf_get_from_drawable(0,buffer_map,
    gtk_colormap->getColormap(),0,0,0,0,width,height);
  ASSERT(image);

  char label[LABEL_LENGTH];
  for (int j=0;j<LABEL_LENGTH;j++) label[j]=0;

  char labelc[SEQUENCE_LENGTH+1];
  labelc[SEQUENCE_LENGTH]='\0';
  unsigned int n=xpm_number;
  for (int i=SEQUENCE_LENGTH-1;i>=0;i--,n/=10) labelc[i]='0'+n %10;

  snprintf(label,LABEL_LENGTH,"%s.%s.png",filename,labelc);
  GError **error=0;
  ASSERT(gdk_pixbuf_save(image,label,"png",error,0));
  g_object_unref(image); image=0;

  xpm_number++;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//finds point where mouse button is pressed
guint GTKWindow::getMouse() {
  mouseX=mouseY=0.;
  getting_mouse=true;
  gtkMain();
  getting_mouse=false;
  return which_button_pressed;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int GTKWindow::getMouse(double &rx, double &ry) {
  guint wbp=getMouse();
  rx=(mouseX-static_cast<gdouble>(xborder))
    /(rwidth-static_cast<gdouble>(xborder));
  ry=static_cast<gdouble>(myY(static_cast<int>(mouseY)+yborder))
    /(rheight-static_cast<gdouble>(yborder));
  return wbp;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//should be called after getMouse, before button is released
void GTKWindow::monitorMouse() {
  monitor_mouse=true;
  gtkMain();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool GTKWindow::buttonIsPressed(double &rx,double &ry) {
  int x=0,y=0;
  monitorMouse();
  rx=(mouseX-static_cast<gdouble>(xborder))
    /(rwidth-static_cast<gdouble>(xborder));
  ry=static_cast<gdouble>(myY(static_cast<int>(mouseY)+yborder))
    /(rheight-static_cast<gdouble>(yborder));
  return monitor_mouse;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::colorRect(const double *ratio,const double *rw,
const double *rh,const double *rx,const double *ry) const {
  colorRect(*ratio,*rw,*rh,*rx,*ry,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::colorRect(const char *color_name,double rw,double rh,
double rx,double ry) {
  colorRect(color_name,rw,rh,rx,ry,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::colorPolygon(int npts,GdkPoint *p,GdkPixmap *drawable,
GdkGC *the_gc) const {
  gdk_draw_polygon(drawable,the_gc,true,p,npts);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::colorPolygon(int color,int npts,GdkPoint *p, 
GdkPixmap *drawable,GdkGC *the_gc) const {
  color=max(0,min(gtk_colormap->getNumberColors(),color));
  setfgColour(gtk_colormap->getColor(color),the_gc);
  gdk_draw_polygon(drawable,the_gc,true,p,npts);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::colorPolygon(int npts,const double *xfrac,
const double *yfrac,GdkPixmap *drawable,GdkGC *the_gc) const {
  GdkPoint *p=OPERATOR_NEW_BRACKET(GdkPoint,npts);
  for (int i=0;i<npts;i++) {
    double rx=min(ONE,max(ZERO,xfrac[i]));
    double ry=min(ONE,max(ZERO,yfrac[i]));
    p[i].x=static_cast<gint>(ixLoc(rx));
    p[i].y=static_cast<gint>(iyLoc(ry));
  }
  colorPolygon(npts,p,drawable,the_gc);
  delete [] p;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::colorPolygon(double ratio,int npts,const double *xfrac,
const double *yfrac,GdkPixmap *drawable,GdkGC *the_gc) const {
  int num_cmap_colors=gtk_colormap->getNumberColors();
  int color=0;
  if (ratio>=ZERO && ratio<=ONE) color=2+int(ratio*(num_cmap_colors-3));
  else if (ratio<ZERO) color=2;
  else color=num_cmap_colors-1;
  setfgColour(gtk_colormap->getColor(color),the_gc);
  colorPolygon(npts,xfrac,yfrac,drawable,the_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::colorPolygon(const double *ratio,const int *npts,
const double *xfrac,const double *yfrac) const {
  colorPolygon(*ratio,*npts,xfrac,yfrac,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::colorPolygon(const int *npts,const double *xfrac,
const double *yfrac) const {
  colorPolygon(*npts,xfrac,yfrac,buffer_map,buffer_gc);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::colorPolygon(int color,int npts,GdkPoint *p) 
const {
  parent->colorPolygon(color,npts,p,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::colorPolygon(const double *ratio,
const int *npts,const double *xfrac,const double *yfrac) const {
  parent->colorPolygon(*ratio,*npts,xfrac,yfrac,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::colorPolygon(const int *npts,
const double *xfrac,const double *yfrac) const {
  parent->colorPolygon(*npts,xfrac,yfrac,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::printOn(ostream& os) const {
  os << "GTKWindow: name = " << name << endl; 

  int num_screens=gdk_display_get_n_screens(display);
  os << "number of screens available = " << num_screens << "\n"
     << "default screen = " << gdk_display_get_default_screen(display)
     << "\n"
     << "screen used by this window = " << screen << endl;
  printScreenInfo(os,display,screen);

  os << "x0 = " << x0 << "\n"
     << "y0 = " << y0 << "\n"
     << "depth = " << depth << endl;
//os << "cmap = " << cmap << "\n"
//   << "xcolormap = " << (XColormap*) xcolormap << endl;
//os << "double_buffer = " << PRINT_bool(double_buffer) << endl;
  VirtualWindow12::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::printScreenInfo(ostream &os,GdkDisplay *d,GdkScreen *s){
  os << "scrn " << gdk_screen_get_number(s) << " is " 
     << gdk_screen_get_width(s)
     << " by " << gdk_screen_get_height(s) << " pixels\n"
     << endl;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKWindow::Buffer::Buffer(GTKWindow *xw) : VirtualWindow12(),parent(xw),
buffer_map(0),buffer_gc(0),x0(0),y0(0) {
  buffer_map=parent->makePixmap();
  buffer_gc=gdk_gc_new(buffer_map);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKWindow::Buffer::~Buffer() { }
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::setLineWidth(const int *lw) const { 
  gdk_gc_set_line_attributes(buffer_gc,(*lw+1)/2,GDK_LINE_SOLID,
    GDK_CAP_ROUND,GDK_JOIN_ROUND);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::putString(const double *px,const double *py,
const char *str,const double * /*angle*/) const {
  GdkFont *text_font=gdk_font_load("10x20");
  gdk_draw_string(buffer_map,text_font,buffer_gc,
    parent->ixLoc(*px),parent->iyLoc(*py),str);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::writeXBorder(const char *str) const {
  GdkFont *text_font=gdk_font_load("10x20");
  gdk_draw_string(buffer_map,text_font,buffer_gc,
    parent->width/2,parent->myY(parent->yborder/2),str);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::writeYBorder(const char *str) const {
  GdkFont *text_font=gdk_font_load("10x20");
  gdk_draw_string(buffer_map,text_font,buffer_gc,
    parent->xborder/2,parent->myY(parent->height/2),str);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::colorRect(const double *ratio,const double *rw,
const double *rh,const double *rx,const double *ry) const {
  parent->colorRect(*ratio,*rw,*rh,*rx,*ry,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::colorRect(const char *color_name,double rw,
double rh,double rx,double ry) {
  parent->colorRect(color_name,rw,rh,rx,ry,buffer_map,buffer_gc);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKWindow::Buffer::printOn(ostream& os) const {
  os << "GTKWindow::Buffer : parent = " << parent << endl;
  VirtualWindow12::printOn(os);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKWindow::QuitButton::QuitButton() : win(0),button(0),label(0) {
  win=gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(win),"Click to Continue");

  button=gtk_button_new();
  gtk_container_add(GTK_CONTAINER(win),button);
  g_signal_connect(G_OBJECT (button),"clicked",
    G_CALLBACK (gtk_main_quit),0);
  gtk_widget_show(button);

  label=gtk_label_new("Continue");
  gtk_container_add(GTK_CONTAINER(button),label);
  gtk_widget_show(label);

  gtk_widget_show(win);
  GTKWindow::gtkMain();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKWindow::QuitButton::~QuitButton() {
  if (win) gtk_widget_destroy(win); 
  win=0;
  button=0;
  label=0;
}
#endif
