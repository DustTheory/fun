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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/GTKGLWindow.C,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#include "GTKGLWindow.H"
#ifdef USE_GTK
#if (SPACEDIM>1)
#include <math.h>

#define NTSC_WIDTH 600
#define NTSC_HEIGHT 420

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#define LABEL_LENGTH 10
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <GL/gl.h>
#include <GL/glu.h>
#include <gdk/gdkglconfig.h>
#include <gtk/gtkwindow.h>

//#include "Debug.H"
#include "GTKWindow.H"
//#include "MyInline.H"
const long maxIntensity = 65535;

#define LABEL_LENGTH 10
#define SEQUENCE_LENGTH 4
#define BUFSIZE 512
#define DEFAULT_FONT "courier 12"

unsigned int GTKGLWindow::StackCount::count=0;
GLenum GTKGLWindow::FeedbackRecord::feedback_type=GL_3D;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::updateMat() {
  if (!mat_valid) {
    glGetDoublev(GL_MODELVIEW_MATRIX, ModelMat);
    glGetDoublev(GL_PROJECTION_MATRIX, ProjMat);
    glGetIntegerv(GL_VIEWPORT, viewport);
    mat_valid=true;
  }
  CHECK_POSITIVE(detMatrix(ModelMat))
  CHECK_POSITIVE(-detMatrix(ProjMat))
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Vector3 GTKGLWindow::clipCoord(const Vector3 &object_coord) {
  updateMat();
  GLdouble winx, winy, winz;
  GLint fake_viewport[4];
    fake_viewport[0]=0; fake_viewport[1]=0;
    fake_viewport[2]=2; fake_viewport[3]=2;
  gluProject(object_coord[0],object_coord[1],object_coord[2],
              ModelMat, ProjMat, fake_viewport, &winx, &winy, &winz );
  return Vector3(winx-ONE,winy-ONE,TWO*winz-ONE);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Vector3 GTKGLWindow::faceLocation(const Pair<GLdouble> &screen_coord) {
  double max_depth=DBL_MAX;
  bool is_interior=false;
  Vector3 best=bound;

  Pair<GLdouble> plo[3],phi[3];
  Vector3 vlo[3],vhi[3];
  Pair<GLdouble> plow=screenCoords(-bound);
  Pair<GLdouble> phigh=screenCoords(bound);
  for (int i=0;i<3;i++) {
    vlo[i]=-bound; vlo[i][i]= bound[i]; plo[i]=screenCoords(vlo[i]);
    vhi[i]= bound; vhi[i][i]=-bound[i]; phi[i]=screenCoords(vhi[i]);
  }
  for (int i=0;i<3;i++) {
    int i1=(i+1)%3,i2=(i+2)%3;
    Pair<GLdouble> tlo=
      closestInteriorPosition(plow,plo[i1],plo[i2],screen_coord);
    Vector3 object_coord=
      -bound+(vlo[i1]+bound)*tlo[0]+(vlo[i2]+bound)*tlo[1];
    Vector3 f=clipCoord(object_coord);
    if (min(tlo[0],tlo[1])>ZERO && max(tlo[0],tlo[1])<ONE) {
      if (!is_interior || f[2]>max_depth) {
        max_depth=f[2];
        best=object_coord;
      }
      is_interior=true;
    } else if (!is_interior && f[2]>max_depth) {
      max_depth=f[2];
      best=object_coord;
    }

    Pair<GLdouble> thi=
      closestInteriorPosition(phigh,phi[i1],phi[i2],screen_coord);
    object_coord=bound+(vhi[i1]-bound)*thi[0]+(vhi[i2]-bound)*thi[1];
    f=clipCoord(object_coord);
    if (min(thi[0],thi[1])>ZERO && max(thi[0],thi[1])<ONE) {
      if (!is_interior || f[2]>max_depth) {
        max_depth=f[2];
        best=object_coord;
      }
      is_interior=true;
    } else if (!is_interior && f[2]>max_depth) {
      max_depth=f[2];
      best=object_coord;
    }
  }
  return best;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::crossHairs(const Vector3 &object_coord) {
  if  (!mouse_list) mouse_list=OPERATOR_NEW MyDisplayList();
  { ListDelimiter ld(mouse_list->getBase(),GL_COMPILE);
    double fudge=radius/30.;
    char label[LABEL_LENGTH];
    Vector3 user_coord=userCoords(object_coord);
    for (int i=0;i<3;i++) {
      if (i==0) setfgColor("red");
      if (i==1) setfgColor("green");
      if (i==2) setfgColor("blue");
      Vector3 low=object_coord; low[i]=-bound[i];
      Vector3 high=object_coord; high[i]=bound[i];
      drawLine(&low,&high);
      for (int j=0;j<LABEL_LENGTH;j++) label[j]=0;
      snprintf(label,LABEL_LENGTH,"%g",user_coord[i]);
      label[LABEL_LENGTH-1]=0;
      Vector3 pos=low; 
      pos[i]=-bound[i]-fudge;
      glRasterPos3f(pos[0],pos[1],pos[2]);
      printString(label);
      pos[i]=bound[i]+fudge;
      glRasterPos3f(pos[0],pos[1],pos[2]);
      printString(label);
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// //usr/share/doc/gtk2-devel-2.4.14/examples/menu/menu.c
gboolean onRealizeGTKGLWindow(GtkWidget *w,gpointer d) {
  GTKGLWindow *glw=reinterpret_cast<GTKGLWindow*>(d);
  glw->ginit();
  return true;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void onClickBreakButtonGTKGLWindow(GtkWidget *w,gpointer d) {
  bool &break_flag=*reinterpret_cast<bool*>(d);
  break_flag=true;
//gtk_main_quit();
}

// /usr/share/doc/gtkglext-devel-1.0.6/examples/logo.c
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGLWindow::GTKGLWindow(const char *na,const Vector3 &b,bool ul,
double rw,double rh,bool rb,bool persp,bool td,bool sbb) : 
VirtualGLWindow(na,b,rb),
window(0),vbox(0),drawing_area(0),gldrawable(0),glcontext(0),
#if (SPACEDIM==3)
slide_bar(0),toggle_button_array(0),
surface_normal_direction(INCREASING_FUNCTION),
#endif
display(0),screen(0),gc(0),image(0),raster_number(0),W(0),H(0),
rwidth(ZERO),rheight(ZERO),radius(0.),active_clip_direction(0),
active_clip_hand(LEFT_HAND),closest_edge(-bound,Vector3(ONE,ZERO,ZERO)),
use_lighting(ul),perspective(persp),three_dimensional(td),
show_bounding_box(sbb),
button_pressed(false),break_flag(false),
finding_mouse(false),cursor(bound),new_model(true),mat_valid(false),
moved(false),cursor_selected(false),
fontbase(0),bounding_box_list(0),mouse_list(0),user_list(0),
best_list(0),best_list_exposed(true),found_escape(false),
double_buffer(true) {
  for (int i=0;i<16;i++) {
    ModelMat[i]=ProjMat[i]=0.;
  }
  for (int i=0;i<16;i+=5) {
    ModelMat[i]=ProjMat[i]=1.;
  }
        
  radius=norm(bound);

#if (SPACEDIM==2)
//y direction up, z direction normal to screen
  Pair<double> zero_pair;
  curquat=trackBall(zero_pair,zero_pair);
#else
//z direction up, y normal to screen
  curquat=Quaternion<double>(ONE,ZERO,ZERO,ONE);
#endif

  GdkGLConfigMode mode=static_cast<GdkGLConfigMode>(
    GDK_GL_MODE_RGB | GDK_GL_MODE_DEPTH | GDK_GL_MODE_DOUBLE);
  GdkGLConfig *glconfig=gdk_gl_config_new_by_mode(mode);
  if (!glconfig) {
    mode=static_cast<GdkGLConfigMode>(
      GDK_GL_MODE_RGB | GDK_GL_MODE_DEPTH);
    glconfig=gdk_gl_config_new_by_mode(mode);
  }
  ASSERT(glconfig);

  GdkDisplayManager *display_manager=gdk_display_manager_get();
  display=gdk_display_manager_get_default_display(display_manager);
  ASSERT(display);
  screen=gdk_display_get_default_screen(display);
  gint dw=gdk_screen_get_width(screen);
  gint dh=gdk_screen_get_height(screen);
  int dm=min(dw,dh);
  rw=min(ONE,rw);
  rh=min(ONE,rh);
  if (rw<=ZERO || rh<=ZERO) {
    W=NTSC_WIDTH;
    H=NTSC_HEIGHT;
  } else {
    W=int(rw*dm);
    H=int(rh*dm);
  }
  rwidth=static_cast<double>(W);
  rheight=static_cast<double>(H);

  window=gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(window),na);
  gtk_container_set_reallocate_redraws(GTK_CONTAINER(window),true);
  g_signal_connect(G_OBJECT(window),"destroy",G_CALLBACK(gtk_main_quit),
    0);
  g_signal_connect(G_OBJECT(window),"delete_event",
    G_CALLBACK(gtk_main_quit),0);

  vbox=gtk_vbox_new(false,0);
  gtk_container_add(GTK_CONTAINER(window),vbox);
  gtk_widget_show(vbox);

  GtkWidget *menubar=gtk_frame_new("");

  GtkWidget *hbuttonbox=gtk_hbutton_box_new();
  gtk_container_set_border_width(GTK_CONTAINER(hbuttonbox),1);
  gtk_container_add(GTK_CONTAINER(menubar),hbuttonbox);
  gtk_button_box_set_layout(GTK_BUTTON_BOX(hbuttonbox),
    GTK_BUTTONBOX_START);
  gtk_button_box_set_spacing(GTK_BUTTON_BOX(hbuttonbox),1);

  GtkWidget *quit_button=gtk_button_new_from_stock(GTK_STOCK_QUIT);
  gtk_container_add(GTK_CONTAINER(hbuttonbox),quit_button);
  g_signal_connect(G_OBJECT(quit_button),"clicked",
    G_CALLBACK(onClickQuitButtonGTKGLWindow),&break_flag);
  gtk_widget_show (quit_button);
  gtk_widget_show (hbuttonbox);
  gtk_widget_show (menubar);

  gtk_box_pack_start_defaults(GTK_BOX(vbox),menubar);

  drawing_area=gtk_drawing_area_new();
  gtk_widget_set_size_request(drawing_area,W,H);
  gtk_widget_set_gl_capability(drawing_area,glconfig,0,true,
    GDK_GL_RGBA_TYPE);

  gtk_widget_add_events(drawing_area, GDK_BUTTON_MOTION_MASK |
    GDK_BUTTON_PRESS_MASK | GDK_BUTTON_PRESS_MASK |
    GDK_BUTTON_RELEASE_MASK | GDK_VISIBILITY_NOTIFY_MASK);

  g_signal_connect_after (G_OBJECT (drawing_area), "realize",
    G_CALLBACK (onRealizeGTKGLWindow),(gpointer) this);
  g_signal_connect (G_OBJECT (drawing_area), "expose_event",
    G_CALLBACK (onExposeEventGTKGLWindow),(gpointer) this);
  g_signal_connect (G_OBJECT (drawing_area), "button_press_event",
    G_CALLBACK (onButtonPressEventGTKGLWindow), (gpointer) this);
  g_signal_connect (G_OBJECT (drawing_area), "button_release_event",
    G_CALLBACK (onButtonReleaseEventGTKGLWindow), (gpointer) this);
  g_signal_connect (G_OBJECT (drawing_area), "motion_notify_event",
    G_CALLBACK (onMotionNotifyEventGTKGLWindow), (gpointer) this);
  g_signal_connect(G_OBJECT (drawing_area), "key_press_event",
    G_CALLBACK (onKeyPressEventGTKGLWindow),
    (gpointer) &found_escape);

  gtk_box_pack_start_defaults(GTK_BOX(vbox),drawing_area);

  gtk_widget_show (drawing_area);
  gtk_widget_show (window);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::ginit() {
  gldrawable=gtk_widget_get_gl_drawable(drawing_area);
  double_buffer=gdk_gl_drawable_is_double_buffered(gldrawable);
  glcontext=gtk_widget_get_gl_context(drawing_area);

// /usr/share/doc/gtkglext-devel-1.0.6/examples/font.c:
  
  ASSERT(gdk_gl_drawable_gl_begin(gldrawable,glcontext));
  fontbase = glGenLists(128);
  PangoFontDescription *font_desc=
    pango_font_description_from_string(DEFAULT_FONT);
  PangoFont *font=gdk_gl_font_use_pango_font(font_desc,0,128,fontbase);
  if (!font) {
    cerr << "Unable to load font: " << DEFAULT_FONT <<endl;
    ASSERT(0);
  }
  pango_font_description_free (font_desc);

  if (use_lighting) {
    glShadeModel(GL_SMOOTH);
    GLfloat light0_specular[] = {1.0f, 1.0f, 1.0f, 1.0};
    GLfloat light0_diffuse[] = {1.0f, 1.0f, 1.0f, 1.0};
    GLfloat position[]={0.0f, 0.0f, 1.0, 0.0f};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, position);
    {
      const GLfloat mat_diffuse1[]={1.0f, 0.5f, 0.5f, 1.0f};
      glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse1);
      const GLfloat mat_diffuse2[]={0.5f, 1.0f, 0.5f, 1.0f};
      glMaterialfv(GL_BACK, GL_DIFFUSE, mat_diffuse2);
    }
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
  } else {
    glShadeModel(GL_FLAT);
  }

  if (three_dimensional) {
    resize();
    glEnable(GL_DEPTH_TEST);
  } else {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-radius,radius,-radius,radius,-radius,radius);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glDisable(GL_DEPTH_TEST);
  }

  VirtualGLWindow::ginit();

  setbgColor(0.7,0.7,0.7);
  makeBoundingBox();

  resize();
  updateMat();
  gdk_gl_drawable_gl_end (gldrawable);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGLWindow::~GTKGLWindow() {
  clearClipPlanes();
#if (SPACEDIM==3)
  if (slide_bar) delete slide_bar; slide_bar=0;
  if (toggle_button_array) delete toggle_button_array; 
    toggle_button_array=0;
#endif
  if (bounding_box_list) {
    delete bounding_box_list; bounding_box_list=0;
  }
  clearDisplayLists();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void onClickQuitButtonGTKGLWindow(GtkWidget *w,gpointer d) {
  GtkWidget *dialog=gtk_dialog_new();
  GtkWidget *question=gtk_label_new ("Do you really want to quit?");
  gtk_box_pack_start(GTK_BOX (GTK_DIALOG (dialog)->vbox),
    question, true, true, 0);

  GtkWidget *no_button=gtk_button_new_from_stock(GTK_STOCK_NO);
  g_signal_connect(G_OBJECT(no_button),"clicked",
    G_CALLBACK(onClickedNoGTKGLWindow),d);
  gtk_box_pack_start(GTK_BOX (GTK_DIALOG (dialog)->action_area),
    no_button,true,true,0);
  gtk_widget_show(no_button);

  GtkWidget *yes_button=gtk_button_new_from_stock(GTK_STOCK_YES);
  g_signal_connect(G_OBJECT(yes_button),"clicked",
    G_CALLBACK(onClickedYesGTKGLWindow),d);
  gtk_box_pack_start(GTK_BOX (GTK_DIALOG (dialog)->action_area),
    yes_button,true,true,0);
  gtk_widget_show(yes_button);
  gtk_widget_show(question);
  gtk_widget_show(dialog);

  GTKWindow::gtkMain();

  if (dialog) gtk_widget_destroy(dialog); dialog=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void onClickedNoGTKGLWindow(GtkWidget *w,gpointer d) {
  bool &break_flag=*reinterpret_cast<bool*>(d);
  break_flag=false;
  gtk_main_quit();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void onClickedYesGTKGLWindow(GtkWidget *w,gpointer d) {
  bool &break_flag=*reinterpret_cast<bool*>(d);
  break_flag=true;
  gtk_main_quit();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::eventLoop(bool loop_only_if_escaped) {
  if (loop_only_if_escaped) {
    while (g_main_pending()) {
      g_main_iteration(false);
      if (found_escape) break;
    }
    if (!found_escape) return;
  }

  break_flag=false;
  while (!break_flag) {
    g_main_iteration(false);
  }
  break_flag=false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gboolean onKeyPressEventGTKGLWindow(GtkWidget *w,GdkEventKey *k,
gpointer d) {
  bool &found_escape=*reinterpret_cast<bool*>(d);
  guint key=k->keyval;
  if (key==GDK_Escape || key==GDK_Q || key==GDK_q) found_escape=true;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int GTKGLWindow::numberSurfaces() const {
#if (SPACEDIM==3)
  if (toggle_button_array!=0) {
    CHECK_SAME(number_surfaces,toggle_button_array->numberSurfaces());
  }
  return number_surfaces;
#else
  return 1;
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Vector3 GTKGLWindow::getMouse() {
  finding_mouse=true;
  GTKWindow::gtkMain();
  finding_mouse=false;
  return cursor;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gboolean onButtonPressEventGTKGLWindow(GtkWidget *w,
GdkEventButton *event,gpointer d) {
  GTKGLWindow *glw=reinterpret_cast<GTKGLWindow*>(d);
  Pair<GLdouble> new_screen_coord=
    glw->xyToScreen(event->x,event->y);
  glw->button_pressed=true;
  glw->moved=false; glw->cursor_selected=false;
  switch (event->button) {
    case 1: // Left 
      break;
    case 2: // Middle
      glw->selectClipPlane(new_screen_coord);
      break;
    case 3: // Right
      glw->markFaceLocation(new_screen_coord);
      break;
    default:
      break;
  }
  glw->old_screen_coord=new_screen_coord;
  glw->expose();
  return false;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gboolean onButtonReleaseEventGTKGLWindow(GtkWidget *w,
GdkEventButton *event,gpointer d) {
  bool break_mouse_flag=false;
  GTKGLWindow *glw=reinterpret_cast<GTKGLWindow*>(d);
  Pair<GLdouble> new_screen_coord=
    glw->xyToScreen(event->x,event->y);

  glw->button_pressed=false;
  glw->moved=false; glw->cursor_selected=false;
  switch (event->button) {
    case 1: // Left 
      break;
    case 2: // Middle
#if (SPACEDIM==3)
      if (glw->plot_obj!=0) {
        if (!glw->plot_obj->drawSurface()) {
          VirtualGLWindow::DrawDelimiter dd(glw,
            glw->clip_plane[glw->active_clip_direction]
                           [glw->active_clip_hand]);
          glw->plot_obj->plot(true,
            glw->clip_plane[glw->active_clip_direction]
                           [glw->active_clip_hand]);
        }
      }
#endif
      break;
    case 3: // Right
      if (glw->finding_mouse) break_mouse_flag=true;
      break;
    default:
      break;
  }
  if (glw->mouse_list) delete glw->mouse_list; glw->mouse_list=0;
  glw->old_screen_coord=new_screen_coord;

  glw->best_list_exposed=true;
  glw->expose();
  if (break_mouse_flag) gtk_main_quit();
  return true;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gboolean onMotionNotifyEventGTKGLWindow(GtkWidget *w,
GdkEventMotion *event,gpointer d) {
  GTKGLWindow *glw=reinterpret_cast<GTKGLWindow*>(d);
  Pair<GLdouble> new_screen_coord=
    glw->xyToScreen(event->x,event->y);
  glw->moved=true;
  if (event->state & GDK_BUTTON1_MASK) { //Left button
    glw->rotateImage(new_screen_coord);
  } else if (event->state & GDK_BUTTON2_MASK) { // Middle button
    glw->moveClipPlane(new_screen_coord);
  } else if (event->state & GDK_BUTTON3_MASK) { // Right button
    glw->moveCursor(new_screen_coord);
  }
  glw->expose();
  return true;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::resize() {
  rwidth=drawing_area->allocation.width;
  rheight=drawing_area->allocation.height;

  makeCurrent();

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  GLfloat min_W_H=min(W,H);
  GLfloat max_W_H=min(W,H);
  GLfloat aspect_ratio=max_W_H/min_W_H;
  if (W <= H) {
    glOrtho(-radius,radius,-radius*aspect_ratio,radius*aspect_ratio,
      -radius,radius);
  } else {
    glOrtho(-radius*aspect_ratio,radius*aspect_ratio,-radius,radius,
      -radius,radius);
  }

  glViewport(-5, -5, W,H);
  updateMat();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::rescale(const Vector3 &l,const Vector3 &h) {
  VirtualGLWindow::rescale(l,h);
  radius=norm(bound);
  makeBoundingBox();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::makeBoundingBox() {
  if  (!bounding_box_list) {
    bounding_box_list=OPERATOR_NEW MyDisplayList();
  }
  { ListDelimiter ld(bounding_box_list->getBase(),GL_COMPILE);
    if (use_lighting) {
      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHT0);
    }
    setfgColor("black");
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    {
      Vector3 low=-bound;
      GLenum mode=GL_RENDER;
      GLenum type=GL_QUADS;
      drawBox(&low,&bound,&mode,&type);
    }
    int three=3;
    setLineWidth(&three);
    for (int i=0;i<3;i++) {
      int i1=(i+1)%3,i2=(i+2)%3;
      if (i==0) setfgColor("red");
      if (i==1) setfgColor("green");
      if (i==2) setfgColor("blue");
      Vector3 low(-bound),high(-bound);
      low[i]=clip_plane[i][LEFT_HAND]->getLocation();
      high[i]=clip_plane[i][RIGHT_HAND]->getLocation();

      drawLine(&low,&high);
      low[i1]=high[i1]=bound[i1];
      drawLine(&low,&high);
      low[i2]=high[i2]=bound[i2];
      drawLine(&low,&high);
      low[i1]=high[i1]=-bound[i1];
      drawLine(&low,&high);
    }
    int one=1;
    setLineWidth(&one);

    float fudge=radius/30.,len=radius/10.,len1=radius/10.+radius/50.;
    Vector3 fudged_origin=-bound;
    for (int i=0;i<3;i++) fudged_origin[i]-=fudge;

    Vector3 vector_end=fudged_origin; vector_end[0]+=len;
    setfgColor("red"); // setfgColor(1.,0.,0.);
    drawLine(&fudged_origin,&vector_end);
    glRasterPos3f(-bound[0]-fudge+len1,-bound[1]-fudge,-bound[2]-fudge);
    printString("x");

    vector_end=fudged_origin; vector_end[1]+=len;
    setfgColor("green"); // setfgColor(0.,1.,0.);
    drawLine(&fudged_origin,&vector_end);
    glRasterPos3f(-bound[0]-fudge,-bound[1]-fudge+len1,-bound[2]-fudge);
    printString("y");

    vector_end=fudged_origin; vector_end[2]+=len;
    setfgColor("blue"); // setfgColor(0.,0.,1.);
    drawLine(&fudged_origin,&vector_end);
    glRasterPos3f(-bound[0]-fudge,-bound[1]-fudge,-bound[2]-fudge+len1);
    printString("z");
    if (use_lighting) {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
    }
  }
}
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGLWindow::SlideBar* GTKGLWindow::slideBar() {
  CHECK_SAME(number_surfaces,1)
  if (!slide_bar) slide_bar=OPERATOR_NEW SlideBar(drawing_area,H); 
  return slide_bar;
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::deleteSlideBar() {
  if (slide_bar) delete slide_bar; slide_bar=0;
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGLWindow::ToggleButtonArray* GTKGLWindow::toggleButtonArray(int ns) {
  CHECK_TEST(ns>1)
  if (toggle_button_array==0) {
    CHECK_SAME(number_surfaces,1)
    number_surfaces=ns;
    toggle_button_array=
      OPERATOR_NEW ToggleButtonArray(vbox,number_surfaces); 
  }
  return toggle_button_array;
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::deleteToggleButtonArray() {
  if (toggle_button_array) delete toggle_button_array; 
  toggle_button_array=0;
  number_surfaces=1;
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::printString(const char *s) {
  glListBase(fontbase);
  glCallLists(strlen(s),GL_UNSIGNED_BYTE,
    reinterpret_cast<const GLubyte *>(s));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::expose() {
  makeCurrent();

  ASSERT(gdk_gl_drawable_gl_begin(gldrawable, glcontext));
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
  disableClipPlanes();

  if (new_model) recalcModelView();
  if (button_pressed || show_bounding_box) bounding_box_list->call();
  if (mouse_list) mouse_list->call();

  enableClipPlanes();
  if ((best_list_exposed || rotate_best) && best_list!=0) { 
    best_list->call();
  } else if (user_list) {
    user_list->call();
  }
  plotClipPlanes();

  if (gdk_gl_drawable_is_double_buffered(gldrawable)) {
    gdk_gl_drawable_swap_buffers(gldrawable);
  } else glFlush();
  gdk_gl_drawable_gl_end(gldrawable);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gboolean onExposeEventGTKGLWindow(GtkWidget *o,GdkEventExpose *event,
gpointer d) {
  GTKGLWindow *gtkw=reinterpret_cast<GTKGLWindow*>(d);
  gtkw->expose();
  return true;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//left button motion
void GTKGLWindow::rotateImage(const Pair<GLdouble> &new_screen_coord) {
  Quaternion<double> rotate_quat=
    trackBall(screenToUnitCircle(old_screen_coord),
              screenToUnitCircle(new_screen_coord));
  curquat=curquat*rotate_quat;
  old_screen_coord=new_screen_coord;
  new_model = true;
  best_list_exposed=rotate_best;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//right button motion
void GTKGLWindow::moveCursor(const Pair<GLdouble> &new_screen_coord) {
  if ((new_screen_coord-old_screen_coord).norm()<4.) return;

  GLdouble max_inner_prod=ZERO;
  unsigned int axis=0;
  for (unsigned int i=0;i<3;i++) {
    GLdouble d=innerProduct(wuvec[i],new_screen_coord-old_screen_coord);
    GLdouble abs_d=d>=ZERO?d:-d;
    if (abs_d>max_inner_prod) {
      max_inner_prod=abs_d; axis=i; 
    }
  }

  Vector3 axis_vector; axis_vector[axis]=ONE;
  Line new_cursor_line=objectCoords(new_screen_coord);
  Line allowed_cursor_line(cursor,axis_vector);
  allowed_cursor_line.nearestPositionTo(new_cursor_line);
  cursor=allowed_cursor_line.getPosition();
  
  if (cursor[axis]<-bound[axis]) cursor[axis]=-bound[axis];
  else if (cursor[axis]>bound[axis]) cursor[axis]=bound[axis];

  crossHairs(cursor);

  old_screen_coord=screenCoords(cursor);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//middle button motion
void GTKGLWindow::moveClipPlane(const Pair<GLdouble> &new_screen_coord){
  Line new_cursor_line=objectCoords(new_screen_coord);
  closest_edge.nearestPositionTo(new_cursor_line);
  cursor=closest_edge.getPosition();
  AxisClipPlane *cp=clip_plane[active_clip_direction][active_clip_hand];

  if (cp->getHand()==LEFT_HAND) {
    cursor[active_clip_direction]=max(-bound[active_clip_direction],
      min(clip_plane[active_clip_direction][RIGHT_HAND]->getLocation(),
          cursor[active_clip_direction]));
  } else {
    cursor[active_clip_direction]=min(bound[active_clip_direction],
      max(clip_plane[active_clip_direction][LEFT_HAND]->getLocation(),
          cursor[active_clip_direction]));
  }

  old_screen_coord=screenCoords(cursor);
  cp->setLocation(cursor[active_clip_direction]);

  makeBoundingBox();
  crossHairs(cursor);
#if (SPACEDIM==3)
  if (plot_obj!=0) {
    if (!plot_obj->drawSurface()) {
      VirtualGLWindow::DrawDelimiter dd(this,cp);
      plot_obj->plot(false,cp);
    }
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::drawLine(const Vector3 *start,const Vector3 *end) const {
  { PrimitiveDelimiter pd(GL_LINES);
    glVertex3f((*start)[0],(*start)[1],(*start)[2]);
    glVertex3f((*end)[0],(*end)[1],(*end)[2]);
  }
}
#if (SPACEDIM==2)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::colorTriangle(const GdkColor *xc,const Vector3 *v0,
const Vector3 *v1,const Vector3 *v2) const {
  setfgColor(*xc);
  setPolygonMode(GL_FRONT,GL_FILL,GL_CCW);
  Vector3 normal=cross((*v1)-(*v0),(*v2)-(*v0));
  double nrm=norm(normal);
  if (nrm>ZERO) normal/=nrm;
  { TriangleDelimiter td;
    if (normal[2]>ZERO) {
      glNormal3d(normal[0],normal[1],normal[2]);
      glVertex3d((*v0)[0],(*v0)[1],(*v0)[2]);
      glVertex3d((*v1)[0],(*v1)[1],(*v1)[2]);
      glVertex3d((*v2)[0],(*v2)[1],(*v2)[2]);
    } else {
      normal=-normal;
      glNormal3d(normal[0],normal[1],normal[2]);
      glVertex3d((*v0)[0],(*v0)[1],(*v0)[2]);
      glVertex3d((*v2)[0],(*v2)[1],(*v2)[2]);
      glVertex3d((*v1)[0],(*v1)[1],(*v1)[2]);
    }
  }
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::colorTriangle(const GdkColor *xc,const Vector3 *v0,
const Vector3 *v1,const Vector3 *v2) const {
  CHECK_TEST(isDrawing())
  setfgColor(*xc);
  setPolygonMode(GL_FRONT,GL_FILL,GL_CCW);
  Vector3 normal=cross((*v1)-(*v0),(*v2)-(*v0));
  double nrm=norm(normal);
  if (nrm>ZERO) normal/=nrm;
  CHECK_TEST((*v0)>=-bound && (*v0)<=bound)
  CHECK_TEST((*v1)>=-bound && (*v1)<=bound)
  CHECK_TEST((*v2)>=-bound && (*v2)<=bound)

  GLfloat emit_max=GLfloat(maxIntensity)*5.f;
  const GLfloat mat_emission1[]=
    {GLfloat(xc->red)/emit_max,GLfloat(xc->green)/emit_max,
     GLfloat(xc->blue)/emit_max,1.f};
  glMaterialfv(GL_FRONT, GL_EMISSION, mat_emission1);
  const GLfloat mat_emission2[]=
   {1.f-mat_emission1[0],1.f-mat_emission1[1],1.f-mat_emission1[2],1.f};
  glMaterialfv(GL_BACK, GL_EMISSION, mat_emission2);

  { TriangleDelimiter td;
    if (surface_normal_direction==INCREASING_FUNCTION) {
      glNormal3d(normal[0],normal[1],normal[2]);
      glVertex3d((*v0)[0],(*v0)[1],(*v0)[2]);
      glVertex3d((*v1)[0],(*v1)[1],(*v1)[2]);
      glVertex3d((*v2)[0],(*v2)[1],(*v2)[2]);
    } else {
      glNormal3d(-normal[0],-normal[1],-normal[2]);
      glVertex3d((*v2)[0],(*v2)[1],(*v2)[2]);
      glVertex3d((*v1)[0],(*v1)[1],(*v1)[2]);
      glVertex3d((*v0)[0],(*v0)[1],(*v0)[2]);
    }
  }
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::drawVector(const Vector3 *base,const Vector3 *head) const
{
  drawLine(base,head);

  double c=cos(2.6),s=sin(2.6);
  Vector3 arrow=*head,dif=(*head)-(*base);
  Vector3 adif(abs(dif[0]),abs(dif[1]),abs(dif[2]));
  int d=0; adif.max(d);
  int d1=(d+1)%3,d2=(d+2)%3;
  arrow[d1]+=(c*dif[d1]-s*dif[d2])*ARROW_SIZE;
  arrow[d2]+=(s*dif[d1]+c*dif[d2])*ARROW_SIZE;
  drawLine(head,&arrow);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::drawVector(const Vector3 *base,const Vector3 *head,
const AxisClipPlane *cp) const {
  unsigned int clip_dir=cp->getDirection();
  GLdouble clip_location=cp->getLocation();
  Vector3 b=*base; b[clip_dir]=clip_location;
  Vector3 h=*head; h[clip_dir]=clip_location;
  drawLine(&b,&h);

  int d1=(clip_dir+1)%3,d2=(clip_dir+2)%3;
  double c=cos(2.6),s=sin(2.6);
  Vector3 a=h,dif=h-b;
  a[d1]+=(c*dif[d1]-s*dif[d2])*ARROW_SIZE;
  a[d2]+=(s*dif[d1]+c*dif[d2])*ARROW_SIZE;
  drawLine(&h,&a);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::plotObject(bool do_best) {
  if (rotate_best) do_best=true;
  if (plot_obj!=0) {
    if (plot_obj->drawSurface()) {
      freeDisplayList(do_best);
      newPage();
      { VirtualGLWindow::DrawDelimiter dd(do_best,this);
        plot_obj->plot(do_best);
      }
      best_list_exposed=do_best;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::colorRect(const GdkColor *xc,const Vector3 *low,
const Vector3 *high,const int *dir) const {
  CHECK_TEST(inside_draw_delimiter)
  setfgColor(*xc);
  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  { QuadDelimiter qd;
    Vector3 v=*low; v[*dir]=HALF*((*low)[*dir]+(*high)[*dir]);
    int dir1=(*dir+1)%3,dir2=(*dir+2)%3;
    glVertex3d(v[0],v[1],v[2]);
    v[dir2]=(*high)[dir2];
    glVertex3d(v[0],v[1],v[2]);
    v[dir1]=(*high)[dir1];
    glVertex3d(v[0],v[1],v[2]);
    v[dir2]=(*low)[dir2];
    glVertex3d(v[0],v[1],v[2]);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::makeCurrent() {
  gdk_gl_drawable_make_current(gldrawable,glcontext);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::recalcModelView() {
  GLfloat m[4][4];

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  curquat.rotationMatrix(m);
  glMultMatrixf(&m[0][0]);
  glTranslatef(ZERO,ZERO,ZERO);
  new_model = false; 
  mat_valid = false;
  updateMat();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// project an object coord into X win coord
Pair<GLdouble> GTKGLWindow::screenCoords(const Vector3 &object_coord) {
  GLdouble winx, winy, winz;
  updateMat();
  gluProject(object_coord[0],object_coord[1],object_coord[2],
              ModelMat, ProjMat, viewport, &winx, &winy, &winz );
  return Pair<GLdouble>(winx,winy);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// convert a X win coord into a line in obj coord
Line GTKGLWindow::objectCoords(Pair<GLdouble> screen_coord) {
  updateMat();

  Vector3 pos0,pos1;
  GLint success=gluUnProject(screen_coord[0],screen_coord[1],ZERO,
    ModelMat,ProjMat,viewport,&pos0[0],&pos0[1],&pos0[2]);
  ASSERT(success);

  success=gluUnProject(screen_coord[0],screen_coord[1],ONE,ModelMat,
    ProjMat,viewport,&pos1[0],&pos1[1],&pos1[2]);
  ASSERT(success);

  return Line(pos0,pos1-pos0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::markFaceLocation(const Pair<double> &new_cursor_coord)
{
  cursor=faceLocation(new_cursor_coord);
  crossHairs(cursor);

  Pair<GLdouble> screen_coord_low=screenCoords(-bound);
  for (int i=0;i<3;i++) { // around the low vertex -bound
    Vector3 neighbor(-bound); neighbor[i]=bound[i];

    Pair<GLdouble> screen_coord=screenCoords(neighbor);
    wuvec[i]=screen_coord-screen_coord_low;
    GLdouble d=wuvec[i].norm();
    if (d>ZERO) wuvec[i] /= d;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::selectClipPlane(
const Pair<GLdouble> &new_screen_coord) {
  Line new_cursor_line=objectCoords(new_screen_coord);
  GLdouble min_dist=DBL_MAX;
  for (unsigned int i=0;i<3;i++) {
    Vector3 axis; axis[i]=ONE;
    Vector3 vertex(-bound);
    int i1=(i+1)%3,i2=(i+2)%3,c1=0;
    for (vertex[i1]=-bound[i1];c1<2;vertex[i1]=-vertex[i1],c1++) {
      int c2=0;
      for (vertex[i2]=-bound[i2];c2<2;vertex[i2]=-vertex[i2],c2++) {
        Line current_edge(vertex,axis);
        double dist=current_edge.nearestPositionTo(new_cursor_line);
        if (dist<min_dist) {
          closest_edge=current_edge;
          active_clip_direction=i;
          min_dist=dist;
        }
      }
    }
  }
  cursor=closest_edge.getPosition();
  if (cursor[active_clip_direction]<-bound[active_clip_direction]) {
    cursor[active_clip_direction]=-bound[active_clip_direction];
  } else if(cursor[active_clip_direction]>bound[active_clip_direction]){
    cursor[active_clip_direction]=bound[active_clip_direction];
  }

  GLdouble middle_pos=
    HALF*(clip_plane[active_clip_direction][LEFT_HAND]->getLocation()
         +clip_plane[active_clip_direction][RIGHT_HAND]->getLocation());
  if (cursor[active_clip_direction]<=middle_pos) {
    active_clip_hand=LEFT_HAND;
  } else active_clip_hand=RIGHT_HAND;
  old_screen_coord=screenCoords(cursor);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::initialize(const VirtualGLWindow::PlotObj *po) {
  plot_obj=po;
#if (SPACEDIM==3)
  if (plot_obj) {
    if (!plot_obj->drawSurface()) {
      disableClipPlanes();
      bool do_best=true;
      for (unsigned int i=0;i<3;i++) {
        { VirtualGLWindow::DrawDelimiter dd(this,
            clip_plane[i][LEFT_HAND]);
          plot_obj->plot(do_best,clip_plane[i][LEFT_HAND]);
        }
        { VirtualGLWindow::DrawDelimiter dd(this,
            clip_plane[i][RIGHT_HAND]);
          plot_obj->plot(do_best,clip_plane[i][RIGHT_HAND]);
        }
      }
    }
  }
#endif
  enableClipPlanes();
//the following is necessary for rotation of the initial image
  if (!rotate_best) plotObject(false);
//the following provides the initial image
  plotObject(true);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::plotClipPlanes() const {
#if (SPACEDIM==3)
  if (plot_obj) {
    if (!plot_obj->drawSurface()) {
      enableClipPlanes();
      for (unsigned int i=0;i<3;i++) {
        clip_plane[i][LEFT_HAND]->callPlotList();
        clip_plane[i][RIGHT_HAND]->callPlotList();
      }
    }
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::terminatePlotObj() {
  if (plot_obj) delete plot_obj;
  plot_obj=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//this must write what openRaster will read
void GTKGLWindow::createRaster(FILE *os) {
  CHECK_TEST(!image);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::writeRaster(FILE *os) {
  CHECK_TEST(!image);
  if (double_buffer) gdk_gl_drawable_swap_buffers(gldrawable);
  ASSERT(gdk_pixbuf_get_from_drawable(image,drawing_area->window,
    gtk_widget_get_colormap(drawing_area),0,0,0,0,W,H));
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
void GTKGLWindow::closeRaster() {
  CHECK_TEST(!image);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::openRaster(FILE *is) {
  CHECK_TEST(!image);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::readRaster(FILE *is) {
  CHECK_TEST(!image);
  char label[LABEL_LENGTH];
  for (int j=0;j<LABEL_LENGTH;j++) label[j]=0;
  snprintf(label,LABEL_LENGTH,"%s.raster.%d",name,raster_number);
  
  GError **error=0;
  image=gdk_pixbuf_new_from_file(label,error);
  ASSERT(image);
  gdk_draw_pixbuf(drawing_area->window,gc,image,0,0,0,0,
    W,H,GDK_RGB_DITHER_NONE,0,0);
  g_object_unref(image); image=0;
  raster_number++;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::printOn(ostream &os) const {
  os << "GTKGLWindow: "
//   << "\n\trwidth,rheight = " << rwidth << " " << rheight
//   << "\n\tcursor = " << cursor
//   << "\n\tbound = " << bound
//   << "\n\tradius = " << radius
//   << "\n\tactive_clip_direction = " << active_clip_direction
//   << "\n\tactive_clip_hand = " << active_clip_hand
     << endl;
//os << "\n\tperspective = " << perspective
//   << "\n\tthree_dimensional = " << three_dimensional
//   << "\n\tbutton_pressed = " << button_pressed
//   << "\n\tnew_model = " << new_model
//   << "\n\tmat_valid = " << mat_valid
//   << "\n\tmoved = " << moved
//   << "\n\tcursor_selected = " << cursor_selected
//   << "\n\tbest_list_exposed = " << best_list_exposed
//   << endl;
//os << "\n\tplot_obj = " << (void*) plot_obj
//   << "\n\tbounding_box_list = " << bounding_box_list
//   << "\n\tmouse_list = " << mouse_list
//   << "\n\tuser_list = " << user_list
//   << "\n\tbest_list = " << best_list
//   << endl;
  VirtualGLWindow::printOn(os);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::drawBox(const Vector3 *low,const Vector3 *high,
const GLenum *mode,const GLenum *type) const {
  CHECK_TEST(*low<*high)
  switch (*type) {
    case GL_POINTS: // normal along radius from center
      { PrimitiveDelimiter pd(*type);
        glNormal3f(-ONE,-ONE,-ONE);
        glVertex3d((*low) [0],(*low) [1],(*low) [2]);
        glNormal3f( ONE,-ONE,-ONE);
        glVertex3d((*high)[0],(*low) [1],(*low) [2]);

        glNormal3f(-ONE, ONE,-ONE);
        glVertex3d((*low) [0],(*high)[1],(*low) [2]);
        glNormal3f( ONE, ONE,-ONE);
        glVertex3d((*high)[0],(*high)[1],(*low) [2]);

        glNormal3f(-ONE,-ONE, ONE);
        glVertex3d((*low) [0],(*low) [1],(*high)[2]);
        glNormal3f( ONE,-ONE, ONE);
        glVertex3d((*high)[0],(*low) [1],(*high)[2]);

        glNormal3f(-ONE, ONE, ONE);
        glVertex3d((*low) [0],(*high)[1],(*high)[2]);
        glNormal3f( ONE, ONE, ONE);
        glVertex3d((*high)[0],(*high)[1],(*high)[2]);
      }
      break;
    case GL_LINES:
      { PrimitiveDelimiter pd(*type);
        { glNormal3f(ZERO,-ONE,-ONE);
          glVertex3d((*low) [0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);

          glNormal3f(ZERO, ONE,-ONE);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);

          glNormal3f(ZERO,-ONE, ONE);
          glVertex3d((*low) [0],(*low) [1],(*high)[2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);

          glNormal3f(ZERO, ONE, ONE);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);
          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
        }
        { glNormal3f(-ONE,ZERO,-ONE);
          glVertex3d((*low) [0],(*low) [1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);

          glNormal3f( ONE,ZERO,-ONE);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);

          glNormal3f(-ONE,ZERO, ONE);
          glVertex3d((*low) [0],(*low) [1],(*high)[2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);

          glNormal3f( ONE,ZERO, ONE);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);
          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
        }
        { glNormal3f(-ONE,-ONE,ZERO);
          glVertex3d((*low) [0],(*low) [1],(*low) [2]);
          glVertex3d((*low) [0],(*low) [1],(*high)[2]);

          glNormal3f( ONE,-ONE,ZERO);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);

          glNormal3f(-ONE, ONE,ZERO);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);

          glNormal3f( ONE, ONE,ZERO);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
        }
      }
      break;
    case GL_POLYGON:
    case GL_TRIANGLE_STRIP:
    case GL_TRIANGLE_FAN:
    case GL_QUADS:
      { QuadDelimiter qd;
        { glNormal3f(-ONE,ZERO,ZERO);
          if (*mode==GL_SELECT) glLoadName(0); // left
          glVertex3d((*low) [0],(*low) [1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);
          glVertex3d((*low) [0],(*low) [1],(*high)[2]);
          if (*mode==GL_FEEDBACK) glPassThrough(1.0); // left->right
          glNormal3f( ONE,ZERO,ZERO);
          if (*mode==GL_SELECT) glLoadName(1); // right
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);
        }
        if (*mode==GL_FEEDBACK) glPassThrough(2.0); // right->bottom
        { glNormal3f(ZERO,-ONE,ZERO);
          if (*mode==GL_SELECT) glLoadName(2); // bottom
          glVertex3d((*low)[0] ,(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);
          glVertex3d((*low)[0] ,(*low) [1],(*high)[2]);
          if (*mode==GL_FEEDBACK) glPassThrough(3.0); // bottom->top
          glNormal3f(ZERO, ONE,ZERO);
          if (*mode==GL_SELECT) glLoadName(3); // top
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);
        }
        if (*mode==GL_FEEDBACK) glPassThrough(4.0); // top->back
        { glNormal3f(ZERO,ZERO,-ONE);
          if (*mode==GL_SELECT) glLoadName(4); // back
          glVertex3d((*low) [0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);
          if (*mode==GL_FEEDBACK) glPassThrough(5.0); // back->front
          glNormal3f(ZERO,ZERO, ONE);
          if (*mode==GL_SELECT) glLoadName(5); // front
          glVertex3d((*low) [0],(*low) [1],(*high)[2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);
          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);
        }
      }
      break;
    case GL_TRIANGLES:
      { PrimitiveDelimiter pd(*type);
        { glNormal3f(-ONE,ZERO,ZERO);
          glVertex3d((*low) [0],(*low) [1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);
          glVertex3d((*low) [0],(*low) [1],(*high)[2]);

          glVertex3d((*low) [0],(*high)[1],(*high)[2]);
          glVertex3d((*low) [0],(*low) [1],(*high)[2]);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);

          glNormal3f( ONE,ZERO,ZERO);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);

          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
        }
        { glNormal3f(ZERO,-ONE,ZERO);
          glVertex3d((*low)[0] ,(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);
          glVertex3d((*low)[0] ,(*low) [1],(*high)[2]);

          glVertex3d((*high)[0],(*low) [1],(*high)[2]);
          glVertex3d((*low)[0] ,(*low) [1],(*high)[2]);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);

          glNormal3f(ZERO, ONE,ZERO);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);

          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);
          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
        }
        { glNormal3f(ZERO,ZERO,-ONE);
          glVertex3d((*low) [0],(*low) [1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);

          glVertex3d((*high)[0],(*high)[1],(*low) [2]);
          glVertex3d((*low) [0],(*high)[1],(*low) [2]);
          glVertex3d((*high)[0],(*low) [1],(*low) [2]);

          glNormal3f(ZERO,ZERO, ONE);
          glVertex3d((*low) [0],(*low) [1],(*high)[2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);

          glVertex3d((*high)[0],(*high)[1],(*high)[2]);
          glVertex3d((*low) [0],(*high)[1],(*high)[2]);
          glVertex3d((*high)[0],(*low) [1],(*high)[2]);
        }
      }
      break;
    case GL_LINE_STRIP:
    case GL_LINE_LOOP: // normal is bogus
      { PrimitiveDelimiter pd(*type);
        glNormal3f(ZERO,ZERO,-ONE);
        glVertex3d((*low) [0],(*low) [1],(*low) [2]);
        glVertex3d((*low) [0],(*high)[1],(*low) [2]);
        glVertex3d((*low) [0],(*high)[1],(*high)[2]);
        glVertex3d((*low) [0],(*low) [1],(*high)[2]);
        glVertex3d((*low) [0],(*low) [1],(*low) [2]);
        glVertex3d((*high)[0],(*low) [1],(*low) [2]);
        glVertex3d((*high)[0],(*low) [1],(*high)[2]);
        glVertex3d((*low) [0],(*low) [1],(*high)[2]);
      }
      { PrimitiveDelimiter pd(*type);
        glNormal3f(ZERO,ZERO, ONE);
        glVertex3d((*high)[0],(*high)[1],(*high)[2]);
        glVertex3d((*high)[0],(*low) [1],(*high)[2]);
        glVertex3d((*high)[0],(*low) [1],(*low) [2]);
        glVertex3d((*high)[0],(*high)[1],(*low) [2]);
        glVertex3d((*high)[0],(*high)[1],(*high)[2]);
        glVertex3d((*low) [0],(*high)[1],(*high)[2]);
        glVertex3d((*low) [0],(*high)[1],(*low) [2]);
        glVertex3d((*high)[0],(*high)[1],(*low) [2]);
      }
      break;
    case GL_QUAD_STRIP: // normal is bogus
      { PrimitiveDelimiter pd(*type);
        glNormal3f(ZERO,ZERO,-ONE);
        glVertex3d((*low) [0],(*low) [1],(*low) [2]);
        glVertex3d((*low) [0],(*low) [1],(*high)[2]);
        glVertex3d((*low) [0],(*high)[1],(*low) [2]);
        glVertex3d((*low) [0],(*high)[1],(*high)[2]);
        glVertex3d((*high)[0],(*high)[1],(*low) [2]);
        glVertex3d((*high)[0],(*high)[1],(*high)[2]);
        glVertex3d((*high)[0],(*low) [1],(*low) [2]);
        glVertex3d((*high)[0],(*low) [1],(*high)[2]);
      }
      break;
    default:
      break;
  }
}



#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGLWindow::SlideBar::SlideBar(GtkWidget *glarea,gint H) : 
var_min(HUGE_VAL),var_max(-HUGE_VAL),frac(HALF),
slide_bar(0) {
  GtkWidget *vbox=gtk_widget_get_parent(glarea);
  slide_bar=gtk_hscale_new_with_range(0.,1.,1.);
  gtk_widget_set_size_request(slide_bar,vbox->allocation.width,
    vbox->allocation.height/16);
  gtk_box_pack_end(GTK_BOX(vbox),slide_bar,false,false,0);
  gtk_widget_show(slide_bar);
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::SlideBar::rescale(double low,double high) {
  CHECK_TEST(low<high)
  short decimal_points=
    static_cast<short>(floor(log10(high-low)));
  double sig_digit=pow(10.,decimal_points);
  int slide_min=0,slide_max=1;
  if (decimal_points<0) {
    slide_max=static_cast<int>(ceil(high/sig_digit));
    slide_min=static_cast<int>(floor(low/sig_digit));
    decimal_points=-decimal_points;
  } else {
    slide_max=static_cast<int>(ceil(high));
    slide_min=static_cast<int>(floor(low));
    decimal_points=0;
  }
  var_min=static_cast<double>(slide_min)*sig_digit;
  var_max=static_cast<double>(slide_max)*sig_digit;
  gtk_range_set_range(GTK_RANGE(slide_bar),slide_min,slide_max);
  gtk_range_set_increments(GTK_RANGE(slide_bar),sig_digit*0.0625,
    sig_digit*0.0625);
  double v=getValue();
  gtk_range_set_value(GTK_RANGE(slide_bar),v);
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::SlideBar::setValue(double v) {
  v=max(var_min,min(var_max,v));
  frac=(v-var_min)/(var_max-var_min);
  gtk_range_set_value(GTK_RANGE(slide_bar),v);
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::SlideBar::setCallback(
void (*callback)(GtkWidget*,gpointer),gpointer user_data) {
  g_signal_connect(G_OBJECT(slide_bar),"value-changed",
    G_CALLBACK(callback),user_data);
}
#endif
#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::SlideBar::printOn(ostream &os) const {
  os << "GTKGLWindow::SlideBar: slide_bar = " << slide_bar << endl;
  os << "\tvar_min = " << var_min << endl;
  os << "\tvar_max = " << var_max << endl;
  os << "\tfrac = " << frac << endl;
}
#endif

#if (SPACEDIM==3)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGLWindow::ToggleButtonArray::ToggleButtonArray(GtkWidget *vbox,
int n) : hbutton_box(0),toggle_button_array(n) {
  hbutton_box=gtk_hbox_new(true,0);
  gtk_box_pack_end(GTK_BOX(vbox),hbutton_box,false,false,0);

  int width=vbox->allocation.width/n;
  int height=vbox->allocation.height/16;
  for (int i=0;i<n;i++) {
    char label[LABEL_LENGTH];
    snprintf(label,LABEL_LENGTH,"%d",i);
    toggle_button_array[i]=gtk_check_button_new_with_label(label);
    gtk_toggle_button_set_active(
      GTK_TOGGLE_BUTTON(toggle_button_array[i]),false);
    gtk_widget_set_size_request(toggle_button_array[i],width,height);
    gtk_box_pack_start(GTK_BOX(hbutton_box),toggle_button_array[i],true,
      false,0);
    gtk_widget_show(toggle_button_array[i]);
  }
  gtk_toggle_button_set_active(
    GTK_TOGGLE_BUTTON(toggle_button_array[0]),true);
  gtk_widget_show(hbutton_box);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::ToggleButtonArray::setCallback(
void (*callback)(GtkWidget*,gpointer),gpointer user_data) {
  int nt=toggle_button_array.getNumber();
  for (int i=0;i<nt;i++) {
    g_signal_connect(G_OBJECT(toggle_button_array[i]),"toggled",
      G_CALLBACK(callback),user_data);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGLWindow::ToggleButtonArray::printOn(ostream &os) const {
  os << "GTKGLWindow::ToggleButtonArray: "
     << "\n\tnumberSurfaces = " << numberSurfaces()
     << endl;
  for (int i=0;i<numberSurfaces();i++) {
    os << "\t drawSurface(" << i << ") = " << drawSurface(i) << endl;
  }
}
#endif

#if (SPACEDIM==3)
#include "NumPtr.C"
INSTANTIATE_NUMPTR(GtkWidget*)
#endif
#endif
#endif
