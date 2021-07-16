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
// "$Header: /home/faculty/johnt/cvs/gui/GTKColorEditor.C,v 1.1 2005/06/08 11:19:01 johnt Exp $"
#ifdef USE_GTK
#include <gtk/gtkcolorseldialog.h>
#include <gtk/gtkdialog.h>
#include <gtk/gtkdrawingarea.h>
#include <gtk/gtkframe.h>
#include <gtk/gtkhscale.h>
#include <gtk/gtkradiobutton.h>
#include <gtk/gtksignal.h>
#include <gtk/gtkstock.h>
#include <gtk/gtktable.h>
#include <gtk/gtktogglebutton.h>
#include <gtk/gtkvbox.h>
#include "GTKColorEditor.H"
#include "ISLList.C"
#include "NumPtr.C"
GtkWidget* GTKColorEditor::dialog=0;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKColorEditor::GTKColorEditor(GTKColormap *map,GtkWidget*) : 
gtkcolormap(map),toggle(map->getNumberColors()),gc(0),form(0),swatch(0),
selected(-1),gui(0) {
  GdkWindow *root_window=gdk_get_default_root_window();
  GdkVisual *visual=gdk_drawable_get_visual(root_window);
  GdkVisualType visual_type=visual->type;
  num_cmap_colors=map->getNumberColors();
  toggle.initialize(0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::selectColorCallback(GtkWidget *w,
gpointer user_data) {
  GTKColorEditor *ce=reinterpret_cast<GTKColorEditor*>(user_data);
  ce->selected==-1;
  for (int i=0;i<ce->toggle.getNumber();i++) {
    if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(ce->toggle[i]))) {
      ce->selected=i;
      break;
    }
  }
  if (ce->warnUserNoColor(ce->swatch)) return;
  const GdkColor &color=ce->gtkcolormap->getColor(ce->selected);
//the following generates a call to changeColorCallback
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(
    GTK_COLOR_SELECTION_DIALOG(ce->swatch)->colorsel),&color);
  gtk_color_selection_set_previous_color(GTK_COLOR_SELECTION(
    GTK_COLOR_SELECTION_DIALOG(ce->swatch)->colorsel),&color);

  GtkRequisition requisition;
  gtk_widget_size_request(ce->drawing_area,&requisition);

  ce->drawRGBBars(ce->selected,ce->selected);
  ce->drawRGBCurves(ce->selected,ce->selected);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::changeColorCallback(GtkWidget *w,
gpointer user_data) {
  GTKColorEditor *ce=reinterpret_cast<GTKColorEditor*>(user_data);
  if (ce->warnUserNoColor(w)) return;
  GdkColor color;
  gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(w),&color);
  ce->colorToggle(ce->toggle[ce->selected],color);
  exposeCallback(GTK_OBJECT(ce->drawing_area),
    reinterpret_cast<GdkEventExpose*>(0),user_data);

  GtkRequisition requisition;
  gtk_widget_size_request(ce->drawing_area,&requisition);

  gint bar_width=static_cast<gint>(
    static_cast<double>(requisition.width-ce->num_cmap_colors-1)
   /static_cast<double>(ce->num_cmap_colors));
  gint bar_height=requisition.height;
  gdk_gc_set_rgb_fg_color(ce->gc,&color);
  gint x=1+ce->selected*(bar_width+1);
  gdk_draw_rectangle(ce->drawing_area->window,ce->gc,TRUE,x,0,
    bar_width,bar_height/2);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::printOn(ostream &os) const {
  os << "GTKColorEditor: "
     << "num_cmap_colors = " << num_cmap_colors
     << endl;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool GTKColorEditor::warnUserNoColor(GtkWidget *parent) const {
  if (selected!=-1) return FALSE;
  if (dialog==0) {
    GtkWidget *dialog=gtk_dialog_new_with_buttons(
      "noColorWarningDialog",0,GTK_DIALOG_DESTROY_WITH_PARENT,
      GTK_STOCK_OK,GTK_RESPONSE_ACCEPT,0);
  }
  gtk_dialog_run(GTK_DIALOG(dialog));
  return FALSE;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::realizeCallback(GtkObject *o,gpointer user_data) {
  GtkWidget *w=reinterpret_cast<GtkWidget*>(o);
  GTKColorEditor* ce=reinterpret_cast<GTKColorEditor*>(user_data);
  ce->gc=gdk_gc_new(w->window);
  ASSERT(ce->gc);
  gdk_window_clear(w->window);
  exposeCallback(GTK_OBJECT(ce->drawing_area),
    reinterpret_cast<GdkEventExpose*>(0),user_data);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::colorToggle(GtkWidget *tog,const GdkColor &bg) 
const {
  gtk_widget_modify_bg(tog,GTK_STATE_NORMAL,&bg);
  gtk_widget_modify_bg(tog,GTK_STATE_ACTIVE,&bg);
  gtk_widget_modify_bg(tog,GTK_STATE_PRELIGHT,&bg);
  gtk_widget_modify_bg(tog,GTK_STATE_SELECTED,&bg);
  gtk_widget_modify_bg(tog,GTK_STATE_INSENSITIVE,&bg);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::configureToggle(GtkWidget *table,GtkWidget *tog,
guint w,guint h,guint rows,int i) const {
  gtk_widget_set_size_request(tog,w,h);
  colorToggle(tog,gtkcolormap->getColor(i));
  guint left_attach=i/rows;
  guint top_attach=i-left_attach*rows;
  gtk_table_attach_defaults(GTK_TABLE(table),tog,
    left_attach,left_attach+1,top_attach,top_attach+1);
  g_signal_connect(G_OBJECT(tog),"toggled",
    G_CALLBACK(selectColorCallback),const_cast<GTKColorEditor*>(this));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::colorEditorCallback(GtkWidget *w,
gpointer user_data) {
  GTKColorEditor* ce=reinterpret_cast<GTKColorEditor*>(user_data);

  if (ce->form) { 
    gtk_widget_map(ce->form); 
  } else {
    ce->form=gtk_window_new(GTK_WINDOW_TOPLEVEL);
    char name[80];
    snprintf(name,80,"GTKColorEditor: %s",
      ce->gtkcolormap->getPalette()->getName());
    gtk_window_set_title(GTK_WINDOW(ce->form),name);

    GtkWidget *vbox=gtk_vbox_new(FALSE,0);
    gtk_container_add(GTK_CONTAINER(ce->form),vbox);
    guint vbox_width=420;

    GtkWidget *frame=gtk_frame_new("colormap colors and rgb graph");
    gtk_container_set_resize_mode(GTK_CONTAINER(frame),
      GTK_RESIZE_IMMEDIATE);
    gtk_container_add(GTK_CONTAINER(vbox),frame);

    ce->drawing_area=gtk_drawing_area_new();
    gtk_widget_set_size_request(ce->drawing_area,vbox_width,100);
    gtk_container_add(GTK_CONTAINER(frame),ce->drawing_area);
    gtk_widget_show(ce->drawing_area);
    gtk_widget_set_events(ce->drawing_area, GDK_EXPOSURE_MASK |
      GDK_BUTTON_MOTION_MASK |
      GDK_BUTTON_PRESS_MASK |
      GDK_BUTTON_RELEASE_MASK);
  //don't know what signal to use for resize
    gtk_signal_connect_after(GTK_OBJECT(ce->drawing_area),"realize",
      G_CALLBACK(realizeCallback),user_data);
    gtk_signal_connect(GTK_OBJECT(ce->drawing_area),"expose_event",
      G_CALLBACK(exposeCallback),user_data);
    gtk_signal_connect(GTK_OBJECT(ce->drawing_area),
      "button_press_event",G_CALLBACK(buttonPressCallback),user_data);
    gtk_widget_show(frame);

    guint cols=6;
    guint rows=ce->num_cmap_colors/cols;
    if (rows*cols<ce->num_cmap_colors) rows++;
    guint button_width=vbox_width/cols;
    guint button_height=20;

    GtkWidget *colors=gtk_frame_new("colormap index selector");
    gtk_widget_set_size_request(colors,button_width*cols,
      button_height*cols);
    gtk_container_add(GTK_CONTAINER(vbox),colors);

    GtkWidget *panel=gtk_table_new(rows,cols,TRUE);
    gtk_container_add(GTK_CONTAINER(colors),panel);


    ce->toggle[0]=gtk_radio_button_new_with_label(0,"0");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ce->toggle[0]),TRUE);
    ce->configureToggle(panel,ce->toggle[0],
      button_width,button_height,rows,0);
    ce->selected=0;
    gtk_widget_show(ce->toggle[0]);
    for (int i=1;i<ce->num_cmap_colors;i++) {
      snprintf(name,80,"%d",i);
      ce->toggle[i]=gtk_radio_button_new_with_label_from_widget(
        GTK_RADIO_BUTTON(ce->toggle[0]),name);
      ce->configureToggle(panel,ce->toggle[i],
        button_width,button_height,rows,i);
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ce->toggle[i]),
        FALSE);
      gtk_widget_show(ce->toggle[i]);
    }
    gtk_widget_show(panel);
    gtk_widget_show(colors);

    ce->swatch=gtk_color_selection_dialog_new("colorSelector");
    g_signal_connect(GTK_COLOR_SELECTION(GTK_COLOR_SELECTION_DIALOG(
       ce->swatch)->colorsel),
      "color-changed",G_CALLBACK(changeColorCallback),user_data);
    g_signal_connect(GTK_COLOR_SELECTION_DIALOG(ce->swatch)->ok_button,
      "clicked",G_CALLBACK(okColorSelectionCallback),user_data);
    g_signal_connect(
      GTK_COLOR_SELECTION_DIALOG(ce->swatch)->cancel_button,
      "clicked",G_CALLBACK(cancelColorSelectionCallback),user_data);
    gtk_widget_show(ce->swatch);
    gtk_widget_show(vbox);
    gtk_widget_show(ce->form);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GdkColor GTKColorEditor::allocNamedColor(char* colorname,
const GdkColor &default_color){
  GdkColor color;
  gint status=gdk_color_parse(colorname,&color);
  bool writeable=TRUE;
  bool best_match=TRUE;
  bool success=gdk_colormap_alloc_color(gtkcolormap->getColormap(),
    &color,writeable,best_match);
  return (success ? color : default_color);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::exposeSwatchCallback(GtkWidget *w,
GdkEventExpose *event,gpointer user_data) {
  GTKColorEditor *ce=reinterpret_cast<GTKColorEditor*>(user_data);
  int i=ce->selected;
  if (i!=-1) {
    GtkRequisition requisition;
    gtk_widget_size_request(w,&requisition);
    gdk_gc_set_foreground(ce->gc,&ce->gtkcolormap->getColor(i));
    gdk_draw_rectangle(w->window,ce->gc,TRUE,0,0,
      requisition.width,requisition.height);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::drawRGBBars(int start,int finish) const {
  GtkRequisition requisition;
  gtk_widget_size_request(drawing_area,&requisition);

  gint bar_width=static_cast<gint>(
    static_cast<double>(requisition.width-num_cmap_colors-1)
   /static_cast<double>(num_cmap_colors));
  gint bar_height=requisition.height;

  int beginning=max(0,start);
  int end=min(num_cmap_colors-1,finish);
  for (int i=beginning,x=1+beginning*(bar_width+1);i<=end;
  i++,x+=bar_width+1) {
    gdk_gc_set_foreground(gc,&gtkcolormap->getColor(i));
    gdk_draw_rectangle(drawing_area->window,gc,TRUE,x,0,
      bar_width,bar_height);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::drawRGBCurves(int start,int finish) const {
  GtkRequisition requisition;
  gtk_widget_size_request(drawing_area,&requisition);

  gint bar_width=static_cast<gint>(
    static_cast<double>(requisition.width-num_cmap_colors-1)
   /static_cast<double>(num_cmap_colors));
  gint bar_height=requisition.height;

  int beginning=(start>0 ? start-1 : 0);
  int end=(finish<num_cmap_colors-1 ? finish+1 : num_cmap_colors-1);
  double r=static_cast<double>(bar_height)/65535.;
  gboolean writeable=TRUE;
  gboolean best_match=TRUE;

#define DRAW_COLOR_CURVE(rgb,name) \
  { \
    GdkColor rgb ## _color; \
    gint status=gdk_color_parse(name,&rgb ## _color); \
    gdk_gc_set_rgb_fg_color(gc,&rgb ## _color); \
    int y=bar_height-static_cast<gint>(r* \
      static_cast<double>(gtkcolormap->getColor(beginning).rgb)); \
    for (int i=beginning+1,x=beginning*(bar_width+1)+bar_width/2; \
    i<=end;i++,x+=bar_width+1) { \
      gint yn=bar_height-static_cast<gint>(r* \
        static_cast<double>(gtkcolormap->getColor(i).rgb)); \
      gdk_draw_line(drawing_area->window,gc,x,y,x+bar_width+1,yn); \
      y=yn; \
    } \
  }

  DRAW_COLOR_CURVE(red,"red")
  DRAW_COLOR_CURVE(green,"green")
  DRAW_COLOR_CURVE(blue,"blue")
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::exposeCallback(GtkObject *o,GdkEventExpose *event,
gpointer user_data) {
  GtkWidget *w=reinterpret_cast<GtkWidget*>(o);
  GTKColorEditor *ce=reinterpret_cast<GTKColorEditor*>(user_data);
//should have w==ce->drawing_area

  GtkRequisition requisition;
  gtk_widget_size_request(w,&requisition);

  gint bar_width=static_cast<gint>(
    static_cast<double>(requisition.width-ce->num_cmap_colors-1)
   /static_cast<double>(ce->num_cmap_colors));
  gint bar_height=requisition.height;

  GdkGCValues values;
  gdk_gc_get_values(ce->gc,&values);
  gdk_gc_set_foreground(ce->gc,&values.background);
  gdk_draw_rectangle(w->window,ce->gc,TRUE,0,0,requisition.width,
    requisition.height);

  ce->drawRGBBars(0,ce->num_cmap_colors);
  ce->drawRGBCurves(0,ce->num_cmap_colors);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::buttonPressCallback(GtkWidget *w,
GdkEventButton *event,gpointer user_data) {
  GTKColorEditor *ce=reinterpret_cast<GTKColorEditor*>(user_data);

  GtkRequisition requisition;
  gtk_widget_size_request(w,&requisition);

  gint bar_width=static_cast<gint>(
    static_cast<double>(requisition.width-ce->num_cmap_colors-1)
   /static_cast<double>(ce->num_cmap_colors));
  int i=static_cast<int>(event->x/static_cast<gdouble>(bar_width+1));

  if (i!=ce->selected) {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(ce->toggle[i]),TRUE);
  }

  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(
    GTK_COLOR_SELECTION_DIALOG(ce->swatch)->colorsel),
    &ce->gtkcolormap->getColor(i));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::okColorSelectionCallback(GtkWidget *w,
gpointer user_data) {
  GTKColorEditor *ce=reinterpret_cast<GTKColorEditor*>(user_data);
  if (ce->warnUserNoColor(ce->swatch)) return;

  GdkColor new_color;
  gtk_color_selection_get_current_color(GTK_COLOR_SELECTION(
    GTK_COLOR_SELECTION_DIALOG(ce->swatch)->colorsel),&new_color);
  ce->gtkcolormap->setColor(ce->selected,new_color);
#ifdef DEBUG
//cout << "\tGdkColormap = " << ce->gtkcolormap->getColormap() << endl;
#endif

  ce->drawRGBBars(ce->selected-1,ce->selected+1);
  ce->drawRGBCurves(ce->selected-1,ce->selected+1);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColorEditor::cancelColorSelectionCallback(GtkWidget *w,
gpointer user_data) {
  GTKColorEditor *ce=reinterpret_cast<GTKColorEditor*>(user_data);
  if (ce->warnUserNoColor(ce->swatch)) return;
  const GdkColor &color=ce->gtkcolormap->getColor(ce->selected);
//the following generates a call to changeColorCallback
  gtk_color_selection_set_current_color(GTK_COLOR_SELECTION(
    GTK_COLOR_SELECTION_DIALOG(ce->swatch)->colorsel),&color);
  gtk_color_selection_set_previous_color(GTK_COLOR_SELECTION(
    GTK_COLOR_SELECTION_DIALOG(ce->swatch)->colorsel),&color);

  ce->drawRGBBars(ce->selected-1,ce->selected+1);
  ce->drawRGBCurves(ce->selected-1,ce->selected+1);
}                                

void GTKColorEditorList::printOn(ostream &os) const {
  os << "\tGTKColorEditorList:" << endl;
}

template class ISLList<GTKColorEditor>;
INSTANTIATE_NUMPTR(GtkWidget*);
#endif
