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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIInputParameter.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "GUIInputParameter.H"
#include "Tracer.H"
#ifdef USE_GTK
#include <gtk/gtkhbox.h>
#include <gtk/gtklabel.h>
#include <gtk/gtkradiobutton.h>
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> GtkWidget* GUIInputParameter<bool>::createWidget(
GtkWidget *parent) {
  hbox=gtk_hbox_new(TRUE,2);
  GtkWidget *label=gtk_label_new(getName());
  gtk_widget_set_size_request(label,50,15);
  gtk_container_add(GTK_CONTAINER(hbox),label);
  gtk_widget_show(label);

  GtkWidget *radio_buttont=gtk_radio_button_new_with_label(0,"True");
  GtkWidget *radio_buttonf=gtk_radio_button_new_with_label_from_widget(
    GTK_RADIO_BUTTON(radio_buttont),"False");
  g_signal_connect(G_OBJECT(radio_buttont),"toggled",
    G_CALLBACK(GTKGUIVirtualInput::onGTKGUIVirtualInputChangedEvent),
    this);
  if (*ptr_data) {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(radio_buttont),TRUE);
  } else {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(radio_buttonf),TRUE);
  }
  gtk_box_pack_start(GTK_BOX(hbox),radio_buttont,TRUE,TRUE,2);
  gtk_box_pack_start(GTK_BOX(hbox),radio_buttonf,TRUE,TRUE,2);
  gtk_container_add(GTK_CONTAINER(parent),hbox);
  gtk_widget_show(radio_buttont);
  gtk_widget_show(radio_buttonf);
  gtk_widget_show(hbox);

  return hbox; 
}
template<> void GUIInputParameter<bool>::reload() {
}

template<> void GUIInputParameter<bool>::valueChanged(GtkWidget *w) {
  GTKGUIVirtualInput::valueChanged(w);
  *ptr_data=gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> GtkWidget* GUIInputParameter<int>::createWidget(
GtkWidget* parent) {
  return GTKGUIVirtualInput::createWidget(parent);
}
template<> void GUIInputParameter<int>::reload() {
  GTKGUIVirtualInput::reload();
}

template<> void GUIInputParameter<int>::valueChanged(GtkWidget* w) {
  GTKGUIVirtualInput::valueChanged(w);
  const char *text=gtk_entry_get_text(GTK_ENTRY(w));
  int t=atoi(text);
  if (t>=lower_bound && t<=upper_bound) *ptr_data=t;
  else gdk_beep();
  writeTo(w);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> GtkWidget* GUIInputParameter<double>::createWidget(
GtkWidget* parent) {
  return GTKGUIVirtualInput::createWidget(parent);
}
template<> void GUIInputParameter<double>::reload() {
  GTKGUIVirtualInput::reload();
}

template<> void GUIInputParameter<double>::valueChanged(GtkWidget* w) {
  GTKGUIVirtualInput::valueChanged(w);
  const char *text=gtk_entry_get_text(GTK_ENTRY(w));
  double t=atof(text);
  if (t>=lower_bound && t<=upper_bound) *ptr_data=t;
  else gdk_beep();
  writeTo(w);
}
#else
#include <Xm/RowColumn.h>
#include <Xm/TextFP.h>
#include <Xm/ToggleB.h>
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> Widget GUIInputParameter<bool>::createWidget(
Widget parent) {
  Widget radiobox = XtVaCreateWidget("radiobox",
                                xmRowColumnWidgetClass,
                                parent,
                                XmNorientation, XmHORIZONTAL,
                                XmNpacking, XmPACK_COLUMN,
                                XmNradioBehavior, True,
                                XmNradioAlwaysOne, True,
                                0 );

  XmString xm_string=XmStringCreateSimple(const_cast<char*>("True"));
  Widget toggle1 = XtVaCreateManagedWidget ( "toggle1",
                               xmToggleButtonWidgetClass,
                               radiobox,
                               XmNlabelType, XmSTRING,
                               XmNlabelString,xm_string,
                               0 );
  XmStringFree(xm_string);

  XtAddCallback(toggle1,
                XmNvalueChangedCallback,
                &GUIVirtualInput::valueChangedCallback,
                reinterpret_cast<XtPointer>(this) );

  xm_string=XmStringCreateSimple(const_cast<char*>("False"));
  Widget toggle2 = XtVaCreateManagedWidget  ( "toggle2",
                               xmToggleButtonWidgetClass,
                               radiobox,
                               XmNlabelType, XmSTRING,
                               XmNlabelString,xm_string,
                               0 );
  XmStringFree(xm_string);

  XtAddCallback(toggle2,
                XmNvalueChangedCallback,
                &GUIVirtualInput::valueChangedCallback,
                reinterpret_cast<XtPointer>(this) );

  if (*ptr_data) XmToggleButtonSetState(toggle1,True,False);
  else XmToggleButtonSetState(toggle2,True,False);

  return radiobox; 
}
template<> void GUIInputParameter<bool>::reload() {
}

template<> void GUIInputParameter<bool>::valueChanged(Widget w) {
  value_changed=TRUE;
  Widget toggle=w;
  if (!XmToggleButtonGetState(toggle)) return;
  XmString xm_string;
  char *label;
  XtVaGetValues(toggle,XmNlabelString,&xm_string,0);
  XmStringGetLtoR(xm_string,XmFONTLIST_DEFAULT_TAG,&label);
  XmStringFree(xm_string);
  if (!strcmp(label,"True")) *ptr_data=1;
  else *ptr_data=0;
  XtFree(label);
  return;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void GUIInputParameter<int>::valueChanged(Widget w) {
  value_changed=TRUE;
  char *text=XmTextFieldGetString(w);
  int t=atoi(text);
  if (t>=lower_bound && t<=upper_bound) *ptr_data=t;
  else XBell(XtDisplay(w),40);
  writeTo(w);
  XtFree(text);
}
template<> void GUIInputParameter<int>::reload() {
  GUIVirtualInput::reload();
}

template<> Widget GUIInputParameter<int>::createWidget(
Widget parent) {
  return GUIVirtualInput::createWidget(parent);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<> void GUIInputParameter<double>::valueChanged(Widget w) {
  value_changed=TRUE;
  char *text=XmTextFieldGetString(w);
  double t=atof(text);
  if (t>=lower_bound && t<=upper_bound) *ptr_data=t;
  else XBell(XtDisplay(w),40);
  writeTo(w);
  XtFree(text);
}
template<> void GUIInputParameter<double>::reload() {
  GUIVirtualInput::reload();
}

template<> Widget GUIInputParameter<double>::createWidget(Widget parent) {
  return GUIVirtualInput::createWidget(parent);
}
#endif
template class GUIInputParameter<bool>;
template class GUIInputParameter<int>;
template class GUIInputParameter<double>;
