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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIEnumInputParameter.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "GUIEnumInputParameter.H"
#include "Tracer.H"
#ifdef USE_GTK
#include <gtk/gtkhbox.h>
#include <gtk/gtklabel.h>
#include <gtk/gtkradiobutton.h>
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GtkWidget* GUIEnumInputParameter::createWidget(GtkWidget *parent) {
  hbox=gtk_hbox_new(TRUE,upper_bound);
  GtkWidget *label=gtk_label_new(getName());
  gtk_widget_set_size_request(label,50,15);
  gtk_container_add(GTK_CONTAINER(hbox),label);
  gtk_widget_show(label);

  GtkWidget *w=gtk_radio_button_new_with_label(0,value_names[0]);
  widgets[0]=w;
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widgets[0]),FALSE);
  gtk_container_add(GTK_CONTAINER(hbox),widgets[0]);
  gtk_widget_show(widgets[0]);
  for (int i=1;i<=upper_bound;i++) {
    widgets[i]=gtk_radio_button_new_with_label_from_widget(
      GTK_RADIO_BUTTON(widgets[0]),value_names[i]);
    g_signal_connect(G_OBJECT(widgets[i]),"toggled",
      G_CALLBACK(GTKGUIVirtualInput::onGTKGUIVirtualInputChangedEvent),
      this);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widgets[i]),FALSE);
    gtk_container_add(GTK_CONTAINER(hbox),widgets[i]);
    gtk_widget_show(widgets[i]);
  }
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widgets[*ptr_data]),
    TRUE);
  gtk_container_add(GTK_CONTAINER(parent),hbox);
  gtk_widget_show(hbox);
  return hbox; 
}

void GUIEnumInputParameter::valueChanged(GtkWidget *w) {
  GTKGUIVirtualInput::valueChanged(w);
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
    *ptr_data=-1;
    for (int i=0;i<=upper_bound;i++) {
      if (w==widgets[i]) { *ptr_data=i; break; }
    }
    CHECK_NONNEGATIVE(*ptr_data);
  } else {
    *ptr_data=0;
    if (w==widgets[0]) *ptr_data=1;
  }
  for (int i=0;i<=upper_bound;i++) {
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widgets[i]),FALSE);
  }
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widgets[*ptr_data]),
    TRUE);
}
#else
#include <Xm/RowColumn.h>
#include <Xm/TextFP.h>
#include <Xm/ToggleB.h>
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Widget GUIEnumInputParameter::createWidget(Widget parent) {
  Widget radiobox = XtVaCreateWidget("radiobox",
                                xmRowColumnWidgetClass,
                                parent,
                                XmNorientation, XmHORIZONTAL,
                                XmNpacking, XmPACK_COLUMN,
                                XmNradioBehavior, True,
                                XmNradioAlwaysOne, True,
                                0 );

  for (int i=0;i<=upper_bound;i++) {
    XmString xm_string=XmStringCreateSimple(value_names[i]);
    widgets[i] = XtVaCreateManagedWidget (value_names[i],
                                 xmToggleButtonWidgetClass,
                                 radiobox,
                                 XmNlabelType, XmSTRING,
                                 XmNlabelString,xm_string,
                                 0 );
    XmStringFree(xm_string);
    XtAddCallback(widgets[i],
                  XmNvalueChangedCallback,
                  &GUIVirtualInput::valueChangedCallback,
                  reinterpret_cast<XtPointer>(this) );
    XmToggleButtonSetState(widgets[i],False,False);
  }
  XmToggleButtonSetState(widgets[*ptr_data],True,False);

  return radiobox; 
}

void GUIEnumInputParameter::valueChanged(Widget w) {
  GUIVirtualInput::valueChanged(w);
  if (XmToggleButtonGetState(w)) {
    *ptr_data=-1;
    for (int i=0;i<=upper_bound;i++) {
      if (w==widgets[i]) { *ptr_data=i; break; }
    }
    CHECK_NONNEGATIVE(*ptr_data);
  }
}
#endif

#include "NumPtr.C"
INSTANTIATE_NUMPTR(char*)
INSTANTIATE_NUMPTR(GUI_WIDGET_POINTER)
