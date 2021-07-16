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
// "$Header: /home/faculty/johnt/cvs/gui/GTKGUIVirtualInput.C,v 1.1 2005/06/08 11:19:02 johnt Exp $"
#include "GTKGUIVirtualInput.H"
#include <gtk/gtkarrow.h>
#include <gtk/gtkbutton.h>
#include <gtk/gtkdialog.h>
#include <gtk/gtkhbox.h>
#include <gtk/gtkmessagedialog.h>
#include <gtk/gtkscrolledwindow.h>
#include <gtk/gtkhseparator.h>
#include <gtk/gtklabel.h>
#include <gtk/gtkvbox.h>
#include "MemoryDebugger.H"
#include "Tracer.H"
#include "NISLList.C"
//::::::::::::::::::::::: VirtualInput :::::::::::::::::::::::::::::::::
GTKGUIVirtualInput::GTKGUIVirtualInput(const char *g) : group_name(0),
hbox(0),group_arrow(0),group_label(0), /* label(0), */ textfield(0),
value_changed(FALSE) {
  group_name=OPERATOR_NEW_BRACKET(char,strlen(g)+1);
  strcpy(group_name,g);
}

GTKGUIVirtualInput::~GTKGUIVirtualInput() {
  delete group_name; group_name=0;
}

void GTKGUIVirtualInput::writeTo(GtkWidget *tf) const {
  gtk_entry_set_text(GTK_ENTRY(tf),getValue());
}

void GTKGUIVirtualInput::onGTKGUIVirtualInputChangedEvent(GtkWidget *w,
gpointer user_data) {
  reinterpret_cast<GTKGUIVirtualInput*>(user_data)->valueChanged(w);
}

void GTKGUIVirtualInput::onGTKGUIVirtualInputFocusOutEvent(GtkWidget *w,
GdkEventFocus*,gpointer user_data) {
  reinterpret_cast<GTKGUIVirtualInput*>(user_data)->valueChanged(w);
}

void GTKGUIVirtualInput::printOn(ostream &os) const {
  os << "GTKGUIVirtualInput: name = " << getName()
     << "\n\tgroup_name = " << group_name << endl;
}

GtkWidget* GTKGUIVirtualInput::createWidget(GtkWidget *parent) {
  hbox=gtk_hbox_new(FALSE,0);

  GtkWidget *label=gtk_label_new(getName());
  gtk_widget_set_size_request(label,50,15);
  gtk_container_add(GTK_CONTAINER(hbox),label);
  gtk_widget_show(label);

  textfield=gtk_entry_new();
  gtk_entry_set_text(GTK_ENTRY(textfield),getValue());
  gtk_container_add(GTK_CONTAINER(hbox),textfield);
  gtk_widget_show(textfield);

  g_signal_connect(G_OBJECT(textfield),"activate",
    G_CALLBACK(onGTKGUIVirtualInputChangedEvent),this);
  g_signal_connect(G_OBJECT(textfield),"focus_out_event",
    G_CALLBACK(onGTKGUIVirtualInputFocusOutEvent),this);

  gtk_container_add(GTK_CONTAINER(parent),hbox);
  gtk_widget_show(hbox);
  return hbox;
}

ostream& operator<<(ostream &os,const GTKGUIVirtualInput& vi) {
  os << "GTKGUIVirtualInput: group_name = " << vi.getGroupName() 
     << ", value = " << vi.getValue();
  return os;
}

void GTKGUIVirtualInput::reload() {
  gtk_entry_set_text(GTK_ENTRY(textfield),getValue());
}

//:::::::::::::: InputParameterList ::::::::::::::::::::::::::::::::::::
GTKGUIInputParameterList::GTKGUIInputParameterList(const char *n) :
InputParameterList(n),widget(0),warn(0), /* needreload(FALSE), */
badparams(FALSE),warn_up(FALSE),createListWinCallback(0),
createListWinCallbackData(0) {
}

void GTKGUIInputParameterList::onGTKGUIInputParameterListDoneClicked(
GtkWidget *w,gpointer user_data) {
  GTKGUIInputParameterList *gui_list=
    reinterpret_cast<GTKGUIInputParameterList*>(user_data);

  bool bad=gui_list->badparams;
  gui_list->badparams=FALSE;

  if (gui_list->createListWinCallback) {
    gui_list->
      createListWinCallback(w,gui_list->createListWinCallbackData);
  }

  if (!gui_list->badparams){
    gtk_widget_hide(gui_list->widget);
    gui_list->badparams=bad;
  } else {
    if (gui_list->warn_up) gui_list->showWarningDialog(0);
  }
}

void GTKGUIInputParameterList::arrowActCallback(GtkWidget *gl,
gpointer user_data) {

  GTKGUIInputParameterList* list=
    reinterpret_cast<GTKGUIInputParameterList*>(user_data);
  
  NISLListNode<VirtualInput> *p=0;
  GTKGUIVirtualInput *gp=0;
  for (p=list->first();p;p=list->next(p)) {
    gp=dynamic_cast<GTKGUIVirtualInput*>(p->getData());
    if (gp->groupLabelWidgetIs(gl)) { break; }
  }

  ASSERT(gp!=0);
  GtkWidget *ga=gp->getGroupArrow();
  GValue value={0,};
  gint adir=-1;
  g_value_init(&value,G_TYPE_INT);
  g_object_get_property(G_OBJECT(ga),"arrow-type",&value);
  adir=g_value_get_int(&value);
  g_value_unset(&value);

  if (adir==GTK_ARROW_RIGHT) { // to show all params in the group
    gtk_arrow_set(GTK_ARROW(ga),GTK_ARROW_DOWN,GTK_SHADOW_NONE);
    gtk_widget_show(gp->getParent());
  } else {
    gtk_arrow_set(GTK_ARROW(ga),GTK_ARROW_RIGHT,GTK_SHADOW_NONE);
    gtk_widget_hide(gp->getParent());
  }
}

void GTKGUIInputParameterList::createListWin(
void (*callback)(GtkWidget*,gpointer),gpointer data,
const char *call_back_button_label) {
  CHECK_TEST(widget==0);
  if (widget) {
    for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
      dynamic_cast<GTKGUIVirtualInput*>(p->getData())->reload();
    }
    gtk_widget_show(widget);
    return;
  }

  widget=gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(widget),getName());

  GtkWidget *vbox=gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(widget),vbox);

  GtkWidget *hbox=gtk_hbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(vbox),hbox);

  GtkWidget *done_button=gtk_button_new_with_label("Done");
  g_signal_connect(G_OBJECT(done_button),"clicked",
    G_CALLBACK(onGTKGUIInputParameterListDoneClicked),this);
  gtk_container_add(GTK_CONTAINER(hbox),done_button);
  gtk_widget_show(done_button);

  if (callback) {
    createListWinCallback=callback;
    createListWinCallbackData=data;
    GtkWidget *call_back_button=0;
    if (call_back_button_label) {
      call_back_button=
        gtk_button_new_with_label(call_back_button_label);
    } else {
      call_back_button=gtk_button_new_with_label("Verify");
    }
    g_signal_connect(G_OBJECT(call_back_button),"clicked",
      G_CALLBACK(callback),data);
    gtk_container_add(GTK_CONTAINER(hbox),call_back_button);
    gtk_widget_show(call_back_button);
  }
  gtk_widget_show(hbox);

  GtkWidget *separator=gtk_hseparator_new();
  gtk_container_add(GTK_CONTAINER(vbox),separator);
  gtk_widget_show(separator);

  GtkWidget *scrolled_window=gtk_scrolled_window_new(0,0);
  gtk_widget_set_size_request(scrolled_window,200,150);
  gtk_container_add(GTK_CONTAINER(vbox),scrolled_window);

  GtkWidget *viewport=gtk_viewport_new(0,0);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrolled_window),
    GTK_POLICY_AUTOMATIC,GTK_POLICY_AUTOMATIC);
  gtk_scrolled_window_add_with_viewport(
    GTK_SCROLLED_WINDOW(scrolled_window),viewport);

  GtkWidget *viewport_vbox=gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(viewport),viewport_vbox);

  char group_label[80];
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    bool done=FALSE;
    GTKGUIVirtualInput *gp=
      dynamic_cast<GTKGUIVirtualInput*>(p->getData());
    CHECK_POINTER(gp);
    for (NISLListNode<VirtualInput> *q=first();q!=p;q=next(q)) {
      GTKGUIVirtualInput *gq=
        dynamic_cast<GTKGUIVirtualInput*>(q->getData());
      if (gq->groupNameSameAs(gp)) { done=TRUE; break; }
    }
    if (done) continue;
    char *group_name=gp->getGroupName();
    if (!strcmp(group_name,"")) strcpy(group_label,"Misc Params");
    else strcpy(group_label,group_name);

    GtkWidget *hbox2=gtk_hbox_new(FALSE,0);
    gtk_container_add(GTK_CONTAINER(viewport_vbox),hbox2);

    GtkWidget *gq_group_arrow=
      gtk_arrow_new(GTK_ARROW_RIGHT,GTK_SHADOW_OUT);
    gtk_widget_set_size_request(gq_group_arrow,15,15);
    gtk_box_pack_start(GTK_BOX(hbox2),gq_group_arrow,FALSE,FALSE,0);
    gtk_widget_show(gq_group_arrow);

    GtkWidget *gq_group_label=gtk_button_new_with_label(group_label);
    gtk_button_set_alignment(GTK_BUTTON(gq_group_label),0.,0.5);
    g_signal_connect(G_OBJECT(gq_group_label),"clicked",
      G_CALLBACK(&GTKGUIInputParameterList::arrowActCallback),this);
    gtk_box_pack_start(GTK_BOX(hbox2),gq_group_label,FALSE,FALSE,0);
    gtk_widget_show(gq_group_label);
    gtk_widget_show(hbox2);

    GtkWidget *vbox2=gtk_vbox_new(FALSE,0);
    gtk_container_add(GTK_CONTAINER(viewport_vbox),vbox2);

    for (NISLListNode<VirtualInput> *q=p;q;q=next(q)) {
      GTKGUIVirtualInput *gq=
        dynamic_cast<GTKGUIVirtualInput*>(q->getData());
      CHECK_POINTER(gq);
      if (gq->groupNameSameAs(gp)) {
        gq->createWidget(vbox2);
        gq->setGroup(gq_group_arrow,gq_group_label);
      }
    }
  }
  gtk_widget_show(viewport_vbox);
  gtk_widget_show(viewport);
  gtk_widget_show(scrolled_window);
  gtk_widget_show(vbox);
  gtk_widget_show(widget);
}

void GTKGUIInputParameterList::reload() {
  CHECK_TEST(widget!=0);
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    dynamic_cast<GTKGUIVirtualInput*>(p->getData())->reload();
  }
  gtk_widget_show(widget);
}

void GTKGUIInputParameterList::showWarningDialog(const char* str) {
  if (str) badparams=TRUE;
  if (!widget) return;

  if (warn) {
    OBSOLETE("case not programmed");
  } else {
    warn=gtk_message_dialog_new(0,GTK_DIALOG_DESTROY_WITH_PARENT,
      GTK_MESSAGE_ERROR,GTK_BUTTONS_CLOSE,str);
  }
  warn_up=TRUE;
  gtk_dialog_run (GTK_DIALOG (warn));
  gtk_widget_destroy (warn); warn=0;
}

bool GTKGUIInputParameterList::getValueChanged() {
  if (!widget) { return TRUE; }
  
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    if (dynamic_cast<GTKGUIVirtualInput*>(p->getData())->valueChanged())
      return TRUE;
  }
  return FALSE;
}

void GTKGUIInputParameterList::unsetValueChanged() {
  if (!widget) return;
  for (NISLListNode<VirtualInput> *p=first();p;p=next(p)) {
    dynamic_cast<GTKGUIVirtualInput*>(p->getData())->
      unsetValueChanged();
  }
}

void GTKGUIInputParameterList::printOn(ostream &os) const {
  os << "GTKGUIInputParameterList\n" 
     << "\tbadparams = " << badparams
//   << "\tneedreload = " << needreload
     << "\twarn_up = " << warn_up << endl;
  InputParameterList::printOn(os);
}
