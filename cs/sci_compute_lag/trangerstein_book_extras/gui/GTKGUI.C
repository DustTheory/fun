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
// "$Header: /home/faculty/johnt/cvs/gui/GTKGUI.C,v 1.1 2005/06/08 11:19:01 johnt Exp $"
#ifdef USE_GTK
#include "GTKGUI.H"
#include <gtk/gtkbutton.h>
#include <gtk/gtkdialog.h>
#include <gtk/gtkentry.h>
#include <gtk/gtkframe.h>
#include <gtk/gtkhbbox.h>
#include <gtk/gtkhbox.h>
#include <gtk/gtklabel.h>
#include <gtk/gtkmenubar.h>
#include <gtk/gtkmenuitem.h>
#include <gtk/gtkmessagedialog.h>
#include <gtk/gtkstock.h>
#include <gtk/gtktogglebutton.h>
#include <gtk/gtkvbox.h>
#include <gtk/gtkwindow.h>
#include "Const.H"
#include "InputParameter.H"
#include "GTKGUIVirtualInput.H"
#include "ClassThread.H"
#include "Tracer.H"

//this code does not handle
//  2) changing number of palette entries
//  3) disabling palette and plot var widgets for restarts

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGUI::GTKGUI(char* n,char* dn,GTKGUIInputParameterList *ipl,
void (*m)(bool),void (*cm)(),void (*cu)(),void (*sd)(),bool ut) :
GUIBase(dn,m,cm,cu,sd,ut),program_toplevel(0),menu_bar(0),quit_menu(0),
vbox(0),busy(0),view_menu(0),main_list(ipl) {
  int argc=0;
  char **argv=0;
  GdkDisplay *display=0;
  if (display_name) GdkDisplay *display=gdk_display_open(display_name);
  if (display==0) {
    GdkDisplayManager *display_manager=gdk_display_manager_get();
    display=gdk_display_manager_get_default_display(display_manager);
  }
  ASSERT(display!=0);

  program_toplevel=gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(program_toplevel),n);

  vbox=gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(program_toplevel),vbox);

  menu_bar=gtk_menu_bar_new();
  gtk_box_pack_start(GTK_BOX(vbox),menu_bar,FALSE,FALSE,2);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGUI::~GTKGUI() {
  if (quit_menu) delete quit_menu; quit_menu=0;
  if (view_menu) delete view_menu; view_menu=0;
  if (busy) gtk_widget_destroy(busy); busy=0;
  if (vbox) gtk_widget_destroy(vbox); vbox=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//XEvent types are defined in X.h
void GTKGUI::eventLoop() {
  while (!quit_called) {
    g_main_iteration(FALSE);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GtkWidget* GTKGUI::createMainWindow1(int argc, char *argv[]) {
  char cmdline[LENGTH_COMMENT];
  cmdline[0]='\0';
  strcat(cmdline,argv[0]);
  for (int i=1;i<argc;i++) {
    strcat(cmdline," ");
    strcat(cmdline,argv[i]);
  }

  GtkWidget *textfield= gtk_entry_new();
  gtk_entry_set_text(GTK_ENTRY(textfield),cmdline);
  gtk_container_add(GTK_CONTAINER(vbox),textfield);
  gtk_widget_show(textfield);
  return vbox;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::createMainWindow2(GtkWidget*) {
  GtkWidget *button=gtk_button_new_with_label("Start Run Now");
  g_signal_connect(G_OBJECT(button),"clicked",
    G_CALLBACK(&GTKGUI::runMainCallback),this);
  gtk_container_add(GTK_CONTAINER(vbox),button);
  gtk_widget_show(button);

  gtk_widget_show(menu_bar);
  gtk_widget_show(vbox);
  gtk_widget_show(program_toplevel);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGUI::PullDownMenu::PullDownMenu(GtkWidget *parent_menu_bar,
const char *name,bool radio) {
  GtkWidget *menu_bar_item=gtk_menu_item_new_with_label(name);
  menu=gtk_menu_new();
  gtk_widget_show(menu_bar_item);
  gtk_menu_item_set_submenu(GTK_MENU_ITEM(menu_bar_item),menu);
  gtk_menu_shell_append(GTK_MENU_SHELL(parent_menu_bar),menu_bar_item);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKGUI::PullDownMenu::~PullDownMenu() {
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GtkWidget* GTKGUI::PullDownMenu::createLabelbutton(const char *name) {
  GtkWidget *label=gtk_label_new(name);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu),label);
  gtk_widget_show(label);
  return label;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GtkWidget* GTKGUI::PullDownMenu::createPushButton(const char *name,
void (*callback)(GtkWidget*,gpointer),gpointer user_data) {
  GtkWidget *push=gtk_menu_item_new_with_label(name);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu),push);
  g_signal_connect(G_OBJECT(push),"activate",
    G_CALLBACK(callback),user_data);
  gtk_widget_show(push);
  return push;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GtkWidget* GTKGUI::PullDownMenu::createToggleButton(const char *name,
void (*callback)(GtkWidget*,gpointer),gpointer user_data,bool set) {
  GtkWidget *tog=gtk_toggle_button_new_with_label(name);
  gtk_menu_shell_append(GTK_MENU_SHELL(menu),tog);
  g_signal_connect(G_OBJECT(tog),"toggled",
    G_CALLBACK(callback),user_data);
  gtk_widget_show(tog);
  return tog;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::createFileMenu() {
  quit_menu=OPERATOR_NEW PullDownMenu(menu_bar,"Quit");
  quit_menu->createPushButton("Quit",&GTKGUI::quitCallback,this);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::quitCallback(GtkWidget *w,gpointer user_data) {
  reinterpret_cast<GTKGUI*>(user_data)->quit(w);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::quit(GtkWidget *w) {
  if (shutdown) shutdown();
  quit_called=TRUE;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::createViewMenu() {
  view_menu=OPERATOR_NEW PullDownMenu(menu_bar,"View");
  view_menu->createPushButton("Main",&GTKGUI::readMainMenuCallback,
    this);
  main_list->
    createListWin(GTKGUI::checkMainInputCallback,this,"Verify");
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::readMainMenuCallback(GtkWidget*,gpointer user_data) {
  GTKGUI *gui=reinterpret_cast<GTKGUI*>(user_data);
  gui->main_list->reload();
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::showWarningDialog(const char* str) {
  if (program_toplevel==0) return;
  GtkWidget *warn=gtk_message_dialog_new(0,
      GTK_DIALOG_DESTROY_WITH_PARENT,GTK_MESSAGE_ERROR,
      GTK_BUTTONS_CLOSE,str);
  gtk_dialog_run(GTK_DIALOG(warn));
  gtk_widget_destroy(warn);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#ifndef NO_THREAD
//called from runMainCallback
void GTKGUI::showBusyDialog(const char *str,const char *str2) {
  if (!use_thread) return;
        
  busy=gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(busy),"Busy Dialog");

  GtkWidget *busy_vbox=gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(busy),busy_vbox);

  if (run_main_thread->isAlive()) {
    GtkWidget *busy_hbox=gtk_hbox_new(FALSE,0);
    gtk_container_add(GTK_CONTAINER(busy_vbox),busy_hbox);
    GtkWidget *busy_stop=gtk_button_new_with_label(str2);
    g_signal_connect(G_OBJECT(busy_stop),"clicked",
      G_CALLBACK(cancelBusyCallback),this);
    gtk_container_add(GTK_CONTAINER(busy_hbox),busy_stop);
    gtk_widget_show(busy_stop);
    gtk_widget_show(busy_hbox);
  } else {
    busy_label=gtk_label_new(str);
    gtk_container_add(GTK_CONTAINER(busy_vbox),busy_label);
    gtk_widget_show(busy_label);
  }

  gtk_widget_show(busy_vbox);
  gtk_widget_show(busy);
}
#endif
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//clicked "Start Run Now" on button in createMainWindow
void GTKGUI::runMainCallback(GtkWidget*,gpointer user_data) {
  GTKGUI *gui=reinterpret_cast<GTKGUI*>(user_data);

  gui->main_list->unsetBadparams();
  gui->checkMainInput();
  if (gui->main_list->getBadparams()) {
    gui->showWarningDialog("Parameter Error: Run Cancelled");
    return;
  }

#ifndef NO_THREAD
  if (gui->use_thread) {
    gui->showBusyDialog("Running...","Stop and Cleanup");
  }
#endif
  gui->runMainProc();
  g_idle_add(&GTKGUI::workCallback,gui);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
gboolean GTKGUI::workCallback(gpointer user_data) {
  GTKGUI *gui=reinterpret_cast<GTKGUI*>(user_data);
  if (!gui) return FALSE;
  gboolean run_finished=gui->mainFinished();
  if (!gui->use_thread) {
    gui->cleanup();
  } else if (run_finished) {
    if (gui->busy) gtk_widget_destroy(gui->busy); gui->busy=0;
    GtkWidget *dialog=gtk_dialog_new_with_buttons("Cleanup",0,
      GTK_DIALOG_DESTROY_WITH_PARENT,"Cleanup",GTK_RESPONSE_ACCEPT,
      0);
    GtkWidget *dialog_label=gtk_label_new("Cleanup");
    g_signal_connect(dialog,"response",
      G_CALLBACK(GTKGUI::onGTKGUIDialogResponse),gui);
    gtk_container_add(GTK_CONTAINER(GTK_DIALOG(dialog)->vbox),
      dialog_label);
    gtk_dialog_run(GTK_DIALOG(dialog));
    gtk_widget_destroy(dialog);
  }
  return !run_finished;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::okBusyCallback(GtkWidget *w,gpointer user_data) {
  if (!user_data) return;
  GTKGUI *gui=reinterpret_cast<GTKGUI*>(user_data);
#ifndef NO_THREAD
  if (!gui->mainRunning()) cancelBusyCallback(w,user_data);
  else if (gui->run_main_thread->isSuspended()) {
    gtk_label_set_label(GTK_LABEL(gui->busy_label),"Stop");
    gui->resumeRun();
  } else {
    gtk_label_set_label(GTK_LABEL(gui->busy_label),"Resume");
    gui->suspendRun();
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::cancelBusyCallback(GtkWidget*,gpointer user_data) {
  if (!user_data) return;
  GTKGUI *gui=reinterpret_cast<GTKGUI*>(user_data);
#ifndef NO_THREAD
  if (gui->run_main_thread) {
    ASSERT(gui->run_main_thread->finished());
    gui->run_main_thread->resetStatus();
  }
  gui->cleanup();
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKGUI::printOn(ostream &os) const {
  os << "GTKGUI: main_list = " << endl;
  os << *main_list;
  GUIBase::printOn(os);
}
#endif
