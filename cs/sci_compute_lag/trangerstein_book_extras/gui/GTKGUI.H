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
#ifndef GTK_GUI_H
#define GTK_GUI_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GTKGUI.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#ifdef USE_GTK
#include <iostream>
#include <gtk/gtkwidget.h>
#include "ClassThread.H"
#include "GUIBase.H"
#include "Tracer.H"

class GTKGUIInputParameterList;
class GTKGUI : public GUIBase {
  public:
    class PullDownMenu;
  private:
    static gboolean workCallback(gpointer);
    GtkWidget *program_toplevel;
    GtkWidget *vbox;
    GtkWidget *menu_bar;
    PullDownMenu *quit_menu;
    GtkWidget *busy,*busy_label;
  protected:
    GTKGUIInputParameterList *main_list;
    PullDownMenu *view_menu;

#ifndef NO_THREAD
    void showBusyDialog(const char*,const char*);
#endif
    void hideBusyDialog() const { 
      if (busy) gtk_widget_hide(busy); 
    }
    void showWarningDialog(const char* str);
    
    GtkWidget* createPulldownMenu(GtkWidget*,const char*,bool radio=FALSE);

    static void cancelBusyCallback(GtkWidget*,gpointer);
    static void onGTKGUIDialogResponse(GtkWidget *dialog,gint arg1,
    gpointer user_data) {
      cancelBusyCallback(dialog,user_data);
    }
    static void quitCallback(GtkWidget*,gpointer);
    virtual void quit(GtkWidget*);
    static void okBusyCallback(GtkWidget*,gpointer);
    static void runMainCallback(GtkWidget*,gpointer); 
    static void readMainMenuCallback(GtkWidget*,gpointer);
    static void checkMainInputCallback(GtkWidget*,gpointer user_data) {
      GTKGUI *gui=reinterpret_cast<GTKGUI*>(user_data);
      if (gui->checkMainInput) gui->checkMainInput();
    }
    virtual GtkWidget* createMainWindow1(int,char*[]);
    virtual void createMainWindow2(GtkWidget*);
  public:
    GTKGUI(char*,char*,GTKGUIInputParameterList*,void (*)(bool),
      void (*)()=0,void (*)()=0,void (*)()=0,bool=FALSE);
    virtual ~GTKGUI();
    virtual void createMainWindow(int argc,char *argv[]) {
      createMainWindow2(createMainWindow1(argc,argv));
    }
    void createFileMenu();
    void createViewMenu();
    GtkWidget* getTopLevelWidget() const { return program_toplevel; }
    void eventLoop();
    virtual void printOn(ostream &os = cerr) const;
};

class GTKGUI::PullDownMenu {
  private:
    GtkWidget *menu;
  public:
    PullDownMenu(GtkWidget*,const char*,bool=FALSE);
    ~PullDownMenu();
    GtkWidget* createLabelbutton(const char*);
    GtkWidget* createPushButton(const char*,
      void (*)(GtkWidget*,gpointer),gpointer);
    GtkWidget* createToggleButton(const char*,
      void (*)(GtkWidget*,gpointer),gpointer,bool);
    GtkWidget* getWidget() { return menu; }
    const GtkWidget* getWidget() const { return menu; }
    void show() { gtk_widget_show(menu); }
};
#endif
#endif
