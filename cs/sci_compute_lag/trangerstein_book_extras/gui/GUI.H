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
#ifndef GUI_H
#define GUI_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUI.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#ifndef USE_GTK
#include <iostream>
extern "C" {
#include <Xm/ToggleB.h>
}
#include "ClassThread.H"
#include "GUIBase.H"

class GUIInputParameterList;
class GUI : public GUIBase {
  private:
    XtAppContext program_app;
    static bool workCallback(XtPointer);
    Widget _menu_bar;
    Widget _file_menu;
    Widget _program_toplevel;
    Widget _busy,_warn;
  protected:
    Widget view_menu;
    Widget _program_window;
    GUIInputParameterList *main_list;

    void showBusyDialog(const char* str);
    void hideBusyDialog() const {if(_busy) XtUnmanageChild(_busy);}
    void showWarningDialog(const char* str);
//  menu util functions
    
    Widget createPulldownMenu(Widget ,const char * ,bool radio=FALSE);
    Widget createOptionMenu(Widget parent,char *name);

//  menu callbacks
#define CALLBACK_MACRO(c,r)                                            \
    static void r ## Callback(Widget w,XtPointer client_data,          \
    XtPointer call_data) {                                             \
      reinterpret_cast<c*>(client_data)->r(w,call_data);               \
    }
    virtual void quit(Widget,XtPointer) {
      if (shutdown) shutdown();
      quit_called=TRUE; 
    }
    static void runMainCallback(Widget,XtPointer,XtPointer); 
    static void cancelBusyCallback(Widget,XtPointer,XtPointer);
    virtual void cancelBusy();
    static void okBusyCallback(Widget,XtPointer,XtPointer);
    virtual void okBusy();

//  main window
    void createBboard(Widget _program_window);
    
//  main window callbacks
    virtual void cmdAct(Widget,XtPointer){;}
    void cancel(Widget w,XtPointer) { XtDestroyWidget(w); }
    static void readMainMenuCallback(Widget,XtPointer,XtPointer);
    static void checkMainInputCallback(Widget,XtPointer client_data,
    XtPointer) {
      reinterpret_cast<GUI*>(client_data)->checkMainInput();
    }

    void destroy(Widget,XtPointer client_data,XtPointer){
      XtUnmanageChild(reinterpret_cast<Widget>(client_data));
    }
    Widget createMainWindow1(int,char*[]);
    void createMainWindow2(Widget);
CALLBACK_MACRO(GUI,cmdAct)
CALLBACK_MACRO(GUI,cancel)
CALLBACK_MACRO(GUI,quit)
  public:
    GUI(char*,char*,GUIInputParameterList*,void (*)(bool),
      void (*)()=0,void (*)()=0,void (*)()=0,bool=FALSE);
    virtual ~GUI();
    virtual void createMainWindow(int,char*[]);
    void createFileMenu();
    void createViewMenu();
    void createSeparator(Widget parent);
    void createLabelbutton(Widget,const char*);
    Widget createPushbutton(Widget,const char*,XtCallbackProc,
      XtPointer);
    Widget createTogglebutton(Widget,const char *,XtCallbackProc ,
                              XtPointer,bool=FALSE);

    Widget getTopLevelWidget() const { return _program_toplevel; }
    void eventLoop();
    virtual void printOn(ostream &os = cerr) const;
};
#endif
#endif
