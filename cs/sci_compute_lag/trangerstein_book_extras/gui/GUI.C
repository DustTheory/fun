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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUI.C,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "GUI.H"
#ifndef USE_GTK
#include <stdlib.h> //for getenv
#include <unistd.h>// for sleep,pause
#include <X11/Intrinsic.h>
#include <Xm/ArrowB.h>
#include <Xm/BulletinB.h>
#include <Xm/CascadeB.h>
#include <Xm/DrawingA.h>
#include <Xm/DrawnB.h>
#include <Xm/FileSB.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/MainW.h>
#include <Xm/MessageB.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/Scale.h>
#include <Xm/ScrolledW.h>
#include <Xm/Text.h>
#include <Xm/TextFP.h>
#include <Xm/ToggleB.h>
#include <Xm/Form.h>
#include <Xm/Separator.h>

#include "Const.H"
#include "InputParameter.H"
#include "GUIVirtualInput.H"
#include "Palette.H"
#include "VC.H"
#include "XTWindow.H"
#include "ClassThread.H"
#include "Tracer.H"

//this code does not handle
//  2) changing number of palette entries
//  3) disabling palette and plot var widgets for restarts

//sgiMode makes the checkmarks for the togglebuttons
static const char* fallback_resources[]={
  "!*fontList: lucidasans-14",
  "*background: LightGrey",
  "*selectColor: red",
  "!*sgiMode: True",
  "! color editor labels:",
  "*selectorLabel*labelString: Palette",
  "*swatchLabel*labelString:   Current Color",
  "*redSlider*troughColor:   red",
  "*blueSlider*troughColor:  blue",
  "*greenSlider*troughColor: green",
  "*noColorWarningDialog*messageString: You must select a color\n\
    from the palette first.",
  "*togglePlotVarWarningDialog*messageString: Cannot change plot var in restart\n\
   Must select New Run in program window",
  0};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GUI::GUI(char* n,char* dn,GUIInputParameterList *ipl,
void (*m)(bool),void (*cm)(),void (*cu)(),void (*sd)(),
bool ut) : GUIBase(dn,m,cm,cu,sd,ut),program_app(0),
_menu_bar(0),_file_menu(0),_program_toplevel(0),_busy(0),_warn(0),
view_menu(0),_program_window(0),main_list(ipl) {
  int argc=0;
  char **argv=0;
  if (!display_name) {
    _program_toplevel=XtAppInitialize(&program_app,const_cast<char*>("Amr"),
      0,0,&argc,argv,const_cast<char**>(fallback_resources),0,0);
  } else {
    XtToolkitInitialize();
    program_app=XtCreateApplicationContext();
    XtAppSetFallbackResources(program_app,
      const_cast<char**>(fallback_resources));
    Display* display=XtOpenDisplay(program_app,
      const_cast<char*>(display_name),0,const_cast<char*>("Amr"),0,0,
      &argc,argv);
    ASSERT(display!=0);
  }
  XtVaSetValues(_program_toplevel,XmNtitle,n,0 );

  _program_window=XtVaCreateManagedWidget("Program",
     xmMainWindowWidgetClass,_program_toplevel,0);
  _menu_bar=
    XmCreateMenuBar(_program_window,const_cast<char*>("menu_bar"),0,0);
  XtManageChild(_menu_bar);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GUI::~GUI() {
  if (_busy) XtDestroyWidget(XtParent(_busy)); _busy=0;
  if (_program_window) XtDestroyWidget(_program_window);
    _program_window=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//XEvent types are defined in X.h
void GUI::eventLoop() {
  XEvent event;
  while (!quit_called) {
    XtAppNextEvent(program_app,&event);
    XtDispatchEvent(&event);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Widget GUI::createMainWindow1(int argc, char *argv[]) {
  Widget bboard=XtVaCreateWidget("bboard",
                                 xmBulletinBoardWidgetClass,
                                 _program_window, 
                                 XmNresizePolicy, XmRESIZE_GROW, 
                                 0); 
  Widget textfield=XtVaCreateManagedWidget("textfield",
                                           xmTextFieldWidgetClass,
                                           bboard, 
                                           XmNx, 30, 
                                           XmNy, 75, 
                                           XmNwidth, 360, 
                                           XmNeditable, False,
                                           0); 
  XtAddCallback(textfield,
                XmNactivateCallback,
                &GUI::cmdActCallback,
                reinterpret_cast<XtPointer>(this) ); 
  XtAddCallback(textfield,
                XmNlosingFocusCallback,
                &GUI::cmdActCallback,
                reinterpret_cast<XtPointer>(this) ); 
  char cmdline[LENGTH_COMMENT];
  cmdline[0]='\0';
  strcat(cmdline,argv[0]);
  for (int i=1;i<argc;i++) {
    strcat(cmdline," ");
    strcat(cmdline,argv[i]);
  }
  XmTextFieldSetString(textfield,cmdline);
  return bboard;
}
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::createMainWindow2(Widget bboard) {
  XmString xm_string
    =XmStringCreateSimple(const_cast<char*>("Start Run Now"));
  Widget button=XtVaCreateManagedWidget("button",
                                        xmPushButtonWidgetClass,
                                        bboard, 
                                        XmNlabelType, XmSTRING, 
                                        XmNlabelString,xm_string,
                                        XmNx, 130, 
                                        XmNy, 125, 
                                        0 ); 
  XmStringFree(xm_string);
  XtAddCallback(button,
                XmNactivateCallback,
                &GUI::runMainCallback,
                reinterpret_cast<XtPointer>(this) ); 
  XtManageChild(bboard);
  XtManageChild(_program_window);
  XtRealizeWidget(_program_toplevel);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::createMainWindow(int argc, char *argv[]) {
  Widget bboard=createMainWindow1(argc,argv);
  createMainWindow2(bboard);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::createFileMenu() {
  _file_menu=createPulldownMenu(_menu_bar,"Quit");
  createPushbutton(_file_menu,"Quit",&GUI::quitCallback,
                   reinterpret_cast<XtPointer>(this));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//callback for view_menu PushButton
void GUI::readMainMenuCallback(Widget,XtPointer client_data,
XtPointer) {
  GUI *gui=reinterpret_cast<GUI*>(client_data);
  gui->main_list->reloadListWin();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::createViewMenu() {
  view_menu=createPulldownMenu(_menu_bar,"View");
  createPushbutton(view_menu,"Main",&GUI::readMainMenuCallback,this);
  char str[7]="Verify";
  main_list->createListWin(_program_toplevel,
    &GUI::checkMainInputCallback,
    reinterpret_cast<XtPointer>(this),str);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::runMainCallback (Widget,XtPointer client_data,XtPointer){
  GUI *gui=reinterpret_cast<GUI*>(client_data);
  GUIInputParameterList::unsetBadparams();
  gui->checkMainInput();
  if (GUIInputParameterList::getBadparams()) {
    gui->showWarningDialog("Parameter Error: Run Cancelled");
    return;
  }
  if (gui->use_thread) {
    gui->showBusyDialog("Running...");
//  for reruns, should have showBusyDialog("Running...",OK,"Stop");
    XtManageChild(XmMessageBoxGetChild(gui->_busy,
                                       XmDIALOG_OK_BUTTON));
    XmString xcl=XmStringCreate(const_cast<char*>("Stop"),
      XmFONTLIST_DEFAULT_TAG);
    XtVaSetValues(gui->_busy, XmNokLabelString, xcl, 0);
    XmStringFree(xcl);
    XtUnmanageChild(XmMessageBoxGetChild(gui->_busy,
                                         XmDIALOG_CANCEL_BUTTON));
  }
  gui->runMainProc();
  XtAppAddWorkProc(gui->program_app,
                   reinterpret_cast<XtWorkProc>(GUI::workCallback),
                   reinterpret_cast<XtPointer>(gui));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bool GUI::workCallback(XtPointer client_data) {
  GUI *gui=reinterpret_cast<GUI*>(client_data);
  if (!gui) return TRUE;
  bool run_finished=!gui->mainRunning() && gui->run_main_called;
  if (!gui->use_thread) {
    if (run_finished) gui->cleanup();
    return run_finished;
  }
  if (run_finished) {
    XtManageChild(XmMessageBoxGetChild(gui->_busy,
                                       XmDIALOG_CANCEL_BUTTON));
    XmString xcl=XmStringCreate(const_cast<char*>("Cleanup"),
      XmFONTLIST_DEFAULT_TAG);
    XtVaSetValues(gui->_busy, XmNcancelLabelString, xcl, 0);
    XmStringFree(xcl);
    XtUnmanageChild(XmMessageBoxGetChild(gui->_busy,
                                         XmDIALOG_OK_BUTTON));
    gui->hideBusyDialog();
    gui->showBusyDialog("Run finished.");
    return TRUE;
  }
  return FALSE;
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//called from runMainCallback
void GUI::showBusyDialog(const char* str) {
  if (!_program_window) return;
  if (!use_thread) return;
        
  XmString xstr=
    XmStringCreateLtoR(const_cast<char*>(str),XmSTRING_DEFAULT_CHARSET);
  if (_busy) {
#ifndef NO_THREAD
    CHECK_POINTER(run_main_thread)
    if (run_main_thread->isAlive()) {
      XtManageChild(XmMessageBoxGetChild(_busy,XmDIALOG_OK_BUTTON));
      XtUnmanageChild(XmMessageBoxGetChild(_busy,
                                           XmDIALOG_CANCEL_BUTTON));
    }
#endif
    XtVaSetValues(_busy,XmNmessageString, xstr, 0);
    XMapRaised(XtDisplay(_busy),XtWindow(_busy));
  } else {
    Arg args[10];
    int n=0;
    Widget  secondary_shell = XtVaCreatePopupShell("secondary_shell",
      topLevelShellWidgetClass, _program_toplevel, 
      XtNinput, False, 
      XmNallowResize, True,
      XmNmappedWhenManaged, False,
      XmNdeleteResponse, XmDESTROY,
      0);

    XmString xok=XmStringCreate(const_cast<char*>("Stop"),
      XmFONTLIST_DEFAULT_TAG);
    XmString xcl=XmStringCreate(const_cast<char*>("Kill"),
      XmFONTLIST_DEFAULT_TAG);
    XtSetArg(args[n], XmNmessageString, xstr); n++;
    XtSetArg(args[n], XmNokLabelString, xok); n++;
    XtSetArg(args[n], XmNcancelLabelString, xcl); n++;
    XtSetArg(args[n], XmNautoUnmanage, False); n++;
    XtSetArg(args[n], XmNnoResize, True); n++;
    XtSetArg(args[n], XmNdeleteResponse, XmDO_NOTHING); n++;
    _busy=XmCreateWorkingDialog(secondary_shell,const_cast<char*>("busy"),
      args,n);
    XmStringFree(xok);XmStringFree(xcl);
    
    XtUnmanageChild(XmMessageBoxGetChild(_busy,XmDIALOG_HELP_BUTTON));
    XtUnmanageChild(XmMessageBoxGetChild(_busy,XmDIALOG_CANCEL_BUTTON));

    XtAddCallback(_busy,XmNokCallback,okBusyCallback,(XtPointer)this);
    XtAddCallback(_busy,XmNcancelCallback,cancelBusyCallback,
                  reinterpret_cast<XtPointer>(this));
  }
  XmStringFree(xstr);
  XtManageChild(_busy);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::okBusyCallback(Widget,XtPointer client_data,XtPointer) {
  GUI *gui=reinterpret_cast<GUI*>(client_data);
  if (!gui) return;
  gui->okBusy();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::cancelBusyCallback(Widget,XtPointer client_data,
XtPointer) {
  GUI *gui=reinterpret_cast<GUI*>(client_data);
  if (!gui) return;
  gui->cancelBusy();
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::showWarningDialog(const char* str) {
  if (!_program_window) return;
  XmString xstr=XmStringCreateLtoR(const_cast<char*>(str),
    XmSTRING_DEFAULT_CHARSET);
  if (_warn) {
    XtVaSetValues(_warn,XmNmessageString, xstr, 0);
    XtManageChild(_warn); 
  } else {
    Arg args[10];
    int n=0;
    XtSetArg(args[n], XmNmessageString, xstr); n++;
    _warn=XmCreateWarningDialog(_program_window,const_cast<char*>("busy"),
      args,n);
  }
  XmStringFree(xstr);
  XtManageChild(_warn);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Widget GUI::createPushbutton(Widget parent,const char *name,
XtCallbackProc callback,XtPointer client_data) {
  Arg args[1];
  Cardinal n=0;
  Widget push=XmCreatePushButton(parent,const_cast<char*>(name),args,n);
  XtAddCallback(push,XmNactivateCallback,callback,client_data);
  XtManageChild(push);
  return push;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::createSeparator(Widget parent) {
  Arg args[1];
  Cardinal n=0;
  XtManageChild(XmCreateSeparator(parent,const_cast<char*>("seq"),args,n));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::createLabelbutton(Widget parent,const char* name) {
  XmString xstr=XmStringCreateSimple(const_cast<char*>(name));
  XtVaCreateManagedWidget("label",
                            xmLabelWidgetClass,parent,
                            XmNalignment,XmALIGNMENT_BEGINNING,
                            XmNlabelType,XmSTRING,
                            XmNlabelString,xstr,
                            0);
  XmStringFree(xstr);
  XmCreateSeparator(parent,const_cast<char*>("seq"),0,0);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Widget GUI::createTogglebutton(Widget parent,const char *name,
XtCallbackProc callback,XtPointer client_data,bool set) {
  XmString st1;
  bool radio=FALSE;
  Widget tog = XtVaCreateManagedWidget(radio?"Radio":name, 
                          xmToggleButtonWidgetClass, parent,
          XmNlabelString, st1=
            XmStringCreateSimple(const_cast<char*>(name)),
          XmNindicatorOn,True, XmNvisibleWhenOff,True,
          XmNset, set, 
          0);
  XmStringFree(st1);
  XtAddCallback(tog,XmNvalueChangedCallback,callback,client_data);
  return tog;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Widget GUI::createPulldownMenu(Widget parent,const char *name,
bool radio) {
  Arg args[5];
  Cardinal n=0;
#ifdef XmNtearOffModel
  XtSetArg(args[n],XmNtearOffModel,XmTEAR_OFF_ENABLED);n++;
#endif
  if (radio) {
    XtSetArg(args[n],XmNradioBehavior,True);n++;
    XtSetArg(args[n],XmNradioAlwaysOne,True);n++;
  }
  Widget menu=XmCreatePulldownMenu(parent,const_cast<char*>(name),args,n);
  
  n=0;
  XtSetArg(args[n],XmNsubMenuId,menu); n++;
  Widget cascade=
    XmCreateCascadeButton(parent,const_cast<char*>(name),args,n);
  XtManageChild(cascade);
  return menu;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Widget GUI::createOptionMenu(Widget parent,char *name) {
  Arg args[5];
  Cardinal n=0;
#ifdef XmNtearOffModel
  XtSetArg(args[n],XmNtearOffModel,XmTEAR_OFF_ENABLED);n++;
#endif
  Widget menu=XmCreateOptionMenu(parent,name,args,n);
  
  Widget pane=XmCreatePulldownMenu(parent,const_cast<char*>("pane"),0,0);
  XtVaSetValues ( menu, XmNsubMenuId, pane, 0 ); 
  XtManageChild ( menu );

  return pane;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::okBusy() { // called by the stop/resume button
#ifndef NO_THREAD
  XmString xok;
  CHECK_POINTER(run_main_thread)
  if (!mainRunning()) cancelBusy();
  else if (run_main_thread->isSuspended()) {
    xok=XmStringCreate(const_cast<char*>("Stop"),XmFONTLIST_DEFAULT_TAG);
    XtVaSetValues(_busy, XmNokLabelString, xok, 0);
    resumeRun();
  } else {
    xok=XmStringCreate(const_cast<char*>("Resume"),XmFONTLIST_DEFAULT_TAG);
    XtVaSetValues(_busy, XmNokLabelString, xok, 0);
    suspendRun();
  }
  XmStringFree(xok);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::cancelBusy() { // called by the cleanup button
#ifndef NO_THREAD
  CHECK_POINTER(run_main_thread)
  if (run_main_thread) {
    ASSERT(run_main_thread->finished());
    run_main_thread->resetStatus();
  }

  XmString xok=
    XmStringCreate(const_cast<char*>("Cleanup."),XmFONTLIST_DEFAULT_TAG);
  XtVaSetValues(_busy, XmNokLabelString,xok,0);
  XmStringFree(xok);
  XtUnmanageChild(_busy);
  cleanup();
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GUI::printOn(ostream &os) const {
  os << "GUI: main_list = " << endl;
  os << *main_list;
  GUIBase::printOn(os);
}
#endif
