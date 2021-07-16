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
#ifndef GUI_VIRTUAL_INPUT_H
#define GUI_VIRTUAL_INPUT_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIVirtualInput.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
extern "C" {
#include <X11/Intrinsic.h>
}
#include "VirtualInput.H"

class GUIVirtualInput {
  private:
    char *group_name;

    Widget group_arrow;
    Widget group_label;
    Widget label;
    Widget textfield;
  protected:
    bool value_changed; // indicate if the textfield gets new value

    virtual void valueChanged(Widget) { value_changed=TRUE; }
    virtual void writeTo(Widget) const ;
  public:
    GUIVirtualInput(const char *g="");
    virtual ~GUIVirtualInput();

    virtual const char* getValue() const =0;
    virtual char* getName() const =0;
    virtual bool isNamed(const char*) const =0;
    virtual ifstream& formattedRead(ifstream &is) =0;
    virtual ofstream& unformattedWrite(ofstream &os) =0;
    virtual ifstream& unformattedRead(ifstream &is) =0;

    virtual void printOn(ostream & =cerr) const ;

    char* getGroupName() const {return group_name;}
    static void valueChangedCallback(Widget,XtPointer,XtPointer) ;
    bool valueChanged() {return value_changed;}
    void unsetValueChanged() {value_changed=FALSE;}
    bool isGroupArrowWidget(Widget w) { return group_arrow==w; }
    bool isGroupName(const char *g) { 
      return strcmp(group_name,g)==0; 
    }
    bool sameGroupNameAs(const GUIVirtualInput *p) { 
      return strcmp(group_name,p->group_name)==0; 
    }
    Widget getTextField() { return textfield; }
    Widget getGroupLabel() { return group_label; }
    void setGroupArrow(Widget w) { group_arrow=w; }
    void setGroupLabel(Widget w) { group_label=w; }
    void setLabel(Widget w) { label=w; }
    void setTextField(Widget w) { textfield=w; }
    void manage(unsigned char);
    void setTopWidget(Widget);
    
//  textfield or radiobox,unmanaged
    virtual Widget createWidget(Widget parent);
    virtual void reload(); //sync textfield with actual data
};
ostream& operator<<(ostream&,const GUIVirtualInput&);

class GUIInputParameterList : public InputParameterList {
  private:
    Widget form_dialog_widget,warn;
    static bool badparams;
    bool warn_up; // set when showWarningDialog() is first called.
    void (*verify)(Widget,XtPointer,XtPointer);
    XtPointer verify_data;
  public:
    GUIInputParameterList(const char*);
    virtual ~GUIInputParameterList() {;}

    Widget getWidget() { return form_dialog_widget; }
    
    void createListWin(Widget parent,
      void (*)(Widget,XtPointer,XtPointer)=0,XtPointer=0,const char* =0);
    void showWarningDialog(const char*);

    void reloadListWin(); // actually reload all
    
    bool getValueChanged();
    void unsetValueChanged();
    static bool getBadparams() {return badparams;}
    static void unsetBadparams() {badparams=False;}

    static void hideCallback(Widget,XtPointer client_data,XtPointer);
    void hide();
    static void reloadCallback(Widget,XtPointer client_data,XtPointer);
    static void arrowActCallback(Widget,XtPointer,XtPointer);

    virtual void printOn(ostream&) const ;
};
#endif
