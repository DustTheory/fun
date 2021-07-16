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
#ifndef GTK_GUI_VIRTUAL_INPUT_H
#define GTK_GUI_VIRTUAL_INPUT_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GTKGUIVirtualInput.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include <gtk/gtkcontainer.h>
#include <gtk/gtkentry.h>
#include "VirtualInput.H"

class GTKGUIVirtualInput {
  private:
    char *group_name;

    GtkWidget *group_arrow;
    GtkWidget *group_label;
    GtkWidget *textfield; // gtk_entry_new
    bool value_changed; // indicate if the textfield gets new value
  protected:
    GtkWidget *hbox;

    virtual void valueChanged(GtkWidget*) { value_changed=TRUE; }
    virtual void writeTo(GtkWidget*) const ;
  public:
    GTKGUIVirtualInput(const char *g="");
    virtual ~GTKGUIVirtualInput();

    virtual const char* getValue() const =0;
    virtual char* getName() const =0;
    virtual bool isNamed(const char*) const =0;
    virtual ifstream& formattedRead(ifstream &is) =0;
    virtual ofstream& unformattedWrite(ofstream &os) =0;
    virtual ifstream& unformattedRead(ifstream &is) =0;

    virtual void printOn(ostream & =cerr) const ;

    char* getGroupName() const {return group_name;}
    static void onGTKGUIVirtualInputChangedEvent(GtkWidget*,gpointer);
    static void onGTKGUIVirtualInputFocusOutEvent(GtkWidget*,
      GdkEventFocus*,gpointer);
    bool valueChanged() {return value_changed;}
    void unsetValueChanged() {value_changed=FALSE;}
    bool groupLabelWidgetIs(GtkWidget* w) { return group_label==w; }
    bool groupNameSameAs(const GTKGUIVirtualInput *p) { 
      return strcmp(group_name,p->group_name)==0; 
    }
    GtkWidget* getParent() { return gtk_widget_get_parent(hbox); }
    GtkWidget* getGroupArrow() { return group_arrow; }
    void setGroup(GtkWidget *a,GtkWidget *l) {
      group_arrow=a; 
      group_label=l;
    }
    
//  textfield or radiobox,unmanaged
    virtual GtkWidget* createWidget(GtkWidget*);
    virtual void reload(); //sync textfield with actual data
};
ostream& operator<<(ostream&,const GTKGUIVirtualInput&);

class GTKGUIInputParameterList : public InputParameterList {
  private:
    GtkWidget *widget,*warn;
    bool badparams;
    bool warn_up; // set when showWarningDialog() is first called.
    void (*createListWinCallback)(GtkWidget*,gpointer);
    gpointer createListWinCallbackData;
  public:
    GTKGUIInputParameterList(const char*);
    virtual ~GTKGUIInputParameterList() {;}
    
    void createListWin(void (*)(GtkWidget*,gpointer)=0,gpointer=0,
      const char* =0);
    void showWarningDialog(const char*);
    
    bool getValueChanged();
    void unsetValueChanged();
    bool getBadparams() {return badparams;}
    void unsetBadparams() {badparams=FALSE;}

    static void onGTKGUIInputParameterListDoneClicked(GtkWidget*,
      gpointer);
    void reload();
    static void arrowActCallback(GtkWidget*,gpointer);

    virtual void printOn(ostream&) const ;
};
#endif
