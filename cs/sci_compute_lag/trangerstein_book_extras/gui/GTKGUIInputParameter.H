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
#ifndef GTK_GUI_INPUT_PARAMETER_H
#define GTK_GUI_INPUT_PARAMETER_H "$Header: /home/faculty/johnt/cvs/gui/GTKGUIInputParameter.H,v 1.2 2007/01/10 13:21:58 johnt Exp $"
#include "InputParameter.H"
#include "GTKGUIVirtualInput.H"

template<class T> class GTKGUIInputParameter : 
public GTKGUIVirtualInput,public InputParameter<T> {
  protected:
    virtual void valueChanged(GtkWidget*) ;
  public:
    GTKGUIInputParameter(T &p,const char *n,T l,T u) : 
      GTKGUIVirtualInput(),InputParameter<T>(p,n,l,u) {;}
    GTKGUIInputParameter(T &p,const char *n,T l,T u,const char *g) : 
      GTKGUIVirtualInput(g),InputParameter<T>(p,n,l,u,g) {;}
    virtual ~GTKGUIInputParameter() {;}

    virtual char* getName() const { 
      return InputParameter<T>::getName();
    }
    virtual bool isNamed(const char *str) const { 
      return InputParameter<T>::isNamed(str);
    }
    virtual char* getValue() const { 
      return InputParameter<T>::getValue();
    }
    virtual ifstream& formattedRead(ifstream &fs) {
      return InputParameter<T>::formattedRead(fs);
    }
    virtual ifstream& unformattedRead(ifstream &fs) {
      return InputParameter<T>::unformattedRead(fs);
    }
    virtual ofstream& unformattedWrite(ofstream &fs) {
      return InputParameter<T>::unformattedWrite(fs);
    }
    virtual GtkWidget* createWidget(GtkWidget *parent);
    virtual void reload() { GUIVirtualInput::reload(); }
    virtual void printOn(ostream &os=cerr) const {
      os << "GUIInputParameter:" << endl;
      GUIVirtualInput::printOn(os);
      InputParameter<T>::printOn(os);
    }
};

#endif
