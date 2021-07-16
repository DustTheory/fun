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
#ifndef GUI_INPUT_PARAMETER_H
#define GUI_INPUT_PARAMETER_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIInputParameter.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "InputParameter.H"
#ifdef USE_GTK
#include "GTKGUIVirtualInput.H"
#define GUI_INPUT_PARAMETER_BASE GTKGUIVirtualInput
#define GUI_INPUT_PARAMETER_LIST_TYPE GTKGUIInputParameterList
#else
#include "GUIVirtualInput.H"
#define GUI_INPUT_PARAMETER_BASE GUIVirtualInput
#define GUI_INPUT_PARAMETER_LIST_TYPE GUIInputParameterList
#endif

template<class T> class GUIInputParameter : 
public GUI_INPUT_PARAMETER_BASE,public InputParameter<T> {
  protected:
    virtual void valueChanged(GUI_WIDGET_POINTER) ;
  public:
    GUIInputParameter(T &p,const char *n,T l,T u) : 
      GUI_INPUT_PARAMETER_BASE(),InputParameter<T>(p,n,l,u) {;}
    GUIInputParameter(T &p,const char *n,T l,T u,const char *g) : 
      GUI_INPUT_PARAMETER_BASE(g),InputParameter<T>(p,n,l,u,g) {;}
    virtual ~GUIInputParameter() {;}

    virtual char* getName() const { 
      return InputParameter<T>::getName();
    }
    virtual bool isNamed(const char *str) const { 
      return InputParameter<T>::isNamed(str);
    }
    virtual const char* getValue() const { 
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
    virtual GUI_WIDGET_POINTER createWidget(GUI_WIDGET_POINTER parent);
    virtual void reload() { GUI_INPUT_PARAMETER_BASE::reload(); }
    virtual void printOn(ostream &os=cerr) const {
      os << "GUIInputParameter:" << endl;
      GUI_INPUT_PARAMETER_BASE::printOn(os);
      InputParameter<T>::printOn(os);
    }
};

#endif
