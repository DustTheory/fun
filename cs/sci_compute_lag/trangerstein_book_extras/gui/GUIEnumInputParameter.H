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
#ifndef GUI_ENUM_INPUT_PARAMETER_H
#define GUI_ENUM_INPUT_PARAMETER_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIEnumInputParameter.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"
#include "GUIInputParameter.H"
#include "NumPtr.H"

#ifndef SET_VALUE_NAMES
#define SET_VALUE_NAMES(m) \
  value_names[m]=OPERATOR_NEW_BRACKET(char,strlen(label ## m)+1); \
  strcpy(value_names[m],label ## m);
#endif
class GUIEnumInputParameter : public GUIInputParameter<int> {
  private:
    NumPtr<char*> value_names;
    NumPtr<GUI_WIDGET_POINTER> widgets;
  protected:
    virtual void valueChanged(GUI_WIDGET_POINTER) ;
  public:
    GUIEnumInputParameter(int &p,const char *n,const char *label0,
    const char *label1) : GUIInputParameter<int>(p,n,0,1),
    value_names(2),widgets(2) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
    }
    GUIEnumInputParameter(const char *g,int &p,const char *n,
    const char *label0,const char *label1) : 
    GUIInputParameter<int>(p,n,0,1,g),value_names(2),widgets(2) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
    }
    GUIEnumInputParameter(int &p,const char *n,const char *label0,
    const char *label1,const char *label2) : 
    GUIInputParameter<int>(p,n,0,2),value_names(3),widgets(3) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
    }
    GUIEnumInputParameter(const char *g,int &p,const char *n,
    const char *label0,const char *label1,const char *label2) : 
    GUIInputParameter<int>(p,n,0,2,g),value_names(3),widgets(3) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
    }
    GUIEnumInputParameter(int &p,const char *n,const char *label0,
    const char *label1,const char *label2,const char *label3) : 
    GUIInputParameter<int>(p,n,0,3),value_names(4),widgets(4) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
    }
    GUIEnumInputParameter(const char *g,int &p,const char *n,
    const char *label0,const char *label1,const char *label2,
    const char *label3) : GUIInputParameter<int>(p,n,0,3,g),
    value_names(4),widgets(4) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
    }
    GUIEnumInputParameter(int &p,const char *n,const char *label0,
    const char *label1,const char *label2,const char *label3,
    const char *label4) : GUIInputParameter<int>(p,n,0,4),
    value_names(5),widgets(5) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
    }
    GUIEnumInputParameter(const char *g,int &p,const char *n,
    const char *label0,const char *label1,const char *label2,
    const char *label3,const char *label4) : 
    GUIInputParameter<int>(p,n,0,4,g),value_names(5),widgets(5) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
    }
    GUIEnumInputParameter(int &p,const char *n,const char *label0,
    const char *label1,const char *label2,const char *label3,
    const char *label4,const char *label5) : 
    GUIInputParameter<int>(p,n,0,5),value_names(6),widgets(6) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
      SET_VALUE_NAMES(5);
    }
    GUIEnumInputParameter(const char *g,int &p,const char *n,
    const char *label0,const char *label1,const char *label2,
    const char *label3,const char *label4,const char *label5) : 
    GUIInputParameter<int>(p,n,0,5,g),value_names(6),widgets(6) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
      SET_VALUE_NAMES(5);
    }
    GUIEnumInputParameter(int &p,const char *n,const char *label0,
    const char *label1,const char *label2,const char *label3,
    const char *label4,const char *label5,const char *label6) : 
    GUIInputParameter<int>(p,n,0,6),value_names(7),widgets(7) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
      SET_VALUE_NAMES(5);
      SET_VALUE_NAMES(6);
    }
    GUIEnumInputParameter(const char *g,int &p,const char *n,
    const char *label0,const char *label1,const char *label2,
    const char *label3,const char *label4,const char *label5,
    const char *label6) : GUIInputParameter<int>(p,n,0,6,g),
    value_names(7),widgets(7) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
      SET_VALUE_NAMES(5);
      SET_VALUE_NAMES(6);
    }
    GUIEnumInputParameter(int &p,const char *n,const char *label0,
    const char *label1,const char *label2,const char *label3,
    const char *label4,const char *label5,const char *label6,
    const char *label7) : GUIInputParameter<int>(p,n,0,7),
    value_names(8),widgets(8) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
      SET_VALUE_NAMES(5);
      SET_VALUE_NAMES(6);
      SET_VALUE_NAMES(7);
    }
    GUIEnumInputParameter(const char *g,int &p,const char *n,
    const char *label0,const char *label1,const char *label2,
    const char *label3,const char *label4,const char *label5,
    const char *label6,const char *label7) : 
    GUIInputParameter<int>(p,n,0,7,g),
    value_names(8),widgets(8) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
      SET_VALUE_NAMES(5);
      SET_VALUE_NAMES(6);
      SET_VALUE_NAMES(7);
    }
    GUIEnumInputParameter(int &p,const char *n,const char *label0,
    const char *label1,const char *label2,const char *label3,
    const char *label4,const char *label5,const char *label6,
    const char *label7,const char *label8) : 
    GUIInputParameter<int>(p,n,0,8),value_names(9),widgets(9) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
      SET_VALUE_NAMES(5);
      SET_VALUE_NAMES(6);
      SET_VALUE_NAMES(7);
      SET_VALUE_NAMES(8);
    }
    GUIEnumInputParameter(const char *g,int &p,const char *n,
    const char *label0,const char *label1,const char *label2,
    const char *label3,const char *label4,const char *label5,
    const char *label6,const char *label7,const char *label8) : 
    GUIInputParameter<int>(p,n,0,8,g),value_names(9),widgets(9) {
      SET_VALUE_NAMES(0);
      SET_VALUE_NAMES(1);
      SET_VALUE_NAMES(2);
      SET_VALUE_NAMES(3);
      SET_VALUE_NAMES(4);
      SET_VALUE_NAMES(5);
      SET_VALUE_NAMES(6);
      SET_VALUE_NAMES(7);
      SET_VALUE_NAMES(8);
    }
    virtual ~GUIEnumInputParameter() {
      for (int i=0;i<=upper_bound;i++) delete [] value_names[i];
    }

    virtual GUI_WIDGET_POINTER createWidget(GUI_WIDGET_POINTER parent);
    virtual void printOn(ostream &os=cerr) const {
      os << "GUIInputParameter:" << endl;
      GUIInputParameter<int>::printOn(os);
    }
};

#endif
