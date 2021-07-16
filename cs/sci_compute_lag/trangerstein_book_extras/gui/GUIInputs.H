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
#ifndef GUI_INPUTS_H
#define GUI_INPUTS_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/GUIInputs.H,v 1.1 2009/08/20 17:32:35 johnt Exp $"

#ifdef USE_GTK
#include "GTKGUIVirtualInput.H"
#define GUI_INPUTS_BASE GTKGUIVirtualInput
#else
#include "GUIVirtualInput.H"
#define GUI_INPUTS_BASE GUIVirtualInput
#endif

class GUIInputString : public GUI_INPUTS_BASE, public InputString {
  private:
    virtual void valueChanged(GUI_WIDGET_POINTER) ;
  public:
    GUIInputString(char *&p,const char *n,int m) : GUI_INPUTS_BASE(),
      InputString(p,n,m) {;}
    GUIInputString(char *&p,const char *n,int m,const char *g) : 
      GUI_INPUTS_BASE(g),InputString(p,n,m) {;}
    ~GUIInputString() {;}

    virtual const char* getValue() const { return InputString::getValue(); }
    virtual char* getName() const { return InputString::getName(); }
    virtual bool isNamed(const char *n) const { 
      return InputString::isNamed(n);
    }
    virtual ifstream& formattedRead(ifstream &is) {
      return InputString::formattedRead(is);
    }
    virtual ofstream& unformattedWrite(ofstream &os) {
      return InputString::unformattedWrite(os);
    }
    virtual ifstream& unformattedRead(ifstream &is) {
      return InputString::unformattedRead(is);
    }
};

class GUIInputIFStream : public GUI_INPUTS_BASE, public InputIFStream {
  private:
    virtual void valueChanged(GUI_WIDGET_POINTER) ;
  public:
    GUIInputIFStream(ifstream &fs,const char *n,int m) : GUI_INPUTS_BASE(),
      InputIFStream(fs,n,m) {;}
    GUIInputIFStream(ifstream &fs,const char *n,int m,const char *g) :
      GUI_INPUTS_BASE(g),InputIFStream(fs,n,m) {;}

    virtual const char* getValue() const { 
      return InputIFStream::getValue(); 
    }
    virtual char* getName() const { return InputIFStream::getName(); }
    virtual bool isNamed(const char *n) const { 
      return InputIFStream::isNamed(n);
    }
    virtual ifstream& formattedRead(ifstream &is) {
      return InputIFStream::formattedRead(is);
    }
    virtual ofstream& unformattedWrite(ofstream &os) {
      return InputIFStream::unformattedWrite(os);
    }
    virtual ifstream& unformattedRead(ifstream &is) {
      return InputIFStream::unformattedRead(is);
    }
};

class GUIInputOFStream : public GUI_INPUTS_BASE, public InputOFStream {
  private:
    virtual void valueChanged(GUI_WIDGET_POINTER) ;
  public:
    GUIInputOFStream(ofstream &fs,const char *n,int m) : GUI_INPUTS_BASE(),
      InputOFStream(fs,n,m) {;}
    GUIInputOFStream(ofstream &fs,const char *n,int m,const char *g) :
      GUI_INPUTS_BASE(g),InputOFStream(fs,n,m) {;}

    virtual const char* getValue() const { 
      return InputOFStream::getValue(); 
    }
    virtual char* getName() const { return InputOFStream::getName(); }
    virtual bool isNamed(const char *n) const { 
      return InputOFStream::isNamed(n);
    }
    virtual ifstream& formattedRead(ifstream &is) {
      return InputOFStream::formattedRead(is);
    }
    virtual ofstream& unformattedWrite(ofstream &os) {
      return InputOFStream::unformattedWrite(os);
    }
    virtual ifstream& unformattedRead(ifstream &is) {
      return InputOFStream::unformattedRead(is);
    }
};
#endif
