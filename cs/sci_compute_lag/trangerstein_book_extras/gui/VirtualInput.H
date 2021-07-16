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
#ifndef VIRTUAL_INPUT_H
#define VIRTUAL_INPUT_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/VirtualInput.H,v 1.1 2009/08/20 17:32:36 johnt Exp $"
#include <iostream>
#include <fstream>
#include <string.h>
#include "Errors.H"
#include "NISLList.H"

#ifdef USE_GTK
#include <gtk/gtkwidget.h>
#define GUI_WIDGET_POINTER GtkWidget*
#else
#include <X11/Intrinsic.h>
#define GUI_WIDGET_POINTER Widget
#endif

class VirtualInput {
  private:
    char *name;
  public:
    VirtualInput(const char *n,const char* =0);
    virtual ~VirtualInput();

    char *getName() const { return name; }
    bool isNamed(const char *n) const { return strcmp(name,n)==0; }

    virtual const char* getValue() const { return ""; }
    virtual ifstream& formattedRead(ifstream &is) { return is; }
    virtual ofstream& unformattedWrite(ofstream &os) { return os; }
    virtual ifstream& unformattedRead(ifstream &is) { return is; }
    virtual void printOn(ostream & =cerr) const ;
};
ostream& operator<<(ostream&,const VirtualInput&);

class InputString : public VirtualInput {
  private:
  protected:
    char **ptr_string,*copy;
    size_t length_string,maxlen;
  public:
    InputString(char*&,const char*,int,const char* =0);
    ~InputString();

    virtual const char* getValue() const { return ptr_string[0]; }
    virtual ifstream& formattedRead(ifstream&);
    virtual ofstream& unformattedWrite(ofstream&);
    virtual ifstream& unformattedRead(ifstream&);
    virtual void printOn(ostream & =cerr) const ;
};

class InputIFStream : public VirtualInput {
  private:
    ifstream *stream;
    char *string;
    size_t maxlen;
  protected:
    void changeString(const char*);
  public:
    InputIFStream(ifstream&,const char*,int,const char* =0);
    ~InputIFStream();

    virtual const char *getValue() const { return string; }
    virtual ifstream& formattedRead(ifstream&);
    virtual ofstream& unformattedWrite(ofstream&);
    virtual ifstream& unformattedRead(ifstream&);
    virtual ifstream& reOpen();
    virtual void printOn(ostream & =cerr) const ;
};


class InputOFStream : public VirtualInput {
  private:
    ofstream *stream;
    char *string;
    size_t maxlen;
  protected:
    void changeString(const char*);
  public:
    InputOFStream(ofstream&,const char*,int,const char* =0);
    ~InputOFStream();

    virtual const char *getValue() const { return string; }
    virtual ifstream& formattedRead(ifstream&);
    virtual ofstream& unformattedWrite(ofstream&);
    virtual ifstream& unformattedRead(ifstream&);
    virtual ofstream& reOpen();
    virtual void printOn(ostream & =cerr) const ;
};

class InputParameterList : public NISLList<VirtualInput> {
  private:
    char *name;
  public:
    InputParameterList(const char*);
    virtual ~InputParameterList();

    char* getName() { return name; }
    ifstream& formattedRead(ifstream&,const char *);
    ofstream& unformattedWrite(ofstream&);
    ifstream& unformattedRead(ifstream&);

    virtual void showWarningDialog(const char*) { OBSOLETE(0);}
    virtual void printOn(ostream&) const ;
};
ostream& operator<<(ostream&,const InputParameterList&);
#endif
