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
#ifndef VIRTUAL_WINDOW_H
#define VIRTUAL_WINDOW_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/VirtualWindow.H,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#include <stdio.h>
#include <iostream>
#ifdef USE_GTK
#include <gdk/gdkcolor.h>
#endif
#include <X11/Xlib.h>
#include "Const.H"
//#include "Errors.H"

/*
//#include "MemoryDebugger.H"
# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif
*/

class VirtualWindow {
  public:
#ifdef USE_GTK
    typedef GdkColor COLOR_TYPE;
#else 
    typedef XColor COLOR_TYPE;
#endif
  protected:
    VirtualWindow(const VirtualWindow&) { 
      OBSOLETE("copy constructor disallowed");
    }
  public:
    VirtualWindow() {}
    virtual ~VirtualWindow() {}
    virtual VirtualWindow* getWindow() = 0;
    virtual const VirtualWindow* getWindow() const = 0;
    virtual const char* getName() const = 0;
    virtual int getWidth() const = 0;
    virtual int getHeight() const = 0;
    virtual void newPage() = 0 ;
    virtual void setLineWidth(const int*) const=0 ;
    virtual void setfgColor(const char*) const=0;
    virtual void setbgColor(const char*) const=0;
    virtual void flush(bool=true) = 0 ;
    virtual void printOn(ostream& os = cerr) const {
      os << "VirtualWindow:" << getName() << endl;
    }
};

class VirtualWindow12 : public VirtualWindow {
  protected:
    VirtualWindow12(const VirtualWindow12&) { 
      OBSOLETE("copy constructor disallowed");
    }
  public:
    VirtualWindow12() {}
    virtual ~VirtualWindow12() {}
    virtual const char* getDisplayName() const = 0 ;
    virtual void setFont(const char*) const = 0 ;
    virtual bool monochrome() const = 0 ;
    virtual void setfgColor(const int*,const int*) const=0 ;
    virtual void putString(const double*,const double*,const char*,
      const double*) const=0;
    virtual void writeXBorder(const char *str) const=0;
    virtual void writeYBorder(const char *str) const=0;
    virtual void movePen(const double*,const double*) const = 0 ;
    virtual void drawLine(const double*,const double*,const bool* =0)
      const=0;
    virtual bool buttonIsPressed(double&,double&) = 0 ;
    virtual int getMouse(double&,double&) = 0 ;
    virtual void colorPolygon(const double*,const int *,const double*,
      const double*) const=0;
    virtual void colorPolygon(const int*,const double*,const double*)
      const=0;
    virtual void colorRect(const double*,const double*,const double*,
      const double*,const double*) const=0;
    virtual int pickColor(double) const=0 ;
    virtual void setfgColor(const char*) const = 0;
    virtual void setfgColor(const double*) const = 0;
    virtual void setbgColor(const char*) const = 0;
    virtual void setbgColor(const double*) const = 0;
    virtual void writePaletteName(FILE*) = 0 ;
    virtual void createRaster(FILE*) = 0 ;
    virtual void writeRaster(FILE*) = 0 ;
    virtual void closeRaster() = 0 ;
    virtual void openRaster(FILE*) = 0 ;
    virtual void readRaster(FILE*) = 0 ;
    virtual void printOn(ostream& os = cerr) const {
      os << "VirtualWindow12: display_name=" << getDisplayName()
         << endl;
      VirtualWindow::printOn(os);
    }
};
#endif
