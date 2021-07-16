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
#ifndef GDKCOLOR_ARRAY_H
#define GDKCOLOR_ARRAY_H "$Header: /home/faculty/johnt/cvs/graphics/GDKColorArray.H,v 1.1 2005/06/08 11:04:34 johnt Exp $"
//#include <X11/Xlib.h>
#include "Const.H"
#include "NumPtr.H"
#include <gdk/gdkcolor.h>

class Palette;
class GDKColorArray : public NumPtr<GdkColor> {
  public:
    GDKColorArray(const Palette&);
    GdkColor getColor(double frac) const {
      if (frac>=ZERO && frac<=ONE) {
        return (*this)[int(frac*double(getNumber()-3))+2];
      }
      return (*this)[0];
    }
    void printOn(ostream &os=cerr) const {
      os << "GDKColorArray:" << endl;
      NumPtr<GdkColor>::printOn(os);
    }
};
#endif
