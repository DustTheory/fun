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
// "$Header: /home/faculty/johnt/cvs/graphics/GDKColorArray.C,v 1.1 2005/06/08 11:04:34 johnt Exp $"
#include "GDKColorArray.H"
#include "Palette.H"
#include <gdk/gdkdisplaymanager.h>
#include <gdk/gdkscreen.h>

GDKColorArray::GDKColorArray(const Palette &pal) : 
NumPtr<GdkColor>(pal.getMapIndex(pal.getNumNamedColors()-1)+1) {
  int num_named_colors=pal.getNumNamedColors();
  GdkDisplayManager *display_manager=gdk_display_manager_get();
  GdkDisplay *display=
    gdk_display_manager_get_default_display(display_manager);
  GdkScreen *screen=gdk_display_get_default_screen(display);
  GdkColormap *cmap=gdk_screen_get_default_colormap(screen);

  GdkColor actual;
  bool success=gdk_color_parse(pal.getColorName(0),&actual);
  (*this)[0]=actual;

  success&=gdk_color_parse(pal.getColorName(1),&actual);
  (*this)[1]=actual;

  for (int i=2;i<num_named_colors;i++) {
    success&=gdk_color_parse(pal.getColorName(i),&actual);
    (*this)[pal.getMapIndex(i)]=actual;
    GdkColor low=(*this)[pal.getMapIndex(i-1)];
    GdkColor high=(*this)[pal.getMapIndex(i)];
    int n=pal.getMapIndex(i)-pal.getMapIndex(i-1);
    for (int j=pal.getMapIndex(i-1)+1;j<pal.getMapIndex(i);j++) {
      float frac=float(j-pal.getMapIndex(i-1))/float(n);
      actual.red=low.red
        +(unsigned short)(float(high.red-low.red)*frac);
      actual.green=low.green
        +(unsigned short)(float(high.green-low.green)*frac);
      actual.blue=low.blue
        +(unsigned short)(float(high.blue-low.blue)*frac);
      (*this)[j]=actual;
    }
  }
}

#include "NumPtr.C"
INSTANTIATE_NUMPTR(GdkColor)
