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
#ifndef GTKCOLORMAP_H
#define GTKCOLORMAP_H "$Header: /home/faculty/johnt/cvs/graphics/GTKColormap.H,v 1.1 2005/06/08 11:04:35 johnt Exp $"
#include <iostream>
#include <gdk/gdkcolor.h>
#include "VC.H"

class GTKColormap : public VirtualColormap {
  private:
    GdkColor *cmap_colors;
    GdkColormap *cmap;
    GdkDisplay *display;
    GdkScreen *screen;

    GTKColormap( const GTKColormap&);
  public:
    GTKColormap(Palette*,const char* =0,GdkDisplay* =0);
    virtual ~GTKColormap();
    void installNewVisual(GdkWindow*,GdkVisual*);
    GdkColor allocColor(const char*);
    void setColor(int i,const GdkColor&);
    GdkDisplay* getDisplay() { return display; }
    GdkScreen* getScreen() {return screen;}
    GdkColormap* getColormap() { return cmap; }
    const GdkColormap* getColormap() const { return cmap; }
    GdkVisual* getVisual() { return gdk_colormap_get_visual(cmap); }
    const GdkColor& getColor(int i) const { return cmap_colors[i]; }
    GdkColor& getColor(int i) { return cmap_colors[i]; }
    virtual void printOn(ostream& os = cerr) const;
};
#endif
