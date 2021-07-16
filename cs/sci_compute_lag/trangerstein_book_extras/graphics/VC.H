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
#ifndef VIRTUAL_COLORMAP_H
#define VIRTUAL_COLORMAP_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/VC.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#include <iostream>
//#include "Errors.H"
#include "Palette.H"

class VirtualColormap {
  private:
    int num_cmap_colors;
    Palette*const palette;
  protected:
    VirtualColormap( const VirtualColormap&) : num_cmap_colors(0),
    palette(0) { 
      OBSOLETE("copy constructor disallowed"); 
    }
  public:
    VirtualColormap(Palette *pal) : num_cmap_colors(0),palette(pal) {
      CHECK_POINTER(pal)
      num_cmap_colors=pal->getMapIndex(pal->getNumNamedColors()-1)+1;
    }
    virtual ~VirtualColormap() {;}
    int getNumberColors() const { return num_cmap_colors; }
    Palette* getPalette() const { return palette; }

    void allocColor() { num_cmap_colors++; }

    virtual void printOn(ostream &os = cerr) const {
      os << "VirtualColormap: " << endl;
      os << "\nnum_cmap_colors = " << num_cmap_colors
         << "\npalette = " << palette << endl;
    }
};
#endif
