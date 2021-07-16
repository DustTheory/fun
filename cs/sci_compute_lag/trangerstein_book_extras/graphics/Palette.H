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
#ifndef PALETTE_H
#define PALETTE_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/Palette.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#include <fstream>

#include "ISLList.H"

class Palette : public ISLListNode {
  private:
    char *name;
    int max_colors,num_named_colors;
    int *map_index;
    char **color_name;

    Palette(const Palette&) { OBSOLETE("copy constructor disallowed"); }
    void defaultColors();
  public:
    Palette();
    Palette(const char*);
    virtual ~Palette();
    bool isOk() const ;
    void insertEntry(int,const char*);
    int findEntry(const char*) const ;

    const char* getName() const { return name; }
    int getMaxColors() const { return max_colors; }
    int getNumNamedColors() const { return num_named_colors; }
    const char* getColorName(int i=0) const { return color_name[i]; }

    int getMapIndex(int i=0) const {
      CHECK_BOUNDS(i,0,max_colors)
      return map_index[i];
    }

    virtual void printOn(ostream& os = cerr) const;
};

class PaletteList : public ISLList<Palette> {
  private:
    PaletteList(const PaletteList&) { 
      OBSOLETE("copy constructor disallowed"); 
    }
  public:
    PaletteList() {;}
    virtual ~PaletteList();
    void makePalette(const char*,ifstream&,Palette*&);
    const Palette* findByName(const char*) const ;
    virtual void printOn(ostream& os = cerr) const;
};

#endif
