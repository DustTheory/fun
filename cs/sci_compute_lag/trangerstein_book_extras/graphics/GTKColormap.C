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
// "$Header: /home/faculty/johnt/cvs/graphics/GTKColormap.C,v 1.1 2005/06/08 11:04:35 johnt Exp $"
#include "GTKColormap.H"
#include <gdk/gdkdisplaymanager.h>
#include <gdk/gdkdrawable.h>
#include <gdk/gdkscreen.h>
#include <stdlib.h>
//#include "MemoryDebugger.H"
# ifndef OPERATOR_NEW
# define OPERATOR_NEW new
# endif
# ifndef OPERATOR_NEW_BRACKET
# define OPERATOR_NEW_BRACKET(T,n) new T[n]
# endif
//#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKColormap::GTKColormap(Palette *pal,const char *dn,
GdkDisplay *dpy) :
VirtualColormap(pal),cmap_colors(0),cmap(0),display(dpy),screen(0) {
  int num_named_colors=pal->getNumNamedColors();

//need background color and two Palette colors
  CHECK_TEST(num_named_colors>=3);
//background color should be map_index[0]
  CHECK_SAME(pal->getMapIndex(0),0);
//first Palette color should be map_index[1]
  CHECK_SAME(pal->getMapIndex(1),1);
#ifdef DEBUG
  for (int i=2;i<num_named_colors;i++) {
    CHECK_TEST(pal->getMapIndex(i-1)<pal->getMapIndex(i));
  }
#endif

  if (!display) {
    if (dn) display=gdk_display_open(dn);
  }
  if (!display) {
    GdkDisplayManager *display_manager=gdk_display_manager_get();
    display=gdk_display_manager_get_default_display(display_manager);
  }
  ASSERT(display);

  screen=gdk_display_get_default_screen(display);
  cmap=gdk_screen_get_default_colormap(screen);
  cmap_colors=OPERATOR_NEW GdkColor[getNumberColors()];
  
  bool success=true;
  for (int i=0;i<num_named_colors;i++) {
    success&=gdk_color_parse(pal->getColorName(i),
                             &cmap_colors[pal->getMapIndex(i)]);
  }
  ASSERT(success);

  for (int i=2;i<num_named_colors;i++) {
    int map_lo=pal->getMapIndex(i-1);
    int map_hi=pal->getMapIndex(i);
    float denom=static_cast<float>(map_hi - map_lo);
    float lo_red=static_cast<float>(cmap_colors[map_lo].red);
    float hi_red=static_cast<float>(cmap_colors[map_hi].red);
    float lo_green=static_cast<float>(cmap_colors[map_lo].green);
    float hi_green=static_cast<float>(cmap_colors[map_hi].green);
    float lo_blue=static_cast<float>(cmap_colors[map_lo].blue);
    float hi_blue=static_cast<float>(cmap_colors[map_hi].blue);
    for (int j=map_lo+1;j<map_hi;j++) {
      float lo_ratio=static_cast<float>(map_hi - j)/denom;
      float hi_ratio=static_cast<float>(j - map_lo)/denom;
      cmap_colors[j].red=
        static_cast<guint16>(lo_red*lo_ratio+hi_red*hi_ratio);
      cmap_colors[j].green=
        static_cast<guint16>(lo_green*lo_ratio+hi_green*hi_ratio);
      cmap_colors[j].blue=
        static_cast<guint16>(lo_blue*lo_ratio+hi_blue*hi_ratio);
    }
  }
  
  success=true;
  gboolean writeable=true;
  gboolean best_match=true;
  for (int i=0;i<getNumberColors();i++) {
    success &=gdk_colormap_alloc_color(cmap,&cmap_colors[i],writeable,
      best_match);
  }
  ASSERT(success);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GTKColormap::~GTKColormap() {
  for (int i=getNumberColors()-1;i>=0;i--) {
    gdk_colormap_free_colors(cmap,&cmap_colors[i],1);
  }
  delete [] cmap_colors;
  cmap_colors=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColormap::setColor(int i,const GdkColor &color) {
  cmap_colors[i].red=color.red;
  cmap_colors[i].green=color.green;
  cmap_colors[i].blue=color.blue;
  gdk_colormap_free_colors(cmap,&cmap_colors[i],1);
  gboolean writeable=true;
  gboolean best_match=true;
  ASSERT(gdk_colormap_alloc_color(cmap,&cmap_colors[i],writeable,best_match));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColormap::installNewVisual(GdkWindow *w,GdkVisual *v) {
  gint depth=gdk_drawable_get_depth(w);
  gboolean allocate=true;
  cmap=gdk_colormap_new(v,allocate);

  gboolean writeable=true;
  gboolean best_match=true;
  for (int i=0;i<getNumberColors();i++) {
    ASSERT(gdk_colormap_alloc_color(cmap,cmap_colors+i,writeable,best_match));
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
GdkColor GTKColormap::allocColor(const char *colorname) {
  int color_index=getPalette()->findEntry(colorname);
  if (color_index<0) {
    getPalette()->insertEntry(getNumberColors(),colorname);
    GdkColor color;
    bool success=gdk_color_parse(colorname,&color);
    ASSERT(success);
    success=gdk_colormap_alloc_color(cmap,&color,true,true);
    ASSERT(success);
    VirtualColormap::allocColor();
    return color;
  } else return cmap_colors[color_index];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void GTKColormap::printOn(ostream &os) const {
  os << "GTKColormap: " << endl;
  os << "\n\tnum_cmap_colors = " << getNumberColors()
     << "\n\tcmap = " << cmap 
     << "\n\tpalette = " << getPalette() << endl;

  for (int i=0;i<getNumberColors();i++) {
    cout << "\tcolor[" << i << "] = " << cmap_colors[i].red
         << " " << cmap_colors[i].green 
         << " " << cmap_colors[i].blue << endl;
  }
  VirtualColormap::printOn(os);
}
