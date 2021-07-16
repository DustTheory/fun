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
#ifndef COLOR_EDITOR_H
#define COLOR_EDITOR_H "$Header: /home/faculty/johnt/cvs/deal_new/gui/ColorEditor.H,v 1.1 2009/08/20 17:32:34 johnt Exp $"

#ifndef USE_GTK
#include <X11/Intrinsic.h>
#include "GUI.H"
#include "GUICallback.H"
#include "ISLList.H"
#include "Palette.H"
#include "XColormap.H"

class RealFlowVariable;
class XPaletteGUI : public GUI {
  private:
    Widget colormap_menu;
  public:
    XPaletteGUI(char *n,char *dn,GUIInputParameterList *ipl,
      void (*m)(bool),void (cm)()=0,void (cu)()=0,void (*sd)()=0,
      bool ut=FALSE) : 
      GUI(n,dn,ipl,m,cm,cu,sd,ut), colormap_menu(0) {;}
    Widget getColormapMenu() const { return colormap_menu; }
    void createViewMenu(bool edit_colormaps) {
      GUI::createViewMenu();
      if (edit_colormaps) {
        colormap_menu=createPulldownMenu(view_menu,"Colormaps");
      }
    }
    virtual void printOn(ostream &os=cerr) const {
      os << "XPaletteGUI: colormap_menu = " << colormap_menu << endl;
      GUI::printOn(os);
    }
};

class XColormap;
class XColorEditor : public ISLListNode {
  public:
    typedef XColormap COLOR_MAP_TYPE;
  private:
    XColormap* xcolormap;
//    don't use the display from XColormap class, 
//    To not have any effect on XWindows.C, use separate event queues.
    Display *display; 
    int screen; //screen from XColormap class, has same hardware display
    Colormap cmap; // cmap from XColormap class

    Widget parent;
    static Widget dialog;
    Widget swatch, current_toggle, form; 
    Widget red_slider,green_slider,blue_slider;
    GC gc;
    int selected;

    Pixel allocNamedColor(char*,Pixel);
    Widget makeSlider(char*,Widget,XtCallbackProc);
    bool warnUserNoColor(Widget);

    static void redSliderMovedCallback(Widget,XtPointer,XtPointer);
    static void blueSliderMovedCallback(Widget,XtPointer,XtPointer);
    static void greenSliderMovedCallback(Widget,XtPointer,XtPointer);

    static void resizeCallback(Widget,XtPointer,XtPointer);
    static void exposeCallback(Widget,XtPointer,XtPointer);
    void expose(Widget,XtPointer);
    static void inputCallback(Widget,XtPointer,XtPointer);
    void input(Widget,XtPointer);

    static void selectColorCallback(Widget,XtPointer,XtPointer);

    static void inputEventHandler(Widget,XtPointer,XEvent*,Boolean*);
  public:
    XColorEditor(XColormap*,Widget);
    static void colorEditorCallback(Widget,XtPointer,XtPointer);

    Palette* getPalette() { return xcolormap->getPalette();}
    void popdown() { XtUnmanageChild(form); }

    void printOn(ostream &os=cerr) const;
};

class XColorEditorList : public ISLList<XColorEditor> {
  private:
  protected:
    XColorEditorList(const XColorEditorList&);
  public:
    XColorEditorList() {;}
    virtual ~XColorEditorList() {;}
    virtual void printOn(ostream&) const;
};
#endif
#endif
