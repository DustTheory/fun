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
#ifndef GTK_COLOR_EDITOR_H
#define GTK_COLOR_EDITOR_H "$Header: /home/faculty/johnt/cvs/gui/GTKColorEditor.H,v 1.1 2005/06/08 11:19:01 johnt Exp $"

#ifdef USE_GTK

#include "GTKColormap.H"
#include "GTKGUI.H"
#include "GUICallback.H"
#include "ISLList.H"
#include "NumPtr.H"
#include "Palette.H"

class GTKPaletteGUI : public GTKGUI {
  private:
    PullDownMenu *colormap_menu;
  public:
    GTKPaletteGUI(char *n,char *dn,GTKGUIInputParameterList *ipl,
      void (*m)(bool),void (cm)()=0,void (cu)()=0,
      void (*sd)()=0,bool ut=FALSE) : 
      GTKGUI(n,dn,ipl,m,cm,cu,sd,ut), colormap_menu(0) {;}
    ~GTKPaletteGUI() { 
      if (colormap_menu) delete colormap_menu; colormap_menu=0;
    }
    PullDownMenu* getColormapPulldown() const { return colormap_menu; }
    GtkWidget* getColormapMenu() const { 
      return (colormap_menu ? colormap_menu->getWidget() : 0);
    }
    virtual void createViewMenu(bool edit_colormaps) {
      GTKGUI::createViewMenu();
      if (edit_colormaps) {
        colormap_menu=OPERATOR_NEW PullDownMenu(view_menu->getWidget(),
          "Colormaps");
      }
    }
    virtual void printOn(ostream &os=cerr) const {
      os << "GTKPaletteGUI: colormap_menu = " << colormap_menu << endl;
      GTKGUI::printOn(os);
    }
};

class GTKColormap;
class GTKColorEditor : public ISLListNode {
  public:
    typedef GTKColormap COLOR_MAP_TYPE;
  private:
    GTKColormap *gtkcolormap;
    int num_cmap_colors;

    GtkWidget *form; 
    GtkWidget *drawing_area;
    GtkWidget *swatch;
    NumPtr<GtkWidget*> toggle;
    GdkGC *gc;
    int selected;
    GTKPaletteGUI* gui;
    static GtkWidget *dialog;

    GdkColor allocNamedColor(char*,const GdkColor&);
    void colorToggle(GtkWidget*,const GdkColor&) const;
    void configureToggle(GtkWidget*,GtkWidget*,guint,guint,guint,int)
      const;
    static void okColorSelectionCallback(GtkWidget*,gpointer);
    static void cancelColorSelectionCallback(GtkWidget*,gpointer);
    bool warnUserNoColor(GtkWidget*) const;

    void drawRGBBars(int,int) const;
    void drawRGBCurves(int,int) const;
    static void realizeCallback(GtkObject*,gpointer);
    static void exposeCallback(GtkObject*,GdkEventExpose*,gpointer);
    static void exposeSwatchCallback(GtkWidget*,GdkEventExpose*,
      gpointer);
    static void buttonPressCallback(GtkWidget*,GdkEventButton*,
      gpointer);

    static void selectColorCallback(GtkWidget*,gpointer);
    static void changeColorCallback(GtkWidget*,gpointer);
  public:
    GTKColorEditor(GTKColormap*,GtkWidget*);
    virtual ~GTKColorEditor() {;}
    static void colorEditorCallback(GtkWidget*,gpointer);

    Palette* getPalette() { return gtkcolormap->getPalette();}
    void popdown() { gtk_widget_unmap(form); }
    void printOn(ostream &os=cerr) const;
};

class GTKColorEditorList : public ISLList<GTKColorEditor> {
  private:
  protected:
    GTKColorEditorList(const GTKColorEditorList&);
  public:
    GTKColorEditorList() {;}
    virtual ~GTKColorEditorList() {;}
    virtual void printOn(ostream&) const;
};
#endif
#endif
