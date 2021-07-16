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
#ifndef SURFACE_H
#define SURFACE_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/Surface.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#if (SPACEDIM==3)

#include <iostream>
#include <X11/Intrinsic.h>
#include <GL/glx.h>
#include "VolGraphTool.H"

template<class T> class Surface : 
public VolGraphTool::WINDOW_TYPE::PlotObj {
  private:
    VolGraphTool *graph_tool;
    T *dataset;
    double var_min;
    double var_max;
    double iso_value;
    bool draw_surface;
    VolGraphTool::WINDOW_TYPE::SlideBar *slide_bar;
    VolGraphTool::WINDOW_TYPE::ToggleButtonArray 
      *toggle_button_array;

    void setGraphToolIsoFrac() const {
      if (slide_bar!=0) {
        double frac=ZERO;
        if (var_max>var_min) {
          frac=(iso_value-var_min)/(var_max-var_min);
        }
        graph_tool->setIsoSurfaceFrac(frac);
      }
    }
  public:
//  cf DataSet::DataSet
    Surface(T *t,VolGraphTool *g) : graph_tool(g),dataset(t),
    slide_bar(0),toggle_button_array(0),draw_surface(true) {
      var_min=dataset->getVarMin();
      var_max=dataset->getVarMax();

      VolGraphTool::WINDOW_TYPE *glw=
        dynamic_cast<VolGraphTool::WINDOW_TYPE*>(g->getWindow());
      CHECK_POINTER(glw)
      slide_bar=glw->slideBar();
      CHECK_TEST(slide_bar!=0)

      graph_tool->setIsoSurfaceFrac(ZERO);
      iso_value=var_min;
      slide_bar->rescale(var_min,var_max);
      slide_bar->setValue(iso_value);
      slide_bar->setCallback(&Surface::scaleCallback,
        reinterpret_cast<XtPointer>(this));
    }
    Surface(T *t,VolGraphTool *g,int ns) : graph_tool(g),dataset(t),
    slide_bar(0),toggle_button_array(0),draw_surface(true) {
      var_min=dataset->getVarMin();
      var_max=dataset->getVarMax();

      VolGraphTool::WINDOW_TYPE *glw=
        dynamic_cast<VolGraphTool::WINDOW_TYPE*>(g->getWindow());
      CHECK_POINTER(glw)
      if (ns<=1) {
        slide_bar=glw->slideBar();
        CHECK_TEST(slide_bar!=0)
        graph_tool->setIsoSurfaceFrac(ZERO);
        iso_value=var_min;
        slide_bar->rescale(var_min,var_max);
        slide_bar->setValue(iso_value);
        slide_bar->setCallback(&Surface::scaleCallback,this);
      } else {
        toggle_button_array=glw->toggleButtonArray(ns);
        CHECK_TEST(toggle_button_array!=0)
        toggle_button_array->setCallback(&Surface::toggleCallback,
          reinterpret_cast<XtPointer>(this));
      }
    }
    ~Surface() {
#ifndef USE_GTK
      if (slide_bar!=0) {
        slide_bar->removeCallback(
          reinterpret_cast<XtCallbackProc>(&Surface::scaleCallback),
          this);
      }
      if (toggle_button_array!=0) {
        toggle_button_array->removeCallback(this,
          reinterpret_cast<XtCallbackProc>(&Surface::toggleCallback));
      }
#endif
    }
    void setIsoValue(double value) {
      CHECK_POINTER(slide_bar)
      iso_value=min(var_max,max(var_min,value));
      slide_bar->setValue(iso_value);
      setGraphToolIsoFrac();
      VolGraphTool::WINDOW_TYPE *glw=
        dynamic_cast<VolGraphTool::WINDOW_TYPE*>(
        graph_tool->getWindow());
      CHECK_POINTER(glw)
      glw->clearDisplayLists();
      glw->plotObject(true);
      graph_tool->expose();
    }
    VolGraphTool* getGraphTool() const { return graph_tool; }

    virtual bool drawSurface() const { return draw_surface; }
    virtual void plot(bool do_best) const { 
      if (slide_bar!=0) {
        setGraphToolIsoFrac();
        dataset->plot(do_best,iso_value);
      } else if (toggle_button_array!=0) {
        int ns=toggle_button_array->numberSurfaces();
        double dv=(var_max-var_min)/static_cast<double>(ns+1);
        for (int i=0;i<ns;i++) {
          if (toggle_button_array->drawSurface(i)) {
            double val=var_min+dv*static_cast<double>(i+1);
            dataset->plot(do_best,val);
          }
        }
      }
    }
    virtual void plot(bool do_best,const AxisClipPlane *cp) const {
      dataset->plot(do_best,cp);
    }
#ifdef USE_GTK
    void scale(GtkRange *widget) {    
      CHECK_POINTER(slide_bar)
      iso_value=gtk_range_get_value(widget);
      setGraphToolIsoFrac();
      VolGraphTool::WINDOW_TYPE *glw=
        dynamic_cast<VolGraphTool::WINDOW_TYPE*>(
        graph_tool->getWindow());
      CHECK_POINTER(glw)
      glw->clearDisplayLists();
      glw->plotObject(false);
      glw->plotObject(true);
      glw->expose();
    }
    static void scaleCallback(GtkWidget *widget,gpointer user_data) {
      reinterpret_cast<Surface*>(user_data)->
        scale(reinterpret_cast<GtkRange*>(widget));
    }
#else
    void scale(XmScaleCallbackStruct *scale_struct) {    
      CHECK_POINTER(slide_bar)
      iso_value=slide_bar->getValue(scale_struct->value);
      setGraphToolIsoFrac();
      GLWindow *glw=dynamic_cast<GLWindow*>(graph_tool->getWindow());
      CHECK_POINTER(glw)
      glw->clearDisplayLists();
      glw->plotObject(false);
      glw->plotObject(true);
      glw->expose();
    }
    static void scaleCallback(Widget,XtPointer client_data,
    XtPointer scale_struct) {
      reinterpret_cast<Surface*>(client_data)->
        scale(reinterpret_cast<XmScaleCallbackStruct*>(scale_struct));
    }
#endif
#ifdef USE_GTK
    void toggle(GtkToggleButton *widget,int) {
      VolGraphTool::WINDOW_TYPE *glw=
        dynamic_cast<VolGraphTool::WINDOW_TYPE*>(
        graph_tool->getWindow());
      CHECK_POINTER(glw)
      glw->clearDisplayLists();
      glw->plotObject(false);
      glw->plotObject(true);
      glw->expose();
    }
    static void toggleCallback(GtkWidget *widget,gpointer user_data) {
      reinterpret_cast<Surface*>(user_data)->
        toggle(reinterpret_cast<GtkToggleButton*>(widget),0);
    }
#else
    void toggle(int ns,XmToggleButtonCallbackStruct *toggle_struct) {
      CHECK_POINTER(toggle_button_array)
      bool set=toggle_struct->set;
      GLWindow *glw=dynamic_cast<GLWindow*>(graph_tool->getWindow());
      CHECK_POINTER(glw)
      glw->clearDisplayLists();
      glw->plotObject(false);
      glw->plotObject(true);
      glw->expose();
    }
    static void toggleCallback(Widget w,XtPointer client_data,
    XtPointer toggle_struct) {
      ToggleButtonArrayCallbackStructure<Surface> *cd=
        reinterpret_cast< ToggleButtonArrayCallbackStructure<Surface>* >
        (client_data);
      cd->getData()->toggle(cd->getToggleNumber(),
        reinterpret_cast<XmToggleButtonCallbackStruct*>(toggle_struct));
    }
#endif
    virtual void printOn(ostream &os=cout) const {
      os << "Surface: " << this << endl;
      os << "\tdataset = " << dataset << endl;
      dataset->printOn(os);
    }
};
#endif

#endif
