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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/areaGraphToolColorGUI.C,v 1.1 2009/08/20 17:32:37 johnt Exp $"
#include <float.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include "AreaGraphTool.H"
#include "GUIInputParameter.H"
#include "GUIInputs.H"
#include "Palette.H"
#include "Tracer.H"
#include "VirtualInput.H"

#ifdef USE_GTK
#include "GTKColorEditor.H"
#else
#include "ColorEditor.H"
#endif

int npts=2;
Palette *pal=0;
AreaGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE *cmap=0;
AreaGraphTool *gt=0;

GUI_INPUT_PARAMETER_LIST_TYPE *main_list=0;
bool skip_gui=FALSE;
char *display_name=0;
char *sim_name=0;
double winsize=0.5;

extern void processCommandLine(int,char *argv[],char*,bool&,
  ifstream&);
extern void readMainInput(ifstream&,InputParameterList*,bool&);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//make list of input parameters to be read by GUI
void makeMainList(GUI_INPUT_PARAMETER_LIST_TYPE *&main_list) {
  main_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main");
  { const char *group="Function";
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(npts,
      "npts",2,INT_MAX,group));
  }
  { const char *group="Graphics";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(winsize,
      "winsize",0.,1.,group));
    main_list->append(OPERATOR_NEW GUIInputString(sim_name,"sim_name",
      LENGTH_NAME,group));
    if (!sim_name) {
      sim_name=OPERATOR_NEW_BRACKET(char,1);
      sim_name[0]='\0';
    }
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//this is called from the GUI because we passed the location of this
//  function to the GUI constructor
//use this to check for inconsistencies between different input params
void checkMainInput() {
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//use this for things that need to be done at the end of each call to
//  runMain by the GUI
void cleanup() {
//things allocated in runMain:
  if (gt) delete gt; gt=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//this is called from the GUI because we passed the location of this
//  function to the GUI constructor
void runMain(bool /*called_before*/) {
//open graphics output window

  double rlo[2]; rlo[0]=0.; rlo[1]=0.;
  double rhi[2]; rhi[0]=1.; rhi[1]=1.;
  double zmin=HUGE_VAL,zmax=-HUGE_VAL;

//  define the surface via a uniform lattice
  double *z=OPERATOR_NEW_BRACKET(double,npts*npts);
  double dx=(rhi[0]-rlo[0])/double(npts);
  double dy=(rhi[1]-rlo[1])/double(npts);
  for (int j=0;j<npts;j++) {
    double y=rlo[1]+dy*double(j);
    for (int i=0;i<npts;i++) {
      double x=rlo[0]+dx*double(i);
      double zz=x*x+y*y;                      // definition of surface
      z[i+j*npts]=zz;
      if (zmin>zz) zmin=zz;
      if (zmax<zz) zmax=zz;
    }
  }
  CHECK_TEST(zmax>zmin);

//open graphics output window
  gt=OPERATOR_NEW AreaGraphTool(sim_name,"x","y",rlo,rhi,cmap,0,
    winsize);
  gt->drawAxes();

//color the image
  double xpos[4],ypos[4];
  for (int j=0;j<npts;j++) {
    ypos[1]=ypos[0]=rlo[1]+dy*double(j);
    ypos[3]=ypos[2]=ypos[1]+dy;
    for (int i=0;i<npts;i++) {
      xpos[3]=xpos[0]=rlo[0]+dx*double(i);
      xpos[2]=xpos[1]=xpos[0]+dx;
      double ratio=(z[i+j*npts]-zmin)/(zmax-zmin);
      ratio=max(0.,min(1.,ratio));
      gt->colorQuad(ratio,xpos,ypos);
    }
  }
  gt->flush();

  if (z) delete [] z;
#ifndef NO_THREAD
  if (skip_gui) {
#endif
    AreaGraphTool::WINDOW_TYPE::QuitButton qb;
    cleanup();
#ifndef NO_THREAD
  }
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//use this for things that need to be done at the end of GUI operations
void shutdown() {
//things allocated in processCommandLine:
  if (display_name) delete display_name; display_name=0;
  if (sim_name) delete sim_name; sim_name=0;

//things allocated in makeMainList:
  while (main_list->notEmpty()) {
    delete main_list->delAfter(0);
  }
  delete main_list; main_list=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int main(int argc,char* argv[]) {
  cout << boolalpha;
  { 
#ifdef MEM_DEBUG
    MemoryDebugger md(1);
#endif
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    ifstream in_file;
    processCommandLine(argc,argv,display_name,skip_gui,in_file);
    makeMainList(main_list);
    if (in_file) {
      readMainInput(in_file,main_list,skip_gui);
      if (skip_gui) checkMainInput();
    }

//  open graphical user interface; turn control of events to it
    pal=OPERATOR_NEW Palette();
    cmap=OPERATOR_NEW AreaGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE(pal);
    if (skip_gui) {
      runMain(0);
      cleanup();
      shutdown();
    } else {
#ifdef USE_GTK
      GTKPaletteGUI gui(argv[0],display_name,main_list,&runMain,
        &checkMainInput,&cleanup,&shutdown,TRUE);
#else
      XPaletteGUI gui(argv[0],display_name,main_list,&runMain,
        &checkMainInput,&cleanup,&shutdown,TRUE);
#endif
      gui.createFileMenu();
      bool edit_colormaps=TRUE;
      gui.createViewMenu(edit_colormaps);
#ifdef USE_GTK
      GTKGUI::PullDownMenu *colormap_menu=gui.getColormapPulldown();
      GTKColorEditor *ce=0;
      if (colormap_menu) {
        ce=OPERATOR_NEW GTKColorEditor(cmap,0);
        colormap_menu->createPushButton("Default Colormap",
          &GTKColorEditor::colorEditorCallback,ce);
      }
#else
      GUI_WIDGET_POINTER colormap_menu=gui.getColormapMenu();
      XColorEditor *ce=0;
      if (colormap_menu) {
        ce=OPERATOR_NEW XColorEditor(cmap,colormap_menu);
        gui.createPushbutton(colormap_menu,"Default Colormap",
          &XColorEditor::colorEditorCallback,ce);
      }
#endif
      gui.createHelpMenu();
      gui.createMainWindow(argc,argv);
      gui.eventLoop();
      delete ce;
    }
    if (cmap) delete cmap; cmap=0;
    if (pal) delete pal; pal=0;
  }
  return EXIT_SUCCESS;
}
