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
// "$Header: /home/faculty/johnt/cvs/deal_new/gui/ivpGUI.C,v 1.1 2009/08/20 17:32:37 johnt Exp $"
#include <float.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include "GUIInputParameter.H"
#include "Debug.H"
#include "GUIInputs.H"
#include "Palette.H"
#include "Tracer.H"
#include "VirtualInput.H"
#include "XYGraphTool.H"

extern double tmax;
extern double x_init;
extern int max_steps;
extern Palette *pal;
extern XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE *cmap;
extern XYGraphTool *gt;

extern GUI_INPUT_PARAMETER_LIST_TYPE *main_list;
extern bool skip_gui;
extern char *display_name;
extern char *sim_name;
extern double winsize;
extern "C" {
  void F77NAME(integrate)(const int &max_steps,const double &tmax,
                  double *soln,double *time);
}
struct odepar_common {
  double rate;
};
extern odepar_common F77NAME(odepar);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//make list of input parameters to be read by GUI
void makeMainList(GUI_INPUT_PARAMETER_LIST_TYPE *&main_list) {
  odepar_.rate=1.;
  main_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main");
  { const char *group="Problem Parameters";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(odepar_.rate,
      "rate",-DBL_MAX,DBL_MAX,group));
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(x_init,
      "initial_x",-DBL_MAX,DBL_MAX,group));
  }
  { const char *group="Timestep Control";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(tmax,"tmax",
      ZERO,DBL_MAX,group));
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(max_steps,
      "max_steps",1,INT_MAX,group));
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
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//this is called from the GUI because we passed the location of this
//  function to the GUI constructor
//use this to check for inconsistencies between different input params
void checkMainInput() {
  if (0.>odepar_.rate && -odepar_.rate*tmax>=2.*double(max_steps)) {
    main_list->showWarningDialog("checkMainInput:\nunstable timestep");
    skip_gui=FALSE;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//this is called from the GUI because we passed the location of this
//  function to the GUI constructor
void runMain(bool /*called_before*/) {
//open graphics output window
  pal=OPERATOR_NEW Palette();
  cmap=OPERATOR_NEW XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE(pal);
  const char *xlabel="t";
  const char *ylabel="x";
  gt=OPERATOR_NEW XYGraphTool(sim_name,xlabel,ylabel,0.,1.,0.,1.,cmap,
    display_name,winsize);

  double *t_array=OPERATOR_NEW_BRACKET(double,max_steps+1);
  double *x_array=OPERATOR_NEW_BRACKET(double,max_steps+1);

#ifdef INDEF
//initialize arrays to help IEEE exception handling catch unassigned
//for double's, can set = HUGE_VAL (which = __infinity)
  for (int i=0;i<=max_steps;i++) {
    t_array[i]=x_array[i]=HUGE_VAL;
  }
#endif
  t_array[0]=0.;  x_array[0]=x_init;

//TimedObject integrate_timing("integrate");
  { //TRACER_CALL(tracer,"integrate block");
    //Timer timer(&integrate_timing);
    integrate_(max_steps,tmax, x_array,t_array);
  }
//integrate_timing.printOn(cout);

//find min,max data values
  double tlo=t_array[0];
  double thi=t_array[max_steps];
  double xlo=x_array[0];
  double xhi=x_array[0];
  { for (int i=1;i<=max_steps;i++) {
    double x=x_array[i];
    if (x<xlo) xlo=x;
    if (x>xhi) xhi=x;
  } }

//initialize graphics display
  gt->newPage();
  gt->rescale(tlo,thi,xlo,xhi);
  gt->setbgColor("white");
  gt->setfgColor("black");
  gt->drawAxes();

//draw numerical solution
  gt->setfgColor("blue");
  gt->movePen(t_array[0],x_array[0]);
  { for (int i=1;i<=max_steps;i++) {
    gt->drawLine(t_array[i],x_array[i]);
  } }
  gt->putString((t_array[0]+t_array[max_steps])*0.5,
    0.9*x_array[0]+0.1*x_array[max_steps],"numerical solution");

//draw analytical solution
  gt->setfgColor("red");
  gt->movePen(t_array[0],x_array[0]);
  { for (int i=1;i<=max_steps;i++) {
    gt->drawLine(t_array[i],x_array[0]*exp(odepar_.rate * t_array[i]));
  } }
  gt->putString((t_array[0]+t_array[max_steps])*0.5,x_array[max_steps],
                "analytical solution");

//force X to perform the requests
  gt->flush();
#ifndef NO_THREAD
  if (skip_gui) {
#endif
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
#ifndef NO_THREAD
  }
#endif

  delete [] t_array;
  delete [] x_array;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//use this for things that need to be done at the end of each call to
//  runMain by the GUI
void cleanup() {
//things allocated in runMain:
  if (gt) delete gt; gt=0;
  if (cmap) delete cmap; cmap=0;
  if (pal) delete pal; pal=0;
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
