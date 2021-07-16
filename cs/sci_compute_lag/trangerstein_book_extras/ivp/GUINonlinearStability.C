#include <cmath>
#include <iostream>
#include <limits> // for numeric_limits
#include <math.h> // for M_PI
#include <rfftw.h>
#include <stdlib.h>

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "XYGraphTool.H"

#ifdef USE_GTK
#include <gtk/gtkmain.h>
#include "GTKGUI.H"
#include "GTKGUIVirtualInput.H"
#include "GTKWindow.H"
#else
#include "GUI.H"
#include "GUIVirtualInput.H"
#endif

GUI_INPUT_PARAMETER_LIST_TYPE *main_list=0;
bool skip_gui=FALSE;
char *display_name=0;
double winsize=0.5;

enum SCHEME{EULER,MODIFIED_MIDPOINT,ADAMS_BASHFORTH2,MODIFIED_EULER};
SCHEME scheme=EULER;
int ischeme=static_cast<int>(scheme);
const char *scheme_name[4]={
  "euler","modified_midpoint","adams-bashforth2","modified euler"};

int burnin=100;
int nplot=10;
int nrandom=10;
double dtmax=1.;
double ymin=0.;
double ymax=1.;
double yold=numeric_limits<double>::infinity();

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double f(double y) { return y*(1.-y); }
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void processCommandLine(int argc,char *argv[],char *display_name,
bool &skip_gui,ifstream &in_file) {
//TRACER_CALL(t,"processCommandLine");
  if (argc<2) skip_gui=FALSE;
  else {
    in_file.open(argv[1],ios::in);
    CHECK_TEST(!in_file.fail());
    if (in_file) {
      int i=2;
      while (i<argc) {
        if (strcmp(argv[i],"-d")==0) {
          display_name=OPERATOR_NEW_BRACKET(char,strlen(argv[++i])+1);
          strcpy(display_name,argv[i]);
        } else {
          cout << " >>>> error...invalid command line option:"
               << argv[i] << endl;
          abort();
        }
        i++;
      }
    } else {
      cerr << "\ncannot open " << argv[1] << " for input" << endl;
      exit(1);
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void readMainInput(ifstream &in_file,InputParameterList *main_list,
bool &skip_gui) {
//TRACER_CALL(t,"readMainInput");
  bool found_main=FALSE;
  char name[LENGTH_NAME], comment[LENGTH_COMMENT];

  in_file.clear(ios::goodbit);
  in_file.seekg(0,ios::beg);
  while ( in_file >> setw(LENGTH_NAME) >> name ) {
    if ( strcmp(name,"Main") == 0 ) {
      in_file.getline( comment, LENGTH_COMMENT);
      found_main=TRUE;
      while ( in_file >> setw(LENGTH_NAME) >> name ) {
#ifdef DEBUG
//      cout << "\tname = " << name << endl;
#endif
        if ( strcmp(name,"end") == 0 ) break;
        else if (strcmp(name,"skip_gui") == 0) {
          int iskip_gui; // type bool not read correctly, so use int
          in_file >> iskip_gui;
          skip_gui= iskip_gui!=0;
        }
        else main_list->formattedRead(in_file,name);
        in_file.getline( comment, LENGTH_COMMENT);
#ifdef DEBUG
//      cout << "\tcomment = " << comment << endl;
#endif
      }
    } else in_file.getline(comment,LENGTH_COMMENT);
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
  if ( !found_main ) skip_gui=FALSE;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void makeMainList(GUI_INPUT_PARAMETER_LIST_TYPE *&main_list) {
//TRACER_CALL(t,"makeMainList");
  main_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Main");

  { const char *group="Numerical Method Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ischeme,
      "scheme",scheme_name[0],scheme_name[1],scheme_name[2],
      scheme_name[3]));
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(burnin,
      "burnin",1,numeric_limits<int>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nplot,
      "nplot",1,numeric_limits<int>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nrandom,
      "nrandom",1,numeric_limits<int>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(dtmax,
      "dtmax",0.,numeric_limits<double>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(ymin,
      "ymin",-numeric_limits<double>::max(),numeric_limits<double>::max(),
      group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(ymax,
      "ymax",-numeric_limits<double>::max(),numeric_limits<double>::max(),
      group) );
  }

  { const char *group="Graphics";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(winsize,
      "winsize",0.,1.,group));
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(t,"checkMainInput");
  scheme=static_cast<SCHEME>(ischeme);
  if (ymin>=ymax) ymax=ymin+1.;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double method(double y,double dt) {
  switch (scheme) {
    case EULER:
      return y+dt*f(y);
    case MODIFIED_MIDPOINT: {
      double ynew=yold+2.*dt*f(y);
      yold=y;
      return ynew;
    }
    case ADAMS_BASHFORTH2: {
      double ynew=y+0.5*dt*(3.*f(y)-f(yold));
      yold=y;
      return ynew;
    }
    case MODIFIED_EULER:
      return y+dt*f(y+0.5*dt*f(y));
    default:
      OBSOLETE("unknown scheme");
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("fixed points","dt","y",0.,dtmax,ymin,ymax,&cmap,
    0,winsize);
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.setfgColor("blue");

  double dy=ymax-ymin;
  double ddt=gt.getPixelWidth();
  for (double dt=ddt;dt<dtmax;dt+=ddt) {
    for (int j=0;j<nrandom;j++) {
      yold=ymin+drand48()*dy;
      double y=ymin+drand48()*dy;
      for (int n=0;n<burnin;n++) {
        if (abs(y)*DBL_EPSILON>1.) break;
        y=method(y,dt);
      }
      for (int n=0;n<nplot;n++) {
        if (abs(y)*DBL_EPSILON>1.) break;
        y=method(y,dt);
        gt.drawPlus(dt,y,ddt);
      }
    }
  }
  gt.flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void cleanup() {
//TRACER_CALL(t,"cleanup");
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//use this for things that need to be done at the end of GUI operations
void shutdown() {
//TRACER_CALL(t,"shutdown");
//things allocated in processCommandLine:
  if (display_name) delete display_name; display_name=0;

//things allocated in makeMainList:
  while (main_list->notEmpty()) {
    delete main_list->delAfter(0);
  }
  delete main_list; main_list=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int main(int argc,char* argv[]) {
//cout << "\tin main" << endl;
#ifdef DEBUG
//setTraps();
#endif
#if MEM_DEBUG
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

  if (skip_gui) {
    runMain(0);
    cleanup();
    shutdown();
  } else {
#ifdef USE_GTK
    GTKGUI gui(argv[0],display_name,main_list,&runMain,
      &checkMainInput,&cleanup,&shutdown,TRUE);
#else
    GUI gui(argv[0],display_name,main_list,&runMain,&checkMainInput,
      &cleanup,&shutdown,TRUE);
#endif
    gui.createFileMenu();
    gui.createViewMenu();
    gui.createHelpMenu();
    gui.createMainWindow(argc,argv);
    gui.eventLoop();
  }
  return EXIT_SUCCESS;
}
template int std::__cmath_power<int>(int,unsigned);
template class InputParameter<bool>;
