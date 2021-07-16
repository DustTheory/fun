#include <cmath>
#include <iostream>
#include <limits>
#include <stdlib.h>

#include <cvode.h>
#include <dense.h>
#include <cvdense.h>

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
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

enum LINEAR_MULTISTEP_METHOD{ADAMS_MULTISTEP_METHOD,BDF_MULTISTEP_METHOD};
LINEAR_MULTISTEP_METHOD lmm=ADAMS_MULTISTEP_METHOD;
int ilmm=static_cast<int>(lmm);
const char *linear_multistep_method_name[2]={"ADAMS_MULTISTEP_METHOD","BDF_MULTISTEP_METHOD"};

enum ITERATIVE_CORRECTION_METHOD{FUNCTIONAL_ITERATION,NEWTON_ITERATION};
ITERATIVE_CORRECTION_METHOD iterative_correction_method=NEWTON_ITERATION;
int iiterative_correction_method=
  static_cast<int>(iterative_correction_method);
const char *iterative_correction_method_name[2]=
  {"functional iteration","Newton"};

enum ABSOLUTE_ERROR_TOLERANCE_TYPE{SCALAR_ABSOLUTE_ERROR_TOLERANCE,
  VECTOR_ABSOLUTE_ERROR_TOLERANCE};
enum TASK_TYPE{NORMAL_TASK,ONE_STEP_TASK};

enum ERROR_FLAG{SINGLE_RELATIVE_SINGLE_ABSOLUTE,
  SINGLE_RELATIVE_VECTOR_ABSOLUTE,
  VECTOR_RELATIVE_SINGLE_ABSOLUTE,
  VECTOR_RELATIVE_VECTOR_ABSOLUTE};
ERROR_FLAG error_flag=SINGLE_RELATIVE_SINGLE_ABSOLUTE;
int ierror_flag=static_cast<int>(error_flag);
const char *error_flag_name[4]={
  "single relative single absolute","single relative vector absolute",
  "vector relative single absolute","vector relative vector absolute"};

double tstart=0.;
double tstop=1.;
double rtol=1.e-4;
double abtol[3]={1.e-15,1.e-15,1.e-15};
int number_coarse_steps=2;
double rate=1.;
double y_init=1.;

/*
typedef void (*RhsFn)(int N,double t,N_Vector y, N_Vector ydot,
  void *f_data);
typedef void (*CVDenseJacFn)(int N,DenseMat J,
  void (*f)(int,double,N_Vector,N_Vector,void*), void *f_data,
  double t,N_Vector y,N_Vector fy,N_Vector ewt,double h,double uround,
  void *jac_data,long int *nfePtr,N_Vector vtemp1,N_Vector vtemp2,
  N_Vector vtemp3);
*/
extern "C" {
/*
//from CVODE/examples/cvdx.c:
  void f(int N,double t,N_Vector y,N_Vector ydot,void *f_data) {
    double y0=y->data[0];
    double y1=y->data[1];
    double y2=y->data[2];
    double yd0=ydot->data[0]=-0.04*y0+1.e4*y1*y2;
    double yd2=ydot->data[2]=3.e7*y1*y1;
    ydot->data[1]=-yd0-yd2;
  }
  void Jac(int N,DenseMat J,RhsFn f,void *f_data,double t,N_Vector y,
  N_Vector fy,N_Vector ewt,double h,double uround,void *jac_data,
  long int *nfePtr,N_Vector vtemp1,N_Vector vtemp2,N_Vector vtemp3) {
    double y0=y->data[0];
    double y1=y->data[1];
    double y2=y->data[2];
    J->data[0][0]=-0.04;
    J->data[1][0]=1.e4*y2;
    J->data[2][0]=1.e4*y1;
    J->data[0][2]=0.;
    J->data[1][2]=6.e7*y1;
    J->data[2][2]=0.;
    J->data[0][1]=0.04;
    J->data[1][1]=-J->data[1][0]-J->data[1][2];
    J->data[2][1]=-J->data[2][0];
  }
*/
//from Hull, Fellen, Enright and Sedgwick problem B3 (nonlinear reaction):
  void f(int N,double t,N_Vector y,N_Vector ydot,void *f_data) {
    double y0=y->data[0];
    double y1=y->data[1];
    double y2=y->data[2];
    ydot->data[0]=-y0;
    double yd2=ydot->data[2]=y1*y1;
    ydot->data[1]=y0-yd2;
  }
  void Jac(int N,DenseMat J,RhsFn f,void *f_data,double t,N_Vector y,
  N_Vector fy,N_Vector ewt,double h,double uround,void *jac_data,
  long int *nfePtr,N_Vector vtemp1,N_Vector vtemp2,N_Vector vtemp3) {
    double y0=y->data[0];
    double y1=y->data[1];
    double y2=y->data[2];
    J->data[0][0]=-1.;
    J->data[1][0]=0.;
    J->data[2][0]=0.;
    J->data[0][2]=0.;
    J->data[1][2]=2.*y1;
    J->data[2][2]=0.;
    J->data[0][1]=1.;
    J->data[1][1]=-J->data[1][2];
    J->data[2][1]=0.;
  }
}
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

  { const char *group="Initial Value Problem Parameters";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(rate,
      "rate",-numeric_limits<double>::max(),numeric_limits<double>::max(),
      group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(tstop,
      "tstop",0.,numeric_limits<double>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(
      y_init,"y_init",-numeric_limits<double>::max(),
      numeric_limits<double>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      number_coarse_steps,"number_coarse_steps",2,INT_MAX,group) );
  }

  { const char *group="Numerical Method Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ilmm,
      "linear_mutistep_method",linear_multistep_method_name[0],
      linear_multistep_method_name[1]));
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iiterative_correction_method,"iterative correction method",
      iterative_correction_method_name[0],
      iterative_correction_method_name[1],
      iterative_correction_method_name[2],
      iterative_correction_method_name[3],
      iterative_correction_method_name[4],
      iterative_correction_method_name[5]));
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ierror_flag,
      "error flag",error_flag_name[0],error_flag_name[1],
      error_flag_name[2],error_flag_name[3]));
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(rtol,
      "rtol",numeric_limits<double>::epsilon(),
      numeric_limits<double>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(abtol[0],
      "atol[0]",numeric_limits<double>::epsilon(),
      numeric_limits<double>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(abtol[1],
      "atol[1]",numeric_limits<double>::epsilon(),
      numeric_limits<double>::max(),group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(abtol[2],
      "atol[2]",numeric_limits<double>::epsilon(),
      numeric_limits<double>::max(),group) );
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
  if (tstop<tstart) tstop=tstart+1.;
  if (abs(y_init)<=0.) y_init=1.;
  lmm=static_cast<LINEAR_MULTISTEP_METHOD>(ilmm);
  iterative_correction_method=
    static_cast<ITERATIVE_CORRECTION_METHOD>(iiterative_correction_method);
  error_flag=static_cast<ERROR_FLAG>(ierror_flag);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
  N_Vector y=N_VNew(3,0);
  y->data[0]=1.;
  y->data[1]=0.;
  y->data[2]=0.;
  N_Vector abstol=N_VNew(3,0);
  abstol->data[0]=abtol[0];
  abstol->data[1]=abtol[1];
  abstol->data[2]=abtol[2];
  int neq=3;
  double t0=0.;
  int itol=static_cast<int>(VECTOR_ABSOLUTE_ERROR_TOLERANCE);
  void *f_data=0;
  FILE *errfp=0;
  boole optIn=0;
//int OPT_SIZE=40;
  long int iopt[OPT_SIZE];
  double ropt[OPT_SIZE];
  void *machEnv=0;
  void *p=CVodeMalloc(neq,f,t0,y,lmm,iterative_correction_method,itol,
    &rtol,reinterpret_cast<void*>(abstol),f_data,errfp,optIn,iopt,ropt,
    machEnv);
  CHECK_POINTER(p);
  CVodeMemRec *cvode_mem=reinterpret_cast<CVodeMemRec*>(p);

  CVDense(cvode_mem,Jac,0);

  double tout=40.;
  double log10=log(10.);

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
//XYGraphTool gt("y versus t","t","y",t0,tout,0.,1.,&cmap,0,winsize);
  XYGraphTool gt("h versus t","t","h",t0,tout,-10.,0.,&cmap,0,winsize);
  gt.newPage();
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();

  double t=t0;
  int itask=static_cast<int>(ONE_STEP_TASK);
  while (t<tout) {
    double tsave=t;
    double ysave[3]={y->data[0],y->data[1],y->data[2]};
    int flag=CVode(cvode_mem,tout,y,&t,itask);
    double logh=log(cvode_mem->cv_h)/log10;
//  cout << "\tt,order,log(h) = " << t << " " << cvode_mem->cv_q << " "
//    << logh << endl;
    gt.setLineWidth(4*cvode_mem->cv_q);
//  gt.setfgColor("red");
//  gt.movePen(tsave,ysave[0]);
//  gt.drawLine(t,y->data[0]);
//  gt.setfgColor("green");
//  gt.movePen(tsave,ysave[1]);
//  gt.drawLine(t,y->data[1]);
//  gt.setfgColor("blue");
//  gt.movePen(tsave,ysave[2]);
//  gt.drawLine(t,y->data[2]);
//  gt.setLineWidth(1);
//  gt.setfgColor("black");
//  gt.movePen(tsave,ysave[0]);
//  gt.drawLine(t,y->data[0]);
//  gt.movePen(tsave,ysave[1]);
//  gt.drawLine(t,y->data[1]);
//  gt.movePen(tsave,ysave[2]);
//  gt.drawLine(t,y->data[2]);

    gt.movePen(tsave,logh);
    gt.drawLine(t,logh);
  }
  gt.flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

  N_VFree(y);
  N_VFree(abstol);
  CVodeFree(cvode_mem);
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
