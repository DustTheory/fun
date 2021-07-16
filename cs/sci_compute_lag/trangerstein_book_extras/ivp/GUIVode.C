#include <cmath>
#include <iostream>
#include <limits>
#include <stdlib.h>

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
//#include "TimedObject.H"
//#include "Tracer.H"
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

enum HOW_STIFF{NOT_STIFF,STIFF};
HOW_STIFF how_stiff=STIFF;
int ihow_stiff=static_cast<int>(how_stiff);
const char *how_stiff_name[2]={"not stiff (ADAMS)","stiff (BDF)"};

enum ITERATIVE_CORRECTION_METHOD{FUNCTIONAL_ITERATION,EXACT_JACOBIAN,
  FINITE_DIFFERENCE_JACOBIAN,DIAGONAL_JACOBIAN,EXACT_BANDED_JACOBIAN,
  FINITE_DIFFERENCE_BANDED_JACOBIAN};
ITERATIVE_CORRECTION_METHOD iterative_correction_method=
  FUNCTIONAL_ITERATION;
int iiterative_correction_method=
  static_cast<int>(iterative_correction_method);
const char *iterative_correction_method_name[6]={
  "functional iteration","exact jacobian","finite difference jacobian",
  "diagonal jacobian","exact banded jacobian",
  "finite difference banded jacobian"};

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
double rtol=sqrt(numeric_limits<double>::epsilon());
int number_coarse_steps=2;
double rate=1.;
double y_init=1.;

extern "C" {
  void F77NAME(dvode)(
    void (*f)(const int &neq,const double &t,const double *y,
              double *ydot,const double *rpar,const int *ipar),
    const int &neq,double *y,double &t,const double &tout,const int &itol,
    const double &rtol,const double &atol,const int &itask,int &istate,
    const int &iopt,double *rwork,const int &lrw,int *iwork,const int &liw,
    void (*jac)(const int &neq,const double &t,const double *y,
                const int &ml,const int &mu,double *pd,const int &nrowpd,
                const double *rpar,const int *ipar),
    const int &mf,const double *rpar,const int *ipar);
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void f(const int &neq,const double &t,const double *y,double *ydot,
const double *rpar,const int *ipar) { 
  ydot[0]=rate*y[0]; 
}
void jac(const int &neq,const double &t,const double *y,const int &ml,
const int &mu,double *pd,const int &nrowpd,const double *rpar,
const int *ipar) { 
  pd[0]=rate; 
}
double exact(const double &t) { return y_init*exp(rate*t); }
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
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ihow_stiff,
      "how stiff",how_stiff_name[0],how_stiff_name[1]));
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
      "relative error bound",numeric_limits<double>::epsilon(),
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
  how_stiff=static_cast<HOW_STIFF>(ihow_stiff);
  iterative_correction_method=
    static_cast<ITERATIVE_CORRECTION_METHOD>(iiterative_correction_method);
  error_flag=static_cast<ERROR_FLAG>(ierror_flag);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
  double emax=-HUGE_VAL;
  double emin=HUGE_VAL;
  double ymax=max(exact(tstart),exact(tstop));
  double ymin=min(exact(tstart),exact(tstop));
  double log10=log(10.);

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("y versus t","t","y",tstart,tstop,ymin,ymax,&cmap,0,
    winsize);
  gt.newPage();
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();

  double *errors=OPERATOR_NEW_BRACKET(double,number_coarse_steps);
  double *times=OPERATOR_NEW_BRACKET(double,number_coarse_steps);

  int neq=1;
  double *y=OPERATOR_NEW_BRACKET(double,1);
  int itol=1;
  double atol=0.;
  int itask=1;
  int istate=1;
  int iopt=0;
    int jsv=1;
    int meth=static_cast<int>(how_stiff)+1;
    int miter=static_cast<int>(iterative_correction_method);
    int maxord=( meth==1 ? 12 : 5 );
    int ml=0;
    int mu=0;
  int mf=jsv*(10*meth+miter);
    int lwm=0;
    switch (miter) {
      case 0:
        lwm=0;
      case 1:
      case 2:
        lwm=(mf>0 ? 2*neq*neq+2 : neq*neq+2);
        break;
      case 3:
        lwm=neq+2;
        break;
      case 4:
      case 5:
        lwm=(mf>0 ? (3*ml+2*mu+2)*neq+2 : (2*ml+mu+1)*neq+2 );
        break;
    }
  int lrw=20+neq*(maxord+1)+3*neq+lwm;
  double *rwork=OPERATOR_NEW_BRACKET(double,lrw);
  int liw=(miter==0 || miter==3 ? 30 : 30+neq);
  int *iwork=OPERATOR_NEW_BRACKET(int,liw);
  double *rpar=OPERATOR_NEW_BRACKET(double,1);
  int *ipar=OPERATOR_NEW_BRACKET(int,1);

  double dt=(tstop-tstart)/static_cast<double>(number_coarse_steps);
  double h0=dt;
  double t=tstart;
  double tout=tstart;
  y[0]=y_init;
  for (int i=0;i<number_coarse_steps-1;i++) {
    tout=tout+dt;
    F77NAME(dvode)(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,
      lrw,iwork,liw,jac,mf,rpar,ipar);
    cout << "\ty(" << t << ") = " << y[0] << endl;
    ASSERT(istate==2);
    times[i]=t;
    errors[i]=log(abs(y[0]-exact(t)))/log10;
    emax=max(emax,errors[i]);
    emin=min(emin,errors[i]);
    if (i==0) gt.movePen(t,y[0]);
    else gt.drawLine(t,y[0]);
  }
  istate=2;
  tout=tstop;
  F77NAME(dvode)(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,
    iwork,liw,jac,mf,rpar,ipar);
  cout << "\ty(" << t << ") = " << y[0] << endl;
  ASSERT(istate==2);
  times[number_coarse_steps-1]=t;
  errors[number_coarse_steps-1]=log(abs(y[0]-exact(t)))/log10;
  emax=max(emax,errors[number_coarse_steps-1]);
  emin=min(emin,errors[number_coarse_steps-1]);
  gt.drawLine(t,y[0]);
  gt.flush();
  {
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }

  XYGraphTool gte("error versus t","t","log_10(error)",tstart,tstop,
    emin,emax,&cmap,0,winsize);
  gte.setbgColor("white");
  gte.setfgColor("black");
  gte.drawAxes();
  gte.movePen(times[0],errors[0]);
  for (int i=1;i<number_coarse_steps;i++) {
    gte.drawLine(times[i],errors[i]);
  }
  gte.flush();
  {
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }

  delete [] ipar; ipar=0;
  delete [] rpar; rpar=0;
  delete [] iwork; iwork=0;
  delete [] rwork; rwork=0;
  delete [] y; y=0;
  delete [] times; times=0;
  delete [] errors; errors=0;
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
