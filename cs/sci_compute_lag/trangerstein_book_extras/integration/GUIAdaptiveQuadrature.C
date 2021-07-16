#include <cmath>
#include <iostream>
#include <math.h> // for HUGE_VAL,M_PI
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

enum SCHEME{RIEMANN,MIDPOINT,TRAPEZOIDAL};
SCHEME scheme=MIDPOINT;
int ischeme=static_cast<int>(scheme);
int test_problem=0;

int order=2;
int two_to_order=4;
double tolerance=sqrt(DBL_EPSILON);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  double F77NAME(midpoint)(const double &a,const double &b,const int &n,
    double (*f)(const double &x));
  double F77NAME(riemann)(const double &a,const double &b,const int &n,
    double (*f)(const double &x));
  double F77NAME(trapezoidal)(const double &a,const double &b,
    const int &n,double (*f)(const double &x));

  double f(const double &x) {
    switch (test_problem) {
      case 1:
        return exp(2.*x)*sin(3.*x); // a=1,b=3
        break;
      case 2:
        return exp(3.*x)*sin(2.*x); // a=1,b=3
        break;
      case 3:
        return 2.*x*cos(2.*x)-__cmath_power(x-2.,2); // a=0,b=5
        break;
      case 4:
        return 4.*x*cos(2.*x)-__cmath_power(x-2.,2); // a=0,b=5
        break;
      case 5:
        return sin(1./x); // a=0.1,b=2
        break;
      case 6:
        return cos(1./x); // a=0.1,b=2
        break;
      case 7:
        return sqrt(x); // a=0,b=1
        break;
      case 8:
        return sqrt(1.-x); // a=0,b=1
        break;
      case 9:
        return (x<1. ? exp(0.25*log(1.-x)) : 0.); // a=0,b=1
        break;
      case 0:
      default:
        return exp(x);
        break;
    }
  }
  double absf(const double &x) { return abs(f(x)); }
}
double lower_bound() {
  switch (test_problem) {
    case 1:
    case 2:
      return 1.;
      break;
    case 5:
    case 6:
      return 0.1;
      break;
    case 0:
    case 3:
    case 4:
    case 7:
    case 8:
    case 9:
    default:
      return 0.;
      break;
  }
}
double upper_bound() {
  switch (test_problem) {
    case 1:
    case 2:
      return 3.;
      break;
    case 3:
    case 4:
      return 5.;
      break;
    case 5:
    case 6:
      return 2.;
      break;
    case 0:
    case 7:
    case 8:
    case 9:
    default:
      return 1.;
      break;
  }
}
double (*quadrature_rule)(const double &a,const double &b,
  const int &n,double (*f)(const double &x)) = 0;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double adapt(const double &a,const double &b,
double (*f)(const double &x),double abs_integral_times_tolerance,
XYGraphTool &gt) {
  double coarse_quadrature=quadrature_rule(a,b,1,f);
  double midpoint=0.5*(b+a);
  double fine_quadrature_left=quadrature_rule(a,midpoint,1,f);
  double fine_quadrature_right=quadrature_rule(midpoint,b,1,f);
  double answer=fine_quadrature_left+fine_quadrature_right;
  double error=abs(answer-coarse_quadrature)
              /static_cast<double>(two_to_order-1);
  gt.movePen(a,0.);
  gt.drawLine(a,f(a));
  gt.movePen(midpoint,0.);
  gt.drawLine(midpoint,f(midpoint));
  gt.movePen(b,0.);
  gt.drawLine(b,f(b));
  if (error<=abs_integral_times_tolerance) return answer;
  else {
    double bound=0.5*abs_integral_times_tolerance;
    return adapt(a,midpoint,f,bound,gt)+adapt(midpoint,b,f,bound,gt);
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

  { const char *group="Numerical Method Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ischeme,
      "scheme","riemann","midpoint","trapezoidal") );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(test_problem,
      "test_problem",0,9,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(tolerance,
      "tolerance",DBL_MIN,1.,group) );
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
  switch (scheme) {
    case RIEMANN:
      quadrature_rule=F77NAME(riemann);
      order=1;
      break;
    case TRAPEZOIDAL:
      quadrature_rule=F77NAME(trapezoidal);
      order=2;
      break;
    case MIDPOINT:
    default:
      quadrature_rule=F77NAME(midpoint);
      order=2;
      break;
  }
  two_to_order=__cmath_power(2,order);;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
  double a=lower_bound();
  double b=upper_bound();
  double dx=(b-a)*.0009765625;
  double ymin=HUGE_VAL;
  double ymax=-HUGE_VAL;
  double x=a;
  for (int i=0;i<=1024;i++,x+=dx) {
    double y=f(x);
    ymin=min(ymin,y);
    ymax=max(ymax,y);
  }

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("f vs x","x","f",a,b,ymin,ymax,&cmap,0,winsize);
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.setfgColor("blue");
  x=a;
  for (int i=0;i<=1024;i++,x+=dx) {
    double y=f(x);
    if (i==0) gt.movePen(x,y);
    else gt.drawLine(x,y);
  }
  gt.setfgColor("green");
  double abs_integral_times_tolerance=
    quadrature_rule(a,b,16,absf)*tolerance;
  adapt(a,b,f,abs_integral_times_tolerance,gt);
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
template double std::__cmath_power<double>(double,unsigned);
#include "InputParameter.C"
template class InputParameter<bool>;
