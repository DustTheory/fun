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
#include "Types.H"
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

enum FUNCTION{EXPONENTIAL,SQUARE_ROOT};
FUNCTION function=EXPONENTIAL;
int ifunction=static_cast<int>(function);

int nmax=2;
int order=2;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  double F77NAME(midpoint)(const double &a,const double &b,const int &n,
    double (*f)(const double &x));
  double F77NAME(riemann)(const double &a,const double &b,const int &n,
    double (*f)(const double &x));
  double F77NAME(trapezoidal)(const double &a,const double &b,
    const int &n,double (*f)(const double &x));

  double f(const double &x) {
    switch (function) {
      case SQUARE_ROOT:
        return sqrt(abs(x));
        break;
      case EXPONENTIAL:
      default:
        return exp(x);
        break;
    }
  }
  double quad(const double &a,const double &b) {
    switch (function) {
      case SQUARE_ROOT:
        return 2.*(sqrt(abs(b*b*b))-sqrt(abs(a*a*a)))/3.;
        break;
      case EXPONENTIAL:
      default:
        return exp(b)-exp(a);
        break;
    }
  }
}
double (*quadrature_rule)(const double &a,const double &b,
  const int &n,double (*f)(const double &x)) = 0;
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
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      ifunction,"function","exponential","square root") );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nmax,
      "max power of 2",2,48,group) );
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
  function=static_cast<FUNCTION>(ifunction);
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
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");

  cout << setprecision(16);

  double *log_error=OPERATOR_NEW_BRACKET(double,nmax);
  double *log_h=OPERATOR_NEW_BRACKET(double,nmax);

  double log10=log(10.);
  double exact=quad(0.,1.);
  int n=1;
  for (int i=0;i<nmax;i++,n*=2) {
    double error=exact-quadrature_rule(0.,1.,n,f);
#ifdef DEBUG
    cout << "error[" << n << "] = " << error << endl;
#endif
    error=abs(error);
    if (error<=0.) error=abs(exact)*DBL_EPSILON;
    log_error[i]=log(error)/log10;
    log_h[i]=log(n)/log10;
  }
  double xmin=HUGE_VAL;
  double xmax=-HUGE_VAL;
  double ymin=HUGE_VAL;
  double ymax=-HUGE_VAL;
  for (int i=0;i<nmax;i++) {
    ymin=min(ymin,log_error[i]);
    ymax=max(ymax,log_error[i]);
    xmin=min(xmin,log_h[i]);
    xmax=max(xmax,log_h[i]);
  }

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("log_10(quadrature error) vs log_10(n)","log_10(n)",
    "log_10(error)",xmin,xmax,ymin,ymax,&cmap,0,winsize);
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.setfgColor("blue");
  for (int i=0;i<nmax;i++) {
    if (i==0) gt.movePen(log_h[i],log_error[i]);
    else gt.drawLine(log_h[i],log_error[i]);
  }
  gt.setfgColor("green");
  double dx=(xmax-xmin)/static_cast<double>(2*nmax);
  for (int i=0;i<nmax;i++) {
    gt.drawPlus(log_h[i],log_error[i],dx);
  }
  {
  XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }

  {
    double *table=OPERATOR_NEW_BRACKET(double,nmax);
    int n=1;
    double p=pow(2.,order);
    for (int i=1;i<=nmax;i++,n*=2) {
      table[nmax-i]=quadrature_rule(0.,1.,n,f);
//    cout << "\n\ti,table[" << nmax-i << "] = " << i << " "
//         << table[nmax-i] << " " << exact - table[nmax-i] << endl;
      double pj=p;
      for (int j=nmax+1-i;j<nmax;j++,pj*=p) {
        table[j]=table[j-1]+(table[j-1]-table[j])/(pj-1.);
//      cout << "\ttable[ " << j << "]  = " << table[j] << " "
//           << exact - table[j] << endl;
      }
      double error=abs(exact-table[nmax-1]);
      if (error<=0.) error=abs(exact)*DBL_EPSILON;
      log_error[i-1]=log(error)/log10;
      log_h[i-1]=log(n)/log10;
    }
    delete [] table;
  }

  xmin=HUGE_VAL;
  xmax=-HUGE_VAL;
  ymin=HUGE_VAL;
  ymax=-HUGE_VAL;
  for (int i=0;i<nmax;i++) {
    ymin=min(ymin,log_error[i]);
    ymax=max(ymax,log_error[i]);
    xmin=min(xmin,log_h[i]);
    xmax=max(xmax,log_h[i]);
  }

  XYGraphTool gt2("log_10(extrapolant error) vs log_10(n)","log_10(n)",
    "log_10(error)",xmin,xmax,ymin,ymax,&cmap,0,winsize);
  gt2.setbgColor("white");
  gt2.setfgColor("black");
  gt2.drawAxes();
  gt2.setfgColor("blue");
  for (int i=0;i<nmax;i++) {
    if (i==0) gt2.movePen(log_h[i],log_error[i]);
    else gt2.drawLine(log_h[i],log_error[i]);
  }
  gt2.setfgColor("green");
  dx=(xmax-xmin)/static_cast<double>(2*nmax);
  for (int i=0;i<nmax;i++) {
    gt2.drawPlus(log_h[i],log_error[i],dx);
  }
  {
  XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }

  delete [] log_error;
  delete [] log_h;
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
#include "InputParameter.C"
template class InputParameter<bool>;
