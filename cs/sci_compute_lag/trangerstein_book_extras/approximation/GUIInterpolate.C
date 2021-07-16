#include <iostream>
#include <limits>
#include <math.h> // for M_PI
#include <stdlib.h>

//using namespace std;
//#define RUNGE
//#define CHEBYSHEV

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
#include "GTKColormap.H"
#include "GTKGUI.H"
#include "GTKGUIVirtualInput.H"
#include "GTKWindow.H"
#else
#include "GUI.H"
#include "GUIVirtualInput.H"
#endif

extern "C" {
  void addpoint_(int&,const double&,const double&,double*,double*);
  void divdif_(const int&,const double*,const double*,double*);
  double lagrangepoly_(const int&,const double&,const double*,
    const double*);
  double newtonpoly_(const int&,const double*,const double&,
    const double*);
}

GUI_INPUT_PARAMETER_LIST_TYPE *main_list=0;
bool skip_gui=FALSE;
char *display_name=0;
double winsize=0.5;

enum POINTS{EQUIDISTANT,CHEBYSHEV};
POINTS points=EQUIDISTANT;
int ipoints=static_cast<int>(points);
const char *points_name[2]={"equidistant","Chebyshev"};

enum PROBLEM{SINE,RUNGE};
PROBLEM problem=SINE;
int iproblem=static_cast<int>(problem);
const char *problem_name[2]={"sine","Runge"};

int nmax=3;

double fcn(double t) { 
  return (problem==SINE ? sin(t) : 1./(1.+25.*t*t));
}
void setField(char* field,int field_width,int val) {
  field[field_width]='\0';
  int i=field_width-1;
  int n=val;
  for (;i>=0;i--,n/=10) field[i]='0'+n %10;
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
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ipoints,
      "points",points_name[0],points_name[1]) );
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,iproblem,
      "problem",problem_name[0],problem_name[1]) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nmax,
      "max number interpolation points",2,30,group) );
  }

  { const char *group="Graphics";
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(winsize,
      "winsize",0.,1.,group));
//  main_list->append(OPERATOR_NEW GUIInputParameter<int>(
//    nplot,"nplot",2,INT_MAX,group) );
  }
#ifdef DEBUG
//main_list->printOn(cout);
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void checkMainInput() {
//TRACER_CALL(t,"checkMainInput");
  problem=static_cast<PROBLEM>(iproblem);
  points=static_cast<POINTS>(ipoints);
//if (problem==SINE) nmax=min(nmax,30);
//else nmax=min(nmax,8);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
  char index_field[3];
  char filename[LENGTH_NAME];

  double xmin=(problem==SINE ? 0. : -1.);
  double xmax=(problem==SINE ? M_PI : 1.);
  double ymin=0.;
  double ymax=1.;
  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  char title[80];
  snprintf(title,80,"interpolate %s function using %s points",
    problem_name[problem],points_name[points]);  
  XYGraphTool *gt=OPERATOR_NEW XYGraphTool(title,"x","y",xmin,xmax,
    ymin,ymax,&cmap,0,0.5);
  gt->setbgColor("white");

  int nmin=2;
//int nmax=(problem==SINE ? 30 : 8);
//larger values lead to noticeable effects from
//  rounding errors in divided difs

  double *difs=OPERATOR_NEW_BRACKET(double,nmax);
  double *errors=OPERATOR_NEW_BRACKET(double,nmax);
  double *x=OPERATOR_NEW_BRACKET(double,nmax);
  double *y=OPERATOR_NEW_BRACKET(double,nmax);

  double errors_min=DBL_MAX;
  double errors_max=-DBL_MAX;
  double log10=log(10.);
  for (int i=0;i<nmax;i++) errors[i]=numeric_limits<double>::infinity();
  TimedObject tl("Lagrange interpolation");
  TimedObject tn("Newton interpolation");
  for (int n=nmin;n<nmax;n++) {
#ifdef INDEF
    for (int i=0;i<=n;i++) {
      difs[i]=x[i]=y[i]=numeric_limits<double>::infinity();
    }
#endif

    double dx=(xmax-xmin)/static_cast<double>(n);
    double dxc=0.5*M_PI/static_cast<double>(n+1);
    double dxw=0.5*(xmax-xmin);
    for (int j=0;j<=n;j++) {
      x[j]=xmin+(points==EQUIDISTANT 
        ? dx*static_cast<double>(j)
        : dxw*(1.-cos(dxc*static_cast<double>(2*j+1))));
      y[j]=fcn(x[j]);
//    cout << "\tx,y[" << j << "] = " << x[j] << " " << y[j] << endl;
    }
    divdif_(n,x,y,difs);
    dx*=0.01;
    double t=xmin;
    if ( /* points==EQUIDISTANT || */ n==nmax-1) {
      ymin=numeric_limits<double>::infinity();
      ymax=-ymin;
      for (int j=0;j<=100*n;j++) {
        double t=xmin+dx*static_cast<double>(j);
        double z=newtonpoly_(n,difs,t,x);
        ymin=min(ymin,z);
        ymax=max(ymax,z);
      }
      gt->newPage();
      gt->rescale(xmin,xmax,ymin,ymax);
      gt->setfgColor("black");
      gt->drawAxes();
      gt->setfgColor("blue");
      gt->movePen(t,fcn(t));
      for (int j=0;j<=100*n;j++) {
        double t=xmin+dx*static_cast<double>(j);
        double z=fcn(t);
        gt->drawLine(t,z);
      }
    }
    double maxerr=0.;
    /* if (points==EQUIDISTANT) */ gt->setfgColor("red");
    // else gt->setfgColor(n-nmin+1,nmax-nmin);
    t=xmin;
    if (n==nmax-1) gt->movePen(t,newtonpoly_(n,difs,t,x));
    for (int j=0;j<=100*n;j++) {
      double t=xmin+dx*static_cast<double>(j);
      double z=newtonpoly_(n,difs,t,x);
      if (n==nmax-1) gt->drawLine(t,z);
      maxerr=max(maxerr,abs(z-fcn(t)));
//    cout << t << " " << z << endl;
    }
    errors[n]=log(maxerr)/log10;
    errors_min=min(errors_min,errors[n]);
    errors_max=max(errors_max,errors[n]);
    //if (points==EQUIDISTANT) gt->flush();
    /* else */ {
      if (n==nmax-1) gt->flush();
//    setField(index_field,3,n);
//    snprintf(filename,LENGTH_NAME,"step_%s",index_field);
//    gt->writeGIF(filename);
    }
    cout << "errors[" << n << "] = " << maxerr << endl;

    { Timer zn(&tn);
      for (int j=0;j<=100*n;j++) {
        double t=xmin+dx*static_cast<double>(j);
        double z=newtonpoly_(n,difs,t,x);
      }
    }
    { Timer zl(&tl);
      for (int j=0;j<=100*n;j++) {
        double t=xmin+dx*static_cast<double>(j);
        double z=lagrangepoly_(n, t,x,y);
      }
    }
  }
  cout << "Lagrange interpolation took " << tl.totalRunTime() 
       << " seconds" << endl;
  cout << "Newton interpolation took " << tn.totalRunTime() 
       << " seconds" << endl;

  snprintf(title,80,"error in interpolating %s using %s points",
    problem_name[problem], points_name[points]);  
  XYGraphTool *gt2=OPERATOR_NEW XYGraphTool(title,"log_10(number points)",
    "log_10(error)",log(static_cast<double>(nmin))/log10,
    log(static_cast<double>(nmax))/log10,errors_min,errors_max,&cmap,
    0,0.5);
  gt2->setbgColor("white");
  gt2->setfgColor("black");
  gt2->drawAxes();
  gt2->setfgColor("blue");
  double t=log(static_cast<double>(nmin))/log10;
  gt2->movePen(t,errors[nmin]);
  for (int n=nmin+1;n<nmax;n++) {
    t=log(static_cast<double>(n))/log10;
    gt2->drawLine(t,errors[n]);
  }
  gt2->flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

  delete [] y;
  delete [] x;
  delete [] errors;
  delete [] difs;
  delete gt; gt=0;
  delete gt2; gt2=0;
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
template class InputParameter<bool>;
