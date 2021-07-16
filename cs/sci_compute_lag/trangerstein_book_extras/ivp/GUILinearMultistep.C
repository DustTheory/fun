#include <cmath>
#include <iostream>
#include <limits>
//#include <math.h>
//#include <rfftw.h>
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

enum SCHEME{ADAMS_BASHFORTH,ADAMS_MOULTON,EXPLICIT_MIDPOINT,
  BACKWARD_DIFFERENTIATION_FORMULA};
SCHEME scheme=ADAMS_BASHFORTH;
int ischeme=static_cast<int>(scheme);

int number_refinements=0;
int coarse_number_steps=1;
int ratio=2;
double rate=1.;
double tstop=1.;
double y_init=1.;

struct multistep_common {
  double lm_gamma[10],tolerance;
  int maxit,order,step;
};
extern multistep_common F77NAME(multistep);

extern "C" {
  void F77_NAME(adams_bashforth_startup)();
  void F77_NAME(adams_bashforth)(const double &dt,const double &t,
    double &y,double (*f)(const double &t,const double &y),
    double *ftable);
  void F77_NAME(adams_moulton_startup)();
  void F77_NAME(adams_moulton)(const double &dt,const double &t,
    double &y,double (*dfdy)(const double &t,const double &y),
    double (*f)(const double &t,const double &y),double *ftable,
    double *ftable_old);
  void F77_NAME(bdf_startup)();
  void F77NAME(bdf)(const double &dt,const double &t,
    double &y,double (*dfdy)(const double &t,const double &y),
    double (*f)(const double &t,const double &y),double *ytable,
    double *ytable_old);
  void F77_NAME(explicit_midpoint)(const double &dt,const double &t,
    double &y,double &yold,
    double (*f)(const double &t,const double &y));
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  double f(const double &t,const double &y) { return rate*y; }
  double dfdy(const double &t,const double &y) { return rate; }
  double exact(const double &t) { return y_init*exp(rate*t); }
//double f(const double &t,const double &y) { return t; }
//double dfdy(const double &t,const double &y) { return 0.; }
//double exact(const double &t) { return y_init+t*t*0.5; }
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
  }

  { const char *group="Numerical Method Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ischeme,
      "scheme","adams-bashforth","adams-moulton","explicit midpoint",
      "bdf") );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      F77NAME(multistep).order,"order",1,10,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      coarse_number_steps,"coarse_number_steps",1,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      number_refinements,"number_refinements",0,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(ratio,
      "ratio",2,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      F77NAME(multistep).maxit,"newton iterations",1,20,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(
      F77NAME(multistep).tolerance,"newton tolerance",
      numeric_limits<double>::epsilon(),numeric_limits<double>::max(),
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
  if (coarse_number_steps<=F77NAME(multistep).order
  && scheme!=EXPLICIT_MIDPOINT) {
//  program resets initial values to exact results for mesh refinement 
    coarse_number_steps=F77NAME(multistep).order+1;
  }
  if (scheme==EXPLICIT_MIDPOINT && coarse_number_steps%2 !=0) {
//  explicit midpoint must run an even number of steps
    coarse_number_steps++;
  }
  if (scheme==BACKWARD_DIFFERENTIATION_FORMULA) {
//  bdf unstable for order > 6
    F77NAME(multistep).order=min(6,F77NAME(multistep).order);
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
#ifdef DEBUG
  cout << "\torder = " << F77NAME(multistep).order << endl;
#endif
  double *log_error=OPERATOR_NEW_BRACKET(double,number_refinements+1);
  double *log_n=OPERATOR_NEW_BRACKET(double,number_refinements+1);
  double log10=log(10.);
  double emax=-numeric_limits<double>::max();
  double emin=numeric_limits<double>::max();

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("y versus t","t","y",0.,1.,0.,1.,&cmap,0,winsize);
  gt.setbgColor("white");

  int n=coarse_number_steps;
  double *ftable=0;
  double *ftable_old=0;
  switch (scheme) {
    case ADAMS_MOULTON:
      ftable=OPERATOR_NEW_BRACKET(double,F77NAME(multistep).order);
      ftable_old=OPERATOR_NEW_BRACKET(double,F77NAME(multistep).order);
      break;
    case EXPLICIT_MIDPOINT:
      break;
    case BACKWARD_DIFFERENTIATION_FORMULA:
      ftable=OPERATOR_NEW_BRACKET(double,F77NAME(multistep).order+1);
      ftable_old=
        OPERATOR_NEW_BRACKET(double,F77NAME(multistep).order+1);
      break;
    case ADAMS_BASHFORTH:
    default:
      ftable=OPERATOR_NEW_BRACKET(double,F77NAME(multistep).order);
  }

  for (int i=0;i<=number_refinements;i++,n*=ratio) {
    F77NAME(multistep).step=0;
    switch (scheme) {
      case ADAMS_MOULTON:
        F77_NAME(adams_moulton_startup)();
        break;
      case EXPLICIT_MIDPOINT:
        break;
      case BACKWARD_DIFFERENTIATION_FORMULA:
        F77_NAME(bdf_startup)();
        break;
      case ADAMS_BASHFORTH:
      default:
        F77_NAME(adams_bashforth_startup)();
        break;
    }
    double t=0.;
    double y=y_init;
    double dt=tstop/static_cast<double>(n);
    double tmin=t;
    double tmax=t;
    double ymin=y_init;
    double ymax=y_init;
    double *tarray=OPERATOR_NEW_BRACKET(double,n+1);
    double *yarray=OPERATOR_NEW_BRACKET(double,n+1);
    tarray[0]=t;
    yarray[0]=y;

#ifdef DEBUG
    cout << "\tnumber steps = " << n << endl;
//  cout << "\tt,y[" << F77NAME(multistep).step << "] = " 
//       << t << " " << y << endl;
#endif
    double yold=y;
    for (;F77NAME(multistep).step<n;t+=dt) {
//    TRACER_CALL(z,"timestep loop");
      switch (scheme) {
        case ADAMS_MOULTON: {
//        TRACER_CALL(z,"ADAMS_MOULTON");
          F77_NAME(adams_moulton)(dt,t,y,dfdy,f,ftable,ftable_old);
          break;
        }
        case EXPLICIT_MIDPOINT:
          if (F77NAME(multistep).step==0) {
//          forward Euler startp
            y=yold+dt*f(t,yold);
            F77NAME(multistep).step++;
          } else {
            F77_NAME(explicit_midpoint)(dt,t,y,yold,f);
          }
          break;
        case BACKWARD_DIFFERENTIATION_FORMULA:
          F77NAME(bdf)(dt,t,y,dfdy,f,ftable,ftable_old);
          break;
        case ADAMS_BASHFORTH:
        default:
          F77_NAME(adams_bashforth)(dt,t,y,f,ftable);
          break;
      }
      if (F77NAME(multistep).step<=F77NAME(multistep).order
      && scheme!=EXPLICIT_MIDPOINT) {
//      fix startup values to get correct order
        y=exact(t+dt);
        switch (scheme) {
          case ADAMS_MOULTON: {
            double fnew=f(t+dt,y);
            ftable[F77NAME(multistep).order-1]=fnew;
            double temp=ftable_old[F77NAME(multistep).order-1];
            for (int i=1;i<F77NAME(multistep).step;i++) {
              fnew=ftable[F77NAME(multistep).order-i]-temp;
              temp=ftable_old[F77NAME(multistep).order-1-i];
              ftable[F77NAME(multistep).order-1-i]=fnew;
            }
            break;
          }
          case BACKWARD_DIFFERENTIATION_FORMULA: {
            double temp=ftable_old[F77NAME(multistep).order];
            ftable[F77NAME(multistep).order]=y;
            for (int i=1;i<=F77NAME(multistep).step;i++) {
              double ytnew=ftable[F77NAME(multistep).order+1-i]-temp;
              temp=ftable_old[F77NAME(multistep).order-i];
              ftable[F77NAME(multistep).order-i]=ytnew;
            }
            break;
          }
        }
      }
      if (scheme==EXPLICIT_MIDPOINT && F77NAME(multistep).step==n-1) {
//      smoothing operation
        y=0.5*(y+yold+dt*f(t+dt,y));
      }
      tarray[F77NAME(multistep).step]=t+dt;
      yarray[F77NAME(multistep).step]=y;
      ymin=min(ymin,y);
      ymax=max(ymax,y);
#ifdef DEBUG
//    cout << "\tt,y[" << F77NAME(multistep).step << "] = " 
//         << t << " " << y << endl;
//    cout << "\ttmin,tmax = " << tmin << " " << tmax << endl;
#endif
    }
    tmax=t;
    gt.newPage();
    gt.rescale(tmin,tmax,ymin,ymax);
    gt.setfgColor("black");
    gt.drawAxes();
    gt.setfgColor("blue");
    gt.movePen(tarray[0],yarray[0]);
    for (int j=1;j<=n;j++) {
      gt.drawLine(tarray[j],yarray[j]);
    }
    gt.setfgColor("red");
    gt.movePen(tarray[0],exact(tarray[0]));
    for (int j=1;j<=n;j++) {
      gt.drawLine(tarray[j],exact(tarray[j]));
    }
    gt.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    double answer=exact(tarray[n]);
    double error=abs(answer-yarray[n]);
    if (error<=0.) error=abs(answer)*numeric_limits<double>::epsilon();
#ifdef DEBUG
//  cout << "\ty[" << n << "] = " << yarray[n] << " " << answer 
//       << " " << answer-yarray[n] << endl;
//  cout << "\terror[" << n << "] = " << error << endl;
#endif
    log_error[i]=log(error)/log10;
    log_n[i]=log(static_cast<double>(n))/log10;
    emax=max(emax,log_error[i]);
    emin=min(emin,log_error[i]);

    delete [] yarray;
    delete [] tarray;
  }
  if (ftable!=0) delete [] ftable;
  if (ftable_old!=0) delete [] ftable_old;

  if (number_refinements>0) {
    XYGraphTool gte("log_10(error) versus log_10(nsteps)","log_10(nsteps)",
      "log_10(error)",log_n[0],log_n[number_refinements],emin,emax,&cmap,
      0,winsize);
    gte.setbgColor("white");
    gte.setfgColor("black");
    gte.drawAxes();
    gte.setfgColor("blue");
    gte.movePen(log_n[0],log_error[0]);
    for (int i=1;i<=number_refinements;i++) {
      gte.drawLine(log_n[i],log_error[i]);
    }
    gte.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }
  delete [] log_error;
  delete [] log_n;
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
