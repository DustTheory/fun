#include <cmath>
#include <iostream>
#include <limits>
#include <stdlib.h>

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "Quadrature.H"
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

int maxit=20;
int order=1;
int number_refinements=0;
int coarse_number_steps=1;
int ratio=2;
double rate=1.;
double tstop=1.;
double y_init=1.;

extern "C" {
  void F77_NAME(deferred_correction)(const int &order,
    const double *nodes,const double *integrals,const double &dt,
    const double &t,double &y,
    double (*f)(const double &t,const double &y),double *fwork,
    double *fnewwork,double *twork,double *ywork,double *ynewwork);
  void F77_NAME(lagrange_polynomial_integrals)(const int &order,
    const double *nodes,const double *weights, double *integrals);
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
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      maxit,"maxit",1,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      order,"order",1,10,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      coarse_number_steps,"coarse_number_steps",1,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      number_refinements,"number_refinements",0,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(ratio,
      "ratio",2,INT_MAX,group) );
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
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
#ifdef DEBUG
//cout << "\torder = " << order << endl;
//cout << "\ttstop = " << tstop << endl;
#endif
  double *nodes=OPERATOR_NEW_BRACKET(double,order+1);
  double *integrals=OPERATOR_NEW_BRACKET(double,(order)*(order+1));
  double *weights=OPERATOR_NEW_BRACKET(double,order+1);
  double *fwork=OPERATOR_NEW_BRACKET(double,order+1);
  double *fnewwork=OPERATOR_NEW_BRACKET(double,order+1);
  double *twork=OPERATOR_NEW_BRACKET(double,order+1);
  double *ywork=OPERATOR_NEW_BRACKET(double,order+1);
  double *ynewwork=OPERATOR_NEW_BRACKET(double,order+1);
#ifdef INDEF
  for (int i=0;i<=order;i++) nodes[i]=numeric_limits<double>::max();
  for (int i=0;i<order*(order+1);i++) {
    integrals[i]=numeric_limits<double>::max();
  }
  for (int i=0;i<order+1;i++) {
    weights[i]=numeric_limits<double>::max();
    fwork[i]=numeric_limits<double>::max();
    fnewwork[i]=numeric_limits<double>::max();
    twork[i]=numeric_limits<double>::max();
    ywork[i]=numeric_limits<double>::max();
    ynewwork[i]=numeric_limits<double>::max();
  }
#endif
  LobattoQuadrature<1> lq(order+1); // quadrature on [0,1]
  for (int i=0;i<=order;i++) { // map quadrature to [-1,1]
    nodes[i]=2.*lq.point(i)[0]-1.;
    weights[i]=2.*lq.weight(i);
  }
#ifdef DEBUG
//for (int i=0;i<=order;i++) {
//  cout << "\tnodes,weights[ " << i << " ] = "
//       << nodes[i] << " " << weights[i] << endl;
//}
#endif
  F77_NAME(lagrange_polynomial_integrals)(order,nodes,weights,integrals);
#ifdef DEBUG
//for (int i=0;i<order;i++) {
//  cout << "\tintegrals from -1 to " << nodes[i+1] << " = ";
//  for (int j=0;j<=order;j++) cout << integrals[i+j*order] << " ";
//  cout << endl;
//}
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

#ifdef DEBUG
  cout << "\tnumber_refinements = " << number_refinements << endl;
#endif
  for (int i=0;i<=number_refinements;i++,n*=ratio) {
    double t=0.;
    double y=y_init;
    double dt=tstop/static_cast<double>(n);
#ifdef DEBUG
//  cout << "\tdt = " << dt << endl;
#endif
    double tmin=t;
    double tmax=t;
    double ymin=y_init;
    double ymax=y_init;
    double *tarray=OPERATOR_NEW_BRACKET(double,n+1);
    double *yarray=OPERATOR_NEW_BRACKET(double,n+1);
    tarray[0]=t;
    yarray[0]=y;

    double yold=y;
    for (int step=0;step<n;step++,t+=dt) {
//    TRACER_CALL(z,"timestep loop");
      F77_NAME(deferred_correction)(order,nodes,integrals,dt,t,y,f, 
        fwork,fnewwork,twork,ywork,ynewwork);
      tarray[step+1]=t+dt;
      yarray[step+1]=y;
      ymin=min(ymin,y);
      ymax=max(ymax,y);
#ifdef DEBUG
//    cout << "\tt,y[" << step << "] = " 
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
//  cout << "\tt[" << n << "] = " << tarray[n] << endl;
    cout << "\terror[" << n << "] = " << error << endl;
#endif
    log_error[i]=log(error)/log10;
    log_n[i]=log(static_cast<double>(n))/log10;
    emax=max(emax,log_error[i]);
    emin=min(emin,log_error[i]);

    delete [] yarray;
    delete [] tarray;
  }

  if (number_refinements>0) {
    XYGraphTool gte("log_10(error) versus log_10(nsteps)",
      "log_10(nsteps)","log_10(error)",log_n[0],log_n[number_refinements],
      emin,emax,&cmap,0,winsize);
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
  delete [] nodes;
  delete [] integrals;
  delete [] weights;
  delete [] fwork;
  delete [] fnewwork;
  delete [] twork;
  delete [] ywork;
  delete [] ynewwork;
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
//F77NAME(machine).roundoff=numeric_limits<double>::epsilon();
//F77NAME(machine).small=numeric_limits<double>::min();
//F77NAME(machine).huge=numeric_limits<double>::max();
//F77NAME(machine).undefind=numeric_limits<double>::max();

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
