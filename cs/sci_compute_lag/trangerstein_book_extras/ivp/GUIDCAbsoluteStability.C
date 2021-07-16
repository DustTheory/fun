#include <cmath>
#include <iostream>
#include <math.h> // for HUGE_VAL,M_PI
#include <rfftw.h>
#include <stdlib.h>

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "Quadrature.H"
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

enum SCHEME{EULER,BACKWARD_EULER};
SCHEME scheme=EULER;
int ischeme=static_cast<int>(scheme);

int nplot=16384;
int order=1;
//int number_refinements=0;
//int coarse_number_steps=1;
//int ratio=2;
complex<double> rate(1.,0.);
//double tstop=1.;
complex<double> z(1.,0.);
double *nodes=0;
double *integrals=0;
double *weights=0;
complex<double> *fwork=0;
complex<double> *fnewwork=0;
double *twork=0;
complex<double> *ywork=0;
complex<double> *ynewwork=0;

extern "C" {
  void F77_NAME(deferred_correctionbe)(const int &order,
    const double *nodes,const double *integrals,const double &dt,
    const double &t,complex<double> &y,
    complex<double> (*f)(const double &t,const complex<double> &y),
    complex<double> (*dfdy)(const double &t,const complex<double> &y),
    complex<double> *fwork,
    complex<double> *fnewwork,double *twork,complex<double> *ywork,
    complex<double> *ynewwork);
  void F77_NAME(deferred_correctione)(const int &order,
    const double *nodes,const double *integrals,const double &dt,
    const double &t,complex<double> &y,
    complex<double> (*f)(const double &t,const complex<double> &y),
    complex<double> *fwork,
    complex<double> *fnewwork,double *twork,complex<double> *ywork,
    complex<double> *ynewwork);
  void F77NAME(hybrd)(
    void (*fcn)(const int &n,double *x,double *fvec,int &iflag),
    const int &n,double *x,double *fvec,const double &xtol,
    const int &maxfev,const int &ml,const int &mu,const double &epsfcn,
    double *diag,const int &mode,const double &factor,const int &nprint,
    int &info,int &nfev,double *fjac,const int &ldfjac,double *r,
    const int &lr,double *qtf,double *wa1,double *wa2,double *wa3,
    double *wa4);
  void F77_NAME(lagrange_polynomial_integrals)(const int &order,
    const double *nodes,const double *weights, double *integrals);
}

//struct machine_common {
//  double roundoff,small,huge,undefind;
//};
//extern machine_common F77NAME(machine);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
extern "C" {
  complex<double> f(const double &t,const complex<double> &y) {
    return rate*y;
  }
  complex<double> dfdy(const double &t,const complex<double> &y) {
    return rate;
  }
//complex<double> exact(const double &t) { return exp(rate*t); }
  void fcn(const int &n,double *x,double *fvec,int &iflag) {
    rate=complex<double>(x[0],x[1]);
    complex<double> y(1.,0.);
    double dt=1.;
    double t=0.;
    switch (scheme) {
      case BACKWARD_EULER: {
        F77_NAME(deferred_correctionbe)(order,nodes,integrals,dt,t,y,f, 
          dfdy, fwork,fnewwork,twork,ywork,ynewwork);
        break;
      }
      case EULER:
      default: {
        F77_NAME(deferred_correctione)(order,nodes,integrals,dt,t,y,f, 
          fwork,fnewwork,twork,ywork,ynewwork);
        break;
      }
    }
    fvec[0]=y.real()-z.real();
    fvec[1]=y.imag()-z.imag();
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

//{ const char *group="Initial Value Problem Parameters";
//  main_list->append(OPERATOR_NEW GUIInputParameter<double>(rate,
//    "rate",-DBL_MAX,DBL_MAX,group) );
//  main_list->append(OPERATOR_NEW GUIInputParameter<double>(tstop,
//    "tstop",0.,DBL_MAX,group) );
//  main_list->append(OPERATOR_NEW GUIInputParameter<double>(
//    y_init,"y_init",-DBL_MAX,DBL_MAX,group) );
//}

  { const char *group="Numerical Method Parameters";
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ischeme,
      "scheme","euler","backward euler"));
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nplot,
      "nplot",2,INT_MAX,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(
      order,"order",1,10,group) );
//  main_list->append(OPERATOR_NEW GUIInputParameter<int>(
//    coarse_number_steps,"coarse_number_steps",1,INT_MAX,group) );
//  main_list->append(OPERATOR_NEW GUIInputParameter<int>(
//    number_refinements,"number_refinements",0,INT_MAX,group) );
//  main_list->append(OPERATOR_NEW GUIInputParameter<int>(ratio,
//    "ratio",2,INT_MAX,group) );
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
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void runMain(bool /*called_before*/) {
//TRACER_CALL(t,"runMain");
  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("absolute stability boundary","Real Part(lambda)",
    "Imaginary Part(lambda)",0.,1.,0.,1.,&cmap,0,winsize);
  gt.setbgColor("white");

#ifdef DEBUG
//cout << "\torder = " << order << endl;
#endif
  nodes=OPERATOR_NEW_BRACKET(double,order+1);
  integrals=OPERATOR_NEW_BRACKET(double,(order)*(order+1));
  weights=OPERATOR_NEW_BRACKET(double,order+1);
  fwork=OPERATOR_NEW_BRACKET(complex<double>,order+1);
  fnewwork=OPERATOR_NEW_BRACKET(complex<double>,order+1);
  twork=OPERATOR_NEW_BRACKET(double,order+1);
  ywork=OPERATOR_NEW_BRACKET(complex<double>,order+1);
  ynewwork=OPERATOR_NEW_BRACKET(complex<double>,order+1);
#ifdef INDEF
  for (int i=0;i<=order;i++) nodes[i]=HUGE_VAL;
  for (int i=0;i<order*(order+1);i++) integrals[i]=HUGE_VAL;
  for (int i=0;i<order+1;i++) {
    weights[i]=HUGE_VAL;
    fwork[i]=complex<double>(HUGE_VAL,HUGE_VAL);
    fnewwork[i]=complex<double>(HUGE_VAL,HUGE_VAL);
    twork[i]=HUGE_VAL;
    ywork[i]=complex<double>(HUGE_VAL,HUGE_VAL);
    ynewwork[i]=complex<double>(HUGE_VAL,HUGE_VAL);
  }
#endif
  LobattoQuadrature<1> lq(order+1); // quadrature on [0,1]
  for (int i=0;i<=order;i++) { // map quadrature to [-1,1]
    nodes[i]=2.*lq.point(i)[0]-1.;
    weights[i]=2.*lq.weight(i);
  }
  F77_NAME(lagrange_polynomial_integrals)(order,nodes,weights,integrals);

  double dt=1.;
  double dtheta=32.*M_PI/static_cast<double>(nplot);
  complex<double> *point=OPERATOR_NEW_BRACKET(complex<double>,nplot+1);
  complex<double> *point9=OPERATOR_NEW_BRACKET(complex<double>,nplot+1);
  complex<double> *point1=OPERATOR_NEW_BRACKET(complex<double>,nplot+1);
  double *xpoint=OPERATOR_NEW_BRACKET(double,nplot+1);
  double *ypoint=OPERATOR_NEW_BRACKET(double,nplot+1);

  double x[2]={HUGE_VAL,HUGE_VAL};
  double fvec[2]={HUGE_VAL,HUGE_VAL};
  double xtol=sqrt(DBL_EPSILON);
  int maxfev=25;
  int ml=2;
  int mu=2;
  double epsfcn=sqrt(DBL_EPSILON);
  double diag[2];
  int mode=1;
  double factor=100.;
  int nprint=0;
  int info=0;
  int nfev=0;
  double fjac[4];
  int ldfjac=2;
  double r[3];
  int lr=3;
  double qtf[2];
  double wa1[2];
  double wa2[2];
  double wa3[2];
  double wa4[2];

  double theta=0.;
  double xmin=HUGE_VAL;
  double xmax=-HUGE_VAL;
  double ymin=HUGE_VAL;
  double ymax=-HUGE_VAL;
  point[0]=complex<double>(0.,0.);
  for (int j=0;j<=nplot;j++,theta+=dtheta) {
#ifdef INDEF
    for (int i=0;i<order+1;i++) {
      fwork[i]=complex<double>(HUGE_VAL,HUGE_VAL);
      fnewwork[i]=complex<double>(HUGE_VAL,HUGE_VAL);
      twork[i]=HUGE_VAL;
      ywork[i]=complex<double>(HUGE_VAL,HUGE_VAL);
      ynewwork[i]=complex<double>(HUGE_VAL,HUGE_VAL);
    }
#endif
    if (j>0) point[j]=point[j-1];
    complex<double> znew(cos(theta),sin(theta));
    z=znew;

    x[0]=point[j].real();
    x[1]=point[j].imag();
    F77NAME(hybrd)(fcn,2,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,mode,factor,
      nprint,info,nfev,fjac,ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4);
    point[j]=complex<double>(x[0],x[1]);
    xmin=min(xmin,x[0]);
    xmax=max(xmax,x[0]);
    ymin=min(ymin,x[1]);
    ymax=max(ymax,x[1]);
  }
  point9[0]=complex<double>(0.,0.);
  for (int j=0;j<=nplot;j++,theta+=dtheta) {
    if (j>0) point9[j]=point9[j-1];
    complex<double> znew(cos(theta),sin(theta));
    z=0.9*znew;

    x[0]=point9[j].real();
    x[1]=point9[j].imag();
    F77NAME(hybrd)(fcn,2,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,mode,factor,
      nprint,info,nfev,fjac,ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4);
    point9[j]=complex<double>(x[0],x[1]);
    xmin=min(xmin,x[0]);
    xmax=max(xmax,x[0]);
    ymin=min(ymin,x[1]);
    ymax=max(ymax,x[1]);
  }
  point1[0]=complex<double>(0.,0.);
  for (int j=0;j<=nplot;j++,theta+=dtheta) {
    if (j>0) point1[j]=point1[j-1];
    complex<double> znew(cos(theta),sin(theta));
    z=1.1*znew;

    x[0]=point1[j].real();
    x[1]=point1[j].imag();
    F77NAME(hybrd)(fcn,2,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,mode,factor,
      nprint,info,nfev,fjac,ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4);
    point1[j]=complex<double>(x[0],x[1]);
    xmin=min(xmin,x[0]);
    xmax=max(xmax,x[0]);
    ymin=min(ymin,x[1]);
    ymax=max(ymax,x[1]);
  }
  gt.newPage();
  gt.rescale(xmin,xmax,ymin,ymax);
  gt.setfgColor("blue");
  gt.movePen(point[0].real(),point[0].imag());
  for (int j=1;j<=nplot;j++) {
    gt.drawLine(point[j].real(),point[j].imag());
  }
  gt.setfgColor("green");
  gt.movePen(point9[0].real(),point9[0].imag());
  for (int j=1;j<=nplot;j++) {
    gt.drawLine(point9[j].real(),point9[j].imag());
  }
  gt.setfgColor("red");
  gt.movePen(point1[0].real(),point1[0].imag());
  for (int j=1;j<=nplot;j++) {
    gt.drawLine(point1[j].real(),point1[j].imag());
  }
//gt.colorPolygon(nplot+1,xpoint,ypoint);
  gt.setfgColor("black");
  gt.drawAxes();
  gt.flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

  delete [] nodes;
  delete [] integrals;
  delete [] weights;
  delete [] fwork;
  delete [] fnewwork;
  delete [] twork;
  delete [] ywork;
  delete [] ynewwork;
  delete [] point1; point1=0;
  delete [] point9; point9=0;
  delete [] point; point=0;
  delete [] xpoint; xpoint=0;
  delete [] ypoint; ypoint=0;
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
//F77NAME(machine).roundoff=DBL_EPSILON;
//F77NAME(machine).small=DBL_MIN;
//F77NAME(machine).huge=DBL_MAX;
//F77NAME(machine).undefind=HUGE_VAL;

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
