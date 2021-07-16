#include <iostream>
#include <limits>
#include <math.h> // for HUGE_VAL,M_PI
#include <stdlib.h>

using namespace std;

#include "Arch.H"
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
  double F77_NAME(c0_spline_eval)(const int &element,const int &nelements,
    const int &order,const double &pt,const double *x,const double *xi,
    const double *y,double *p,double *q);
  void F77_NAME(c1_quadratic_coefs)(const int &nelements,const double *x,
    const double *y,const double *yp,double *z);
  double F77_NAME(c1_quadratic_eval)(const int &element,
    const int &nelements,const double &pt,const double *x,const double *y,
    double *z);
  double F77_NAME(c1_spline_eval)(const int &element,const int &nelements,
    const int &order,const double &pt,const double &sigma,const double *x,
    const double *xi,const double *y,const double *yp);
  void F77_NAME(c2_cubic_coefs)(const int &nelements,double *d,double *e,
    const double *x,const double *y,const double *yp,double *z);
  double F77_NAME(c2_cubic_eval)(const int &element,const int &nelements,
    const double &pt,const double *x,const double *y,double *z);
  void F77_NAME(c2_quartic_coefs)(const int &nelements,const double *x,
    const double *y,const double *yp,double *ypp,double *z);
  double F77_NAME(c2_quartic_eval)(const int &element,const int &nelements,
    const double &pt,const double *x,const double *y,const double *yp,
    double *z);
  double F77_NAME(c2_spline_eval)(const int &element,const int &nelements,
    const int &order,const double &pt,const double &sigma,
    const double &tau,const double *x,const double *xi,const double *y,
    const double *yp,const double *ypp);
}

GUI_INPUT_PARAMETER_LIST_TYPE *main_list=0;
bool skip_gui=FALSE;
char *display_name=0;
double winsize=0.5;

enum TEST_FUNCTION{PIECEWISE_CONSTANT,PIECEWISE_LINEAR,
  PIECEWISE_QUADRATIC,SMOOTH};
TEST_FUNCTION test_function=SMOOTH;
int itest_function=test_function;
enum POINTS{EQUIDISTANT,RANDOM};
POINTS points=EQUIDISTANT;
int ipoints=static_cast<int>(points);

int nmax=3;
int order=1;
int continuity=0;
double jump0=sqrt(0.2);
double jump1=1.;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double fcn(double t) { 
  switch (test_function) {
    case PIECEWISE_CONSTANT:
      return (t<jump0 ? 0. : 1.);
    case PIECEWISE_LINEAR:
      return (t<jump0 ? jump1*t/jump0 :
              1.-(1.-jump1)*(1.-t)/(1.-jump0));
    case PIECEWISE_QUADRATIC: {
      double den=1./(jump0*(1.-jump0));
      double c=(jump0-jump1)*den;
      double bl=(jump1-jump0*jump0)*den;
      double br=(jump1+jump0*(-2.+jump0))*den;
      return (t<jump0 ? t*(bl+c*t) :
              1.+(1.-t)*(br+c-c*t));
    }
    case SMOOTH:
//    return sin(M_PI_2*t);
      return 1./(1.+25.*t*t); // runge example
    default:
      abort();
      return HUGE_VAL;
  }
}
double fcn_deriv(double t) { 
  switch (test_function) {
    case PIECEWISE_CONSTANT:
      return 0.;
    case PIECEWISE_LINEAR:
      return (t<jump0 ? jump1/jump0 : (1.-jump1)/(1.-jump0));
    case PIECEWISE_QUADRATIC: {
      double den=1./(jump0*(1.-jump0));
      double c=(jump0-jump1)*den;
      double bl=(jump1-jump0*jump0)*den;
      double br=(jump1+jump0*(-2.+jump0))*den;
      return (t<jump0 ? bl+2.*c*t : -br-2.*c*(1.-t));
    }
    case SMOOTH:
//    return M_PI_2*cos(M_PI_2*t);
      return -50.*t/__cmath_power(1.+25.*t*t,2);
    default:
      abort();
      return HUGE_VAL;
  }
}
double fcn_2nd_deriv(double t) { 
  switch (test_function) {
    case PIECEWISE_CONSTANT:
    case PIECEWISE_LINEAR:
      return 0.;
    case PIECEWISE_QUADRATIC: {
      return 2.*(jump0-jump1)/(jump0*(1.-jump0));;
    }
    case SMOOTH: {
//    return M_PI_2*M_PI_2*sin(M_PI_2*t);
      double den=1.+25.*t*t;
      return 50.*(75.*t*t-1.)/(den*den*den);
    }
    default:
      abort();
      return HUGE_VAL;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//selection sort:
void sort(const int &n,double *x) {
  int j=0;
  int jnext=1;
  for (;jnext<=n;j=jnext,jnext++) {
    int jsmall=j;
    for (int k=jnext;k<=n;k++) {
      if (x[jsmall]>x[k]) jsmall=k;
    }
    if (j!=jsmall) {
      double t=x[j];
      x[j]=x[jsmall];
      x[jsmall]=t;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int bin_search(int high,int low, double t,const double *x) {
  for (; high-low>1;) {
    int i=(high+low)/2;
    if (t<x[i]) high=i;
    else low=i;
  }
  return low;
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
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      itest_function,"test function","piecewise constant",
      "piecewise linear","piecewise quadratic","smooth") );
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,ipoints,
      "points","equidistant","random") );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nmax,
      "max number elements",2,100,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(order,
      "order",1,11,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(continuity,
      "continuity",0,3,group) );
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(
      jump0,"discontinuity x location",0.,1.,group));
    main_list->append(OPERATOR_NEW GUIInputParameter<double>(
      jump1,"discontinuity y location",0.,1.,group));
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
  test_function=static_cast<TEST_FUNCTION>(itest_function);
  points=static_cast<POINTS>(ipoints);
  continuity=min(continuity,min(2,order-1));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void runMain(bool /*called_before*/) {
  TRACER_CALL(zz,"runMain");
  double xmin=-1.;
  double xmax=1.;
  double ymin=-1.;
  double ymax=1.;
  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("interpolation","x","y",xmin,xmax,ymin,ymax,&cmap,0,0.5);
  gt.setbgColor("white");

  int nmin=2;

  double *errors_alloc=OPERATOR_NEW_BRACKET(double,nmax-nmin+1);
  double *errors=errors_alloc-nmin;
  double *x=OPERATOR_NEW_BRACKET(double,nmax+1);
  double *y=
    OPERATOR_NEW_BRACKET(double,nmax*max(1,order-2*continuity)+1);
  double *yp=OPERATOR_NEW_BRACKET(double,
    ((order>continuity+1 && continuity>0) ? nmax+1 : continuity));
  double *ypp=OPERATOR_NEW_BRACKET(double,
    ((order>2*continuity && continuity>1) ? nmax+1 : 1 ));
  double *xi=OPERATOR_NEW_BRACKET(double,order-2*continuity+1);
  double *z=OPERATOR_NEW_BRACKET(double,nmax+1);

  double errors_min=numeric_limits<double>::max();
  double errors_max=-numeric_limits<double>::max();
  double log10=log(10.);
  for (int i=nmin;i<=nmax;i++) {
    errors[i]=numeric_limits<double>::max();
  }
  TimedObject tl("Lagrange interpolation");
  TimedObject tn("Newton interpolation");
  for (int n=nmin;n<=nmax;n++) {
#ifdef INDEF
    for (int i=0;i<=n;i++) {
      x[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<=n*(order-continuity)+continuity;i++) {
      y[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;
    i<=((order>continuity+1 && continuity>0) ? n : continuity-1);i++) {
      yp[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<=((order>2*continuity && continuity>1) ? n : 0 );i++)
    {
      ypp[i]=numeric_limits<double>::infinity();
    }
    for (int i=0;i<=n;i++) {
      z[i]=numeric_limits<double>::infinity();
    }
#endif
//  generate mesh
    int ninterior=order-2*continuity;
    if (points==EQUIDISTANT) {
      double dx=(xmax-xmin)/static_cast<double>(n);
      x[0]=xmin;
      for (int i=1;i<n;i++) x[i]=x[i-1]+dx;
      x[n]=xmax;
      if (ninterior>0) {
        double dxi=1./static_cast<double>(ninterior);
        xi[0]=0.;
        for (int j=1;j<ninterior;j++) xi[j]=xi[j-1]+dxi;
        xi[ninterior]=1.;
      }
    } else {
      x[0]=xmin;
      for (int j=1;j<n;j++) {
        x[j]=xmin+(xmax-xmin)
            *static_cast<double>(rand())/static_cast<double>(RAND_MAX);
      }
      x[n]=xmax;
      sort(n,x);
      if (ninterior>0) {
        for (int j=1;j<=ninterior/2;j++) {
          xi[j]=0.5*static_cast<double>(rand())
               /static_cast<double>(RAND_MAX);
        }
        if (ninterior%2==0) xi[ninterior/2]=0.5;
        xi[0]=0.;
        sort(order/2-continuity,xi+1);
//      place symmetrically:
        for (int j=order/2-continuity+1;j<=order-2*continuity;j++) {
          xi[j]=1.-xi[order-2*continuity-j];
        }
      }
    }
#ifdef DEBUG
//  for (int j=0;j<=n;j++) {
//    cout << "x[" << j << "] = " << x[j] << endl;
//  }
//  if (ninterior>0) {
//    for (int j=0;j<=n;j++) {
//      cout << "xi[" << j << "] = " << xi[j] << endl;
//    }
//  }
#endif
    int nint=max(1,ninterior);
    for (int j=0;j<=n;j++) {
      y[j*nint]=fcn(x[j]);
    }
    if (ninterior>0) {
      for (int j=0;j<n;j++) {
        double dx=x[j+1]-x[j];
        int i=j*ninterior;
        for (int k=1;k<ninterior;k++) {
          y[i+k]=fcn(x[j]+dx*xi[k]);
        }
      }
    }
#ifdef DEBUG
//  for (int i=0;i<=n*max(1,order-2*continuity);j++) {
//    cout << "y[" << i << "] = " << y[i] << endl;
//  }
#endif
    if (continuity>0) {
      if (order>continuity+1) {
        for (int j=0;j<=n;j++) yp[j]=fcn_deriv(x[j]);
#ifdef DEBUG
//      for (int j=0;j<=n;j++) {
//        cout << "yp[" << j << "] = " << yp[j] << endl;
//      }
#endif
      } else {
        yp[0]=fcn_deriv(xmin);
        if (continuity==2) yp[1]=fcn_deriv(xmax);
#ifdef DEBUG
//      cout << "yp[0] = " << yp[0] << endl;
//      if (continuity==2) cout << "yp[1] = " << yp[1] << endl;
#endif
      }
    }
    if (continuity>1) {
      if (order>2*continuity) {
        for (int j=0;j<=n;j++) ypp[j]=fcn_2nd_deriv(x[j]);
#ifdef DEBUG
//      for (int j=0;j<=n;j++) {
//        cout << "ypp[" << j << "] = " << ypp[j] << endl;
//      }
#endif
      } else {
        ypp[0]=fcn_2nd_deriv(xmin);
#ifdef DEBUG
//      cout << "ypp[0] = " << ypp[0] << endl;
#endif
      }
    }
    double t=xmin;
    double yy=fcn(t);
    double ylo=yy;
    double yhi=yy;
    double dx=(xmax-xmin)/128.;
    for (int j=0;j<=128*n;j++) {
      t=xmin+dx*static_cast<double>(j);
      yy=fcn(t);
      ylo=min(ylo,yy);
      yhi=max(yhi,yy);
    }

    gt.newPage();
    gt.rescale(xmin,xmax,ylo,yhi);
    gt.setfgColor("black");
    gt.drawAxes();
    gt.setfgColor("blue");
    gt.movePen(xmin,fcn(xmin));
    for (int j=0;j<=128*n;j++) {
      t=xmin+dx*static_cast<double>(j);
      gt.drawLine(t,fcn(t));
    }
    double maxerr=0.;
    gt.setfgColor("red");
    if (continuity==0) {
      double *p=OPERATOR_NEW_BRACKET(double,order+1);
      double *q=OPERATOR_NEW_BRACKET(double,order);
      int element=0;
      double yy=F77_NAME(c0_spline_eval)(element,n,order,xmin,x,xi,y,p,q);
      gt.movePen(xmin,yy);
      double dx=(xmax-xmin)/static_cast<double>(128*n);
      double t=xmin+dx;
      for (int j=1;j<=128*n;j++,t+=dx) {
        if (j==128*n) t=xmax;
        element=bin_search(n,element, t,x);
        yy=F77_NAME(c0_spline_eval)(element,n,order,t,x,xi,y,p,q);
        gt.drawLine(t,yy);
        maxerr=max(maxerr,abs(yy-fcn(t)));
      }
      delete [] p; p=0;
      delete [] q; q=0;
    } else if (continuity==1) {
      if (order==2) {
        F77_NAME(c1_quadratic_coefs)(n,x,y,yp,z);
        int element=0;
        double yy=F77_NAME(c1_quadratic_eval)(element,n,xmin,x,y,z);
        gt.movePen(xmin,yy);
        double dx=(xmax-xmin)/static_cast<double>(128*n);
        double t=xmin+dx;
        for (int j=1;j<=128*n;j++,t+=dx) {
          if (j==128*n) t=xmax;
          element=bin_search(n,element, t,x);
          yy=F77_NAME(c1_quadratic_eval)(element,n,t,x,y,z);
          gt.drawLine(t,yy);
          maxerr=max(maxerr,abs(yy-fcn(t)));
        }
      } else {
        int element=0;
        double sum=2.;
        for (int i=1;i<=order-3;i++) sum+=1./xi[i];
        double yy=
          F77_NAME(c1_spline_eval)(element,n,order,xmin,sum,x,xi,y,yp);
        gt.movePen(xmin,yy);
        double dx=(xmax-xmin)/static_cast<double>(128*n);
        double t=xmin+dx;
        for (int j=1;j<=128*n;j++,t+=dx) {
          if (j==128*n) t=xmax;
          element=bin_search(n,element, t,x);
          yy=F77_NAME(c1_spline_eval)(element,n,order,t,sum,x,xi,y,yp);
          gt.drawLine(t,yy);
          maxerr=max(maxerr,abs(yy-fcn(t)));
        }
      }
    } else {
      if (order==3) {
        double *d=OPERATOR_NEW_BRACKET(double,n+1);
        double *e=OPERATOR_NEW_BRACKET(double,n);
        F77_NAME(c2_cubic_coefs)(n,d,e,x,y,yp,z);
        delete [] d; d=0;
        delete [] e; e=0;
        int element=0;
        double yy=F77_NAME(c2_cubic_eval)(element,n,xmin,x,y,z);
        gt.movePen(xmin,yy);
        double dx=(xmax-xmin)/static_cast<double>(128*n);
        double t=xmin+dx;
        for (int j=1;j<=128*n;j++,t+=dx) {
          if (j==128*n) t=xmax;
          element=bin_search(n,element, t,x);
          double yy=F77_NAME(c2_cubic_eval)(element,n,t,x,y,z);
          gt.drawLine(t,yy);
          maxerr=max(maxerr,abs(yy-fcn(t)));
        }
      } else if (order==4) {
        F77_NAME(c2_quartic_coefs)(n,x,y,yp,ypp,z);
        int element=0;
        double yy=F77_NAME(c2_quartic_eval)(element,n,xmin,x,y,yp,z);
        gt.movePen(xmin,yy);
        double dx=(xmax-xmin)/static_cast<double>(128*n);
        double t=xmin+dx;
        for (int j=1;j<=128*n;j++,t+=dx) {
          if (j==128*n) t=xmax;
          element=bin_search(n,element, t,x);
          double yy=F77_NAME(c2_quartic_eval)(element,n,t,x,y,yp,z);
          gt.drawLine(t,yy);
          maxerr=max(maxerr,abs(yy-fcn(t)));
        }
      } else {
        int element=0;
        double sigma=0.;
        double tau=0.;
        for (int i=1;i<=order-5;i++) {
          double sum=0.;
          for (int k=1;k<=order-5;k++) {
            if (k!=i) sum+=1./xi[k];
          }
          sigma+=1./xi[i];
          tau+=sum/xi[i];
        }
        double yy=F77_NAME(c2_spline_eval)(element,n,order,xmin,sigma,tau,
          x,xi,y,yp,ypp);
        gt.movePen(xmin,yy);
        double dx=(xmax-xmin)/static_cast<double>(128*n);
        double t=xmin+dx;
        for (int j=1;j<=128*n;j++,t+=dx) {
          if (j==128*n) t=xmax;
          element=bin_search(n,element, t,x);
          yy=F77_NAME(c2_spline_eval)(element,n,order,t,sigma,tau,
            x,xi,y,yp,ypp);
          gt.drawLine(t,yy);
          maxerr=max(maxerr,abs(yy-fcn(t)));
        }
      }
    }
    maxerr=max(maxerr,numeric_limits<double>::epsilon());
    errors[n]=log(maxerr)/log10;
    errors_min=min(errors_min,errors[n]);
    errors_max=max(errors_max,errors[n]);
    gt.flush();
    cout << "errors[" << n << "] = " << maxerr << endl;

    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }

  if (nmax>nmin) {
    XYGraphTool gt2("log error vs log number cells",
      "log_10(number cells)","log_10(error)",
      log(static_cast<double>(nmin))/log10,
      log(static_cast<double>(nmax))/log10,errors_min,errors_max,&cmap,
      0,0.5);
    gt2.setbgColor("white");
    gt2.setfgColor("black");
    gt2.drawAxes();
    gt2.setfgColor("blue");
    double t=log(static_cast<double>(nmin))/log10;
    gt2.movePen(t,errors[nmin]);
    for (int n=nmin+1;n<nmax;n++) {
      t=log(static_cast<double>(n))/log10;
      gt2.drawLine(t,errors[n]);
    }
    gt2.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  }

  delete [] z; z=0;
  delete [] ypp; ypp=0;
  delete [] yp; yp=0;
  delete [] y; y=0;
  delete [] xi; xi=0;
  delete [] x;
  delete [] errors_alloc;
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
