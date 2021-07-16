#include <cmath>
#include <iostream>
#include <math.h> // for M_PI
#include <limits>
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

extern "C" {
  void dfftb_(const int&,double*,const double*);
  void dfftf_(const int&,double*,const double*);
  void dffti_(const int&,double*);
}

GUI_INPUT_PARAMETER_LIST_TYPE *main_list=0;
bool skip_gui=FALSE;
char *display_name=0;
double winsize=0.5;

enum FFT_CODE{FFTPACK,FFTWREAL,FFTWCOMPLEX};
FFT_CODE code=FFTWCOMPLEX;
int icode=static_cast<int>(code);
enum PROBLEM{GAUSSIAN,PIECEWISE_CONSTANT};
PROBLEM problem=GAUSSIAN;
int iproblem=static_cast<int>(problem);

int nmax=3;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double fcn(double x) { 
  switch (problem) {
    case GAUSSIAN:
      x-=M_PI;
      return exp(-2.*x*x); // exercise 9.4 #1 
//    return 1.;
//    return cos(x);
//    return sin(x);
//    return cos(2.*x);
    case PIECEWISE_CONSTANT: {
      double pi3=M_PI/3.;
      return (x<pi3 || 5.*pi3 < x ? 0. : (x<M_PI ? 1. : -1.));
    }
    default:
      abort();
      return 0.;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
double fftpack_poly(double theta,int n,double *fhat) {
  double fp=0.;
  double angle=theta;
  for (int i=1;i<n-1;i+=2,angle+=theta) {
    fp+=fhat[i]*cos(angle)-fhat[i+1]*sin(angle);
  }
  fp*=2.;
  fp+=fhat[0];
  if (n%2==0) fp+=fhat[n-1]*cos(angle);
  return fp/static_cast<double>(n);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,icode,
      "code","fft pack","fftw real","fftw complex") );
    main_list->append(OPERATOR_NEW GUIEnumInputParameter(group,iproblem,
      "problem","gaussian","piecewise constant") );
    main_list->append(OPERATOR_NEW GUIInputParameter<int>(nmax,
      "max power of 2",3,13,group) );
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
  problem=static_cast<PROBLEM>(iproblem);
  code=static_cast<FFT_CODE>(icode);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

void runMain(bool /*called_before*/) {
//TRACER_CALL(z,"runMain");
  char index_field[3];
  char filename[LENGTH_NAME];

  double xmin=0.;
  double xmax=2.*M_PI;
  double ymin=-numeric_limits<double>::infinity();
  switch (problem) {
    case GAUSSIAN:
      ymin=0.;
      break;
    case PIECEWISE_CONSTANT:
      ymin=-1.;
      break;
    default:
      abort();
  }
  double ymax=1.;
  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("trig polynomial","x","y",xmin,xmax,ymin,ymax,&cmap,0,0.5);
  gt.setbgColor("white");

  int nmin=2;
  int power2=__cmath_power(2,nmax);
#ifdef DEBUG
//cout << "\tnmax = " << nmax << endl;
//cout << "\tpower2 = " << power2 << endl;
//cout << "\tcode = " << code << endl;
#endif
  fftw_plan cplan_backward;
  fftw_plan cplan_forward;
  rfftw_plan rplan_backward;
  rfftw_plan rplan_forward;

  double *df=0;
  double *DFHAT=0;
  double *dwsave=0;
  double *dwsave2=0;
  fftw_complex *cf=0;
  fftw_complex *cfhat=0;
  fftw_complex *CF=0;
  fftw_complex *CFHAT=0;
  fftw_real *rf=0;
  fftw_real *rfhat=0;
  fftw_real *RF=0;
  fftw_real *RFHAT=0;
  switch (code) {
    case FFTPACK:
      df=OPERATOR_NEW_BRACKET(double,power2);
      DFHAT=OPERATOR_NEW_BRACKET(double,power2);
      dwsave=OPERATOR_NEW_BRACKET(double,2*power2+15);
      dwsave2=OPERATOR_NEW_BRACKET(double,2*power2+15);
      dffti_(power2,dwsave2);
      break;
    case FFTWCOMPLEX :
      cf=OPERATOR_NEW_BRACKET(fftw_complex,power2);
      cfhat=OPERATOR_NEW_BRACKET(fftw_complex,power2);
      CF=OPERATOR_NEW_BRACKET(fftw_complex,power2);
      CFHAT=OPERATOR_NEW_BRACKET(fftw_complex,power2);
      cplan_backward=
        fftw_create_plan(power2,FFTW_BACKWARD,FFTW_ESTIMATE);
      break;
    case FFTWREAL:
      rf=OPERATOR_NEW_BRACKET(fftw_real,power2);
      rfhat=OPERATOR_NEW_BRACKET(fftw_real,power2);
      RF=OPERATOR_NEW_BRACKET(fftw_real,power2);
      RFHAT=OPERATOR_NEW_BRACKET(fftw_real,power2);
      rplan_backward=
        rfftw_create_plan(power2,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE);
      break;
    default:
      abort();
  }
  double *errors=OPERATOR_NEW_BRACKET(double,nmax+1);

  double errors_min=numeric_limits<double>::max();
  double errors_max=-numeric_limits<double>::max();
  double log10=log(10.);
  for (int i=0;i<=nmax;i++) {
    errors[i]=numeric_limits<double>::infinity();
  }
  for (int n=4,nm=nmin;nm<=nmax;n*=2,nm++) {
//  TRACER_CALL(t,"loop");
#ifdef DEBUG
//  cout << "\n\tn = " << n << endl;
#endif
#ifdef INDEF
    switch (code) {
      case FFTPACK: {
        for (int i=0;i<power2;i++) df[i]=numeric_limits<double>::infinity();
        for (int i=0;i<power2;i++) {
          DFHAT[i]=numeric_limits<double>::infinity();
        }
        break;
      }
      case FFTWCOMPLEX: {
        for (int i=0;i<power2;i++) {
          cf[i].re=cf[i].im=numeric_limits<double>::infinity();
        }
        for (int i=0;i<power2;i++) {
          CFHAT[i].re=CFHAT[i].im=numeric_limits<double>::infinity();
        }
        break;
      }
      case FFTWREAL: {
        for (int i=0;i<power2;i++) {
          rf[i]=numeric_limits<double>::infinity();
          rfhat[i]=numeric_limits<double>::infinity();
        }
        for (int i=0;i<power2;i++) {
          RFHAT[i]=numeric_limits<double>::infinity();
        }
        break;
      }
      default:
        abort();
    }
# endif

    double dx=(xmax-xmin)/static_cast<double>(n);
    double x=xmin;
    double ymax=-numeric_limits<double>::max();
    double ymin=numeric_limits<double>::max();
    switch (code) {
      case FFTPACK: {
//      TRACER_CALL(t,"runMain FFTPACK");
        for (int j=0;j<n;j++,x+=dx) df[j]=fcn(x);
#ifdef DEBUG
//      for (int j=0;j<n;j++) {
//        cout << "f[" << j << "] = " << df[j] << endl;
//      }
#endif
        dffti_(n,dwsave); // prime factorization of n + fourier modes
        dfftf_(n,df,dwsave); // df replaced by fourier coefficients
#ifdef DEBUG
//      for (int j=0;j<n;j++) {
//        cout << "fhat[" << j << "] = " << df[j] << endl;
//      }
#endif
//      fftpack orders the fourier modes as follows:
//      cos(x*0),cos(x),sin(x),cos(x*2),sin(x*2),...,cos(x*n/2)
        double factor=1./static_cast<double>(n);
        for (int j=0;j<n;j++) {
          DFHAT[j]=df[j]*factor;
        }
        DFHAT[n-1]*=0.5;
        for (int j=n;j<power2;j++) DFHAT[j]=0.;
#ifdef DEBUG
//      for (int j=0;j<power2;j++) {
//        cout << "FHAT[" << j << "] = " << DFHAT[j] << endl;
//      }
#endif
        dfftb_(power2,DFHAT,dwsave2); // DFHAT replaced by inverse F.T.
        for (int j=0;j<power2;j++) {
          double y=DFHAT[j];
          ymax=max(ymax,y);
          ymin=min(ymin,y);
        }
#ifdef DEBUG
//      for (int j=0;j<power2;j++) {
//        cout << "F[" << j << "] = " << DFHAT[j] << endl;
//      }
#endif
        break;
      }
      case FFTWCOMPLEX: {
//      TRACER_CALL(t,"runMain FFTWCOMPLEX");
        for (int j=0;j<n;j++,x+=dx) {
          cf[j].re=fcn(x);
          cf[j].im=0.;
        }
        cplan_forward=
          fftw_create_plan(n,FFTW_FORWARD,FFTW_ESTIMATE);
        fftw_one(cplan_forward,cf,cfhat);
        fftw_real factor=1./static_cast<fftw_real>(n);
        for (int j=0;j<=n/2;j++) {
          CFHAT[j].re=cfhat[j].re*factor;
          CFHAT[j].im=cfhat[j].im*factor;
        }
        for (int j=n/2+1;j<=power2-n/2;j++) CFHAT[j].re=CFHAT[j].im=0.;
        for (int j=1;j<n/2;j++) {
          CFHAT[power2-j].re=cfhat[n-j].re*factor;
          CFHAT[power2-j].im=cfhat[n-j].im*factor;
        }
        fftw_one(cplan_backward,CFHAT,CF);
        for (int j=0;j<power2;j++) {
          double y=CF[j].re;
          ymax=max(ymax,y);
          ymin=min(ymin,y);
        }
        break;
      }
      case FFTWREAL: {
//      TRACER_CALL(t,"runMain FFTWREAL");
        for (int j=0;j<n;j++,x+=dx) rf[j]=fcn(x);
        rplan_forward=
          rfftw_create_plan(n,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);
        rfftw_one(rplan_forward,rf,rfhat);
        fftw_real factor=1./static_cast<fftw_real>(n);
//      fftw orders the fourier modes as follows:
//      cos(x*0),cos(x)/2,cos(x*2)/2,..,cos(x*[n/2-1])/2,cos(x*n/2),
//        -sin(x*[n/2-1])/2,...,-sin(x*2)/2,-sin(x)/2
        for (int j=0;j<n/2;j++) {
          RFHAT[j]=rfhat[j]*factor;
        }
        {
          int j=n/2;
          RFHAT[j]=rfhat[j]*factor*0.5;
        }
        for (int j=n/2+1;j<=power2-n/2;j++) RFHAT[j]=0;
        for (int j=1;j<n/2;j++) {
          RFHAT[power2-j]=rfhat[n-j]*factor;
        }
        rfftw_one(rplan_backward,RFHAT,RF);
        for (int j=0;j<power2;j++) {
          double y=RF[j];
          ymax=max(ymax,y);
          ymin=min(ymin,y);
        }
        break;
      }
      default:
        abort();
    }

    dx=(xmax-xmin)/static_cast<double>(power2);
    x=xmin;
    gt.newPage();
    gt.rescale(xmin,xmax,ymin,ymax);
    gt.setfgColor("black");
    gt.drawAxes();
    gt.setfgColor("red");
    for (int j=0;j<power2;j++,x+=dx) {
      double y=fcn(x);
      if (j==0) gt.movePen(x,y);
      else gt.drawLine(x,y);
    }
    double sum_err=0.;
//  gt.setfgColor(n-nmin+1,nmax-nmin);
    gt.setfgColor("blue");
    x=xmin;
    for (int j=0;j<power2;j++,x+=dx) {
      double y;
      switch (code) {
        case FFTPACK:
          y=DFHAT[j];
          break;
        case FFTWCOMPLEX:
          y=CF[j].re;
          break;
        case FFTWREAL:
          y=RF[j];
          break;
        default:
          abort();
      }
      if (j==0) gt.movePen(x,y);
      else gt.drawLine(x,y);
      sum_err+=abs(y-fcn(x));
//    cout << t << " " << z << endl;
    }
    sum_err=max(sum_err/power2,numeric_limits<double>::epsilon());
    gt.setfgColor("green");
    x=xmin;
    dx=(xmax-xmin)/static_cast<double>(n);
    for (int j=0;j<n;j++,x+=dx) {
      double y=fcn(x);
      gt.drawPlus(x,y,dx*0.25);
    }
    errors[nm]=log(sum_err)/log10;
    errors_min=min(errors_min,errors[nm]);
    errors_max=max(errors_max,errors[nm]);
    gt.flush();
    cout << "errors[" << n << "] = " << sum_err << endl;
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
    if (code==FFTWCOMPLEX) fftw_destroy_plan(cplan_forward);
    if (code==FFTWREAL) fftw_destroy_plan(rplan_forward);
  }

  XYGraphTool gt2("log(L1 error) vs log(n)","log_10(n)","log_10(error)",
    static_cast<double>(nmin)*log(2.)/log10,
    static_cast<double>(nmax)*log(2.)/log10,errors_min,errors_max,
    &cmap, 0,0.5);
  gt2.setbgColor("white");
  gt2.setfgColor("black");
  gt2.drawAxes();
  gt2.setfgColor("blue");
  double t=static_cast<double>(nmin)*log(2.)/log10;
  gt2.movePen(t,errors[nmin]);
  for (int n=nmin+1;n<=nmax;n++) {
    t=static_cast<double>(n)*log(2.)/log10;
    gt2.drawLine(t,errors[n]);
  }
  gt2.flush();
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

  if (df!=0) delete [] df;
  if (DFHAT!=0) delete [] DFHAT;
  if (dwsave!=0) delete [] dwsave;
  if (dwsave2!=0) delete [] dwsave2;
  if (cf!=0) delete [] cf;
  if (cfhat!=0) delete [] cfhat;
  if (CF!=0) delete [] CF;
  if (CFHAT!=0) delete [] CFHAT;
  if (rf!=0) delete [] rf;
  if (rfhat!=0) delete [] rfhat;
  if (RF!=0) delete [] RF;
  if (RFHAT!=0) delete [] RFHAT;
  if (code==FFTWCOMPLEX) fftw_destroy_plan(cplan_backward);
  if (code==FFTWREAL) fftw_destroy_plan(rplan_backward);
  delete [] errors;
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
