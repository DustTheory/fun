#include <iostream>
#include <math.h> // for HUGE_VAL,M_PI
#include <stdlib.h>

using namespace std;

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "GUIInputs.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "TimedObject.H"
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

enum TestFunction {piecewise_constant,piecewise_linear,
  piecewise_quadratic,smooth};
enum InterpolationMesh {equally_spaced,randomly_spaced};
TestFunction test_function=smooth;
InterpolationMesh interpolation_mesh=equally_spaced;
int itest_function=test_function;
int iinterpolation_mesh=interpolation_mesh;

double jump0=HUGE_VAL;
double jump1=HUGE_VAL;

char* display_name=0;
bool skip_gui=false;
int nmin=2;
int nmax=2;
GUI_INPUT_PARAMETER_LIST_TYPE *param_list=0;
Palette *pal=0;
XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE *cmap=0;
XYGraphTool *gt=0;
XYGraphTool *gt2=0;
extern "C" {
  void F77NAME(dptsv)(const int &n,const int &nrhs,double *d,double *e,
    double *b,const int &ldb,int &info);
}

double fcn(double t) { 
  switch (test_function) {
    case piecewise_constant:
      return (t<jump0 ? 0. : 1.);
    case piecewise_linear:
      return (t<jump0 ? jump1*t/jump0 : 
              1.-(1.-jump1)*(1.-t)/(1.-jump0));
    case piecewise_quadratic: {
      double den=1./(jump0*(1.-jump0));
      double c=(jump0-jump1)*den;
      double bl=(jump1-jump0*jump0)*den;
      double br=(jump1+jump0*(-2.+jump0))*den;
      return (t<jump0 ? t*(bl+c*t) : 
              1.+(1.-t)*(br+c-c*t));
    }
    case smooth:
//    return sin(M_PI_2*t);
      return 1./(1.+25.*t*t); // runge example
    default:
      abort();
      return HUGE_VAL;
  }
}

int bin_search(int high,int low, double t,const double *x) {
  for (; high-low>1;) {
    int i=(high+low)/2;
    if (t<x[i]) high=i;
    else low=i;
  }
  return low;
}

void setField(char* field,INTEGER field_width,INTEGER val) {
  field[field_width]='\0';
  int i=field_width-1;
  int n=val;
  for (;i>=0;i--,n/=10) field[i]='0'+n %10;
}

void checkMainInput() {
  if (nmax<=nmin) nmax=nmin+1;
}

void runMain(BOOLEAN /*called_before*/) {
  test_function=static_cast<TestFunction>(itest_function);
  char filebase[LENGTH_NAME];
  char index_field[3];
  char filename[LENGTH_NAME];
  double xmin=0.;
  double xmax=1.;
  double ymin=0.;
  double ymax=1.;
  switch (test_function) {
    case piecewise_constant:
      snprintf(filebase,LENGTH_NAME,"piecewise_constant");
      break;
    case piecewise_linear:
      snprintf(filebase,LENGTH_NAME,"piecewise_linear");
      break;
    case piecewise_quadratic:
      snprintf(filebase,LENGTH_NAME,"piecewise_quadratic");
      break;
    case smooth:
      snprintf(filebase,LENGTH_NAME,"smooth");
      break;
    default:
      abort();
      break;
  }
  if (pal==0) pal=OPERATOR_NEW Palette();
  if (cmap==0) {
    cmap=OPERATOR_NEW XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE(pal);
  }
  if (gt==0) {
    gt=OPERATOR_NEW XYGraphTool("function and spline","x","y",
      xmin,xmax,ymin,ymax,cmap,0,0.5);
  }
  gt->setbgColor("white");

  double *a=OPERATOR_NEW_BRACKET(double,nmax);
  double *b=OPERATOR_NEW_BRACKET(double,nmax);
  double *c=OPERATOR_NEW_BRACKET(double,nmax);
  double *coef=OPERATOR_NEW_BRACKET(double,nmax);
  double *diag=OPERATOR_NEW_BRACKET(double,nmax);
  double *dif=OPERATOR_NEW_BRACKET(double,nmax);
  double *subdiag=OPERATOR_NEW_BRACKET(double,nmax);
  double *errors=OPERATOR_NEW_BRACKET(double,nmax);
  double *x=OPERATOR_NEW_BRACKET(double,nmax);
  double *y=OPERATOR_NEW_BRACKET(double,nmax);

  double errors_min=DBL_MAX;
  double errors_max=-DBL_MAX;
  double log10=log(10.);
  for (int i=0;i<nmax;i++) {
    errors[i]=HUGE_VAL;
  }
  for (int n=nmin;n<nmax;n++) {
#ifdef INDEF
    for (int i=0;i<=n;i++) {
      a[i]=b[i]=c[i]=coef[i]=diag[i]=dif[i]=subdiag[i]=x[i]=y[i]=
        HUGE_VAL;
    }
#endif

//  generate interpolation mesh and array of function values
    double dx=(xmax-xmin)/static_cast<double>(n);
    x[0]=xmin;
    y[0]=fcn(x[0]);
    for (int j=1;j<n;j++) {
      switch (interpolation_mesh) {
        case equally_spaced:
          x[j]=x[j-1]+dx;
          break;
        case randomly_spaced:
          x[j]=x[j-1]+drand48()*(1.-x[j-1]);
          break;
        default:
          abort();
      }
      coef[j]=y[j]=fcn(x[j]);
    }
    x[n]=xmax;
    y[n]=fcn(x[n]);

//  prepare right-hand side for smooth spline linear system
    for (int j=0;j<n;j++) {
      subdiag[j]=x[j+1]-x[j];
      dif[j]=(y[j+1]-y[j])/subdiag[j];
    }
    for (int j=0;j<n;j++) {
      diag[j]=2.*(subdiag[j]+subdiag[j+1]);
    }
    for (int j=1;j<n;j++) {
      coef[j]=6.*(dif[j]-dif[j-1]);
    }

//  call lapack to solve symmetric positive-definite linear system
    int nm1=n-1;
    int nrhs=1;
    int info=INT_MAX;
    F77NAME(dptsv)(nm1,nrhs,diag,subdiag+1,coef+1,n,info);
    coef[0]=coef[n]=0.;

    for (int j=0;j<n;j++) {
      double h=x[j+1]-x[j];
      a[j]=(coef[j+1]-coef[j])/(6.*h);
      b[j]=0.5*coef[j];
      c[j]=dif[j]-h*(coef[j+1]+2.*coef[j])/6.;
    }

//  find spline bounds before plotting
    double ylo=HUGE_VAL;
    double yhi=-HUGE_VAL;
    dx=(xmax-xmin)*0.001;
    double t=xmin;
    for (int j=0;j<=1000;j++,t+=dx) {
      double z=fcn(t);
      ylo=min(ylo,z);
      yhi=max(yhi,z);
      int i=bin_search(n,0,t,x);
      double dt=t-x[i];
      z=y[i]+dt*(c[i]+dt*(b[i]+dt*a[i]));
      ylo=min(ylo,z);
      yhi=max(yhi,z);
    }

//  plot function
    t=xmin;
    dx=(xmax-xmin)*0.001;
    gt->newPage();
    gt->rescale(xmin,xmax,ylo,yhi);
    gt->setfgColor("black");
    gt->drawAxes();
    gt->setfgColor("blue");
    gt->movePen(t,fcn(t));
    for (int j=0;j<=1000;j++,t+=dx) {
      gt->drawLine(t,fcn(t));
    }
//  plot spline interpolant
    gt->setfgColor("red");
    t=xmin;
    double maxerr=0.;
    for (int j=0;j<=1000;j++,t+=dx) {
      int i=bin_search(n,0,t,x);
      double dt=t-x[i];
      double s=y[i]+dt*(c[i]+dt*(b[i]+dt*a[i]));
      if (j==0) gt->movePen(t,s);
      else gt->drawLine(t,s);
      maxerr=max(maxerr,abs(s-fcn(t)));
    }
    errors[n]=log(maxerr)/log10;
    errors_min=min(errors_min,errors[n]);
    errors_max=max(errors_max,errors[n]);
    gt->flush();
//  setField(index_field,3,n);
//  snprintf(filename,LENGTH_NAME,"%s_%s",filebase,index_field);
//  gt->writeXPM(filename);
    cout << "errors[" << n << "] = " << maxerr << endl;
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

  }

//plot errors in mesh refinement study
  if (gt2==0) {
    gt2=OPERATOR_NEW XYGraphTool("log_10(error) vs log_10(n points)",
      "log_10(n points)","log_10(error)",
      log(static_cast<double>(nmin))/log10,
      log(static_cast<double>(nmax))/log10,errors_min,errors_max,cmap,
      0,0.5);
  }
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

  delete [] a;
  delete [] b;
  delete [] c;
  delete [] coef;
  delete [] diag;
  delete [] dif;
  delete [] subdiag;
  delete [] errors;
  delete [] x;
  delete [] y;
}

void cleanup() {
  if (gt) delete gt; gt=0;
  if (gt2) delete gt2; gt2=0;
  if (cmap) delete cmap; cmap=0;
  if (pal) delete pal; pal=0;
}

void shutdown() {
  while (param_list->notEmpty()) delete param_list->delAfter(0);
  delete param_list; param_list=0;
}

int main(int argc,char* argv[]) {
#ifdef DEBUG
  setTraps();
#endif
  {
#ifdef MEM_DEBUG
    MemoryDebugger md(1);
#endif
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

  jump0=drand48();
  jump1=drand48();
  param_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Param List");
  { char *group="Problem Parameters";
    param_list->append(OPERATOR_NEW GUIInputParameter<int>(
      nmin,"minimum number of interpolation points",2,INT_MAX,group));
    param_list->append(OPERATOR_NEW GUIInputParameter<int>(
      nmax,"maximum number of interpolation points",2,INT_MAX,group));
    param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      itest_function,"test function","piecewise constant",
      "piecewise linear","piecewise quadratic","smooth") );
    param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iinterpolation_mesh,"mesh points","equidistant","random") );
    param_list->append(OPERATOR_NEW GUIInputParameter<double>(
      jump0,"discontinuity x location",0.,1.,group));
    param_list->append(OPERATOR_NEW GUIInputParameter<double>(
      jump1,"discontinuity y location",0.,1.,group));
  }

  if (argc>1) {
      ifstream in_file;
      in_file.open(argv[1],ios::in);
      if (in_file) {
        istream &infile(in_file);
        infile.clear(ios::goodbit);
        infile.seekg(0,ios::beg);
        char name[LENGTH_NAME], comment[LENGTH_COMMENT];
        BOOLEAN found_main=FALSE;
        while ( in_file >> setw(LENGTH_NAME) >> name ) {
          if ( strcmp(name,"Main") == 0 ) {
            in_file.getline( comment, LENGTH_COMMENT);
            found_main=TRUE;
            while ( in_file >> setw(LENGTH_NAME) >> name ) {
              if ( strcmp(name,"end") == 0 ) break;
              else if (strcmp(name,"skip_gui") == 0) { 
                int iskip_gui; // type bool not read correctly
                in_file >> iskip_gui;
                skip_gui= iskip_gui!=0;
              }
              else param_list->formattedRead(in_file,name);
              in_file.getline( comment, LENGTH_COMMENT);
            }
          } else in_file.getline(comment,LENGTH_COMMENT);
        }
        if ( !found_main ) skip_gui=FALSE;
      }
      if (skip_gui) checkMainInput();
    } else skip_gui=FALSE;

    if (skip_gui) {
      runMain(FALSE);
      cleanup();
      shutdown();
    } else {
#ifdef USE_GTK
      GTKGUI gui(argv[0],display_name,param_list,&runMain,
        &checkMainInput,&cleanup,&shutdown,TRUE);
#else
      GUI gui(argv[0],display_name,param_list,&runMain,
        &checkMainInput,&cleanup,&shutdown,TRUE);
#endif
      gui.createFileMenu();
      gui.createViewMenu();
      gui.createHelpMenu();
      gui.createMainWindow(argc,argv);
      gui.eventLoop();
    }

  }
  return EXIT_SUCCESS;
}

template class InputParameter<bool>;
template class InputParameter<double>;
template class InputParameter<int>;
