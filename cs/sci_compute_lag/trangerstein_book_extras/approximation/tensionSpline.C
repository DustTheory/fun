#include <iostream>
//#include <math.h> // for HUGE_VAL,M_PI
#include <stdlib.h>

using namespace std;

#include "Debug.H"
#include "GUIEnumInputParameter.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
//#include "XColormap.H"
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
double tension=1.;

double jump0=sqrt(0.2);
double jump1=0.5;

char* display_name=0;
bool skip_gui=false;
int n=2;
int nplot=1000;
GUI_INPUT_PARAMETER_LIST_TYPE *param_list=0;
Palette *pal=0;
XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE *cmap=0;
XYGraphTool *gt=0;
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
    case smooth: {
//    return sin(M_PI_2*t);
//    return sinh(tension*t);
      double x=2.*t-1.;
      return 1./(1.+25.*x*x);
    }
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

void setField(char* field,int field_width,int val) {
  field[field_width]='\0';
  int i=field_width-1;
  int n=val;
  for (;i>=0;i--,n/=10) field[i]='0'+n %10;
}

void checkMainInput() { 
}

void runMain(bool /*called_before*/) {
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
    gt=OPERATOR_NEW XYGraphTool("tension spline","x","y",
      xmin,xmax,ymin,ymax,cmap,0,0.5);
  }
  gt->setbgColor("white");

  double *alpha=OPERATOR_NEW_BRACKET(double,n);
  double *beta=OPERATOR_NEW_BRACKET(double,n);
  double *gamma=OPERATOR_NEW_BRACKET(double,n);
  double *h=OPERATOR_NEW_BRACKET(double,n);

  double *diag=OPERATOR_NEW_BRACKET(double,n-1);
  double *subdiag=OPERATOR_NEW_BRACKET(double,n-1);

  double *x=OPERATOR_NEW_BRACKET(double,n+1);
  double *y=OPERATOR_NEW_BRACKET(double,n+1);
  double *z=OPERATOR_NEW_BRACKET(double,n+1);

//double log10=log(10.);
#ifdef INDEF
  for (int i=0;i<n-1;i++) {
    diag[i]=subdiag[i]=HUGE_VAL;
  }
  for (int i=0;i<n;i++) {
    alpha[i]=beta[i]=gamma[i]=h[i]=HUGE_VAL;
  }
  for (int i=0;i<=n;i++) {
    x[i]=y[i]=z[i]=HUGE_VAL;
  }
//for (int i=0;i<=n;i++) {
//  coef[i]=diag[i]=dif[i]=subdiag[i]=HUGE_VAL;
//}
#endif

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
    y[j]=fcn(x[j]);
  }
  x[n]=xmax;
  y[n]=fcn(x[n]);
#ifdef DEBUG
//for (int j=0;j<=n;j++) {
//  cout << "\tx,y[" << n << "] = " << x[n] << " " << y[n] << endl;
//}
#endif

  double tension2=tension*tension;
  for (int j=0;j<n;j++) {
    h[j]=x[j+1]-x[j];
    double th=tension*h[j];
    double thoversinh=th/sinh(th);
    double hinv=1./h[j];
    alpha[j]=(1.-thoversinh)*hinv;
    beta[j]=(cosh(th)*thoversinh-1.)*hinv;
    gamma[j]=tension2*(y[j+1]-y[j])*hinv;
#ifdef DEBUG
//  cout << "\talpha,beta,gamma[" << j << "] = " << alpha[j] << " " 
//       << beta[j] << " " << gamma[j] << endl;
#endif
  }

  for (int j=0;j<n-1;j++) {
    subdiag[j]=alpha[j];
    diag[j]=beta[j]+beta[j+1];
#ifdef DEBUG
//  cout << "\tdiag,subdiag[" << j << "] = " << diag[j] << " " 
//       << subdiag[j] << endl;
#endif
  }
  for (int j=1;j<n;j++) {
    z[j]=gamma[j]-gamma[j-1];
#ifdef DEBUG
//  cout << "\tz[" << j << "] = " << z[j] << endl;
#endif
  }

  int nm1=n-1;
  int nrhs=1;
  int info=INT_MAX;
  F77NAME(dptsv)(nm1,nrhs,diag,subdiag+1,z+1,n,info);
  z[0]=0.;
  z[n]=0.;
#ifdef DEBUG
//  for (int j=0;j<=n;j++) {
//    cout << "\tcoef[" << j << "] = " << coef[j] << endl;
//  }
#endif

  double ylo=HUGE_VAL;
  double yhi=-HUGE_VAL;
  dx=(xmax-xmin)/static_cast<double>(nplot);
  double t=xmin;
  double t2=tension*tension;
  for (int j=0;j<=nplot;j++,t+=dx) {
    double s=fcn(t);
    ylo=min(ylo,s);
    yhi=max(yhi,s);

    int i=bin_search(n,0,t,x);
    double dtm=t-x[i];
    double dtp=x[i+1]-t;
    s=(z[i]*sinh(tension*dtp)
     +z[i+1]*sinh(tension*dtm))
     /(t2*sinh(tension*h[i]))
     +(y[i]-z[i]/t2)*dtp/h[i]
     +(y[i+1]-z[i+1]/t2)*dtm/h[i];
    ylo=min(ylo,s);
    yhi=max(yhi,s);
  }

  t=xmin;
  dx=(xmax-xmin)/static_cast<double>(nplot);
  gt->newPage();
  gt->rescale(xmin,xmax,ylo,yhi);
  gt->setfgColor("black");
  gt->drawAxes();
  gt->setfgColor("blue");
  gt->movePen(t,fcn(t));
  for (int j=0;j<=nplot;j++,t+=dx) {
    gt->drawLine(t,fcn(t));
  }
  gt->setfgColor("red");
  t=xmin;
  double maxerr=0.;
  for (int j=0;j<=nplot;j++,t+=dx) {
    int i=bin_search(n,0,t,x);
    double dtm=t-x[i];
    double dtp=x[i+1]-t;
    double s=(z[i]*sinh(tension*dtp)
             +z[i+1]*sinh(tension*dtm))
            /(t2*sinh(tension*h[i]))
            +(y[i]-z[i]/t2)*dtp/h[i]
            +(y[i+1]-z[i+1]/t2)*dtm/h[i];
#ifdef DEBUG
//  cout << "s,y,i[" << j << "] = " << s << " " << fcn(t) << " " << i
//       << endl;
#endif
    if (j==0) gt->movePen(t,s);
    else gt->drawLine(t,s);
  }
  dx=(xmax-xmin)/static_cast<double>(n);
  for (int j=0;j<=n;j++) {
    gt->drawPlus(x[j],y[j],0.5*dx);
  }
  gt->flush();
//setField(index_field,3,n);
//snprintf(filename,LENGTH_NAME,"%s_%s",filebase,index_field);
//gt->writeXPM(filename);
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

  delete [] alpha;
  delete [] beta;
  delete [] gamma;
  delete [] h;
  delete [] diag;
  delete [] subdiag;
  delete [] x;
  delete [] y;
  delete [] z;
}

void cleanup() { }

void shutdown() {
  if (gt) delete gt; gt=0;
  if (cmap) delete cmap; cmap=0;
  if (pal) delete pal; pal=0;
  while (param_list->notEmpty()) delete param_list->delAfter(0);
  delete param_list; param_list=0;
}

int main(int argc,char* argv[]) {
#ifdef DEBUG
//setTraps();
#endif
#ifdef MEM_DEBUG
  MemoryDebugger md(1);
#endif
  {
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

  param_list=OPERATOR_NEW GUI_INPUT_PARAMETER_LIST_TYPE("Param List");
  { const char *group="Problem Parameters";
    param_list->append(OPERATOR_NEW GUIInputParameter<int>(
      n,"interpolation points",3,INT_MAX,group));
    param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      itest_function,"test function","piecewise constant",
      "piecewise linear","piecewise quadratic","smooth") );
    param_list->append(OPERATOR_NEW GUIEnumInputParameter(group,
      iinterpolation_mesh,"mesh points","equidistant","random") );
    param_list->append(OPERATOR_NEW GUIInputParameter<double>(
      tension,"tension",0.,DBL_MAX,group));
    param_list->append(OPERATOR_NEW GUIInputParameter<double>(
      jump0,"discontinuity x location",0.,1.,group));
    param_list->append(OPERATOR_NEW GUIInputParameter<double>(
      jump1,"discontinuity y location",0.,1.,group));
    param_list->append(OPERATOR_NEW GUIInputParameter<int>(
      nplot,"plotting points",2,INT_MAX,group));
  }

  if (argc>1) {
      ifstream in_file;
      in_file.open(argv[1],ios::in);
      if (in_file) {
        istream &infile(in_file);
        infile.clear(ios::goodbit);
        infile.seekg(0,ios::beg);
        char name[LENGTH_NAME], comment[LENGTH_COMMENT];
        bool found_main=FALSE;
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
