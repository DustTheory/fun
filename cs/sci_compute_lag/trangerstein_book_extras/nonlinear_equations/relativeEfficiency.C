#include <float.h>
#include <fstream>
#include <limits>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
#include <iostream>
#include <math.h>
#include <sys/times.h>
#include <sys/param.h>
#include <unistd.h>

#include "MemoryDebugger.H"
#include "Palette.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XYGraphTool.H"

#define npts 100

inline double f(double x,double &scale) { 
  double term1=x*x;
  double term2=4.*sin(x);
  scale=2.*max(abs(term1),abs(term2));
  return term1-term2;
}
inline double fp(double x) { return 2.*x-4.*cos(x); }
inline double fmodified(double x,double &scale) {
  double term1=0.25*x;
  double term2=(abs(x)==0. ? 1. : sin(x)/x);
  scale=2.*max(abs(term1),abs(term2));
  return term1-term2;
}
double A=-4.,B=4.;

//inline double f(double x,double &scale) { 
//  double term=atan(x);
//  scale=abs(term);
//  return term;
//}
//inline double fp(double x) { return 1./(1.+x*x); }
//double A=-4.,B=4.;

//inline double f(double x,double &scale) { 
//  double term=exp(x)*cos(x);
//  scale=abs(term);
//  return term;
//}
//inline double fp(double x) { return exp(x)*(cos(x)-sin(x)); }
//double A=-4.,B=4.;

int main(int argc,char *argv[]) {
  cout << boolalpha;
  MemoryDebugger md(1);
  {
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    double *xarray=OPERATOR_NEW_BRACKET(double,npts+1);
    double *farray=OPERATOR_NEW_BRACKET(double,npts+1);
    double *fmarray=OPERATOR_NEW_BRACKET(double,npts+1);
    double scale;
    double dx=(B-A)/double(npts);
    xarray[0]=A; farray[0]=f(A,scale); fmarray[0]=fmodified(A,scale);
    double fmax=max(farray[0],fmarray[0]);
    double fmin=min(farray[0],fmarray[0]);
    for (int i=1;i<=npts;i++) {
      xarray[i]=xarray[i-1]+dx;
      farray[i]=f(xarray[i],scale);
      fmarray[i]=fmodified(xarray[i],scale);
      fmin=min(fmin,min(farray[i],fmarray[i]));
      fmax=max(fmax,max(farray[i],fmarray[i]));
    }

    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    {
      XYGraphTool gtf("functions","x","f",A,B,fmin,fmax,&cmap,0,winsize);
      gtf.newPage();
      gtf.setbgColor("white");
      gtf.setfgColor("black");
      gtf.drawAxes();
      gtf.setfgColor("blue");
      gtf.movePen(xarray[0],farray[0]);
      for (int i=1;i<=npts;i++) gtf.drawLine(xarray[i],farray[i]);
      gtf.setfgColor("green");
      gtf.movePen(xarray[0],fmarray[0]);
      for (int i=1;i<=npts;i++) gtf.drawLine(xarray[i],fmarray[i]);
      gtf.flush();
      cout << "\toriginal function in blue, modified function in green"
           << endl;
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }

    int nsteps=48; // number of bits in mantissa for double
    int function_count = 1000000; // make functions take measureable time
    double rt_epsilon=sqrt(numeric_limits<double>::epsilon());
    clock_t start;
    struct tms usage;

//  Newton's method:
    double *logx_error_newton=OPERATOR_NEW_BRACKET(double,nsteps+1);
    double *logf_newton=OPERATOR_NEW_BRACKET(double,nsteps+1);
    double *elapsed_newton=OPERATOR_NEW_BRACKET(double,nsteps+1);
    double x=1.2,fx=f(x,scale),fpx=fp(x),dxold=-fx/fpx;
    elapsed_newton[0]=0.;
    xarray[0]=x;
    farray[0]=fx;
    int nsteps_newton=nsteps;

    times(&usage);
    start=usage.tms_utime;
    for (int i=1;i<=nsteps_newton;i++) {
      double dx=-fx/fpx;
      if (abs(dx)<=rt_epsilon*(abs(x)+scale/abs(fpx)) &&
      abs(dx)>=abs(dxold)) { // maximum attainable accuracy
        nsteps_newton=i-1;
        break; 
      }
      x=x+dx;
      dxold=dx;
      for (int k=0;k<function_count;k++) {
        fx=f(x,scale);
        fpx=fp(x);
      }
      xarray[i]=x;
      farray[i]=fx;
      times(&usage);
      elapsed_newton[i]=1.e-6*static_cast<double>(usage.tms_utime-start)
                /static_cast<double>(sysconf(_SC_CLK_TCK));
    }

//  secant method:
    int nsteps_secant=nsteps+1;
    double *xarray_secant=OPERATOR_NEW_BRACKET(double,nsteps_secant+1);
    double *farray_secant=OPERATOR_NEW_BRACKET(double,nsteps_secant+1);
    double *logx_error_secant=OPERATOR_NEW_BRACKET(double,nsteps_secant+1);
    double *logf_secant=OPERATOR_NEW_BRACKET(double,nsteps_secant+1);
    double *elapsed_secant=OPERATOR_NEW_BRACKET(double,nsteps_secant+1);
    double xold=1.2,fxold=f(xold,scale);
    elapsed_secant[0]=0.;
    xarray_secant[0]=xold;
    farray_secant[0]=fxold;
    times(&usage);
    start=usage.tms_utime;
    x=2.6;
    dxold=x-xold;
    for (int k=0;k<function_count;k++) fx=f(x,scale);
    times(&usage);
    elapsed_secant[1]=1.e-6*static_cast<double>(usage.tms_utime-start)
              /static_cast<double>(sysconf(_SC_CLK_TCK));
    xarray_secant[1]=x;
    farray_secant[1]=fx;

    for (int i=2;i<=nsteps_secant;i++) {
      if (x==xold || fx==fxold) { // prevent no change or divide by zero
        nsteps_secant=i-1;
        break;
      }
      double dx=-fx*(x-xold)/(fx-fxold);
      double xnew=x+dx;
      if (abs(fx)<rt_epsilon*scale && abs(dx)>=abs(dxold)) {
//      max attainable accuracy
        nsteps_secant=i-1;
        break;
      }
      xold=x;
      fxold=fx;
      dxold=dx;
      x=xnew;
      for (int k=0;k<function_count;k++) fx=f(x,scale);
      xarray_secant[i]=x;
      farray_secant[i]=fx;
      times(&usage);
      elapsed_secant[i]=1.e-6*static_cast<double>(usage.tms_utime-start)
                  /static_cast<double>(sysconf(_SC_CLK_TCK));
    }

//  bisection:
    int nsteps_bisection=nsteps+1;
    double *xarray_bisection=
      OPERATOR_NEW_BRACKET(double,nsteps_bisection+1);
    double *farray_bisection=
      OPERATOR_NEW_BRACKET(double,nsteps_bisection+1);
    double *logx_error_bisection=
      OPERATOR_NEW_BRACKET(double,nsteps_bisection+1);
    double *logf_bisection=OPERATOR_NEW_BRACKET(double,nsteps_bisection+1);
    double *elapsed_bisection=
      OPERATOR_NEW_BRACKET(double,nsteps_bisection+1);
    xold=1.2;
    fxold=f(xold,scale);
    elapsed_bisection[0]=0.;
    xarray_bisection[0]=xold;
    farray_bisection[0]=fxold;
    times(&usage);
    start=usage.tms_utime;
    x=2.6;
    for (int k=0;k<function_count;k++) fx=f(x,scale);
    times(&usage);
    elapsed_bisection[1]=1.e-6*static_cast<double>(usage.tms_utime-start)
              /static_cast<double>(sysconf(_SC_CLK_TCK));
    xarray_bisection[1]=x;
    farray_bisection[1]=fx;

    for (int i=2;i<=nsteps_bisection;i++) {
      double xnew=0.5*(x+xold);
      if (xnew<=xold || xnew>=x) {
//      max attainable accuracy
        nsteps_bisection=i-1;
        break;
      }
      double fxnew;
      for (int k=0;k<function_count;k++) fxnew=f(xnew,scale);
      if (fxnew*fxold>=0.) {
        xold=xnew;
        fxold=fxnew;
      } else {
        x=xnew;
        fx=fxnew;
      }
      xarray_bisection[i]=xnew;
      farray_bisection[i]=fxnew;
      times(&usage);
      elapsed_bisection[i]=1.e-6*static_cast<double>(usage.tms_utime-start)
                  /static_cast<double>(sysconf(_SC_CLK_TCK));
    }

    double logfmin=numeric_limits<double>::max();
    double logfmax=-numeric_limits<double>::max();
    double logxmin=numeric_limits<double>::max();
    double logxmax=-numeric_limits<double>::max();
    for (int i=0;i<=nsteps;i++) {
      logf_newton[i]=log10(max(numeric_limits<double>::epsilon(),
                               abs(farray[i])));
      logfmin=min(logfmin,logf_newton[i]);
      logfmax=max(logfmax,logf_newton[i]);

      logx_error_newton[i]=log10(max(numeric_limits<double>::epsilon(),
                              abs(xarray[i]-x)));
      logxmin=min(logxmin,logx_error_newton[i]);
      logxmax=max(logxmax,logx_error_newton[i]);
    }
    for (int i=0;i<=nsteps_secant;i++) {
      logf_secant[i]=log10(max(numeric_limits<double>::epsilon(),
                         abs(farray_secant[i])));
      logfmin=min(logfmin,logf_secant[i]);
      logfmax=max(logfmax,logf_secant[i]);

      logx_error_secant[i]=log10(max(numeric_limits<double>::epsilon(),
        abs(xarray_secant[i]-x)));
      logxmin=min(logxmin,logx_error_secant[i]);
      logxmax=max(logxmax,logx_error_secant[i]);
    }
    for (int i=0;i<=nsteps_bisection;i++) {
      logf_bisection[i]=log10(max(numeric_limits<double>::epsilon(),
                         abs(farray_bisection[i])));
      logfmin=min(logfmin,logf_bisection[i]);
      logfmax=max(logfmax,logf_bisection[i]);

      logx_error_bisection[i]=log10(max(numeric_limits<double>::epsilon(),
        abs(xarray_bisection[i]-x)));
      logxmin=min(logxmin,logx_error_bisection[i]);
      logxmax=max(logxmax,logx_error_bisection[i]);
    }

    double tmax=max(elapsed_newton[nsteps_newton],
      max(elapsed_secant[nsteps_secant],
          elapsed_bisection[nsteps_bisection]));
    XYGraphTool gt("log f vs time","time","log_10 ( | f | )",
      0.,tmax,logfmin,logfmax,&cmap,0,winsize);
    XYGraphTool gt2("log (x-z) vs time","time","log_10 ( | x-z | )",
      0.,tmax,logfmin,logfmax,&cmap,0,winsize);

    {
      gt.newPage();
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      gt.movePen(elapsed_newton[0],logf_newton[0]);
      for (int i=1;i<=nsteps_newton;i++) {
        gt.drawLine(elapsed_newton[i],logf_newton[i]);
      }
      gt.setfgColor("green");
      gt.movePen(elapsed_secant[0],logf_secant[0]);
      for (int i=1;i<=nsteps_secant;i++) {
        gt.drawLine(elapsed_secant[i],logf_secant[i]);
      }
      gt.flush();
      gt.setfgColor("red");
      gt.movePen(elapsed_bisection[0],logf_bisection[0]);
      for (int i=1;i<=nsteps_bisection;i++) {
        gt.drawLine(elapsed_bisection[i],logf_bisection[i]);
      }
      gt.flush();

      gt2.newPage();
      gt2.setbgColor("white");
      gt2.setfgColor("black");
      gt2.drawAxes();
      gt2.setfgColor("blue");
      gt2.movePen(elapsed_newton[0],logx_error_newton[0]);
      for (int i=1;i<=nsteps_newton;i++) {
        gt2.drawLine(elapsed_newton[i],logx_error_newton[i]);
      }
      gt2.setfgColor("green");
      gt2.movePen(elapsed_secant[0],logx_error_secant[0]);
      for (int i=1;i<=nsteps_secant;i++) {
        gt2.drawLine(elapsed_secant[i],logx_error_secant[i]);
      }
      gt2.setfgColor("red");
      gt2.movePen(elapsed_bisection[0],logx_error_bisection[0]);
      for (int i=1;i<=nsteps_bisection;i++) {
        gt2.drawLine(elapsed_bisection[i],logx_error_bisection[i]);
      }
      gt2.flush();
      cout << "\n\tNewton took " << nsteps_newton << " steps" << endl;
      cout << "\tsecant took " << nsteps_secant << " steps" << endl;
      cout << "\tbisection took " << nsteps_bisection << " steps" << endl;
      cout << "\tNewton in blue, secant in green, bisection in red" << endl;

      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }

    delete [] xarray; xarray=0;
    delete [] farray; farray=0;
    delete [] fmarray; fmarray=0;
    delete [] logx_error_newton; logx_error_newton=0;
    delete [] logf_newton; logf_newton=0;
    delete [] elapsed_newton; elapsed_newton=0;
    delete [] xarray_secant; xarray_secant=0;
    delete [] farray_secant; farray_secant=0;
    delete [] logx_error_secant; logx_error_secant=0;
    delete [] logf_secant; logf_secant=0;
    delete [] elapsed_secant; elapsed_secant=0;
    delete [] xarray_bisection; xarray_bisection=0;
    delete [] farray_bisection; farray_bisection=0;
    delete [] logx_error_bisection; logx_error_bisection=0;
    delete [] logf_bisection; logf_bisection=0;
    delete [] elapsed_bisection; elapsed_bisection=0;
  } // pal,cmap,gt go out of scope here 
  return EXIT_SUCCESS;
}
