#include <float.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
#include <iostream>
#include <math.h>
#include <unistd.h>

#include "MemoryDebugger.H"
#include "Palette.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XYGraphTool.H"

#define npts 100

inline double f(double x) { return -10.*pow(1.-x,9)-2.*x; }
double A=-10.,B=10.;

//inline double f(double x) { return x*x-4.*sin(x); }
double f_gsl(double x,void *p) { return f(x); }
//double A=-4.,B=4.;

//inline double f(double x) { return atan(x); }
//double A=-4.,B=4.;

//inline double f(double x) { return exp(x)*cos(x); }
//double A=-4.,B=4.;

void computePlotPoints(double a,double b,double &fmin,double &fmax,
double *xarray,double *farray) {
  double dx=(b-a)*0.0625;
  a-=dx; b+=dx;

  dx=(b-a)/double(npts);
  xarray[0]=a; farray[0]=f(a);
  fmax=fmin=f(a);
  for (int i=1;i<=npts;i++) {
    xarray[i]=xarray[i-1]+dx; farray[i]=f(xarray[i]);
    fmin=min(fmin,farray[i]);
    fmax=max(fmax,farray[i]);
  }
}

void plotPoints(const double *xarray,const double *farray,
XYGraphTool &gt) {
  gt.movePen(xarray[0],farray[0]);
  for (int i=1;i<=npts;i++) {
    gt.drawLine(xarray[i],farray[i]);
  }
}

int main(int argc,char *argv[]) {
  cout << boolalpha;
  { MemoryDebugger md(1);
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

    int nsteps = numeric_limits<double>::digits;
    double *logf = OPERATOR_NEW_BRACKET(double,nsteps);
    double logfmin=numeric_limits<double>::max();
    double logfmax=-numeric_limits<double>::max();

    double x=B;
    if (A>B) { x=A; A=B; B=x; }
    cout << "a,b = " << A << " " << B << endl;

    double dx=(B-A)*0.0078125;
    double xarray[npts+1],farray[npts+1];
    double fmin,fmax;
    computePlotPoints(A,B,fmin,fmax,xarray,farray);

//  setup interactive graphics
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("bisection","x","f",A,B,fmin,fmax,&cmap,0,winsize);

//  initialize graphics display
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();
//  draw function
    gt.setfgColor("blue");
    plotPoints(xarray,farray,gt);
//  force X to perform the requests
    gt.flush();

    double a,b,y;
    cout << "click mouse in window to choose starting point" << endl;
    int button=gt.getMouse(a,y);
    gt.setfgColor("red");
    gt.movePen(a,gt.getLowY());
    gt.drawLine(a,gt.getHighY());
    cout << "click mouse in window again to choose second starting point"
         << endl;
    button=gt.getMouse(b,y);
    gt.movePen(b,gt.getLowY());
    gt.drawLine(b,gt.getHighY());
    if (a>b) { x=a; a=b; b=x; }
    double fa=f(a),fb=f(b);

//  GSL iteration
    gsl_root_fsolver *solver_gsl
      = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
    gsl_function fcn_struct_gsl;
      fcn_struct_gsl.function = &f_gsl;
      fcn_struct_gsl.params = 0;
    int status=gsl_root_fsolver_set(solver_gsl,&fcn_struct_gsl,a,b);
    int step=0;
    double epsabs=0.;
    double epsrel=numeric_limits<double>::epsilon();
    do {
      step++;
      status = gsl_root_fsolver_iterate(solver_gsl);
      double r = gsl_root_fsolver_root(solver_gsl);
      double x_lo = gsl_root_fsolver_x_lower(solver_gsl);
      double x_hi = gsl_root_fsolver_x_upper(solver_gsl);
      status = gsl_root_test_interval(x_lo, x_hi, epsabs, epsrel);
//    cout << "\tstep,x,x_lo,x_hi = " << step << " " << r << " " 
//         << x_lo << " " << x_hi << endl;
    } while (status == GSL_CONTINUE && step < nsteps);
    cout << "GSL root = " << gsl_root_fsolver_root(solver_gsl) 
         << " found in " << step << " iterations" << endl;
    gsl_root_fsolver_free(solver_gsl); solver_gsl = 0;

//  loop over iterations
    for (int i=0;i<nsteps;i++) {
      x=0.5*(a+b);
      double fx=f(x);
      cout << "step = " << i << " f(" << x << ") = " << fx << endl;

      logf[i]=log10(max(numeric_limits<double>::epsilon(),abs(fx)));
      logfmin=min(logfmin,logf[i]);
      logfmax=max(logfmax,logf[i]);

      gt.newPage();
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      plotPoints(xarray,farray,gt);

//    mark points
      gt.setfgColor("red");
      gt.movePen(a,gt.getLowY());
      gt.drawLine(a,gt.getHighY());
      gt.movePen(b,gt.getLowY());
      gt.drawLine(b,gt.getHighY());
      gt.setfgColor("green");
      gt.drawCross(x,fx,dx);
      gt.flush();

      XYGraphTool::WINDOW_TYPE::QuitButton qb;

      if (a>=x || x>=b) {
        nsteps=i+1;
        break;
      }
      if (fx*fa>0.) { a=x; fa=fx; }
      else { b=x; fb=fx; }
    }
    cout << setprecision(16);
    cout << "final solution = " << x << endl;
    XYGraphTool gt2("error","step","log_10(f(x_i))",
      0.,static_cast<double>(nsteps-1),logfmin,logfmax,&cmap,0,winsize);
    gt2.newPage();
    gt2.setbgColor("white");
    gt2.setfgColor("black");
    gt2.drawAxes();
    gt2.setfgColor("blue");
    gt2.movePen(0.,logf[0]);
    for (int i=1;i<nsteps;i++) {
      gt2.drawLine(static_cast<double>(i),logf[i]);
    }
    gt2.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    delete [] logf;
  } // pal,cmap,gt go out of scope here 
  return EXIT_SUCCESS;
}
