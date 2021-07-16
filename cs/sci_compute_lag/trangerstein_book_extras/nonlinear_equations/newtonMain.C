#include <float.h>
#include <fstream>
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

inline double f(double x) { return x*x-4.*sin(x); }
inline double fscale(double x) { return max(x*x,4.*sin(x)); }
inline double fp(double x) { return 2.*x-4.*cos(x); }
double f_gsl(double x,void *p) { return f(x); }
double fp_gsl(double x,void *p) { return fp(x); }
void fdf_gsl(double x,void *p,double *y,double *dy) {
  *y=f(x);
  *dy=fp(x);
}
double A=-4.,B=4.;

//inline double f(double x) { return atan(x); }
//inline double fp(double x) { return 1./(1.+x*x); }
//double A=-4.,B=4.;

//inline double f(double x) { return exp(x)*cos(x); }
//inline double fp(double x) { return exp(x)*(cos(x)-sin(x)); }
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

void plotTangentLine(double a,double b,double x,double fx,double fpx,
XYGraphTool &gt) {
//TRACER_CALL(t,"plotTangentLine")
//cout << "\tx = " << x << endl;
//cout << "\tfx = " << fx << endl;
//cout << "\tfpx = " << fpx << endl;
  double xlo=gt.getLowX();
  double xhi=gt.getHighX();
  double ylo=gt.getLowY();
  double yhi=gt.getHighY();
  gt.movePen(x,fx);
  double y=fx+fpx*(xlo-x);
  if (y<ylo) {
    gt.drawLine((ylo-fx)/fpx+x,ylo);
  } else if (y>yhi) {
    gt.drawLine((yhi-fx)/fpx+x,yhi);
  } else {
    gt.drawLine(xlo,y);
  }
  gt.movePen(x,fx);
  y=fx+fpx*(xhi-x);
  if (y<ylo) {
    gt.drawLine((ylo-fx)/fpx+x,ylo);
  } else if (y>yhi) {
    gt.drawLine((yhi-fx)/fpx+x,yhi);
  } else {
    gt.drawLine(xhi,y);
  }
}

int main(int argc,char *argv[]) {
  cout << boolalpha;
  { MemoryDebugger md(1);
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

    int nsteps=7;
    double *logf = OPERATOR_NEW_BRACKET(double,nsteps);
    double logfmin=numeric_limits<double>::max();
    double logfmax=-numeric_limits<double>::max();

    double x;
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
    XYGraphTool gt("newton","x","f",A,B,fmin,fmax,&cmap,0,winsize);

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

    double y;
    cout << "click mouse in window to choose starting guess" <<endl;
    int button=gt.getMouse(x,y);

//  GSL iteration
    gsl_root_fdfsolver *solver_gsl
      = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
    gsl_function_fdf fcn_struct_gsl;
      fcn_struct_gsl.f = &f_gsl; 
      fcn_struct_gsl.df = &fp_gsl; 
      fcn_struct_gsl.fdf = &fdf_gsl; 
      fcn_struct_gsl.params = 0;
    int status=gsl_root_fdfsolver_set(solver_gsl,&fcn_struct_gsl,x);
    int step=0;
    double epsabs=0.;
    double epsrel=numeric_limits<double>::epsilon();
    do {
      step++;
      status=gsl_root_fdfsolver_iterate(solver_gsl);
      double xold=x;
      double xnew=gsl_root_fdfsolver_root(solver_gsl);
      status=gsl_root_test_delta(xnew,xold,epsabs,epsrel);
//    cout << "\tstep,xnew,xold = " << step << " " << xnew << " " 
//         << xold << endl;
    } while (status == GSL_CONTINUE && step < nsteps);
    double z=gsl_root_fdfsolver_root(solver_gsl);
    cout << "GSL root = " << z << " found in " << step << " iterations"
         << endl;
    gsl_root_fdfsolver_free(solver_gsl); solver_gsl = 0;

//  loop over iterations
    for (int i=0;i<nsteps;i++) {
      double fx=f(x), fpx=fp(x);
//    cout << "f(" << x << ") = " << fx << ", f'(x) = " << fpx << endl;

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
      plotTangentLine(xarray[0],xarray[npts],x,fx,fpx,gt);
      gt.setfgColor("green");
      gt.drawCross(x,fx,dx);
      gt.flush();

      XYGraphTool::WINDOW_TYPE::QuitButton qb;

      double s=fx/fpx;
      cout << "i,s,f,error = " << i << " " << s/abs(x) << " "
           << fx/fscale(x) << " " << abs(x-z)/abs(z) << endl;
      x=x-s;
    }
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
