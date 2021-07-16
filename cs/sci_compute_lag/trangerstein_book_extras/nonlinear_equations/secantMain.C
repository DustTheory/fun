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
double f_gsl(double x,void *p) { return f(x); }
struct secant_params {
  double x0;
};
double fp_gsl(double x,void *p) { 
  struct secant_params *sp=reinterpret_cast<secant_params*>(p);
  return (f(x)-f(sp->x0))/(x-sp->x0);
}
void fdf_gsl(double x,void *p,double *y,double *dy) {
  *y=f(x);
  struct secant_params *sp=reinterpret_cast<secant_params*>(p);
  *dy=(f(x)-f(sp->x0))/(x-sp->x0);
}
double A=-4.,B=4.;

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

void plotSecantLine(double x1,double f1,double x2,double f2,
XYGraphTool &gt ) {
  if ( abs( f2 - f1 ) > 0. ) {
    double fpx = ( f2 - f1 ) / ( x2 - x1 );
    double secantY_at_xlo = f1 - fpx * ( x1 - gt.getLowX() );
    double secantY_at_xhi = f1 + fpx * ( gt.getHighX() - x1 );
    double secantX_at_ylo = x1 - ( f1 - gt.getLowY() ) / fpx;
    double secantX_at_yhi = x1 + ( gt.getHighY() - f1 ) / fpx;
    if ( gt.getLowY() <= secantY_at_xlo &&
    secantY_at_xlo <= gt.getHighY() ) {
      if ( gt.getLowX() <= secantX_at_ylo &&
      secantX_at_ylo <= gt.getHighX() ) {
        gt.movePen( gt.getLowX(), secantY_at_xlo );
        gt.drawLine( secantX_at_ylo, gt.getLowY() );
      } else {
        if ( gt.getLowY() <= secantY_at_xhi &&
        secantY_at_xhi <= gt.getHighY() ) {
          gt.movePen( gt.getLowX(), secantY_at_xlo );
          gt.drawLine( gt.getHighX(), secantY_at_xhi );
        } else if ( gt.getLowX() <= secantX_at_yhi &&
        secantX_at_yhi <= gt.getHighX() ) {
          gt.movePen( gt.getLowX(), secantY_at_xlo );
          gt.drawLine( secantX_at_yhi, gt.getHighY() );
        }
      }
    } else if ( gt.getLowX() <= secantX_at_ylo &&
    secantX_at_ylo <= gt.getHighX() ) {
      if ( gt.getLowY() <= secantY_at_xhi &&
      secantY_at_xhi <= gt.getHighY() ) {
        gt.movePen( secantX_at_ylo, gt.getLowY() );
        gt.drawLine( gt.getHighX(), secantY_at_xhi );
      } else if ( gt.getLowX() <= secantX_at_yhi &&
      secantX_at_yhi <= gt.getHighX() ) {
        gt.movePen( secantX_at_ylo, gt.getLowY() );
        gt.drawLine( secantX_at_yhi, gt.getHighY() );
      }
    } else if ( gt.getLowY() <= secantY_at_xhi &&
    secantY_at_xhi <= gt.getHighY() ) {
      if ( gt.getLowX() <= secantX_at_yhi &&
      secantX_at_yhi <= gt.getHighX() ) {
        gt.movePen( gt.getHighX(), secantY_at_xhi );
        gt.drawLine( secantX_at_yhi, gt.getHighY() );
      }
    }
  } else {
    if ( gt.getLowY() <= f1 && f1 <= gt.getHighY() ) {
      gt.movePen( gt.getLowX(), f1 );
      gt.drawLine( gt.getHighX(), f1 );
    }
  }
}

int main(int argc,char *argv[]) {
  cout << boolalpha;
  { MemoryDebugger md(1);
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

    int nsteps=15;
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

    double oldx,y;
    cout << "click mouse in window to choose starting point" << endl;
    int button=gt.getMouse(oldx,y);
    cout << "click mouse in window again to choose second starting point"
         << endl;
    button=gt.getMouse(x,y);
    double foldx=f(oldx);
    double fx=f(x);

//  GSL iteration
    gsl_root_fdfsolver *solver_gsl
      = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_secant);
    struct secant_params sp;
      sp.x0 = oldx;
    gsl_function_fdf fcn_struct_gsl;
      fcn_struct_gsl.f = &f_gsl;
      fcn_struct_gsl.df = &fp_gsl;
      fcn_struct_gsl.fdf = &fdf_gsl;
      fcn_struct_gsl.params = &sp;
    int status=gsl_root_fdfsolver_set(solver_gsl,&fcn_struct_gsl,x);
    int step=0;
    double epsabs=0.;
    double epsrel=numeric_limits<double>::epsilon();
    double xold=x;
    do {
      step++;
//    first step is Newton; later steps are secant
      status=gsl_root_fdfsolver_iterate(solver_gsl);
      double xnew=gsl_root_fdfsolver_root(solver_gsl);
      cout << "\tstep,root = " << step << " " << xnew << endl;
      if (xnew==xold) break;
      status=gsl_root_test_delta(xnew,xold,epsabs,epsrel);
      if (status==GSL_CONTINUE) {
        status=gsl_root_test_residual(f(xnew),epsrel);
      }
      xold=xnew;
    } while (status == GSL_CONTINUE && step < nsteps);
    double z=gsl_root_fdfsolver_root(solver_gsl);
    cout << "GSL root = " << z << " found in " << step << " iterations"
         << endl;
    gsl_root_fdfsolver_free(solver_gsl); solver_gsl = 0;

//  loop over iterations
    for (int i=0;i<nsteps;i++) {
      double slope=(fx-foldx)/(x-oldx);
      double newx=x-fx/slope;
      double fnewx=f(newx);
      cout << "f(" << newx << ") = " << fnewx << endl;

      logf[i]=log10(max(numeric_limits<double>::epsilon(),abs(fnewx)));
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
      gt.drawCross(oldx,foldx,dx);
      gt.drawCross(x,fx,dx);
      gt.setfgColor("green");
      gt.drawCross(newx,fnewx,dx);
      gt.flush();

      XYGraphTool::WINDOW_TYPE::QuitButton qb;

      oldx=x; foldx=fx;
      x=newx; fx=fnewx;
      if (x==oldx || fx==foldx) {
        nsteps = i+1;
        break;
      }
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
