#include <float.h>
#include <fstream>
#include <limits>
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
inline double fp(double x) { return 2.*x-4.*cos(x); }
inline double fscale(double x) { return max(x*x,4.*sin(x)); }
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

    int nsteps=20;
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

    double x0,y;
    cout << "click mouse in window to choose starting guess" <<endl;
    int button=gt.getMouse(x,y);

//  fixed slope
    x0=x;
    double fpx=fp(x);
    for (int i=0;i<nsteps;i++) {
      double fx=f(x);

      logf[i]=log10(max(numeric_limits<double>::epsilon(),abs(fx)));
      logfmin=min(logfmin,logf[i]);
      logfmax=max(logfmax,logf[i]);

      double s=fx/fpx;
//    cout << "i,s,f,error = " << i << " " << s/abs(x) << " "
//         << fx/fscale(x) << " " << abs(x-z)/abs(z) << endl;
      x=x-s;
    }
    XYGraphTool gt2("fixed slope error","step","log_10(f(x_i))",
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

//  finite difference slope
    XYGraphTool gt3("fixed increment errors","step","log_10(f(x_i))",
      0.,static_cast<double>(nsteps-1),-16.,0.,&cmap,0,winsize);
    gt3.newPage();
    gt3.setbgColor("white");
    gt3.setfgColor("black");
    gt3.drawAxes();
    for (int ic=1;ic<=7;ic++) {
      double h=pow(0.1,ic);
      x=x0;
      for (int i=0;i<nsteps;i++) {
        double fx=f(x);
        double fpx=(f(x+h)-fx)/h;

        logf[i]=log10(max(numeric_limits<double>::epsilon(),abs(fx)));
        logfmin=min(logfmin,logf[i]);
        logfmax=max(logfmax,logf[i]);

        double s=fx/fpx;
        x=x-s;
      }
      gt3.setfgColor(ic,8);
      gt3.movePen(0.,logf[0]);
      for (int i=1;i<nsteps;i++) {
        gt3.drawLine(static_cast<double>(i),logf[i]);
      }
    }
    gt3.flush();

//  vanishing increment
    XYGraphTool gt4("vanishing increment error","step","log_10(f(x_i))",
      0.,static_cast<double>(nsteps-1),-16.,0.,&cmap,0,winsize);
    x=x0;
    double fpx0=fp(x0);
    for (int i=0;i<nsteps;i++) {
      double fx=f(x);
      double h=fx/fpx0;
      double fpx=(f(x+h)-fx)/h;

      logf[i]=log10(max(numeric_limits<double>::epsilon(),abs(fx)));
      logfmin=min(logfmin,logf[i]);
      logfmax=max(logfmax,logf[i]);

      double s=(abs(fpx)>0. ? fx/fpx : 0.);
      x=x-s;
    }
    gt4.newPage();
    gt4.setbgColor("white");
    gt4.setfgColor("black");
    gt4.drawAxes();
    gt4.setfgColor("blue");
    gt4.movePen(0.,logf[0]);
    for (int i=1;i<nsteps;i++) {
      gt4.drawLine(static_cast<double>(i),logf[i]);
    }
    gt4.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    delete [] logf;
  } // pal,cmap,gt go out of scope here 
  return EXIT_SUCCESS;
}
