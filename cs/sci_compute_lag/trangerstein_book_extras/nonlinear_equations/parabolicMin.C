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

extern "C" {
  double fmin_(const double &ax,const double &bx,double (*)(const double&),
    const double&);
}

#define npts 100

//double f(const double &x) { return x*x-4.*sin(x); }
double f(const double &x) { 
  double xm1=x-1.;
  return M_PI+xm1*xm1*(20.+xm1*(8.+xm1));
}
double A=-4.,B=4.;

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

void plotParabolic(double x0,double f0,double x1,double f1,
double x2,double f2,XYGraphTool &gt ) {
  if (x1==x0 || x1==x2 || x2==x0) return;
  double s1=(f1-f0)/(x1-x0);
  double s2=(f2-f1)/(x2-x1);
  double ss=(s2-s1)/(x2-x0);
  double dx = ( gt.getHighX() - gt.getLowX() ) / npts;
  bool currently_drawing=false;
  double x=gt.getLowX();
  for (double x=gt.getLowX();x<=gt.getHighX();x+=dx) {
    double y=((x-x1)*ss+s2)*(x-x2)+f2;
    if (y>=gt.getLowY() && y<=gt.getHighY()) {
      if (currently_drawing) gt.drawLine(x,y);
      else {
        gt.movePen(x,y);
        currently_drawing=true;
      }
    } else currently_drawing=false;
  }
}

int main(int argc,char *argv[]) {
  cout << boolalpha;
  { MemoryDebugger md(1);
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
//  double z=fmin_(A,B,f,numeric_limits<double>::epsilon());
    double z=1.;

    int nsteps=25;
    double *loge = OPERATOR_NEW_BRACKET(double,nsteps);
    double logemin=numeric_limits<double>::max();
    double logemax=-numeric_limits<double>::max();

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
    XYGraphTool gt("objective","x","f",A,B,fmin,fmax,&cmap,0,winsize);

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

    double x2,x1,f2,f1,fx,y;
    cout << "click mouse in window to choose starting point" << endl;
    int button=gt.getMouse(x2,y);
    f2=f(x2);
    cout << "click mouse in window again to choose second starting point"
         << endl;
    button=gt.getMouse(x1,y);
    f1=f(x1);
    if (f2<f1) {
      x=x2; x2=x1; x1=x;
      y=f2; f2=f1; f1=y;
    }
    cout << "click mouse in window again to choose second starting point"
         << endl;
    button=gt.getMouse(x,y);
    fx=f(x);
    if (fx>f2) {
      double t=x; x=x1; x1=x2; x2=t;
      y=fx; fx=f1; f1=f2; f2=y;
    } else if (fx>f1) {
      double t=x; x=x1; x1=t;
      y=fx; fx=f1; f1=y;
    }
    double s12=(f2-f1)/(x2-x1);
    cout << setprecision(16);

//  loop over iterations
    for (int i=0;i<nsteps;i++) {
      double s01=(f1-fx)/(x1-x);
      double ss=(s01-s12)/(x-x2);
      if (ss==0.) {
        nsteps = i;
        break;
      }
      double newx=0.5*(x+x1-s01/ss);
      double fnewx=f(newx);
      cout << "f(" << newx << ") = " << fnewx << endl;

      loge[i]=log10(max(numeric_limits<double>::epsilon(),abs(newx-z)));
      logemin=min(logemin,loge[i]);
      logemax=max(logemax,loge[i]);

      gt.newPage();
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      plotPoints(xarray,farray,gt);
      gt.setfgColor("brown");
      plotParabolic(x2,f2,x1,f1,x,fx,gt);

//    mark points
      gt.setfgColor("red");
      gt.drawCross(x2,f2,dx);
      gt.drawCross(x1,f1,dx);
      gt.drawCross(x,fx,dx);
      gt.setfgColor("green");
      gt.drawCross(newx,fnewx,dx);
      gt.flush();
      XYGraphTool::WINDOW_TYPE::QuitButton qb;

      x2=x1; f2=f1;
      x1=x; f1=fx;
      s12=s01;
      x=newx; fx=fnewx;
      if (x1==x || x==x2) {
        nsteps = i+1;
        break;
      }
    }
    XYGraphTool gt2("error","step","log_10(x_i-z))",
      0.,static_cast<double>(nsteps-1),logemin,logemax,&cmap,0,winsize);
    gt2.newPage();
    gt2.setbgColor("white");
    gt2.setfgColor("black");
    gt2.drawAxes();
    gt2.setfgColor("blue");
    gt2.movePen(0.,loge[0]);
    for (int i=1;i<nsteps;i++) {
      gt2.drawLine(static_cast<double>(i),loge[i]);
    }
    gt2.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    delete [] loge;
  } // pal,cmap,gt go out of scope here 
  return EXIT_SUCCESS;
}
