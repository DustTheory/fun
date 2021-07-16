#include <float.h>
#include <fstream.h>
#include <iostream.h>
#include <math.h>
#include "setTraps.H"
#include "MemoryDebugger.H"
#include "Palette.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XColormap.H"
#include "XYGraphTool.H"
#include <unistd.h>
#define npts 100

inline double f(double x) { return x*x-4.*sin(x); }
//inline double f(double x) { return atan(x); }

//take a=-3.5, b=2.5
//inline double f(double x) { return sin(x)-cos(2.*x); }

//take a=0., b=5. : secant fails
//inline double f(double x) { return exp(-x*x)-1.e-4; }

//take a=0., b=9.999e-3 : secant fails
//inline double f(double x) { double y=1.-100.*x; return x/(y*y); }

//take a=0., b=10.
//inline double f(double x) { return -3.9+x*(12.+x*(-9.+x*2.)); }

void computePlotPoints(double a,double b,double &fmin,double &fmax,
double *xarray,double *farray) {
  double dx=(b-a)*0.125;
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

void plotSecantLine(double a,double b,double x1,double x2,double fx1,
double fx2,XYGraphTool &gt) {
  double slope=(fx2-fx1)/(x2-x1);
  gt.movePen(a,fx1+slope*(a-x1));
  gt.drawLine(b,fx1+slope*(b-x1));
}

int main(int /*argc*/,char* /*argv[]*/) {
#ifdef DEBUG
  setTraps();
#endif

  { MemoryDebugger md(1); {

    int nsteps;
    cout << "enter number steps" << endl;
    cin >> nsteps;
    if (nsteps<0) nsteps=0;
    cout << "number steps = " << nsteps << endl;

    double a,b,x;
    cout << "enter a" << endl;
    cin >> a;
    cout << "enter b" << endl;
    cin >> b;
    if (a>b) { x=a; a=b; b=x; }
    cout << "a,b = " << a << " " << b << endl;

    double xold=a;
    double xnew=b;
    double fold=f(xold);
    double fnew=f(xnew);
    assert(fold*fnew<=0.);

    if (abs(fnew)>abs(fold)) {
      double t=xold; xold=xnew; xnew=t;
      t=fold; fold=fnew; fnew=t;
    }
    cout << "\txold,fold = " << xold << " " << fold << endl;
    cout << "\txnew,fnew = " << xnew << " " << fnew << endl;

    double dx=(b-a)*0.0625;
    double xarray[npts+1],farray[npts+1];
    double fmin,fmax;
    computePlotPoints(a,b,fmin,fmax,xarray,farray);

//  setup interactive graphics
    Palette pal;
    XColormap cmap(&pal);
    REAL winsize=0.5;
    XYGraphTool gt("newton",a,b,fmin,fmax,&cmap,NULL,winsize);

//  loop over iterations
    for (int i=0;i<nsteps;i++) {
      double slopei=(xnew-xold)/(fnew-fold);
      double x=xnew-fnew*slopei;
      double fx=f(x);
      cout << "\tx,fx = " << x << " " << fx << endl;

//    initialize graphics display
      gt.newPage();
      gt.setbgColor("black");
      gt.setfgColor("white");
      gt.drawAxes();

//    draw function
      gt.setfgColor("blue");
//    a=min(xold,xnew); b=max(xold,xnew);
//    dx=(b-a)*0.0625; a-=dx; b+=dx;
//    computePlotPoints(a,b,fmin,fmax,xarray,farray);
//    gt.rescale(xarray[0],xarray[npts],fmin,fmax);
      plotPoints(xarray,farray,gt);

//    mark points
      gt.setfgColor("red");
      plotSecantLine(xarray[0],xarray[npts],xold,xnew,fold,fnew,gt);
      gt.setfgColor("green");
      gt.drawCross(x,fx,dx);

//    force X to perform the requests
      gt.flush();

      wait();
      if (i==0) wait();
      if (fx*fold>=0.) {
	xold=x; fold=fx;
      } else {
	xnew=x; fnew=fx;
      }
    }
  } // pal,cmap,gt go out of scope here 
  }
  return EXIT_SUCCESS;
}

#ifdef __GNUC__
extern "C" {
void MAIN__(void) {;}
}
#endif
