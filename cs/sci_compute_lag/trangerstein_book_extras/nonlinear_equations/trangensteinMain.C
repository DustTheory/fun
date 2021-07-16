#include <float.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include "SetTraps.H"
#include "MemoryDebugger.H"
#include "Palette.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XYGraphTool.H"
#include <unistd.h>
#define npts 100

//inline double f(double x) { return x*x-4.*sin(x); }
//inline double typf(double x) { return max(x*x,4.*sin(x)); }

  inline double f(double x) { return atan(x); }
  inline double typf(double x) { return M_PI_4; }

//take a=-3.5, b=2.5
//inline double f(double x) { return sin(x)-cos(2.*x); }
//inline double typf(double x) { return sin(x); }

//take a=0., b=5. : secant fails
//inline double f(double x) { return exp(-x*x)-1.e-4; }
//inline double typf(double x) { return 1.e-4; }

//take a=0., b=9.999e-3 : secant fails
//inline double f(double x) { double y=1.-100.*x; return x/(y*y)-1.; }
//inline double typf(double x) { return 1.; }

//take a=0., b=10.
//  also try a=.32, b=2.8
//inline double f(double x) { return -3.9+x*(12.+x*(-9.+x*2.)); }
//inline double typf(double x) { return 3.9; }

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

void plotMullerTraub(double *xarray,double x0,double x1,double x2,
double f0,double f1,double f2,XYGraphTool &gt) {
  double slope10=(f1-f0)/(x1-x0);
  double slope21=(f2-f1)/(x2-x1);
  double seconddif=(slope21-slope10)/(x2-x0);
  gt.movePen(xarray[0],
             f2+(xarray[0]-x2)*(slope21+(xarray[0]-x1)*seconddif));
  for (int i=1;i<=npts;i++) {
    gt.drawLine(xarray[i],
               f2+(xarray[i]-x2)*(slope21+(xarray[i]-x1)*seconddif));
  }
}

void plotRationalPoly(double *xarray,double x0,double x1,double x2,
double f0,double f1,double f2,XYGraphTool &gt) {
  double slope10=(f1-f0)/(x1-x0);
  double slope21=(f2-f1)/(x2-x1);
  double beta=slope21/slope10;
  double B=1.-beta;
  double C=x2*beta-x0;
  double den=f2-beta*f0;
  double A=x2-f2*(x2-x0)/den;
//cout << "\n\tzero of rational polynomial = " << A << endl;
  B /= den;
  C /= den;
  gt.movePen(xarray[0],(xarray[0]-A)/(B*xarray[0]+C));
  for (int i=1;i<=npts;i++) {
    gt.drawLine(xarray[i],(xarray[i]-A)/(B*xarray[i]+C));
  }
}

int main(int argc,char *argv[]) {
#ifdef DEBUG
//setTraps();
#endif

  { MemoryDebugger md(1); {

    double conv=DBL_EPSILON;
    double mindx=0.125;
//  allow muller-traub to work unimpeded within 4 steps of the solution
    double term=0.5*(1.+sqrt(5.));
    term*=term; term*=term;
    double cutoff=exp( log(conv)/term );

    int nsteps=20;
//  cout << "enter number steps" << endl;
//  cin >> nsteps;
//  if (nsteps<0) nsteps=0;
//  cout << "number steps = " << nsteps << endl;

    double a,b,x;
    cout << "enter graph lower x coordinate" << endl;
    cin >> a;
    cout << "enter graph upper x coordinate" << endl;
    cin >> b;
    if (a>b) { x=a; a=b; b=x; }
    cout << "a,b = " << a << " " << b << endl;

    double dx=(b-a)*0.0625;
    double xarray[npts+1],farray[npts+1];
    double fmin,fmax;
    computePlotPoints(a,b,fmin,fmax,xarray,farray);

//  setup interactive graphics
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("newton","x","f",a,b,fmin,fmax,&cmap,NULL,winsize);

//  initialize graphics display
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

//  draw function
    gt.setfgColor("blue");
    plotPoints(xarray,farray,gt);
    gt.flush();

    cout << "click mouse in graph tool to select x0" << endl;
    double xold=a;
    double fold=0.;
    gt.getMouse(xold,fold);
    fold=f(xold);
    gt.setfgColor("green");
    gt.drawCross(xold,fold,dx);
    gt.flush();

    cout << "click mouse in graph tool to select x1" << endl;
    double xnew=b;
    double fnew=0.;
    gt.getMouse(xnew,fnew);
    fnew=f(xnew);
    gt.drawCross(xnew,fnew,dx);
    gt.flush();

    CHECK_TEST(fold*fnew<=0.);

//  if (abs(fnew)>abs(fold)) {
//    double t=xold; xold=xnew; xnew=t;
//    t=fold; fold=fnew; fnew=t;
//  }
    cout << "\txold,fold = " << xold << " " << fold << endl;
    cout << "\txnew,fnew = " << xnew << " " << fnew << endl;

//  loop over iterations
    double fakef=fold;
    bool previous_rational=false;
    for (int i=0;i<nsteps;i++) {
      double slopei=(xnew-xold)/(fnew-fakef);
      double x=xnew-fnew*slopei;
      double fx=f(x);
      double fsize=typf(x);
      cout << "\tit,x,fx = " << i << " " << x << " " << fx << " " 
           << endl;
      double differ=x-xnew;
//    if (abs(fx)<=conv*fsize) break;

//    initialize graphics display
      gt.newPage();
      gt.setbgColor("white");
      gt.setfgColor("black");
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
      if (fx*fnew<=0.) {
        plotSecantLine(xarray[0],xarray[npts],x,xnew,fx,fnew,gt);
      } else {
        plotSecantLine(xarray[0],xarray[npts],x,xold,fx,fold,gt);
      }
      gt.setfgColor("yellow");
      plotRationalPoly(xarray,xold,xnew,x,fold,fnew,fx,gt);
      gt.setfgColor("cyan");
      plotMullerTraub(xarray,xold,xnew,x,fold,fnew,fx,gt);
      gt.setfgColor("green");
      gt.drawCross(xold,fold,dx);
      gt.drawCross(xnew,fnew,dx);
      gt.drawCross(x,fx,dx);

//    force X to perform the requests
      gt.flush();

      if (i==0) wait();
      if (fx*fold<=0.) {
//      rational polynomial interpolation
        double slopen=(fx-fnew)/differ;
        double beta=slopei*slopen;
        if (0.<beta && beta<1. && !previous_rational) {
          fakef=fold*beta;
          previous_rational=true;
          cout << "try rational interpolation" << endl;
        } else {
//      quadratic interpolation
          double seconddif=(slopen-(fnew-fold)/(xnew-xold))/(x-xold);
          double omega=slopen+seconddif*(x-xnew);
          double radical=sqrt(omega*omega-4.*fx*seconddif);
          double fprime=(fx*(x-xold)>0. ? omega+radical :omega-radical);
          fakef=fx-0.5*(x-xold)*fprime;
          previous_rational=false;
          cout << "try quadratic interpolation" << endl;
        }
        double temp=-fakef/fx;
        if ((temp<mindx) || (temp>1./mindx) && abs(fx)>cutoff*fsize) {
          fakef=-fx;
          previous_rational=false;
          cout << "temp = " << temp << " ==> use bisection" << endl;
        } else {
          cout << "use higher-order interpolation" << endl;
        }
      } else {
        xold=xnew; fakef=fnew; fold=fnew; 
        previous_rational=false;
        cout << "secant" << endl;
      }

//    if (abs(xnew-xold)<=conv*max(abs(xold),abs(xnew))) break;
      xnew=x; fnew=fx;

//    cout << "xold,fold = " << xold << " " << fold << endl;
//    cout << "xnew,fnew = " << xnew << " " << fnew << endl;
//    cout << "fakef = " << fakef << endl;
      wait();
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
