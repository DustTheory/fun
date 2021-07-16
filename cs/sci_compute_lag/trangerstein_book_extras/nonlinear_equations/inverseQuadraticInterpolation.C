#include <float.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
  
#include "SetTraps.H"
#include "MemoryDebugger.H"
#include "Palette.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XYGraphTool.H"
#define npts 100

inline double f(double x) { return x*x-4.*sin(x); }
double A=-4.,B=4.;

//inline double f(double x) { return atan(x); }
//double A=-4.,B=4.;

//inline double f(double x) { return sin(x)-cos(2.*x); }
//double A=-3.5, B=2.5;

//inline double f(double x) { return exp(-x*x)-1.e-4; }
//double A=0., B=5.; // secant fails

//inline double f(double x) { return -3.9+x*(12.+x*(-9.+x*2.)); }
//double A=0., B=10.;

void computePlotPoints(double a,double b,double &fmin,double &fmax,
double *xarray,double *farray) {
//double dx=(b-a)*0.125;
//a-=dx; b+=dx;

  double dx=(b-a)/static_cast<double>(npts);
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

void plotInverseQuadratic(double x1,double f1,double x2,double f2,
double x3,double f3,XYGraphTool &gt ) {
//TRACER_CALL(t,"plotInverseQuadratic");
  double dy=(gt.getHighY()-gt.getLowY())/static_cast<double>(npts);
//cout << "\tx1,f1,interp = " << x1 << " " << f1 << " "
//     << (x1-ar)/(br*x1+cr) << endl;
//cout << "\tx2,f2,interp = " << x2 << " " << f2 << " "
//     << (x2-ar)/(br*x2+cr) << endl;
//cout << "\tx3,f3,interp = " << x3 << " " << f3 << " "
//     << (x3-ar)/(br*x3+cr) << endl;
  gt.movePen(x3,f3);
  for (double y=f3-dy;y>=gt.getLowY();y-=dy) {
    double r1=(y-f1)/(f2-f3);
    double r2=(y-f2)/(f3-f1);
    double r3=(y-f3)/(f1-f2);
    double x=-x1*r2*r3-x2*r3*r1-x3*r1*r2;
    if (x<=gt.getHighX() && x>=gt.getLowX()) gt.drawLine(x,y);
    else gt.movePen(x,y);
  }
  gt.movePen(x3,f3);
  for (double y=f3+dy;y<=gt.getHighY();y+=dy) {
    double r1=(y-f1)/(f2-f3);
    double r2=(y-f2)/(f3-f1);
    double r3=(y-f3)/(f1-f2);
    double x=-x1*r2*r3-x2*r3*r1-x3*r1*r2;
    if (x<=gt.getHighX() && x>=gt.getLowX()) gt.drawLine(x,y);
    else gt.movePen(x,y);
  }
}

int main(int argc,char *argv[]) {
  cout << boolalpha;
  { MemoryDebugger md(1); {
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    double ninth=1./9.;
    int nsteps=50;
    double dx=(B-A)*0.015625;

    double xarray[npts+1],farray[npts+1];
    double fmin,fmax;
    computePlotPoints(A,B,fmin,fmax,xarray,farray);

//  setup interactive graphics
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("rational interpolation","x","f",A,B,fmin,fmax,&cmap,0,
      winsize);
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();
    gt.setfgColor("blue");
    plotPoints(xarray,farray,gt);
    gt.flush();

    double xold,xnew,x;
    cout << "click mouse in window to choose starting point" << endl;
    int button=gt.getMouse(xold,x);
    gt.setfgColor("green");
    gt.movePen(xold,gt.getLowY());
    gt.drawLine(xold,gt.getHighY());
    gt.flush();
    cout << "click mouse in window again to choose second starting point"
         << " with function value of opposite sign" << endl;
    button=gt.getMouse(xnew,x);

    double fold=f(xold);
    double fnew=f(xnew);
    CHECK_TEST(fold*fnew<=0.);

    if (abs(fnew)>abs(fold)) {
      x=fold; fold=fnew; fnew=x;
      x=xold; xold=xnew; xnew=x;
    }
    cout << "\txold,fold = " << xold << " " << fold << endl;
    cout << "\txnew,fnew = " << xnew << " " << fnew << endl;
    cout << "\tsecant step" << endl;

    double zeta=max(abs(xold),abs(xnew));
    double slopei=(xnew-xold)/(fnew-fold);
    x=xnew-fnew*slopei;
    double xdif=abs(x-xnew);
    if (xdif<=numeric_limits<double>::epsilon()*zeta) {
      cout << "\tsmall relative error in solution" << endl;
      return EXIT_SUCCESS;
    }
    if (xdif<=sqrt(numeric_limits<double>::epsilon())*zeta &&
    xdif >= abs(xnew-xold)) {
      cout << "\tmaximum attainable accuracy" << endl;
      return EXIT_SUCCESS;
    }

//  loop over iterations
    bool special=false;
    for (int i=0;i<nsteps;i++) {
      double fx=f(x);
      if (fx==0.) {
        cout << "\tzero function value" << endl;
        break;
      }
      cout << "\tx,fx = " << x << " " << fx << endl;
//    cout << "\txold,xnew,x,fx = " << xold << " " << xnew << " " << x
//         << " " << fx << endl;

//    initialize graphics display
      gt.newPage();
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();

//    draw function
      gt.setfgColor("blue");
      plotPoints(xarray,farray,gt);
      gt.setfgColor("green");
      gt.drawCross(x,fx,dx);
      gt.movePen(xold,gt.getLowY());
      gt.drawLine(xold,gt.getHighY());
      gt.movePen(xnew,gt.getLowY());
      gt.drawLine(xnew,gt.getHighY());
      gt.setfgColor("red");

      if (fx*fnew<0.) {
        cout << "\tsecant step" << endl;
	xold=xnew; fold=fnew; 
        xnew=x; fnew=fx;
        plotSecantLine(xold,fold,xnew,fnew,gt);
        zeta=max(abs(xold),abs(xnew));
        slopei=(xnew-xold)/(fnew-fold);
        x=xnew-fnew*slopei;
        xdif=abs(x-xnew);
        if (abs(x-xnew)<=numeric_limits<double>::epsilon()*zeta) {
          cout << "\tsmall relative error in solution" << endl;
          break;
        }
        if (xdif<=sqrt(numeric_limits<double>::epsilon())*zeta &&
        xdif >= abs(xnew-xold)) {
          cout << "\tmaximum attainable accuracy" << endl;
          break;
        }
        special=false;
      } else {
        plotInverseQuadratic(xold,fold,xnew,fnew,x,fx,gt);
        if (fx!=fnew) {
          double rhoa=fold/(fnew-fx);
          double rhob=fnew/(fx-fold);
          double rhoc=fx/(fold-fnew);
          double z=-xold*rhoc*rhob-xnew*rhoa*rhoc-x*rhoa*rhob;
          double rho=(z-min(xold,xnew))/abs(xnew-xold);
          xnew=x; fnew=fx;
          if ((0.<rho && rho<1.) && (!special || (0.1<rho && rho<0.9))) {
            cout << "\tinverse quadratic interpolant" << endl;
            x=z;
          } else {
            cout << "\tout of bounds ==> bisection" << endl;
            x=(xold+xnew)*0.5;
          }
        } else {
          cout << "\tinterpolant undefined ==> bisection" << endl;
          xnew=x; fnew=fx;
          x=(xold+xnew)*0.5;
        }
        special=true;
      }
      gt.flush();
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
      if (x==xold || x==xnew) {
        cout << "\tno change in solution" << endl;
        break;
      }
    }
    cout << setprecision(16);
    cout << "\tx = " << x << endl;
  } // pal,cmap,gt go out of scope here 
  }
  return EXIT_SUCCESS;
}
