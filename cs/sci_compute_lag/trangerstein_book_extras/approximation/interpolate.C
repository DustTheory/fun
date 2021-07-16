#include <iostream>
#include <math.h> // for HUGE_VAL,M_PI
#include <stdlib.h>

using namespace std;
//#define RUNGE
//#define CHEBYSHEV

//#include "Debug.H"
#include "MemoryDebugger.H"
//#include "SetTraps.H"
#include "TimedObject.H"
#include "GTKColormap.H"
#include "XYGraphTool.H"
extern "C" {
  void addpoint_(int&,const double&,const double&,double*,double*);
  void divdif_(const int&,const double*,const double*,double*);
  double lagrangepoly_(const int&,const double&,const double*,
    const double*);
  double newtonpoly_(const int&,const double*,const double&,
    const double*);
}

double fcn(double t) { 
#ifndef RUNGE
  return sin(t); // exercise 9.2 #3 
#else
  return 1./(1.+25.*t*t); // runge example
#endif
}

int main(int argc,char *argv[]) {
  {
#ifdef MEM_DEBUG
    MemoryDebugger md(1);
#endif
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    char index_field[3];
    char filename[LENGTH_NAME];

#ifndef RUNGE
    double xmin=0.;
    double xmax=M_PI;
#else
    double xmin=-1.;
    double xmax=1.;
#endif
    double ymin=0.;
    double ymax=1.;
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    XYGraphTool gt(
#ifdef RUNGE
      "Runge function",
#else
      "sine function",
#endif
      "x","y",xmin,xmax,ymin,ymax,&cmap,0,0.5);
    gt.setbgColor("white");

    int nmin=2;
#ifndef RUNGE
    int nmax=30; //larger values lead to noticeable effects from
                 //rounding errors in divided difs
#else
    int nmax=12;
#endif

    double *difs=OPERATOR_NEW_BRACKET(double,nmax);
    double *errors=OPERATOR_NEW_BRACKET(double,nmax);
    double *x=OPERATOR_NEW_BRACKET(double,nmax);
    double *y=OPERATOR_NEW_BRACKET(double,nmax);

    double errors_min=DBL_MAX;
    double errors_max=-DBL_MAX;
    double log10=log(10.);
    for (int i=0;i<nmax;i++) {
      errors[i]=HUGE_VAL;
    }
    TimedObject tl("Lagrange interpolation");
    TimedObject tn("Newton interpolation");
    for (int n=nmin;n<nmax;n++) {
#ifdef INDEF
      for (int i=0;i<=n;i++) {
        difs[i]=x[i]=y[i]=HUGE_VAL;
      }
#endif

      double dx=(xmax-xmin)/static_cast<double>(n);
#ifdef CHEBYSHEV
      double dxc=0.5*M_PI/static_cast<double>(n+1);
      double dxw=0.5*(xmax-xmin);
#endif
      for (int j=0;j<=n;j++) {
#ifndef CHEBYSHEV
//      equally-spaced interpolation points
        x[j]=xmin+dx*static_cast<double>(j);
#else
//      chebyshev interpolation points on [0,M_PI]:
        x[j]=xmin+dxw*(1.-cos(dxc*static_cast<double>(2*j+1)));
#endif
        y[j]=fcn(x[j]);
//      cout << "\tx,y[" << j << "] = " << x[j] << " " << y[j] << endl;
      }
      ymin=HUGE_VAL;
      ymax=-HUGE_VAL;
      divdif_(n,x,y,difs);
      dx*=0.01;
      for (int j=0;j<=100*n;j++) {
        double t=xmin+dx*static_cast<double>(j);
        double z=newtonpoly_(n,difs,t,x);
        double ft=fcn(t);
        ymin=min(ymin,min(ft,z));
        ymax=max(ymax,max(ft,z));
      }
      gt.rescale(xmin,xmax,ymin,ymax);

      double t=xmin;
#ifdef RUNGE
//    if (n==nmin) {
#endif
        gt.newPage();
        gt.setfgColor("black");
        gt.drawAxes();
        gt.setfgColor("blue");
        gt.movePen(t,fcn(t));
        for (int j=0;j<=100*n;j++) {
          double t=xmin+dx*static_cast<double>(j);
          double z=fcn(t);
          gt.drawLine(t,z);
        }
#ifdef RUNGE
//    }
#endif
      double maxerr=0.;
#ifndef RUNGE
      gt.setfgColor("red");
#else
      gt.setfgColor(n-nmin+1,nmax-nmin);
#endif
      t=xmin;
      gt.movePen(t,newtonpoly_(n,difs,t,x));
      for (int j=0;j<=100*n;j++) {
        double t=xmin+dx*static_cast<double>(j);
        double z=newtonpoly_(n,difs,t,x);
        gt.drawLine(t,z);
        maxerr=max(maxerr,abs(z-fcn(t)));
//      cout << t << " " << z << endl;
      }
      errors[n]=log(maxerr)/log10;
      errors_min=min(errors_min,errors[n]);
      errors_max=max(errors_max,errors[n]);
      gt.flush();
      cout << "errors[" << n << "] = " << maxerr << endl;

      { Timer zn(&tn);
        for (int j=0;j<=100*n;j++) {
          double t=xmin+dx*static_cast<double>(j);
          double z=newtonpoly_(n,difs,t,x);
        }
      }
      { Timer zl(&tl);
        for (int j=0;j<=100*n;j++) {
          double t=xmin+dx*static_cast<double>(j);
          double z=lagrangepoly_(n, t,x,y);
        }
      }
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }
    cout << "Lagrange interpolation took " << tl.totalRunTime() 
         << " seconds" << endl;
    cout << "Newton interpolation took " << tn.totalRunTime() 
         << " seconds" << endl;

    XYGraphTool gt2("interpolation errors",
      "log_10(number interpolation points)","log_10(error)",
      log(static_cast<double>(nmin))/log10,
      log(static_cast<double>(nmax))/log10,errors_min,errors_max,&cmap,
      0,0.5);
    gt2.setbgColor("white");
    gt2.setfgColor("black");
    gt2.drawAxes();
    gt2.setfgColor("blue");
    double t=log(static_cast<double>(nmin))/log10;
    gt2.movePen(t,errors[nmin]);
    for (int n=nmin+1;n<nmax;n++) {
      t=log(static_cast<double>(n))/log10;
      gt2.drawLine(t,errors[n]);
    }
    gt2.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    delete [] y;
    delete [] x;
    delete [] errors;
    delete [] difs;
  }
}
