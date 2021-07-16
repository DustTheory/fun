#include <iostream>
#include <limits>
//#include <math.h> // for HUGE_VAL,M_PI
#include <stdlib.h>

using namespace std;

#include "Arch.H"
#include "Debug.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "GTKColormap.H"
#include "XYGraphTool.H"
extern "C" {
  double F77_NAME(c0_spline_eval)(const int &element,const int &nelements,
    const int &order,const double &pt,const double *x,const double *xi,
    const double *y,double *p,double *q);
  void F77_NAME(c1_quadratic_coefs)(const int &nelements,const double *x,
    const double *y,const double *yp,double *z);
  double F77_NAME(c1_quadratic_eval)(const int &element,
    const int &nelements,const double &pt,const double *x,const double *y,
    double *z);
  double F77_NAME(c1_spline_eval)(const int &element,const int &nelements,
    const int &order,const double &pt,const double &sigma,const double *x,
    const double *xi,const double *y,const double *yp);
  void F77_NAME(c2_cubic_coefs)(const int &nelements,double *d,double *e,
    const double *x,const double *y,const double *yp,double *z);
  double F77_NAME(c2_cubic_eval)(const int &element,const int &nelements,
    const double &pt,const double *x,const double *y,double *z);
  void F77_NAME(c2_quartic_coefs)(const int &nelements,const double *x,
    const double *y,const double *yp,double *ypp,double *z);
  double F77_NAME(c2_quartic_eval)(const int &element,const int &nelements,
    const double &pt,const double *x,const double *y,const double *yp,
    double *z);
  double F77_NAME(c2_spline_eval)(const int &element,const int &nelements,
    const int &order,const double &pt,const double &sigma,
    const double &tau,const double *x,const double *xi,const double *y,
    const double *yp,const double *ypp);
}

//
double fcn(double t) { 
  return 1./(1.+25.*t*t); // runge example
}
double fcn_deriv(double t) { 
  double den=1.+25.*t*t;
  return -50.*t/(den*den);
}
double fcn_2nd_deriv(double t) { 
  double den=1.+25.*t*t;
  return 50.*(75.*t*t-1.)/(den*den*den);
}
//
/*
double fcn(double t) { 
  return pow(t,7);
}
double fcn_deriv(double t) { 
  return 7.*pow(t,6);
}
double fcn_2nd_deriv(double t) { 
  return 42.*pow(t,5);
}
*/
//selection sort:
void sort(const int &n,double *x) {
  int j=0;
  int jnext=1;
  for (;jnext<=n;j=jnext,jnext++) {
    int jsmall=j;
    for (int k=jnext;k<=n;k++) {
      if (x[jsmall]>x[k]) jsmall=k;
    }
    if (j!=jsmall) {
      double t=x[j];
      x[j]=x[jsmall];
      x[jsmall]=t;
    }
  }
}

int bin_search(int high,int low, double t,const double *x) {
  for (; high-low>1;) {
    int i=(high+low)/2;
    if (t<x[i]) high=i;
    else low=i;
  }
  return low;
}

int main(int argc,char *argv[]) {
#ifdef DEBUG
//setTraps();
#endif
  {
    MemoryDebugger md(1);
    {
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

    double xmin=-1.;
    double xmax=1.;
    double ymin=-1.;
    double ymax=1.;
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    XYGraphTool gt("interpolation","x","y",xmin,xmax,ymin,ymax,&cmap,0,0.5);
    gt.setbgColor("white");

    int nmin=2;
    int nmax=13;

    int order=1;
    int continuity=0;
#ifdef DEBUG
    cout << "\torder = " << order << " continuity = " << continuity << endl;
#endif

    double *errors_alloc=OPERATOR_NEW_BRACKET(double,nmax-nmin+1);
    double *errors=errors_alloc-nmin;
    double *x=OPERATOR_NEW_BRACKET(double,nmax+1);
    double *y=OPERATOR_NEW_BRACKET(double,nmax*max(1,order-2*continuity)+1);
    double *yp=OPERATOR_NEW_BRACKET(double,
        ((order>continuity+1 && continuity>0) ? nmax+1 : continuity));
    double *ypp=OPERATOR_NEW_BRACKET(double,
        ((order>2*continuity && continuity>1) ? nmax+1 : 1 ));
    double *xi=OPERATOR_NEW_BRACKET(double,order-2*continuity+1);
    double *z=OPERATOR_NEW_BRACKET(double,nmax+1);

    double errors_min=numeric_limits<double>::max();
    double errors_max=-numeric_limits<double>::max();
    double log10=log(10.);
    for (int i=nmin;i<=nmax;i++) {
      errors[i]=numeric_limits<double>::max();
    }
    TimedObject tl("Lagrange interpolation");
    TimedObject tn("Newton interpolation");
    for (int n=nmin;n<nmax;n++) {
#ifdef INDEF
      for (int i=0;i<=n;i++) {
        x[i]=numeric_limits<double>::infinity();
      }
      for (int i=0;i<=n*max(1,order-2*continuity)+1;i++) {
        y[i]=numeric_limits<double>::infinity();
      }
#endif

//    generate mesh
      double dx=(xmax-xmin)/static_cast<double>(n);
      x[0]=xmin;
      for (int i=1;i<n;i++) x[i]=x[i-1]+dx;
      x[n]=xmax;
#ifdef DEBUG
//    for (int i=0;i<=n;i++) {
//      cout << "x[" << i << "] = " << x[i] << endl;
//    }
#endif
/*
//    random mesh nodes
      x[0]=xmin;
      for (int j=1;j<n;j++) {
        x[j]=xmin+(xmax-xmin)
            *static_cast<double>(rand())/static_cast<double>(RAND_MAX));
      }
      x[n]=xmax;
      sort(n,x);
*/
      int ninterior=order-2*continuity;
      if (ninterior>0) {
        double dxi=1./static_cast<double>(ninterior);
        xi[0]=0.;
        for (int j=1;j<ninterior;j++) xi[j]=xi[j-1]+dxi;
        xi[ninterior]=1.;
/*
//      random interior points
        for (int j=1;j<=ninterior/2;j++) {
          xi[j]=0.5*static_cast<double>(rand())
                   /static_cast<double>(RAND_MAX);
        }
        if (ninterior%2==0) xi[ninterior/2]=0.5;
        xi[0]=0.;
        sort(order/2-continuity,xi+1);
//      place symmetrically:
        for (int j=order/2-continuity+1,j<=order-2*continuity;j++) {
          xi[j]=1.-xi[order-2*continuity-j];
        }
*/
#ifdef DEBUG
//      for (int i=0;i<=ninterior;i++) {
//        cout << "xi[" << i << "] = " << xi[i] << endl;
//      }
#endif
      }

      int nint=max(1,ninterior);
      for (int j=0;j<=n;j++) {
        y[j*nint]=fcn(x[j]);
      }
      if (ninterior>0) {
        for (int j=0;j<n;j++) {
          double dx=x[j+1]-x[j];
          int i=j*ninterior;
          for (int k=1;k<ninterior;k++) {
            y[i+k]=fcn(x[j]+dx*xi[k]);
          }
        }
      }
#ifdef DEBUG
//    for (int j=0;j<=n*max(1,order-2*continuity);j++) {
//      cout << "y[" << j << "] = " << y[j] << endl;
//    }
#endif
      if (continuity>0) {
        if (order>continuity+1) {
          for (int j=0;j<=n;j++) yp[j]=fcn_deriv(x[j]);
#ifdef DEBUG
          for (int j=0;j<=n;j++) {
            cout << "yp[" << j << "] = " << yp[j] << endl;
          }
#endif
        } else {
          yp[0]=fcn_deriv(xmin);
          if (continuity==2) yp[1]=fcn_deriv(xmax);
#ifdef DEBUG
          cout << "yp[0] = " << yp[0] << endl;
          if (continuity==2) cout << "yp[1] = " << yp[1] << endl;
#endif
        }
      }
      if (continuity>1) {
        if (order>2*continuity) {
          for (int j=0;j<=n;j++) ypp[j]=fcn_2nd_deriv(x[j]);
#ifdef DEBUG
          for (int j=0;j<=n;j++) {
            cout << "ypp[" << j << "] = " << ypp[j] << endl;
          }
#endif
        } else {
          ypp[0]=fcn_2nd_deriv(xmin);
#ifdef DEBUG
          cout << "ypp[0] = " << ypp[0] << endl;
#endif
        }
      }
      gt.newPage();
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      gt.movePen(xmin,fcn(xmin));
      for (int j=0;j<=128*n;j++) {
        int i=j/128;
        double dx=(x[i+1]-x[i])/128.;
        double t=x[i]+dx*static_cast<double>(j-i*128);
        gt.drawLine(t,fcn(t));
      }
      double maxerr=0.;
      gt.setfgColor("red");
      if (continuity==0) {
//      TRACER_CALL(tt,"C0 spline");
        double *p=OPERATOR_NEW_BRACKET(double,order+1);
        double *q=OPERATOR_NEW_BRACKET(double,order);
        int element=0;
        double yy=
          F77_NAME(c0_spline_eval)(element,n,order,xmin,x,xi,y,p,q);  
        gt.movePen(xmin,yy);
#ifdef DEBUG
//      cout << xmin << " " << yy << endl;
#endif
        double dx=(xmax-xmin)/static_cast<double>(128*n);
        double t=xmin+dx;
        for (int j=1;j<=128*n;j++,t+=dx) {
          if (j==128*n) t=xmax;
          element=bin_search(n,element, t,x);
          yy=F77_NAME(c0_spline_eval)(element,n,order,t,x,xi,y,p,q);  
          gt.drawLine(t,yy);
          maxerr=max(maxerr,abs(yy-fcn(t)));
#ifdef DEBUG
//        cout << t << " " << yy << endl;
#endif
        }
        delete [] p; p=0;
        delete [] q; q=0;
      } else if (continuity==1) {
        if (order==2) {
//        TRACER_CALL(tt,"C1 quadratic");
          F77_NAME(c1_quadratic_coefs)(n,x,y,yp,z);
#ifdef DEBUG
//        cout << "\tafter c1_quadratic_coefs" << endl;
//        mem_check();
#endif
          int element=0;
          double yy=F77_NAME(c1_quadratic_eval)(element,n,xmin,x,y,z);  
#ifdef DEBUG
//        cout << "\tafter c1_quadratic_eval" << endl;
//        mem_check();
#endif
          gt.movePen(xmin,yy);
          double dx=(xmax-xmin)/static_cast<double>(128*n);
          double t=xmin+dx;
          for (int j=1;j<=128*n;j++,t+=dx) {
            if (j==128*n) t=xmax;
            element=bin_search(n,element, t,x);
            yy=F77_NAME(c1_quadratic_eval)(element,n,t,x,y,z);  
#ifdef DEBUG
//          cout << "\tafter c1_quadratic_eval" << endl;
//          mem_check();
#endif
            gt.drawLine(t,yy);
            maxerr=max(maxerr,abs(yy-fcn(t)));
  //        cout << t << " " << yy << endl;
          }
        } else {
//        TRACER_CALL(tt,"C1 spline");
          int element=0;
          double sum=2.;
          for (int i=1;i<=order-3;i++) sum+=1./xi[i];
          double yy=
            F77_NAME(c1_spline_eval)(element,n,order,xmin,sum,x,xi,y,yp);  
          gt.movePen(xmin,yy);
          double dx=(xmax-xmin)/static_cast<double>(128*n);
          double t=xmin+dx;
          for (int j=1;j<=128*n;j++,t+=dx) {
            if (j==128*n) t=xmax;
            element=bin_search(n,element, t,x);
            yy=F77_NAME(c1_spline_eval)(element,n,order,t,sum,x,xi,y,yp);  
            gt.drawLine(t,yy);
            maxerr=max(maxerr,abs(yy-fcn(t)));
  //        cout << t << " " << yy << endl;
          }
        }
      } else {
        if (order==3) {
//        TRACER_CALL(tt,"C2 cubic");
          double *d=OPERATOR_NEW_BRACKET(double,n+1);
          double *e=OPERATOR_NEW_BRACKET(double,n);
          F77_NAME(c2_cubic_coefs)(n,d,e,x,y,yp,z);
          delete [] d; d=0;
          delete [] e; e=0;
          int element=0;
          double yy=F77_NAME(c2_cubic_eval)(element,n,xmin,x,y,z);  
          gt.movePen(xmin,yy);
          double dx=(xmax-xmin)/static_cast<double>(128*n);
          double t=xmin+dx;
          for (int j=1;j<=128*n;j++,t+=dx) {
            if (j==128*n) t=xmax;
            element=bin_search(n,element, t,x);
            double yy=F77_NAME(c2_cubic_eval)(element,n,t,x,y,z);  
            gt.drawLine(t,yy);
            maxerr=max(maxerr,abs(yy-fcn(t)));
  //        cout << t << " " << yy << endl;
          }
        } else if (order==4) {
//        TRACER_CALL(tt,"C2 quartic");
          F77_NAME(c2_quartic_coefs)(n,x,y,yp,ypp,z);
          int element=0;
          double yy=F77_NAME(c2_quartic_eval)(element,n,xmin,x,y,yp,z);  
          gt.movePen(xmin,yy);
          double dx=(xmax-xmin)/static_cast<double>(128*n);
          double t=xmin+dx;
          for (int j=1;j<=128*n;j++,t+=dx) {
            if (j==128*n) t=xmax;
            element=bin_search(n,element, t,x);
            double yy=F77_NAME(c2_quartic_eval)(element,n,t,x,y,yp,z);  
            gt.drawLine(t,yy);
            maxerr=max(maxerr,abs(yy-fcn(t)));
  //        cout << t << " " << yy << endl;
          }
        } else {
//        TRACER_CALL(tt,"C2 spline");
          int element=0;
          double sigma=0.;
          double tau=0.;
          for (int i=1;i<=order-5;i++) {
            double sum=0.;
            for (int k=1;k<=order-5;k++) {
              if (k!=i) sum+=1./xi[k];
            }
            sigma+=1./xi[i];
            tau+=sum/xi[i];
          }
          double yy=F77_NAME(c2_spline_eval)(element,n,order,xmin,sigma,tau,
            x,xi,y,yp,ypp);  
          gt.movePen(xmin,yy);
          double dx=(xmax-xmin)/static_cast<double>(128*n);
          double t=xmin+dx;
          for (int j=1;j<=128*n;j++,t+=dx) {
            if (j==128*n) t=xmax;
            element=bin_search(n,element, t,x);
            yy=F77_NAME(c2_spline_eval)(element,n,order,t,sigma,tau,
              x,xi,y,yp,ypp);  
            gt.drawLine(t,yy);
            maxerr=max(maxerr,abs(yy-fcn(t)));
  //        cout << t << " " << yy << endl;
          }
        }
      }
      maxerr=max(maxerr,numeric_limits<double>::epsilon());
      errors[n]=log(maxerr)/log10;
      errors_min=min(errors_min,errors[n]);
      errors_max=max(errors_max,errors[n]);
      gt.flush();
#ifdef DEBUG
//    cout << "errors[" << n << "] = " << maxerr << endl;
#endif

      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }
#ifdef DEBUG
//  cout << "\tbefore gt2" << endl;
//  md.printOn(cout);
#endif

    {
    XYGraphTool gt2("errors","log_10(number cells)","log_10(error)",
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
    }

#ifdef DEBUG
//  cout << "\tafter gt2" << endl;
//  md.printOn(cout);
#endif
    delete [] z; z=0;
    delete [] ypp; ypp=0;
    delete [] yp; yp=0;
    delete [] y; y=0;
    delete [] xi; xi=0;
    delete [] x; x=0;
    delete [] errors_alloc; errors_alloc=0;
  }
  }
}
