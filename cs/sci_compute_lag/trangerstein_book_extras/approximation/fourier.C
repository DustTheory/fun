#include <cmath>
#include <iostream>
#include <limits>
#include <math.h> // for M_PI
#include <rfftw.h>
#include <stdlib.h>

using namespace std;
//#define GAUSSIAN
#define NPTS 1024

//#define USE_COMPLEX
//#define USE_FFTPACK
#define USE_FFTW

#ifdef USE_FFTPACK
# undef USE_COMPLEX
#endif

#include "Debug.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "XYGraphTool.H"
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#include "GTKColormap.H"
#endif
#ifdef USE_FFTPACK
extern "C" {
  void dfftb_(const int&,double*,const double*);
  void dfftf_(const int&,double*,const double*);
  void dffti_(const int&,double*);
}
#endif

double fcn(double x) { 
#ifdef GAUSSIAN
  x-=M_PI;
  return exp(-2.*x*x); // exercise 9.2 #3 
#else
  double pi3=M_PI/3.;
  return (x<pi3 || 5.*pi3 < x ? 0. : (x<M_PI ? 1. : -1.));

//return cos(2.*x);
#endif
}
#ifdef USE_FFTPACK
double fftpack_poly(double theta,int n,double *fhat) {
  double fp=0.;
  double angle=theta;
  for (int i=1;i<n-1;i+=2,angle+=theta) {
    fp+=fhat[i]*cos(angle)-fhat[i+1]*sin(angle);
  }
  fp*=2.;
  fp+=fhat[0];
  if (n%2==0) fp+=fhat[n-1]*cos(angle);
  return fp/static_cast<double>(n);
}
#endif
void setField(char* field,int field_width,int val) {
  field[field_width]='\0';
  int i=field_width-1;
  int n=val;
  for (;i>=0;i--,n/=10) field[i]='0'+n %10;
}

int main(int argc,char** argv) {
#ifdef DEBUG
//setTraps();
#endif
#ifdef MEM_DEBUG
  MemoryDebugger md(1);
#endif
  {
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    char index_field[3];
    char filename[LENGTH_NAME];

    double xmin=0.;
    double xmax=2.*M_PI;
#ifdef GAUSSIAN
    double ymin=0.;
#else
    double ymin=-1.;
#endif
    double ymax=1.;
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    XYGraphTool gt("interpolation","x","y",xmin,xmax,ymin,ymax,&cmap,0,0.5);
    gt.setbgColor("white");

    int nmin=2;
    int nmax=8;
    int power2=__cmath_power(2,nmax);

#ifdef USE_FFTPACK
    double *f=OPERATOR_NEW_BRACKET(double,power2);
    double *FHAT=OPERATOR_NEW_BRACKET(double,NPTS);
    double *wsave=OPERATOR_NEW_BRACKET(double,2*power2+15);
    double *wsave2=OPERATOR_NEW_BRACKET(double,2*NPTS+15);
#endif
#ifdef USE_FFTW
# ifdef USE_COMPLEX
    fftw_complex *f=OPERATOR_NEW_BRACKET(fftw_complex,power2);
    fftw_complex *fhat=OPERATOR_NEW_BRACKET(fftw_complex,power2);
    fftw_complex *F=OPERATOR_NEW_BRACKET(fftw_complex,NPTS);
    fftw_complex *FHAT=OPERATOR_NEW_BRACKET(fftw_complex,NPTS);
# else
    fftw_real *f=OPERATOR_NEW_BRACKET(fftw_real,power2);
    fftw_real *fhat=OPERATOR_NEW_BRACKET(fftw_real,power2);
    fftw_real *F=OPERATOR_NEW_BRACKET(fftw_real,NPTS);
    fftw_real *FHAT=OPERATOR_NEW_BRACKET(fftw_real,NPTS);
# endif
#endif
    double *errors=OPERATOR_NEW_BRACKET(double,nmax+1);

    double errors_min=numeric_limits<double>::max();
    double errors_max=-numeric_limits<double>::max();
    double log10=log(10.);
    for (int i=0;i<=nmax;i++) {
      errors[i]=numeric_limits<double>::infinity();
    }
#ifdef USE_FFTPACK
    dffti_(NPTS,wsave2);
#endif
#ifdef USE_FFTW
# ifdef USE_COMPLEX
    fftw_plan plan_backward=
      fftw_create_plan(NPTS,FFTW_BACKWARD,FFTW_ESTIMATE);
# else
    rfftw_plan plan_backward=
      rfftw_create_plan(NPTS,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE);
# endif
#endif
    for (int n=4,nm=nmin;nm<=nmax;n*=2,nm++) {
//    TRACER_CALL(t,"loop");
#ifdef INDEF
      for (int i=0;i<power2;i++) {
# ifdef USE_FFTPACK
        f[i]=numeric_limits<double>::infinity();
# endif
# ifdef USE_FFTW
#   ifdef USE_COMPLEX
        f[i].re=f[i].im=numeric_limits<double>::infinity();
#   else
        f[i]=numeric_limits<double>::infinity();
        fhat[i]=numeric_limits<double>::infinity();
#   endif
# endif
      }
      for (int i=0;i<NPTS;i++) {
# ifdef USE_FFTPACK
        FHAT[i]=numeric_limits<double>::infinity();
# endif
# ifdef USE_FFTW
#   ifdef USE_COMPLEX
        FHAT[i].re=FHAT[i].im=numeric_limits<double>::infinity();
#   else
        FHAT[i]=numeric_limits<double>::infinity();
#   endif
# endif
      }
#endif

      double dx=(xmax-xmin)/static_cast<double>(n);
      double x=xmin;
      for (int j=0;j<n;j++,x+=dx) {
#ifdef USE_FFTPACK
        f[j]=fcn(x);
//      cout << "\tt,f[" << j << "] = " << x << " " << f[j] << endl;
#endif
#ifdef USE_FFTW
# ifdef USE_COMPLEX
        f[j].re=fcn(x);
        f[j].im=0.;
# else
        f[j]=fcn(x);
//      cout << "\tt,f[" << j << "] = " << x << " " << f[j] << endl;
# endif
#endif
      }
#ifdef USE_FFTPACK
      dffti_(n,wsave);
      dfftf_(n,f,wsave);
# ifdef DEBUG
//    for (int j=0;j<n;j++) {
//      cout << "f[" << j << "] = " << f[j] << endl;
//    }
# endif
#endif
#ifdef USE_FFTW
# ifdef USE_COMPLEX
      fftw_plan plan_forward=
        fftw_create_plan(n,FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_one(plan_forward,f,fhat);
#   ifdef DEBUG
//    for (int j=0;j<n;j++) {
//      cout << "fhat[" << j << "] = " << fhat[j].re << " + i * "
//           << fhat[j].im << endl;
//    }
#   endif
# else
      rfftw_plan plan_forward=
        rfftw_create_plan(n,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);
      rfftw_one(plan_forward,f,fhat);
#   ifdef DEBUG
//    for (int j=0;j<n;j++) {
//      cout << "fhat[" << j << "] = " << fhat[j] << endl;
//    }
#   endif
# endif
#endif

#ifdef USE_FFTPACK
//    fftpack orders the fourier modes as follows:
//    cos(x*0),cos(x),sin(x),cos(x*2),sin(x*2),...,cos(x*n/2)
      double factor=1./static_cast<double>(n);
      for (int j=0;j<n;j++) {
        FHAT[j]=fhat[j]*factor;
      }
      for (int j=n;j<NPTS;j++) FHAT[j]=0.;
# ifdef DEBUG
//    for (int j=0;j<NPTS;j++) {
//      cout << "FHAT[" << j << "] = " << FHAT[j] << endl;
//    }
# endif
      dfftb_(NPTS,FHAT,wsave2);
#endif
#ifdef USE_FFTW
      fftw_real factor=1./static_cast<fftw_real>(n);
# ifdef USE_COMPLEX
      for (int j=0;j<=n/2;j++) {
        FHAT[j].re=fhat[j].re*factor;
        FHAT[j].im=fhat[j].im*factor;
      }
      for (int j=n/2+1;j<=NPTS-n/2;j++) FHAT[j].re=FHAT[j].im=0.;
      for (int j=1;j<n/2;j++) {
        FHAT[NPTS-j].re=fhat[n-j].re*factor;
        FHAT[NPTS-j].im=fhat[n-j].im*factor;
      }
#   ifdef DEBUG
//    for (int j=0;j<NPTS;j++) {
//      cout << "FHAT[" << j << "] = " << FHAT[j].re << " + i * "
//           << FHAT[j].im << endl;
//    }
#   endif
      fftw_one(plan_backward,FHAT,F);
#   ifdef DEBUG
//    for (int j=0;j<NPTS;j++) {
//      cout << "F[" << j << "] = " << F[j].re << " + i * "
//           << F[j].im << endl;
//    }
#   endif
# else
//    fftw orders the fourier modes as follows:
//    cos(x*0),cos(x)/2,cos(x*2)/2,...,cos(x*[n/2-1])/2,cos(x*n/2),
//      -sin(x*[n/2-1])/2,...,-sin(x*2)/2,-sin(x)/2
      for (int j=0;j<n/2;j++) {
        FHAT[j]=fhat[j]*factor;
      }
      {
        int j=n/2;
        FHAT[j]=fhat[j]*factor*0.5;
      }
      for (int j=n/2+1;j<=NPTS-n/2;j++) FHAT[j]=0;
      for (int j=1;j<n/2;j++) {
        FHAT[NPTS-j]=fhat[n-j]*factor;
      }
      rfftw_one(plan_backward,FHAT,F);
#   ifdef DEBUG
//    for (int j=0;j<NPTS;j++) {
//      cout << "F[" << j << "] = " << F[j] << endl;
//    }
#   endif
# endif
#endif

      dx=(xmax-xmin)/static_cast<double>(NPTS);
//    if (n==4) {
        x=xmin;
        gt.newPage();
        gt.setfgColor("black");
        gt.drawAxes();
        gt.setfgColor("red");
        for (int j=0;j<NPTS;j++,x+=dx) {
          double y=fcn(x);
          if (j==0) gt.movePen(x,y);
          else gt.drawLine(x,y);
        }
//    }
      double maxerr=0.;
//    gt.setfgColor(n-nmin+1,nmax-nmin);
      gt.setfgColor("blue");
      x=xmin;
      for (int j=0;j<NPTS;j++,x+=dx) {
#ifdef USE_FFTPACK
//      double y=fftpack_poly(x,n,f);
        double y=FHAT[j];
#endif
#ifdef USE_FFTW
# ifdef USE_COMPLEX
        fftw_real y=F[j].re;
# else
//      fftw_real y=fftw_poly(x,n,fhat);
        fftw_real y=F[j];
# endif
#endif
        if (j==0) gt.movePen(x,y);
        else gt.drawLine(x,y);
        maxerr=max(maxerr,abs(y-fcn(x)));
//      cout << t << " " << z << endl;
      }
      maxerr=max(maxerr,DBL_EPSILON);
      gt.setfgColor("green");
      x=xmin;
      dx=(xmax-xmin)/static_cast<double>(n);
      for (int j=0;j<n;j++,x+=dx) {
        double y=fcn(x);
        gt.drawPlus(x,y,dx*0.25);
      }
      errors[nm]=log(maxerr)/log10;
      errors_min=min(errors_min,errors[nm]);
      errors_max=max(errors_max,errors[nm]);
      gt.flush();
      cout << "errors[" << n << "] = " << maxerr << endl;
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
#ifdef USE_FFTW
      fftw_destroy_plan(plan_forward);
#endif
    }

    XYGraphTool gt2("errors","log10 2^n","log10 error",
      static_cast<double>(nmin)*log(2.)/log10,
      static_cast<double>(nmax)*log(2.)/log10,errors_min,errors_max,
      &cmap, 0,0.5);
    gt2.setbgColor("white");
    gt2.setfgColor("black");
    gt2.drawAxes();
    gt2.setfgColor("blue");
    double t=static_cast<double>(nmin)*log(2.)/log10;
    gt2.movePen(t,errors[nmin]);
    for (int n=nmin+1;n<nmax;n++) {
      t=static_cast<double>(n)*log(2.)/log10;
      gt2.drawLine(t,errors[n]);
    }
    gt2.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

#ifdef USE_FFTPACK
    delete [] f;
    delete [] fpoly;
    delete [] wsave;
    delete [] wsave2;
#endif
#ifdef USE_FFTW
    delete [] f;
    delete [] fhat;
    delete [] F;
    delete [] FHAT;
    fftw_destroy_plan(plan_backward);
#endif
    delete [] errors;
  }
}
template int std::__cmath_power<int>(int,unsigned);
