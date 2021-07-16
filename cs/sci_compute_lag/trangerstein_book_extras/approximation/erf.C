#include <iostream>
#include <limits>
#include <math.h> // for M_PI
#include <stdlib.h>

using namespace std;

//#include "Debug.H"
#include "MemoryDebugger.H"
//#include "SetTraps.H"
//#include "TimedObject.H"
#include "GTKColormap.H"
#include "XYGraphTool.H"

double computePoints(int n,int w,double xmax,double &ymin,double &ymax,
double *x,double *e,double *f) {
  double dx=xmax / w;
  double xi=0.;
  ymin=numeric_limits<double>::max();
  ymax=-numeric_limits<double>::max();
  for (int i=0;i<=w;i++,xi+=dx) {
    double sum=0.;
    double sign=1.;
    double x2=xi*xi;
    double xpower=xi;
    double factorial=1.;
    for (int j=0;j<=n;j++) {
      sum+=sign*xpower/factorial;
      sign=-sign;
      xpower*=x2;
      factorial*=static_cast<double>((2*j+2)*(2*j+3));
    }
    x[i]=xi;
    e[i]=erf(xi);
    f[i]=sum*2./sqrt(M_PI);
    ymin=min(ymin,min(e[i],f[i]));
    ymax=max(ymax,max(e[i],f[i]));
  }
}

int main(int argc,char *argv[]) {
#ifdef MEM_DEBUG
  MemoryDebugger md(1);
#endif
#ifdef USE_GTK
  GTKWindow::gtkInit(argc,argv);
#endif
  int n=4;
  int w=512;
  double xmax=4.;
  double *x=OPERATOR_NEW_BRACKET(double,w+1);
  double *e=OPERATOR_NEW_BRACKET(double,w+1);
  double *f=OPERATOR_NEW_BRACKET(double,w+1);
  double ymin=numeric_limits<double>::max();
  double ymax=-numeric_limits<double>::max();
  computePoints(n,w,xmax,ymin,ymax,x,e,f);

  Palette pal;
  XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
  XYGraphTool gt("error function","x","y",0.,xmax,ymin,ymax,&cmap,0,0.5);
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.setfgColor("red");
  gt.movePen(x[0],e[0]);
  for (int i=1;i<=w;i++) {
    gt.drawLine(x[i],e[i]);
  }
  gt.setfgColor("blue");
  gt.movePen(x[0],f[0]);
  for (int i=1;i<=w;i++) {
    gt.drawLine(x[i],f[i]);
  }
  XYGraphTool::WINDOW_TYPE::QuitButton qb;

  delete [] x;
  delete [] e;
  delete [] f;
}
