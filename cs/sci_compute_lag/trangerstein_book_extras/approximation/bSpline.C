#include <iostream>
#include <limits>
//#include <math.h> // for HUGE_VAL,M_PI
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

int main(int argc,char *argv[]) {
#ifdef MEM_DEBUG
  MemoryDebugger md(1);
#endif
  {
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif
    int order=1;
    int n=128;
    double dx=1./static_cast<double>(n);
    double *bold=OPERATOR_NEW_BRACKET(double,n+1);
    for (int i=0;i<n;i++) bold[i]=1.;
    bold[n]=0.;

    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    {
      XYGraphTool gt("cardinal B-spline","x","B_k",0.,1.,0.,1.,&cmap,0,0.5);
      gt.setbgColor("white");
      gt.newPage();
      gt.setfgColor( "black" );
      gt.drawAxes();
      gt.setfgColor( "blue" );
      gt.movePen(0.,bold[0]);
      for (int i=1;i<=n;i++) gt.drawLine(i*dx,bold[i]);
      cout << "order = " << order << endl;
      gt.flush();
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
    }

    double *bnew=0;
    while (true) {
      order++;
      bnew=OPERATOR_NEW_BRACKET(double,order*n+1);
      for (int i=0;i<n;i++) {
        double t=static_cast<double>(i)*dx;
        bnew[i]=bold[i]*t/static_cast<double>(order-1);
      }
      for (int i=n;i<=(order-1)*n;i++) {
        double t=static_cast<double>(i)*dx;
        bnew[i]=(bold[i]*t+bold[i-n]*static_cast<double>(order-t))
               /static_cast<double>(order-1);
      }
      for (int i=(order-1)*n+1;i<=order*n;i++) {
        double t=static_cast<double>(i)*dx;
        bnew[i]=bold[i-n]*static_cast<double>(order-t)
               /static_cast<double>(order-1);
      }
      XYGraphTool gt("cardinal B-spline","x","B_k",
        0.,static_cast<double>(order),0.,1.,&cmap,0,0.5);
      gt.newPage();
      gt.rescale(0.,static_cast<double>(order),0.,1.);
      gt.setfgColor("black");
      gt.drawAxes();
      gt.setfgColor("blue");
      gt.movePen(0.,bnew[0]);
      for (int i=1;i<=order*n;i++) {
        gt.drawLine(static_cast<double>(i)*dx,bnew[i]);
      }
      gt.flush();
      cout << "order = " << order << endl;
      XYGraphTool::WINDOW_TYPE::QuitButton qb;
      delete [] bold;
      bold=bnew;
    }
    delete [] bnew; bnew=0;
  }
}
