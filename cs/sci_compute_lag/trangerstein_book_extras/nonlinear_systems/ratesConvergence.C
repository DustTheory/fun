#include <float.h>
#include <fstream>
#include <limits>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_roots.h>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
#include <iostream>
//#include <math.h>
//#include <unistd.h>

#include "MemoryDebugger.H"
#include "Palette.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "Tracer.H"
#include "VGT.H"
#include "XYGraphTool.H"

int main(int argc,char *argv[]) {
  cout << boolalpha;
  { MemoryDebugger md(1);
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

    int nlinear=505;
    double w[nlinear + 1],xlinear[nlinear+1],ylinear[nlinear+1];
    double zold = 1.;
    double znew = 0.5;
    w[ 0 ] = 0;
    xlinear[ 0 ] = - log( zold );
    ylinear[ 0 ] = - log( znew );
    double xmin = xlinear[ 0 ];
    double ymin = ylinear[ 0 ];
    double xmax = xmin;
    double ymax = ymin;
    for (int k = 1; k <= nlinear; k++) {
      zold = znew;
      znew = zold * 0.5;
      w[ k ] = static_cast<double>(k);
      xlinear[ k ] = ylinear[ k - 1 ];
      ylinear[ k ] = - log( znew );
      xmin = min( xmin, xlinear[ k ] );
      ymin = min( ymin, ylinear[ k ] );
      xmax = max( xmax, xlinear[ k ] );
      ymax = max( ymax, ylinear[ k ] );
    }

    int nsuper=97;
    double xsuper[nsuper+1],ysuper[nsuper+1];
    zold = 2.;
    znew = 2.;
    xsuper[ 0 ] = - log( zold );
    ysuper[ 0 ] = - log( znew );
    xmin = min( xmin, xsuper[ 0 ] );
    ymin = min( ymin, ysuper[ 0 ] );
    xmax = max( xmax, xsuper[ 0 ] );
    ymax = max( ymax, ysuper[ 0 ] );
    for (int k = 1; k <= nsuper; k++) {
      zold = znew;
      znew = zold / static_cast<double>(k+1);
      xsuper[ k ] = ysuper[ k - 1 ];
      ysuper[ k ] = - log( znew );
      xmin = min( xmin, xsuper[ k ] );
      ymin = min( ymin, ysuper[ k ] );
      xmax = max( xmax, xsuper[ k ] );
      ymax = max( ymax, ysuper[ k ] );
    }

    int nquad=8;
    double xquad[nquad+1],yquad[nquad+1];
    zold = 0.5;
    znew = 0.25;
    xquad[ 0 ] = - log( zold );
    yquad[ 0 ] = - log( znew );
    xmin = min( xmin, xquad[ 0 ] );
    ymin = min( ymin, yquad[ 0 ] );
    xmax = max( xmax, xquad[ 0 ] );
    ymax = max( ymax, yquad[ 0 ] );
    for (int k = 1; k <= nquad; k++) {
      zold = znew;
      znew = zold * zold;
      xquad[ k ] = yquad[ k - 1 ];
      yquad[ k ] = - log( znew );
      xmin = min( xmin, xquad[ k ] );
      ymin = min( ymin, yquad[ k ] );
      xmax = max( xmax, xquad[ k ] );
      ymax = max( ymax, yquad[ k ] );
    }
#ifdef DEBUG
//  for (int k = 0; k <= nlinear; k++ ) {
//    cout << "k,w,xlinear,ylinear = " << k << " " << w[k] << " " << xlinear[k] << " "
//         << ylinear[k] << endl;
//  }
//  cout << "xmin,xmax = " << xmin << " " << xmax << endl;
//  cout << "ymin,ymax = " << ymin << " " << ymax << endl;
#endif

    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    double winsize=0.5;
    XYGraphTool gt("-log(x_k - x) vs -log(x_{k-1} - x)",
      "-log(x_{k-1} - x)","-log(x_k - x)",xmin,xmax, ymin,ymax,&cmap,0,
      winsize);
    gt.newPage();
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();
    gt.setfgColor("red");
    gt.movePen( xlinear[0], ylinear[0] );
    for (int k = 1; k <= nlinear; k ++ ) {
      gt.drawLine( xlinear[k], ylinear[k] );
    }
    gt.setfgColor("green");
    gt.movePen( xsuper[0], ysuper[0] );
    for (int k = 1; k <= nsuper; k ++ ) {
      gt.drawLine( xsuper[k], ysuper[k] );
    }
    gt.setfgColor("blue");
    gt.movePen( xquad[0], yquad[0] );
    for (int k = 1; k <= nquad; k ++ ) {
      gt.drawLine( xquad[k], yquad[k] );
    }
    gt.flush();

    XYGraphTool gt2("log(x_k - x) vs k","k","log(x_k-x)",
      w[1],w[nquad],-ymax,-ymin,&cmap,0,winsize);
    gt2.newPage();
    gt2.setbgColor("white");
    gt2.setfgColor("black");
    gt2.drawAxes();
    gt2.setfgColor("red");
    gt2.movePen( w[1], -ylinear[1] );
    for (int k = 2; k <= nquad; k ++ ) {
      gt2.drawLine( w[k], -ylinear[k] );
    }
    gt2.setfgColor("green");
    gt2.movePen( w[1], -ysuper[1] );
    for (int k = 2; k <= nquad; k ++ ) {
      gt2.drawLine( w[k], -ysuper[k] );
    }
    gt2.setfgColor("blue");
    gt2.movePen( w[1], -yquad[1] );
    for (int k = 2; k <= nquad; k ++ ) {
      gt2.drawLine( w[k], -yquad[k] );
    }
    gt2.flush();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;
  } // pal,cmap,gt go out of scope here 
  return EXIT_SUCCESS;
}
