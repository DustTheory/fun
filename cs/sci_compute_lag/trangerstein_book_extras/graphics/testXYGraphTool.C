//**********************************************************************
// Copyright 2006 John A. Trangenstein
//
// This software is made available for research and instructional use 
// only. 
// You may copy and use this software without charge for these 
// non-commercial purposes, provided that the copyright notice and 
// associated text is reproduced on all copies.  
// For all other uses (including distribution of modified versions), 
// please contact the author at
//   John A. Trangenstein
//   Department of Mathematics
//   Duke University
//   Durham, NC 27708-0320
//   USA
// or
//   johnt@math.duke.edu
// 
// This software is made available "as is" without any assurance that it
// is completely correct, or that it will work for your purposes.  
// Use the software at your own risk.
//**********************************************************************
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/testXYGraphTool.C,v 1.1 2009/08/20 17:31:49 johnt Exp $"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <gtk/gtk.h>
#include "MemoryDebugger.H"
#include "Palette.H"
//#include "Tracer.H"
#include "XYGraphTool.H"
#include <unistd.h>

int main(int argc,char* argv[]) {
  cout << boolalpha;
  {
#ifdef MEM_DEBUG
//  MemoryDebugger md(1);
#endif
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

//  rate in ode dx/dt = rate * x
    double rate=1.;

//  max time and timesteps
    double tmax=2.;
    int nsteps = 10;

//  window size, as fraction of min(screen width, screen height)
    double winsize=0.5;

//  time and solution arrays
    double *t_array=OPERATOR_NEW_BRACKET(double,nsteps+1);
    double *x_array=OPERATOR_NEW_BRACKET(double,nsteps+1);

//  initial condition
    t_array[0]=0.;  
    x_array[0]=1.;

    double dt=tmax/double(nsteps);
    for (int i=1;i<=nsteps;i++) {
      x_array[i] = x_array[i-1]+dt*rate*x_array[i-1]; // forward euler
      t_array[i] = t_array[i-1]+dt;
    }

//find min,max data values
    double tlo=t_array[0];
    double thi=t_array[nsteps];
    double xlo=x_array[0];
    double xhi=x_array[0];
    for (int i=1;i<=nsteps;i++) {
      double x=x_array[i];
      if (x<xlo) xlo=x;
      if (x>xhi) xhi=x;
    }

//open graphics output window
    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    XYGraphTool gt("Forward Euler for dx/dt=rate*x","t","x",0.,1.,0.,1.,
      &cmap,0,winsize);

//initialize graphics display
    gt.newPage();
    gt.rescale(tlo,thi,xlo,xhi); // XYGraphTool constructor could set
    gt.setbgColor("white");
    gt.setfgColor("black");
    gt.drawAxes();

//draw numerical solution
    gt.setfgColor("blue");
    gt.movePen(t_array[0],x_array[0]);
    for (int i=1;i<=nsteps;i++) {
      gt.drawLine(t_array[i],x_array[i]);
    }
    gt.putString((t_array[0]+t_array[nsteps])*0.5,
      0.9*x_array[0]+0.1*x_array[nsteps],"numerical solution");

//draw analytical solution
    gt.setfgColor("red");
    gt.movePen(t_array[0],x_array[0]);
    for (int i=1;i<=nsteps;i++) {
      gt.drawLine(t_array[i],x_array[0]*exp(rate * t_array[i]));
    }
    gt.putString((t_array[0]+t_array[nsteps])*0.5,
                  x_array[nsteps],"analytical solution");

//force X to perform the requests
    gt.flush();

//  pause so that we can view the results
//  gt.eventLoop();
    XYGraphTool::WINDOW_TYPE::QuitButton qb;

    delete [] t_array;
    delete [] x_array;

  }
  return EXIT_SUCCESS;
}
