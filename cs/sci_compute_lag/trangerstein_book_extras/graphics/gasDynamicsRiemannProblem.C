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
// "$Header: /home/faculty/johnt/cvs/deal_new/graphics/gasDynamicsRiemannProblem.C,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#include <algorithm>
#include <float.h>
#include <fstream>
#ifdef USE_GTK
#include <gtk/gtkmain.h>
#endif
#include <iomanip>
#include <iostream>
#include <math.h>
//#include "Debug.H"
//#include "MemoryDebugger.H"
#include "Palette.H"
//#include "Tracer.H"
#include "VGT.H"
#include "XYGraphTool.H"

#include <unistd.h>

#define NPLOT 100
#define ITMAX 10

void plotSlowRarefaction(XYGraphTool &gt,double gamma,double cleft,
double pleft,double vleft,double pstar) {
//TRACER_CALL(t,"plotSlowRarefaction");
  double rarefaction_exponent=0.5*(1.-1./gamma);
  gt.movePen(vleft,pleft);
  gt.setfgColor("blue");
  for (int i=1;i<=NPLOT;i++) {
    double p=pleft-(pleft-pstar)*double(i)/double(NPLOT);
    double c=cleft*pow(p/pleft,rarefaction_exponent);
    double v=vleft-2.*(c-cleft)/(gamma-1.);
    gt.drawLine(v,p);
  }
}
void plotSlowShock(XYGraphTool &gt,double gamma,double cleft,
double pleft,double vleft,double pstar) {
//TRACER_CALL(t,"plotSlowShock");
  double factor=0.5*(gamma+1.)/gamma;
  gt.movePen(vleft,pleft);
  gt.setfgColor("cyan");
//double rholeft=gamma*pleft/(cleft*cleft);
//double etotleft=pleft/(gamma-1.)+0.5*rholeft*vleft*vleft;
//double efluxleft=(etotleft+pleft)*vleft;
  for (int i=1;i<=NPLOT;i++) {
    double p=pleft+(pstar-pleft)*double(i)/double(NPLOT);
    double z=(p-pleft)/pleft;
    double v=vleft-(z*cleft/gamma)/sqrt(1.+factor*z);
    gt.drawLine(v,p);

//  double rho=rholeft*(1.+z*factor)/(1.+0.5*z*(gamma-1.)/gamma);
//  cout << "\ni = " << i << endl;
//  cout << "\tmass shock speed = " 
//       << (v*rho-vleft*rholeft)/(rho-rholeft) << endl;
//  cout << "\tmomentum shock speed = " 
//       << (p+rho*v*v-pleft-rholeft*vleft*vleft)/(rho*v-rholeft*vleft)
//       << endl;
//  double etot=p/(gamma-1.)+0.5*rho*v*v;
//  double eflux=(etot+p)*v;
//  cout << "\tenergy shock speed = " 
//       << (eflux-efluxleft)/(etot-etotleft) << endl;
  }
}
void plotFastRarefaction(XYGraphTool &gt,double gamma,double cright,
double pright,double vright,double pstar) {
//TRACER_CALL(t,"plotFastRarefaction");
  double rarefaction_exponent=0.5*(1.-1./gamma);
  gt.movePen(vright,pright);
  gt.setfgColor("red");
  for (int i=1;i<=NPLOT;i++) {
    double p=pright-(pright-pstar)*double(i)/double(NPLOT);
    double c=cright*pow(p/pright,rarefaction_exponent);
    double v=vright+2.*(c-cright)/(gamma-1.);
    gt.drawLine(v,p);
  }
}
void plotFastShock(XYGraphTool &gt,double gamma,double cright,
double pright,double vright,double pstar) {
//TRACER_CALL(t,"plotFastShock");
  double factor=0.5*(gamma+1.)/gamma;
  gt.movePen(vright,pright);
  gt.setfgColor("magenta");
  double rhoright=gamma*pright/(cright*cright);
  double etotright=pright/(gamma-1.)+0.5*rhoright*vright*vright;
  double efluxright=(etotright+pright)*vright;
  for (int i=1;i<=NPLOT;i++) {
    double p=pright+(pstar-pright)*double(i)/double(NPLOT);
    double z=(p-pright)/pright;
    double v=vright+(z*cright/gamma)/sqrt(1.+factor*z);
    gt.drawLine(v,p);

//  double rho=rhoright*(1.+z*factor)/(1.+0.5*z*(gamma-1.)/gamma);
//  cout << "\ni = " << i << endl;
//  cout << "\tmass shock speed = " 
//       << (v*rho-vright*rhoright)/(rho-rhoright) << endl;
//  cout << "\tmomentum shock speed = " 
//       << (p+rho*v*v-pright-rhoright*vright*vright)
//          /(rho*v-rhoright*vright) << endl;
//  double etot=p/(gamma-1.)+0.5*rho*v*v;
//  double eflux=(etot+p)*v;
//  cout << "\tenergy shock speed = " 
//       << (eflux-efluxright)/(etot-etotright) << endl;
  }
}
void solveShockShock(XYGraphTool &gt,double gamma,double pleft,
double rholeft,double vleft,double pright,double rhoright,double vright,
double &p,double &rhol,double &rhor,double &v) {
//TRACER_CALL(t,"solveShockShock");
  double cleft=sqrt(gamma*pleft/rholeft);
  double impedleft=cleft*rholeft;
  double cright=sqrt(gamma*pright/rhoright);
  double impedright=cright*rhoright;
  double shock_factor=0.5*(gamma+1.)/gamma;

  p=(cright*pright+cleft*pleft-cleft*cright*(vright-vleft))
   /(cright+cleft);
  v=0.5*(vleft+vright);

  for (int it=0;it<ITMAX;it++) {
//  cout << "\n\tit = " << it << endl;
//  cout << "v,p = " << v << " " << p << endl;
    double zleft=p/pleft-1.;
    double wleft=impedleft*sqrt(max(1.-shock_factor,
                                    1.+shock_factor*zleft));
    double zetaleft=2.*wleft*wleft*wleft
                   /(wleft*wleft+impedleft*impedleft);
    double vl=vleft+(pleft-p)/wleft;

    double zright=p/pright-1.;
    double wright=impedright*sqrt(max(1.-shock_factor,
                                      1.+shock_factor*zright));
    double zetaright=2.*wright*wright*wright
                    /(wright*wright+impedright*impedright);
    double vr=vright-(pright-p)/wright;

    p=p-(vr-vl)*zetaright*zetaleft/(zetaright+zetaleft);
    v=(zetaright*vr+zetaleft*vl)/(zetaright+zetaleft);
  }
//cout << "\n\tshock-shock intersection:" << endl;
//cout << "v,p = " << v << " " << p << endl;
  double zleft=p/pleft-1.;
  rhol=rholeft*(1.+zleft/(gamma+0.5*zleft*(gamma-1.)));
  double zright=p/pright-1.;
  rhor=rhoright*(1.+zright/(gamma+0.5*zright*(gamma-1.)));
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.flush();
  plotSlowShock(gt,gamma,cleft,pleft,vleft,p);
  plotFastShock(gt,gamma,cright,pright,vright,p);
  gt.flush();
}
void solveRarefactionShock(XYGraphTool &gt,double gamma,double pleft,
double rholeft,double vleft,double pright,double rhoright,double vright,
double &p,double &rhol,double &rhor,double &v) {
//TRACER_CALL(t,"solveRarefactionShock");
  double cleft=sqrt(gamma*pleft/rholeft);
  double cright=sqrt(gamma*pright/rhoright);
  double impedright=cright*rhoright;
  double rarefaction_exponent=0.5*(1.-1./gamma);
  double shock_factor=0.5*(gamma+1.)/gamma;
  p=min(pleft,max(pright,p));
//cout << "\trarefaction_exponent = " << rarefaction_exponent << endl;
  for (int it=0;it<ITMAX;it++) {
//  cout << "\tit,v,p = " << it << " " << v << " " << p << endl;
    double zright=p/pright-1.;
    double wright=impedright*sqrt(1.+shock_factor*zright);
    double coverw=impedright/wright;
    double zetaright=2.*wright/(1.+coverw*coverw);
    double vr=vright-(pright-p)/wright;
    rhor=rhoright*(1.+zright/(gamma+0.5*zright*(gamma-1.)));

    double wleft=cleft*pow(p/pleft,rarefaction_exponent);
    double vl=vleft-2.*(wleft-cleft)/(gamma-1.);
//  cout << "\twleft,vl = " << wleft << " " << vl << endl;
    rhol=rholeft*pow(p/pleft,1./gamma);

    double dvldp=-1./sqrt(gamma*p*rhol);
    double dvrdp=1./zetaright;
    p=p-(vr-vl)/(dvrdp-dvldp);
    v=0.5*(vl+vr);
  }
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.flush();
  plotSlowRarefaction(gt,gamma,cleft,pleft,vleft,p);
  plotFastShock(gt,gamma,cright,pright,vright,p);
  gt.flush();
}
void solveShockRarefaction(XYGraphTool &gt,double gamma,double pleft,
double rholeft,double vleft,double pright,double rhoright,double vright,
double &p,double &rhol,double &rhor,double &v) {
//TRACER_CALL(t,"solveShockRarefaction");
//cout << "\tpleft,vleft = " << pleft << " " << vleft << endl;
//cout << "\tpright,vright = " << pright << " " << vright << endl;
  double cleft=sqrt(gamma*pleft/rholeft);
  double impedleft=cleft*rholeft;
  double cright=sqrt(gamma*pright/rhoright);
  double rarefaction_exponent=0.5*(1.-1./gamma);
  double shock_factor=0.5*(gamma+1.)/gamma;

  for (int it=0;it<ITMAX;it++) {
//  cout << "\n\tit = " << it << endl;
//  cout << "v,p = " << v << " " << p << endl;
    double zleft=p/pleft-1.;
    double wleft=impedleft*sqrt(1.+shock_factor*zleft);
    double coverw=impedleft/wleft;
    double zetaleft=2.*wleft/(1.+coverw*coverw);
    double vl=vleft+(pleft-p)/wleft;
    rhol=rholeft*(1.+zleft/(gamma+0.5*zleft*(gamma-1.)));

    double wright=cright*pow(p/pright,rarefaction_exponent);
    double vr=vright+2.*(wright-cright)/(gamma-1.);
    rhor=rhoright*pow(p/pright,1./gamma);

    double dvrdp=1./sqrt(gamma*p*rhor);
    double dvldp=-1./zetaleft;
    p=p-(vr-vl)/(dvrdp-dvldp);
    v=0.5*(vl+vr);
  }
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.flush();
  plotSlowShock(gt,gamma,cleft,pleft,vleft,p);
  plotFastRarefaction(gt,gamma,cright,pright,vright,p);
  gt.flush();
}
void solveRarefactionRarefaction(XYGraphTool &gt,double gamma,
double pleft,double rholeft,double vleft,double pright,double rhoright,
double vright,double &p,double &rhol,double &rhor,double &v) {
//TRACER_CALL(t,"solveRarefactionRarefaction");
//cout << "vleft,pleft = " << vleft << " " << pleft << endl;
//cout << "vright,pright = " << vright << " " << pright << endl;
  double cleft=sqrt(gamma*pleft/rholeft);
  double cright=sqrt(gamma*pright/rhoright);
  double rarefaction_exponent=0.5*(1.-1./gamma);

  double vlmax=vleft+2.*cleft/(gamma-1.);
  double vrmin=vright-2.*cright/(gamma-1.);
//cout << "\tvlmax,vrmin = " << vlmax << " " << vrmin << endl;
  double c=0.;
  if (vlmax>vrmin) {
    if (p<0.) p=0.5*min(pleft,pright);
    for (int it=0;it<ITMAX;it++) {
//    cout << "\tv,p = " << v << " " << p << endl;
      double wright=cright*pow(p/pright,rarefaction_exponent);
      double vr=vright+2.*(wright-cright)/(gamma-1.);
      rhor=rhoright*pow(p/pright,1./gamma);

      double wleft=cleft*pow(p/pleft,rarefaction_exponent);
      double vl=vleft-2.*(wleft-cleft)/(gamma-1.);
      rhol=rholeft*pow(p/pleft,1./gamma);

      double dvldp=-1./sqrt(gamma*p*rhol);
      double dvrdp=1./sqrt(gamma*p*rhor);
      double pstar=p-(vr-vl)/(dvrdp-dvldp);
      if (pstar<0.) p=0.5*p;
      else p=pstar;
      v=0.5*(vl+vr);
    }
  } else {
    p=0.;
    v=0.5*(vlmax+vrmin);
  }
  gt.setbgColor("white");
  gt.setfgColor("black");
  gt.drawAxes();
  gt.flush();
  plotSlowRarefaction(gt,gamma,cleft,pleft,vleft,p);
  plotFastRarefaction(gt,gamma,cright,pright,vright,p);
  gt.flush();
}
void solveRiemannProblem(XYGraphTool &gt,double gamma,double pleft,
double rholeft,double vleft,double pright,double rhoright,double vright,
double &p,double &rhol,double &rhor,double &v) {
//TRACER_CALL(t,"solveRiemannProblem");
  double rarefaction_exponent=0.5*(1.-1./gamma);
  double shock_factor=0.5*(gamma+1.)/gamma;
  pleft=max(DBL_EPSILON,pleft);
  double cleft=sqrt(gamma*pleft/rholeft);
  double impedleft=cleft*rholeft;
  pright=max(DBL_EPSILON,pright);
  double cright=sqrt(gamma*pright/rhoright);
  double impedright=cright*rhoright;
  p=(impedleft*pright+impedright*pleft
    -impedleft*impedright*(vright-vleft))/(impedright+impedleft);
  p=max(0.5*min(pleft,pright),p);
  v=(impedleft*vleft+impedright*vright-pright+pleft)
   /(impedleft+impedright);
  if (pright>pleft) {
    if (vright>vleft) {
      double c=cright*pow(pleft/pright,rarefaction_exponent);
      v=vright+2.*(c-cright)/(gamma-1.);
      if (vleft>v) {
        solveShockRarefaction(gt,gamma,pleft,rholeft,vleft,pright,
          rhoright,vright,p,rhol,rhor,v);
      } else {
        solveRarefactionRarefaction(gt,gamma,pleft,rholeft,vleft,pright,
          rhoright,vright,p,rhol,rhor,v);
      }
    } else {
      double z=pright/pleft-1.;
      double w=impedleft*sqrt(1.+z*shock_factor);
      v=vleft-(pright-pleft)/w;
      if (vright>v) {
        solveShockRarefaction(gt,gamma,pleft,rholeft,vleft,pright,
          rhoright,vright,p,rhol,rhor,v);
      } else {
        solveShockShock(gt,gamma,pleft,rholeft,vleft,pright,
          rhoright,vright,p,rhol,rhor,v);
      }
    }
  } else {
    if (vright>vleft) {
      double c=cleft*pow(pright/pleft,rarefaction_exponent);
      v=vleft-2.*(c-cleft)/(gamma-1.);
      if (vright>v) {
        solveRarefactionRarefaction(gt,gamma,pleft,rholeft,vleft,pright,
          rhoright,vright,p,rhol,rhor,v);
      } else {
        solveRarefactionShock(gt,gamma,pleft,rholeft,vleft,pright,
          rhoright,vright,p,rhol,rhor,v);
      }
    } else {
      double z=pleft/pright-1.;
      double w=impedright*sqrt(1.+z*shock_factor);
      v=vright+(pleft-pright)/w;
      if (vleft<v) {
        solveRarefactionShock(gt,gamma,pleft,rholeft,vleft,pright,
          rhoright,vright,p,rhol,rhor,v);
      } else {
        solveShockShock(gt,gamma,pleft,rholeft,vleft,pright,rhoright,
          vright,p,rhol,rhor,v);
      }
    }
  }
}

int main(int argc,char* argv[]) {
  cout << boolalpha;
//cout << "\tin main" << endl;
  {
#ifdef _MEM_CHECK_
    MemoryDebugger md(1);
#endif
#ifdef USE_GTK
    GTKWindow::gtkInit(argc,argv);
#endif

    double winsize=0.5;
    double gamma=1.4;
    double vlo=0.,vhi=1.;
    double plo=0.,phi=1.;
    double rholeft=1.,rhoright=1.;

    Palette pal;
    XYGraphTool::WINDOW_TYPE::COLOR_MAP_TYPE cmap(&pal);
    XYGraphTool gt(
      "Gas Dynamics Riemann Problem","Velocity","Pressure",vlo,vhi,
      plo,phi,&cmap,NULL,winsize);

    while (1) {
      gt.setbgColor("white");
      gt.setfgColor("black");
      gt.drawAxes();
      gt.flush();

      double vleft=vlo,pleft=plo;
      cout << "use mouse to select left state for Riemann problem"
           << endl;
      gt.getMouse(vleft,pleft);
      if (pleft<DBL_EPSILON) pleft=DBL_EPSILON;
      double cleft=sqrt(gamma*pleft/rholeft);
//    cout << "vleft,pleft,rholeft = " 
//         << vleft << " " << pleft << " " << rholeft << endl;
//    cout << "\tcleft = " << cleft << endl;
      plotSlowShock(gt,gamma,cleft,pleft,vleft,phi);
      plotSlowRarefaction(gt,gamma,cleft,pleft,vleft,plo);
      gt.flush();

      double vright=vlo,pright=plo;
      cout << "use mouse to select right state for Riemann problem"
           << endl;
      gt.getMouse(vright,pright);
      if (pright<DBL_EPSILON) pright=DBL_EPSILON;
      double cright=sqrt(gamma*pright/rhoright);
//    cout << "vright,pright,rhoright = " 
//         << vright << " " << pright << " " << rhoright << endl;
//    cout << "\tcright = " << cright << endl;

      double p=pright;
      double v=vright;
      double rhol=rholeft;
      double rhor=rhoright;
      solveRiemannProblem(gt,gamma,pleft,rholeft,vleft,pright,
        rhoright,vright,p,rhol,rhor,v);

      while (gt.buttonIsPressed(vright,pright)) {
//      TRACER_CALL(t,"buttonIsPressed");
        gt.setbgColor("white");
        gt.setfgColor("black");
        gt.drawAxes();
//      cout << "vright,pright,rhoright = " 
//           << vright << " " << pright << " " << rhoright << endl;
        solveRiemannProblem(gt,gamma,pleft,rholeft,vleft,pright,
          rhoright,vright,p,rhol,rhor,v);
      }

      if (pleft<p) {
        if (p>pright) {
//        TRACER_CALL(t,"shock-shock");
          double speedmin=v-sqrt(gamma*p/rhol);
          double speedmax=v+sqrt(gamma*p/rhor);
          double slow_shock=(v*rhol-vleft*rholeft)/(rhol-rholeft);
          double fast_shock=(v*rhor-vright*rhoright)/(rhor-rhoright);
          XYGraphTool gtv("Velocity vs. Wavespeed","Wavespeed",
            "Velocity",speedmin,speedmax,vright,vleft,&cmap,NULL,
            winsize);
          gtv.setbgColor("white");
          gtv.setfgColor("black");
          gtv.drawAxes();
          gtv.movePen(speedmin,vleft);
          gtv.setfgColor("orange");
          gtv.drawLine(slow_shock,vleft);
          gtv.setfgColor("cyan");
          gtv.drawLine(slow_shock,v);
          gtv.setfgColor("orange");
          gtv.drawLine(fast_shock,v);
          gtv.setfgColor("magenta");
          gtv.drawLine(fast_shock,vright);
          gtv.setfgColor("orange");
          gtv.drawLine(speedmax,vright);
          gtv.flush();

          XYGraphTool gtp("Pressure vs. Wavespeed","Wavespeed",
            "Pressure",speedmin,speedmax,min(pleft,pright),p,&cmap,NULL,
            winsize);
          gtp.setbgColor("white");
          gtp.setfgColor("black");
          gtp.drawAxes();
          gtp.movePen(speedmin,pleft);
          gtp.setfgColor("orange");
          gtp.drawLine(slow_shock,pleft);
          gtp.setfgColor("cyan");
          gtp.drawLine(slow_shock,p);
          gtp.setfgColor("orange");
          gtp.drawLine(fast_shock,p);
          gtp.setfgColor("magenta");
          gtp.drawLine(fast_shock,pright);
          gtp.setfgColor("orange");
          gtp.drawLine(speedmax,pright);
          gtp.flush();

          XYGraphTool gtrho("Density vs. Wavespeed","Wavespeed",
            "Density",speedmin,speedmax,min(rholeft,rhoright),
            max(rhol,rhor),&cmap,NULL,winsize);
          gtrho.setbgColor("white");
          gtrho.setfgColor("black");
          gtrho.drawAxes();
          gtrho.movePen(speedmin,rholeft);
          gtrho.setfgColor("orange");
          gtrho.drawLine(slow_shock,rholeft);
          gtrho.setfgColor("cyan");
          gtrho.drawLine(slow_shock,rhol);
          gtrho.setfgColor("orange");
          gtrho.drawLine(v,rhol);
          gtrho.setfgColor("green");
          gtrho.drawLine(v,rhor);
          gtrho.setfgColor("orange");
          gtrho.drawLine(fast_shock,rhor);
          gtrho.setfgColor("magenta");
          gtrho.drawLine(fast_shock,rhoright);
          gtrho.setfgColor("orange");
          gtrho.drawLine(speedmax,rhoright);
          gtrho.flush();

          gt.flush();
          XYGraphTool::WINDOW_TYPE::QuitButton qb;
        } else {
//        TRACER_CALL(t,"shock-rarefaction");
          double speedmin=v-sqrt(gamma*p/rhol);
          double speedmax=max(vleft-cleft,vright+cright);
          double slow_shock=(v*rhol-vleft*rholeft)/(rhol-rholeft);

          XYGraphTool gtv("Velocity vs. Wavespeed","Wavespeed",
            "Velocity",speedmin,speedmax,v,max(vleft,vright),&cmap,NULL,
            winsize);
          gtv.setbgColor("white");
          gtv.setfgColor("black");
          gtv.drawAxes();
          gtv.movePen(speedmin,vleft);
          gtv.setfgColor("orange");
          gtv.drawLine(slow_shock,vleft);
          gtv.setfgColor("cyan");
          gtv.drawLine(slow_shock,v);
          gtv.setfgColor("orange");
          double c=sqrt(gamma*p/rhor);
          gtv.drawLine(v+c,v);

          XYGraphTool gtp("Pressure vs. Wavespeed","Wavespeed",
            "Pressure",speedmin,speedmax,pleft,pright,&cmap,NULL,
            winsize);
          gtp.setbgColor("white");
          gtp.setfgColor("black");
          gtp.drawAxes();
          gtp.movePen(speedmin,pleft);
          gtp.setfgColor("orange");
          gtp.drawLine(slow_shock,pleft);
          gtp.setfgColor("cyan");
          gtp.drawLine(slow_shock,p);
          gtp.setfgColor("orange");
          gtp.drawLine(v+c,p);


          XYGraphTool gtrho("Density vs. Wavespeed","Wavespeed",
            "Density",speedmin,speedmax,min(rholeft,rhor),
            max(rhol,rhoright),&cmap,NULL,winsize);
          gtrho.setbgColor("white");
          gtrho.setfgColor("black");
          gtrho.drawAxes();
          gtrho.movePen(speedmin,rholeft);
          gtrho.setfgColor("orange");
          gtrho.drawLine(slow_shock,rholeft);
          gtrho.setfgColor("cyan");
          gtrho.drawLine(slow_shock,rhol);
          gtrho.setfgColor("orange");
          gtrho.drawLine(v,rhol);
          gtrho.setfgColor("green");
          gtrho.drawLine(v,rhor);
          gtrho.setfgColor("orange");
          gtrho.drawLine(v+c,rhor);

          gtv.setfgColor("red");
          gtp.setfgColor("red");
          gtrho.setfgColor("red");
          double rarefaction_exponent=0.5*(1.-1./gamma);
          for (int i=1;i<=NPLOT;i++) {
            double pi=pright-(pright-p)*double(NPLOT-i)/double(NPLOT);
            double ci=cright*pow(pi/pright,rarefaction_exponent);
            double vi=vright+2.*(ci-cright)/(gamma-1.);
            double rhoi=rhoright*pow(pi/pright,1./gamma);
            gtv.drawLine(vi+ci,vi);
            gtp.drawLine(vi+ci,pi);
            gtrho.drawLine(vi+ci,rhoi);
          }
          gtv.setfgColor("orange");
          gtv.drawLine(speedmax,vright);
          gtv.flush();
          gtp.setfgColor("orange");
          gtp.drawLine(speedmax,pright);
          gtp.flush();
          gtrho.setfgColor("orange");
          gtrho.drawLine(speedmax,rhoright);
          gtrho.flush();

          gt.flush();
          XYGraphTool::WINDOW_TYPE::QuitButton qb;
        }
      } else {
        if (p>pright) {
//        TRACER_CALL(t,"rarefaction-shock");
          double speedmax=
            max(v-sqrt(gamma*p/rhol),v+sqrt(gamma*p/rhor));
          double speedmin=min(vleft-cleft,vright+cright);
          double fast_shock=(v*rhor-vright*rhoright)/(rhor-rhoright);

          XYGraphTool gtv("Velocity vs. Wavespeed","Wavespeed",
            "Velocity",speedmin,speedmax,
            min(vleft,vright),v,&cmap,NULL,winsize);
          gtv.setbgColor("white");
          gtv.setfgColor("black");
          gtv.drawAxes();
          gtv.movePen(speedmin,vleft);
          gtv.setfgColor("orange");
          gtv.drawLine(vleft-cleft,vleft);

          XYGraphTool gtp("Pressure vs. Wavespeed","Wavespeed",
            "Pressure",speedmin,speedmax,pright,pleft,&cmap,NULL,
            winsize);
          gtp.setbgColor("white");
          gtp.setfgColor("black");
          gtp.drawAxes();
          gtp.movePen(speedmin,pleft);
          gtp.setfgColor("orange");
          gtp.drawLine(vleft-cleft,pleft);

          XYGraphTool gtrho("Density vs. Wavespeed","Wavespeed",
            "Density",speedmin,speedmax,min(rhol,rhoright),
            max(rholeft,rhor),&cmap,NULL,winsize);
          gtrho.setbgColor("white");
          gtrho.setfgColor("black");
          gtrho.drawAxes();
          gtrho.movePen(speedmin,rholeft);
          gtrho.setfgColor("orange");
          gtrho.drawLine(vleft-cleft,rholeft);

          gtv.setfgColor("blue");
          gtp.setfgColor("blue");
          gtrho.setfgColor("blue");
          double rarefaction_exponent=0.5*(1.-1./gamma);
          for (int i=1;i<=NPLOT;i++) {
            double pi=pleft-(pleft-p)*double(i)/double(NPLOT);
            double ci=cleft*pow(pi/pleft,rarefaction_exponent);
            double vi=vleft-2.*(ci-cleft)/(gamma-1.);
            double rhoi=rholeft*pow(pi/pleft,1./gamma);
            gtv.drawLine(vi-ci,vi);
            gtp.drawLine(vi-ci,pi);
            gtrho.drawLine(vi-ci,rhoi);
          }

          gtv.setfgColor("orange");
          gtv.drawLine(fast_shock,v);
          gtv.setfgColor("magenta");
          gtv.drawLine(fast_shock,vright);
          gtv.setfgColor("orange");
          gtv.drawLine(speedmax,vright);
          gtv.flush();

          gtp.setfgColor("orange");
          gtp.drawLine(fast_shock,p);
          gtp.setfgColor("magenta");
          gtp.drawLine(fast_shock,pright);
          gtp.setfgColor("orange");
          gtp.drawLine(speedmax,pright);
          gtp.flush();

          gtrho.setfgColor("orange");
          gtrho.drawLine(v,rhol);
          gtrho.setfgColor("green");
          gtrho.drawLine(v,rhor);
          gtrho.setfgColor("orange");
          gtrho.drawLine(fast_shock,rhor);
          gtrho.setfgColor("magenta");
          gtrho.drawLine(fast_shock,rhoright);
          gtrho.setfgColor("orange");
          gtrho.drawLine(speedmax,rhoright);
          gtrho.flush();

          gt.flush();
          XYGraphTool::WINDOW_TYPE::QuitButton qb;
        } else {
//        TRACER_CALL(t,"rarefaction-rarefaction");
          double speedmax=vright+cright;
          double speedmin=vleft-cleft;

          XYGraphTool gtv("Velocity vs. Wavespeed","Wavespeed",
            "Velocity",speedmin,speedmax,vleft,vright,&cmap,NULL,
            winsize);
          gtv.setbgColor("white");
          gtv.setfgColor("black");
          gtv.drawAxes();
          gtv.movePen(speedmin,vleft);
          gtv.setfgColor("orange");
          gtv.drawLine(vleft-cleft,vleft);

          XYGraphTool gtp("Pressure vs. Wavespeed","Wavespeed",
            "Pressure",speedmin,speedmax,p,max(pleft,pright),&cmap,NULL,
            winsize);
          gtp.setbgColor("white");
          gtp.setfgColor("black");
          gtp.drawAxes();
          gtp.movePen(speedmin,pleft);
          gtp.setfgColor("orange");
          gtp.drawLine(vleft-cleft,pleft);

          XYGraphTool gtrho("Density vs. Wavespeed","Wavespeed",
            "Density",speedmin,speedmax,min(rhol,rhor),
            max(rholeft,rhoright),&cmap,NULL,winsize);
          gtrho.setbgColor("white");
          gtrho.setfgColor("black");
          gtrho.drawAxes();
          gtrho.movePen(speedmin,rholeft);
          gtrho.setfgColor("orange");
          gtrho.drawLine(vleft-cleft,rholeft);

          gtv.setfgColor("blue");
          gtp.setfgColor("blue");
          gtrho.setfgColor("blue");
          double rarefaction_exponent=0.5*(1.-1./gamma);
          for (int i=1;i<=NPLOT;i++) {
            double pi=pleft-(pleft-p)*double(i)/double(NPLOT);
            double ci=cleft*pow(pi/pleft,rarefaction_exponent);
            double vi=vleft-2.*(ci-cleft)/(gamma-1.);
            double rhoi=rholeft*pow(pi/pleft,1./gamma);
            gtv.drawLine(vi-ci,vi);
            gtp.drawLine(vi-ci,pi);
            gtrho.drawLine(vi-ci,rhoi);
          }

          double vlmax=vleft+2.*cleft/(gamma-1.);
          double vrmin=vright-2.*cright/(gamma-1.);
          if (vlmax>vrmin) {
            gtv.setfgColor("orange");
            double c=sqrt(gamma*p/rhor);
            gtv.drawLine(v+c,v);
            gtp.setfgColor("orange");
            gtp.drawLine(v+c,p);
            gtrho.setfgColor("orange");
            gtrho.drawLine(v,rhol);
            gtrho.setfgColor("green");
            gtrho.drawLine(v,rhor);
            gtrho.setfgColor("orange");
            gtrho.drawLine(v+c,rhor);
          } else {
            gtv.setfgColor("orange");
            gtv.drawLine(vrmin,vrmin);
            gtp.setfgColor("orange");
            gtp.drawLine(vrmin,p);
            gtrho.setfgColor("orange");
            gtrho.drawLine(vrmin,0.);
          }

          gtv.setfgColor("red");
          gtp.setfgColor("red");
          gtrho.setfgColor("red");
          for (int i=1;i<=NPLOT;i++) {
            double pi=pright-(pright-p)*double(NPLOT-i)/double(NPLOT);
            double ci=cright*pow(pi/pright,rarefaction_exponent);
            double vi=vright+2.*(ci-cright)/(gamma-1.);
            double rhoi=rhoright*pow(pi/pright,1./gamma);
//          cout << "\tvi,pi = " << vi << " " << pi << endl;
            gtv.drawLine(vi+ci,vi);
            gtp.drawLine(vi+ci,pi);
            gtrho.drawLine(vi+ci,rhoi);
          }
          gtv.setfgColor("orange");
          gtv.drawLine(speedmax,vright);
          gtv.flush();
          gtp.setfgColor("orange");
          gtp.drawLine(speedmax,pright);
          gtp.flush();
          gtrho.setfgColor("orange");
          gtrho.drawLine(speedmax,rhoright);
          gtrho.flush();

          gt.flush();
          XYGraphTool::WINDOW_TYPE::QuitButton qb;
        }
      }

    } // pal,cmap,gt go out of scope here 
  } // scope of MemoryDebugger
  return EXIT_SUCCESS;
}
