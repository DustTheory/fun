#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <sys/times.h>
#include <sys/param.h>
using namespace std;

double rate=1.;
inline double f(double u,double t) { return rate*u; }
inline double exact(double u0,double t) { return u0*exp(rate*t); }
inline double euler(double u,double t,double dt) { return u+dt*f(u,t); }
inline double modeuler(double u,double t,double dt) { 
  double dt2=0.5*dt;
  return u+dt*f(euler(u,t,dt2),t+dt2);
}

int main(int /*argc*/,char** /*argv*/) {
  int nlevels=11;
  int repeat=100000;
  int nsteps;
  double rate,tmax,u0;

//cout << "enter rate" << endl;
//cin >> rate;
//cout << "enter initial u" << endl;
//cin >> u0;
//cout << "enter tmax" << endl;
//cin >> tmax;
//cout << "enter number steps" << endl;
//cin >> nsteps;
  rate=1.;
  u0=1.;
  tmax=1.;
  nsteps=16;
  if (nsteps<=0) nsteps=1;
//cout << "rate,u0,tmax,nsteps = " << rate << " " << u0 << " " << tmax
//  << " " << nsteps << endl;

  struct tms buffer;
  ofstream out_efficiency,out_order;
  out_efficiency.open("euler_efficiency",ios::out);
  out_order.open("euler_order",ios::out);
  int ns=nsteps;
  for (int j=1;j<nlevels;j++) {
    double dt=tmax/static_cast<double>(nsteps);
    times(&buffer);
    clock_t times_start=buffer.tms_utime;
  
    double u;
    for (int r=0;r<repeat;r++) {
      double t=0.;
      u=u0;
      for (int i=0;i<nsteps;i++) {
        u=euler(u,t,dt);
        t=t+dt;
      }
    }

    times(&buffer);
    clock_t times_stop=buffer.tms_utime;
    double elapsed=static_cast<double>(times_stop-times_start)
      / ( static_cast<double>(HZ) * static_cast<double>(repeat) );
    elapsed=max(elapsed,numeric_limits<double>::epsilon());

    double error=abs(exact(u0,tmax)-u);
    out_efficiency << log10(elapsed) << " " << log10(error) << endl;
    out_order << -log10(dt) << " " << log10(error) << endl;
    cout << "j,time,error = " << j << " " << elapsed << " " << error << endl;
    nsteps*=2;
  }
  out_efficiency << endl;
  out_order << endl;

  nsteps=ns;
  for (int j=1;j<nlevels;j++) {
    double dt=tmax/static_cast<double>(nsteps);
    times(&buffer);
    clock_t times_start=buffer.tms_utime;
  
    double u;
    for (int r=0;r<repeat;r++) {
      double t=0.;
      u=u0;
      t=0.;
      for (int i=0;i<nsteps;i++) {
        u=modeuler(u,t,dt);
        t=t+dt;
      }
    }

    times(&buffer);
    clock_t times_stop=buffer.tms_utime;
    double elapsed=static_cast<double>(times_stop-times_start)
      / ( static_cast<double>(HZ) * static_cast<double>(repeat) );
    elapsed=max(elapsed,numeric_limits<double>::epsilon());

    double error=abs(exact(u0,tmax)-u);
    out_efficiency << log10(elapsed) << " " << log10(error) << endl;
    out_order << -log10(dt) << " " << log10(error) << endl;
    cout << "j,time,error = " << j << " " << elapsed << " " << error << endl;
    nsteps*=2;
  }
  out_efficiency.close();
  out_order.close();

  return EXIT_SUCCESS;
}
