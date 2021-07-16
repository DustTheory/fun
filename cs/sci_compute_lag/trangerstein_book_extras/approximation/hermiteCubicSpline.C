#include <iostream>
#include <math.h> // for HUGE_VAL,M_PI
#include <stdlib.h>

using namespace std;

#include "Arch.H"
#include "Debug.H"
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "TimedObject.H"
#include "XColormap.H"
#include "XYGraphTool.H"

//#define NSPLINE 20
#define NSPLINE 5
struct spline_common {
  double potential_min_spline,potential_max_spline,dpotential_spline,
       potential_spline[NSPLINE],
       ak1_spline[NSPLINE],ak1_derivative_spline[NSPLINE],
       bk1_spline[NSPLINE],bk1_derivative_spline[NSPLINE];
};
extern spline_common F77_NAME(spline_common);

struct time_independent_k_common {
  double time_independent_k_mult,time_independent_k_reversal_potential,
       expmultak1,expshiftak1,multak1,
       expmultbk1d,expmultbk1n1,expmultbk1n2,
         expshiftbk1d,expshiftbk1n1,expshiftbk1n2,
         multbk1n1;
};
extern time_independent_k_common F77_NAME(time_independent_k_common);

extern "C" {
  void F77_NAME(spline_setup)();
  void F77_NAME(spline_values)(const double&,double&,double&);
  void F77_NAME(spline_derivative_values)(const double&,double&,
    double&);
  double F77_NAME(time_independent_potassium_ak1)(const double&);
  double F77_NAME(time_independent_potassium_ak1_derivative)(
    const double&);
  double F77_NAME(time_independent_potassium_bk1)(const double&);
  double F77_NAME(time_independent_potassium_bk1_derivative)(
    const double&);
}

int main(int /*argc*/,char** /*argv*/) {
#ifdef DEBUG
  setTraps();
#endif
  F77_NAME(spline_common).potential_min_spline=-100.;
  F77_NAME(spline_common).potential_max_spline=50.;
  F77_NAME(time_independent_k_common).time_independent_k_mult=0.6047;
  F77_NAME(time_independent_k_common).
    time_independent_k_reversal_potential=-87.26;
  F77_NAME(time_independent_k_common).expmultak1=.2385;
  F77_NAME(time_independent_k_common).expshiftak1=-59.215;
  F77_NAME(time_independent_k_common).multak1=1.02;
  F77_NAME(time_independent_k_common).expmultbk1d=-.5143;
  F77_NAME(time_independent_k_common).expmultbk1n1=.08032;
  F77_NAME(time_independent_k_common).expmultbk1n2=.06175;
  F77_NAME(time_independent_k_common).expshiftbk1d=4.753;
  F77_NAME(time_independent_k_common).expshiftbk1n1=5.476;
  F77_NAME(time_independent_k_common).expshiftbk1n2=-594.31;
  F77_NAME(time_independent_k_common).multbk1n1=.49124;

  F77_NAME(spline_setup)();
  double vmin=F77_NAME(spline_common).potential_min_spline;
  double vmax=F77_NAME(spline_common).potential_max_spline;
  double ak1min=HUGE_VAL;
  double ak1max=-HUGE_VAL;
  double bk1min=HUGE_VAL;
  double bk1max=-HUGE_VAL;

  double dv=F77_NAME(spline_common).dpotential_spline/128.;
  double vm=vmin;
  for (int i=0;i<=NSPLINE*128;i++,vm+=dv) {
    double ak1=F77_NAME(time_independent_potassium_ak1)(vm);
    ak1min=min(ak1min,ak1);
    ak1max=max(ak1max,ak1);
    double bk1=F77_NAME(time_independent_potassium_bk1)(vm);
    bk1min=min(bk1min,bk1);
    bk1max=max(bk1max,bk1);
  }
  Palette pal;
  XColormap cmap(&pal);
  XYGraphTool gtak1("ak1",vmin,vmax,ak1min,ak1max,&cmap,0,0.5);
  gtak1.setbgColor("white");
  gtak1.newPage();
  gtak1.setfgColor("black");
  gtak1.drawAxes();
  gtak1.setfgColor("blue");

  XYGraphTool gtbk1("bk1",vmin,vmax,bk1min,bk1max,&cmap,0,0.5);
  gtbk1.setbgColor("white");
  gtbk1.newPage();
  gtbk1.setfgColor("black");
  gtbk1.drawAxes();
  gtbk1.setfgColor("blue");

  vm=vmin;
  gtak1.movePen(vm,F77_NAME(time_independent_potassium_ak1)(vm));
  gtbk1.movePen(vm,F77_NAME(time_independent_potassium_bk1)(vm));
  vm+=dv;
  for (int i=1;i<=NSPLINE*128;i++,vm+=dv) {
    gtak1.drawLine(vm,F77_NAME(time_independent_potassium_ak1)(vm));
    gtbk1.drawLine(vm,F77_NAME(time_independent_potassium_bk1)(vm));
  }

  gtak1.setfgColor("red");
  gtbk1.setfgColor("red");
  vm=vmin;
  double ak1=HUGE_VAL;
  double bk1=HUGE_VAL;
  F77_NAME(spline_values)(vm, ak1,bk1);
  gtak1.movePen(vm,ak1);
  gtbk1.movePen(vm,bk1);
  vm+=dv;
  for (int i=1;i<=NSPLINE*128;i++,vm+=dv) {
    F77_NAME(spline_values)(vm, ak1,bk1);
    gtak1.drawLine(vm,ak1);
    gtbk1.drawLine(vm,bk1);
  }
  gtak1.flush();
  gtbk1.flush();

  TimedObject tf("Function evaluation");
  { Timer zf(&tf);
    for (int j=0;j<1000;j++) {
      vm=vmin;
      for (int i=0;i<=NSPLINE*128;i++,vm+=dv) {
        double ak1=F77_NAME(time_independent_potassium_ak1)(vm);
        double bk1=F77_NAME(time_independent_potassium_bk1)(vm);
      }
    }
  }

  TimedObject ts("Spline evaluation");
  { Timer zs(&ts);
    for (int j=0;j<1000;j++) {
      vm=vmin;
      for (int i=0;i<=NSPLINE*128;i++,vm+=dv) {
        F77_NAME(spline_values)(vm, ak1,bk1);
      }
    }
  }

  cout << "Function evaluation took " << tf.totalRunTime() 
       << " seconds" << endl;
  cout << "Spline evaluation took " << ts.totalRunTime() 
       << " seconds" << endl;
  wait();

}
