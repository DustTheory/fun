#include <cmath>
#include <fstream>
#include <iostream>
#include "MemoryDebugger.H"
#include "TimedObject.H"
#include "Tracer.H"
   
extern "C" {
  void integrate_(const int &nsteps,const double &tmax,
                  double *soln,double *time); 
}  
struct odepar_common {
  double rate;
}; 
extern odepar_common odepar_;
   
#define MAX_STEPS 100000
   
int main(int /*argc*/,char** /*argv*/) {
  MemoryDebugger md(1);
  cout << boolalpha;
  { Tracer constructor_tracer("constructor block");
    int nsteps; 
    double rate,tmax,u; 
    cout << "enter rate " << endl;
    cin >> rate;
    cout << "rate = " << rate << endl;
    cout << "enter initial u" << endl;
    cin >> u;
    cout << "enter tmax" << endl;
    cin >> tmax;
    cout << "enter number steps" << endl;
    cin >> nsteps; 
    nsteps=max(0,min(MAX_STEPS,nsteps));
    cout << "number steps = " << nsteps << endl;
   
    double *t_array=OPERATOR_NEW_BRACKET(double,nsteps);
    double *u_array=OPERATOR_NEW_BRACKET(double,nsteps);
   
    t_array[0]=0.;
   
    TimedObject integrate_timing("integrate");
    { Tracer tracer("integrate block");
      Timer timer(&integrate_timing);
      integrate_(nsteps,tmax, u_array,t_array);
    }
    integrate_timing.printOn(cout);
   
    TimedObject formatted_write_timing("formatted write");
    { Timer timer(&formatted_write_timing);
      ofstream os;
      os.open("mixed_main_output");
      for (int i=0;i<=nsteps;i++) {
        os << t_array[i] << " " << u_array[i] << endl;
      }       
    }
    formatted_write_timing.printOn(cout);

    delete [] t_array; t_array=0;
    delete [] u_array; u_array=0;
  }
  return EXIT_SUCCESS;
}
