#include <fstream> // for ofstream
#include <iostream> // for cin, cout
#include <limits> // for numeric_limits
#include "MemoryDebugger.H"
#include "SetTraps.H"
#include "Tracer.H"
#include "TimedObject.H"

extern "C" {
  void integrate_(const int &nsteps,const double &tmax,
                  double *soln,double *time);
}
struct odepar_common {
  double rate;
};
extern odepar_common odepar_;

int main(int /*argc*/,char** /* argv */) {
  cout << boolalpha;
  MemoryDebugger md(1);
  { Tracer constructor_block_tracer("constructor block");
    int nsteps;
    double tmax,x;
    cout << "enter rate" << endl;
    cin >> odepar_.rate;
    cout << "rate = " << odepar_.rate << endl;
    cout << "enter initial x" << endl;
    cin >> x;
    cout << "enter tmax" << endl;
    cin >> tmax;
    cout << "enter number steps" << endl;
    cin >> nsteps;
    if (nsteps<0) nsteps=0;
    cout << "number steps = " << nsteps << endl;

    double *t_array=OPERATOR_NEW_BRACKET(double,nsteps+1);
    double *x_array=OPERATOR_NEW_BRACKET(double,nsteps+1);

//  Initialize arrays to help IEEE exception handling catch unassigned
    for (int i=0;i<=nsteps;i++) {
      t_array[i]=x_array[i]=numeric_limits<double>::infinity();
    }
    t_array[0]=0.;
    x_array[0]=x;
    TimedObject integrate_timing("integrate");
    { Tracer tracer("integrate block"); 
      Timer timer(&integrate_timing);
      integrate_(nsteps,tmax, x_array,t_array);
    }
    integrate_timing.printOn(cout);

    TimedObject formatted_write_timing("formatted write");
    { Timer timer(&formatted_write_timing);
      ofstream os;
      os.open("mixed_main_output_formatted");
      for (int i=0;i<=nsteps;i++) {
        os << t_array[i] << " " << x_array[i] << endl;
      }
    }
    formatted_write_timing.printOn(cout);

    TimedObject unformatted_write_timing("unformatted write");
    { Timer timer(&unformatted_write_timing);
      ofstream os;
      os.open("mixed_main_output_unformatted");
      int nchar=(nsteps+1)*sizeof(double);
      os.write( (char*) t_array, nchar );
      os.write( (char*) x_array, nchar );
    }
    unformatted_write_timing.printOn(cout);

    delete [] t_array;
    delete [] x_array;
  }
  return EXIT_SUCCESS;
}
