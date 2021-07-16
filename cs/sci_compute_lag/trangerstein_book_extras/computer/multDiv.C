#include <cstdlib> // for drand48, EXIT_SUCCESS
#include <iostream>
#include "TimedObject.H"
using namespace std;

int main(int /*argc*/,char** /*argv*/) {
  cout << boolalpha;
  int nmin=262144,nmax=16777216;

  { TimedObject mult_timing("mult");
    cout << "timing for multiplication" << endl;
    double t=0.1;
    for (int n=nmin;n<=nmax;n*=2) {
      double *A=new double[n];
      for (int i=0;i<n;i++) A[i]=drand48();

      { Timer timer(&mult_timing);
        for (int i=0;i<n;i++) A[i]*=t;
      }
      cout << "n = " << n << endl;
      mult_timing.printOn(cout);
      delete [] A;
    }
  }

  { TimedObject div_timing("div");
    cout << "\n\ntiming for division" << endl;
    double t=10.;
    for (int n=nmin;n<=nmax;n*=2) {
      double *A=new double[n];
      for (int i=0;i<n;i++) A[i]=drand48();

      { Timer timer(&div_timing);
        for (int i=0;i<n;i++) A[i]/=t;
      }
      cout << "n = " << n << endl;
      div_timing.printOn(cout);
      delete [] A;
    }
  }

  { TimedObject mult2_timing("mult2");
    cout << "\n\ntiming for multiplication" << endl;
    double t=0.5;
    for (int n=nmin;n<=nmax;n*=2) {
      double *A=new double[n];
      for (int i=0;i<n;i++) A[i]=drand48();

      { Timer timer(&mult2_timing);
        for (int i=0;i<n;i++) A[i]*=t;
      }
      cout << "n = " << n << endl;
      mult2_timing.printOn(cout);
      delete [] A;
    }
  }

  { TimedObject div2_timing("div2");
    cout << "\n\ntiming for division" << endl;
    double t=2.;
    for (int n=nmin;n<=nmax;n*=2) {
      double *A=new double[n];
      for (int i=0;i<n;i++) A[i]=drand48();

      { Timer timer(&div2_timing);
        for (int i=0;i<n;i++) A[i]/=t;
      }
      cout << "n = " << n << endl;
      div2_timing.printOn(cout);
      delete [] A;
    }
  }
  return EXIT_SUCCESS;
}
