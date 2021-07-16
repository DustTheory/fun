#include <cstdlib> // for EXIT_SUCCESS
#include <iostream>
#include "TimedObject.H"
using namespace std;

extern "C" {
  void fcn1_(const int &n,double *A);
  void fcn2_(const int &n,double *A);
  void fcn3_(double &A);
  double fcn4_(const double &x);
}

int main(int /*argc*/,char** /*argv*/) {
  cout << boolalpha;
  int nmin=262144,nmax=16777216;

  { TimedObject fcn1_timing("fcn1");
    cout << "timing for fcn1" << endl;
    for (int n=nmin;n<=nmax;n*=2) {
      double *A=new double[n];

      { Timer timer(&fcn1_timing);
        fcn1_(n,A);
      }
      cout << "n = " << n << endl;
      fcn1_timing.printOn(cout);
      delete [] A;
    }
  }

  { TimedObject fcn2_timing("fcn2");
    cout << "\n\ntiming for fcn2" << endl;
    for (int n=nmin;n<=nmax;n*=2) {
      double *A=new double[n];

      { Timer timer(&fcn2_timing);
        fcn2_(n,A);
      }
      cout << "n = " << n << endl;
      fcn2_timing.printOn(cout);
      delete [] A;
    }
  }

  { TimedObject fcn3_timing("fcn3");
    cout << "\n\ntiming for fcn3" << endl;
    for (int n=nmin;n<=nmax;n*=2) {
      double *A=new double[n];

      { Timer timer(&fcn3_timing);
        for (int i=0;i<n;i++) {
  	  fcn3_(A[i]);
        }
      }
      cout << "n = " << n << endl;
      fcn3_timing.printOn(cout);
      delete [] A;
    }
  }

  { TimedObject fcn4_timing("fcn4");
    cout << "\n\ntiming for fcn4" << endl;
    for (int n=nmin;n<=nmax;n*=2) {
      double *A=new double[n];

      { Timer timer(&fcn4_timing);
        for (int i=0;i<n;i++) {
  	  A[i]=fcn4_(2.2);
        }
      }
      cout << "n = " << n << endl;
      fcn4_timing.printOn(cout);
      delete [] A;
    }
  }
  return EXIT_SUCCESS;
}
