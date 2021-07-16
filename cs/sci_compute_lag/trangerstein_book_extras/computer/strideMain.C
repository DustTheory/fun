#include <cstdlib> // for EXIT_SUCCESS
#include <iostream> // for cout
#include "PermutationGenerator.H"
#include "TimedObject.H"
using namespace std;

extern "C" {
  void stride1_(const int &m,const int &n,double *A);
  void stridem_(const int &m,const int &n,double *A);
  void strider_(const int &m,const int &n,const int *colperm,
    const int *rowperm,double *A);
}

int main(int /*argc*/,char** /*argv*/) {
  cout << boolalpha;
  { TimedObject stride1_timing("stride1");
    cout << "timing for stride1" << endl;
    for (int n=128;n<=8192;n*=2) {
      double *A=new double[n*n];

      { Timer timer(&stride1_timing);
        stride1_(n,n,A);
      }
      cout << "n = " << n << endl;
      stride1_timing.printOn(cout);
      delete [] A;
    }
  }

  { TimedObject stridem_timing("stridem");
    cout << "\n\ntiming for stridem" << endl;
    for (int n=128;n<=8192;n*=2) {
      double *A=new double[n*n];

      { Timer timer(&stridem_timing);
        stridem_(n,n,A);
      }
      cout << "n = " << n << endl;
      stridem_timing.printOn(cout);
      delete [] A;
    }
  }

  { TimedObject strider_timing("strider");
    cout << "\n\ntiming for strider" << endl;
    for (int n=128;n<=8192;n*=2) {
      NumPtr<int> column_perm(n);
      RandomPermutation(column_perm);
      NumPtr<int> row_perm(n);
      RandomPermutation(row_perm);
      double *A=new double[n*n];

      { Timer timer(&strider_timing);
        strider_(n,n,column_perm.getData(),row_perm.getData(),A);
      }
      cout << "n = " << n << endl;
      strider_timing.printOn(cout);
      delete [] A;
    }
  }
  return EXIT_SUCCESS;
}

#include "NumPtr.C"
INSTANTIATE_NUMPTR(int);
