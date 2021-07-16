#include <cstdlib> // for EXIT_SUCCESS
#include <cstring> // for memcpy
#include <iostream> // for cout
#include <vector> // for vector
#include "TimedObject.H"
using namespace std;

extern "C" {
  void dcopy_(const int &n,const double *dx,const int &incx,
    double *dy,const int &incy);
  void dlacpy_(const char &uplo,const int &m,const int &n,const double *A,
    const int &lda,double *B,const int &ldb);
  void f90copy_(const int &m,const int &n,const double *A,double *B);
}

int main(int /*argc*/,char** /*argv*/) {
  cout << boolalpha;
  int repeat=64;
  {
    for (int n=512;n<=8192;n*=2) {
      double *A=new double[n*n];
      double *B=new double[n*n];
      for (int i=0;i<n*n;i++) A[i]=static_cast<double>(i);
      TimedObject memcpy_timing("memcpy");
      {
        Timer timer(&memcpy_timing);
        for (int r=0;r<repeat;r++) {
          memcpy(B,A,n*n*sizeof(double));
        }
      }
      cout << "n = " << n << endl;
      memcpy_timing.printOn(cout);
      delete [] A;
      delete [] B;
    }
  }
  { 
    cout << endl;
    for (int n=512;n<=8192;n*=2) {
      double *A=new double[n*n];
      double *B=new double[n*n];
      for (int i=0;i<n*n;i++) A[i]=static_cast<double>(i);
      TimedObject dcopy_timing("dcopy");
      {
        Timer timer(&dcopy_timing);
        for (int r=0;r<repeat;r++) {
          dcopy_(n*n,A,1,B,1);
        }
      }
      cout << "n = " << n << endl;
      dcopy_timing.printOn(cout);
      delete [] A;
      delete [] B;
    }
  }
  { 
    cout << endl;
    for (int n=512;n<=8192;n*=2) {
      double *A=new double[n*n];
      double *B=new double[n*n];
      for (int i=0;i<n*n;i++) A[i]=static_cast<double>(i);
      TimedObject dlacpy_timing("dlacpy");
      {
        Timer timer(&dlacpy_timing);
        for (int r=0;r<repeat;r++) {
          dlacpy_('A',n,n,A,n,B,n);
        }
      }
      cout << "n = " << n << endl;
      dlacpy_timing.printOn(cout);
      delete [] A;
      delete [] B;
    }
  }
  { 
    cout << endl;
    for (int n=512;n<=8192;n*=2) {
      vector<double> A(n*n);
      vector<double> B(n*n);
      for (int i=0;i<n*n;i++) A[i]=static_cast<double>(i);
      TimedObject stl_copy_timing("STL copy");
      {
        Timer timer(&stl_copy_timing);
        for (int r=0;r<repeat;r++) {
          copy(A.begin(),A.end(),B.begin());
        }
      }
      cout << "n = " << n << endl;
      stl_copy_timing.printOn(cout);
    }
  }
  { 
    cout << endl;
    for (int n=512;n<=8192;n*=2) {
      double *A=new double[n*n];
      double *B=new double[n*n];
      for (int i=0;i<n*n;i++) A[i]=static_cast<double>(i);
      TimedObject f90_copy_timing("F90 copy");
      {
        Timer timer(&f90_copy_timing);
        for (int r=0;r<repeat;r++) {
          f90copy_(n,n,A,B);
        }
      }
      cout << "n = " << n << endl;
      f90_copy_timing.printOn(cout);
      delete [] A;
      delete [] B;
    }
  }
  return EXIT_SUCCESS;
}
