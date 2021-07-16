#include <cstdlib> // drand48
#include <iostream> // for cout
#include <sys/times.h> // for function times and struct tms
#include <unistd.h> // for function sysconf
using namespace std; // otherwise, use std::cout below
extern "C" {
//C++ allows multiple functions with the same name but different arguments
//It does this by mangling the function name with the list of argument
//  types to create a unique internal function name.
//To call functions from Fortran or C, which do not mangle function names,
//  we declare the functions to C++ inside the extern "C" block:
  void daxpy_(const int &n,const double &da,const double *dx,
    const int &incx,double *dy,const int &incy);
  double ddot_(const int &n,const double *dx,const int &incx,
    const double *dy,const int &incy);
  void dgemv_(const char &trans,const int &m,const int &n,
    const double &alpha,const double *A,const int &lda,const double *x,
    const int &incx,const double &beta,double *y,const int &incy);
}
int main(int argc,char **argv) {
// argc is number of words in the command line, and
// argv is an array of character strings of the words on the command line


  int m=512;
  int n=m;
  double x[n],y[m];
// other data types:
//   float complex c;
//   double complex z;
  double A[n][m]; // A[column][row]
  double *A_dynamic=new double[m*n]; // dynamic memory allocation
//we recommend dynamic memory allocation for large arrays

//fill array A with random numbers
  for (int j=0;j<n;j++) {
    for (int i=0;i<m;i++) {
      A[j][i]=drand48(); // in stdlib.h
      A_dynamic[i+j*m]=drand48();
    }
    x[j]=drand48();
  }
  for (int i=0;i<m;i++) y[i]=0.;

//times(&usage) will return
// usage.tms_utime = User CPU time
// usage.tms_stime = System CPU time
// usage.tms_cutime = User CPU time of dead children
// usage.tms_cstime = System CPU time of dead children
// all are measured in CLK_TCKths of a second
  struct tms usage; // in sys/time.h
  clock_t start; // clock_t in sys/time.h
  times(&usage);
  start=usage.tms_utime;
  for (int k=0;k<4096;k++) {
    for (int j=0;j<n;j++) {
      double xj=x[j];
      double *Aj=A[j]; // save array access time
      for (int i=0;i<m;i++) {
        y[i]+=xj*Aj[i];
//      alternatively,
//      y[i] += x[j] * A[j][i];
//      or
//      y[i] += A_dynamic[i+j*m] * x[j];
      }
    }
  }
  times(&usage);
  double elapsed=static_cast<double>(usage.tms_utime-start)
                /(static_cast<double>(sysconf(_SC_CLK_TCK))*4096.);
  cout << "time for double loop j over i = " << elapsed << endl;

  times(&usage);
  start=usage.tms_utime;
  for (int k=0;k<4096;k++) {
    for (int j=0;j<n;j++) {
      daxpy_(m,x[j],&A[j][0],1,y,1);
//    note that we did not need to copy the j'th column of A to a vector */
//    alternatively, we could use A_dynamic:
//    daxpy_(m,x[j],&A_dynamic[j*m],1,y,1);
    }
  }
  times(&usage);
  elapsed=static_cast<double>(usage.tms_utime-start)
         /(static_cast<double>(sysconf(_SC_CLK_TCK))*4096.);
  cout << "time for loop j over daxpy = " << elapsed << endl;

  times(&usage);
  start=usage.tms_utime;
  for (int k=0;k<4096;k++) {
    for (int i=0;i<m;i++) {
      y[i]=ddot_(n,&A[0][i],m,x,1);
//    note that we did not need to copy the i'th row of A to a vector */
//    alternatively, we could use A_dynamic:
//    y[i]=ddot_(n,&A_dynamic[i],m,x,incx);
    }
  }
  times(&usage);
  elapsed=static_cast<double>(usage.tms_utime-start)
         /(static_cast<double>(sysconf(_SC_CLK_TCK))*4096.);
  cout << "time for loop i over ddot = " << elapsed << endl;

  times(&usage);
  start=usage.tms_utime;
  for (int k=0;k<4096;k++) {
    dgemv_('N',m,n,1.,A[0],m,x,1,0.,y,1);
//  alternatively, we could use A_dynamic:
//  dgemv_('N',m,n,1.,A_dynamic,m,x,1,0.,y,1);
  }
  times(&usage);
  elapsed=static_cast<double>(usage.tms_utime-start)
         /(static_cast<double>(sysconf(_SC_CLK_TCK))*4096.);
  cout << "time for dgemv = " << elapsed << endl;

  delete [] A_dynamic; A_dynamic=0; /* deallocation */
}
