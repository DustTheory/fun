#include <cstdlib> // for EXIT_SUCCESS
#include <iostream> // for cout
#include <cstdlib> // for drand48
#include "TimedObject.H" // for TimedObject, Timer
using namespace std;

extern "C" {
  void dgemv_simple__(const int &m,const int &n,double *A,const int &lda,
    const double *x,double *y);
  void dgemv_reversed__(const int &m,const int &n,double *A,const int &lda,
    const double *x,double *y);
  void dgemv_daxpy__(const int &m,const int &n,double *A,const int &lda,
    const double *x,double *y);
  void dgemv_ddot__(const int &m,const int &n,double *A,const int &lda,
    const double *x,double *y);
  void dgemv_(const char&,const int &m,const int &n,const double &alpha,
    double *A,const int &lda,const double *x,const int &incx,
    const double &beta,double *y,const int &incy);
  void daxpy_(const int &n,const double &a,const double *x,const int &incx,
    double *y,const int &incy);
  double ddot_(const int &n,const double *x,const int &incx,
    double *y,const int &incy);
}

int main(int /*argc*/,char** /*argv*/) {
  int n=4096;
//int n=512;
//double AA[3][2]; // 3 arrays of length 2 stored row-wise
//for (int i=0;i<3;i++) {
//  for (int j=0;j<2;j++) {
//    cout << "&AA[" << i << "][" << j << "] = " << &AA[i][j] << endl;
//  }
//}
  double *A=new double[n*n];
  double *x=new double[n];
  double *y=new double[n];
  for (int i=0;i<n*n;i++) A[i]=drand48();
  for (int i=0;i<n;i++) {
    x[i]=drand48();
  }
//
  for (int i=0;i<n;i++) y[i]=0.;
  {
    TimedObject ti("dgemv_simple");
    { Timer timer(&ti);
      dgemv_simple__(n,n,A,n,x,y);
    }
    ti.printOn(cout);
    for (int i=0;i<n;i++) y[i]=0.;
    { Timer timer(&ti);
      for (int j=0;j<n;j++) {
        double xj=x[j];
        double *colj=A+j*n;
        for (int i=0;i<n;i++) y[i]+=colj[i]*xj;
      }
    }
    ti.printOn(cout);
  }
//
//
  for (int i=0;i<n;i++) y[i]=0.;
  {
    TimedObject ti("dgemv_reversed");
    { Timer timer(&ti);
      dgemv_reversed__(n,n,A,n,x,y);
    }
    ti.printOn(cout);
    for (int i=0;i<n;i++) y[i]=0.;
    { Timer timer(&ti);
      for (int i=0;i<n;i++) {
        double sum=0.;
        double *rowi=A+i;
        for (int j=0;j<n;j++) sum+=rowi[j*n]*x[j];
        y[i]+=sum;
      }
    }
    ti.printOn(cout);
  }
//
/*
  for (int i=0;i<n;i++) y[i]=0.;
  {
    TimedObject ti("dgemv_daxpy");
    { Timer timer(&ti);
      dgemv_daxpy__(n,n,A,n,x,y);
    }
    ti.printOn(cout);
    for (int i=0;i<n;i++) y[i]=0.;
    { Timer timer(&ti);
      for (int j=0;j<n;j++) {
        daxpy_(n,x[j],A+j*n,1,y,1);
      }
    }
    ti.printOn(cout);
  }
  for (int i=0;i<n;i++) y[i]=0.;
  {
    TimedObject ti("dgemv_ddot");
    { Timer timer(&ti);
      dgemv_ddot__(n,n,A,n,x,y);
    }
    ti.printOn(cout);
    for (int i=0;i<n;i++) y[i]=0.;
    { Timer timer(&ti);
      for (int i=0;i<n;i++) {
        y[i]+=ddot_(n,A+i,n,x,1);
      }
    }
    ti.printOn(cout);
  }
  for (int i=0;i<n;i++) y[i]=0.;
  {
    TimedObject ti("dgemv");
    { Timer timer(&ti);
      dgemv_('N',n,n,1.,A,n,x,1,0.,y,1);
    }
    ti.printOn(cout);
  }
*/
  delete [] A;
  delete [] x;
  delete [] y;

  return EXIT_SUCCESS;
}
