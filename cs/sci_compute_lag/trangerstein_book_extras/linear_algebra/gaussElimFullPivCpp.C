#include <cstdlib> // drand48
#include <cstring> // memcpy
#include <iostream> // for cout

#include "GaussianFactorization.H"
#include "SquareMatrix.H"
using namespace std;
extern "C" {
  void dgemv_(const char&,const int &m,const int &n,const double &alpha,
    const double *A,const int &lda,const double *x,const int &incx,
    const double &beta,double *y,const int &incy);
  void dgesc2_(const int &n,double *A,const int &lda,double *rhs,
    int *ipiv,int *jpiv,double &scale);
  void dgetc2_(const int &n,double *A,const int &lda,int *ipiv,int *jpiv,
    int &info);
}
int main(int argc,char **argv) {
  {
    cout << "\ncalling LAPACK directly:" << endl;
    int n=10;
    double A[n][n];
    for (int j=0;j<n;j++) {
      for (int i=0;i<=j;i++) A[j][i]=drand48();
    }

    double LU[n][n];
    memcpy(LU[0],A[0],n*n*sizeof(double));
    int info,ipiv[n],jpiv[n];
    dgetc2_(n,LU[0],n,ipiv,jpiv,info);
    cout << "\tinfo = " << info << endl;
    cout << "\tipiv = ";
    for (int i=0;i<n-1;i++) cout << ipiv[i] << " ";
    cout << endl;
    cout << "\tjpiv = ";
    for (int j=0;j<n-1;j++) cout << jpiv[j] << " ";
    cout << endl;

    double b[n];
    for (int i=0;i<n;i++) b[i]=drand48();

    double x[n];
    memcpy(x,b,n*sizeof(double));
    double s;
    dgesc2_(n,LU[0],n,x,ipiv,jpiv,s);
    cout << "\tscale = " << s << endl;
    cout << "\tx = ";
    for (int j=0;j<n;j++) cout << x[j] << " ";
    cout << endl;

    dgemv_('N',n,n,-1.,A[0],n,x,1,1.,b,1);
    cout << "\terror = ";
    for (int i=0;i<n;i++) cout << b[i] << " ";
    cout << endl;
  }

  {
    cout << "\ncalling lapack++:" << endl;
    int n=10;
    SquareMatrix<double,double> A(n);
    for (int j=0;j<n;j++) {
      for (int i=0;i<=j;i++) A(i,j)=drand48();
    }

    GaussianFactorization<double,double,SquareMatrix<double,double> >
      gf(A,Factorization::PIVOT_ROWS_AND_COLUMNS);

    Vector<double,double> b(n);
    for (int i=0;i<n;i++) b[i]=drand48();

    Vector<double,double> x(n);
    gf.solve(b,x);
    cout << "\tx = ";
    x.printOn(cout);

    A.gemv(1.,x,-1.,b);
    cout << "\terror = ";
    b.printOn(cout);
  }
}
