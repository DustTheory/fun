#include <cstdlib> // drand48
#include <cstring> // memcpy
#include <iostream> // for cout

#include "GaussianFactorization.H"
#include "SquareMatrix.H"
using namespace std;
extern "C" {
  void dgesvx_(const char &fact,const char &trans,const int &n,
    const int &nrhs,const double *A,const int &lda,double *AF,
    const int &ldaf,int *ipiv,char &equed,double *r,double *c,
    const double *b,const int &ldb,double *x,const int &ldx,
    double &rcond,double *ferr,double *berr,double *work,int *iwork,
    int &info);
}
int main(int argc,char **argv) {
  int n=1024;
  double **A=new double*[n];
  double *space=new double[n*n];
  for (int j=0;j<n;j++) {
    A[j]=space;
    space+=n;
  }
  {
    cout << "\ncalling LAPACK directly:" << endl;
    double **LU=new double*[n];
    space=new double[n*n];
    for (int j=0;j<n;j++) {
      LU[j]=space;
      space+=n;
    }
    for (int j=0;j<n;j++) {
      for (int i=0;i<n;i++) A[j][i]=drand48();
    }
    double *b=new double[n];
    double *x=new double[n];
    double *r=new double[n];
    double *c=new double[n];
    double *work=new double[4*n];
    int *ipiv=new int[n];
    int *iwork=new int[n];

    for (int i=0;i<n;i++) b[i]=drand48();
    char equed;
    double rcond,ferr[1],berr[1];
    int info;
    dgesvx_('E','N',n,1,A[0],n,LU[0],n,ipiv,equed,r,c,b,n,x,n,
      rcond,ferr,berr,work,iwork,info);
    cout << "\tinfo,equed,rcond,ferr,berr = " << info << " " << equed
         << " " << rcond << " " << ferr[0] << " " << berr[0] << endl;

    delete [] iwork; iwork=0;
    delete [] ipiv; ipiv=0;
    delete [] work; work=0;
    delete [] c; c=0;
    delete [] r; r=0;
    delete [] x; x=0;
    delete [] b; b=0;
    delete [] LU[0]; LU[0]=0;
    delete [] LU; LU=0;
  }

  {
    cout << "\ncalling lapack++:" << endl;
    int n=1024;
    SquareMatrix<double,double> S(n);
    memcpy(S.addr(),A[0],n*n*sizeof(double));
    Vector<double,double> b(n);
    for (int i=0;i<n;i++) b[i]=drand48();

    GaussianFactorization<double,double,SquareMatrix<double,double> >
      gf(S,Factorization::PIVOT_ROWS,
        Factorization::EQUILIBRATE_ROWS_AND_COLUMNS);

    Vector<double,double> x(n);
    gf.solve(b,x);
    double rcond=gf.reciprocalConditionNumber(Factorization::INFINITY_NORM);
    double berr,ferr;
    gf.improve(b,x,berr,ferr);
    cout << "rcond,ferr,berr = " << rcond << " " << ferr << " " << berr
         << endl;
  }
  delete [] A[0]; A[0]=0;
  delete [] A; A=0;
}
