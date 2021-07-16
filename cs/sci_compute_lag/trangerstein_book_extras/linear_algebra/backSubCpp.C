#include <cstdlib> // drand48
#include <cstring> // memcpy
#include <iostream> // for cout
using namespace std;
extern "C" {
  void saxpy_(const int &n,const float &sa,const float *sx,const int &incx,
    float *sy,const int &incy);
  void strmv_(const char &uplo,const char &trans,const char &diag,
    const int &n,const float *A,const int &lda,const float *x,
    const int &incx);
  void strsv_(const char &uplo,const char &trans,const char &diag,
    const int &n,const float *A,const int &lda,const float *x,
    const int &incx);
}
int main(int argc,char **argv) {
  int n=10;
  float R[n][n];
  for (int j=0;j<n;j++) {
    for (int i=0;i<=j;i++) R[j][i]=drand48();
  }

  float y[n];
  for (int i=0;i<n;i++) y[i]=drand48();
  cout << "y = ";
  for (int i=0;i<n;i++) cout << y[i] << " ";
  cout << endl;

  float x[n];
  memcpy(x,y,n*sizeof(float));
  strsv_('U','N','N',n,R[0],n,x,1);
  cout << "x = ";
  for (int i=0;i<n;i++) cout << x[i] << " ";
  cout << endl;

  float Rx[n];
  memcpy(Rx,x,n*sizeof(float));
  strmv_('U','N','N',n,R[0],n,Rx,1);
  saxpy_(n,-1.,Rx,1,y,1);
  cout << "error = ";
  for (int i=0;i<n;i++) cout << y[i] << " ";
  cout << endl;
}
