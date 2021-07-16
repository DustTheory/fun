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
  int m=10;
  float L[m][m];
  for (int j=0;j<m;j++) {
    L[j][j]=1.;
    for (int i=j+1;i<m;i++) L[j][i]=drand48();
  }

  float b[m];
  for (int i=0;i<m;i++) b[i]=drand48();
  cout << "b = ";
  for (int i=0;i<m;i++) cout << b[i] << " ";
  cout << endl;

  float y[m];
  memcpy(y,b,m*sizeof(float));
  strsv_('L','N','U',m,L[0],m,y,1);
  cout << "y = ";
  for (int i=0;i<m;i++) cout << y[i] << " ";
  cout << endl;

  float Ly[m];
  memcpy(Ly,y,m*sizeof(float));
  strmv_('L','N','U',m,L[0],m,Ly,1);
  saxpy_(m,-1.,Ly,1,b,1);
  cout << "error = ";
  for (int i=0;i<m;i++) cout << b[i] << " ";
  cout << endl;
}
