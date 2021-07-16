//triangular systems overwrite A:
//  left-triangular coefficients stored in A below the diagonal
//  right-triangular coefficients stored on and above the diagonal
void gaussElimNoPiv(int n, float **A) {
  for (int k=0;k<n;k++) {
    double diag=A[k][k];
    if (abs(diag)<=0.) return;
    int kp1=k+1;
    diag=1./diag;
    for (int i=kp1;i<n;i++) {
      A[k][i] *= diag;
      for (int j=kp1;j<n;j++) {
        A[j][i] -= A[k][i] * A[j][k];
      }
    }
  }
}
extern "C" {
  void sscal_(const int &n,const float &sa,float *sx,const int &incx);
  void sger_(const int &m,const int &n,const float &alpha,const float *x,
    const int &incx,const float *y,const int &incy,float *A,
    const int &lda);
}
void gaussElimNoPivWithBlas(int n, float **A,int lda) {
  for (int k=0;k<n;k++) {
    double diag=A[k][k];
    if (abs(diag)<=0.) return;
    int kp1=k+1;
    sscal_(n-k-1,1./diag,&A[k][kp1],1);
    sger_(n-k-1,n-k-1,-1.,&A[k][kp1],1,&A[kp1][k],lda,&A[kp1][kp1],lda);
  }
}
