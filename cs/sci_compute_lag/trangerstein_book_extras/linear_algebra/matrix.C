#include "matrix.H"
#include "vector.H"
extern "C" {
  void strsv_simple__(const bool&,const bool&,const int&,const float*,
    const int&,float*);
  void gauss_elim__(int&,float*,int*,int*,int&);
}
void Matrix::fill(float alpha) {
  for (int i=0;i<sz;i++) data[i]=alpha;
}
void Matrix::copy(const Matrix &G) {
  setSizes(G.sizes[0],G.sizes[1]);
  for (int j=0;j<sizes[1];j++) {
    for (int i=0;i<sizes[0];i++) (*this)(i,j)=G(i,j);
  }
}
void Matrix::copyInto(Matrix &G) const {
  int Gm=G.size(0),Am=size(0);
  int m=(Am<Gm ? Am : Gm);
  int n=(size(1)<G.size(1) ? size(1) : G.size(1));
  const float *src=data; float *dest=G.data; int len=m*sizeof(float);
  for (int j=0;j<n;j++,src+=Am,dest+=Gm) { memcpy(dest,src,len); }
}
Matrix& Matrix::operator+=(const Matrix &G) {
  ASSERT(sameSizeAs(G));
  for (int j=0;j<sizes[1];j++) {
    for (int i=0;i<sizes[0];i++) (*this)(i,j)+=G(i,j);
  }
  return *this;
}
Matrix& Matrix::operator-=(const Matrix &G) {
  ASSERT(sameSizeAs(G));
  for (int j=0;j<sizes[1];j++) {
    for (int i=0;i<sizes[0];i++) (*this)(i,j)-=G(i,j);
  }
  return *this;
}
Matrix& Matrix::operator*=(float alpha) {
  for (int ij=0;ij<sz;ij++) data[ij]*=alpha; return *this;
}
Matrix& Matrix::operator/=(float alpha) {
  ASSERT(alpha!=0.); float denom=1./alpha;
  for (int ij=0;ij<sz;ij++) data[ij]*=denom; return *this;
}
Matrix* Matrix::operator+(const Matrix &G) const {
  ASSERT(sameSizeAs(G));
  Matrix *sum=new Matrix(*this); (*sum) += G; return sum;
}
Matrix* Matrix::operator-(const Matrix &G) const {
  ASSERT(sameSizeAs(G));
  Matrix *dif=new Matrix(*this); (*dif) -= G; return dif;
}
Vector* Matrix::operator*(const Vector &x) const {
  int m=sizes[0],n=sizes[1]; ASSERT(n==x.size());
  Vector *Ax=new Vector(m); *Ax=0.;
  for (int j=0;j<n;j++) {
    float xj=x[j]; for (int i=0;i<m;i++) (*Ax)[i] += (*this)(i,j)*xj;
  }
}
Matrix* Matrix::operator*(const Matrix &X) const {
  int m=sizes[0],n=sizes[1],p=X.sizes[1];
  ASSERT(n==X.sizes[0]);
  Matrix *AX=new Matrix(m,p); AX->fill(0.);
  for (int k=0;k<p;k++) {
    for (int j=0;j<n;j++) {
      float Xjk=X(j,k);
      for (int i=0;i<m;i++) (*AX)(i,k) += (*this)(i,j)*Xjk;
    }
  }
  return AX;
}
Matrix* Matrix::transpose() const {
  int m=size(0),n=size(1);
  Matrix *X=new Matrix(n,m);
  for (int j=0;j<n;j++) {
    for (int i=0;i<m;i++) (*X)(j,i)=operator()(i,j);
  }
  return X;
}
void Matrix::printOn(ostream& s) const {
  s << "Matrix(" << size(0) << " x " << size(1) << ")\n" ;
  for (int i=0; i< size(0); i++) {
    for (int j=0; j< size(1); j++) s << operator()(i,j) << "  ";
    s << "\n";
  }
}
Vector* Matrix::solve(const Vector &b) const {
  int m=sizes[0], n=sizes[1]; ASSERT(m==b.size());
  Matrix Af(m,n); Af.copy(*this);
  Vector *x=new Vector(n); Vector bf(b.size()); bf.copy(b);

  int *ipiv=new int[m]; int *jpiv=new int[n]; int npivots;
  gauss_elim__(m,Af.data,ipiv,jpiv,npivots);
  { for (int i=0;i<npivots;i++) {
    int ip=ipiv[i]-1;
    if (ip!=i) { float t=bf[i]; bf[i]=bf[ip]; bf[ip]=t; }
  } }

  bool upper=false,unit=true;
  strsv_simple__(upper,unit,npivots,Af.data,m,&bf[0]);
  { for (int i=0;i<npivots;i++) (*x)[i]=bf[i]; }
  { for (int i=npivots;i<n;i++) (*x)[i]=0.; }

  upper=1; unit=0;
  strsv_simple__(upper,unit,npivots,Af.data,m,&(*x)[0]);
  { for (int j=npivots-1;j>=0;j--) {
    int jp=jpiv[j]-1;
    if (jp!=j) { float t=(*x)[j]; (*x)[j]=(*x)[jp]; (*x)[jp]=t; }
  } }
  delete [] ipiv; delete [] jpiv;
  return x;
}
