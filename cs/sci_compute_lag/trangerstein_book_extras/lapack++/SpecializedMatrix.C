#include "SpecializedMatrix.H"
#include <limits>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> CompanionMatrix<R,C>::CompanionMatrix(
const Vector<R,C> &p) : UpperHessenbergMatrix<R,C>(p.size()-1,
Matrix<R,C>::zero_) {
  int n=this->size(0);
  for (int j=0;j<n;j++) (*this)(0,j)=-p[j+1]/p[0];
  for (int j=0;j<n-1;j++) (*this)(j+1,j)=Matrix<R,C>::one_;
}

template<typename R,typename C> void CompanionMatrix<R,C>::printOn(
ostream& s) const {
  s << "CompanionMatrix" << endl;
  UpperHessenbergMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> HadamardMatrix<R,C>::HadamardMatrix(int n) :
SquareMatrix<R,C>(n,Matrix<R,C>::one_) {
  int k=0;
  int sz=-1;
  if (n%20 == 0) {
    int pof2=n/10;
    while (pof2%2==0) {
      pof2/=2;
      k++;
    }
    CHECK_TEST(pof2==1);
    (*this)(1,1)=Matrix<R,C>::mone_;
    (*this)(2,1)=Matrix<R,C>::mone_;
    (*this)(5,1)=Matrix<R,C>::mone_;
    (*this)(6,1)=Matrix<R,C>::mone_;
    (*this)(7,1)=Matrix<R,C>::mone_;
    (*this)(8,1)=Matrix<R,C>::mone_;
    (*this)(10,1)=Matrix<R,C>::mone_;
    (*this)(12,1)=Matrix<R,C>::mone_;
    (*this)(17,1)=Matrix<R,C>::mone_;
    (*this)(18,1)=Matrix<R,C>::mone_;
    (*this)(1,2)=Matrix<R,C>::mone_;
    (*this)(4,2)=Matrix<R,C>::mone_;
    (*this)(5,2)=Matrix<R,C>::mone_;
    (*this)(6,2)=Matrix<R,C>::mone_;
    (*this)(7,2)=Matrix<R,C>::mone_;
    (*this)(9,2)=Matrix<R,C>::mone_;
    (*this)(11,2)=Matrix<R,C>::mone_;
    (*this)(16,2)=Matrix<R,C>::mone_;
    (*this)(17,2)=Matrix<R,C>::mone_;
    (*this)(19,2)=Matrix<R,C>::mone_;
    (*this)(3,3)=Matrix<R,C>::mone_;
    (*this)(4,3)=Matrix<R,C>::mone_;
    (*this)(5,3)=Matrix<R,C>::mone_;
    (*this)(6,3)=Matrix<R,C>::mone_;
    (*this)(8,3)=Matrix<R,C>::mone_;
    (*this)(10,3)=Matrix<R,C>::mone_;
    (*this)(15,3)=Matrix<R,C>::mone_;
    (*this)(16,3)=Matrix<R,C>::mone_;
    (*this)(18,3)=Matrix<R,C>::mone_;
    (*this)(19,3)=Matrix<R,C>::mone_;
    (*this)(2,4)=Matrix<R,C>::mone_;
    (*this)(3,4)=Matrix<R,C>::mone_;
    (*this)(4,4)=Matrix<R,C>::mone_;
    (*this)(5,4)=Matrix<R,C>::mone_;
    (*this)(7,4)=Matrix<R,C>::mone_;
    (*this)(9,4)=Matrix<R,C>::mone_;
    (*this)(14,4)=Matrix<R,C>::mone_;
    (*this)(15,4)=Matrix<R,C>::mone_;
    (*this)(17,4)=Matrix<R,C>::mone_;
    (*this)(18,4)=Matrix<R,C>::mone_;
    (*this)(1,5)=Matrix<R,C>::mone_;
    (*this)(2,5)=Matrix<R,C>::mone_;
    (*this)(3,5)=Matrix<R,C>::mone_;
    (*this)(4,5)=Matrix<R,C>::mone_;
    (*this)(6,5)=Matrix<R,C>::mone_;
    (*this)(8,5)=Matrix<R,C>::mone_;
    (*this)(13,5)=Matrix<R,C>::mone_;
    (*this)(14,5)=Matrix<R,C>::mone_;
    (*this)(16,5)=Matrix<R,C>::mone_;
    (*this)(17,5)=Matrix<R,C>::mone_;
    (*this)(1,6)=Matrix<R,C>::mone_;
    (*this)(2,6)=Matrix<R,C>::mone_;
    (*this)(3,6)=Matrix<R,C>::mone_;
    (*this)(5,6)=Matrix<R,C>::mone_;
    (*this)(7,6)=Matrix<R,C>::mone_;
    (*this)(12,6)=Matrix<R,C>::mone_;
    (*this)(13,6)=Matrix<R,C>::mone_;
    (*this)(15,6)=Matrix<R,C>::mone_;
    (*this)(16,6)=Matrix<R,C>::mone_;
    (*this)(19,6)=Matrix<R,C>::mone_;
    (*this)(1,7)=Matrix<R,C>::mone_;
    (*this)(2,7)=Matrix<R,C>::mone_;
    (*this)(4,7)=Matrix<R,C>::mone_;
    (*this)(6,7)=Matrix<R,C>::mone_;
    (*this)(11,7)=Matrix<R,C>::mone_;
    (*this)(12,7)=Matrix<R,C>::mone_;
    (*this)(14,7)=Matrix<R,C>::mone_;
    (*this)(15,7)=Matrix<R,C>::mone_;
    (*this)(18,7)=Matrix<R,C>::mone_;
    (*this)(19,7)=Matrix<R,C>::mone_;
    (*this)(1,8)=Matrix<R,C>::mone_;
    (*this)(3,8)=Matrix<R,C>::mone_;
    (*this)(5,8)=Matrix<R,C>::mone_;
    (*this)(10,8)=Matrix<R,C>::mone_;
    (*this)(11,8)=Matrix<R,C>::mone_;
    (*this)(13,8)=Matrix<R,C>::mone_;
    (*this)(14,8)=Matrix<R,C>::mone_;
    (*this)(17,8)=Matrix<R,C>::mone_;
    (*this)(18,8)=Matrix<R,C>::mone_;
    (*this)(19,8)=Matrix<R,C>::mone_;
    (*this)(2,9)=Matrix<R,C>::mone_;
    (*this)(4,9)=Matrix<R,C>::mone_;
    (*this)(9,9)=Matrix<R,C>::mone_;
    (*this)(10,9)=Matrix<R,C>::mone_;
    (*this)(12,9)=Matrix<R,C>::mone_;
    (*this)(13,9)=Matrix<R,C>::mone_;
    (*this)(16,9)=Matrix<R,C>::mone_;
    (*this)(17,9)=Matrix<R,C>::mone_;
    (*this)(18,9)=Matrix<R,C>::mone_;
    (*this)(19,9)=Matrix<R,C>::mone_;
    (*this)(1,10)=Matrix<R,C>::mone_;
    (*this)(3,10)=Matrix<R,C>::mone_;
    (*this)(8,10)=Matrix<R,C>::mone_;
    (*this)(9,10)=Matrix<R,C>::mone_;
    (*this)(11,10)=Matrix<R,C>::mone_;
    (*this)(12,10)=Matrix<R,C>::mone_;
    (*this)(15,10)=Matrix<R,C>::mone_;
    (*this)(16,10)=Matrix<R,C>::mone_;
    (*this)(17,10)=Matrix<R,C>::mone_;
    (*this)(18,10)=Matrix<R,C>::mone_;
    (*this)(2,11)=Matrix<R,C>::mone_;
    (*this)(7,11)=Matrix<R,C>::mone_;
    (*this)(8,11)=Matrix<R,C>::mone_;
    (*this)(10,11)=Matrix<R,C>::mone_;
    (*this)(11,11)=Matrix<R,C>::mone_;
    (*this)(14,11)=Matrix<R,C>::mone_;
    (*this)(15,11)=Matrix<R,C>::mone_;
    (*this)(16,11)=Matrix<R,C>::mone_;
    (*this)(17,11)=Matrix<R,C>::mone_;
    (*this)(19,11)=Matrix<R,C>::mone_;
    (*this)(1,12)=Matrix<R,C>::mone_;
    (*this)(6,12)=Matrix<R,C>::mone_;
    (*this)(7,12)=Matrix<R,C>::mone_;
    (*this)(9,12)=Matrix<R,C>::mone_;
    (*this)(10,12)=Matrix<R,C>::mone_;
    (*this)(13,12)=Matrix<R,C>::mone_;
    (*this)(14,12)=Matrix<R,C>::mone_;
    (*this)(15,12)=Matrix<R,C>::mone_;
    (*this)(16,12)=Matrix<R,C>::mone_;
    (*this)(18,12)=Matrix<R,C>::mone_;
    (*this)(5,13)=Matrix<R,C>::mone_;
    (*this)(6,13)=Matrix<R,C>::mone_;
    (*this)(8,13)=Matrix<R,C>::mone_;
    (*this)(9,13)=Matrix<R,C>::mone_;
    (*this)(12,13)=Matrix<R,C>::mone_;
    (*this)(13,13)=Matrix<R,C>::mone_;
    (*this)(14,13)=Matrix<R,C>::mone_;
    (*this)(15,13)=Matrix<R,C>::mone_;
    (*this)(17,13)=Matrix<R,C>::mone_;
    (*this)(19,13)=Matrix<R,C>::mone_;
    (*this)(4,14)=Matrix<R,C>::mone_;
    (*this)(5,14)=Matrix<R,C>::mone_;
    (*this)(7,14)=Matrix<R,C>::mone_;
    (*this)(8,14)=Matrix<R,C>::mone_;
    (*this)(11,14)=Matrix<R,C>::mone_;
    (*this)(12,14)=Matrix<R,C>::mone_;
    (*this)(13,14)=Matrix<R,C>::mone_;
    (*this)(14,14)=Matrix<R,C>::mone_;
    (*this)(16,14)=Matrix<R,C>::mone_;
    (*this)(18,14)=Matrix<R,C>::mone_;
    (*this)(3,15)=Matrix<R,C>::mone_;
    (*this)(4,15)=Matrix<R,C>::mone_;
    (*this)(6,15)=Matrix<R,C>::mone_;
    (*this)(7,15)=Matrix<R,C>::mone_;
    (*this)(10,15)=Matrix<R,C>::mone_;
    (*this)(11,15)=Matrix<R,C>::mone_;
    (*this)(12,15)=Matrix<R,C>::mone_;
    (*this)(13,15)=Matrix<R,C>::mone_;
    (*this)(15,15)=Matrix<R,C>::mone_;
    (*this)(17,15)=Matrix<R,C>::mone_;
    (*this)(2,16)=Matrix<R,C>::mone_;
    (*this)(3,16)=Matrix<R,C>::mone_;
    (*this)(5,16)=Matrix<R,C>::mone_;
    (*this)(6,16)=Matrix<R,C>::mone_;
    (*this)(9,16)=Matrix<R,C>::mone_;
    (*this)(10,16)=Matrix<R,C>::mone_;
    (*this)(11,16)=Matrix<R,C>::mone_;
    (*this)(12,16)=Matrix<R,C>::mone_;
    (*this)(14,16)=Matrix<R,C>::mone_;
    (*this)(16,16)=Matrix<R,C>::mone_;
    (*this)(1,17)=Matrix<R,C>::mone_;
    (*this)(2,17)=Matrix<R,C>::mone_;
    (*this)(4,17)=Matrix<R,C>::mone_;
    (*this)(5,17)=Matrix<R,C>::mone_;
    (*this)(8,17)=Matrix<R,C>::mone_;
    (*this)(9,17)=Matrix<R,C>::mone_;
    (*this)(10,17)=Matrix<R,C>::mone_;
    (*this)(11,17)=Matrix<R,C>::mone_;
    (*this)(13,17)=Matrix<R,C>::mone_;
    (*this)(15,17)=Matrix<R,C>::mone_;
    (*this)(1,18)=Matrix<R,C>::mone_;
    (*this)(3,18)=Matrix<R,C>::mone_;
    (*this)(4,18)=Matrix<R,C>::mone_;
    (*this)(7,18)=Matrix<R,C>::mone_;
    (*this)(8,18)=Matrix<R,C>::mone_;
    (*this)(9,18)=Matrix<R,C>::mone_;
    (*this)(10,18)=Matrix<R,C>::mone_;
    (*this)(12,18)=Matrix<R,C>::mone_;
    (*this)(14,18)=Matrix<R,C>::mone_;
    (*this)(19,18)=Matrix<R,C>::mone_;
    (*this)(2,19)=Matrix<R,C>::mone_;
    (*this)(3,19)=Matrix<R,C>::mone_;
    (*this)(6,19)=Matrix<R,C>::mone_;
    (*this)(7,19)=Matrix<R,C>::mone_;
    (*this)(8,19)=Matrix<R,C>::mone_;
    (*this)(9,19)=Matrix<R,C>::mone_;
    (*this)(11,19)=Matrix<R,C>::mone_;
    (*this)(13,19)=Matrix<R,C>::mone_;
    (*this)(18,19)=Matrix<R,C>::mone_;
    (*this)(19,19)=Matrix<R,C>::mone_;
    sz=20;
  } else if (n%12 == 0) {
    int pof2=n/6;
    while (pof2%2==0) {
      pof2/=2;
      k++;
    }
    CHECK_TEST(pof2==1);
    (*this)(1,1)=Matrix<R,C>::mone_;
    (*this)(2,1)=Matrix<R,C>::mone_;
    (*this)(4,1)=Matrix<R,C>::mone_;
    (*this)(5,1)=Matrix<R,C>::mone_;
    (*this)(6,1)=Matrix<R,C>::mone_;
    (*this)(10,1)=Matrix<R,C>::mone_;
    (*this)(2,2)=Matrix<R,C>::mone_;
    (*this)(3,2)=Matrix<R,C>::mone_;
    (*this)(5,2)=Matrix<R,C>::mone_;
    (*this)(6,2)=Matrix<R,C>::mone_;
    (*this)(7,2)=Matrix<R,C>::mone_;
    (*this)(11,2)=Matrix<R,C>::mone_;
    (*this)(1,3)=Matrix<R,C>::mone_;
    (*this)(3,3)=Matrix<R,C>::mone_;
    (*this)(4,3)=Matrix<R,C>::mone_;
    (*this)(6,3)=Matrix<R,C>::mone_;
    (*this)(7,3)=Matrix<R,C>::mone_;
    (*this)(8,3)=Matrix<R,C>::mone_;
    (*this)(2,4)=Matrix<R,C>::mone_;
    (*this)(4,4)=Matrix<R,C>::mone_;
    (*this)(5,4)=Matrix<R,C>::mone_;
    (*this)(7,4)=Matrix<R,C>::mone_;
    (*this)(8,4)=Matrix<R,C>::mone_;
    (*this)(9,4)=Matrix<R,C>::mone_;
    (*this)(3,5)=Matrix<R,C>::mone_;
    (*this)(5,5)=Matrix<R,C>::mone_;
    (*this)(6,5)=Matrix<R,C>::mone_;
    (*this)(8,5)=Matrix<R,C>::mone_;
    (*this)(9,5)=Matrix<R,C>::mone_;
    (*this)(10,5)=Matrix<R,C>::mone_;
    (*this)(4,6)=Matrix<R,C>::mone_;
    (*this)(6,6)=Matrix<R,C>::mone_;
    (*this)(7,6)=Matrix<R,C>::mone_;
    (*this)(9,6)=Matrix<R,C>::mone_;
    (*this)(10,6)=Matrix<R,C>::mone_;
    (*this)(11,6)=Matrix<R,C>::mone_;
    (*this)(1,7)=Matrix<R,C>::mone_;
    (*this)(5,7)=Matrix<R,C>::mone_;
    (*this)(7,7)=Matrix<R,C>::mone_;
    (*this)(8,7)=Matrix<R,C>::mone_;
    (*this)(10,7)=Matrix<R,C>::mone_;
    (*this)(11,7)=Matrix<R,C>::mone_;
    (*this)(1,8)=Matrix<R,C>::mone_;
    (*this)(2,8)=Matrix<R,C>::mone_;
    (*this)(6,8)=Matrix<R,C>::mone_;
    (*this)(8,8)=Matrix<R,C>::mone_;
    (*this)(9,8)=Matrix<R,C>::mone_;
    (*this)(11,8)=Matrix<R,C>::mone_;
    (*this)(1,9)=Matrix<R,C>::mone_;
    (*this)(2,9)=Matrix<R,C>::mone_;
    (*this)(3,9)=Matrix<R,C>::mone_;
    (*this)(7,9)=Matrix<R,C>::mone_;
    (*this)(9,9)=Matrix<R,C>::mone_;
    (*this)(10,9)=Matrix<R,C>::mone_;
    (*this)(2,10)=Matrix<R,C>::mone_;
    (*this)(3,10)=Matrix<R,C>::mone_;
    (*this)(4,10)=Matrix<R,C>::mone_;
    (*this)(8,10)=Matrix<R,C>::mone_;
    (*this)(10,10)=Matrix<R,C>::mone_;
    (*this)(11,10)=Matrix<R,C>::mone_;
    (*this)(1,11)=Matrix<R,C>::mone_;
    (*this)(3,11)=Matrix<R,C>::mone_;
    (*this)(4,11)=Matrix<R,C>::mone_;
    (*this)(5,11)=Matrix<R,C>::mone_;
    (*this)(9,11)=Matrix<R,C>::mone_;
    (*this)(11,11)=Matrix<R,C>::mone_;
    sz=12;
  } else {
    int pof2=n;
    while (pof2%2==0) {
      pof2/=2;
      k++;
    }
    CHECK_TEST(pof2==1);
    (*this)(1,1)=Matrix<R,C>::mone_;
    sz=2;
  }
  for (int kk=2;kk<=k;kk++,sz*=2) {
    for (int j=0;j<sz;j++) {
      memcpy(this->addr(sz,j),this->addr(0,j),sz*sizeof(C));
      memcpy(this->addr(0,j+sz),this->addr(0,j),sz*sizeof(C));
      for (int i=0;i<sz;i++) {
//      (*this)(i+sz,j)=(*this)(i,j);
//      (*this)(i,j+sz)=(*this)(i,j);
        (*this)(i+sz,j+sz)=-(*this)(i,j);
      }
    }
  }
}

template<typename R,typename C> void HadamardMatrix<R,C>::printOn(
ostream& s) const {
  s << "HadamardMatrix" << endl;
  SquareMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> HankelMatrix<R,C>::HankelMatrix(
const Vector<R,C> &d) : SymmetricMatrix<R,C>((d.size()+1)/2) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    memcpy(this->addr(j,j),d.addr(j+j),(n-j)*sizeof(C));
//  for (int i=j;i<n;i++) (*this)(i,j)=d[i+j];
  }
}

template<typename R,typename C> HankelMatrix<R,C>::HankelMatrix(
const Vector<R,C> &c,const Vector<R,C> &r) : 
SymmetricMatrix<R,C>(c.size()) {
  int n=this->size(0);
  CHECK_SAME(n,r.size());
  for (int j=0;j<n;j++) {
    if (2*j<n) {
      memcpy(this->addr(j,j),c.addr(2*j),(n-2*j)*sizeof(C));
      if (j>0) memcpy(this->addr(n-j,j),r.addr(1),j*sizeof(C));
    } else memcpy(this->addr(j,j),r.addr(2*j-n+1),(n-j)*sizeof(C));
//  for (int i=j;i<n;i++) {
//    int k=i+j;
//    (*this)(i,j)=(k<n ? c[k] : r[k-n+1]);
//  }
  }
}

template<typename R,typename C> void HankelMatrix<R,C>::printOn(
ostream& s) const {
  s << "HankelMatrix" << endl;
  SquareMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> HilbertMatrix<R,C>::HilbertMatrix(
int m,int n) : Matrix<R,C>(m,n) {
  Vector<R,C> temp(m+n-1);
  for (int k=0;k<m+n-1;k++) temp[k]=Matrix<R,C>::one_/static_cast<R>(k+1);
  for (int j=0;j<n;j++) {
    memcpy(this->addr(0,j),temp.addr(j),m*sizeof(C));
//  for (int i=0;i<n;i++) (*this)(i,j)=temp[i+j];
  }
}

template<typename R,typename C> void HilbertMatrix<R,C>::printOn(
ostream& s) const {
  s << "HilbertMatrix" << endl;
  Matrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> KahanMatrix<R,C>::KahanMatrix(int n,C s) :
UpperTriangularMatrix<R,C>(n) {
  R rone=static_cast<R>(1.);
  s=Matrix<R,C>::one_*min(rone,abs(s));
  R a=abs(s);
  C c=Matrix<R,C>::one_*sqrt(rone-a*a);
  C diag=Matrix<R,C>::one_;
  C offdiag=-c;
  for (int i=0;i<n;i++,diag*=s,offdiag*=s) {
    (*this)(i,i)=diag;
    for (int j=i+1;j<n;j++) (*this)(i,j)=offdiag;
  }
}

template<typename R,typename C> void KahanMatrix<R,C>::printOn(
ostream& s) const {
  s << "KahanMatrix" << endl;
  UpperTriangularMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> LauchliMatrix<R,C>::LauchliMatrix(int n) :
Matrix<R,C>(n+1,n,Matrix<R,C>::zero_) {
  R eps=sqrt(numeric_limits<R>::epsilon());
  C small=Matrix<R,C>::one_*eps;
  for (int j=0;j<n;j++) {
    (*this)(0,j)=Matrix<R,C>::one_;
    (*this)(j+1,j)=small;
  }
}

template<typename R,typename C> void LauchliMatrix<R,C>::printOn(
ostream& s) const {
  s << "LauchliMatrix" << endl;
  Matrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> PascalMatrix<R,C>::PascalMatrix(int n) :
SymmetricMatrix<R,C>(n,Matrix<R,C>::one_) {
  for (int k=1;k<2*n-1;k++) {
    for (int i=(k-1)/2+1;i<=min(k-1,n-1);i++) {
      int j=k-i;
      if (i-1>=j) (*this)(i,j)=(*this)(i,j-1)+(*this)(i-1,j);
      else (*this)(i,j)=(*this)(i,j-1)+(*this)(j,i-1);
    }
  }
}

template<typename R,typename C> void PascalMatrix<R,C>::printOn(
ostream& s) const {
  s << "PascalMatrix" << endl;
  SymmetricMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> RosserMatrix<R,C>::RosserMatrix() :
SymmetricMatrix<R,C>(8,Matrix<R,C>::one_) {
  (*this)(0,0)*=611.;
  (*this)(1,0)*=196.;
  (*this)(2,0)*=-192.;
  (*this)(3,0)*=407.;
  (*this)(4,0)*=-8.;
  (*this)(5,0)*=-52.;
  (*this)(6,0)*=-49.;
  (*this)(7,0)*=29.;
  (*this)(1,1)*=899.;
  (*this)(2,1)*=113.;
  (*this)(3,1)*=-192.;
  (*this)(4,1)*=-71.;
  (*this)(5,1)*=-43;
  (*this)(6,1)*=-8.;
  (*this)(7,1)*=-44.;
  (*this)(2,2)*=899.;
  (*this)(3,2)*=196.;
  (*this)(4,2)*=61.;
  (*this)(5,2)*=49;
  (*this)(6,2)*=8.;
  (*this)(7,2)*=52.;
  (*this)(3,3)*=611.;
  (*this)(4,3)*=8.;
  (*this)(5,3)*=44;
  (*this)(6,3)*=59.;
  (*this)(7,3)*=-23.;
  (*this)(4,4)*=411.;
  (*this)(5,4)*=-599;
  (*this)(6,4)*=208.;
  (*this)(7,4)*=208.;
  (*this)(5,5)*=411;
  (*this)(6,5)*=208.;
  (*this)(7,5)*=208.;
  (*this)(6,6)*=99.;
  (*this)(7,6)*=-911.;
  (*this)(7,7)*=99.;
}

template<typename R,typename C> void RosserMatrix<R,C>::printOn(
ostream& s) const {
  s << "RosserMatrix" << endl;
  SymmetricMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> ToeplitzMatrix<R,C>::ToeplitzMatrix(
const Vector<R,C> &c,const Vector<R,C> &r) : SquareMatrix<R,C>(c.size()) {
  int n=this->size(0);
  CHECK_SAME(n,r.size());
  Vector<R,C> temp(n);
  for (int j=0;j<n;j++) temp[n-1-j]=r[j];
  for (int j=0;j<n;j++) {
    if (j>0) memcpy(this->addr(0,j),temp.addr(n-j),j*sizeof(C));
    memcpy(this->addr(j,j),c.addr(0),(n-j)*sizeof(C));
  }
}

template<typename R,typename C> void ToeplitzMatrix<R,C>::printOn(
ostream& s) const {
  s << "ToeplitzMatrix" << endl;
  SquareMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C>
SymmetricToeplitzMatrix<R,C>::SymmetricToeplitzMatrix(
const Vector<R,C> &c) : SymmetricMatrix<R,C>(c.size()) {
  int n=this->size(0);
  for (int j=0;j<n;j++) {
    memcpy(this->addr(j,j),c.addr(0),(n-j)*sizeof(C));
  }
}

template<typename R,typename C> void SymmetricToeplitzMatrix<R,C>::printOn(
ostream& s) const {
  s << "SymmetricToeplitzMatrix" << endl;
  SymmetricMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> VandermondeMatrix<R,C>::VandermondeMatrix(
const Vector<R,C> &x) : SquareMatrix<R,C>(x.size()) {
  int n=this->size(0);
  for (int i=0;i<n;i++) (*this)(i,0)=Matrix<R,C>::one_;
  for (int j=1;j<n;j++) {
    for (int i=0;i<n;i++) {
      (*this)(i,j)=(*this)(i,j-1)*x[i];
    }
  }
}

template<typename R,typename C> void VandermondeMatrix<R,C>::printOn(
ostream& s) const {
  s << "VandermondeMatrix" << endl;
  SquareMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C>
WilkinsonMatrix<R,C>::WilkinsonMatrix(int n) :
SymmetricTridiagonalMatrix<R,C>(n,Matrix<R,C>::one_) {
  for (int j=0;j<(n+1)/2;j++) {
    this->diagonalValue(j)=static_cast<R>(0.5)*static_cast<R>(n-1-2*j);
    this->diagonalValue(n-1-j)=this->diagonalValue(j);
  }
}

template<typename R,typename C> void WilkinsonMatrix<R,C>::printOn(
ostream& s) const {
  s << "WilkinsonMatrix" << endl;
  SymmetricTridiagonalMatrix<R,C>::printOn(s);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> void testSpecializedMatrix(R fscalar,
C scalar) {
/*
  {
    Vector<R,C> p(4);
    for (int i=0;i<4;i++) p[i]=scalar*(4-i);
    CompanionMatrix<R,C> M(p);
    M.printOn(cout);
  }
*/
/*
  {
    HadamardMatrix<R,C> H(2);
    H.printOn(cout);
  }
*/
/*
  {
    HadamardMatrix<R,C> H(4);
    H.printOn(cout);
  }
*/
/*
  {
    HadamardMatrix<R,C> H(8);
    SquareMatrix<R,C> P(8,scalar-scalar);
    P.gemm(scalar/scalar,H,H,scalar-scalar,'C','N');
    H.printOn(cout);
    P.printOn(cout);
  }
*/
/*
  {
    HadamardMatrix<R,C> H(12);
    SquareMatrix<R,C> P(12,scalar-scalar);
    P.gemm(scalar/scalar,H,H,scalar-scalar,'C','N');
    H.printOn(cout);
    P.printOn(cout);
  }
*/
/*
  {
    HadamardMatrix<R,C> H(24);
    SquareMatrix<R,C> P(24,scalar-scalar);
    P.gemm(scalar/scalar,H,H,scalar-scalar,'C','N');
    H.printOn(cout);
    P.printOn(cout);
  }
*/
/*
  {
    HadamardMatrix<R,C> H(20);
    SquareMatrix<R,C> P(20,scalar-scalar);
    P.gemm(scalar/scalar,H,H,scalar-scalar,'C','N');
    H.printOn(cout);
    P.printOn(cout);
  }
*/
/*
  {
    HadamardMatrix<R,C> H(40);
    SquareMatrix<R,C> P(40,scalar-scalar);
    P.gemm(scalar/scalar,H,H,scalar-scalar,'C','N');
    H.printOn(cout);
    P.printOn(cout);
  }
*/
/*
  {
    Vector<R,C> d(5);
    for (int i=0;i<5;i++) d[i]=scalar*(i+1);
    HankelMatrix<R,C> H(d);
    H.printOn(cout);
  }
*/
/*
  {
    Vector<R,C> d(6);
    for (int i=0;i<6;i++) d[i]=scalar*(i+1);
    HankelMatrix<R,C> H(d);
    H.printOn(cout);
  }
*/
/*
  {
    Vector<R,C> c(3);
    Vector<R,C> r(3);
    for (int i=0;i<3;i++) {
      c[i]=scalar*(i+1);
      r[i]=scalar*(1-i);
    }
    HankelMatrix<R,C> H(c,r);
    H.printOn(cout);
  }
*/
/*
  {
    HilbertMatrix<R,C> H(5,4);
    H.printOn(cout);
  }
*/
/*
  {
    KahanMatrix<R,C> K(5,scalar/pow(abs(scalar),2));
    K.printOn(cout);
  }
*/
/*
  {
    LauchliMatrix<R,C> L(3);
    L.printOn(cout);
  }
*/
/*
  {
    PascalMatrix<R,C> P(5);
    P.printOn(cout);
  }
*/
/*
  {
    RosserMatrix<R,C> M;
    M.printOn(cout);
  }
*/
/*
  {
    Vector<R,C> c(3);
    Vector<R,C> r(3);
    for (int i=0;i<3;i++) {
      c[i]=scalar*(i+1);
      r[i]=scalar*(1-i);
    }
    ToeplitzMatrix<R,C> T(c,r);
    T.printOn(cout);
  }
*/
/*
  {
    Vector<R,C> c(3);
    for (int i=0;i<3;i++) c[i]=scalar*(i+1);
    SymmetricToeplitzMatrix<R,C> S(c);
    S.printOn(cout);
  }
*/
/*
  {
    Vector<R,C> x(4);
    for (int i=0;i<4;i++) x[i]=scalar*(i+1);
    VandermondeMatrix<R,C> V(x);
    V.printOn(cout);
  }
*/
  {
    WilkinsonMatrix<R,C> W(4);
    W.printOn(cout);
  }
  {
    WilkinsonMatrix<R,C> W(5);
    W.printOn(cout);
  }
}
