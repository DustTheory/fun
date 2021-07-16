#ifndef _SPECIALIZED_MATRIX_H
#define _SPECIALIZED_MATRIX_H
#include "BandMatrix.H"

template<typename R,typename C> class CompanionMatrix :
public UpperHessenbergMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit CompanionMatrix(const CompanionMatrix &M) :
    UpperHessenbergMatrix<R,C>(M.size(0)) { 
      UpperHessenbergMatrix<R,C>::copy(M);
    }
//  use copy instead of operator=
    CompanionMatrix& operator=(const CompanionMatrix &M) {
      resize(M.size(0));
      UpperHessenbergMatrix<R,C>::copy(M);
      return *this;
    }
    virtual CompanionMatrix& operator=(C scalar) {
      UpperHessenbergMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
//  CompanionMatrix() : UpperHessenbergMatrix<R,C>() { }
//  explicit CompanionMatrix(int n) : UpperHessenbergMatrix<R,C>(n) {}
    explicit CompanionMatrix(const Vector<R,C> &p);
    virtual ~CompanionMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class HadamardMatrix :
public SquareMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit HadamardMatrix(const HadamardMatrix &H) :
    SquareMatrix<R,C>(H.size(0)) { 
      SquareMatrix<R,C>::copy(H);
    }
//  use copy instead of operator=
    HadamardMatrix& operator=(const HadamardMatrix &H) {
      resize(H.size(0));
      SquareMatrix<R,C>::copy(H);
      return *this;
    }
    virtual HadamardMatrix& operator=(C scalar) {
      SquareMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
//  HadamardMatrix() : SquareMatrix<R,C>() { }
    explicit HadamardMatrix(int n);
    virtual ~HadamardMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class HankelMatrix :
public SymmetricMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit HankelMatrix(const HankelMatrix &H) :
    SymmetricMatrix<R,C>(H.size(0)) { 
      SymmetricMatrix<R,C>::copy(H);
    }
//  use copy instead of operator=
    HankelMatrix& operator=(const HankelMatrix &H) {
      resize(H.size(0));
      SymmetricMatrix<R,C>::copy(H);
      return *this;
    }
    virtual HankelMatrix& operator=(C scalar) {
      SymmetricMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
//  HankelMatrix() : SymmetricMatrix<R,C>() { }
//  explicit HankelMatrix(int n) : SymmetricMatrix<R,C>(n) {}
    explicit HankelMatrix(const Vector<R,C> &d);
    explicit HankelMatrix(const Vector<R,C> &c,const Vector<R,C> &r);
    virtual ~HankelMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class HilbertMatrix :
public Matrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit HilbertMatrix(const HilbertMatrix &H) :
    Matrix<R,C>(H.size(0),H.size(1)) { 
      Matrix<R,C>::copy(H);
    }
//  use copy instead of operator=
    HilbertMatrix& operator=(const HilbertMatrix &H) {
      resize(H.size(0),H.size(1));
      Matrix<R,C>::copy(H);
      return *this;
    }
    virtual HilbertMatrix& operator=(C scalar) {
      Matrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
//  HilbertMatrix() : Matrix<R,C>() { }
    explicit HilbertMatrix(int m,int n);
    virtual ~HilbertMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class KahanMatrix :
public UpperTriangularMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit KahanMatrix(const KahanMatrix &K) :
    UpperTriangularMatrix<R,C>(K.size(0)) { 
      UpperTriangularMatrix<R,C>::copy(K);
    }
//  use copy instead of operator=
    KahanMatrix& operator=(const KahanMatrix &K) {
      resize(K.size(0),K.size(1));
      UpperTriangularMatrix<R,C>::copy(K);
      return *this;
    }
    virtual KahanMatrix& operator=(C scalar) {
      Matrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
//  KahanMatrix() : UpperTriangularMatrix<R,C>() { }
    explicit KahanMatrix(int n,C s);
    virtual ~KahanMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class LauchliMatrix : public Matrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit LauchliMatrix(const LauchliMatrix &L) :
    Matrix<R,C>(L.size(0),L.size(1)) { 
      Matrix<R,C>::copy(L);
    }
//  use copy instead of operator=
    LauchliMatrix& operator=(const LauchliMatrix &L) {
      resize(L.size(0),L.size(1));
      Matrix<R,C>::copy(L);
      return *this;
    }
    virtual LauchliMatrix& operator=(C scalar) {
      Matrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
//  LauchliMatrix() : Matrix<R,C>() { }
    explicit LauchliMatrix(int n);
    virtual ~LauchliMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class PascalMatrix :
public SymmetricMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit PascalMatrix(const PascalMatrix &P) :
    SymmetricMatrix<R,C>(P.size(0)) { 
      SymmetricMatrix<R,C>::copy(P);
    }
//  use copy instead of operator=
    PascalMatrix& operator=(const PascalMatrix &P) {
      resize(P.size(0));
      SymmetricMatrix<R,C>::copy(P);
      return *this;
    }
    virtual PascalMatrix& operator=(C scalar) {
      SymmetricMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
//  PascalMatrix() : Matrix<R,C>() { }
    explicit PascalMatrix(int n);
    virtual ~PascalMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class RosserMatrix :
public SymmetricMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit RosserMatrix(const RosserMatrix &M) :
    SymmetricMatrix<R,C>(M.size(0)) { 
      SymmetricMatrix<R,C>::copy(M);
    }
//  use copy instead of operator=
    RosserMatrix& operator=(const RosserMatrix &M) {
      resize(M.size(0));
      SymmetricMatrix<R,C>::copy(M);
      return *this;
    }
    virtual RosserMatrix& operator=(C scalar) {
      SymmetricMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
//  RosserMatrix() : Matrix<R,C>() { }
    explicit RosserMatrix();
    virtual ~RosserMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class ToeplitzMatrix :
public SquareMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit ToeplitzMatrix(const ToeplitzMatrix &T) :
    SquareMatrix<R,C>(T.size(0)) { 
      SquareMatrix<R,C>::copy(T);
    }
//  use copy instead of operator=
    ToeplitzMatrix& operator=(const ToeplitzMatrix &T) {
      resize(T.size(0));
      SquareMatrix<R,C>::copy(T);
      return *this;
    }
    virtual ToeplitzMatrix& operator=(C scalar) {
      SquareMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
    explicit ToeplitzMatrix(const Vector<R,C> &c,const Vector<R,C> &r);
    virtual ~ToeplitzMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class SymmetricToeplitzMatrix :
public SymmetricMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit SymmetricToeplitzMatrix(const SymmetricToeplitzMatrix &S) :
    SymmetricMatrix<R,C>(S.size(0)) { 
      SymmetricMatrix<R,C>::copy(S);
    }
//  use copy instead of operator=
    SymmetricToeplitzMatrix& operator=(const SymmetricToeplitzMatrix &S) {
      resize(S.size(0));
      SymmetricMatrix<R,C>::copy(S);
      return *this;
    }
    virtual SymmetricToeplitzMatrix& operator=(C scalar) {
      SymmetricMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
    explicit SymmetricToeplitzMatrix(const Vector<R,C> &c);
    virtual ~SymmetricToeplitzMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class VandermondeMatrix :
public SquareMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit VandermondeMatrix(const VandermondeMatrix &V) :
    SquareMatrix<R,C>(V.size(0)) { 
      SquareMatrix<R,C>::copy(V);
    }
//  use copy instead of operator=
    VandermondeMatrix& operator=(const VandermondeMatrix &V) {
      resize(V.size(0));
      SquareMatrix<R,C>::copy(V);
      return *this;
    }
    virtual VandermondeMatrix& operator=(C scalar) {
      SquareMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
    explicit VandermondeMatrix(const Vector<R,C> &x);
    virtual ~VandermondeMatrix() {}
    virtual void printOn(ostream&) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class WilkinsonMatrix :
public SymmetricTridiagonalMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit WilkinsonMatrix(const WilkinsonMatrix &V) :
    SymmetricTridiagonalMatrix<R,C>(V.size(0)) { 
      SymmetricTridiagonalMatrix<R,C>::copy(V);
    }
//  use copy instead of operator=
    WilkinsonMatrix& operator=(const WilkinsonMatrix &V) {
      resize(V.size(0));
      SymmetricTridiagonalMatrix<R,C>::copy(V);
      return *this;
    }
    virtual WilkinsonMatrix& operator=(C scalar) {
      SymmetricTridiagonalMatrix<R,C>::operator=(scalar);
      return *this;
    }
  public:
    explicit WilkinsonMatrix(int n);
    virtual ~WilkinsonMatrix() {}
    virtual void printOn(ostream&) const ;
};
template<typename R,typename C> void testSpecializedMatrix(R,C);

//replaces LaPack gmc.H gmd.H gmf.H gmi.H gmli.H
// Modified from gmd.H by John Trangenstein, 11/9/96 and 5/14/12
#endif
