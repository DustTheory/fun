#ifndef _TRAPEZOIDAL_MATRIX_H_
#define _TRAPEZOIDAL_MATRIX_H_
#include "SquareMatrix.H"

//this pure virtual class prevents various Matrix functions and
//  introduces outofbounds_ and safety_
template<typename R,typename C> class TrapezoidalMatrix : public Matrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit TrapezoidalMatrix(const TrapezoidalMatrix &T) :
      Matrix<R,C>(T.size(0),T.size(1)) { assert(0);}
// Use copy instead
    TrapezoidalMatrix& operator=(const TrapezoidalMatrix &T) {// not virtual
      Matrix<R,C>::resize(T.size(0),T.size(1)); copy(T); return *this;
    }
//  prevent Matrix functions
    void interchangeColumns(int,int) {
      OBSOLETE("inappropriate for this class");
    }
    void interchangeRows(int,int) {
      OBSOLETE("inappropriate for this class");
    }
    void gemv(char trans,C alpha,const Vector<R,C> &x,C beta,
    Vector<R,C> &y) const {
      OBSOLETE("inappropriate for this class");
    }
    void ger(C alpha,const Vector<R,C> &x,const Vector<R,C> &y) {
      OBSOLETE("inappropriate for this class");
    }
    void gerc(C alpha,const Vector<R,C> &x,const Vector<R,C> &y) {
      OBSOLETE("inappropriate for this class");
    }
    void geru(C alpha,const Vector<R,C> &x,const Vector<R,C> &y) {
      OBSOLETE("inappropriate for this class");
    }
    void gemm(char transa,char transb,C alpha,const Matrix<R,C> &A,
    const Matrix<R,C> &B,C beta) {
      OBSOLETE("inappropriate for this class");
    }
    R equilibrate(Vector<R,R> &r,Vector<R,R> &c,R &rowcnd,R &colcnd)
    const {
      OBSOLETE("inappropriate for this class");
    }
    virtual void copyFrom(char,int,int,const Matrix<R,C>&) {
      OBSOLETE("inappropriate for this class");
    }
    void set(char uplo,C diag,C offdiag) {
      OBSOLETE("inappropriate for this class");
    }
    int lastNonzeroColumn() const {
      OBSOLETE("inappropriate for this class");
    }
    int lastNonzeroRow() const {
      OBSOLETE("inappropriate for this class");
    }
  protected:
    static C outofbounds_;
    static C safety_;
  public:
    TrapezoidalMatrix() : Matrix<R,C>() {;}
    explicit TrapezoidalMatrix(int m,int n) : Matrix<R,C>(m,n) {;}
    virtual ~TrapezoidalMatrix() {}

//  pure virtuals: must be defined by derived classes
    virtual Matrix<R,C>* makeMatrix() const = 0;
    virtual TrapezoidalMatrix& operator=(C) = 0;
    virtual C* addr() = 0;
    virtual const C* addr() const = 0;
    virtual C* addr(int,int) = 0;
    virtual const C* addr(int,int) const = 0;
    virtual C& operator()(int,int) = 0;
    virtual C operator()(int,int) const =0;
    virtual void copy(const Matrix<R,C>&) = 0;
    virtual Vector<R,C>* trmv(const Vector<R,C>&,char='N') const = 0;
//  BLAS routines
    virtual Matrix<R,C>* trmm(const Matrix<R,C>&,char='L',char='N')
      const = 0;
//  LaPack routines
    virtual void copyFrom(int,int,const Matrix<R,C>&) = 0;
    virtual void scale(R denominator,R numerator) = 0;
    virtual void set(C offdiag,C diag) = 0;
    virtual R normFrobenius() const = 0;
    virtual R normInfinity() const = 0;
    virtual R normMaxEntry() const = 0;
    virtual R normOne() const = 0;
    virtual void solve(const Vector<R,C>&,Vector<R,C>&,char='T') const = 0;
    virtual void solve(const Matrix<R,C>&,Matrix<R,C>&,char='L',char='T')
      const = 0;

    virtual void printOn(ostream&) const = 0;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class UpperTrapezoidalMatrix;

template<typename R,typename C> class LowerTrapezoidalMatrix : 
public TrapezoidalMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit LowerTrapezoidalMatrix(const LowerTrapezoidalMatrix &L) : 
    TrapezoidalMatrix<R,C>(L.size(0),L.size(1)) {
      copyFrom(L.size(0),L.size(1),L);
    }
// Use copy instead
    LowerTrapezoidalMatrix& operator=(const LowerTrapezoidalMatrix &L) {
      int m=L.size(0),n=L.size(1);
      Matrix<R,C>::resize(m,n);
      copyFrom(m,n,L);
      return *this;
    }
//  prevent Matrix function
    virtual void copyFrom(char,int m,int n,const Matrix<R,C> &L) {
      OBSOLETE("inappropriate for this class");
    }
  protected:
    virtual void fillWith(C d) { Matrix<R,C>::set('L',d,d); }
  public:
    LowerTrapezoidalMatrix() : TrapezoidalMatrix<R,C>() {}
    explicit LowerTrapezoidalMatrix(int m,int n) :
    TrapezoidalMatrix<R,C>(m,n) {
      CHECK_TEST(m>=n);
    }
    explicit LowerTrapezoidalMatrix(int m,int n,C d) : 
    TrapezoidalMatrix<R,C>(m,n) {
      CHECK_TEST(m>=n);
      fillWith(d);
    }
    virtual ~LowerTrapezoidalMatrix() {}
    virtual Matrix<R,C>* makeMatrix() const;

    virtual LowerTrapezoidalMatrix& operator=(C scalar) {
      fillWith(scalar);
      return *this; 
    }

    virtual C* addr() { return Matrix<R,C>::addr(); }
    virtual const C* addr() const { return Matrix<R,C>::addr(); }
    virtual C* addr(int i,int j) {
      CHECK_TEST(i>=j);
      return Matrix<R,C>::addr(i,j);
    }
    virtual const C* addr(int i,int j) const {
      CHECK_TEST(i>=j);
      return Matrix<R,C>::addr(i,j);
    }
    virtual C& operator()(int i,int j) {
      if (i>=j) return Matrix<R,C>::operator()(i,j);
      this->safety_=this->outofbounds_; 
      return this->safety_; 
    }
    virtual C operator()(int i,int j) const {
      if (i>=j) return Matrix<R,C>::operator()(i,j);
      return this->outofbounds_;
    }
    virtual void copy(const Matrix<R,C> &L) {
      int m=this->size(0),n=this->size(1);
      CHECK_SAME(m,L.size(0));
      CHECK_SAME(n,L.size(1));
      copyFrom(m,n,L);
    }
    LowerTrapezoidalMatrix& operator+=(const LowerTrapezoidalMatrix&);
    LowerTrapezoidalMatrix& operator-=(const LowerTrapezoidalMatrix&);
    virtual LowerTrapezoidalMatrix& operator*=(C);
    virtual LowerTrapezoidalMatrix& operator/=(C);

//  these return Matrix* to remind us that the return value
//    must be deleted by the calling procedure, to avoid memory leaks
    virtual LowerTrapezoidalMatrix* operator+(
      const LowerTrapezoidalMatrix &L) const;
    virtual Matrix<R,C>* operator+(const Matrix<R,C> &M) const;

    virtual LowerTrapezoidalMatrix* operator-(
      const LowerTrapezoidalMatrix &L) const;
    virtual Matrix<R,C>* operator-(const Matrix<R,C> &M) const;

    virtual LowerTrapezoidalMatrix* operator*(C c) const {
      LowerTrapezoidalMatrix *L=OPERATOR_NEW LowerTrapezoidalMatrix(*this);
      *L*=c;
      return L;
    }
    virtual LowerTrapezoidalMatrix* operator/(C c) const {
      LowerTrapezoidalMatrix *L=OPERATOR_NEW LowerTrapezoidalMatrix(*this);
      *L/=c;
      return L;
    }
    virtual LowerTrapezoidalMatrix* operator*(
      const LowerTrapezoidalMatrix<R,C>&) const;
    virtual Matrix<R,C>* operator*(const Matrix<R,C>&) const;
    virtual Vector<R,C>* operator*(const Vector<R,C>&) const;

//  use trmv or trmm instead of transpose:
//  virtual UpperTrapezoidalMatrix<R,C>* transpose() const;
//  virtual UpperTrapezoidalMatrix<R,C>* conjugateTranspose() const;

//  BLAS routines
    virtual Vector<R,C>* trmv(const Vector<R,C>&,char trans='N') const;
     // op(THIS)*x
    virtual Matrix<R,C>* trmm(const Matrix<R,C>&,char side='L',
      char trans='N') const; // op(THIS) * B or B * op(THIS)
//  LaPack routines
    virtual void copyFrom(int m,int n,const Matrix<R,C> &L) {
      Matrix<R,C>::copyFrom('L',m,n,L);
    }
    virtual void scale(R denominator,R numerator) {
      Matrix<R,C>::scale('L',0,0,denominator,numerator);
    }
    virtual void set(C offdiag,C diag) {
      Matrix<R,C>::set('L',offdiag,diag);
    }
    virtual R normFrobenius() const;
    virtual R normInfinity() const;
    virtual R normMaxEntry() const;
    virtual R normOne() const;

    void solve(const Matrix<R,C>&,Matrix<R,C>&,char side='L',char trans='N')
      const;
    void solve(const Vector<R,C>&,Vector<R,C>&,char trans='N') const;

    virtual void printOn(ostream&) const ;
};
//  these return Matrix* to remind us that the return value
//    must be deleted by the calling procedure, to avoid memory leaks
template<typename R,typename C> Matrix<R,C>* operator+(
  const Matrix<R,C> &M,const LowerTrapezoidalMatrix<R,C> &L);
template<typename R,typename C> Matrix<R,C>* operator-(
  const Matrix<R,C> &M,const LowerTrapezoidalMatrix<R,C> &L);
template<typename R,typename C> Matrix<R,C>* operator*(
  const Matrix<R,C>&,const LowerTrapezoidalMatrix<R,C>&);

template<typename R,typename C> void testLowerTrapezoidalMatrix(R,C);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class LowerTriangularMatrix : 
public LowerTrapezoidalMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit LowerTriangularMatrix(const LowerTriangularMatrix &L) : 
    LowerTrapezoidalMatrix<R,C>(L.size(0),L.size(0)) { 
      LowerTrapezoidalMatrix<R,C>::copy(L);
    }
// Use copy instead
    LowerTriangularMatrix& operator=(const LowerTriangularMatrix &L) {
      Matrix<R,C>::resize(L.size(0),L.size(0)); 
      LowerTrapezoidalMatrix<R,C>::copy(L); 
      return *this;
    }
  protected:
  public:
    LowerTriangularMatrix() : LowerTrapezoidalMatrix<R,C>() {;}
    explicit LowerTriangularMatrix(int n) :
      LowerTrapezoidalMatrix<R,C>(n,n) {}
    explicit LowerTriangularMatrix(int n,C d) : 
      LowerTrapezoidalMatrix<R,C>(n,n,d) {}
    virtual ~LowerTriangularMatrix() {}

    virtual LowerTriangularMatrix& operator=(C scalar) {
      LowerTrapezoidalMatrix<R,C>::operator=(scalar);
      return *this; 
    }

    R reciprocalConditionNumber(char) const;

    virtual void printOn(ostream&) const ;
};
template<typename R,typename C> void testLowerTriangularMatrix(R,C);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class UnitUpperTrapezoidalMatrix;

template<typename R,typename C> class UnitLowerTrapezoidalMatrix : 
public LowerTrapezoidalMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit
      UnitLowerTrapezoidalMatrix(const UnitLowerTrapezoidalMatrix &L) : 
      LowerTrapezoidalMatrix<R,C>(L.size(0),L.size(1)) { copy(L);}
//  use copy instead
    UnitLowerTrapezoidalMatrix& operator=(
    const UnitLowerTrapezoidalMatrix &L) {
      int m=L.size(0),n=L.size(1);
      Matrix<R,C>::resize(m,n);
      copyFrom(m,n,L);
      return *this;
    }
//  prevent LowerTrapezoidalMatrix functions
    virtual void copyFrom(char,int m,int n,const Matrix<R,C> &L) {
      OBSOLETE("inappropriate for this class");
    }
    virtual UnitLowerTrapezoidalMatrix& operator*=(C) {
      OBSOLETE("inappropriate for this class"); return *this;
    }
    virtual UnitLowerTrapezoidalMatrix& operator/=(C) {
      OBSOLETE("inappropriate for this class"); return *this;
    }
    virtual void scale(R denominator,R numerator) {
      OBSOLETE("inappropriate for this class");
    }
    virtual void set(C offdiag,C diag) {
      OBSOLETE("inappropriate for this class");
    }
  protected:
    virtual void fillWith(C);
  public:
    UnitLowerTrapezoidalMatrix() : LowerTrapezoidalMatrix<R,C>() {;}
    explicit UnitLowerTrapezoidalMatrix(int m,int n) : 
      LowerTrapezoidalMatrix<R,C>(m,n) {}
    explicit UnitLowerTrapezoidalMatrix(int m,int n,C d) :
      LowerTrapezoidalMatrix<R,C>(m,n) { fillWith(d); }
    virtual ~UnitLowerTrapezoidalMatrix() {}

    virtual UnitLowerTrapezoidalMatrix& operator=(C scalar) { 
      fillWith(scalar); return *this; 
    }
    virtual void copy(const Matrix<R,C> &L) {
      int m=this->size(0),n=this->size(1);
      CHECK_SAME(m,L.size(0));
      CHECK_SAME(n,L.size(1));
      copyFrom(m,n,L);
    }

    virtual C* addr() { return Matrix<R,C>::addr(); }
    virtual const C* addr() const { return Matrix<R,C>::addr(); }
    virtual C* addr(int i,int j) {
      CHECK_TEST(i>=j); // calling routine must not use diagonal entries
      return Matrix<R,C>::addr(i,j);
    }
    virtual const C* addr(int i,int j) const {
      CHECK_TEST(i>=j); // calling routine must not use diagonal entries
      return Matrix<R,C>::addr(i,j);
    }
    virtual C& operator()(int i,int j) {
      if (i>j) return Matrix<R,C>::operator()(i,j);
      if (i==j) {
        this->safety_=this->one_; 
        return this->safety_; 
      }
      this->safety_=this->outofbounds_; 
      return this->safety_; 
    }
    virtual C operator()(int i,int j) const {
      if (i>j) return Matrix<R,C>::operator()(i,j);
      if (i==j) return this->one_;
      return this->outofbounds_;
    }
    virtual LowerTrapezoidalMatrix<R,C>* operator+(
      const LowerTrapezoidalMatrix<R,C>&) const;
    virtual Matrix<R,C>* operator+(const Matrix<R,C>&) const;
    virtual LowerTrapezoidalMatrix<R,C>* operator-(
      const LowerTrapezoidalMatrix<R,C>&) const;
    virtual Matrix<R,C>* operator-(const Matrix<R,C>&) const;
    virtual LowerTrapezoidalMatrix<R,C>* operator*(C) const;
    virtual LowerTrapezoidalMatrix<R,C>* operator/(C) const;
    virtual LowerTrapezoidalMatrix<R,C>* operator*(
      const LowerTrapezoidalMatrix<R,C>&) const;
    virtual Matrix<R,C>* operator*(const Matrix<R,C>&) const;
    virtual Vector<R,C>* operator*(const Vector<R,C>&) const;

//  use trmv or trmm instead of transpose:
//  virtual UnitUpperTrapezoidalMatrix<R,C>* transpose() const;
//  virtual UnitUpperTrapezoidalMatrix<R,C>* conjugateTranspose() const;

//  BLAS routines
    virtual Vector<R,C>* trmv(const Vector<R,C>&,char='N') const;
     // op(THIS) * x
    virtual Matrix<R,C>* trmm(const Matrix<R,C>&,char side='L',
      char trans='N') const; // op(THIS) * B or B * op(THIS)
//  LaPack routines
    virtual void copyFrom(int,int,const Matrix<R,C>&);
    virtual R normFrobenius() const;
    virtual R normInfinity() const;
    virtual R normMaxEntry() const;
    virtual R normOne() const;

    void solve(const Matrix<R,C>&,Matrix<R,C>&,char side='L',char trans='N')
      const;
    void solve(const Vector<R,C>&,Vector<R,C>&,char trans='N') const;
    virtual void printOn(ostream&) const ;
};
template<typename R,typename C> void testUnitLowerTrapezoidalMatrix(R,C);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//although a UnitLowerTriangularMatrix is both a LowerTriangularMatrix
//  and a UnitLowerTrapezoidalMatrix, giving this class multiple
//  inheritance would cause programming problems
//we include this class so that reciprocalConditionNumber is defined
//  without accessing diagonal entries.
template<typename R,typename C> class UnitLowerTriangularMatrix : 
public UnitLowerTrapezoidalMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit UnitLowerTriangularMatrix(
    const UnitLowerTriangularMatrix &L) : 
    UnitLowerTrapezoidalMatrix<R,C>(L.size(0),L.size(0)) { 
      UnitLowerTrapezoidalMatrix<R,C>::copy(L);
    }
// Use copy instead
    UnitLowerTriangularMatrix& operator=(
    const UnitLowerTriangularMatrix &L) {
      Matrix<R,C>::resize(L.size(0),L.size(0)); 
      UnitLowerTrapezoidalMatrix<R,C>::copy(L); 
      return *this;
    }
  protected:
  public:
    UnitLowerTriangularMatrix() : UnitLowerTrapezoidalMatrix<R,C>() {;}
    explicit UnitLowerTriangularMatrix(int n) :
      UnitLowerTrapezoidalMatrix<R,C>(n,n) {}
    explicit UnitLowerTriangularMatrix(int n,C d) : 
      UnitLowerTrapezoidalMatrix<R,C>(n,n,d) {}
    virtual ~UnitLowerTriangularMatrix() {}

    virtual UnitLowerTriangularMatrix& operator=(C scalar) {
      UnitLowerTrapezoidalMatrix<R,C>::operator=(scalar);
      return *this; 
    }

    R reciprocalConditionNumber(char) const;

    virtual void printOn(ostream&) const ;
};
template<typename R,typename C> void testUnitLowerTriangularMatrix(R,C);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class UpperTrapezoidalMatrix : 
public TrapezoidalMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit UpperTrapezoidalMatrix(const UpperTrapezoidalMatrix &U) : 
    TrapezoidalMatrix<R,C>(U.size(0),U.size(1)) { 
      copyFrom(U.size(0),U.size(1),U);
    }
// Use copy instead
    UpperTrapezoidalMatrix& operator=(const UpperTrapezoidalMatrix &U) {
      int m=U.size(0),n=U.size(1);
      Matrix<R,C>::resize(m,n);
      copyFrom(m,n,U);
      return *this;
    }
//  prevent Matrix functions
    virtual void copyFrom(char,int,int,const Matrix<R,C>&) {
      OBSOLETE("inappropriate for this class");
    }
  protected:
    virtual void fillWith(C d) { Matrix<R,C>::set('U',d,d); }
  public:
    UpperTrapezoidalMatrix() : TrapezoidalMatrix<R,C>() {}
    explicit UpperTrapezoidalMatrix(int m,int n) :
    TrapezoidalMatrix<R,C>(m,n) {
      CHECK_TEST(m<=n);
    }
    explicit UpperTrapezoidalMatrix(int m,int n,C d) : 
    TrapezoidalMatrix<R,C>(m,n) {
      CHECK_TEST(m<=n);
      fillWith(d);
    }
    virtual ~UpperTrapezoidalMatrix() {}
    virtual Matrix<R,C>* makeMatrix() const;

    virtual UpperTrapezoidalMatrix& operator=(C scalar) { 
      fillWith(scalar);
      return *this;
    }

    virtual C* addr() { return Matrix<R,C>::addr(); }
    virtual const C* addr() const { return Matrix<R,C>::addr(); }
    virtual C* addr(int i,int j) {
      CHECK_TEST(i<=j);
      return Matrix<R,C>::addr(i,j);
    }
    virtual const C* addr(int i,int j) const {
      CHECK_TEST(i<=j);
      return Matrix<R,C>::addr(i,j);
    }
    virtual C& operator()(int i,int j) {
      if (i<=j) return Matrix<R,C>::operator()(i,j);
      this->safety_=this->outofbounds_;
      return this->safety_; 
    }
    virtual C operator()(int i,int j) const {
      if (i<=j) return Matrix<R,C>::operator()(i,j);
      return this->outofbounds_;
    }
    virtual void copy(const Matrix<R,C> &U) {
      int m=this->size(0),n=this->size(1);
      CHECK_SAME(m,U.size(0));
      CHECK_SAME(n,U.size(1));
      copyFrom(m,n,U);
    }

    UpperTrapezoidalMatrix& operator+=(const UpperTrapezoidalMatrix<R,C>&);
    UpperTrapezoidalMatrix& operator-=(const UpperTrapezoidalMatrix<R,C>&);
    virtual UpperTrapezoidalMatrix& operator*=(C);
    virtual UpperTrapezoidalMatrix& operator/=(C);

    virtual UpperTrapezoidalMatrix* operator+(
      const UpperTrapezoidalMatrix<R,C> &U) const;
    virtual SquareMatrix<R,C>* operator+(
      const LowerTrapezoidalMatrix<R,C> &L) const;
    virtual Matrix<R,C>* operator+(const Matrix<R,C> &M) const;

    virtual UpperTrapezoidalMatrix* operator-(
      const UpperTrapezoidalMatrix<R,C> &U) const;
    virtual SquareMatrix<R,C>* operator-(
      const LowerTrapezoidalMatrix<R,C> &L) const;
    virtual Matrix<R,C>* operator-(const Matrix<R,C> &M) const;

    virtual UpperTrapezoidalMatrix* operator*(C z) const {
      UpperTrapezoidalMatrix *U=
        OPERATOR_NEW UpperTrapezoidalMatrix(*this);
      *U*=z;
      return U;
    }
    virtual UpperTrapezoidalMatrix* operator/(C z) const {
      UpperTrapezoidalMatrix *U=
        OPERATOR_NEW UpperTrapezoidalMatrix(*this);
      *U/=z;
      return U;
    }

    virtual UpperTrapezoidalMatrix<R,C>* operator*(
      const UpperTrapezoidalMatrix<R,C>&) const;
    virtual Matrix<R,C>* operator*(
      const LowerTrapezoidalMatrix<R,C>&) const;
    virtual Matrix<R,C>* operator*(const Matrix<R,C>&) const;
    virtual Vector<R,C>* operator*(const Vector<R,C>&) const;

//  use trmv or trmm instead of transpose:
//  virtual LowerTrapezoidalMatrix<R,C>* transpose() const;
//  virtual LowerTrapezoidalMatrix<R,C>* conjugateTranspose() const;

//  BLAS routines
    virtual Vector<R,C>* trmv(const Vector<R,C>&,char trans='N') const;
     // op(THIS) * x
    virtual Matrix<R,C>* trmm(const Matrix<R,C>&,char side='L',
      char trans='N') const; // OP(THIS) * B or B * op(THIS)
//  LaPack routines
    virtual void copyFrom(int m,int n,const Matrix<R,C> &U) {
      Matrix<R,C>::copyFrom('U',m,n,U);
    }
    virtual void scale(R denominator,R numerator) {
      Matrix<R,C>::scale('U',0,0,denominator,numerator);
    }
    virtual void set(C offdiag,C diag) {
      Matrix<R,C>::set('U',offdiag,diag);
    }
    virtual R normFrobenius() const;
    virtual R normInfinity() const;
    virtual R normMaxEntry() const;
    virtual R normOne() const;
    virtual void solve(const Matrix<R,C>&,Matrix<R,C>&,char side='L',
      char trans='N') const;
    virtual void solve(const Vector<R,C>&,Vector<R,C>&,char trans='N')
      const;

    virtual void printOn(ostream&) const ;
};

template<typename R,typename C> SquareMatrix<R,C>* operator+(
  const UnitLowerTrapezoidalMatrix<R,C> &L,
  const UpperTrapezoidalMatrix<R,C> &U);
template<typename R,typename C> SquareMatrix<R,C>* operator+(
  const LowerTrapezoidalMatrix<R,C> &L,
  const UpperTrapezoidalMatrix<R,C> &U);
template<typename R,typename C> Matrix<R,C>* operator+(
  const Matrix<R,C> &M,const UpperTrapezoidalMatrix<R,C> &U);

template<typename R,typename C> SquareMatrix<R,C>* operator-(
  const UnitLowerTrapezoidalMatrix<R,C> &L,
  const UpperTrapezoidalMatrix<R,C> &U);
template<typename R,typename C> SquareMatrix<R,C>* operator-(
  const LowerTrapezoidalMatrix<R,C> &L,
  const UpperTrapezoidalMatrix<R,C> &U);
template<typename R,typename C> Matrix<R,C>* operator-(
  const Matrix<R,C> &M,const UpperTrapezoidalMatrix<R,C> &U);

template<typename R,typename C> Matrix<R,C>* operator*(
  const UnitLowerTrapezoidalMatrix<R,C>&,
  const UpperTrapezoidalMatrix<R,C>&);
template<typename R,typename C> Matrix<R,C>* operator*(
  const LowerTrapezoidalMatrix<R,C>&,const UpperTrapezoidalMatrix<R,C>&);
template<typename R,typename C> Matrix<R,C>* operator*(
  const Matrix<R,C>&,const UpperTrapezoidalMatrix<R,C>&);

template<typename R,typename C> void testUpperTrapezoidalMatrix(R,C);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class UpperTriangularMatrix : 
public UpperTrapezoidalMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit UpperTriangularMatrix(const UpperTriangularMatrix &U) : 
    UpperTrapezoidalMatrix<R,C>(U.size(0),U.size(0)) { 
      UpperTrapezoidalMatrix<R,C>::copy(U);
    }
// Use copy instead
    UpperTriangularMatrix& operator=(const UpperTriangularMatrix &U) {
      Matrix<R,C>::resize(U.size(0),U.size(0)); 
      UpperTrapezoidalMatrix<R,C>::copy(U); 
      return *this;
    }
  protected:
  public:
    UpperTriangularMatrix() : UpperTrapezoidalMatrix<R,C>() {;}
    explicit UpperTriangularMatrix(int n) :
      UpperTrapezoidalMatrix<R,C>(n,n) {}
    explicit UpperTriangularMatrix(int n,C d) : 
      UpperTrapezoidalMatrix<R,C>(n,n,d) {}
    virtual ~UpperTriangularMatrix() {}

    virtual UpperTriangularMatrix& operator=(C scalar) {
      UpperTrapezoidalMatrix<R,C>::operator=(scalar);
      return *this;
    }
    R reciprocalConditionNumber(char) const;
//  UpperTriangularMatrix* inverse() const;

    virtual void printOn(ostream&) const ;
};
template<typename R,typename C> void testUpperTriangularMatrix(R,C);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class UnitUpperTrapezoidalMatrix : 
public UpperTrapezoidalMatrix<R,C> {
  private:
    explicit
      UnitUpperTrapezoidalMatrix(const UnitUpperTrapezoidalMatrix &U) : 
        UpperTrapezoidalMatrix<R,C>(U.size(0),U.size(1)) { copy(U);}
//  use copy instead
    UnitUpperTrapezoidalMatrix& operator=(
    const UnitUpperTrapezoidalMatrix &U) {
      int m=U.size(0),n=U.size(1);
      Matrix<R,C>::resize(m,n);
      copyFrom(m,n,U);
      return *this;
    }
//  prevent UpperTrapezoidalMatrix functions
    virtual UnitUpperTrapezoidalMatrix& operator*=(C) {
      OBSOLETE("inappropriate for this class"); return *this;
    }
    virtual UnitUpperTrapezoidalMatrix& operator/=(C) {
      OBSOLETE("inappropriate for this class"); return *this;
    }
    virtual void scale(R denominator,R numerator) {
      OBSOLETE("inappropriate for this class");
    }
    virtual void set(C offdiag,C diag) {
      OBSOLETE("inappropriate for this class");
    }
    virtual void copyFrom(char,int,int,const Matrix<R,C>&) {
      OBSOLETE("inappropriate for this class");
    }
  protected:
    virtual void fillWith(C);
  public:
    UnitUpperTrapezoidalMatrix() : UpperTrapezoidalMatrix<R,C>() {;}
    explicit UnitUpperTrapezoidalMatrix(int m,int n) : 
      UpperTrapezoidalMatrix<R,C>(m,n) {}
    explicit UnitUpperTrapezoidalMatrix(int m,int n,C d) :
      UpperTrapezoidalMatrix<R,C>(m,n) { fillWith(d); }
    virtual ~UnitUpperTrapezoidalMatrix() {}

    virtual UnitUpperTrapezoidalMatrix& operator=(C scalar) { 
      fillWith(scalar); return *this; 
    }
    virtual void copy(const Matrix<R,C> &U) {
      int m=this->size(0),n=this->size(1);
      CHECK_SAME(m,U.size(0));
      CHECK_SAME(n,U.size(1));
      copyFrom(m,n,U);
    }
    virtual void copyFrom(int,int,const Matrix<R,C>&);

    virtual C* addr() { return Matrix<R,C>::addr(); }
    virtual const C* addr() const { return Matrix<R,C>::addr(); }
    virtual C* addr(int i,int j) {
      CHECK_TEST(i<=j);// calling routine must not use diagonal entries
      return Matrix<R,C>::addr(i,j);
    }
    virtual const C* addr(int i,int j) const {
      CHECK_TEST(i<=j);// calling routine must not use diagonal entries
      return Matrix<R,C>::addr(i,j);
    }
    virtual C& operator()(int i,int j) {
      if (i<j) return Matrix<R,C>::operator()(i,j);
      if (i==j) {
        this->safety_=this->one_; 
        return this->safety_; 
      }
      this->safety_=this->outofbounds_; 
      return this->safety_; 
    }
    virtual C operator()(int i,int j) const {
      if (i<j) return Matrix<R,C>::operator()(i,j);
      if (i==j) return this->one_;
      return this->outofbounds_;
    }
    virtual UpperTrapezoidalMatrix<R,C>* operator+(
      const UpperTrapezoidalMatrix<R,C>&) const;
    virtual SquareMatrix<R,C>* operator+(
      const LowerTrapezoidalMatrix<R,C> &L) const;
    virtual Matrix<R,C>* operator+(const Matrix<R,C>&) const;

    virtual UpperTrapezoidalMatrix<R,C>* operator-(
      const UpperTrapezoidalMatrix<R,C>&) const;
    virtual SquareMatrix<R,C>* operator-(
      const LowerTrapezoidalMatrix<R,C> &L) const;
    virtual Matrix<R,C>* operator-(const Matrix<R,C>&) const;

    virtual UpperTrapezoidalMatrix<R,C>* operator*(C) const;
    virtual UpperTrapezoidalMatrix<R,C>* operator/(C) const;

    virtual UpperTrapezoidalMatrix<R,C>* operator*(
      const UpperTrapezoidalMatrix<R,C>&) const ;
    virtual Matrix<R,C>* operator*(const LowerTrapezoidalMatrix<R,C>&)
      const;
    virtual Matrix<R,C>* operator*(const Matrix<R,C>&) const;
    virtual Vector<R,C>* operator*(const Vector<R,C>&) const;

//  use trmv or trmm instead of transpose:
//  virtual UnitLowerTrapezoidalMatrix<R,C>* transpose() const;
//  virtual UnitLowerTrapezoidalMatrix<R,C>* conjugateTranspose() const;

//  BLAS routines
    virtual Vector<R,C>* trmv(const Vector<R,C>&,char trans='N') const;
     // op(THIS) * x
    virtual Matrix<R,C>* trmm(const Matrix<R,C>&,char side='L',
      char trans='N') const; // op(THIS) * B or B * op(THIS)
//  LaPack routines
    virtual R normFrobenius() const;
    virtual R normInfinity() const;
    virtual R normMaxEntry() const;
    virtual R normOne() const;
    virtual void solve(const Matrix<R,C>&,Matrix<R,C>&,char side='L',
      char trans='N') const;
    virtual void solve(const Vector<R,C>&,Vector<R,C>&,char trans='N')
      const;

    virtual void printOn(ostream&) const ;
};
template<typename R,typename C> void testUnitUpperTrapezoidalMatrix(R,C);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename R,typename C> class UnitUpperTriangularMatrix : 
public UnitUpperTrapezoidalMatrix<R,C> {
  private:
//  prevent public use of copy constructor
    explicit UnitUpperTriangularMatrix(const UnitUpperTriangularMatrix &U) :
    UnitUpperTrapezoidalMatrix<R,C>(U.size(0),U.size(0)) { 
      UnitUpperTrapezoidalMatrix<R,C>::copy(U);
    }
// Use copy instead
    UnitUpperTriangularMatrix& operator=(
    const UnitUpperTriangularMatrix &U) {
      Matrix<R,C>::resize(U.size(0),U.size(0)); 
      UnitUpperTrapezoidalMatrix<R,C>::copy(U); 
      return *this;
    }
  protected:
  public:
    UnitUpperTriangularMatrix() : UnitUpperTrapezoidalMatrix<R,C>() {;}
    explicit UnitUpperTriangularMatrix(int n) :
      UnitUpperTrapezoidalMatrix<R,C>(n,n) {}
    explicit UnitUpperTriangularMatrix(int n,C d) : 
      UnitUpperTrapezoidalMatrix<R,C>(n,n,d) {}
    virtual ~UnitUpperTriangularMatrix() {}

    virtual UnitUpperTriangularMatrix& operator=(C scalar) {
      UnitUpperTrapezoidalMatrix<R,C>::operator=(scalar);
      return *this;
    }
    R reciprocalConditionNumber(char) const;
//  UnitUpperTriangularMatrix* inverse() const;

    virtual void printOn(ostream&) const ;
};
template<typename R,typename C> void testUnitUpperTriangularMatrix(R,C);
#endif 
