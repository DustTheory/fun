#ifndef _LINEAR_PROGRAM_H
#define _LINEAR_PROGRAM_H
#include "GramSchmidtQRFactorization.H"
#include "LaEnum.H"
#include "SquareMatrix.H"

//this is a base class for the classes below
//it holds information common to both treatments of the simplex method
template<class T> class VirtualLinearProgram {
  private:
  protected:
    Matrix<T,T> *A; // in original order
    Vector<T,T> *b; // in original order
    Vector<T,T> *c; // in original order
    Vector<T,T> *x; // in original order
    Vector<T,T> *y; // in original order
    STATUS_OPTION current_status;

    static T zero_;
    static T one_ ;
    static T huge_ ;
//  static int initialized_ ;
  public:
    VirtualLinearProgram(const Matrix<T,T> &Ain,const Vector<T,T> &bin,
      const Vector<T,T> &c);
    virtual ~VirtualLinearProgram() {
      delete A; A=0;
      delete b; b=0;
      delete c; c=0;
      delete x; x=0;
      delete y; y=0;
    }
    const Vector<T,T>& currentPrimalSolution() const { return *x; }
    const Vector<T,T>& currentDualSolution() const { return *y; }
    STATUS_OPTION currentStatus() const { return current_status; }
    virtual T currentValue() const;

    virtual STATUS_OPTION findBasicFeasibleGuess() = 0;
    virtual STATUS_OPTION simplexStep() = 0;

    virtual void printOn(ostream& = cerr) const ;
};
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//Standard form of the linear program:
// solve min c^t x subject to Ax=b, x>=0
template<class T> class SFLinearProgram : public VirtualLinearProgram<T> {
  protected:
    Vector<T,T> *r; // dual slack vector
    Vector<T,T> *h; // used to find basic var to leave
    GramSchmidtQRFactorization<T,T> *QR;
    SquareMatrix<T,T> *Abasic;
    Vector<T,T> *xbasic;
    Vector<T,T> *cbasic;
    Vector<T,T> *resid;
    int *column;    // the current internal order of cols
    int new_basic;

    int smallestDualSlack();
  public:
    SFLinearProgram(const Matrix<T,T> &Ain,const Vector<T,T> &bin,
      const Vector<T,T> &c_trans_in);
    ~SFLinearProgram() {
      if (r!=0) delete r; r=0; 
      if (h!=0) delete h; h=0;
      if (QR!=0) delete QR; QR=0; 
      if (Abasic!=0) delete Abasic; Abasic=0; 
      if (xbasic!=0) delete xbasic; xbasic=0; 
      if (cbasic!=0) delete cbasic; cbasic=0; 
      if (resid!=0) delete resid; resid=0; 
      if (column != 0 ) delete column; column=0;
    }
//  T currentValue() const;

    void specifyBasicVariables(int *basic_indices);
    STATUS_OPTION findBasicFeasibleGuess();
    STATUS_OPTION simplexStep();

    void costBounds(Vector<T,T> &lower,Vector<T,T> &upper) const ;
    T costSensitivity(int j,Vector<T,T> &dydcj) const ;

    void constraintBounds(Vector<T,T> &lower,Vector<T,T> &upper) const ;
    T constraintSensitivity(int i,Vector<T,T> &dxdbi) const ;

//  void arrayBounds(int i,int j,T &lower,T &upper) const ;
//  T arraySensitivity(int i,int j,Matrix<T,T> &primal_derivative,
//    Matrix<T,T> &dual_derivative) const ;

    void printOn(ostream& = cerr) const ;
};
//template<class T> ostream& operator<<(ostream &s,
//const SFLinearProgram<T> &LP){
//  LP.printOn(s); return s;
//}
template<class T> void testSFLinearProgram(T);
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Dual simplex algorithm
// solve min c^t x subject to A x >= b  , x >= 0
//    or max y^t b subject to y^t A <= c^t, y >= 0
template<class T> class LinearProgram : 
public VirtualLinearProgram<T> {
  private:
  protected:
    GramSchmidtQRFactorization<T,T> *QR;
    Vector<T,T> *bbasic; // in permuted order
    Vector<T,T> *cbasic; // in permuted order
    Vector<T,T> *xbasic; // in permuted order
    Vector<T,T> *ybasic; // in permuted order
    Vector<T,T> *r; // dual slack vector, in permuted order
    Vector<T,T> *s; // primal slack vector, in permuted order
    int *column;
    int *row;
    int basic_number;

    void computeSolution();
    STATUS_OPTION primalSimplexStep();
    STATUS_OPTION dualSimplexStep();
    void drop(int,int);
    void add(int,int);
    void switchRows(int,int);
    void switchColumns(int,int);
    int basicPrimalPivot(const Vector<T,T>&,T&);
    int nonBasicPrimalPivot(const Vector<T,T>&,T&);
    int basicDualPivot(Vector<T,T>&,T&);
    int nonBasicDualPivot(const Vector<T,T>&,T&);
  public:
    LinearProgram(const Matrix<T,T> &Ain,const Vector<T,T> &bin,
      const Vector<T,T> &cin);
    LinearProgram(const LinearProgram &LP,
      const Matrix<T,T> &new_rows_or_cols);
    ~LinearProgram() {
      if (QR!=0) delete QR; QR=0;
      if (bbasic!=0) delete bbasic; bbasic=0;
      if (cbasic!=0) delete cbasic; cbasic=0;
      if (xbasic!=0) delete xbasic; xbasic=0;
      if (ybasic!=0) delete ybasic; ybasic=0;
      if (r!=0) delete r; r=0;
      if (s!=0) delete s; s=0;
      if (column!=0) delete column; column=0;
      if (row!=0) delete row; row=0;
    }

    void specifyBasicVariables(int basic_no,int *basic_row_index,
      int *basic_column_index);
    STATUS_OPTION findBasicFeasibleGuess();
    STATUS_OPTION findDualBasicFeasibleGuess();
    STATUS_OPTION findPrimalBasicFeasibleGuess();
    STATUS_OPTION simplexStep();

//  void costBounds(Vector<T,T> &lower,Vector<T,T> &upper) const;
//  T costSensitivity(int j,Vector<T,T> &dydcj) const;
//  void constraintBounds(Vector<T,T> &lower,Vector<T,T> &upper) const;
//  T constraintSensitivity(int i,Vector<T,T> &dxdbi) const;

//  void arrayBounds(int i,int j,T &lower,T &upper) const ;
//  T arraySensitivity(int i,int j,Matrix<T,T> &primal_derivative,
//    Matrix<T,T> &dual_derivative) const ;

    void printOn(ostream& = cerr) const ;
};
//template<class T> inline ostream& operator<<(ostream &s,
//const LinearProgram<T> &LP){
//  LP.printOn(s); return s;
//}
template<class T> void testLinearProgram(T);

//template<> status_option LinearProgram<float>::findBasicFeasibleGuess();
//template<> status_option 
//  LinearProgram<float>::findPrimalBasicFeasibleGuess();
//template<> status_option 
//  LinearProgram<float>::findDualBasicFeasibleGuess();
//template<> void LinearProgram<float>::computeSolution();
//template<> int LinearProgram<float>::basicPrimalPivot(Matrix<float>&,
//  float&);
//template<> int LinearProgram<float>::basicDualPivot(Matrix<float>&,
//  float&);
//template<> void LinearProgram<float>::add(int,int);
//template<> void LinearProgram<float>::switchRows(int,int);
//template<> void LinearProgram<float>::switchColumns(int,int);
#endif
