#include "CholeskyFactorization.H"
#include "BandMatrix.H"
#include "SymmetricMatrix.H"

template<typename F,typename Z,typename M> void
CholeskyFactorization<F,Z,M>::printOn(ostream &os) const {
  os << "CholeskyFactorization: equilibrate option = " << equ_op << endl;
  if (L != 0) {
    os << "\treciprocal pivot growth = " << rpvgrw << endl;
    os << "\tL:" << endl;
    L->printOn(cout);
  }
  if (s != 0) {
    os << "\tscond = " << scond << endl;
    os << "\tscale factors:" << endl;
    s->printOn(cout);
  }
}

template<typename F,typename Z> void
testCholeskyFactorization(F fscalar,Z scalar) {
  SymmetricPositiveMatrix<F,Z> A(3);
  Vector<F,Z> b(3);
  Vector<F,Z> p(3,1.e-10);
  A(0,0)= 30.*fscalar;
  A(1,0)=static_cast<F>(-13.)*scalar;
    A(1,1)= 13.*fscalar;
  A(2,0)=static_cast<F>(20.)*scalar;
    A(2,1)=static_cast<F>(-8.)*scalar;
    A(2,2)=16.*fscalar;
  b[0]=   scalar;
  b[1]=static_cast<F>(2.)*scalar;
  b[2]=static_cast<F>(3.)*scalar;
/*
  {
    Vector<F,Z> x(3);
    Vector<F,Z> r(3);
    for (Factorization::EQUILIBRATE_OPTION
    eo=Factorization::NO_EQUILIBRATION;;
    eo=Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      CholeskyFactorization<F,Z,SymmetricPositiveMatrix<F,Z> > CF(A,eo);
      cout << "\neo = " << eo << endl;
      CF.printOn(cout);
      cout << "reciprocal condition number = "
           << CF.reciprocalConditionNumber() << endl;
      cout << "reciprocal pivot growth = " << CF.reciprocalPivotGrowth()
           << endl;

      CF.solve(b,x);
      cout << "b = " << endl;
      b.printOn(cout);
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      A.symv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
      cout << "r = " << endl;
      r.printOn(cout);

      F berr,ferr;
      x+=p;
      CF.improve(b,x,berr,ferr);
      cout << "berr,ferr = " << berr << " " << ferr << endl;
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      A.symv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
      cout << "r = " << endl;
      r.printOn(cout);

      if (eo==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) break;
    }
  }
*/
/*
  {
    Matrix<F,Z> BL(3,2),BR(2,3);
    Matrix<F,Z> XL(3,2),XR(2,3);
    Matrix<F,Z> PL(3,2,1.e-10),PR(2,3,1.e-10);
    Matrix<F,Z> RL(3,2),RR(2,3);
    BL(0,0)=BL(0,1)=BR(0,0)=BR(1,0)=b[0];
    BL(1,0)=BL(1,1)=BR(0,1)=BR(1,1)=b[1];
    BL(2,0)=BL(2,1)=BR(0,2)=BR(1,2)=b[2];
//  for (int j=0;j<3;j++) {
//    A(0,j)*=3.e15;
//  }
//  BL(0,0)*=3.e15;
//  BR(0,0)*=3.e15;
//  for (int i=0;i<3;i++) {
//    A(i,1)/=6.e12;
//    A(i,2)/=9.e12;
//  }
//  A.printOn(cout);
//  BL.printOn(cout);
//  BR.printOn(cout);

    Vector<F,F> berr(2),ferr(2);
    for (Factorization::EQUILIBRATE_OPTION
    eo=Factorization::NO_EQUILIBRATION;;
    eo=Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      CholeskyFactorization<F,Z,SymmetricPositiveMatrix<F,Z> > CF(A,eo);
      cout << "\neo = " << eo << endl;
      CF.printOn(cout);
      cout << "reciprocal condition number = "
           << CF.reciprocalConditionNumber() << endl;
      cout << "reciprocal pivot growth = " << CF.reciprocalPivotGrowth()
           << endl;

      for (Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE;;
      so=Factorization::RIGHT_SIDE) {
//      Factorization::SIDE_OPTION so=Factorization::RIGHT_SIDE;
        Matrix<F,Z> &BB=(so==Factorization::LEFT_SIDE ? BL : BR);
        Matrix<F,Z> &XX=(so==Factorization::LEFT_SIDE ? XL : XR);
        Matrix<F,Z> &R2=(so==Factorization::LEFT_SIDE ? RL : RR);
        cout << "\nso = " << so << endl;
        CF.solve(BB,XX,so);
        cout << "BB = " << endl;
        BB.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          A.symm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L');
        } else {
          A.symm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);

        XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
        CF.improve(BB,XX,berr,ferr,so);
        cout << "berr = " << endl;
        berr.printOn(cout);
        cout << "ferr = " << endl;
        ferr.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          A.symm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L');
        } else {
          A.symm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);
        if (so==Factorization::RIGHT_SIDE) break;
      }
      if (eo==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) break;
    }
  }
*/
  SymmetricPositiveTridiagonalMatrix<F,Z> T(3);
  T.diagonalValue(0)=3.*fscalar; T.lowerDiagonalValue(0)=-scalar;
  T.diagonalValue(1)=4.*fscalar;
    T.lowerDiagonalValue(1)=static_cast<F>(-2.)*scalar;
  T.diagonalValue(2)=5.*fscalar;
//T.printOn(cout);
///.printOn(cout);
/*
  {
    Vector<F,Z> x(3);
    Vector<F,Z> p(3,1.e-10);
    Vector<F,Z> r(3);
    Factorization::EQUILIBRATE_OPTION eo=Factorization::NO_EQUILIBRATION;
      CholeskyFactorization<F,Z,SymmetricPositiveTridiagonalMatrix<F,Z> >
        CF(T,eo);
      cout << "\neo = " << eo << endl;
      CF.printOn(cout);
      cout << "reciprocal condition number = "
           << CF.reciprocalConditionNumber() << endl;
      F berr,ferr;

      CF.solve(b,x);
      cout << "b = " << endl;
      b.printOn(cout);
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      T.stmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
      cout << "r = " << endl;
      r.printOn(cout);

      x+=p;
      CF.improve(b,x,berr,ferr);
      cout << "berr,ferr = " << berr << " " << ferr << endl;
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      T.stmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
      cout << "r = " << endl;
      r.printOn(cout);
  }
*/
/*
  {
    Matrix<F,Z> BL(3,2),BR(2,3);
    Matrix<F,Z> XL(3,2),XR(2,3);
    Matrix<F,Z> PL(3,2,1.e-10),PR(2,3,1.e-10);
    Matrix<F,Z> RL(3,2),RR(2,3);
    BL(0,0)=BL(0,1)=BR(0,0)=BR(1,0)=b[0];
    BL(1,0)=BL(1,1)=BR(0,1)=BR(1,1)=b[1];
    BL(2,0)=BL(2,1)=BR(0,2)=BR(1,2)=b[2];
    T.printOn(cout);
    BL.printOn(cout);
    BR.printOn(cout);

    Vector<F,F> berr(2),ferr(2);
    Factorization::EQUILIBRATE_OPTION eo=Factorization::NO_EQUILIBRATION;
      CholeskyFactorization<F,Z,SymmetricPositiveTridiagonalMatrix<F,Z> >
        CF(T,eo);
      cout << "\neo = " << eo << endl;
      CF.printOn(cout);
      cout << "reciprocal condition number = "
           << CF.reciprocalConditionNumber() << endl;

      for (Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE;;
      so=Factorization::RIGHT_SIDE) {
//      Factorization::SIDE_OPTION so=Factorization::RIGHT_SIDE;
        Matrix<F,Z> &BB=(so==Factorization::LEFT_SIDE ? BL : BR);
        Matrix<F,Z> &XX=(so==Factorization::LEFT_SIDE ? XL : XR);
        Matrix<F,Z> &R2=(so==Factorization::LEFT_SIDE ? RL : RR);
        cout << "\nso = " << so << endl;

        CF.solve(BB,XX,so);
        cout << "BB = " << endl;
        BB.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          T.stmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L');
        } else {
          T.stmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);

        XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
        CF.improve(BB,XX,berr,ferr,so);
        cout << "berr = " << endl;
        berr.printOn(cout);
        cout << "ferr = " << endl;
        ferr.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          T.stmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L');
        } else {
          T.stmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);
        if (so==Factorization::RIGHT_SIDE) break;
      }
  }
*/
  SymmetricPositiveBandMatrix<F,Z> M(3,1);
  M(0,0)=3.*fscalar;
  M(1,0)=   -scalar;M(1,1)= 4.*fscalar;
                    M(2,1)=static_cast<F>(-2.)*scalar;M(2,2)= 5.*fscalar;
/*
  {
    Vector<F,Z> x(3);
    Vector<F,Z> r(3);
    for (Factorization::EQUILIBRATE_OPTION
    eo=Factorization::NO_EQUILIBRATION;;
    eo=Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      CholeskyFactorization<F,Z,SymmetricPositiveBandMatrix<F,Z> >
        CF(M,eo);
      cout << "\neo = " << eo << endl;
      cout << "reciprocal condition number = "
           << CF.reciprocalConditionNumber() << endl;
      cout << "reciprocal pivot growth = " << CF.reciprocalPivotGrowth()
           << endl;

      F berr,ferr;

      CF.solve(b,x);
      cout << "b = " << endl;
      b.printOn(cout);
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      M.sbmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
      cout << "r = " << endl;
      r.printOn(cout);

      x+=p;
      CF.improve(b,x,berr,ferr);
      cout << "berr,ferr = " << berr << " " << ferr << endl;
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      M.sbmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
      cout << "r = " << endl;
      r.printOn(cout);
      if (eo==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) break;
    }
  }
*/
//
  {
    Matrix<F,Z> BL(3,2),BR(2,3);
    Matrix<F,Z> XL(3,2),XR(2,3);
    Matrix<F,Z> PL(3,2,1.e-10),PR(2,3,1.e-10);
    Matrix<F,Z> RL(3,2),RR(2,3);
    BL(0,0)=BL(0,1)=BR(0,0)=BR(1,0)=b[0];
    BL(1,0)=BL(1,1)=BR(0,1)=BR(1,1)=b[1];
    BL(2,0)=BL(2,1)=BR(0,2)=BR(1,2)=b[2];
//  M.printOn(cout);
//  BL.printOn(cout);
//  BR.printOn(cout);

    Vector<F,F> berr(2),ferr(2);
    for (Factorization::EQUILIBRATE_OPTION
    eo=Factorization::NO_EQUILIBRATION;;
    eo=Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      CholeskyFactorization<F,Z,SymmetricPositiveBandMatrix<F,Z> >
        CF(M,eo);
      cout << "\neo = " << eo << endl;
      cout << "reciprocal condition number = "
           << CF.reciprocalConditionNumber() << endl;
      cout << "reciprocal pivot growth = " << CF.reciprocalPivotGrowth()
           << endl;

      for (Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE;;
      so=Factorization::RIGHT_SIDE) {
        Matrix<F,Z> &BB=(so==Factorization::LEFT_SIDE ? BL : BR);
        Matrix<F,Z> &XX=(so==Factorization::LEFT_SIDE ? XL : XR);
        Matrix<F,Z> &R2=(so==Factorization::LEFT_SIDE ? RL : RR);
        cout << "\nso = " << so << endl;
        CF.solve(BB,XX,so);
        cout << "BB = " << endl;
        BB.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          M.sbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L');
        } else {
          M.sbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);

        XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
        CF.improve(BB,XX,berr,ferr,so);
        cout << "berr = " << endl;
        berr.printOn(cout);
        cout << "ferr = " << endl;
        ferr.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          M.sbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L');
        } else {
          M.sbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);
        if (so==Factorization::RIGHT_SIDE) break;
      }
      if (eo==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) break;
    }
  }
//
}
     
// Modified from ltgmd.C by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
