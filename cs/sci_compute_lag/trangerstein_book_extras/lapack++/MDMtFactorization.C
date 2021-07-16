#include "MDMtFactorization.H"
#include "SymmetricMatrix.H"

template<typename F,typename Z> void
MDMtFactorization<F,Z>::printOn(ostream &os) const {
  os << "MDMtFactorization: equilibrate option = " << equ_op << endl;
  if (MD != 0) {
    os << "\treciprocal pivot growth = " << rpvgrw << endl;
    os << "\tL:" << endl;
    MD->printOn(cout);
    os << "\tipiv = ";
    for (int i=0;i<MD->size(0);i++) os << ipiv[i] << " ";
    os << endl;
  }
  if (s != 0) {
    os << "\tscond = " << scond << endl;
    os << "\tscale factors:" << endl;
    s->printOn(cout);
  }
}

template<typename F,typename Z> void
testMDMtFactorization(F fscalar,Z scalar) {
  SymmetricMatrix<F,Z> A(3);
  Vector<F,Z> b(3);
  Vector<F,Z> p(3,1.e-10);
  A(0,0)=-30.*fscalar;
  A(1,0)=static_cast<F>(-13.)* scalar;
    A(1,1)= 13.*fscalar;
  A(2,0)=static_cast<F>(20.)*scalar;
    A(2,1)=static_cast<F>(-8.)*scalar;
    A(2,2)=16.*fscalar;
  Vector<F,Z> exact(3); 
    for (int i=0;i<3;i++) exact[i]=static_cast<F>(i)*Vector<F,Z>::one_;
  b=Vector<F,Z>::zero_;
  A.symv(Vector<F,Z>::one_,exact,Vector<F,Z>::zero_,b);
/*
  {
    Vector<F,Z> x(3);
    Vector<F,Z> r(3);
    for (Factorization::EQUILIBRATE_OPTION
    eo=Factorization::NO_EQUILIBRATION;;
    eo=Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) {
      MDMtFactorization<F,Z> MF(A,eo);
      cout << "\neo = " << eo << endl;
      MF.printOn(cout);
      cout << "reciprocal condition number = "
           << MF.reciprocalConditionNumber() << endl;
      cout << "reciprocal pivot growth = " << MF.reciprocalPivotGrowth()
           << endl;

      MF.solve(b,x);
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
      MF.improve(b,x,berr,ferr);
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
//
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
      MDMtFactorization<F,Z> MF(A,eo);
      cout << "\neo = " << eo << endl;
      MF.printOn(cout);
      cout << "reciprocal condition number = "
           << MF.reciprocalConditionNumber() << endl;
      cout << "reciprocal pivot growth = " << MF.reciprocalPivotGrowth()
           << endl;

      for (Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE;;
      so=Factorization::RIGHT_SIDE) {
//      Factorization::SIDE_OPTION so=Factorization::RIGHT_SIDE;
        Matrix<F,Z> &BB=(so==Factorization::LEFT_SIDE ? BL : BR);
        Matrix<F,Z> &XX=(so==Factorization::LEFT_SIDE ? XL : XR);
        Matrix<F,Z> &R2=(so==Factorization::LEFT_SIDE ? RL : RR);
        cout << "\nso = " << so << endl;
        MF.solve(BB,XX,so);
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
        MF.improve(BB,XX,berr,ferr,so);
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
//
}
     
// Modified from ltgmd.C by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
