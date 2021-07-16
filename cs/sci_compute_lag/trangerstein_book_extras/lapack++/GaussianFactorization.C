#include "GaussianFactorization.H"
#include "BandMatrix.H"
#include "SquareMatrix.H"

template<typename F,typename Z,typename M> void
GaussianFactorization<F,Z,M>::printOn(ostream &s) const {
  s << "GaussianFactorization: pivot option = " << piv_op
    << "\n\tequilibrate option = " << equ_op << endl;
//const UnitLowerTrapezoidalMatrix<F,Z> &L=lower();
//s << "L: \n" << endl;
//L.printOn(s);
//const UpperTrapezoidalMatrix<F,Z> &U=upper();
//s << "U: \n" << endl;
//U.printOn(s);
  if (LU != 0) {
    s << "\treciprocal pivot growth = " << rpvgrw << endl;
    s << "\tLU:" << endl;
    LU->printOn(cout);
  }
  if (r != 0) {
    s << "\trowcnd = " << rowcnd << endl;
    s << "\trow scale factors:" << endl;
    r->printOn(cout);
  }
  if (c != 0) {
    s << "\tcolcnd = " << colcnd << endl;
    s << "\tcolumn scale factors:" << endl;
    c->printOn(cout);
  }
  if (ipiv != 0) {
    s << "\trow pivots: \n";
    for (int i=0;i<LU->size(0)-1;i++) s << ipiv[i] << " ";
    s << endl;
  }
  if (jpiv != 0) {
    s << "\tcolumn pivots: \n";
    for (int j=0;j<LU->size(1)-1;j++) s << jpiv[j] << " ";
    s << endl;
  }
}

template<typename F,typename Z> void
testGaussianFactorization(F fscalar,Z scalar) {
  SquareMatrix<F,Z> A(3);
  Vector<F,Z> b(3);
  Vector<F,Z> p(3,1.e-10);
  A(0,0)=static_cast<F>(1.e-15)*scalar;
    A(0,1)=static_cast<F>(2.)*scalar;
    A(0,2)=static_cast<F>(3.)*scalar; b[0]=   scalar;
  A(1,0)=      -scalar;
    A(1,1)=static_cast<F>(3.)*scalar;
    A(1,2)=static_cast<F>(4.)*scalar; b[1]=static_cast<F>(2.)*scalar;
  A(2,0)=static_cast<F>(-2.)*scalar;
    A(2,1)=   scalar;
    A(2,2)=static_cast<F>(5.)*scalar; b[2]=static_cast<F>(3.)*scalar;
/*
  {
//  Stewart, Example 2.8 p. 123
    Vector<F,Z> x(3);
    Vector<F,Z> r(3);
    for (Factorization::EQUILIBRATE_OPTION
    eo=Factorization::NO_EQUILIBRATION;;eo++) {
      for (Factorization::PIVOT_OPTION po=Factorization::NO_PIVOTING;;po++)
      {
        GaussianFactorization<F,Z,SquareMatrix<F,Z> > GF(A,po,eo);
        cout << "\npo,eo = " << po << " " << eo << endl;
        GF.printOn(cout);
        cout << "reciprocal condition number = "
             << GF.reciprocalConditionNumber() << endl;
        cout << "reciprocal pivot growth = " << GF.reciprocalPivotGrowth()
             << endl;

        GF.solve(b,x);
        cout << "b = " << endl;
        b.printOn(cout);
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        A.gemv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
        cout << "r = " << endl;
        r.printOn(cout);

        F berr,ferr;
        x+=p;
        GF.improve(b,x,berr,ferr);
        cout << "berr,ferr = " << berr << " " << ferr << endl;
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        A.gemv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
        cout << "r = " << endl;
        r.printOn(cout);

        GF.solve(b,x,Factorization::TRANSPOSE);
        cout << "\nTRANSPOSE" << endl;
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        A.gemv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'T');
        cout << "r = " << endl;
        r.printOn(cout);

        x+=p;
        GF.improve(b,x,berr,ferr,Factorization::TRANSPOSE);
        cout << "berr,ferr = " << berr << " " << ferr << endl;
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        A.gemv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'T');
        cout << "r = " << endl;
        r.printOn(cout);

        if (sizeof(F)!=sizeof(Z)) {
          GF.solve(b,x,Factorization::CONJUGATE_TRANSPOSE);
          cout << "\nCONJUGATE_TRANSPOSE" << endl;
          cout << "x = " << endl;
          x.printOn(cout);
          r.copy(b);
          A.gemv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'C');
          cout << "r = " << endl;
          r.printOn(cout);

          x+=p;
          GF.improve(b,x,berr,ferr,Factorization::CONJUGATE_TRANSPOSE);
          cout << "berr,ferr = " << berr << " " << ferr << endl;
          cout << "x = " << endl;
          x.printOn(cout);
          r.copy(b);
          A.gemv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'C');
          cout << "r = " << endl;
          r.printOn(cout);
        }

        if (po==Factorization::PIVOT_ROWS_AND_COLUMNS) break;
      }
      if (eo==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) break;
    }
  }
*/
/*
  {
//  Stewart, p. 158
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
    eo=Factorization::NO_EQUILIBRATION;;eo++) {
      for (Factorization::PIVOT_OPTION po=Factorization::NO_PIVOTING;;po++)
      {
//    Factorization::PIVOT_OPTION po=Factorization::PIVOT_ROWS;
        GaussianFactorization<F,Z,SquareMatrix<F,Z> > GF(A,po,eo);
        cout << "\npo,eo = " << po << " " << eo << endl;
        GF.printOn(cout);
        cout << "reciprocal condition number = "
             << GF.reciprocalConditionNumber() << endl;
        cout << "reciprocal pivot growth = " << GF.reciprocalPivotGrowth()
             << endl;

        for (Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE;;
        so=Factorization::RIGHT_SIDE) {
//        Factorization::SIDE_OPTION so=Factorization::RIGHT_SIDE;
          Matrix<F,Z> &BB=(so==Factorization::LEFT_SIDE ? BL : BR);
          Matrix<F,Z> &XX=(so==Factorization::LEFT_SIDE ? XL : XR);
          Matrix<F,Z> &R2=(so==Factorization::LEFT_SIDE ? RL : RR);
          cout << "\nso = " << so << endl;
          GF.solve(BB,XX,Factorization::NO_TRANSPOSE,so);
          cout << "BB = " << endl;
          BB.printOn(cout);
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            R2.gemm(Matrix<F,Z>::one_,A,XX,Matrix<F,Z>::mone_,'N','N');
          } else {
            R2.gemm(Matrix<F,Z>::one_,XX,A,Matrix<F,Z>::mone_,'N','N');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);

          XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
          GF.improve(BB,XX,berr,ferr,Factorization::NO_TRANSPOSE,so);
          cout << "berr = " << endl;
          berr.printOn(cout);
          cout << "ferr = " << endl;
          ferr.printOn(cout);
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            R2.gemm(Matrix<F,Z>::one_,A,XX,Matrix<F,Z>::mone_,'N','N');
          } else {
            R2.gemm(Matrix<F,Z>::one_,XX,A,Matrix<F,Z>::mone_,'N','N');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);

          GF.solve(BB,XX,Factorization::TRANSPOSE,so);
          cout << "\nTRANSPOSE" << endl;
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            R2.gemm(Matrix<F,Z>::one_,A,XX,Matrix<F,Z>::mone_,'T','N');
          } else {
            R2.gemm(Matrix<F,Z>::one_,XX,A,Matrix<F,Z>::mone_,'N','T');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);

          XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
          GF.improve(BB,XX,berr,ferr,Factorization::TRANSPOSE,so);
          cout << "berr = " << endl;
          berr.printOn(cout);
          cout << "ferr = " << endl;
          ferr.printOn(cout);
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            R2.gemm(Matrix<F,Z>::one_,A,XX,Matrix<F,Z>::mone_,'T','N');
          } else {
            R2.gemm(Matrix<F,Z>::one_,XX,A,Matrix<F,Z>::mone_,'N','T');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);
          
          if (sizeof(F)!=sizeof(Z)) {
            GF.solve(BB,XX,Factorization::CONJUGATE_TRANSPOSE,so);
            cout << "\nCONJUGATE_TRANSPOSE" << endl;
            cout << "XX = " << endl;
            XX.printOn(cout);
            R2.copy(BB);
            if (so==Factorization::LEFT_SIDE) {
              R2.gemm(Matrix<F,Z>::one_,A,XX,Matrix<F,Z>::mone_,'C','N');
            } else {
              R2.gemm(Matrix<F,Z>::one_,XX,A,Matrix<F,Z>::mone_,'N','C');
            }
            cout << "R2 = " << endl;
            R2.printOn(cout);

            XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
            GF.improve(BB,XX,berr,ferr,
              Factorization::CONJUGATE_TRANSPOSE,so);
            cout << "berr = " << endl;
            berr.printOn(cout);
            cout << "ferr = " << endl;
            ferr.printOn(cout);
            cout << "XX = " << endl;
            XX.printOn(cout);
            R2.copy(BB);
            if (so==Factorization::LEFT_SIDE) {
              R2.gemm(Matrix<F,Z>::one_,A,XX,Matrix<F,Z>::mone_,'C','N');
            } else {
              R2.gemm(Matrix<F,Z>::one_,XX,A,Matrix<F,Z>::mone_,'N','C');
            }
            cout << "R2 = " << endl;
            R2.printOn(cout);
          }
          if (so==Factorization::RIGHT_SIDE) break;
        }
        if (po==Factorization::PIVOT_ROWS_AND_COLUMNS) break;
      }
      if (eo==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) break;
    }
  }
*/
  TridiagonalMatrix<F,Z> T(3);
  T(0,0)=static_cast<F>(3.)*scalar;
    T(0,1)=static_cast<F>(-2.)*scalar;
  T(1,0)=  -scalar;
    T(1,1)=static_cast<F>(4.)*scalar;
    T(1,2)=static_cast<F>(-3.)*scalar;
  T(2,1)=static_cast<F>(-2.)*scalar;
    T(2,2)=static_cast<F>(5.)*scalar;
  T.printOn(cout);
  b.printOn(cout);
/*
  {
    Vector<F,Z> x(3);
    Vector<F,Z> p(3,1.e-10);
    Vector<F,Z> r(3);
    Factorization::EQUILIBRATE_OPTION eo=Factorization::NO_EQUILIBRATION;
    for (Factorization::PIVOT_OPTION po=Factorization::NO_PIVOTING;;po++) {
//  Factorization::PIVOT_OPTION po=Factorization::PIVOT_ROWS;
      GaussianFactorization<F,Z,TridiagonalMatrix<F,Z> > GF(T,po,eo);
      cout << "\npo = " << po << endl;
      GF.printOn(cout);
      cout << "reciprocal condition number = "
           << GF.reciprocalConditionNumber() << endl;
      F berr,ferr;

      GF.solve(b,x);
      cout << "b = " << endl;
      b.printOn(cout);
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      T.gtmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
      cout << "r = " << endl;
      r.printOn(cout);

      x+=p;
      GF.improve(b,x,berr,ferr);
      cout << "berr,ferr = " << berr << " " << ferr << endl;
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      T.gtmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
      cout << "r = " << endl;
      r.printOn(cout);

      GF.solve(b,x,Factorization::TRANSPOSE);
      cout << "\nTRANSPOSE" << endl;
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      T.gtmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'T');
      cout << "r = " << endl;
      r.printOn(cout);

      x+=p;
      GF.improve(b,x,berr,ferr,Factorization::TRANSPOSE);
      cout << "berr,ferr = " << berr << " " << ferr << endl;
      cout << "x = " << endl;
      x.printOn(cout);
      r.copy(b);
      T.gtmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'T');
      cout << "r = " << endl;
      r.printOn(cout);

      if (sizeof(F)!=sizeof(Z)) {
        GF.solve(b,x,Factorization::CONJUGATE_TRANSPOSE);
        cout << "\nCONJUGATE_TRANSPOSE" << endl;
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        T.gtmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'C');
        cout << "r = " << endl;
        r.printOn(cout);

        x+=p;
        GF.improve(b,x,berr,ferr,Factorization::CONJUGATE_TRANSPOSE);
        cout << "berr,ferr = " << berr << " " << ferr << endl;
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        T.gtmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'C');
        cout << "r = " << endl;
        r.printOn(cout);
      }

      if (po!=Factorization::NO_PIVOTING) break;
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
    T.printOn(cout);
    BL.printOn(cout);
    BR.printOn(cout);

    Vector<F,F> berr(2),ferr(2);
    Factorization::EQUILIBRATE_OPTION eo=Factorization::NO_EQUILIBRATION;
    for (Factorization::PIVOT_OPTION po=Factorization::NO_PIVOTING;;po++)
    {
//    Factorization::PIVOT_OPTION po=Factorization::PIVOT_ROWS;
      GaussianFactorization<F,Z,TridiagonalMatrix<F,Z> > GF(T,po,eo);
      cout << "\npo,eo = " << po << " " << eo << endl;
      GF.printOn(cout);
      cout << "reciprocal condition number = "
           << GF.reciprocalConditionNumber() << endl;

      for (Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE;;
      so=Factorization::RIGHT_SIDE) {
//      Factorization::SIDE_OPTION so=Factorization::RIGHT_SIDE;
        Matrix<F,Z> &BB=(so==Factorization::LEFT_SIDE ? BL : BR);
        Matrix<F,Z> &XX=(so==Factorization::LEFT_SIDE ? XL : XR);
        Matrix<F,Z> &R2=(so==Factorization::LEFT_SIDE ? RL : RR);
        cout << "\nso = " << so << endl;

        GF.solve(BB,XX,Factorization::NO_TRANSPOSE,so);
        cout << "BB = " << endl;
        BB.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','N');
        } else {
          T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','N');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);

        XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
        GF.improve(BB,XX,berr,ferr,Factorization::NO_TRANSPOSE,so);
        cout << "berr = " << endl;
        berr.printOn(cout);
        cout << "ferr = " << endl;
        ferr.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','N');
        } else {
          T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','N');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);

        GF.solve(BB,XX,Factorization::TRANSPOSE,so);
        cout << "\nTRANSPOSE" << endl;
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','T');
        } else {
          T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','T');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);

        XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
        GF.improve(BB,XX,berr,ferr,Factorization::TRANSPOSE,so);
        cout << "berr = " << endl;
        berr.printOn(cout);
        cout << "ferr = " << endl;
        ferr.printOn(cout);
        cout << "XX = " << endl;
        XX.printOn(cout);
        R2.copy(BB);
        if (so==Factorization::LEFT_SIDE) {
          T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','T');
        } else {
          T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','T');
        }
        cout << "R2 = " << endl;
        R2.printOn(cout);
        
        if (sizeof(F)!=sizeof(Z)) {
          GF.solve(BB,XX,Factorization::CONJUGATE_TRANSPOSE,so);
          cout << "\nCONJUGATE_TRANSPOSE" << endl;
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','C');
          } else {
            T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','C');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);

          XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
          GF.improve(BB,XX,berr,ferr,
            Factorization::CONJUGATE_TRANSPOSE,so);
          cout << "berr = " << endl;
          berr.printOn(cout);
          cout << "ferr = " << endl;
          ferr.printOn(cout);
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','C');
          } else {
            T.gtmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','C');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);
        }

        if (so==Factorization::RIGHT_SIDE) break;
      }
      if (po!=Factorization::NO_PIVOTING) break;
    }
  }
*/
  BandMatrix<F,Z> M(3,1,2);
  M(0,0)=static_cast<F>(3.)*scalar;
    M(0,1)=static_cast<F>(-2.)*scalar;
    M(0,2)=    scalar;
  M(1,0)=-scalar;
    M(1,1)=static_cast<F>(4.)*scalar;
    M(1,2)=static_cast<F>(-3.)*scalar;
  M(2,1)=static_cast<F>(-2.)*scalar;
    M(2,2)=static_cast<F>(5.)*scalar;
/*
  {
    Vector<F,Z> x(3);
    Vector<F,Z> r(3);
    for (Factorization::EQUILIBRATE_OPTION
    eo=Factorization::NO_EQUILIBRATION;;eo++) {
      for (Factorization::PIVOT_OPTION po=Factorization::NO_PIVOTING;;
      po=Factorization::PIVOT_ROWS) {
        GaussianFactorization<F,Z,BandMatrix<F,Z> > GF(M,po,eo);
        cout << "\npo,eo = " << po << " " << eo << endl;
        cout << "reciprocal condition number = "
             << GF.reciprocalConditionNumber() << endl;
        cout << "reciprocal pivot growth = " << GF.reciprocalPivotGrowth()
             << endl;

        F berr,ferr;

        GF.solve(b,x);
        cout << "b = " << endl;
        b.printOn(cout);
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        M.gbmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
        cout << "r = " << endl;
        r.printOn(cout);

        x+=p;
        GF.improve(b,x,berr,ferr);
        cout << "berr,ferr = " << berr << " " << ferr << endl;
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        M.gbmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r);
        cout << "r = " << endl;
        r.printOn(cout);

        GF.solve(b,x,Factorization::TRANSPOSE);
        cout << "\nTRANSPOSE" << endl;
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        M.gbmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'T');
        cout << "r = " << endl;
        r.printOn(cout);

        x+=p;
        GF.improve(b,x,berr,ferr,Factorization::TRANSPOSE);
        cout << "berr,ferr = " << berr << " " << ferr << endl;
        cout << "x = " << endl;
        x.printOn(cout);
        r.copy(b);
        M.gbmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'T');
        cout << "r = " << endl;
        r.printOn(cout);

        if (sizeof(F)!=sizeof(Z)) {
          GF.solve(b,x,Factorization::CONJUGATE_TRANSPOSE);
          cout << "\nCONJUGATE_TRANSPOSE" << endl;
          cout << "x = " << endl;
          x.printOn(cout);
          r.copy(b);
          M.gbmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'C');
          cout << "r = " << endl;
          r.printOn(cout);

          x+=p;
          GF.improve(b,x,berr,ferr,Factorization::CONJUGATE_TRANSPOSE);
          cout << "berr,ferr = " << berr << " " << ferr << endl;
          cout << "x = " << endl;
          x.printOn(cout);
          r.copy(b);
          M.gbmv(Matrix<F,Z>::one_,x,Matrix<F,Z>::mone_,r,'C');
          cout << "r = " << endl;
          r.printOn(cout);
        }

        if (po==Factorization::PIVOT_ROWS) break;
      }
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
    eo=Factorization::NO_EQUILIBRATION;;eo++) {
      for (Factorization::PIVOT_OPTION po=Factorization::NO_PIVOTING;;
      po=Factorization::PIVOT_ROWS) {
        GaussianFactorization<F,Z,BandMatrix<F,Z> > GF(M,po,eo);
        cout << "\npo,eo = " << po << " " << eo << endl;
        cout << "reciprocal condition number = "
             << GF.reciprocalConditionNumber() << endl;
        cout << "reciprocal pivot growth = " << GF.reciprocalPivotGrowth()
             << endl;

        for (Factorization::SIDE_OPTION so=Factorization::LEFT_SIDE;;
        so=Factorization::RIGHT_SIDE) {
          Matrix<F,Z> &BB=(so==Factorization::LEFT_SIDE ? BL : BR);
          Matrix<F,Z> &XX=(so==Factorization::LEFT_SIDE ? XL : XR);
          Matrix<F,Z> &R2=(so==Factorization::LEFT_SIDE ? RL : RR);
          cout << "\nso = " << so << endl;
          GF.solve(BB,XX,Factorization::NO_TRANSPOSE,so);
          cout << "BB = " << endl;
          BB.printOn(cout);
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','N');
          } else {
            M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','N');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);

          XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
          GF.improve(BB,XX,berr,ferr,Factorization::NO_TRANSPOSE,so);
          cout << "berr = " << endl;
          berr.printOn(cout);
          cout << "ferr = " << endl;
          ferr.printOn(cout);
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','N');
          } else {
            M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','N');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);

          GF.solve(BB,XX,Factorization::TRANSPOSE,so);
          cout << "\nTRANSPOSE" << endl;
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','T');
          } else {
            M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','T');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);

          XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
          GF.improve(BB,XX,berr,ferr,Factorization::TRANSPOSE,so);
          cout << "berr = " << endl;
          berr.printOn(cout);
          cout << "ferr = " << endl;
          ferr.printOn(cout);
          cout << "XX = " << endl;
          XX.printOn(cout);
          R2.copy(BB);
          if (so==Factorization::LEFT_SIDE) {
            M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','T');
          } else {
            M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','T');
          }
          cout << "R2 = " << endl;
          R2.printOn(cout);
          
          if (sizeof(F)!=sizeof(Z)) {
            GF.solve(BB,XX,Factorization::CONJUGATE_TRANSPOSE,so);
            cout << "\nCONJUGATE_TRANSPOSE" << endl;
            cout << "XX = " << endl;
            XX.printOn(cout);
            R2.copy(BB);
            if (so==Factorization::LEFT_SIDE) {
              M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','C');
            } else {
              M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','C');
            }
            cout << "R2 = " << endl;
            R2.printOn(cout);

            XX+=(so==Factorization::LEFT_SIDE ? PL : PR);
            GF.improve(BB,XX,berr,ferr,
              Factorization::CONJUGATE_TRANSPOSE,so);
            cout << "berr = " << endl;
            berr.printOn(cout);
            cout << "ferr = " << endl;
            ferr.printOn(cout);
            cout << "XX = " << endl;
            XX.printOn(cout);
            R2.copy(BB);
            if (so==Factorization::LEFT_SIDE) {
              M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'L','C');
            } else {
              M.gbmm(Matrix<F,Z>::one_,XX,Matrix<F,Z>::mone_,R2,'R','C');
            }
            cout << "R2 = " << endl;
            R2.printOn(cout);
          }
          if (so==Factorization::RIGHT_SIDE) break;
        }
        if (po==Factorization::PIVOT_ROWS) break;
      }
      if (eo==Factorization::EQUILIBRATE_ROWS_AND_COLUMNS) break;
    }
  }
//
}
     
// Modified from ltgmd.C by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
