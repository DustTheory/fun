#include "HouseholderQRFactorization.H"
#include "Tracer.H"

template<typename F,typename Z> void
HouseholderQRFactorization<F,Z>::printOn(ostream &s) const {
  s << "HouseholderQRFactorization: \n";
  s << "QR: \n" << endl;
  QR->printOn(cout);
  OrthogonalMatrix<F,Z> *Q=orthogonalPart();
  s << "Q: \n" << endl;
  Q->printOn(cout);
  delete Q;
  TrapezoidalMatrix<F,Z> *R=rightTrapezoidalPart();
  s << "R: \n" << endl;
  R->printOn(cout);
  delete R;
  s << "tau: \n" << endl;
  tau->printOn(cout);
  s << "piv_op = " << piv_op << endl;
  if (piv_op!=Factorization::NO_PIVOTING) {
    s << "jpvt = ";
    for (int j=0;j<QR->size(1)-1;j++) s << jpvt[j] << " ";
    s << endl;
  }
  s << "iascl,ascl,anrm = " << iascl << " " << ascl << " " << anrm << endl;
}

#include "SpecializedMatrix.H"
template<typename F,typename Z> void testHouseholderQRFactorization(
F fscalar,Z scalar) {
  cout << "scalar = " << scalar << endl;
  LauchliMatrix<F,Z> L(3);
  cout << "L:" << endl;
  L.printOn(cout);

  { TRACER_CALL(tt,"HouseholderQRFactorization NO_PIVOTING");
    HouseholderQRFactorization<F,Z> HF(L);
//  cout << "HouseholderQRFactorization(L)" << endl;
//  HF.printOn(cout);

//
    { TRACER_CALL(ttt,"m>n and overdetermined");
//    cout << "\nm>n and overdetermined" << endl;
      Vector<F,Z> b(4,scalar-scalar); b[0]=scalar/scalar;
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> x(3);
      F residual_norm=HF.solve(b,x);
      cout << "x:" << endl;
      x.printOn(cout);
      cout << "residual_norm = " << residual_norm << endl;
      Vector<F,Z> r(4);
      r.copy(b);
      L.gemv(-scalar/scalar,x,scalar/scalar,r,'N');
      cout << "\tr: " << endl;
      r.printOn(cout);
      Vector<F,Z> Atr(3,scalar-scalar);
      L.gemv(scalar/scalar,r,scalar-scalar,Atr,'C');
      cout << "\tA' * r: " << endl;
      Atr.printOn(cout);

      Vector<F,Z> xx(3);
      xx.copy(x);
      HF.improve(b,xx);
      cout << "xx:" << endl;
      xx.printOn(cout);
      r.copy(b);
      L.gemv(-scalar/scalar,xx,scalar/scalar,r,'N');
      cout << "\tr: " << endl;
      r.printOn(cout);
      Atr=scalar-scalar;
      L.gemv(scalar/scalar,r,scalar-scalar,Atr,'C');
      cout << "\tA' * r: " << endl;
      Atr.printOn(cout);
    }
//
//
    { TRACER_CALL(ttt,"m>n and underdetermined");
//    cout << "\nm>n and underdetermined" << endl;
      Vector<F,Z> b(3,scalar/scalar);
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> x(4);
      F residual_norm=HF.solve(b,x,Factorization::TRANSPOSE);
      cout << "x:" << endl;
      x.printOn(cout);
      cout << "residual_norm = " << residual_norm << endl;
      Vector<F,Z> s(3);
      F rn=HF.solve(x,s,Factorization::NO_TRANSPOSE);
      cout << "s:" << endl;
      s.printOn(cout);
      cout << "rn = " << rn << endl;

      Vector<F,Z> xx(4);
      xx.copy(x);
      HF.improve(b,xx,Factorization::TRANSPOSE);
      cout << "xx:" << endl;
      xx.printOn(cout);
      Vector<F,Z> r(3);
      r.copy(b);
      L.gemv(-scalar/scalar,xx,scalar/scalar,r,'C');
      cout << "\tr: " << endl;
      r.printOn(cout);
    }
//

//
    { TRACER_CALL(ttt,"m>n and overdetermined");
//    cout << "\nm>n and overdetermined" << endl;
      Matrix<F,Z> B(4,2,scalar-scalar);
      B(0,0)=scalar/scalar;
      B(1,1)=scalar/scalar;
      B(2,1)=scalar/scalar;
      B(3,1)=scalar/scalar;
      cout << "B:" << endl;
      B.printOn(cout);
      Matrix<F,Z> X(3,2);
      Vector<F,Z> residual_norm(2);
      HF.solve(B,X,residual_norm);
      cout << "X:" << endl;
      X.printOn(cout);
      cout << "residual_norm:" << endl;
      residual_norm.printOn(cout);
      Matrix<F,Z> R(4,2);
      R.copy(B);
      R.gemm(-scalar/scalar,L,X,scalar/scalar,'N','N');
      cout << "\tR: " << endl;
      R.printOn(cout);
      Matrix<F,Z> AtR(3,2,scalar-scalar);
      AtR.gemm(scalar/scalar,L,R,scalar-scalar,'C','N');
      cout << "\tA' * R: " << endl;
      AtR.printOn(cout);

      Matrix<F,Z> XX(3,2);
      XX.copy(X);
      HF.improve(B,XX);
      cout << "XX:" << endl;
      XX.printOn(cout);
      R.copy(B);
      R.gemm(-scalar/scalar,L,XX,scalar/scalar,'N','N');
      cout << "\tR: " << endl;
      R.printOn(cout);
      AtR=scalar-scalar;
      AtR.gemm(scalar/scalar,L,R,scalar-scalar,'C','N');
      cout << "\tA' * R: " << endl;
      AtR.printOn(cout);
    }
//
//
    { TRACER_CALL(ttt,"m>n and underdetermined");
//    cout << "\nm>n and underdetermined" << endl;
      Matrix<F,Z> B(3,2,scalar/scalar);
      B(0,1)=scalar/scalar;
      B(1,1)=-scalar/scalar;
      B(2,1)=scalar/scalar;
      cout << "B:" << endl;
      B.printOn(cout);
      Matrix<F,Z> X(4,2);
      Vector<F,Z> residual_norm(2);
      HF.solve(B,X,residual_norm,Factorization::TRANSPOSE);
      cout << "X:" << endl;
      X.printOn(cout);
      cout << "residual_norm:" << endl;
      residual_norm.printOn(cout);
      Matrix<F,Z> S(3,2);
      Vector<F,Z> rn(2);
      HF.solve(X,S,rn,Factorization::NO_TRANSPOSE);
      cout << "S:" << endl;
      S.printOn(cout);
      cout << "rn: " << endl;
      rn.printOn(cout);

      Matrix<F,Z> XX(4,2);
      XX.copy(X);
      HF.improve(B,XX,Factorization::TRANSPOSE);
      cout << "XX:" << endl;
      XX.printOn(cout);
      Matrix<F,Z> R(3,2);
      R.copy(B);
      R.gemm(-scalar/scalar,L,XX,scalar/scalar,'C','N');
      cout << "\tR: " << endl;
      R.printOn(cout);
    }
//

    Matrix<F,Z> A(3,4);
    for (int j=0;j<4;j++) {
      for (int i=0;i<3;i++) A(i,j)=L(j,i);
    }
//  cout << "A:" << endl;
//  A.printOn(cout);

    HouseholderQRFactorization<F,Z> HF2(A);
//  cout << "HouseholderQRFactorization(A)" << endl;
//  HF2.printOn(cout);

/*
    {
      cout << "\nm<n and underdetermined" << endl;
      Vector<F,Z> b(3,scalar/scalar);
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> x(4);
      F residual_norm=HF2.solve(b,x);
      cout << "x:" << endl;
      x.printOn(cout);
      cout << "residual_norm = " << residual_norm << endl;

      Vector<F,Z> xx(4);
      xx.copy(x);
      HF2.improve(b,xx);
      cout << "xx:" << endl;
      xx.printOn(cout);
    }
*/
/*
    {
      cout << "\nm<n and overdetermined" << endl;
      Vector<F,Z> b(4,scalar-scalar);
      b[0]=scalar/scalar;
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> x(3);
      F residual_norm=HF2.solve(b,x,Factorization::TRANSPOSE);
      cout << "x:" << endl;
      x.printOn(cout);
      cout << "residual_norm = " << residual_norm << endl;

      Vector<F,Z> xx(3);
      xx.copy(x);
      HF2.improve(b,xx,Factorization::TRANSPOSE);
      cout << "xx:" << endl;
      xx.printOn(cout);
    }
*/

/*
    {
      cout << "\nm<n and underdetermined" << endl;
      Matrix<F,Z> B(3,2,scalar/scalar);
      B(0,1)=scalar/scalar;
      B(1,1)=-scalar/scalar;
      B(2,1)=scalar/scalar;
      cout << "B:" << endl;
      B.printOn(cout);
      Matrix<F,Z> X(4,2);
      Vector<F,Z> residual_norm(2);
      HF2.solve(B,X,residual_norm);
      cout << "X:" << endl;
      X.printOn(cout);
      cout << "residual_norm:" << endl;
      residual_norm.printOn(cout);

      Matrix<F,Z> XX(4,2);
      XX.copy(X);
      HF2.improve(B,XX);
      cout << "XX:" << endl;
      XX.printOn(cout);
    }
*/
/*
    {
      cout << "\nm<n and overdetermined" << endl;
      Matrix<F,Z> B(4,2,scalar-scalar);
      B(0,0)=scalar/scalar;
      B(1,1)=scalar/scalar;
      B(2,1)=scalar/scalar;
      B(3,1)=scalar/scalar;
      cout << "B:" << endl;
      B.printOn(cout);
      Matrix<F,Z> X(3,2);
      Vector<F,Z> residual_norm(2);
      HF2.solve(B,X,residual_norm,Factorization::TRANSPOSE);
      cout << "X:" << endl;
      X.printOn(cout);
      cout << "residual_norm:" << endl;
      residual_norm.printOn(cout);

      Matrix<F,Z> XX(3,2);
      XX.copy(X);
      HF2.improve(B,XX,Factorization::TRANSPOSE);
      cout << "XX:" << endl;
      XX.printOn(cout);
    }
*/
  }

  { TRACER_CALL(tt,"HouseholderQRFactorization COLUMN_PIVOTING");
    Matrix<F,Z> LL(4,3);
    LL.copy(L);
    Z one=scalar/scalar;
    Z two=one+one;
    Z four=two+two;
    for (int i=0;i<4;i++) LL(i,1)=two*L(i,1);
    for (int i=0;i<4;i++) LL(i,2)=four*L(i,2);
    HouseholderQRFactorization<F,Z> HF(LL,
      Factorization::PIVOT_COLUMNS);
//  HF.printOn(cout);

//
    { TRACER_CALL(ttt,"m>n and overdetermined");
      Vector<F,Z> b(4,scalar-scalar); b[0]=scalar/scalar;
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> x(3);
      F residual_norm=HF.solve(b,x);
      cout << "x:" << endl;
      x.printOn(cout);
      cout << "residual_norm = " << residual_norm << endl;
      Vector<F,Z> r(4);
      r.copy(b);
      LL.gemv(-scalar/scalar,x,scalar/scalar,r,'N');
      cout << "\tr: " << endl;
      r.printOn(cout);
      Vector<F,Z> Atr(3,scalar-scalar);
      LL.gemv(scalar/scalar,r,scalar-scalar,Atr,'C');
      cout << "\tA' * r: " << endl;
      Atr.printOn(cout);

//    Vector<F,Z> xx(3);
//    xx.copy(x);
//    HF.improve(b,xx);
//    cout << "xx:" << endl;
//    xx.printOn(cout);
//    r.copy(b);
//    LL.gemv(-scalar/scalar,x,scalar/scalar,r,'N');
//    cout << "\tr: " << endl;
//    r.printOn(cout);
//    Atr=scalar-scalar;
//    LL.gemv(scalar/scalar,r,scalar-scalar,Atr,'C');
//    cout << "\tA' * r: " << endl;
//    Atr.printOn(cout);
    }
//
//
    { TRACER_CALL(ttt,"m>n and underdetermined");
      Vector<F,Z> b(3,scalar/scalar);
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> x(4);
      F residual_norm=HF.solve(b,x,Factorization::TRANSPOSE);
      cout << "x:" << endl;
      x.printOn(cout);
      cout << "residual_norm = " << residual_norm << endl;
      Vector<F,Z> s(3);
      F rn=HF.solve(x,s,Factorization::NO_TRANSPOSE);
      cout << "s:" << endl;
      s.printOn(cout);
      cout << "rn = " << rn << endl;

//    Vector<F,Z> xx(4);
//    xx.copy(x);
//    HF.improve(b,xx,Factorization::TRANSPOSE);
//    cout << "xx:" << endl;
//    xx.printOn(cout);
//    Vector<F,Z> r(3);
//    r.copy(b);
//    LL.gemv(-scalar/scalar,x,scalar/scalar,r,'C');
//    cout << "\tr: " << endl;
//    r.printOn(cout);
    }
//

//
    { TRACER_CALL(ttt,"m>n and overdetermined");
      Matrix<F,Z> B(4,2,scalar-scalar);
      B(0,0)=scalar/scalar;
      B(1,1)=scalar/scalar;
      B(2,1)=scalar/scalar;
      B(3,1)=scalar/scalar;
      cout << "B:" << endl;
      B.printOn(cout);
      Matrix<F,Z> X(3,2);
      Vector<F,Z> residual_norm(2);
      HF.solve(B,X,residual_norm);
      cout << "X:" << endl;
      X.printOn(cout);
      cout << "residual_norm:" << endl;
      residual_norm.printOn(cout);
      Matrix<F,Z> R(4,2);
      R.copy(B);
      R.gemm(-scalar/scalar,LL,X,scalar/scalar,'N','N');
      cout << "\tR: " << endl;
      R.printOn(cout);
      Matrix<F,Z> AtR(3,2,scalar-scalar);
      AtR.gemm(scalar/scalar,LL,R,scalar-scalar,'C','N');
      cout << "\tA' * R: " << endl;
      AtR.printOn(cout);

//    Matrix<F,Z> XX(3,2);
//    XX.copy(X);
//    XX(0,0)=scalar/scalar;
//    HF.improve(B,XX);
//    cout << "XX:" << endl;
//    XX.printOn(cout);
//    R.copy(B);
//    R.gemm(-scalar/scalar,LL,XX,scalar/scalar,'N','N');
//    cout << "\tR: " << endl;
//    R.printOn(cout);
//    AtR=scalar-scalar;
//    AtR.gemm(scalar/scalar,LL,R,scalar-scalar,'C','N');
//    cout << "\tA' * R: " << endl;
//    AtR.printOn(cout);
    }
//
//
    { TRACER_CALL(ttt,"m>n and underdetermined");
      Matrix<F,Z> B(3,2,scalar/scalar);
      B(0,1)=scalar/scalar;
      B(1,1)=-scalar/scalar;
      B(2,1)=scalar/scalar;
      cout << "B:" << endl;
      B.printOn(cout);
      Matrix<F,Z> X(4,2);
      Vector<F,Z> residual_norm(2);
      HF.solve(B,X,residual_norm,Factorization::TRANSPOSE);
      cout << "X:" << endl;
      X.printOn(cout);
      cout << "residual_norm:" << endl;
      residual_norm.printOn(cout);
      Matrix<F,Z> S(3,2);
      Vector<F,Z> rn(2);
      HF.solve(X,S,rn,Factorization::NO_TRANSPOSE);
      cout << "S:" << endl;
      S.printOn(cout);
      cout << "rn: " << endl;
      rn.printOn(cout);

//    Matrix<F,Z> XX(4,2);
//    XX.copy(X);
//    HF.improve(B,XX,Factorization::TRANSPOSE);
//    cout << "XX:" << endl;
//    XX.printOn(cout);
//    Matrix<F,Z> R(3,2);
//    R.copy(B);
//    R.gemm(-scalar/scalar,LL,X,scalar/scalar,'C','N');
//    cout << "\tR: " << endl;
//    R.printOn(cout);
    }
//

    Matrix<F,Z> A(3,4);
    for (int j=0;j<4;j++) {
      for (int i=0;i<3;i++) A(i,j)=LL(j,i);
    }
//  cout << "A:" << endl;
//  A.printOn(cout);

    HouseholderQRFactorization<F,Z>
      HF2(A,Factorization::PIVOT_COLUMNS);
//  cout << "HouseholderQRFactorization(A)" << endl;
//  HF2.printOn(cout);

/*
    {
      cout << "\nm<n and underdetermined" << endl;
      Vector<F,Z> b(3,scalar/scalar);
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> x(4);
      F residual_norm=HF2.solve(b,x);
      cout << "x:" << endl;
      x.printOn(cout);
      cout << "residual_norm = " << residual_norm << endl;

      Vector<F,Z> xx(4);
      xx.copy(x);
      HF2.improve(b,xx);
      cout << "xx:" << endl;
      xx.printOn(cout);
    }
*/
/*
    {
      cout << "\nm<n and overdetermined" << endl;
      Vector<F,Z> b(4,scalar-scalar);
      b[0]=scalar/scalar;
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> x(3);
      F residual_norm=HF2.solve(b,x,Factorization::TRANSPOSE);
      cout << "x:" << endl;
      x.printOn(cout);
      cout << "residual_norm = " << residual_norm << endl;

      Vector<F,Z> xx(3);
      xx.copy(x);
      HF2.improve(b,xx,Factorization::TRANSPOSE);
      cout << "xx:" << endl;
      xx.printOn(cout);
    }
*/

/*
    {
      cout << "\nm<n and underdetermined" << endl;
      Matrix<F,Z> B(3,2,scalar/scalar);
      B(0,1)=scalar/scalar;
      B(1,1)=-scalar/scalar;
      B(2,1)=scalar/scalar;
      cout << "B:" << endl;
      B.printOn(cout);
      Matrix<F,Z> X(4,2);
      Vector<F,Z> residual_norm(2);
      HF2.solve(B,X,residual_norm);
      cout << "X:" << endl;
      X.printOn(cout);
      cout << "residual_norm:" << endl;
      residual_norm.printOn(cout);

      Matrix<F,Z> XX(4,2);
      XX.copy(X);
      HF2.improve(B,XX);
      cout << "XX:" << endl;
      XX.printOn(cout);
    }
*/
/*
    {
      cout << "\nm<n and overdetermined" << endl;
      Matrix<F,Z> B(4,2,scalar-scalar);
      B(0,0)=scalar/scalar;
      B(1,1)=scalar/scalar;
      B(2,1)=scalar/scalar;
      B(3,1)=scalar/scalar;
      cout << "B:" << endl;
      B.printOn(cout);
      Matrix<F,Z> X(3,2);
      Vector<F,Z> residual_norm(2);
      HF2.solve(B,X,residual_norm,Factorization::TRANSPOSE);
      cout << "X:" << endl;
      X.printOn(cout);
      cout << "residual_norm:" << endl;
      residual_norm.printOn(cout);

      Matrix<F,Z> XX(3,2);
      XX.copy(X);
      HF2.improve(B,XX,Factorization::TRANSPOSE);
      cout << "XX:" << endl;
      XX.printOn(cout);
    }
*/
  }

}
     
// Modified from ltgmd.C by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
