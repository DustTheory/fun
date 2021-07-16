#include "CompleteOrthogonalDecomposition.H"

template<typename F,typename Z> void
CompleteOrthogonalDecomposition<F,Z>::printOn(ostream &s) const {
  s << "CompleteOrthogonalDecomposition: \n";
  s << "rank = " << the_rank << endl;
  s << "URV': \n" << endl;
  URV->printOn(cout);
  OrthogonalMatrix<F,Z> *U=this->leftOrthogonalPart();
  s << "U: \n" << endl;
  U->printOn(cout);
  delete U;
  OrthogonalMatrix<F,Z> *V=this->rightOrthogonalPartTransposed();
  s << "V': \n" << endl;
  V->printOn(cout);
  delete V;
  UpperTriangularMatrix<F,Z> *R=this->upperTriangularPart();
  s << "R: \n" << endl;
  R->printOn(cout);
  delete R;
  s << "utau: \n" << endl;
  utau->printOn(cout);
  s << "vtau: \n" << endl;
  vtau->printOn(cout);
  s << "jpvt = ";
  for (int i=0;i<URV->size(1);i++) s << jpvt[i] << " ";
  s << endl;
}

#include "SpecializedMatrix.H"
template<typename F,typename Z> void testCompleteOrthogonalDecomposition(
F fscalar,Z scalar) {
  cout << "scalar = " << scalar << endl;
  Z zero=scalar-scalar;
  Z one=scalar/scalar;
  LauchliMatrix<F,Z> L(3);
//cout << "L:" << endl;
//L.printOn(cout);

  F rcond=10.;
//A=U*R*V'
//U=[ 1  e 0 0 ], R=[sqrt{3}]
//  [ e -1 0 0 ]
//  [ 0  0 1 0 ]
//  [ 0  0 0 1 ]
//V'=[ sqrt{1/3}     sqrt{1/3}        sqrt{1/3}    ]
//   [ sqrt{1/3}  (1-sqrt{1/3})/2 -(1+sqrt{1/3})/2 ]
//   [ sqrt{1/3} -(1+sqrt{1/3})/2  (1-sqrt{1/3})/2 ]
  CompleteOrthogonalDecomposition<F,Z> cod(L,rcond);
  cout << "CompleteOrthogonalDecomposition(L)" << endl;
  cod.printOn(cout);

//
  {
//  x = [ 1/3 ], s = [ 1 ] / 3, r = [ 0 ] ( -e / 3 )
//      [ 1/3 ]      [ e ]          [ 1 ]
//      [ 1/3 ]      [ 0 ]          [ 1 ]
//                   [ 0 ]          [ 1 ]
//  A' r = [ 1 ] ( - e^2 / 3 ), x - A' s = [ 1 ] ( - e^2 / 3 )
//         [ 1 ]                           [ 0 ]
//         [ 1 ]                           [ 0 ]
    cout << "\nm>n" << endl;
    Vector<F,Z> b(4,zero); b[0]=one;
    cout << "b:" << endl;
    b.printOn(cout);
    Vector<F,Z> x(3);
    F residual_norm=cod.solve(b,x);
    cout << "x:" << endl;
    x.printOn(cout);
    cout << "residual_norm = " << residual_norm << endl;
    Vector<F,Z> r(4);
    r.copy(b);
    L.gemv(-one,x,one,r,'N');
    cout << "\tr: " << endl;
    r.printOn(cout);
    Vector<F,Z> Atr(3,zero);
    L.gemv(one,r,zero,Atr,'C');
    cout << "\tA' * r: " << endl;
    Atr.printOn(cout);
    Vector<F,Z> s(4);
    F rn=cod.solve(x,s,Factorization::TRANSPOSE);
    cout << "s:" << endl;
    s.printOn(cout);
    Vector<F,Z> dx(3);
    dx.copy(x);
    L.gemv(-one,s,one,dx,'C');
    cout << "\tx-A's:" << endl;
    dx.printOn(cout);
  }
//
//
  {
//  x = [ 1 ], s = [ 1 ] / 3, r = [ 1 ] ( -e^2 )
//      [ e ]      [ 1 ]          [ 0 ]
//      [ 0 ]      [ 1 ]          [ 0 ]
//      [ 0 ]                    
//  A r = [ 1 ] ( - e^2 ), x - A s = [  0 ] ( e / 3 )
//        [ e ]                      [  2 ]
//        [ 0 ]                      [ -1 ]
//        [ 0 ]                      [ -1 ]
    cout << "\nm>n transpose" << endl;
    Vector<F,Z> b(3,one);
    cout << "b:" << endl;
    b.printOn(cout);
    Vector<F,Z> x(4);
    F residual_norm=cod.solve(b,x,Factorization::TRANSPOSE);
    cout << "x:" << endl;
    x.printOn(cout);
    cout << "residual_norm = " << residual_norm << endl;
    Vector<F,Z> r(3);
    r.copy(b);
    L.gemv(-one,x,one,r,'C');
    cout << "\tr: " << endl;
    r.printOn(cout);
    Vector<F,Z> Atr(4,zero);
    L.gemv(one,r,zero,Atr,'N');
    cout << "\tA * r: " << endl;
    Atr.printOn(cout);
    Vector<F,Z> s(3);
    F rn=cod.solve(x,s,Factorization::NO_TRANSPOSE);
    cout << "s:" << endl;
    s.printOn(cout);
    Vector<F,Z> dx(4);
    dx.copy(x);
    L.gemv(-one,s,one,dx,'N');
    cout << "\tx-A s:" << endl;
    dx.printOn(cout);
  }
//

//
  {
//  X = [ 1  e ] / 3 , S = [ 1  e  ] / 3, R = [   0   -e ]
//      [ 1  e ]           [ e e^2 ]          [ -e/3   1 ]
//      [ 1  e ]           [ 0  0  ]          [ -e/3   1 ]
//                         [ 0  0  ]          [ -e/3   1 ]
//  A' R = [ 1 0 ] ( - e^2 / 3 ), x - A' s = [ 1 0 ] ( - e^2 / 3 )
//         [ 1 0 ]                           [ 0 0 ]
//         [ 1 0 ]                           [ 0 0 ]
    cout << "\nm>n and overdetermined" << endl;
    Matrix<F,Z> B(4,2,zero);
    B(0,0)=one;
    B(1,1)=one;
    B(2,1)=one;
    B(3,1)=one;
    cout << "B:" << endl;
    B.printOn(cout);
    Matrix<F,Z> X(3,2);
    Vector<F,F> residual_norm(2);
    cod.solve(B,X,residual_norm);
    cout << "X:" << endl;
    X.printOn(cout);
    cout << "residual_norm:" << endl;
    residual_norm.printOn(cout);
    Matrix<F,Z> R(4,2);
    R.copy(B);
    R.gemm(-one,L,X,one,'N','N');
    cout << "\tR: " << endl;
    R.printOn(cout);
    Matrix<F,Z> AtR(3,2,zero);
    AtR.gemm(one,L,R,zero,'C','N');
    cout << "\tA' * R: " << endl;
    AtR.printOn(cout);
    Matrix<F,Z> S(4,2);
    Vector<F,F> Rn(2);
    cod.solve(X,S,Rn,Factorization::TRANSPOSE);
    cout << "S:" << endl;
    S.printOn(cout);
    Matrix<F,Z> dX(3,2);
    dX.copy(X);
    dX.gemm(-one,L,S,one,'C','N');
    cout << "\tX-A'S:" << endl;
    dX.printOn(cout);
  }
//
//
  {
//  X = [ 1  1/3 ], S = [ 1/3 1/9 ], R = [ -e^2  2/3 ]
//      [ e  e/3 ]      [ 1/3 1/9 ]      [  0   -4/3 ]
//      [ 0   0  ]      [ 1/3 1/9 ]      [  0    2/3 ]
//      [ 0   0  ]                        
//  A R = [ -e^2/3   0  ], x - A s = [   0    0  ]
//        [ -e^3   2e/3 ]            [ 2e/3 2e/9 ]
//        [    0  -4e/3 ]            [ -e/3 -e/9 ]
//        [    0   2e/3 ]            [ -e/3 -e/9 ]
    cout << "\nm>n and underdetermined" << endl;
    Matrix<F,Z> B(3,2,one);
    B(0,1)=one;
    B(1,1)=-one;
    B(2,1)=one;
    cout << "B:" << endl;
    B.printOn(cout);
    Matrix<F,Z> X(4,2);
    Vector<F,F> residual_norm(2);
    cod.solve(B,X,residual_norm,Factorization::TRANSPOSE);
    cout << "X:" << endl;
    X.printOn(cout);
    cout << "residual_norm:" << endl;
    residual_norm.printOn(cout);
    Matrix<F,Z> R(3,2);
    R.copy(B);
    R.gemm(-one,L,X,one,'C','N');
    cout << "\tR: " << endl;
    R.printOn(cout);
    Matrix<F,Z> AtR(4,2,zero);
    AtR.gemm(one,L,R,zero,'N','N');
    cout << "\tA * R: " << endl;
    AtR.printOn(cout);
    Matrix<F,Z> S(3,2);
    Vector<F,F> rn(2);
    cod.solve(X,S,rn,Factorization::NO_TRANSPOSE);
    cout << "S:" << endl;
    S.printOn(cout);
    Matrix<F,Z> dX(4,2);
    dX.copy(X);
    dX.gemm(-one,L,S,one,'N','N');
    cout << "\tX-A S:" << endl;
    dX.printOn(cout);
  }
//

  Matrix<F,Z> A(3,4);
  for (int j=0;j<4;j++) {
    for (int i=0;i<3;i++) A(i,j)=L(j,i);
  }
//cout << "A:" << endl;
//A.printOn(cout);

//A=U*R*V'
//U =[ sqrt{1/3}     sqrt{1/3}        sqrt{1/3}    ], R=[sqrt{3}]
//   [ sqrt{1/3}  (1-sqrt{1/3})/2 -(1+sqrt{1/3})/2 ]
//   [ sqrt{1/3} -(1+sqrt{1/3})/2  (1-sqrt{1/3})/2 ]
//V'=[  1   e/3  e/3  e/3 ]
//   [ e/3  1/3 -2/3 -2/3 ]
//   [ e/3 -2/3  1/3 -2/3 ]
//   [ e/3  1/3 -2/3  1/3 ]
  CompleteOrthogonalDecomposition<F,Z> cod2(A,rcond);
  cout << "CompleteOrthogonalDecomposition(A)" << endl;
  cod2.printOn(cout);

//
  {
//  x = [  1  ], s = [ 1 ] / 3, r = [ 1 ] ( -e^2/3 )
//      [ e/3 ]      [ 1 ]          [ 1 ]
//      [ e/3 ]      [ 1 ]          [ 1 ]
//      [ e/3 ]                    
//  A r = [ 3 ] ( - e^2/3 ), x - A s = [ 0 ]
//        [ e ]                        [ 0 ]
//        [ e ]                        [ 0 ]
//        [ e ]                        [ 0 ]
    cout << "\nm<n and underdetermined" << endl;
    Vector<F,Z> b(3,one);
    cout << "b:" << endl;
    b.printOn(cout);
    Vector<F,Z> x(4);
    F residual_norm=cod2.solve(b,x);
    cout << "x:" << endl;
    x.printOn(cout);
    cout << "residual_norm = " << residual_norm << endl;
    Vector<F,Z> r(3);
    r.copy(b);
    A.gemv(-one,x,one,r,'N');
    cout << "\tr: " << endl;
    r.printOn(cout);
    Vector<F,Z> Atr(4,zero);
    A.gemv(one,r,zero,Atr,'C');
    cout << "\tA' * r: " << endl;
    Atr.printOn(cout);
    Vector<F,Z> s(3);
    F rn=cod2.solve(x,s,Factorization::TRANSPOSE);
    cout << "s:" << endl;
    s.printOn(cout);
    Vector<F,Z> dx(4);
    dx.copy(x);
    A.gemv(-one,s,one,dx,'C');
    cout << "\tx-A's:" << endl;
    dx.printOn(cout);
  }
//
//
  {
//  x = [ 1 ]/3, s = [ 1/3 ], r = [ 0 ] ( -e/3 )
//      [ 1 ]        [ e/9 ]      [ 1 ]
//      [ 1 ]        [ e/9 ]      [ 1 ]
//                   [ e/9 ]      [ 1 ]
//  A r = [ 1 ] ( -e^2/3 ), x - A s = [ 1 ] ( -e^2/9 )
//        [ 1 ]                       [ 0 ]
//        [ 1 ]                       [ 0 ]
    cout << "\nm<n and overdetermined" << endl;
    Vector<F,Z> b(4,zero);
    b[0]=one;
    cout << "b:" << endl;
    b.printOn(cout);
    Vector<F,Z> x(3);
    F residual_norm=cod2.solve(b,x,Factorization::TRANSPOSE);
    cout << "x:" << endl;
    x.printOn(cout);
    cout << "residual_norm = " << residual_norm << endl;
    Vector<F,Z> r(4);
    r.copy(b);
    A.gemv(-one,x,one,r,'C');
    cout << "\tr: " << endl;
    r.printOn(cout);
    Vector<F,Z> Atr(3,zero);
    A.gemv(one,r,zero,Atr,'N');
    cout << "\tA * r: " << endl;
    Atr.printOn(cout);
    Vector<F,Z> s(4);
    F rn=cod2.solve(x,s,Factorization::NO_TRANSPOSE);
    cout << "s:" << endl;
    s.printOn(cout);
    Vector<F,Z> dx(3);
    dx.copy(x);
    A.gemv(-one,s,one,dx,'N');
    cout << "\tx-A s:" << endl;
    dx.printOn(cout);
  }
//

//
  {
//  X = [  1  1/3 ], S = [ 1/3 1/9 ], R = [ e^2/3  2/3 ]
//      [ e/3 e/9 ]      [ 1/3 1/9 ]      [ e^2/3 -4/3 ]
//      [ e/3 e/9 ]      [ 1/3 1/9 ]      [ e^2/3  2/3 ]
//      [ e/3 e/9 ]                    
//  A' R = [ e^2      0  ], x - A s = [ 0  0 ]
//         [ e^3/3  2e/3 ]            [ 0  0 ]
//         [ e^3/3 -4e/3 ]            [ 0  0 ]
//         [ e^3/3  2e/3 ]            [ 0  0 ]
    cout << "\nm<n and underdetermined" << endl;
    Matrix<F,Z> B(3,2,scalar/scalar);
    B(0,1)=scalar/scalar;
    B(1,1)=-scalar/scalar;
    B(2,1)=scalar/scalar;
    cout << "B:" << endl;
    B.printOn(cout);
    Matrix<F,Z> X(4,2);
    Vector<F,F> residual_norm(2);
    cod2.solve(B,X,residual_norm);
    cout << "X:" << endl;
    X.printOn(cout);
    cout << "residual_norm:" << endl;
    residual_norm.printOn(cout);
    Matrix<F,Z> R(3,2);
    R.copy(B);
    R.gemm(-one,A,X,one,'N','N');
    cout << "\tR: " << endl;
    R.printOn(cout);
    Matrix<F,Z> AtR(4,2,zero);
    AtR.gemm(one,A,R,zero,'C','N');
    cout << "\tA'* R: " << endl;
    AtR.printOn(cout);
    Matrix<F,Z> S(3,2);
    Vector<F,F> rn(2);
    cod2.solve(X,S,rn,Factorization::TRANSPOSE);
    cout << "S:" << endl;
    S.printOn(cout);
    Matrix<F,Z> dX(4,2);
    dX.copy(X);
    dX.gemm(-one,A,S,one,'C','N');
    cout << "\tX-A'S:" << endl;
    dX.printOn(cout);
  }
//
//
  {
//  X = [ 1 e ]/3, S = [ 1/3   e/3  ], R = [   0   -e ]
//      [ 1 e ]        [ e/9  e^2/3 ]      [ -e/3   1 ]
//      [ 1 e ]        [ e/9  e^2/3 ]      [ -e/3   1 ]
//                     [ e/9  e^2/3 ]      [ -e/3   1 ]
//  A R = [ -e^2/3 0 ], x - A s = [ -e^2/9 -e^3/9 ]
//        [ -e^2/3 0 ]            [ -e^2/9 -e^3/9 ]
//        [ -e^2/3 0 ]            [ -e^2/9 -e^3/9 ]
    cout << "\nm<n and overdetermined" << endl;
    Matrix<F,Z> B(4,2,scalar-scalar);
    B(0,0)=scalar/scalar;
    B(1,1)=scalar/scalar;
    B(2,1)=scalar/scalar;
    B(3,1)=scalar/scalar;
    cout << "B:" << endl;
    B.printOn(cout);
    Matrix<F,Z> X(3,2);
    Vector<F,F> residual_norm(2);
    cod2.solve(B,X,residual_norm,Factorization::TRANSPOSE);
    cout << "X:" << endl;
    X.printOn(cout);
    cout << "residual_norm:" << endl;
    residual_norm.printOn(cout);
    Matrix<F,Z> R(4,2);
    R.copy(B);
    R.gemm(-one,A,X,one,'C','N');
    cout << "\tR: " << endl;
    R.printOn(cout);
    Matrix<F,Z> AtR(3,2,zero);
    AtR.gemm(one,A,R,zero,'N','N');
    cout << "\tA * R: " << endl;
    AtR.printOn(cout);
    Matrix<F,Z> S(4,2);
    Vector<F,F> Rn(2);
    cod2.solve(X,S,Rn,Factorization::NO_TRANSPOSE);
    cout << "S:" << endl;
    S.printOn(cout);
    Matrix<F,Z> dX(3,2);
    dX.copy(X);
    dX.gemm(-one,A,S,one,'N','N');
    cout << "\tX-A S:" << endl;
    dX.printOn(cout);
  }
//

}
     
// Modified from ltgmd.C by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
