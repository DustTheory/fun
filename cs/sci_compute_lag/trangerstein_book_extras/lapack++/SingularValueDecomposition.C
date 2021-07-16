#include "SingularValueDecomposition.H"

template<typename F,typename Z> void
SingularValueDecomposition<F,Z>::printOn(ostream &os) const {
  os << "SingularValueDecomposition: \n";
  os << "s': \n" << endl;
  s->printOn(cout);
  os << "U': \n" << endl;
  U->printOn(cout);
  os << "Vtranspose': \n" << endl;
  Vtranspose->printOn(cout);
  os << endl;
}

#include "SpecializedMatrix.H"
template<typename F,typename Z> void testSingularValueDecomposition(
F fscalar,Z scalar) {
  cout << "scalar = " << scalar << endl;
  LauchliMatrix<F,Z> L(3);
//cout << "L:" << endl;
//L.printOn(cout);

  SingularValueDecomposition<F,Z> SVD(L);
//cout << "SingularValueDecomposition(L)" << endl;
//SVD.printOn(cout);

  F rcond=0.1;
  F ridge=10.;
//
  {
//  x = [ 1/3, 1/3, 1/3 ]
    cout << "\nm>n" << endl;
    Vector<F,Z> b(4,scalar-scalar); b[0]=scalar/scalar;
    cout << "b:" << endl;
    b.printOn(cout);
    Vector<F,Z> x(3);
    SVD.solve(b,x,rcond);
    cout << "x(solve):" << endl;
    x.printOn(cout);
    SVD.regularize(b,x,ridge);
    cout << "x(regularize):" << endl;
    x.printOn(cout);
  }
//
//
  {
//  x = [ 1, e/3, e/3, e/3 ]
    cout << "\nm>n transpose" << endl;
    Vector<F,Z> b(3,scalar/scalar);
    cout << "b:" << endl;
    b.printOn(cout);
    Vector<F,Z> x(4);
    SVD.solve(b,x,rcond,Factorization::TRANSPOSE);
    cout << "x(solve):" << endl;
    x.printOn(cout);
    SVD.regularize(b,x,ridge,Factorization::TRANSPOSE);
    cout << "x(regularize):" << endl;
    x.printOn(cout);
  }
//

//
  {
//  X = [ 1/3  e/3 ]
//      [ 1/3  e/3 ]
//      [ 1/3  e/3 ]
    cout << "\nm>n and overdetermined" << endl;
    Matrix<F,Z> B(4,2,scalar-scalar);
    B(0,0)=scalar/scalar;
    B(1,1)=scalar/scalar;
    B(2,1)=scalar/scalar;
    B(3,1)=scalar/scalar;
    cout << "B:" << endl;
    B.printOn(cout);
    Matrix<F,Z> X(3,2);
    SVD.solve(B,X,rcond);
    cout << "X(solve):" << endl;
    X.printOn(cout);
    SVD.regularize(B,X,ridge);
    cout << "X(regularize):" << endl;
    X.printOn(cout);
  }
//
//
  {
//  X = [  1   1/3 ]
//      [ e/3  e/3 ]
//      [ e/3   0  ]
//      [ e/3   0  ]
    cout << "\nm>n and underdetermined" << endl;
    Matrix<F,Z> B(3,2,scalar/scalar);
    B(0,1)=scalar/scalar;
    B(1,1)=-scalar/scalar;
    B(2,1)=scalar/scalar;
    cout << "B:" << endl;
    B.printOn(cout);
    Matrix<F,Z> X(4,2);
    Vector<F,Z> residual_norm(2);
    SVD.solve(B,X,rcond,Factorization::TRANSPOSE);
    cout << "X(solve):" << endl;
    X.printOn(cout);
    SVD.regularize(B,X,ridge,Factorization::TRANSPOSE);
    cout << "X(regularize):" << endl;
    X.printOn(cout);
  }
//

  Matrix<F,Z> A(3,4);
  for (int j=0;j<4;j++) {
    for (int i=0;i<3;i++) A(i,j)=L(j,i);
  }
//cout << "A:" << endl;
//A.printOn(cout);

  SingularValueDecomposition<F,Z> SVD2(A);
//cout << "SingularValueDecomposition(A)" << endl;
//SVD2.printOn(cout);

//
  {
//  x = [ 1, e/3, e/3, e/3 ]
    cout << "\nm<n and underdetermined" << endl;
    Vector<F,Z> b(3,scalar/scalar);
    cout << "b:" << endl;
    b.printOn(cout);
    Vector<F,Z> x(4);
    SVD2.solve(b,x,rcond);
    cout << "x(solve):" << endl;
    x.printOn(cout);
    SVD2.regularize(b,x,ridge);
    cout << "x(regularize):" << endl;
    x.printOn(cout);
  }
//
//
  {
//  x = [ 1/3, 1/3, 1/3 ]
    cout << "\nm<n and overdetermined" << endl;
    Vector<F,Z> b(4,scalar-scalar);
    b[0]=scalar/scalar;
    cout << "b:" << endl;
    b.printOn(cout);
    Vector<F,Z> x(3);
    SVD2.solve(b,x,rcond,Factorization::TRANSPOSE);
    cout << "x(solve):" << endl;
    x.printOn(cout);
    SVD2.regularize(b,x,ridge,Factorization::TRANSPOSE);
    cout << "x(regularize):" << endl;
    x.printOn(cout);
  }
//

//
  {
//  X = [  1  1/3 ]
//      [ e/3 e/9 ]
//      [ e/3 e/9 ]
//      [ e/3 e/9 ]
    cout << "\nm<n and underdetermined" << endl;
    Matrix<F,Z> B(3,2,scalar/scalar);
    B(0,1)=scalar/scalar;
    B(1,1)=-scalar/scalar;
    B(2,1)=scalar/scalar;
    cout << "B:" << endl;
    B.printOn(cout);
    Matrix<F,Z> X(4,2);
    SVD2.solve(B,X,rcond);
    cout << "X(solve):" << endl;
    X.printOn(cout);
    SVD2.regularize(B,X,ridge);
    cout << "X(regularize):" << endl;
    X.printOn(cout);
  }
//
//
  {
//  X = [ 1/3 e/3 ]
//      [ 1/3 e/3 ]
//      [ 1/3 e/3 ]
    cout << "\nm<n and overdetermined" << endl;
    Matrix<F,Z> B(4,2,scalar-scalar);
    B(0,0)=scalar/scalar;
    B(1,1)=scalar/scalar;
    B(2,1)=scalar/scalar;
    B(3,1)=scalar/scalar;
    cout << "B:" << endl;
    B.printOn(cout);
    Matrix<F,Z> X(3,2);
    SVD2.solve(B,X,rcond,Factorization::TRANSPOSE);
    cout << "X(solve):" << endl;
    X.printOn(cout);
    SVD2.regularize(B,X,ridge,Factorization::TRANSPOSE);
    cout << "X(regularize):" << endl;
    X.printOn(cout);
  }
//

}
     
// Modified from ltgmd.C by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
