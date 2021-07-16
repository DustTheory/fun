#include "GramSchmidtQRFactorization.H"

template<typename F,typename Z> void
GramSchmidtQRFactorization<F,Z>::printOn(ostream &s) const {
  s << "GramSchmidtQRFactorization: \n";
  s << "Q: " << Q << "\n" << endl;
  Q->printOn(s);
  s << "R: " << R << "\n" << endl;
  R->printOn(s);
}

template<typename F,typename Z> void
testGramSchmidtQRFactorization(F fscalar,Z scalar) {
  cout << "scalar = " << scalar << endl;
  Z zero=scalar-scalar;
  Z one=scalar/scalar;
  Matrix<F,Z> *A=OPERATOR_NEW Matrix<F,Z>(3,2);
  (*A)(0,0)=scalar; (*A)(0,1)=scalar;
  (*A)(1,0)=scalar; (*A)(1,1)=zero;
  (*A)(2,0)=zero;   (*A)(2,1)=scalar;
  cout << "A:" << endl;
  A->printOn(cout);

  GramSchmidtQRFactorization<F,Z> *GSF=
    OPERATOR_NEW GramSchmidtQRFactorization<F,Z>(*A);
//cout << "GramSchmidtQRFactorization(A)" << endl;
//GSF->printOn(cout);

//
  const OrthogonalMatrix<F,Z> &Q=GSF->orthogonalPart();
  const UnitUpperTrapezoidalMatrix<F,Z> *R=GSF->rightTrapezoidalPart();
  {
    Matrix<F,Z> QtQ(2,2,zero);
    QtQ.gemm(one,Q,Q,zero,'C','N');
    cout << "Q^T * Q :" << endl;
    QtQ.printOn(cout);
  }
  {
    Matrix<F,Z> *QRminusA=R->trmm(Q,'R','N');
    (*QRminusA)-=(*A);
    cout << "Q * R - A :" << endl;
    QRminusA->printOn(cout);
    delete QRminusA; QRminusA=0;
  }
  delete R; R=0;
//

//
  {
  Vector<F,Z> b(3);
  b[0]=scalar;
  b[1]=zero;
  b[2]=zero;
  cout << "b:" << endl;
  b.printOn(cout);

  Vector<F,Z> x(2);
  Vector<F,Z> r(3);
  cout << "solveOverdetermined" << endl;
  GSF->solveOverdetermined(b,x,r);
  cout << "x:" << endl;
  x.printOn(cout);
  cout << "residual:" << endl;
  r.printOn(cout);

  Vector<F,Z> rerror(3);
  rerror.copy(r);
  rerror-=b;
  A->gemv(one,x,one,rerror,'N');
  cout << "r-b+A*x:" << endl;
  rerror.printOn(cout);

  Vector<F,Z> zerror(2,zero);
  A->gemv(one,r,zero,zerror,'C');
  cout << "A^H*residual:" << endl;
  zerror.printOn(cout);
  }
//

//
  {
  Matrix<F,Z> B(3,2);
  B(0,0)=scalar; B(0,1)=zero;
  B(1,0)=zero;   B(1,1)=scalar;
  B(2,0)=zero;   B(2,1)=scalar;
  cout << "B:" << endl;
  B.printOn(cout);

  Matrix<F,Z> X(2,2);
  Matrix<F,Z> R(3,2);
  cout << "solveOverdetermined" << endl;
  GSF->solveOverdetermined(B,X,R);
  cout << "X:" << endl;
  X.printOn(cout);
  cout << "residual:" << endl;
  R.printOn(cout);

  Matrix<F,Z> Rerror(3,2);
  Rerror.copy(R);
  Rerror-=B;
  Rerror.gemm(one,*A,X,one,'N','N');
  cout << "R-B+A*x:" << endl;
  Rerror.printOn(cout);

  Matrix<F,Z> Zerror(2,2,zero);
  Zerror.gemm(one,*A,R,zero,'C','N');
  cout << "A^H*residual:" << endl;
  Zerror.printOn(cout);
  }
//

//
  {
  Vector<F,Z> b(2);
  b[0]=scalar;
  b[1]=scalar;
  cout << "b:" << endl;
  b.printOn(cout);

  Vector<F,Z> x(3);
  cout << "solveUnderdetermined" << endl;
  GSF->solveUnderdetermined(b,x);
  cout << "x:" << endl;
  x.printOn(cout);

  Vector<F,Z> berror(2);
  berror.copy(b);
  A->gemv(-one,x,one,berror,'C');
  cout << "b-A*x:" << endl;
  berror.printOn(cout);

  Vector<F,Z> xerror(3);
  Vector<F,Z> s(2);
  GSF->solveOverdetermined(x,s,xerror);
  cout << "xerror:" << endl;
  xerror.printOn(cout);
  }
//

//
  {
  Matrix<F,Z> B(2,3);
  B(0,0)=scalar; B(0,1)=scalar; B(0,2)=zero;
  B(1,0)=scalar; B(1,1)=zero;   B(1,2)=scalar;
  cout << "B:" << endl;
  B.printOn(cout);

  Matrix<F,Z> X(3,3);
  cout << "solveUnderdetermined" << endl;
  GSF->solveUnderdetermined(B,X);
  cout << "X:" << endl;
  X.printOn(cout);

  Matrix<F,Z> Berror(2,3);
  Berror.copy(B);
  Berror.gemm(-one,*A,X,one,'C','N');
  cout << "B-A*X:" << endl;
  Berror.printOn(cout);

  Matrix<F,Z> Xerror(3,3);
  Matrix<F,Z> S(2,3);
  GSF->solveOverdetermined(X,S,Xerror);
  cout << "Xerror:" << endl;
  Xerror.printOn(cout);
  }
//

//
  {
  Vector<F,Z> b(3);
  b[0]=scalar;
  b[1]=zero;
  b[2]=zero;
  cout << "b:" << endl;
  b.printOn(cout);

  Vector<F,Z> x(2,zero);
  Vector<F,Z> r(3,zero);
//GSF->solveOverdetermined(b,x,r);
  cout << "improveOverdetermined" << endl;
  GSF->improveOverdetermined(b,x,r);
  cout << "x:" << endl;
  x.printOn(cout);
  cout << "residual:" << endl;
  r.printOn(cout);

  Vector<F,Z> rerror(3);
  rerror.copy(r);
  rerror-=b;
  A->gemv(one,x,one,rerror,'N');
  cout << "r-b+A*x:" << endl;
  rerror.printOn(cout);

  Vector<F,Z> zerror(2,zero);
  A->gemv(one,r,zero,zerror,'C');
  cout << "A^H*residual:" << endl;
  zerror.printOn(cout);
  }
//

//
  {
  Matrix<F,Z> B(3,2);
  B(0,0)=scalar; B(0,1)=zero;
  B(1,0)=zero;   B(1,1)=scalar;
  B(2,0)=zero;   B(2,1)=scalar;
  cout << "B:" << endl;
  B.printOn(cout);

  Matrix<F,Z> X(2,2,zero);
  Matrix<F,Z> R(3,2,zero);
//GSF->solveOverdetermined(B,X,R);
  cout << "improveOverdetermined" << endl;
  GSF->improveOverdetermined(B,X,R);
  cout << "X:" << endl;
  X.printOn(cout);
  cout << "residual:" << endl;
  R.printOn(cout);

  Matrix<F,Z> Rerror(3,2);
  Rerror.copy(R);
  Rerror-=B;
  Rerror.gemm(one,*A,X,one,'N','N');
  cout << "R-B+A*x:" << endl;
  Rerror.printOn(cout);

  Matrix<F,Z> Zerror(2,2,zero);
  Zerror.gemm(one,*A,R,zero,'C','N');
  cout << "A^H*residual:" << endl;
  Zerror.printOn(cout);
  }
//

//
  {
  Vector<F,Z> b(2);
  b[0]=scalar;
  b[1]=scalar;
  cout << "b:" << endl;
  b.printOn(cout);

  Vector<F,Z> x(3,zero);
  cout << "x:" << endl;
  x.printOn(cout);
  cout << "improveUnderdetermined" << endl;
//GSF->solveUnderdetermined(b,x);
  cout << "x:" << endl;
  x.printOn(cout);
  GSF->improveUnderdetermined(b,x);
  cout << "x:" << endl;
  x.printOn(cout);

  Vector<F,Z> berror(2);
  berror.copy(b);
  A->gemv(-one,x,one,berror,'C');
  cout << "b-A*x:" << endl;
  berror.printOn(cout);

  Vector<F,Z> xerror(3);
  Vector<F,Z> s(2);
  GSF->solveOverdetermined(x,s,xerror);
  cout << "xerror:" << endl;
  xerror.printOn(cout);
  }
//

//
  {
  Matrix<F,Z> B(2,3);
  B(0,0)=scalar; B(0,1)=scalar; B(0,2)=zero;
  B(1,0)=scalar; B(1,1)=zero;   B(1,2)=scalar;
  cout << "B:" << endl;
  B.printOn(cout);

  Matrix<F,Z> X(3,3,zero);
  cout << "improveUnderdetermined" << endl;
//GSF->solveUnderdetermined(B,X);
  GSF->improveUnderdetermined(B,X);
  cout << "X:" << endl;
  X.printOn(cout);

  Matrix<F,Z> Berror(2,3);
  Berror.copy(B);
  Berror.gemm(-one,*A,X,one,'C','N');
  cout << "B-A*X:" << endl;
  Berror.printOn(cout);

  Matrix<F,Z> Xerror(3,3);
  Matrix<F,Z> S(2,3);
  GSF->solveOverdetermined(X,S,Xerror);
  cout << "Xerror:" << endl;
  Xerror.printOn(cout);
  }
//

/*
  {
  GSF->dropColumn(1);
  cout << "\nafter dropColumn(1):" << endl;
//GSF->printOn(cout);
  const OrthogonalMatrix<F,Z> &Q=GSF->orthogonalPart();
  const UpperTrapezoidalMatrix<F,Z> &R=GSF->rightTrapezoidalPart();
  Matrix<F,Z> QtQ(1,1,zero);
  QtQ.gemm(one,Q,Q,zero,'C','N');
  cout << "Q^H * Q :" << endl;
  QtQ.printOn(cout);
  Matrix<F,Z> *QR=R.trmm(Q,'R','N');
  cout << "Q * R :" << endl;
  QR->printOn(cout);
  delete QR; QR=0;
  }
*/

/*
  {
  GSF->addColumn(1,*A);
  cout << "\nafter addColumn(1):" << endl;
//GSF->printOn(cout);
  const OrthogonalMatrix<F,Z> &Q=GSF->orthogonalPart();
  const UpperTrapezoidalMatrix<F,Z> &R=GSF->rightTrapezoidalPart();
  Matrix<F,Z> QtQ(2,2,zero);
  QtQ.gemm(one,Q,Q,zero,'C','N');
  cout << "Q^H * Q :" << endl;
  QtQ.printOn(cout);
  Matrix<F,Z> *QR=R.trmm(Q,'R','N');
  cout << "Q * R :" << endl;
  QR->printOn(cout);
  delete QR; QR=0;
  }
*/

/*
  {
  GSF->dropColumn(0);
  cout << "\nafter dropColumn(0):" << endl;
//GSF->printOn(cout);
  const OrthogonalMatrix<F,Z> &Q=GSF->orthogonalPart();
  const UpperTrapezoidalMatrix<F,Z> &R=GSF->rightTrapezoidalPart();
  Matrix<F,Z> QtQ(1,1,zero);
  QtQ.gemm(one,Q,Q,zero,'C','N');
  cout << "Q^H * Q :" << endl;
  QtQ.printOn(cout);
  Matrix<F,Z> *QR=R.trmm(Q,'R','N');
  cout << "Q * R :" << endl;
  QR->printOn(cout);
  delete QR; QR=0;
  }
*/

/*
  {
  GSF->addColumn(0,*A);
  cout << "\nafter addColumn(0):" << endl;
//GSF->printOn(cout);
  const OrthogonalMatrix<F,Z> &Q=GSF->orthogonalPart();
  const UpperTrapezoidalMatrix<F,Z> &R=GSF->rightTrapezoidalPart();
  Matrix<F,Z> QtQ(2,2,zero);
  QtQ.gemm(one,Q,Q,zero,'C','N');
  cout << "Q^H * Q :" << endl;
  QtQ.printOn(cout);
  Matrix<F,Z> *QR=R.trmm(Q,'R','N');
  cout << "Q * R :" << endl;
  QR->printOn(cout);
  delete QR; QR=0;
  }
*/

/*
  {
  GSF->dropRow(1);
  cout << "\nafter dropRow(1):" << endl;
//GSF->printOn(cout);
  const OrthogonalMatrix<F,Z> &Q=GSF->orthogonalPart();
  const UpperTrapezoidalMatrix<F,Z> &R=GSF->rightTrapezoidalPart();
  Matrix<F,Z> QtQ(2,2,zero);
  QtQ.gemm(one,Q,Q,zero,'C','N');
  cout << "Q^H * Q :" << endl;
  QtQ.printOn(cout);
  Matrix<F,Z> *QR=R.trmm(Q,'R','N');
  cout << "Q * R :" << endl;
  QR->printOn(cout);
  delete QR; QR=0;
  }
*/

/*
  {
  GSF->addRow(1,*A);
  cout << "\nafter addRow(1):" << endl;
//GSF->printOn(cout);
  const OrthogonalMatrix<F,Z> &Q=GSF->orthogonalPart();
  const UpperTrapezoidalMatrix<F,Z> &R=GSF->rightTrapezoidalPart();
  Matrix<F,Z> QtQ(2,2,zero);
  QtQ.gemm(one,Q,Q,zero,'C','N');
  cout << "Q^H * Q :" << endl;
  QtQ.printOn(cout);
  Matrix<F,Z> *QR=R.trmm(Q,'R','N');
  cout << "Q * R :" << endl;
  QR->printOn(cout);
  delete QR; QR=0;
  }
*/

  delete A;
  delete GSF;
}
     
// Modified from ltgmd.C by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
