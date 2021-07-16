#include "BandMatrix.H"
//#include "MyInline.H"

template<typename F,typename Z> BandMatrix<F,Z>*
TridiagonalMatrix<F,Z>::operator*(const TridiagonalMatrix<F,Z> &M) const {
  int n=size(0);
  CHECK_SAME(n,M.size(0));
  int nb=min(n-1,2);
  BandMatrix<F,Z> *B=OPERATOR_NEW BandMatrix<F,Z>(n,nb,nb);
  for (int j=2;j<n;j++) (*B)(j-2,j)=(*U)[j-2]*(*M.U)[j-1];
  for (int j=1;j<n;j++) {
    (*B)(j-1,j)=(*D)[j-1]*(*M.U)[j-1]+(*U)[j-1]*(*M.D)[j];
  }
  for (int j=0;j<n;j++) {
    Z &bjj=(*B)(j,j);
    bjj=(*D)[j]*(*M.D)[j];
    if (j>0) bjj+=(*L)[j-1]*(*M.U)[j-1];
    if (j<n-1) bjj+=(*U)[j]*(*M.L)[j];
  }
  for (int j=0;j<n-1;j++) {
    (*B)(j+1,j)=(*L)[j]*(*M.D)[j]+(*D)[j+1]*(*M.L)[j];
  }
  for (int j=0;j<n-2;j++) (*B)(j+2,j)=(*L)[j+1]*(*M.L)[j];
  return B;
}

template<typename F,typename Z> void TridiagonalMatrix<F,Z>::printOn(
ostream& s) const {
  s << "TridiagonalMatrix(" << size(0) << " x " << size(1) << ")\n" ;
  for (int i=0; i<dim; i++) {
    for (int j=0;j<dim; j++) s << operator()(i,j) << "  ";
    s << "\n";
  }
}

template<typename F,typename Z> void testTridiagonalMatrix(F fscalar,
Z scalar) {
  TridiagonalMatrix<F,Z> *T1=OPERATOR_NEW TridiagonalMatrix<F,Z>;
  cout << "\nafter TridiagonalMatrix()" << endl;
  T1->printOn(cout);
  delete T1; T1=0;

  int n=3;
  cout << "\nn = " << n << endl;
  TridiagonalMatrix<F,Z> *T2=OPERATOR_NEW TridiagonalMatrix<F,Z>(n);;
  cout << "after TridiagonalMatrix(n)" << endl;
  T2->printOn(cout);

  TridiagonalMatrix<F,Z> *T3=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  cout << "\nafter TridiagonalMatrix(n,scalar)" << endl;
  T3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=T3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }

  *T2=scalar;
  cout << "\nafter operator=" << endl;
  T2->printOn(cout);

  T3->resize(n+1);
  cout << "\nafter resize(n)" << endl;
  T3->printOn(cout);

  T3->resize(*T2);
  cout << "\nafter resize(T)" << endl;
  T3->printOn(cout);

  T3->copy(*T2);
  cout << "\nafter copy(T)" << endl;
  T3->printOn(cout);

  cout << "\nT3+T2:" << endl;
  TridiagonalMatrix<F,Z> *T=(*T3)+(*T2);
  T->printOn(cout);
  delete T; T=0;

  cout << "\nTridiagonalMatrix + SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *Sq=(*T2)+(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix + TridiagonalMatrix" << endl;
  Sq=(*Sy)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  cout << "\nTridiagonalMatrix + UnitUpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *U=OPERATOR_NEW 
    UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*U);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitUpperTrapezoidalMatrix + TridiagonalMatrix" << endl;
  Sq=(*U)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete U; U=0;

  cout << "\nTridiagonalMatrix + UpperTrapezoidalMatrix" << endl;
  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*U);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUpperTrapezoidalMatrix + TridiagonalMatrix" << endl;
  Sq=(*U)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete U; U=0;

  cout << "\nTridiagonalMatrix + UnitLowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *L=OPERATOR_NEW 
    UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitLowerTrapezoidalMatrix + TridiagonalMatrix" << endl;
  Sq=(*L)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nTridiagonalMatrix + LowerTrapezoidalMatrix" << endl;
  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTrapezoidalMatrix + TridiagonalMatrix" << endl;
  Sq=(*L)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nTridiagonalMatrix + SquareMatrix" << endl;
  SquareMatrix<F,Z> *SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*T2)+(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix + TridiagonalMatrix" << endl;
  Sq=(*SQ)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nTridiagonalMatrix + Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix + TridiagonalMatrix" << endl;
  Sq=(*M)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nT3-T2:" << endl;
  T=(*T3)-(*T2);
  T->printOn(cout);
  delete T; T=0;

  cout << "\nTridiagonalMatrix - SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sq=(*T2)-(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix - TridiagonalMatrix" << endl;
  Sq=(*Sy)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  cout << "\nTridiagonalMatrix - UnitUpperTrapezoidalMatrix" << endl;
  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*U);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitUpperTrapezoidalMatrix - TridiagonalMatrix" << endl;
  Sq=(*U)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete U; U=0;

  cout << "\nTridiagonalMatrix - UpperTrapezoidalMatrix" << endl;
  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*U);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUpperTrapezoidalMatrix - TridiagonalMatrix" << endl;
  Sq=(*U)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete U; U=0;

  cout << "\nTridiagonalMatrix - UnitLowerTrapezoidalMatrix" << endl;
  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitLowerTrapezoidalMatrix - TridiagonalMatrix" << endl;
  Sq=(*L)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nTridiagonalMatrix - LowerTrapezoidalMatrix" << endl;
  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTrapezoidalMatrix - TridiagonalMatrix" << endl;
  Sq=(*L)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nTridiagonalMatrix - SquareMatrix" << endl;
  SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*T2)-(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix - TridiagonalMatrix" << endl;
  Sq=(*SQ)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nTridiagonalMatrix - Matrix" << endl;
  M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix - TridiagonalMatrix" << endl;
  Sq=(*M)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nT3*scalar :" << endl;
  T=(*T3)*scalar;
  T->printOn(cout);
  delete T; T=0;

  cout << "\nT3/scalar.:" << endl;
  T=(*T3)/scalar;
  T->printOn(cout);
  delete T; T=0;

  cout << "\nTridiagonalMatrix * TridiagonalMatrix" << endl;
  BandMatrix<F,Z> *B=(*T2)*(*T3);
  B->printOn(cout);
  delete B; B=0;
  delete T3; T3=0;

  cout << "\nTridiagonalMatrix * SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sq=(*T2)*(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix * TridiagonalMatrix" << endl;
  Sq=(*Sy)*(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  L= OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nTridiagonalMatrix * UnitLowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,2,scalar);
  cout << "\nTridiagonalMatrix * UnitLowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,1,scalar);
  cout << "\nTridiagonalMatrix * UnitLowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(4,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(5,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L= OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,2,scalar);
  cout << "\nTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,1,scalar);
  cout << "\nTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(4,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(5,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nTridiagonalMatrix * UnitUpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,4,scalar);
  cout << "\nTridiagonalMatrix * UnitUpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,5,scalar);
  cout << "\nTridiagonalMatrix * UnitUpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(2,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(1,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,4,scalar);
  cout << "\nTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,5,scalar);
  cout << "\nTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(2,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(1,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * TridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  SquareMatrix<F,Z> *S=OPERATOR_NEW SquareMatrix<F,Z>(3,scalar);
  cout << "\nTridiagonalMatrix * SquareMatrix" << endl;
  M=(*T2)*(*S);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nSquareMatrix * TridiagonalMatrix" << endl;
  M=(*S)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete S; S=0;

  Matrix<F,Z> *M2=OPERATOR_NEW Matrix<F,Z>(3,2,scalar);
  cout << "\nTridiagonalMatrix * Matrix" << endl;
  M=(*T2)*(*M2);
  M->printOn(cout);
  delete M; M=0;
  delete M2; M2=0;

  M2=OPERATOR_NEW Matrix<F,Z>(4,3,scalar);
  cout << "\nMatrix * TridiagonalMatrix" << endl;
  M=(*M2)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete M2; M2=0;

  Vector<F,Z> *v=OPERATOR_NEW Vector<F,Z>(3,scalar);
  cout << "\nTridiagonalMatrix * Vector" << endl;
  Vector<F,Z> *w=(*T2)*(*v);
  w->printOn(cout);
  delete w; w=0;
  delete v; v=0;

  for (int i=0;i<n;i++) (*T2)(i,i)=static_cast<F>(2.)*scalar;
  for (int i=0;i<n-1;i++) (*T2)(i+1,i)=scalar;
  for (int i=0;i<n-1;i++) (*T2)(i,i+1)=static_cast<F>(3.)*scalar;
  cout <<"\nT2 = " << endl;
  T2->printOn(cout);
/*
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(n,scalar);
  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(n,scalar);
  cout << "\nlagtm" << endl;
  T2->lagtm(scalar,*x,scalar,*b);
  b->printOn(cout);

  *b=scalar;
  cout << "\nlagtm('T')" << endl;
  T2->lagtm(scalar,*x,scalar,*b,'T');
  b->printOn(cout);
  delete b; b=0;
  delete x; x=0;

  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(n,1,scalar);
  Matrix<F,Z> *B=OPERATOR_NEW Matrix<F,Z>(n,1,scalar);
  cout << "\nlagtm" << endl;
  T2->lagtm(scalar,*X,scalar,*B);
  B->printOn(cout);

  *B=scalar;
  cout << "\nlagtm('T')" << endl;
  T2->lagtm(scalar,*X,scalar,*B,'T');
  B->printOn(cout);
  delete B; B=0;
  delete X; X=0;
*/
  cout << "\nnormFrobenius = " << T2->normFrobenius() << endl;
  cout << "normInfinity = " << T2->normInfinity() << endl;
  cout << "normMaxEntry = " << T2->normMaxEntry() << endl;
  cout << "normOne = " << T2->normOne() << endl;
  cout << "reciprocalConditionNumber('I') = "
       << T2->reciprocalConditionNumber('I') << endl;
  cout << "reciprocalConditionNumber('O') = "
       << T2->reciprocalConditionNumber('O') << endl;
  delete T2; T2=0;

  TridiagonalMatrix<F,Z> *TT=OPERATOR_NEW TridiagonalMatrix<F,Z>(3,3);
  (*TT)(0,0)=static_cast<F>(2.)*scalar;
    (*TT)(0,1)=static_cast<F>(-2.)*scalar;
  (*TT)(1,0)=static_cast<F>(-1.)*scalar;
    (*TT)(1,1)=static_cast<F>(3.)*scalar;
    (*TT)(1,2)=static_cast<F>(-3.)*scalar;
  (*TT)(2,1)=static_cast<F>(-2.)*scalar;
    (*TT)(2,2)=static_cast<F>(4.)*scalar;

  cout << "\ngtmv('N')" << endl;
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(n,scalar);
  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(n,scalar);
  TT->gtmv(scalar,*x,scalar,*b,'N');
  b->printOn(cout);
  delete b; b=0;

  cout << "\ngtmv('T')" << endl;
  b=OPERATOR_NEW Vector<F,Z>(n,scalar);
  TT->gtmv(scalar,*x,scalar,*b,'T');
  b->printOn(cout);
  delete b; b=0;
  delete x; x=0;

  cout << "\ngtmm('L','N')" << endl;
  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(n,2,scalar);
  Matrix<F,Z> *B2=OPERATOR_NEW Matrix<F,Z>(n,2,scalar);
  TT->gtmm(scalar,*X,scalar,*B2,'L','N');
  B2->printOn(cout);
  delete B2; B2=0;

  cout << "\ngtmm('L','T')" << endl;
  B2=OPERATOR_NEW Matrix<F,Z>(n,2,scalar);
  TT->gtmm(scalar,*X,scalar,*B2,'L','T');
  B2->printOn(cout);
  delete B2; B2=0;
  delete X; X=0;

  cout << "\ngtmm('R','N')" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(2,n,scalar);
  B2=OPERATOR_NEW Matrix<F,Z>(2,n,scalar);
  TT->gtmm(scalar,*X,scalar,*B2,'R','N');
  B2->printOn(cout);
  delete B2; B2=0;

  cout << "\ngtmm('R','T')" << endl;
  B2=OPERATOR_NEW Matrix<F,Z>(2,n,scalar);
  TT->gtmm(scalar,*X,scalar,*B2,'R','T');
  B2->printOn(cout);
  delete B2; B2=0;
  delete X; X=0;

  b=OPERATOR_NEW Vector<F,Z>(3);
  (*b)[0]=static_cast<F>(0.)*scalar;
  (*b)[1]=static_cast<F>(6.)*scalar;
  (*b)[2]=static_cast<F>(0.)*scalar;
  cout << "\nsquare system of linear equations" << endl;
  TT->printOn(cout);
  b->printOn(cout);
  cout << "solve('N'):" << endl;
  x=OPERATOR_NEW Vector<F,Z>(3);
  TT->solve(*b,*x);
  x->printOn(cout);
  Vector<F,Z> *r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  TT->gtmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r);
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;

  cout << "\nsolve('T'):" << endl;
  TT->solve(*b,*x,'T');
  x->printOn(cout);
  r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  TT->gtmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'T');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;
  delete b; b=0;

  Matrix<F,Z> *BB=OPERATOR_NEW Matrix<F,Z>(3,1);
  (*BB)(0,0)=static_cast<F>(0.)*scalar;
  (*BB)(1,0)=static_cast<F>(6.)*scalar;
  (*BB)(2,0)=static_cast<F>(0.)*scalar;
  cout << "\nsolve('L','N'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(3,1);
  TT->solve(*BB,*X);
  X->printOn(cout);
  Matrix<F,Z> *R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*BB);
  TT->gtmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L','N');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  cout << "\nsolve('L','T'):" << endl;
  TT->solve(*BB,*X,'L','T');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*BB);
  TT->gtmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L','T');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete BB; BB=0;

  BB=OPERATOR_NEW Matrix<F,Z>(1,3);
  (*BB)(0,0)=static_cast<F>(0.)*scalar;
  (*BB)(0,1)=static_cast<F>(6.)*scalar;
  (*BB)(0,2)=static_cast<F>(0.)*scalar;
  cout << "\nsolve('R','N'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  TT->solve(*BB,*X,'R','N');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*BB);
  TT->gtmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R','N');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  cout << "\nsolve('R','T'):" << endl;
  TT->solve(*BB,*X,'R','T');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*BB);
  TT->gtmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R','T');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete BB; BB=0;
  delete TT; TT=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> BandMatrix<F,Z>*
SymmetricTridiagonalMatrix<F,Z>::operator*(
const SymmetricTridiagonalMatrix<F,Z> &S) const {
  int n=size(0);
  CHECK_SAME(n,S.size(0));
  int nb=min(n-1,2);
  BandMatrix<F,Z> *B=OPERATOR_NEW BandMatrix<F,Z>(n,nb,nb);
  for (int j=2;j<n;j++) (*B)(j-2,j)=
    upperDiagonalValue(j-2)*S.upperDiagonalValue(j-1);
  for (int j=1;j<n;j++) {
    (*B)(j-1,j)=diagonalValue(j-1)*S.upperDiagonalValue(j-1)
      +upperDiagonalValue(j-1)*S.diagonalValue(j);
  }
  for (int j=0;j<n;j++) {
    Z &bjj=(*B)(j,j);
    bjj=diagonalValue(j)*S.diagonalValue(j);
    if (j>0) bjj+=lowerDiagonalValue(j-1)*S.upperDiagonalValue(j-1);
    if (j<n-1) bjj+=upperDiagonalValue(j)*S.lowerDiagonalValue(j);
  }
  for (int j=0;j<n-1;j++) {
    (*B)(j+1,j)=lowerDiagonalValue(j)*S.diagonalValue(j)
      +diagonalValue(j+1)*S.lowerDiagonalValue(j);
  }
  for (int j=0;j<n-2;j++) {
    (*B)(j+2,j)=lowerDiagonalValue(j+1)*S.lowerDiagonalValue(j);
  }
  return B;
}

template<typename F,typename Z> BandMatrix<F,Z>*
SymmetricTridiagonalMatrix<F,Z>::operator*(
const TridiagonalMatrix<F,Z> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  int nb=min(n-1,2);
  BandMatrix<F,Z> *B=OPERATOR_NEW BandMatrix<F,Z>(n,nb,nb);
  for (int j=2;j<n;j++) (*B)(j-2,j)=upperDiagonalValue(j-2)*T(j-1,j);
  for (int j=1;j<n;j++) {
    (*B)(j-1,j)=diagonalValue(j-1)*T(j-1,j)
      +upperDiagonalValue(j-1)*T(j,j);
  }
  for (int j=0;j<n;j++) {
    Z &bjj=(*B)(j,j);
    bjj=diagonalValue(j)*T(j,j);
    if (j>0) bjj+=lowerDiagonalValue(j-1)*T(j-1,j);
    if (j<n-1) bjj+=upperDiagonalValue(j)*T(j+1,j);
  }
  for (int j=0;j<n-1;j++) {
    (*B)(j+1,j)=lowerDiagonalValue(j)*T(j,j)+diagonalValue(j+1)*T(j+1,j);
  }
  for (int j=0;j<n-2;j++) (*B)(j+2,j)=lowerDiagonalValue(j+1)*T(j+1,j);
  return B;
}

template<typename F,typename Z> BandMatrix<F,Z>* operator*(
const TridiagonalMatrix<F,Z> &T,
const SymmetricTridiagonalMatrix<F,Z> &St) {
  int n=St.size(0);
  CHECK_SAME(n,T.size(0));
  int nb=min(n-1,2);
  BandMatrix<F,Z> *B=OPERATOR_NEW BandMatrix<F,Z>(n,nb,nb);
  for (int j=2;j<n;j++) (*B)(j-2,j)=T(j-2,j-1)*St.upperDiagonalValue(j-1);
  for (int j=1;j<n;j++) {
    (*B)(j-1,j)=T(j-1,j-1)*St.upperDiagonalValue(j-1)
      +T(j-1,j)*St.diagonalValue(j);
  }
  for (int j=0;j<n;j++) {
    Z &bjj=(*B)(j,j);
    bjj=T(j,j)*St.diagonalValue(j);
    if (j>0) bjj+=T(j,j-1)*St.upperDiagonalValue(j-1);
    if (j<n-1) bjj+=T(j,j+1)*St.lowerDiagonalValue(j);
  }
  for (int j=0;j<n-1;j++) {
    (*B)(j+1,j)=T(j+1,j)*St.diagonalValue(j)
      +T(j+1,j+1)*St.lowerDiagonalValue(j);
  }
  for (int j=0;j<n-2;j++) (*B)(j+2,j)=T(j+2,j+1)*St.lowerDiagonalValue(j);
  return B;
}

template<typename F,typename Z> void
SymmetricTridiagonalMatrix<F,Z>::printOn(ostream& s) const {
  s << "SymmetricTridiagonalMatrix(" << size(0) << " x " << size(1)
    << ")\n" ;
  for (int i=0; i<dim; i++) {
    for (int j=0;j<dim;j++) s << (*this)(i,j) << " ";
    s << endl;
  }
}

template<typename F,typename Z> void
testSymmetricTridiagonalMatrix(F fscalar,Z scalar) {
  SymmetricTridiagonalMatrix<F,Z> *T1=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>;
  cout << "\nafter SymmetricTridiagonalMatrix()" << endl;
  T1->printOn(cout);
  delete T1; T1=0;

  int n=3;
  cout << "\nn = " << n << endl;
  SymmetricTridiagonalMatrix<F,Z> *T2=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n);;
  cout << "after SymmetricTridiagonalMatrix(n)" << endl;
  T2->printOn(cout);

  SymmetricTridiagonalMatrix<F,Z> *T3=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  cout << "\nafter SymmetricTridiagonalMatrix(n,scalar)" << endl;
  T3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=T3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }

  *T2=scalar;
  cout << "\nafter operator=" << endl;
  T2->printOn(cout);

  T3->resize(n+1);
  cout << "\nafter resize(n)" << endl;
  T3->printOn(cout);

  T3->resize(*T2);
  cout << "\nafter resize(T)" << endl;
  T3->printOn(cout);

  T3->copy(*T2);
  cout << "\nafter copy(T)" << endl;
  T3->printOn(cout);

  cout << "\nT3+T2:" << endl;
  SymmetricTridiagonalMatrix<F,Z> *T=(*T3)+(*T2);
  T->printOn(cout);
  delete T; T=0;

  cout << "\nSymmetricTridiagonalMatrix + TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  TridiagonalMatrix<F,Z> *Tr2=(*T2)+(*Tr);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nTridiagonalMatrix + SymmetricTridiagonalMatrix" << endl;
  Tr2=(*Tr)+(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Tr; Tr=0;

  cout << "\nSymmetricTridiagonalMatrix + SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *Sy2=(*T2)+(*Sy);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  cout << "\nSymmetricMatrix + SymmetricTridiagonalMatrix" << endl;
  Sy2=(*Sy)+(*T2);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;
  delete Sy; Sy=0;

  cout << "\nSymmetricTridiagonalMatrix + UnitUpperTrapezoidalMatrix"
       << endl;
  UpperTrapezoidalMatrix<F,Z> *U=OPERATOR_NEW 
    UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  UpperHessenbergMatrix<F,Z> *H=(*T2)+(*U);
  H->printOn(cout);
  delete H; H=0;

  cout << "\nUnitUpperTrapezoidalMatrix + SymmetricTridiagonalMatrix"
       << endl;
  H=(*U)+(*T2);
  H->printOn(cout);
  delete H; H=0;
  delete U; U=0;

  cout << "\nSymmetricTridiagonalMatrix + UpperTrapezoidalMatrix" << endl;
  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  H=(*T2)+(*U);
  H->printOn(cout);
  delete H; H=0;

  cout << "\nUpperTrapezoidalMatrix + SymmetricTridiagonalMatrix" << endl;
  H=(*U)+(*T2);
  H->printOn(cout);
  delete H; H=0;
  delete U; U=0;

  cout << "\nSymmetricTridiagonalMatrix + UnitLowerTrapezoidalMatrix"
       << endl;
  LowerTrapezoidalMatrix<F,Z> *L=OPERATOR_NEW 
    UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *Sq=(*T2)+(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitLowerTrapezoidalMatrix + SymmetricTridiagonalMatrix"
       << endl;
  Sq=(*L)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nSymmetricTridiagonalMatrix + LowerTrapezoidalMatrix" << endl;
  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTrapezoidalMatrix + SymmetricTridiagonalMatrix" << endl;
  Sq=(*L)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nSymmetricTridiagonalMatrix + SquareMatrix" << endl;
  SquareMatrix<F,Z> *SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*T2)+(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix + SymmetricTridiagonalMatrix" << endl;
  Sq=(*SQ)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nSymmetricTridiagonalMatrix + Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix + SymmetricTridiagonalMatrix" << endl;
  Sq=(*M)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nT3-T2:" << endl;
  T=(*T3)-(*T2);
  T->printOn(cout);
  delete T; T=0;

  cout << "\nSymmetricTridiagonalMatrix - TridiagonalMatrix" << endl;
  Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  Tr2=(*T2)-(*Tr);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nTridiagonalMatrix - SymmetricTridiagonalMatrix" << endl;
  Tr2=(*Tr)-(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Tr; Tr=0;

  cout << "\nSymmetricTridiagonalMatrix - SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sq=(*T2)-(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix - SymmetricTridiagonalMatrix" << endl;
  Sq=(*Sy)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  cout << "\nSymmetricTridiagonalMatrix - UnitUpperTrapezoidalMatrix"
       << endl;
  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*U);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitUpperTrapezoidalMatrix - SymmetricTridiagonalMatrix"
       << endl;
  Sq=(*U)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete U; U=0;

  cout << "\nSymmetricTridiagonalMatrix - UpperTrapezoidalMatrix" << endl;
  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*U);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUpperTrapezoidalMatrix - SymmetricTridiagonalMatrix" << endl;
  Sq=(*U)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete U; U=0;

  cout << "\nSymmetricTridiagonalMatrix - UnitLowerTrapezoidalMatrix"
       << endl;
  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitLowerTrapezoidalMatrix - SymmetricTridiagonalMatrix"
       << endl;
  Sq=(*L)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nSymmetricTridiagonalMatrix - LowerTrapezoidalMatrix" << endl;
  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTrapezoidalMatrix - SymmetricTridiagonalMatrix" << endl;
  Sq=(*L)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nSymmetricTridiagonalMatrix - SquareMatrix" << endl;
  SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*T2)-(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix - SymmetricTridiagonalMatrix" << endl;
  Sq=(*SQ)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nSymmetricTridiagonalMatrix - Matrix" << endl;
  M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix - SymmetricTridiagonalMatrix" << endl;
  Sq=(*M)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nT3*scalar :" << endl;
  Tr=(*T3)*scalar;
  Tr->printOn(cout);
  delete Tr; Tr=0;

  cout << "\nT3/scalar.:" << endl;
  Tr=(*T3)/scalar;
  Tr->printOn(cout);
  delete Tr; Tr=0;

  cout << "\nSymmetricTridiagonalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  BandMatrix<F,Z> *B=(*T2)*(*T3);
  B->printOn(cout);
  delete B; B=0;
  delete T3; T3=0;

  cout << "\nSymmetricTridiagonalMatrix * TridiagonalMatrix" << endl;
  Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  B=(*T2)*(*Tr);
  B->printOn(cout);
  delete B; B=0;


  cout << "\nTridiagonalMatrix * SymmetricTridiagonalMatrix" << endl;
  B=(*Tr)*(*T2);
  B->printOn(cout);
  delete B; B=0;
  delete Tr; Tr=0;

  cout << "\nSymmetricTridiagonalMatrix * SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sq=(*T2)*(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix * SymmetricTridiagonalMatrix" << endl;
  Sq=(*Sy)*(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  L= OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitLowerTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,2,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitLowerTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,1,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitLowerTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(4,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(5,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L= OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,2,scalar);
  cout << "\nSymmetricTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,1,scalar);
  cout << "\nSymmetricTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(4,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(5,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitUpperTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,4,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitUpperTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,5,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitUpperTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * SymmetricUnitTridiagonalMatrix"
       << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(2,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(1,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,4,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,5,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(2,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(1,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  SquareMatrix<F,Z> *S=OPERATOR_NEW SquareMatrix<F,Z>(3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * SquareMatrix" << endl;
  M=(*T2)*(*S);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nSquareMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*S)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete S; S=0;

  Matrix<F,Z> *M2=OPERATOR_NEW Matrix<F,Z>(3,2,scalar);
  cout << "\nSymmetricTridiagonalMatrix * Matrix" << endl;
  M=(*T2)*(*M2);
  M->printOn(cout);
  delete M; M=0;
  delete M2; M2=0;

  M2=OPERATOR_NEW Matrix<F,Z>(4,3,scalar);
  cout << "\nMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*M2)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete M2; M2=0;

  Vector<F,Z> *v=OPERATOR_NEW Vector<F,Z>(3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * Vector" << endl;
  Vector<F,Z> *w=(*T2)*(*v);
  w->printOn(cout);
  delete w; w=0;
  delete v; v=0;

  for (int i=0;i<n;i++) T2->diagonalValue(i)=2.*fscalar;
  for (int i=0;i<n-1;i++) T2->lowerDiagonalValue(i)=scalar;
  cout <<"\nT2 = " << endl;
  T2->printOn(cout);
/*
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(n,scalar);
  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(n,scalar);
  cout << "\nlagtm" << endl;
  T2->lagtm(scalar,*x,scalar,*b);
  b->printOn(cout);

  *b=scalar;
  cout << "\nlagtm('T')" << endl;
  T2->lagtm(scalar,*x,scalar,*b,'T');
  b->printOn(cout);
  delete b; b=0;
  delete x; x=0;

  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(n,1,scalar);
  Matrix<F,Z> *B=OPERATOR_NEW Matrix<F,Z>(n,1,scalar);
  cout << "\nlagtm" << endl;
  T2->lagtm(scalar,*X,scalar,*B);
  B->printOn(cout);

  *B=scalar;
  cout << "\nlagtm('T')" << endl;
  T2->lagtm(scalar,*X,scalar,*B,'T');
  B->printOn(cout);
  delete B; B=0;
  delete X; X=0;
*/
  cout << "\nnormFrobenius = " << T2->normFrobenius() << endl;
  cout << "normInfinity = " << T2->normInfinity() << endl;
  cout << "normMaxEntry = " << T2->normMaxEntry() << endl;
  cout << "normOne = " << T2->normOne() << endl;
  cout << "reciprocalConditionNumber('I') = "
       << T2->reciprocalConditionNumber('I') << endl;
  cout << "reciprocalConditionNumber('O') = "
       << T2->reciprocalConditionNumber('O') << endl;
  delete T2; T2=0;

  if (sizeof(Z)==sizeof(F)) {
    SymmetricTridiagonalMatrix<F,F> *TTT=
      OPERATOR_NEW SymmetricTridiagonalMatrix<F,F>(3,3);
    TTT->diagonalValue(0)=2.*fscalar;
    TTT->lowerDiagonalValue(0)=-1.*fscalar;
      TTT->diagonalValue(1)=3.*fscalar;
      TTT->lowerDiagonalValue(1)=-2.*fscalar;
        TTT->diagonalValue(2)=4.*fscalar;
    cout << "\neigenvalues" << endl;
    OrthogonalMatrix<F,F> *Q=OPERATOR_NEW OrthogonalMatrix<F,F>(3,3);
    Vector<F,F> *lambda=eigenvalues(*TTT,Q);
    lambda->printOn(cout);
    Q->printOn(cout);
    SquareMatrix<F,F> *RR=OPERATOR_NEW SquareMatrix<F,F>(3);
    for (int j=0;j<3;j++) {
      for (int i=0;i<3;i++) {
        (*RR)(i,j)=TTT->diagonalValue(i)*(*Q)(i,j)-(*Q)(i,j)*(*lambda)[j];
        if (i>0) (*RR)(i,j)+=TTT->lowerDiagonalValue(i-1)*(*Q)(i-1,j);
        if (i<2) (*RR)(i,j)+=TTT->upperDiagonalValue(i)*(*Q)(i+1,j);
      }
    }
    cout << "\nTTT * Q - Q * lambda" << endl;
    RR->printOn(cout);
    RR->gemm(Vector<F,F>::one_,*Q,*Q,Vector<F,F>::zero_,'C','N');
    cout << "\nQ^H * Q" << endl;
    RR->printOn(cout);
    delete RR; RR=0;
    delete lambda; lambda=0;
    delete Q; Q=0;
    delete TTT; TTT=0;
  }

  SymmetricTridiagonalMatrix<F,Z> *TT=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(3,3);
  TT->diagonalValue(0)=2.*fscalar;
  TT->lowerDiagonalValue(0)=static_cast<F>(-1.)*scalar;
    TT->diagonalValue(1)=static_cast<F>(3.)*fscalar;
    TT->lowerDiagonalValue(1)=static_cast<F>(-2.)*scalar;
      TT->diagonalValue(2)=static_cast<F>(4.)*fscalar;

  cout << "\nstmv" << endl;
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(n,scalar);
  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(n,scalar);
  TT->stmv(scalar,*x,scalar,*b);
  b->printOn(cout);
  delete b; b=0;
  delete x; x=0;

  cout << "\nstmm('L')" << endl;
  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(n,2,scalar);
  Matrix<F,Z> *B2=OPERATOR_NEW Matrix<F,Z>(n,2,scalar);
  TT->stmm(scalar,*X,scalar,*B2,'L');
  B2->printOn(cout);
  delete B2; B2=0;
  delete X; X=0;

  cout << "\nstmm('R')" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(2,n,scalar);
  B2=OPERATOR_NEW Matrix<F,Z>(2,n,scalar);
  TT->stmm(scalar,*X,scalar,*B2,'R');
  B2->printOn(cout);
  delete B2; B2=0;
  delete X; X=0;

  b=OPERATOR_NEW Vector<F,Z>(3);
  (*b)[0]=static_cast<F>(0.)*scalar;
  (*b)[1]=static_cast<F>(6.)*scalar;
  (*b)[2]=static_cast<F>(0.)*scalar;
  cout << "\nsquare system of linear equations" << endl;
  TT->printOn(cout);
  b->printOn(cout);
  cout << "solve:" << endl;
  x=OPERATOR_NEW Vector<F,Z>(3);
  TT->solve(*b,*x);
  x->printOn(cout);
  Vector<F,Z> *r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  TT->stmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r);
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;
  delete b; b=0;

  Matrix<F,Z> *BB=OPERATOR_NEW Matrix<F,Z>(3,1);
  (*BB)(0,0)=static_cast<F>(0.)*scalar;
  (*BB)(1,0)=static_cast<F>(6.)*scalar;
  (*BB)(2,0)=static_cast<F>(0.)*scalar;
  cout << "\nsolve('L'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(3,1);
  TT->solve(*BB,*X);
  X->printOn(cout);
  Matrix<F,Z> *R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*BB);
  TT->stmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete BB; BB=0;

  BB=OPERATOR_NEW Matrix<F,Z>(1,3);
  (*BB)(0,0)=static_cast<F>(0.)*scalar;
  (*BB)(0,1)=static_cast<F>(6.)*scalar;
  (*BB)(0,2)=static_cast<F>(0.)*scalar;
  cout << "\nsolve('R'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  TT->solve(*BB,*X,'R');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*BB);
  TT->stmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete BB; BB=0;
  delete TT; TT=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
SymmetricPositiveTridiagonalMatrix<F,Z>::printOn(ostream& s) const {
  int n=this->size(0);
  s << "SymmetricPositiveTridiagonalMatrix(" << n << " x " << n
    << ")\n" ;
  for (int i=0; i<n; i++) {
    for (int j=0;j<n;j++) s << (*this)(i,j) << " ";
    s << endl;
  }
}

template<typename F,typename Z> void
testSymmetricPositiveTridiagonalMatrix(F fscalar,Z scalar) {
  SymmetricPositiveTridiagonalMatrix<F,Z> *T1=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>;
  cout << "\nafter SymmetricPositiveTridiagonalMatrix()" << endl;
  T1->printOn(cout);
  delete T1; T1=0;

  int n=3;
  cout << "\nn = " << n << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *T2=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n);
  cout << "\nafter SymmetricPositiveTridiagonalMatrix(int)" << endl;
  T2->printOn(cout);

  SymmetricPositiveTridiagonalMatrix<F,Z> *T3=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  cout << "\nafter SymmetricPositiveTridiagonalMatrix(int,scalar)" << endl;
  T3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=T3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }

  *T2=scalar;
  cout << "\nafter T2=scalar" << endl;
  T2->printOn(cout);

  T3->resize(n+1);
  cout << "\nafter resize(n)" << endl;
  T3->printOn(cout);
  
  T3->resize(*T2);
  cout << "\nafter resize(T)" << endl;
  T3->printOn(cout);
  
  T3->copy(*T2);
  cout << "\nafter copy(T)" << endl;
  T3->printOn(cout);

  cout << "\nT3+T2:" << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *T=(*T3)+(*T2);
  T->printOn(cout);
  delete T; T=0;

  cout
    << "\nSymmetricPositiveTridiagonalMatrix + SymmetricTridiagonalMatrix"
    << endl;
  SymmetricTridiagonalMatrix<F,Z> *Str=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  SymmetricTridiagonalMatrix<F,Z> *Str2=(*T2)+(*Str);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  cout
    << "\nSymmetricTridiagonalMatrix + SymmetricPositiveTridiagonalMatrix"
    << endl;
  Str2=(*Str)+(*T2);
  Str2->printOn(cout);
  delete Str2; Str2=0;
  delete Str; Str=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + TridiagonalMatrix"
       << endl;
  TridiagonalMatrix<F,Z> *Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  TridiagonalMatrix<F,Z> *Tr2=(*T2)+(*Tr);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nTridiagonalMatrix + SymmetricPositiveTridiagonalMatrix"
       << endl;
  Tr2=(*Tr)+(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Tr; Tr=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *Sy2=(*T2)+(*Sy);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  cout << "\nSymmetricMatrix + SymmetricPositiveTridiagonalMatrix" << endl;
  Sy2=(*Sy)+(*T2);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;
  delete Sy; Sy=0;

  cout
    << "\nSymmetricPositiveTridiagonalMatrix + UnitUpperTrapezoidalMatrix"
    << endl;
  UpperTrapezoidalMatrix<F,Z> *U=OPERATOR_NEW 
    UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  UpperHessenbergMatrix<F,Z> *H=(*T2)+(*U);
  H->printOn(cout);
  delete H; H=0;

  cout
    << "\nUnitUpperTrapezoidalMatrix + SymmetricPositiveTridiagonalMatrix"
    << endl;
  H=(*U)+(*T2);
  H->printOn(cout);
  delete H; H=0;
  delete U; U=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + UpperTrapezoidalMatrix"
       << endl;
  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  H=(*T2)+(*U);
  H->printOn(cout);
  delete H; H=0;

  cout << "\nUpperTrapezoidalMatrix + SymmetricPositiveTridiagonalMatrix"
       << endl;
  H=(*U)+(*T2);
  H->printOn(cout);
  delete H; H=0;
  delete U; U=0;

  cout
    << "\nSymmetricPositiveTridiagonalMatrix + UnitLowerTrapezoidalMatrix"
    << endl;
  LowerTrapezoidalMatrix<F,Z> *L=OPERATOR_NEW 
    UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *Sq=(*T2)+(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout
    << "\nUnitLowerTrapezoidalMatrix + SymmetricPositiveTridiagonalMatrix"
    << endl;
  Sq=(*L)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + LowerTrapezoidalMatrix"
       << endl;
  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTrapezoidalMatrix + SymmetricPositiveTridiagonalMatrix"
       << endl;
  Sq=(*L)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + SquareMatrix" << endl;
  SquareMatrix<F,Z> *SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*T2)+(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix + SymmetricPositiveTridiagonalMatrix" << endl;
  Sq=(*SQ)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix + SymmetricPositiveTridiagonalMatrix" << endl;
  Sq=(*M)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nT3-T2:" << endl;
  Str=(*T3)+(*T2);
  Str->printOn(cout);
  delete Str; Str=0;

  cout
    << "\nSymmetricPositiveTridiagonalMatrix - SymmetricTridiagonalMatrix"
    << endl;
  Str=OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  Str2=(*T2)-(*Str);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  cout
    << "\nSymmetricTridiagonalMatrix - SymmetricPositiveTridiagonalMatrix"
    << endl;
  Str2=(*Str)-(*T2);
  Str2->printOn(cout);
  delete Str2; Str2=0;
  delete Str; Str=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - TridiagonalMatrix"
       << endl;
  Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  Tr2=(*T2)-(*Tr);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nTridiagonalMatrix - SymmetricPositiveTridiagonalMatrix"
       << endl;
  Tr2=(*Tr)-(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Tr; Tr=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sy2=(*T2)-(*Sy);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  cout << "\nSymmetricMatrix - SymmetricPositiveTridiagonalMatrix" << endl;
  Sy2=(*Sy)-(*T2);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;
  delete Sy; Sy=0;

  cout
    << "\nSymmetricPositiveTridiagonalMatrix - UnitUpperTrapezoidalMatrix"
    << endl;
  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  H=(*T2)-(*U);
  H->printOn(cout);
  delete H; H=0;

  cout
    << "\nUnitUpperTrapezoidalMatrix - SymmetricPositiveTridiagonalMatrix"
    << endl;
  H=(*U)-(*T2);
  H->printOn(cout);
  delete H; H=0;
  delete U; U=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - UpperTrapezoidalMatrix"
       << endl;
  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  H=(*T2)-(*U);
  H->printOn(cout);
  delete H; H=0;

  cout << "\nUpperTrapezoidalMatrix - SymmetricPositiveTridiagonalMatrix"
       << endl;
  H=(*U)-(*T2);
  H->printOn(cout);
  delete H; H=0;
  delete U; U=0;

  cout
    << "\nSymmetricPositiveTridiagonalMatrix - UnitLowerTrapezoidalMatrix"
    << endl;
  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout
    << "\nUnitLowerTrapezoidalMatrix - SymmetricPositiveTridiagonalMatrix"
    << endl;
  Sq=(*L)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - LowerTrapezoidalMatrix"
       << endl;
  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*L);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTrapezoidalMatrix - SymmetricPositiveTridiagonalMatrix"
       << endl;
  Sq=(*L)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete L; L=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - SquareMatrix" << endl;
  SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*T2)-(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix - SymmetricPositiveTridiagonalMatrix" << endl;
  Sq=(*SQ)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - Matrix" << endl;
  M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix - SymmetricPositiveTridiagonalMatrix" << endl;
  Sq=(*M)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nT3*scalar :" << endl;
  Tr=(*T3)*scalar;
  Tr->printOn(cout);
  delete Tr; Tr=0;
  
  cout << "\nT3/scalar.:" << endl;
  Tr=(*T3)/scalar;
  Tr->printOn(cout);
  delete Tr; Tr=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix * SymmetricPositiveTridiagonalMatrix" << endl;
  BandMatrix<F,Z> *B=(*T2)*(*T3);
  B->printOn(cout);
  delete B; B=0;
  delete T3; T3=0;

  cout
    << "\nSymmetricPositiveTridiagonalMatrix * SymmetricTridiagonalMatrix"
    << endl;
  Str=OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  B=(*T2)*(*Str);
  B->printOn(cout);
  delete B; B=0;

  cout
    << "\nSymmetricTridiagonalMatrix * SymmetricPositiveTridiagonalMatrix"
    << endl;
  B=(*Str)*(*T2);
  B->printOn(cout);
  delete B; B=0;
  delete Str; Str=0;

  cout << "\nSymmetricTridiagonalMatrix * TridiagonalMatrix" << endl;
  Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  B=(*T2)*(*Tr);
  B->printOn(cout);
  delete B; B=0;


  cout << "\nTridiagonalMatrix * SymmetricTridiagonalMatrix" << endl;
  B=(*Tr)*(*T2);
  B->printOn(cout);
  delete B; B=0;
  delete Tr; Tr=0;

  cout << "\nSymmetricTridiagonalMatrix * SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sq=(*T2)*(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix * SymmetricTridiagonalMatrix" << endl;
  Sq=(*Sy)*(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  L= OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitLowerTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,2,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitLowerTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,1,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitLowerTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(4,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(5,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L= OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,2,scalar);
  cout << "\nSymmetricTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,1,scalar);
  cout << "\nSymmetricTridiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  M=(*T2)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(4,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(5,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*L)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitUpperTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,4,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitUpperTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,5,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UnitUpperTrapezoidalMatrix"
       << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * SymmetricUnitTridiagonalMatrix"
       << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(2,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(1,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix"
       << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,4,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,5,scalar);
  cout << "\nSymmetricTridiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  M=(*T2)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(2,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(1,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*U)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  SquareMatrix<F,Z> *S=OPERATOR_NEW SquareMatrix<F,Z>(3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * SquareMatrix" << endl;
  M=(*T2)*(*S);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nSquareMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*S)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete S; S=0;

  Matrix<F,Z> *M2=OPERATOR_NEW Matrix<F,Z>(3,2,scalar);
  cout << "\nSymmetricTridiagonalMatrix * Matrix" << endl;
  M=(*T2)*(*M2);
  M->printOn(cout);
  delete M; M=0;
  delete M2; M2=0;

  M2=OPERATOR_NEW Matrix<F,Z>(4,3,scalar);
  cout << "\nMatrix * SymmetricTridiagonalMatrix" << endl;
  M=(*M2)*(*T2);
  M->printOn(cout);
  delete M; M=0;
  delete M2; M2=0;

  Vector<F,Z> *v=OPERATOR_NEW Vector<F,Z>(3,scalar);
  cout << "\nSymmetricTridiagonalMatrix * Vector" << endl;
  Vector<F,Z> *w=(*T2)*(*v);
  w->printOn(cout);
  delete w; w=0;
  delete v; v=0;

  for (int i=0;i<n;i++) T2->diagonalValue(i)=2.*fscalar;
  for (int i=0;i<n-1;i++) T2->lowerDiagonalValue(i)=scalar;
  cout <<"\nT2 = " << endl;
  T2->printOn(cout);
  cout << "\nnormFrobenius = " << T2->normFrobenius() << endl;
  cout << "normInfinity = " << T2->normInfinity() << endl;
  cout << "normMaxEntry = " << T2->normMaxEntry() << endl;
  cout << "normOne = " << T2->normOne() << endl;
  cout << "reciprocalConditionNumber('I') = "
       << T2->reciprocalConditionNumber('I') << endl;
  cout << "reciprocalConditionNumber('O') = "
       << T2->reciprocalConditionNumber('O') << endl;
  delete T2; T2=0;

  if (sizeof(Z)==sizeof(F)) {
    SymmetricPositiveTridiagonalMatrix<F,F> *TTT=
      OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,F>(3,3);
    TTT->diagonalValue(0)= 2.*fscalar;
    TTT->lowerDiagonalValue(0)=-1.*fscalar;
      TTT->diagonalValue(1)= 3.*fscalar;
      TTT->lowerDiagonalValue(1)=-2.*fscalar;
        TTT->diagonalValue(2)=4.*fscalar;
    OrthogonalMatrix<F,F> *Q=OPERATOR_NEW OrthogonalMatrix<F,F>(3,3);
    Vector<F,F> *lambda=eigenvalues(*TTT,Q);
    cout << "\nafter eigenvalues" << endl;
    lambda->printOn(cout);
    Q->printOn(cout);
    SquareMatrix<F,F> *RR=OPERATOR_NEW SquareMatrix<F,F>(3);
    for (int j=0;j<3;j++) {
      for (int i=0;i<3;i++) {
        (*RR)(i,j)=TTT->diagonalValue(i)*(*Q)(i,j)-(*Q)(i,j)*(*lambda)[j];
        if (i>0) (*RR)(i,j)+=TTT->lowerDiagonalValue(i-1)*(*Q)(i-1,j);
        if (i<2) (*RR)(i,j)+=TTT->upperDiagonalValue(i)*(*Q)(i+1,j);
      }
    }
    cout << "\nTTT * Q - Q * lambda" << endl;
    RR->printOn(cout);
    delete RR; RR=0;
    delete Q; Q=0;
    delete lambda; lambda=0;
    delete TTT; TTT=0;
  }

  SymmetricPositiveTridiagonalMatrix<F,Z> *TT=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(3,3);
  TT->diagonalValue(0)= static_cast<F>(2.)*fscalar;
  TT->lowerDiagonalValue(0)=static_cast<F>(-1.)*scalar;
    TT->diagonalValue(1)=static_cast<F>(3.)*fscalar;
    TT->lowerDiagonalValue(1)=static_cast<F>(-2.)*scalar;
      TT->diagonalValue(2)=static_cast<F>(4.)*fscalar;
  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(3);
  (*b)[0]=static_cast<F>(0.)*scalar;
  (*b)[1]=static_cast<F>(6.)*scalar;
  (*b)[2]=static_cast<F>(0.)*scalar;
  cout << "\nsquare system of linear equations" << endl;
  TT->printOn(cout);
  b->printOn(cout);
  cout << "solve:" << endl;
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(3);
  TT->solve(*b,*x);
  x->printOn(cout);
  Vector<F,Z> *r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  TT->stmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r);
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;
  delete b; b=0;
 
  Matrix<F,Z> *BB=OPERATOR_NEW Matrix<F,Z>(3,1);
  (*BB)(0,0)=static_cast<F>(0.)*scalar;
  (*BB)(1,0)=static_cast<F>(6.)*scalar;
  (*BB)(2,0)=static_cast<F>(0.)*scalar;
  cout << "\nsolve('L'):" << endl;
  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(3,1);
  TT->solve(*BB,*X);
  X->printOn(cout);
  Matrix<F,Z> *R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*BB);
  TT->stmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete BB; BB=0;

  BB=OPERATOR_NEW Matrix<F,Z>(1,3);
  (*BB)(0,0)=static_cast<F>(0.)*scalar;
  (*BB)(0,1)=static_cast<F>(6.)*scalar;
  (*BB)(0,2)=static_cast<F>(0.)*scalar;
  cout << "\nsolve('R'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  TT->solve(*BB,*X,'R');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*BB);
  TT->stmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete BB; BB=0;
  delete TT; TT=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> DiagonalMatrix<F,Z>*
DiagonalMatrix<F,Z>::operator*(
const DiagonalMatrix<F,Z> &M) const {
  CHECK_SAME(dim,M.dim);
  DiagonalMatrix<F,Z> *P=OPERATOR_NEW DiagonalMatrix<F,Z>(dim);
  for (int j=0;j<dim;j++) {
    P->D->operator[](j)=D->operator[](j)*M.D->operator[](j);
  }
  return P;
}

template<typename F,typename Z> TridiagonalMatrix<F,Z>* 
DiagonalMatrix<F,Z>::operator*(const TridiagonalMatrix<F,Z> &T) const {
  int n=size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<F,Z> *S=OPERATOR_NEW TridiagonalMatrix<F,Z>(n);
  S->copy(T);
  for (int i=0;i<n;i++) {
    Z di=D->operator[](i);
    if (i>0) (*S)(i,i-1)*=di;
    (*S)(i,i)*=di;
    if (i<n-1) (*S)(i,i+1)*=di;
  }
  return S;
}

template<typename F,typename Z> TridiagonalMatrix<F,Z>* operator*(
const TridiagonalMatrix<F,Z> &T,const DiagonalMatrix<F,Z> &A) {
  int n=A.size(0);
  CHECK_SAME(n,T.size(0));
  TridiagonalMatrix<F,Z> *S=OPERATOR_NEW TridiagonalMatrix<F,Z>(n);
  S->copy(T);
  for (int j=0;j<n;j++) {
    Z dj=A[j];
    if (j>0) (*S)(j-1,j)*=dj;
    (*S)(j,j)*=dj;
    if (j<n-1) (*S)(j+1,j)*=dj;
  }
  return S;
}

template<typename F,typename Z> UpperTrapezoidalMatrix<F,Z>*
DiagonalMatrix<F,Z>::operator*(const UpperTrapezoidalMatrix<F,Z> &U)
const {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(m,size(0));
  UpperTrapezoidalMatrix<F,Z> *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(m,n);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<F,Z>*>(&U)==0);
  for (int i=0;i<m;i++) {
    Z di=D->operator[](i);
    if (U_non_unit) {
      const Z *Uij=U.addr(i,i);
      Z *Sij=S->addr(i,i);
      for (int j=i;j<n;j++,Uij+=m,Sij+=m) *Sij=di*(*Uij);
    } else {
      (*S)(i,i)=di;
      if (i+1<n) {
        const Z *Uij=U.addr(i,i+1);
        Z *Sij=S->addr(i,i+1);
        for (int j=i+1;j<n;j++,Uij+=m,Sij+=m) *Sij=di*(*Uij);
      }
    }
  }
  return S;
}

template<typename F,typename Z> UpperTrapezoidalMatrix<F,Z>* operator*(
const UpperTrapezoidalMatrix<F,Z> &U,const DiagonalMatrix<F,Z> &A) {
  int m=U.size(0),n=U.size(1);
  CHECK_SAME(n,A.size(0));
  UpperTrapezoidalMatrix<F,Z> *S=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(m,n);
  bool U_non_unit=(dynamic_cast<const
    UnitUpperTrapezoidalMatrix<F,Z>*>(&U)==0);
  for (int j=0;j<n;j++) {
    Z dj=A[j];
    const Z *Uij=U.addr(0,j);
    Z *Sij=S->addr(0,j);
    if (U_non_unit) {
      for (int i=0;i<=min(m-1,j);i++,Uij++,Sij++) *Sij=(*Uij)*dj;
    } else {
      for (int i=0;i<min(m,j);i++,Uij++,Sij++) *Sij=(*Uij)*dj;
      if (j<m) (*S)(j,j)=dj;
    }
  }
  return S;
}

template<typename F,typename Z> LowerTrapezoidalMatrix<F,Z>*
DiagonalMatrix<F,Z>::operator*(const LowerTrapezoidalMatrix<F,Z> &L)
const {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(m,size(0));
  LowerTrapezoidalMatrix<F,Z> *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(m,n);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<F,Z>*>(&L)==0);
  for (int i=0;i<m;i++) {
    Z di=D->operator[](i);
    const Z *Lij=L.addr(i,0);
    Z *Sij=S->addr(i,0);
    if (L_non_unit) {
      for (int j=0;j<=min(i,n-1);j++,Lij+=m,Sij+=m) *Sij=di*(*Lij);
    } else {
      for (int j=0;j<min(i,n);j++,Lij+=m,Sij+=m) *Sij=di*(*Lij);
      if (i<n) (*S)(i,i)=di;
    }
  }
  return S;
}

template<typename F,typename Z> LowerTrapezoidalMatrix<F,Z>* operator*(
const LowerTrapezoidalMatrix<F,Z> &L,const DiagonalMatrix<F,Z> &A) {
  int m=L.size(0),n=L.size(1);
  CHECK_SAME(n,A.size(0));
  LowerTrapezoidalMatrix<F,Z> *S=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(m,n);
  bool L_non_unit=(dynamic_cast<const
    UnitLowerTrapezoidalMatrix<F,Z>*>(&L)==0);
  for (int j=0;j<n;j++) {
    Z dj=A[j];
    if (L_non_unit) {
      const Z *Lij=L.addr(j,j);
      Z *Sij=S->addr(j,j);
      for (int i=j;i<m;i++,Lij++,Sij++) *Sij=(*Lij)*dj;
    } else {
      (*S)(j,j)=dj;
      if (j+1<m) {
        const Z *Lij=L.addr(j+1,j);
        Z *Sij=S->addr(j+1,j);
        for (int i=j+1;i<m;i++,Lij++,Sij++) *Sij=(*Lij)*dj;
      }
    }
  }
  return S;
}

template<typename F,typename Z> SquareMatrix<F,Z>*
DiagonalMatrix<F,Z>::operator*(
const SquareMatrix<F,Z> &M) const {
  CHECK_SAME(dim,M.size(0));
  SquareMatrix<F,Z> *S=
    OPERATOR_NEW SquareMatrix<F,Z>(dim);
  for (int i=0;i<dim;i++) {
    Z di=D->operator[](i);
    const Z *Mij=M.addr(i,0);
    Z *Sij=S->addr(i,0);
    for (int j=0;j<dim;j++,Mij+=dim,Sij+=dim) *Sij=di*(*Mij);
  }
  return S;
}

template<typename F,typename Z> SquareMatrix<F,Z>* operator*(
const SquareMatrix<F,Z> &M,
const DiagonalMatrix<F,Z> &A) {
  int n=A.size(0);
  CHECK_SAME(n,M.size(0));
  SquareMatrix<F,Z> *S=
    OPERATOR_NEW SquareMatrix<F,Z>(n);
  for (int j=0;j<n;j++) {
    Z dj=A[j];
    const Z *Mij=M.addr(0,j);
    Z *Sij=S->addr(0,j);
    for (int i=0;i<n;i++,Mij++,Sij++) *Sij=(*Mij)*dj;
  }
  return S;
}

template<typename F,typename Z> Matrix<F,Z>*
DiagonalMatrix<F,Z>::operator*(const Matrix<F,Z> &M) const {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(m,dim);
  Matrix<F,Z> *S=OPERATOR_NEW Matrix<F,Z>(m,n);
  for (int i=0;i<m;i++) {
    Z di=D->operator[](i);
    const Z *Mij=M.addr(i,0);
    Z *Sij=S->addr(i,0);
    for (int j=0;j<n;j++,Mij+=m,Sij+=m) *Sij=di*(*Mij);
  }
  return S;
}

template<typename F,typename Z> Matrix<F,Z>* operator*(
const Matrix<F,Z> &M,const DiagonalMatrix<F,Z> &A) {
  int m=M.size(0),n=M.size(1);
  CHECK_SAME(n,A.size(0));
  Matrix<F,Z> *S=OPERATOR_NEW Matrix<F,Z>(m,n);
  for (int j=0;j<n;j++) {
    Z dj=A[j];
    const Z *Mij=M.addr(0,j);
    Z *Sij=S->addr(0,j);
    for (int i=0;i<m;i++,Mij++,Sij++) *Sij=(*Mij)*dj;
  }
  return S;
}

template<typename F,typename Z> Vector<F,Z>*
DiagonalMatrix<F,Z>::operator*(const Vector<F,Z> &v) const {
  CHECK_SAME(dim,v.size());
  Vector<F,Z> *s=OPERATOR_NEW Vector<F,Z>(dim);
  const Z *di=D->addr();
  const Z *vi=v.addr();
  Z *si=s->addr();
  for (int i=0;i<dim;i++,di++,vi++,si++) *si=(*di)*(*vi);
  return s;
}

template<typename F,typename Z> void DiagonalMatrix<F,Z>::dmv(Z alpha,
const Vector<F,Z> &x,Z beta,Vector<F,Z> &b)
const {
  CHECK_SAME(dim,x.size());
  CHECK_SAME(dim,b.size());
  const Z *di=addr();
  const Z *xi=x.addr();
  Z *bi=b.addr();
  for (int i=0;i<dim;i++,di++,xi++,bi++) {
    *bi=(*di)*(*xi)*alpha+(*bi)*beta; 
  }
}

template<typename F,typename Z> void DiagonalMatrix<F,Z>::dmm(Z alpha,
const Matrix<F,Z> &X,Z beta,Matrix<F,Z> &B,char side) const {
  int m=X.size(0),n=X.size(1);
  CHECK_SAME(m,B.size(0));
  CHECK_SAME(n,B.size(1));
  if (side=='L' || side=='l') {
    CHECK_SAME(m,dim);
    for (int j=0;j<n;j++) {
      const Z *Di=addr();
      const Z *Xij=X.addr(0,j);
      Z *Bij=B.addr(0,j);
      for (int i=0;i<m;i++,Di++,Xij++,Bij++) {
        *Bij=(*Di)*(*Xij)*alpha+(*Bij)*beta; 
      }
    }
  } else {
    CHECK_SAME(n,dim);
    const Z *Dj=addr();
    for (int j=0;j<n;j++,Dj++) {
      const Z *Xij=X.addr(0,j);
      Z *Bij=B.addr(0,j);
      Z t=(*Dj)*alpha;
      for (int i=0;i<m;i++,Xij++,Bij++) {
        *Bij=(*Xij)*t+(*Bij)*beta; 
      }
    }
  }
}

template<typename F,typename Z> void DiagonalMatrix<F,Z>::solve(
const Vector<F,Z> &b,Vector<F,Z> &x) const {
  CHECK_SAME(dim,b.size());
  CHECK_SAME(dim,x.size());
  for (int i=0;i<dim;i++) {
    Z di=D->operator[](i);
    CHECK_NONZERO(abs(di));
    x[i]=b[i]/di;
  }
}

template<typename F,typename Z> void DiagonalMatrix<F,Z>::solve(
const Matrix<F,Z> &B,Matrix<F,Z> &X,char side) const {
  int m=B.size(0),n=B.size(1);
  if (side=='L' || side=='l') {
    CHECK_SAME(m,dim);
    for (int i=0;i<m;i++) {
      Z di=D->operator[](i);
      CHECK_NONZERO(abs(di));
      const Z *Bij=B.addr(i,0);
      Z *Xij=X.addr(i,0);
      for (int j=0;j<n;j++,Bij+=m,Xij+=m) *Xij=(*Bij)/di;
    }
  } else {
    CHECK_SAME(n,dim);
    for (int j=0;j<n;j++) {
      Z dj=D->operator[](j);
      CHECK_NONZERO(abs(dj));
      const Z *Bij=B.addr(0,j);
      Z *Xij=X.addr(0,j);
      for (int i=0;i<m;i++,Bij++,Xij++) *Xij=(*Bij)/dj;
    }
  }
}

template<typename F,typename Z> void
DiagonalMatrix<F,Z>::printOn(ostream& s) const {
  int n=this->size(0);
  s << "DiagonalMatrix(" << n << " x " << n << ")\n" ;
  for (int i=0;i<n;i++) {
    for (int j=0;j<i;j++) s << Matrix<F,Z>::zero_ << " ";
    s << (*D)[i] << "  ";
    for (int j=i+1;j<n;j++) s << Matrix<F,Z>::zero_ << " ";
    s << endl;
  }
}

template<typename F,typename Z> void testDiagonalMatrix(F fscalar,
Z scalar) {
  DiagonalMatrix<F,Z> *T1=OPERATOR_NEW DiagonalMatrix<F,Z>;
  cout << "\nafter DiagonalMatrix()" << endl;
  T1->printOn(cout);
  delete T1; T1=0;

  int n=3;
  cout << "\nn = " << n << endl;
  DiagonalMatrix<F,Z> *T2=OPERATOR_NEW DiagonalMatrix<F,Z>(n);
  cout << "\nafter DiagonalMatrix(int)" << endl;
  T2->printOn(cout);

  DiagonalMatrix<F,Z> *T3=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  cout << "\nafter DiagonalMatrix(int,scalar)" << endl;
  T3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=T3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }

  *T2=scalar;
  cout << "\nafter T2=scalar" << endl;
  T2->printOn(cout);

  T3->resize(n+1);
  cout << "\nafter resize(n)" << endl;
  T3->printOn(cout);
  
  T3->resize(*T2);
  cout << "\nafter resize(T)" << endl;
  T3->printOn(cout);
  
  T3->copy(*T2);
  cout << "\nafter copy(T)" << endl;
  T3->printOn(cout);

  cout << "\nT3+T2:" << endl;
  DiagonalMatrix<F,Z> *T=(*T3)+(*T2);
  T->printOn(cout);
  delete T; T=0;

  cout << "\nDiagonalMatrix + SymmetricPositiveTridiagonalMatrix" << endl;
  DiagonalMatrix<F,F> *DR=OPERATOR_NEW DiagonalMatrix<F,F>(n,fscalar);
  SymmetricPositiveTridiagonalMatrix<F,Z> *Sptr=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  SymmetricTridiagonalMatrix<F,Z> *Str2=(*DR)+(*Sptr);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + DiagonalMatrix" << endl;
  Str2=(*Sptr)+(*DR);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  if (sizeof(F)<sizeof(Z)) {
    cout << "\nDiagonalMatrix + SymmetricPositiveTridiagonalMatrix"
         << endl;
    void *V=(*T2)+(*Sptr);
    TridiagonalMatrix<F,Z> *TDM=
      reinterpret_cast<TridiagonalMatrix<F,Z>*>(V);
    TDM->printOn(cout);
    delete TDM; TDM=0;

    cout << "\nSymmetricPositiveTridiagonalMatrix + DiagonalMatrix"
         << endl;
    V=(*Sptr)+(*T2);
    TDM=reinterpret_cast<TridiagonalMatrix<F,Z>*>(V);
    TDM->printOn(cout);
    delete TDM; TDM=0;
  }
  delete Sptr; Sptr=0;

  cout << "\nDiagonalMatrix + SymmetricTridiagonalMatrix" << endl;
  SymmetricTridiagonalMatrix<F,Z> *Str=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  Str2=(*DR)+(*Str);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  cout << "\nSymmetricTridiagonalMatrix + DiagonalMatrix" << endl;
  Str2=(*Str)+(*DR);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  if (sizeof(F)<sizeof(Z)) {
    cout << "\nDiagonalMatrix + SymmetricTridiagonalMatrix" << endl;
    void *V=(*T2)+(*Str);
    TridiagonalMatrix<F,Z> *TDM=
      reinterpret_cast<TridiagonalMatrix<F,Z>*>(V);
    TDM->printOn(cout);
    delete TDM; TDM=0;

    cout << "\nSymmetricTridiagonalMatrix + DiagonalMatrix" << endl;
    V=(*Str)+(*T2);
    TDM=reinterpret_cast<TridiagonalMatrix<F,Z>*>(V);
    TDM->printOn(cout);
    delete TDM; TDM=0;
  }
  delete Str; Str=0;

  cout << "\nDiagonalMatrix + TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  TridiagonalMatrix<F,Z> *Tr2=(*T2)+(*Tr);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nTridiagonalMatrix + DiagonalMatrix" << endl;
  Tr2=(*Tr)+(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Tr; Tr=0;

  cout << "\nDiagonalMatrix + SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *Syp=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *Sy2=(*DR)+(*Syp);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  cout << "\nSymmetricPositive Matrix + DiagonalMatrix" << endl;
  Sy2=(*Syp)+(*DR);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  if (sizeof(F)<sizeof(Z)) {
    cout << "\nDiagonalMatrix + SymmetricPositiveMatrix" << endl;
    void *V=(*T2)+(*Syp);
    SquareMatrix<F,Z> *SM=reinterpret_cast<SquareMatrix<F,Z>*>(V);
    SM->printOn(cout);
    delete SM; SM=0;

    cout << "\nSymmetricPositiveMatrix + DiagonalMatrix" << endl;
    V=(*Syp)+(*T2);
    SM=reinterpret_cast<SquareMatrix<F,Z>*>(V);
    SM->printOn(cout);
    delete SM; SM=0;
  }
  delete Syp; Syp=0;

  cout << "\nDiagonalMatrix + SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sy2=(*DR)+(*Sy);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  cout << "\nSymmetricMatrix + DiagonalMatrix" << endl;
  Sy2=(*Sy)+(*DR);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  if (sizeof(F)<sizeof(Z)) {
    cout << "\nDiagonalMatrix + SymmetricMatrix" << endl;
    void *V=(*T2)+(*Sy);
    SquareMatrix<F,Z> *SM=reinterpret_cast<SquareMatrix<F,Z>*>(V);
    SM->printOn(cout);
    delete SM; SM=0;

    cout << "\nSymmetricMatrix + DiagonalMatrix" << endl;
    V=(*Sy)+(*T2);
    SM=reinterpret_cast<SquareMatrix<F,Z>*>(V);
    SM->printOn(cout);
    delete SM; SM=0;
  }
  delete Sy; Sy=0;

  cout << "\nDiagonalMatrix + OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Orth=
    OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  Matrix<F,Z> *M2=(*T2)+(*Orth);
  M2->printOn(cout);
  delete M2; M2=0;

  cout << "\nOrthogonalMatrix + DiagonalMatrix" << endl;
  M2=(*Orth)+(*T2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Orth; Orth=0;

  cout << "\nDiagonalMatrix + UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *Uutg=OPERATOR_NEW 
    UnitUpperTriangularMatrix<F,Z>(n,scalar);
  UpperTriangularMatrix<F,Z> *Utg2=(*T2)+(*Uutg);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;

  cout << "\nUnitUpperTriangularMatrix + DiagonalMatrix" << endl;
  Utg2=(*Uutg)+(*T2);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;
  delete Uutg; Uutg=0;

  cout << "\nDiagonalMatrix + UnitUpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *Uutp=OPERATOR_NEW 
    UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Utg2=(*T2)+(*Uutp);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;

  cout << "\nUnitUpperTrapezoidalMatrix + DiagonalMatrix" << endl;
  Utg2=(*Uutp)+(*T2);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;
  delete Uutp; Uutp=0;

  cout << "\nDiagonalMatrix + UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *Utg=OPERATOR_NEW 
    UpperTriangularMatrix<F,Z>(n,scalar);
  Utg2=(*T2)+(*Utg);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;

  cout << "\nUpperTriangularMatrix + DiagonalMatrix" << endl;
  Utg2=(*Utg)+(*T2);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;
  delete Utg; Utg=0;

  cout << "\nDiagonalMatrix + UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *Utp=OPERATOR_NEW 
    UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Utg2=(*T2)+(*Utp);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;

  cout << "\nUpperTrapezoidalMatrix + DiagonalMatrix" << endl;
  Utg2=(*Utp)+(*T2);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;
  delete Utp; Utp=0;

  cout << "\nDiagonalMatrix + UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *Ultg=OPERATOR_NEW 
    UnitLowerTriangularMatrix<F,Z>(n,scalar);
  LowerTriangularMatrix<F,Z> *Ltg2=(*T2)+(*Ultg);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;

  cout << "\nUnitLowerTriangularMatrix + DiagonalMatrix" << endl;
  Ltg2=(*Ultg)+(*T2);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;
  delete Ultg; Ultg=0;

  cout << "\nDiagonalMatrix + UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *Ultp=OPERATOR_NEW 
    UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Ltg2=(*T2)+(*Ultp);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;

  cout << "\nUnitLowerTrapezoidalMatrix + DiagonalMatrix" << endl;
  Ltg2=(*Ultp)+(*T2);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;
  delete Ultp; Ultp=0;

  cout << "\nDiagonalMatrix + LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *Ltg=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  Ltg2=(*T2)+(*Ltg);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;

  cout << "\nLowerTriangularMatrix + DiagonalMatrix" << endl;
  Ltg2=(*Ltg)+(*T2);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;
  delete Ltg; Ltg=0;

  cout << "\nDiagonalMatrix + LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *Ltp=OPERATOR_NEW 
    LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Ltg2=(*T2)+(*Ltp);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;

  cout << "\nLowerTrapezoidalMatrix + DiagonalMatrix" << endl;
  Ltg2=(*Ltp)+(*T2);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;
  delete Ltp; Ltp=0;

  cout << "\nDiagonalMatrix + SquareMatrix" << endl;
  SquareMatrix<F,Z> *SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *Sq=(*T2)+(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix + DiagonalMatrix" << endl;
  Sq=(*SQ)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nDiagonalMatrix + Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*T2)+(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix + DiagonalMatrix" << endl;
  Sq=(*M)+(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nDiagonalMatrix - SymmetricPositiveTridiagonalMatrix" << endl;
  Sptr=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  Str2=(*DR)-(*Sptr);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - DiagonalMatrix" << endl;
  Str2=(*Sptr)-(*DR);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  if (sizeof(F)<sizeof(Z)) {
    cout << "\nDiagonalMatrix - SymmetricPositiveTridiagonalMatrix"
         << endl;
    void *V=(*T2)-(*Sptr);
    TridiagonalMatrix<F,Z> *TDM=
      reinterpret_cast<TridiagonalMatrix<F,Z>*>(V);
    TDM->printOn(cout);
    delete TDM; TDM=0;

    cout << "\nSymmetricPositiveTridiagonalMatrix - DiagonalMatrix"
         << endl;
    V=(*Sptr)-(*T2);
    TDM=reinterpret_cast<TridiagonalMatrix<F,Z>*>(V);
    TDM->printOn(cout);
    delete TDM; TDM=0;
  }
  delete Sptr; Sptr=0;

  cout << "\nDiagonalMatrix - SymmetricTridiagonalMatrix" << endl;
  Str=OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  Str2=(*DR)-(*Str);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  cout << "\nSymmetricTridiagonalMatrix - DiagonalMatrix" << endl;
  Str2=(*Str)-(*DR);
  Str2->printOn(cout);
  delete Str2; Str2=0;

  if (sizeof(F)<sizeof(Z)) {
    cout << "\nDiagonalMatrix - SymmetricTridiagonalMatrix" << endl;
    void *V=(*T2)-(*Str);
    TridiagonalMatrix<F,Z> *TDM=
      reinterpret_cast<TridiagonalMatrix<F,Z>*>(V);
    TDM->printOn(cout);
    delete TDM; TDM=0;

    cout << "\nSymmetricTridiagonalMatrix - DiagonalMatrix" << endl;
    V=(*Str)-(*T2);
    TDM=reinterpret_cast<TridiagonalMatrix<F,Z>*>(V);
    TDM->printOn(cout);
    delete TDM; TDM=0;
  }
  delete Str; Str=0;

  cout << "\nDiagonalMatrix - TridiagonalMatrix" << endl;
  Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  Tr2=(*T2)-(*Tr);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nTridiagonalMatrix - DiagonalMatrix" << endl;
  Tr2=(*Tr)-(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Tr; Tr=0;

  cout << "\nDiagonalMatrix - SymmetricPositiveMatrix" << endl;
  Syp=OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  Sy2=(*DR)-(*Syp);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  cout << "\nSymmetricPositive Matrix - DiagonalMatrix" << endl;
  Sy2=(*Syp)-(*DR);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  if (sizeof(F)<sizeof(Z)) {
    cout << "\nDiagonalMatrix - SymmetricPositiveMatrix" << endl;
    void *V=(*T2)-(*Syp);
    SquareMatrix<F,Z> *SM=reinterpret_cast<SquareMatrix<F,Z>*>(V);
    SM->printOn(cout);
    delete SM; SM=0;

    cout << "\nSymmetricPositiveMatrix - DiagonalMatrix" << endl;
    V=(*Syp)-(*T2);
    SM=reinterpret_cast<SquareMatrix<F,Z>*>(V);
    SM->printOn(cout);
    delete SM; SM=0;
  }
  delete Syp; Syp=0;

  cout << "\nDiagonalMatrix - SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sy2=(*DR)-(*Sy);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  cout << "\nSymmetricMatrix - DiagonalMatrix" << endl;
  Sy2=(*Sy)-(*DR);
  Sy2->printOn(cout);
  delete Sy2; Sy2=0;

  if (sizeof(F)<sizeof(Z)) {
    cout << "\nDiagonalMatrix - SymmetricMatrix" << endl;
    void *V=(*T2)-(*Sy);
    SquareMatrix<F,Z> *SM=reinterpret_cast<SquareMatrix<F,Z>*>(V);
    SM->printOn(cout);
    delete SM; SM=0;

    cout << "\nSymmetricMatrix - DiagonalMatrix" << endl;
    V=(*Sy)-(*T2);
    SM=reinterpret_cast<SquareMatrix<F,Z>*>(V);
    SM->printOn(cout);
    delete SM; SM=0;
  }
  delete Sy; Sy=0;
  delete DR; DR=0;

  cout << "\nDiagonalMatrix - OrthogonalMatrix" << endl;
  Orth=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  M2=(*T2)-(*Orth);
  M2->printOn(cout);
  delete M2; M2=0;

  cout << "\nOrthogonalMatrix - DiagonalMatrix" << endl;
  M2=(*Orth)-(*T2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Orth; Orth=0;

  cout << "\nDiagonalMatrix - UnitUpperTriangularMatrix" << endl;
  Uutg=OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  Utg2=(*T2)-(*Uutg);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;

  cout << "\nUnitUpperTriangularMatrix - DiagonalMatrix" << endl;
  Utg2=(*Uutg)-(*T2);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;
  delete Uutg; Uutg=0;

  cout << "\nDiagonalMatrix - UnitUpperTrapezoidalMatrix" << endl;
  Uutp=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Utg2=(*T2)-(*Uutp);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;

  cout << "\nUnitUpperTrapezoidalMatrix - DiagonalMatrix" << endl;
  Utg2=(*Uutp)-(*T2);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;
  delete Uutp; Uutp=0;

  cout << "\nDiagonalMatrix - UpperTriangularMatrix" << endl;
  Utg=OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  Utg2=(*T2)-(*Utg);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;

  cout << "\nUpperTriangularMatrix - DiagonalMatrix" << endl;
  Utg2=(*Utg)-(*T2);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;
  delete Utg; Utg=0;

  cout << "\nDiagonalMatrix - UpperTrapezoidalMatrix" << endl;
  Utp=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  Utg2=(*T2)-(*Utp);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;

  cout << "\nUnitUpperTrapezoidalMatrix - DiagonalMatrix" << endl;
  Utg2=(*Utp)-(*T2);
  Utg2->printOn(cout);
  delete Utg2; Utg2=0;
  delete Utp; Utp=0;

  cout << "\nDiagonalMatrix - UnitLowerTriangularMatrix" << endl;
  Ultg=OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  Ltg2=(*T2)-(*Ultg);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;

  cout << "\nUnitLowerTriangularMatrix - DiagonalMatrix" << endl;
  Ltg2=(*Ultg)-(*T2);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;
  delete Ultg; Ultg=0;

  cout << "\nDiagonalMatrix - UnitLowerTrapezoidalMatrix" << endl;
  Ultp=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Ltg2=(*T2)-(*Ultp);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;

  cout << "\nUnitLowerTrapezoidalMatrix - DiagonalMatrix" << endl;
  Ltg2=(*Ultp)-(*T2);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;
  delete Ultp; Ultp=0;

  cout << "\nDiagonalMatrix - LowerTriangularMatrix" << endl;
  Ltg=OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  Ltg2=(*T2)-(*Ltg);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;

  cout << "\nLowerTriangularMatrix - DiagonalMatrix" << endl;
  Ltg2=(*Ltg)-(*T2);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;
  delete Ltg; Ltg=0;

  cout << "\nDiagonalMatrix - LowerTrapezoidalMatrix" << endl;
  Ltp=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Ltg2=(*T2)-(*Ltp);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;

  cout << "\nUnitLowerTrapezoidalMatrix - DiagonalMatrix" << endl;
  Ltg2=(*Ltp)-(*T2);
  Ltg2->printOn(cout);
  delete Ltg2; Ltg2=0;
  delete Ltp; Ltp=0;

  cout << "\nDiagonalMatrix - SquareMatrix" << endl;
  SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*T2)-(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix - DiagonalMatrix" << endl;
  Sq=(*SQ)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nDiagonalMatrix - Matrix" << endl;
  M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*T2)-(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix - DiagonalMatrix" << endl;
  Sq=(*M)-(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nT3*scalar :" << endl;
  DiagonalMatrix<F,Z> *D=(*T3)*scalar;
  D->printOn(cout);
  delete D; D=0;
  
  cout << "\nT3/scalar.:" << endl;
  D=(*T3)/scalar;
  D->printOn(cout);
  delete D; D=0;

  cout << "\nDiagonalMatrix * DiagonalMatrix" << endl;
  D=(*T2)*(*T3);
  D->printOn(cout);
  delete D; D=0;
  delete T3; T3=0;

  cout << "\nDiagonalMatrix * SymmetricPositiveTridiagonalMatrix" << endl;
  Sptr=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  Tr2=(*T2)*(*Sptr);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix * DiagonalMatrix" << endl;
  Tr2=(*Sptr)*(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Sptr; Sptr=0;

  cout << "\nDiagonalMatrix * SymmetricTridiagonalMatrix" << endl;
  Str=OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  Tr2=(*T2)*(*Str);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nSymmetricTridiagonalMatrix * DiagonalMatrix" << endl;
  Tr2=(*Str)*(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Str; Str=0;

  cout << "\nDiagonalMatrix * TridiagonalMatrix" << endl;
  Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  Tr2=(*T2)*(*Tr);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;

  cout << "\nTridiagonalMatrix * DiagonalMatrix" << endl;
  Tr2=(*Tr)*(*T2);
  Tr2->printOn(cout);
  delete Tr2; Tr2=0;
  delete Tr; Tr=0;

  cout << "\nDiagonalMatrix * SymmetricPositiveMatrix" << endl;
  Syp=OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  Sq=(*T2)*(*Syp);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricPositiveMatrix * DiagonalMatrix" << endl;
  Sq=(*Syp)*(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Syp; Syp=0;

  cout << "\nDiagonalMatrix * SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sq=(*T2)*(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix * DiagonalMatrix" << endl;
  Sq=(*Sy)*(*T2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  Ultg= OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(3,scalar);
  cout << "\nDiagonalMatrix * UnitLowerTriangularMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *Ltp2=(*T2)*(*Ultg);
  Ltp2->printOn(cout);
  delete Ltp2; Ltp2=0;

  cout << "\nUnitLowerTriangularMatrix * DiagonalMatrix" << endl;
  Ltp2=(*Ultg)*(*T2);
  Ltp2->printOn(cout);
  delete Ltp2; Ltp2=0;
  delete Ultg; Ultg=0;

  Ultp=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(3,2,scalar);
  cout << "\nDiagonalMatrix * UnitLowerTrapezoidalMatrix" << endl;
  Ltp2=(*T2)*(*Ultp);
  Ltp2->printOn(cout);
  delete Ltp2; Ltp2=0;
  delete Ultp; Ultp=0;

  Ultp=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(4,3,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * DiagonalMatrix" << endl;
  Ltp2=(*Ultp)*(*T2);
  Ltp2->printOn(cout);
  delete Ltp2; Ltp2=0;
  delete Ultp; Ultp=0;

  LowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(3,2,scalar);
  cout << "\nDiagonalMatrix * LowerTrapezoidalMatrix" << endl;
  Ltp2=(*T2)*(*L);
  Ltp2->printOn(cout);
  delete Ltp2; Ltp2=0;
  delete L; L=0;

  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(4,3,scalar);
  cout << "\nLowerTrapezoidalMatrix * DiagonalMatrix" << endl;
  Ltp2=(*L)*(*T2);
  Ltp2->printOn(cout);
  delete Ltp2; Ltp2=0;
  delete L; L=0;

  Uutg=OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(3,scalar);
  cout << "\nDiagonalMatrix * UnitUpperTriangularMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *Utp2=(*T2)*(*Uutg);
  Utp2->printOn(cout);
  delete Utp2; Utp2=0;

  cout << "\nUnitUpperTriangularMatrix * DiagonalMatrix" << endl;
  Utp2=(*Uutg)*(*T2);
  Utp2->printOn(cout);
  delete Utp2; Utp2=0;
  delete Uutg; Uutg=0;

  Uutp=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(3,4,scalar);
  cout << "\nDiagonalMatrix * UnitUpperTrapezoidalMatrix" << endl;
  Utp2=(*T2)*(*Uutp);
  Utp2->printOn(cout);
  delete Utp2; Utp2=0;
  delete Uutp; Uutp=0;

  Uutp=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(2,3,scalar);
  cout << "\nUnitUpperTrapezoidalMatrix * DiagonalMatrix" << endl;
  Utp2=(*Uutp)*(*T2);
  Utp2->printOn(cout);
  delete Utp2; Utp2=0;
  delete Uutp; Uutp=0;

  UpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(3,4,scalar);
  cout << "\nDiagonalMatrix * UpperTrapezoidalMatrix" << endl;
  Utp2=(*T2)*(*U);
  Utp2->printOn(cout);
  delete Utp2; Utp2=0;
  delete U; U=0;

  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(2,3,scalar);
  cout << "\nUpperTrapezoidalMatrix * DiagonalMatrix" << endl;
  Utp2=(*U)*(*T2);
  Utp2->printOn(cout);
  delete Utp2; Utp2=0;
  delete U; U=0;

  Sq=OPERATOR_NEW SquareMatrix<F,Z>(3,scalar);
  cout << "\nDiagonalMatrix * SquareMatrix" << endl;
  SquareMatrix<F,Z> *Sq2=(*T2)*(*Sq);
  Sq2->printOn(cout);
  delete Sq2; Sq2=0;

  cout << "\nSquareMatrix * DiagonalMatrix" << endl;
  Sq2=(*Sq)*(*T2);
  Sq2->printOn(cout);
  delete Sq2; Sq2=0;
  delete Sq; Sq=0;

  M=OPERATOR_NEW Matrix<F,Z>(3,2,scalar);
  cout << "\nDiagonalMatrix * Matrix" << endl;
  M2=(*T2)*(*M);
  M2->printOn(cout);
  delete M2; M2=0;
  delete M; M=0;

  M=OPERATOR_NEW Matrix<F,Z>(4,3,scalar);
  cout << "\nMatrix * DiagonalMatrix" << endl;
  M2=(*M)*(*T2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete M; M=0;

  Vector<F,Z> *v=OPERATOR_NEW Vector<F,Z>(3,scalar);
  cout << "\nDiagonalMatrix * Vector" << endl;
  Vector<F,Z> *w=(*T2)*(*v);
  w->printOn(cout);
  delete w; w=0;
  delete v; v=0;
  delete T2; T2=0;

  DiagonalMatrix<F,Z> *TT=OPERATOR_NEW DiagonalMatrix<F,Z>(3);
  for (int i=0;i<3;i++) (*TT)[i]=F(i+1)*scalar;
  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(3);
  for (int i=0;i<3;i++) (*b)[i]=F((i+1)*(i+1))*scalar;
  cout << "\nsquare system of linear equations" << endl;
  TT->printOn(cout);
  b->printOn(cout);
  cout << "solve:" << endl;
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(3);
  TT->solve(*b,*x);
  x->printOn(cout);
  Vector<F,Z> *r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  TT->dmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r);
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;
  delete b; b=0;
 
  Matrix<F,Z> *BB=OPERATOR_NEW Matrix<F,Z>(3,1);
  for (int i=0;i<3;i++) (*BB)(i,0)=F((i+1)*(i+1))*scalar;
  cout << "\nsolve('L'):" << endl;
  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(3,1);
  TT->solve(*BB,*X);
  X->printOn(cout);
  Matrix<F,Z> *R=OPERATOR_NEW Matrix<F,Z>(BB->size(0),BB->size(1));
  R->copy(*BB);
  TT->dmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete BB; BB=0;

  BB=OPERATOR_NEW Matrix<F,Z>(1,3);
  for (int i=0;i<3;i++) (*BB)(0,i)=F((i+1)*(i+1))*scalar;
  cout << "\nsolve('R'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  TT->solve(*BB,*X,'R');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(BB->size(0),BB->size(1));
  R->copy(*BB);
  TT->dmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete BB; BB=0;
  delete TT; TT=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
UpperHessenbergMatrix<F,Z>::fillWith(Z d) {
  int n=this->size(0);
  Z *colj=addr();
  for (int j=0;j<n;j++,colj+=n) {
    for (int i=0;i<=min(n-1,j+1);i++) colj[i]=d;
  }
}

template<typename F,typename Z> 
UpperHessenbergMatrix<F,Z>::UpperHessenbergMatrix(
const TridiagonalMatrix<F,Z> &T) : SquareMatrix<F,Z>(T.size(0)) {
  int n=this->size(0);
  Z *colj=addr();
  for (int j=0;j<n;j++,colj+=n) {
    for (int i=0;i<j-1;i++) colj[i]=this->zero_;
    if (j>0) colj[j-1]=T(j-1,j);
    colj[j]=T(j,j);
    if (j<n-1) colj[j+1]=T(j+1,j);
  }
}

template<typename F,typename Z> void
UpperHessenbergMatrix<F,Z>::printOn(ostream& s) const {
  int n=this->size(0);
  s << "UpperHessenbergMatrix(" << n << " x " << n << ")\n" ;
  for (int i=0; i<n; i++) {
    for (int j=0;j<n; j++) s << operator()(i,j) << "  ";
    s << endl;
  }
}

template<typename F,typename Z> void
testUpperHessenbergMatrix(F fscalar,Z scalar) {
  UpperHessenbergMatrix<F,Z> *H1=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>;
  cout << "\nafter UpperHessenbergMatrix()" << endl;
  H1->printOn(cout);
  delete H1; H1=0;

  int n=4;
  cout << "\nn = " << n << endl;
  UpperHessenbergMatrix<F,Z> *H2=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n);
  cout << "\nafter UpperHessenbergMatrix(int)" << endl;
  H2->printOn(cout);

  UpperHessenbergMatrix<F,Z> *H3=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  cout << "\nafter UpperHessenbergMatrix(int,scalar)" << endl;
  H3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=H3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }

  TridiagonalMatrix<F,Z> *T=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  UpperHessenbergMatrix<F,Z> *H4=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(*T);
  cout << "\nafter UpperHessenbergMatrix(TridiagonalMatrix)" << endl;
  H4->printOn(cout);
  delete H4; H4=0;
  delete T; T=0;

  *H2=scalar;
  cout << "\nafter H2=scalar" << endl;
  H2->printOn(cout);

  H3->resize(n+1);
  cout << "\nafter resize(n)" << endl;
  H3->printOn(cout);
  
  H3->resize(*H2);
  cout << "\nafter resize(H2)" << endl;
  H3->printOn(cout);
  
  H3->copy(*H2);
  cout << "\nafter copy(H2)" << endl;
  H3->printOn(cout);

  cout << "\nH3+H2:" << endl;
  UpperHessenbergMatrix<F,Z> *HT=(*H3)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUpperHessenbergMatrix + DiagonalMatrix:" << endl;
  DiagonalMatrix<F,Z> *D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  HT=(*H2)+(*D);
  HT->printOn(cout);
  delete HT; HT=0;
  delete D; D=0;

  cout << "\nDiagonalMatrix + UpperHessenbergMatrix:" << endl;
  D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  HT=(*D)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete D; D=0;

  cout << "\nUpperHessenbergMatrix + SymmetricPositiveTridiagonalMatrix"
       << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *Sptr=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  HT=(*H2)+(*Sptr);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + UpperHessenbergMatrix"
       << endl;
  HT=(*Sptr)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Sptr; Sptr=0;

  cout << "\nUpperHessenbergMatrix + SymmetricTridiagonalMatrix" << endl;
  SymmetricTridiagonalMatrix<F,Z> *Str=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  HT=(*H2)+(*Str);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nSymmetricTridiagonalMatrix + UpperHessenbergMatrix" << endl;
  HT=(*Str)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Str; Str=0;

  cout << "\nUpperHessenbergMatrix + TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  HT=(*H2)+(*Tr);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nTridiagonalMatrix + UpperHessenbergMatrix" << endl;
  HT=(*Tr)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Tr; Tr=0;

  cout << "\nUpperHessenbergMatrix + SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *Syp=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *Sq=(*H2)+(*Syp);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricPositive Matrix + UpperHessenbergMatrix" << endl;
  Sq=(*Syp)+(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Syp; Syp=0;

  cout << "\nUpperHessenbergMatrix + SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sq=(*H2)+(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix + UpperHessenbergMatrix" << endl;
  Sq=(*Sy)+(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  cout << "\nUpperHessenbergMatrix + OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Orth=
    OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  Matrix<F,Z> *M2=(*H2)+(*Orth);
  M2->printOn(cout);
  delete M2; M2=0;

  cout << "\nOrthogonalMatrix + UpperHessenbergMatrix" << endl;
  M2=(*Orth)+(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Orth; Orth=0;

  cout << "\nUpperHessenbergMatrix + UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *Uutg=OPERATOR_NEW 
    UnitUpperTriangularMatrix<F,Z>(n,scalar);
  HT=(*H2)+(*Uutg);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUnitUpperTriangularMatrix + UpperHessenbergMatrix" << endl;
  HT=(*Uutg)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Uutg; Uutg=0;

  cout << "\nUpperHessenbergMatrix + UnitUpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *Uutp=OPERATOR_NEW 
    UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  HT=(*H2)+(*Uutp);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUnitUpperTrapezoidalMatrix + UpperHessenbergMatrix" << endl;
  HT=(*Uutp)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Uutp; Uutp=0;

  cout << "\nUpperHessenbergMatrix + UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *Utg=OPERATOR_NEW 
    UpperTriangularMatrix<F,Z>(n,scalar);
  HT=(*H2)+(*Utg);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUpperTriangularMatrix + UpperHessenbergMatrix" << endl;
  HT=(*Utg)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Utg; Utg=0;

  cout << "\nUpperHessenbergMatrix + UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *Utp=OPERATOR_NEW 
    UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  HT=(*H2)+(*Utp);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUpperTrapezoidalMatrix + UpperHessenbergMatrix" << endl;
  HT=(*Utp)+(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Utp; Utp=0;

  cout << "\nUpperHessenbergMatrix + UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *Ultg=OPERATOR_NEW 
    UnitLowerTriangularMatrix<F,Z>(n,scalar);
  Sq=(*H2)+(*Ultg);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitLowerTriangularMatrix + UpperHessenbergMatrix" << endl;
  Sq=(*Ultg)+(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Ultg; Ultg=0;

  cout << "\nUpperHessenbergMatrix + UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *Ultp=OPERATOR_NEW 
    UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*H2)+(*Ultp);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitLowerTrapezoidalMatrix + UpperHessenbergMatrix" << endl;
  Sq=(*Ultp)+(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Ultp; Ultp=0;

  cout << "\nUpperHessenbergMatrix + LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *Ltg=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  Sq=(*H2)+(*Ltg);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTriangularMatrix + UpperHessenbergMatrix" << endl;
  Sq=(*Ltg)+(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Ltg; Ltg=0;

  cout << "\nUpperHessenbergMatrix + LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *Ltp=OPERATOR_NEW 
    LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*H2)+(*Ltp);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTrapezoidalMatrix + UpperHessenbergMatrix" << endl;
  Sq=(*Ltp)+(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Ltp; Ltp=0;

  cout << "\nUpperHessenbergMatrix + SquareMatrix" << endl;
  SquareMatrix<F,Z> *SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*H2)+(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix + UpperHessenbergMatrix" << endl;
  Sq=(*SQ)+(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nUpperHessenbergMatrix + Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*H2)+(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix + UpperHessenbergMatrix" << endl;
  Sq=(*M)+(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nH3-H2:" << endl;
  HT=(*H3)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUpperHessenbergMatrix - DiagonalMatrix" << endl;
  D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  HT=(*H2)-(*D);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nDiagonalMatrix - UpperHessenbergMatrix" << endl;
  HT=(*D)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete D; D=0;

  cout << "\nUpperHessenbergMatrix - SymmetricPositiveTridiagonalMatrix"
       << endl;
  Sptr=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  HT=(*H2)-(*Sptr);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - UpperHessenbergMatrix"
       << endl;
  HT=(*Sptr)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Sptr; Sptr=0;

  cout << "\nUpperHessenbergMatrix - SymmetricTridiagonalMatrix" << endl;
  Str=OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  HT=(*H2)-(*Str);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nSymmetricTridiagonalMatrix - UpperHessenbergMatrix" << endl;
  HT=(*Str)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Str; Str=0;

  cout << "\nUpperHessenbergMatrix - TridiagonalMatrix" << endl;
  Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  HT=(*H2)-(*Tr);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nTridiagonalMatrix - UpperHessenbergMatrix" << endl;
  HT=(*Tr)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Tr; Tr=0;

  cout << "\nUpperHessenbergMatrix - SymmetricPositiveMatrix" << endl;
  Syp=OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  Sq=(*H2)-(*Syp);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricPositive Matrix - UpperHessenbergMatrix" << endl;
  Sq=(*Syp)-(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Syp; Syp=0;

  cout << "\nUpperHessenbergMatrix - SymmetricMatrix" << endl;
  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  Sq=(*H2)-(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix - UpperHessenbergMatrix" << endl;
  Sq=(*Sy)-(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  cout << "\nUpperHessenbergMatrix - OrthogonalMatrix" << endl;
  Orth=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  M2=(*H2)-(*Orth);
  M2->printOn(cout);
  delete M2; M2=0;

  cout << "\nOrthogonalMatrix - UpperHessenbergMatrix" << endl;
  M2=(*Orth)-(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Orth; Orth=0;

  cout << "\nUpperHessenbergMatrix - UnitUpperTriangularMatrix" << endl;
  Uutg=OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  HT=(*H2)-(*Uutg);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUnitUpperTriangularMatrix - UpperHessenbergMatrix" << endl;
  HT=(*Uutg)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Uutg; Uutg=0;

  cout << "\nUpperHessenbergMatrix - UnitUpperTrapezoidalMatrix" << endl;
  Uutp=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  HT=(*H2)-(*Uutp);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUnitUpperTrapezoidalMatrix - UpperHessenbergMatrix" << endl;
  HT=(*Uutp)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Uutp; Uutp=0;

  cout << "\nUpperHessenbergMatrix - UpperTriangularMatrix" << endl;
  Utg=OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  HT=(*H2)-(*Utg);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUpperTriangularMatrix - UpperHessenbergMatrix" << endl;
  HT=(*Utg)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Utg; Utg=0;

  cout << "\nUpperHessenbergMatrix - UpperTrapezoidalMatrix" << endl;
  Utp=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  HT=(*H2)-(*Utp);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nUpperTrapezoidalMatrix - UpperHessenbergMatrix" << endl;
  HT=(*Utp)-(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete Utp; Utp=0;

  cout << "\nUpperHessenbergMatrix - UnitLowerTriangularMatrix" << endl;
  Ultg=OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  Sq=(*H2)-(*Ultg);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitLowerTriangularMatrix - UpperHessenbergMatrix" << endl;
  Sq=(*Ultg)-(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Ultg; Ultg=0;

  cout << "\nUpperHessenbergMatrix - UnitLowerTrapezoidalMatrix" << endl;
  Ultp=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*H2)-(*Ultp);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nUnitLowerTrapezoidalMatrix - UpperHessenbergMatrix" << endl;
  Sq=(*Ultp)-(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Ultp; Ultp=0;

  cout << "\nUpperHessenbergMatrix - LowerTriangularMatrix" << endl;
  Ltg=OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  Sq=(*H2)-(*Ltg);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTriangularMatrix - UpperHessenbergMatrix" << endl;
  Sq=(*Ltg)-(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Ltg; Ltg=0;

  cout << "\nUpperHessenbergMatrix - LowerTrapezoidalMatrix" << endl;
  Ltp=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  Sq=(*H2)-(*Ltp);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nLowerTrapezoidalMatrix - UpperHessenbergMatrix" << endl;
  Sq=(*Ltp)-(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Ltp; Ltp=0;

  cout << "\nUpperHessenbergMatrix - SquareMatrix" << endl;
  SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  Sq=(*H2)-(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix - UpperHessenbergMatrix" << endl;
  Sq=(*SQ)-(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  cout << "\nUpperHessenbergMatrix - Matrix" << endl;
  M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Sq=(*H2)-(*M);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nMatrix - UpperHessenbergMatrix" << endl;
  Sq=(*M)-(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete M; M=0;

  cout << "\nH3*scalar :" << endl;
  HT=(*H3)*scalar;
  HT->printOn(cout);
  delete HT; HT=0;
  
  cout << "\nH3/scalar.:" << endl;
  HT=(*H3)/scalar;
  HT->printOn(cout);
  delete HT; HT=0;
  delete H3; H3=0;

  UpperHessenbergMatrix<F,Z> *H=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * UpperHessenbergMatrix" << endl;
  SquareMatrix<F,Z> *S2=(*H2)*(*H);
  S2->printOn(cout);
  delete S2; S2=0;
  delete H; H=0;

  D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * DiagonalMatrix" << endl;
  HT=(*H2)*(*D);
  HT->printOn(cout);
  delete HT; HT=0;

  cout << "\nDiagonalMatrix * UpperHessenbergMatrix" << endl;
  HT=(*D)*(*H2);
  HT->printOn(cout);
  delete HT; HT=0;
  delete D; D=0;

  Sptr=OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * SymmetricPositiveTridiagonalMatrix"
       << endl;
  Sq=(*H2)*(*Sptr);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix * UpperHessenbergMatrix"
       << endl;
  Sq=(*Sptr)*(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sptr; Sptr=0;

  Str=OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * SymmetricTridiagonalMatrix" << endl;
  Sq=(*H2)*(*Str);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricTridiagonalMatrix * UpperHessenbergMatrix" << endl;
  Sq=(*Str)*(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Str; Str=0;

  Tr=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * TridiagonalMatrix" << endl;
  Sq=(*H2)*(*Tr);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nTridiagonalMatrix * UpperHessenbergMatrix" << endl;
  Sq=(*Tr)*(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Tr; Tr=0;

  Syp=OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * SymmetricPositiveMatrix" << endl;
  Sq=(*H2)*(*Syp);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricPositiveMatrix * UpperHessenbergMatrix" << endl;
  Sq=(*Syp)*(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Syp; Syp=0;

  Sy=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * SymmetricMatrix" << endl;
  Sq=(*H2)*(*Sy);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix * UpperHessenbergMatrix" << endl;
  Sq=(*Sy)*(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete Sy; Sy=0;

  Orth=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  cout << "\nUpperHessenbergMatrix * OrthogonalMatrix" << endl;
  M2=(*H2)*(*Orth);
  M2->printOn(cout);
  delete M2; M2=0;

  cout << "\nOrthogonalMatrix * UpperHessenbergMatrix" << endl;
  M2=(*Orth)*(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Orth; Orth=0;

  Uutg=OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * UnitUpperTriangularMatrix" << endl;
  M2=(*H2)*(*Uutg);
  M2->printOn(cout);
  delete M2; M2=0;

  cout << "\nUnitUpperTriangularMatrix * UpperHessenbergMatrix" << endl;
  M2=(*Uutg)*(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Uutg; Uutg=0;

  Uutp=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n+1,scalar);
  cout << "\nUpperHessenbergMatrix * UnitUpperTrapezoidalMatrix" << endl;
  M2=(*H2)*(*Uutp);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Uutp; Uutp=0;

  cout << "\nUnitUpperTrapezoidalMatrix * UpperHessenbergMatrix" << endl;
  Uutp=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n-1,n,scalar);
  M2=(*Uutp)*(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Uutp; Uutp=0;

  Utp=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n+1,scalar);
  cout << "\nUpperHessenbergMatrix * UpperTrapezoidalMatrix" << endl;
  M2=(*H2)*(*Utp);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Utp; Utp=0;

  Utp=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n-1,n,scalar);
  cout << "\nUpperTrapezoidalMatrix * UpperHessenbergMatrix" << endl;
  M2=(*Utp)*(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Utp; Utp=0;

  Ultg= OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * UnitLowerTriangularMatrix" << endl;
  M2=(*H2)*(*Ultg);
  M2->printOn(cout);
  delete M2; M2=0;

  cout << "\nUnitLowerTriangularMatrix * UpperHessenbergMatrix" << endl;
  M2=(*Ultg)*(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Ultg; Ultg=0;

  Ultp=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n-1,scalar);
  cout << "\nUpperHessenbergMatrix * UnitLowerTrapezoidalMatrix" << endl;
  M2=(*H2)*(*Ultp);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Ultp; Ultp=0;

  Ultp=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n+1,n,scalar);
  cout << "\nUnitLowerTrapezoidalMatrix * UpperHessenbergMatrix" << endl;
  M2=(*Ultp)*(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Ultp; Ultp=0;

  Ltp=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n-1,scalar);
  cout << "\nUpperHessenbergMatrix * LowerTrapezoidalMatrix" << endl;
  M2=(*H2)*(*Ltp);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Ltp; Ltp=0;

  Ltp=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n+1,n,scalar);
  cout << "\nLowerTrapezoidalMatrix * UpperHessenbergMatrix" << endl;
  M2=(*Ltp)*(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete Ltp; Ltp=0;

  SQ=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * SquareMatrix" << endl;
  Sq=(*H2)*(*SQ);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSquareMatrix * UpperHessenbergMatrix" << endl;
  Sq=(*SQ)*(*H2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete SQ; SQ=0;

  M=OPERATOR_NEW Matrix<F,Z>(n,n-1,scalar);
  cout << "\nUpperHessenbergMatrix * Matrix" << endl;
  M2=(*H2)*(*M);
  M2->printOn(cout);
  delete M2; M2=0;
  delete M; M=0;

  M=OPERATOR_NEW Matrix<F,Z>(n+1,n,scalar);
  cout << "\nMatrix * UpperHessenbergMatrix" << endl;
  M2=(*M)*(*H2);
  M2->printOn(cout);
  delete M2; M2=0;
  delete M; M=0;

  Vector<F,Z> *v=OPERATOR_NEW Vector<F,Z>(n,scalar);
  cout << "\nUpperHessenbergMatrix * Vector" << endl;
  Vector<F,Z> *w=(*H2)*(*v);
  w->printOn(cout);
  delete w; w=0;
  delete v; v=0;

//cout << "\ntranspose" << endl;
//SquareMatrix<F,Z> *S=H2->transpose();
//S->printOn(cout);
//delete S; S=0;

  cout << "\ncopyFrom" << endl;
  H1=OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,Vector<F,Z>::one_);
  H1->copyFrom(n-1,*H2);
  H1->printOn(cout);
  delete H1; H1=0;

  for (int j=0;j<n;j++) {
    for (int i=0;i<=min(n-1,j+1);i++) (*H2)(i,j)=F(i+j)*scalar;
  }
  cout <<"\nH2 = " << endl;
  H2->printOn(cout);
  cout << "\nnormFrobenius = " << H2->normFrobenius() << endl;
  cout << "normInfinity = " << H2->normInfinity() << endl;
  cout << "normMaxEntry = " << H2->normMaxEntry() << endl;
  cout << "normOne = " << H2->normOne() << endl;

  SquareMatrix<F,complex<F> > *V=
    OPERATOR_NEW SquareMatrix<F,complex<F> >(n);
  SquareMatrix<F,complex<F> > *U=
    OPERATOR_NEW SquareMatrix<F,complex<F> >(n);
  Vector<F,complex<F> > *lambda=H2->eigenvalues(V,U);
  cout << "\nafter eigenvalues" << endl;
  lambda->printOn(cout);
  V->printOn(cout);
  U->printOn(cout);
  SquareMatrix<F,complex<F> > *R=
    OPERATOR_NEW SquareMatrix<F,complex<F> >(n);
  for (int j=0;j<n;j++) {
    for (int i=0;i<n;i++) {
      (*R)(i,j)=-(*U)(i,j)*(*lambda)[j];
      for (int k=max(0,i-1);k<n;k++) {
        (*R)(i,j)+=(*H2)(i,k)*(*U)(k,j);
      }
    }
  }
  cout << "\nH * U - U * lambda" << endl;
  R->printOn(cout);
  for (int j=0;j<n;j++) {
    for (int i=0;i<n;i++) {
      (*R)(i,j)=-(*lambda)[i]*conj((*V)(j,i));
      for (int k=0;k<=min(n-1,j+1);k++) {
        (*R)(i,j)+=conj((*V)(k,i))*(*H2)(k,j);
      }
    }
  }
  cout << "\nV * H - lambda * V" << endl;
  R->printOn(cout);
  delete R; R=0;
  delete U; U=0;
  delete V; V=0;
  delete lambda; lambda=0;

  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(4);
  (*b)[0]=static_cast<F>(12.)*scalar;
  (*b)[1]=static_cast<F>(20.)*scalar;
  (*b)[2]=static_cast<F>(24.)*scalar;
  (*b)[3]=static_cast<F>(22.)*scalar;
  cout << "\nsquare system of linear equations" << endl;
  H2->printOn(cout);
  b->printOn(cout);
  cout << "solve('N'):" << endl;
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(4);
  H2->solve(*b,*x);
  x->printOn(cout);
  Vector<F,Z> *r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  H2->uhmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'N');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;
  delete b; b=0;

  b=OPERATOR_NEW Vector<F,Z>(4);
  (*b)[0]=static_cast<F>(2.)*scalar;
  (*b)[1]=static_cast<F>(12.)*scalar;
  (*b)[2]=static_cast<F>(28.)*scalar;
  (*b)[3]=static_cast<F>(36.)*scalar;
  cout << "\nsquare system of linear equations" << endl;
  H2->printOn(cout);
  b->printOn(cout);
  cout << "solve('T'):" << endl;
  x=OPERATOR_NEW Vector<F,Z>(4);
  H2->solve(*b,*x,'T');
  x->printOn(cout);
  r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  H2->uhmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'T');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;
  delete b; b=0;

  if (sizeof(F)!=sizeof(Z)) {
    b=OPERATOR_NEW Vector<F,Z>(4);
    (*b)[0]=static_cast<F>(2.)*scalar;
    (*b)[1]=static_cast<F>(12.)*scalar;
    (*b)[2]=static_cast<F>(28.)*scalar;
    (*b)[3]=static_cast<F>(36.)*scalar;
    cout << "\nsquare system of linear equations" << endl;
    H2->printOn(cout);
    b->printOn(cout);
    cout << "solve('C'):" << endl;
    x=OPERATOR_NEW Vector<F,Z>(4);
    H2->solve(*b,*x,'C');
    x->printOn(cout);
    r=OPERATOR_NEW Vector<F,Z>(b->size());
    r->copy(*b);
    H2->uhmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'C');
    cout << "r = " << endl;
    r->printOn(cout);
    delete r; r=0;
    delete x; x=0;
    delete b; b=0;
  }
  delete H2; H2=0;

  UpperHessenbergMatrix<F,Z> *HH=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(3,3);
  (*HH)(0,0)=static_cast<F>(2.)*scalar;
    (*HH)(0,1)=static_cast<F>(-1.)*scalar;
    (*HH)(0,2)=static_cast<F>(0.)*scalar;
  (*HH)(1,0)=static_cast<F>(-1.)*scalar;
    (*HH)(1,1)=static_cast<F>(2.)*scalar;
    (*HH)(1,2)=static_cast<F>(-1.)*scalar;
  (*HH)(2,1)=static_cast<F>(-1.)*scalar;
    (*HH)(2,2)=static_cast<F>(2.)*scalar;
  b=OPERATOR_NEW Vector<F,Z>(3);
  (*b)[0]=static_cast<F>(-6.)*scalar;
  (*b)[1]=static_cast<F>(8.)*scalar;
  (*b)[2]=static_cast<F>(-6.)*scalar;
  cout << "\nsquare system of linear equations" << endl;
  HH->printOn(cout);
  b->printOn(cout);
  cout << "solve('N'):" << endl;
  x=OPERATOR_NEW Vector<F,Z>(3);
  HH->solve(*b,*x);
  x->printOn(cout);
  r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  HH->uhmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'N');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;

  cout << "solve('T'):" << endl;
  x=OPERATOR_NEW Vector<F,Z>(3);
  HH->solve(*b,*x,'T');
  x->printOn(cout);
  r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  HH->uhmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'T');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "solve('C'):" << endl;
    x=OPERATOR_NEW Vector<F,Z>(3);
    HH->solve(*b,*x,'C');
    x->printOn(cout);
    r=OPERATOR_NEW Vector<F,Z>(b->size());
    r->copy(*b);
    HH->uhmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'C');
    cout << "r = " << endl;
    r->printOn(cout);
    delete r; r=0;
    delete x; x=0;
  }
  delete b; b=0;
 
  Matrix<F,Z> *B=OPERATOR_NEW Matrix<F,Z>(3,1);
  (*B)(0,0)=static_cast<F>(-6.)*scalar;
  (*B)(1,0)=static_cast<F>(8.)*scalar;
  (*B)(2,0)=static_cast<F>(-6.)*scalar;
  cout << "\nsolve('L','N'):" << endl;
  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(3,1);
  HH->solve(*B,*X);
  X->printOn(cout);
  Matrix<F,Z> *RR=OPERATOR_NEW Matrix<F,Z>(B->size(0),B->size(1));
  RR->copy(*B);
  HH->uhmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*RR,'L','N');
  cout << "RR = " << endl;
  RR->printOn(cout);
  delete RR; RR=0;
  delete X; X=0;

  cout << "\nsolve('L','T'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(3,1);
  HH->solve(*B,*X,'L','T');
  X->printOn(cout);
  RR=OPERATOR_NEW Matrix<F,Z>(B->size(0),B->size(1));
  RR->copy(*B);
  HH->uhmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*RR,'L','T');
  cout << "RR = " << endl;
  RR->printOn(cout);
  delete RR; RR=0;
  delete X; X=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "\nsolve('L','C'):" << endl;
    X=OPERATOR_NEW Matrix<F,Z>(3,1);
    HH->solve(*B,*X,'L','C');
    X->printOn(cout);
    RR=OPERATOR_NEW Matrix<F,Z>(B->size(0),B->size(1));
    RR->copy(*B);
    HH->uhmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*RR,'L','C');
    cout << "RR = " << endl;
    RR->printOn(cout);
    delete RR; RR=0;
    delete X; X=0;
  }
  delete B; B=0;

  B=OPERATOR_NEW Matrix<F,Z>(1,3);
  (*B)(0,0)=static_cast<F>(-6.)*scalar;
  (*B)(0,1)=static_cast<F>(8.)*scalar;
  (*B)(0,2)=static_cast<F>(-6.)*scalar;
  cout << "\nsolve('R','N'):" << endl;
  HH->printOn(cout);
  B->printOn(cout);
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  HH->solve(*B,*X,'R','N');
  X->printOn(cout);
  RR=OPERATOR_NEW Matrix<F,Z>(B->size(0),B->size(1));
  RR->copy(*B);
  HH->uhmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*RR,'R','N');
  cout << "RR = " << endl;
  RR->printOn(cout);
  delete RR; RR=0;
  delete X; X=0;

  cout << "\nsolve('R','T'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  HH->solve(*B,*X,'R','T');
  X->printOn(cout);
  RR=OPERATOR_NEW Matrix<F,Z>(B->size(0),B->size(1));
  RR->copy(*B);
  HH->uhmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*RR,'R','T');
  cout << "RR = " << endl;
  RR->printOn(cout);
  delete RR; RR=0;
  delete X; X=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "\nsolve('R','C'):" << endl;
    X=OPERATOR_NEW Matrix<F,Z>(1,3);
    HH->solve(*B,*X,'R','C');
    X->printOn(cout);
    RR=OPERATOR_NEW Matrix<F,Z>(B->size(0),B->size(1));
    RR->copy(*B);
    HH->uhmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*RR,'R','C');
    cout << "RR = " << endl;
    RR->printOn(cout);
    delete RR; RR=0;
    delete X; X=0;
  }
  delete B; B=0;
  delete HH; HH=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void BandMatrix<F,Z>::fillWith(Z d) {
  Z *colj=addr();
  for (int j=0;j<dim;j++,colj+=nt) {
    for (int i=max(0,nsup-j);i<min(dim+nsup-j,nt);i++) colj[i]=d;
  }
}

template<typename F,typename Z> void BandMatrix<F,Z>::printOn(
ostream& s) const {
  s << "BandMatrix(" << size(0) << " x " << size(1) << ")\n" ;
  s << "nsub = " << nsub << ", nsup = " << nsup << ", nt = " << nt << endl;
  for (int i=0; i< size(0); i++) {
    for (int j=0; j< size(1); j++) s << operator()(i,j) << "  ";
    s << "\n";
  }
//if (AB) {
//  s << "AB(" << AB->size(0) << " x " << AB->size(1) << ")" << endl;
//  s << "AB:" << endl;
//  AB->printOn(s);
//} else s << "AB = 0" << endl;
}

template<typename F,typename Z> void testBandMatrix(F fscalar,Z scalar) {
  BandMatrix<F,Z> *A1=OPERATOR_NEW BandMatrix<F,Z>;
  cout << "\nafter BandMatrix()" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);

  int n=4,nsub=2,nsup=1;
  cout << "n, nsub, nsup = " << n << " " << nsub << " " << nsup  <<endl;
  BandMatrix<F,Z> *A2=OPERATOR_NEW BandMatrix<F,Z>(n,nsub,nsup);;
  cout << "\nafter BandMatrix(n,nsub,nsup)" << endl;
  cout << "A2:" << endl;
  A2->printOn(cout);

  BandMatrix<F,Z> *A3=OPERATOR_NEW BandMatrix<F,Z>(n,nsub,nsup,scalar);
  cout << "\nafter A3=BandMatrix(n,nsub,nsup,scalar)" << endl;
  cout << "A3:" << endl;
  A3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=A3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }
  delete A3; A3=0;

  *A2=scalar;
  cout << "\nafter A2=scalar" << endl;
  cout << "A2:" << endl;
  A2->printOn(cout);

  A1->resize(n,nsub,nsup);
  cout << "\nafter A1->resize(n,nsub,nsup)" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);
  delete A1; A1=0;

  A1=OPERATOR_NEW BandMatrix<F,Z>();
  A1->resize(*A2);
  cout << "\nafter A1->resize(*A2)" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);

  A1->copy(*A2);
  cout << "\nafter A1->copy(A2)" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);

  cout << "\nA1+A2:" << endl;
  BandMatrix<F,Z> *T=(*A1)+(*A2);
  T->printOn(cout);
  delete T; T=0;

  {
  cout << "\nBandMatrix + UpperHessenbergMatrix" << endl;
  UpperHessenbergMatrix<F,Z> *H=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*H);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperHessenbergMatrix + BandMatrix" << endl;
  S=(*H)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete H; H=0;
  }

  {
  cout << "\nBandMatrix + DiagonalMatrix" << endl;
  DiagonalMatrix<F,Z> *D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)+(*D);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nDiagonalMatrix + BandMatrix" << endl;
  B=(*D)+(*A1);
  B->printOn(cout);
  delete B; B=0;
  delete D; D=0;
  }

  {
  cout << "\nBandMatrix + SymmetricPositiveTridiagonalMatrix" << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *SPT=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)+(*SPT);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + BandMatrix" << endl;
  B=(*SPT)+(*A1);
  B->printOn(cout);
  delete B; B=0;
  delete SPT; SPT=0;
  }

  {
  cout << "\nBandMatrix + SymmetricTridiagonalMatrix" << endl;
  SymmetricTridiagonalMatrix<F,Z> *ST=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)+(*ST);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricTridiagonalMatrix + BandMatrix" << endl;
  B=(*ST)+(*A1);
  B->printOn(cout);
  delete B; B=0;
  delete ST; ST=0;
  }

  {
  cout << "\nBandMatrix + TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *TM=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)+(*TM);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nTridiagonalMatrix + BandMatrix" << endl;
  B=(*TM)+(*A1);
  B->printOn(cout);
  delete B; B=0;
  delete TM; TM=0;
  }

  {
  cout << "\nBandMatrix + SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *SP=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*SP);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSymmetricPositiveMatrix + BandMatrix" << endl;
  S=(*SP)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete SP; SP=0;
  }

  {
  cout << "\nBandMatrix + SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *S=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *Sq=(*A1)+(*S);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix + BandMatrix" << endl;
  Sq=(*S)+(*A1);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete S; S=0;
  }

  {
  cout << "\nBandMatrix + OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  SquareMatrix<F,Z> *S=(*A1)+(*Q);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nOrthogonalMatrix + BandMatrix" << endl;
  S=(*Q)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete Q; Q=0;
  }

  {
  cout << "\nBandMatrix + UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTriangularMatrix + BandMatrix" << endl;
  S=(*U)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix + UnitUpperTrapezoidalMatrix" << endl;
  UnitUpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTrapezoidalMatrix + BandMatrix" << endl;
  S=(*U)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix + UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTriangularMatrix + BandMatrix" << endl;
  S=(*U)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix + UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTrapezoidalMatrix + BandMatrix" << endl;
  S=(*U)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix + UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTriangularMatrix + BandMatrix" << endl;
  S=(*L)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix + UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTrapezoidalMatrix + BandMatrix" << endl;
  S=(*L)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix + LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTriangularMatrix + BandMatrix" << endl;
  S=(*L)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix + LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTrapezoidalMatrix + BandMatrix" << endl;
  S=(*L)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix + SquareMatrix" << endl;
  SquareMatrix<F,Z> *SM=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*SM);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSquareMatrix + BandMatrix" << endl;
  S=(*SM)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete SM; SM=0;
  }

  {
  cout << "\nBandMatrix + Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*M);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nMatrix + BandMatrix" << endl;
  S=(*M)+(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete M; M=0;
  }

  cout << "\nA1-A2:" << endl;
  T=(*A1)-(*A2);
  T->printOn(cout);
  delete T; T=0;

  {
  cout << "\nBandMatrix - UpperHessenbergMatrix" << endl;
  UpperHessenbergMatrix<F,Z> *H=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*H);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperHessenbergMatrix - BandMatrix" << endl;
  S=(*H)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete H; H=0;
  }

  {
  cout << "\nBandMatrix - DiagonalMatrix" << endl;
  DiagonalMatrix<F,Z> *D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)-(*D);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nDiagonalMatrix - BandMatrix" << endl;
  B=(*D)-(*A1);
  B->printOn(cout);
  delete B; B=0;
  delete D; D=0;
  }

  {
  cout << "\nBandMatrix - SymmetricPositiveTridiagonalMatrix" << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *SPT=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)-(*SPT);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - BandMatrix" << endl;
  B=(*SPT)-(*A1);
  B->printOn(cout);
  delete B; B=0;
  delete SPT; SPT=0;
  }

  {
  cout << "\nBandMatrix - SymmetricTridiagonalMatrix" << endl;
  SymmetricTridiagonalMatrix<F,Z> *ST=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)-(*ST);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricTridiagonalMatrix - BandMatrix" << endl;
  B=(*ST)-(*A1);
  B->printOn(cout);
  delete B; B=0;
  delete ST; ST=0;
  }

  {
  cout << "\nBandMatrix - TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *TM=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)-(*TM);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nTridiagonalMatrix - BandMatrix" << endl;
  B=(*TM)-(*A1);
  B->printOn(cout);
  delete B; B=0;
  delete TM; TM=0;
  }

  {
  cout << "\nBandMatrix - SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *SP=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*SP);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSymmetricPositiveMatrix - BandMatrix" << endl;
  S=(*SP)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete SP; SP=0;
  }

  {
  cout << "\nBandMatrix - SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *S=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *Sq=(*A1)-(*S);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix - BandMatrix" << endl;
  Sq=(*S)-(*A1);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete S; S=0;
  }

  {
  cout << "\nBandMatrix - OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  SquareMatrix<F,Z> *S=(*A1)-(*Q);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nOrthogonalMatrix - BandMatrix" << endl;
  S=(*Q)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete Q; Q=0;
  }

  {
  cout << "\nBandMatrix - UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTriangularMatrix - BandMatrix" << endl;
  S=(*U)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix - UnitUpperTrapezoidalMatrix" << endl;
  UnitUpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTrapezoidalMatrix - BandMatrix" << endl;
  S=(*U)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix - UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTriangularMatrix - BandMatrix" << endl;
  S=(*U)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix - UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTrapezoidalMatrix - BandMatrix" << endl;
  S=(*U)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix - UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTriangularMatrix - BandMatrix" << endl;
  S=(*L)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix - UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTrapezoidalMatrix - BandMatrix" << endl;
  S=(*L)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix - LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTriangularMatrix - BandMatrix" << endl;
  S=(*L)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix - LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTrapezoidalMatrix - BandMatrix" << endl;
  S=(*L)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix - SquareMatrix" << endl;
  SquareMatrix<F,Z> *SM=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*SM);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSquareMatrix - BandMatrix" << endl;
  S=(*SM)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete SM; SM=0;
  }

  {
  cout << "\nBandMatrix - Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*M);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nMatrix - BandMatrix" << endl;
  S=(*M)-(*A1);
  S->printOn(cout);
  delete S; S=0;
  delete M; M=0;
  }

  cout << "\nA1*scalar :" << endl;
  T=(*A1)*scalar;
  T->printOn(cout);
  delete T; T=0;

  cout << "\nA1/scalar.:" << endl;
  T=(*A1)/scalar;
  T->printOn(cout);
  delete T; T=0;
  delete A1; A1=0;
  delete A2; A2=0;

  A1=OPERATOR_NEW BandMatrix<F,Z>(n,1,2,scalar);
  A2=OPERATOR_NEW BandMatrix<F,Z>(n,2,1,scalar);
  cout << "\nA1 * A2 = " << endl;
  T=(*A1)*(*A2);
  T->printOn(cout);
  delete T; T=0;

  cout << "\nA2 * A1 = " << endl;
  T=(*A2)*(*A1);
  T->printOn(cout);
  delete T; T=0;

  {
  cout << "\nBandMatrix * UpperHessenbergMatrix" << endl;
  UpperHessenbergMatrix<F,Z> *H=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)*(*H);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperHessenbergMatrix * BandMatrix" << endl;
  S=(*H)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete H; H=0;
  }

  {
  cout << "\nBandMatrix * DiagonalMatrix" << endl;
  DiagonalMatrix<F,Z> *D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)*(*D);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nDiagonalMatrix * BandMatrix" << endl;
  B=(*D)*(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete D; D=0;
  }

  {
  cout << "\nBandMatrix * SymmetricPositiveTridiagonalMatrix" << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *SPT=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)*(*SPT);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix * BandMatrix" << endl;
  B=(*SPT)*(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete SPT; SPT=0;
  }

  {
  cout << "\nBandMatrix * SymmetricTridiagonalMatrix" << endl;
  SymmetricTridiagonalMatrix<F,Z> *ST=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)*(*ST);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricTridiagonalMatrix * BandMatrix" << endl;
  B=(*ST)*(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete ST; ST=0;
  }

  {
  cout << "\nBandMatrix * TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *TM=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)*(*TM);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nTridiagonalMatrix * BandMatrix" << endl;
  B=(*TM)*(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete TM; TM=0;
  }

  {
  cout << "\nBandMatrix * SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *SP=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)*(*SP);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSymmetricPositiveMatrix * BandMatrix" << endl;
  S=(*SP)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SP; SP=0;
  }

  {
  cout << "\nBandMatrix * SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *S=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *Sq=(*A1)*(*S);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix * BandMatrix" << endl;
  Sq=(*S)*(*A2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete S; S=0;
  }

  {
  cout << "\nBandMatrix * OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  Matrix<F,Z> *M=(*A1)*(*Q);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nOrthogonalMatrix * BandMatrix" << endl;
  M=(*Q)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete Q; Q=0;
  }

  {
  cout << "\nBandMatrix * UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  Matrix<F,Z> *M=(*A1)*(*U);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nUnitUpperTriangularMatrix * BandMatrix" << endl;
  M=(*U)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix * UnitUpperTrapezoidalMatrix" << endl;
  UnitUpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n+1,scalar);
  Matrix<F,Z> *M=(*A1)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  cout << "\nUnitUpperTrapezoidalMatrix * BandMatrix" << endl;
  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n-1,n,scalar);
  M=(*U)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix * UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  Matrix<F,Z> *M=(*A1)*(*U);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nUpperTriangularMatrix * BandMatrix" << endl;
  M=(*U)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix * UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n+1,scalar);
  Matrix<F,Z> *M=(*A1)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  cout << "\nUpperTrapezoidalMatrix * BandMatrix" << endl;
  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n-1,n,scalar);
  M=(*U)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;
  }

  {
  cout << "\nBandMatrix * UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  Matrix<F,Z> *M=(*A1)*(*L);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nUnitLowerTriangularMatrix * BandMatrix" << endl;
  M=(*L)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix * UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n-1,scalar);
  Matrix<F,Z> *M=(*A1)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  cout << "\nUnitLowerTrapezoidalMatrix * BandMatrix" << endl;
  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n+1,n,scalar);
  M=(*L)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix * LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  Matrix<F,Z> *M=(*A1)*(*L);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nLowerTriangularMatrix * BandMatrix" << endl;
  M=(*L)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix * LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n-1,scalar);
  Matrix<F,Z> *M=(*A1)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  cout << "\nLowerTrapezoidalMatrix * BandMatrix" << endl;
  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n+1,n,scalar);
  M=(*L)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;
  }

  {
  cout << "\nBandMatrix * SquareMatrix" << endl;
  SquareMatrix<F,Z> *SM=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)*(*SM);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSquareMatrix * BandMatrix" << endl;
  S=(*SM)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SM; SM=0;
  }

  {
  cout << "\nBandMatrix * Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Matrix<F,Z> *S=(*A1)*(*M);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nMatrix * BandMatrix" << endl;
  S=(*M)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete M; M=0;
  }

  Vector<F,Z> *v=OPERATOR_NEW Vector<F,Z>(n,scalar);
  cout << "\nBandMatrix * Vector" << endl;
  Vector<F,Z> *w=(*A1)*(*v);
  w->printOn(cout);
  delete w; w=0;

  w=(*A2)*(*v);
  w->printOn(cout);
  delete w; w=0;
  delete v; v=0;
  delete A1; A1=0;
  delete A2; A2=0;

  BandMatrix<F,Z> *AA=OPERATOR_NEW BandMatrix<F,Z>(3,1,1);
  (*AA)(0,0)=static_cast<F>(2.)*scalar;
    (*AA)(0,1)=static_cast<F>(-3.)*scalar;
  (*AA)(1,0)=static_cast<F>(-1.)*scalar;
    (*AA)(1,1)=static_cast<F>(2.)*scalar;
    (*AA)(1,2)=static_cast<F>(-3.)*scalar;
  (*AA)(2,1)=static_cast<F>(-1.)*scalar;
    (*AA)(2,2)=static_cast<F>(2.)*scalar;
  cout << "\nequilibrate" << endl;
  AA->printOn(cout);
  Vector<F,F> *rr=OPERATOR_NEW Vector<F,F>(3);
  F rowcnd,colcnd;
  Vector<F,F> *c=OPERATOR_NEW Vector<F,F>(3);
  F anorm=AA->equilibrate(*rr,*c,rowcnd,colcnd);
  cout << "rr = " << endl;
  rr->printOn(cout);
  cout << "c = " << endl;
  c->printOn(cout);
  cout << "rowcnd,colcnd = " << rowcnd << " " << colcnd << endl;
  delete rr; rr=0;
  delete c; c=0;

  cout << "\nnormFrobenius = " << AA->normFrobenius() << endl;
  cout << "normInfinity = " << AA->normInfinity() << endl;
  cout << "normMaxEntry = " << AA->normMaxEntry() << endl;
  cout << "normOne = " << AA->normOne() << endl;
  cout << "reciprocalConditionNumber('I') = "
       << AA->reciprocalConditionNumber('I') << endl;
  cout << "reciprocalConditionNumber('O') = "
       << AA->reciprocalConditionNumber('O') << endl;

//cout << "\ntranspose" << endl;
//T=AA->transpose();
//T->printOn(cout);
//delete T; T=0;

  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(3);
  (*b)[0]=0.;
  (*b)[1]=6.;
  (*b)[2]=0.;
  cout << "\nsquare system of linear equations" << endl;
  cout << "A:" << endl;
  AA->printOn(cout);
  cout << "b:" << endl;
  b->printOn(cout);
  cout << "solve('N'):" << endl;
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(3);
  AA->solve(*b,*x);
  cout << "x = " << endl;
  x->printOn(cout);
  Vector<F,Z> *r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  AA->gbmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'N');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;

  cout << "\nsolve('T')" << endl;
  AA->solve(*b,*x,'T');
  x->printOn(cout);
  r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  AA->gbmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'T');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "\nsolve('C')" << endl;
    AA->solve(*b,*x,'C');
    x->printOn(cout);
    r=OPERATOR_NEW Vector<F,Z>(b->size());
    r->copy(*b);
    AA->gbmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'C');
    cout << "r = " << endl;
    r->printOn(cout);
    delete r; r=0;
  }
  delete x; x=0;
  delete b; b=0;

  Matrix<F,Z> *B=OPERATOR_NEW Matrix<F,Z>(3,1);
  (*B)(0,0)=0.;
  (*B)(1,0)=6.;
  (*B)(2,0)=0.;
  cout << "\nsolve('L','N'):" << endl;
  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(3,1);
  AA->solve(*B,*X);
  X->printOn(cout);
  Matrix<F,Z> *R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*B);
  AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L','N');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  cout << "\nsolve('L','T'):" << endl;
  AA->solve(*B,*X,'L','T');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*B);
  AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L','T');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "\nsolve('L','C'):" << endl;
    AA->solve(*B,*X,'L','C');
    X->printOn(cout);
    R=OPERATOR_NEW Matrix<F,Z>(3,1);
    R->copy(*B);
    AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L','C');
    cout << "R = " << endl;
    R->printOn(cout);
    delete R; R=0;
  }
  delete X; X=0;
  delete B; B=0;

  B=OPERATOR_NEW Matrix<F,Z>(1,3);
  (*B)(0,0)=0.;
  (*B)(0,1)=6.;
  (*B)(0,2)=0.;
  cout << "\nsolve('R','N'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  AA->solve(*B,*X,'R','N');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*B);
  AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R','N');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  cout << "\nsolve('R','T'):" << endl;
  AA->solve(*B,*X,'R','T');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*B);
  AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R','T');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "\nsolve('R','C'):" << endl;
    AA->solve(*B,*X,'R','C');
    X->printOn(cout);
    R=OPERATOR_NEW Matrix<F,Z>(1,3);
    R->copy(*B);
    AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R','C');
    cout << "R = " << endl;
    R->printOn(cout);
    delete R; R=0;
  }
  delete X; X=0;
  delete B; B=0;
  delete AA; AA=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void SymmetricBandMatrix<F,Z>::printOn(
ostream& s) const {
  s << "SymmetricBandMatrix(" << size(0) << " x " << size(1) << ")\n" ;
//s << "nsub = " << nsub << ", nsup = " << nsup << ", nt = " << nt << endl;
  for (int i=0; i< size(0); i++) {
    for (int j=0; j< size(1); j++) s << operator()(i,j) << "  ";
    s << "\n";
  }
//if (AB) {
//  s << "AB(" << AB->size(0) << " x " << AB->size(1) << ")" << endl;
//  s << "AB:" << endl;
//  AB->printOn(s);
//} else s << "AB = 0" << endl;
}

template<typename F,typename Z> void testSymmetricBandMatrix(F fscalar,
Z scalar) {
  SymmetricBandMatrix<F,Z> *A1=OPERATOR_NEW SymmetricBandMatrix<F,Z>;
  cout << "\nafter SymmetricBandMatrix()" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);

  int n=4,nsub=2;
  cout << "n, nsub = " << n << " " << nsub << endl;
  SymmetricBandMatrix<F,Z> *A2=
    OPERATOR_NEW SymmetricBandMatrix<F,Z>(n,nsub);;
  cout << "\nafter SymmetricBandMatrix(n,nsub)" << endl;
  cout << "A2:" << endl;
  A2->printOn(cout);

  SymmetricBandMatrix<F,Z> *A3=
    OPERATOR_NEW SymmetricBandMatrix<F,Z>(n,nsub,scalar);
  cout << "\nafter A3=SymmetricBandMatrix(n,nsub,scalar)" << endl;
  cout << "A3:" << endl;
  A3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=A3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }

  {
  cout << "\nmakeBandMatrix" << endl;
  BandMatrix<F,Z> *B=A3->makeBandMatrix();
  B->printOn(cout);
  delete B; B=0;
  }

  {
  cout << "\nmakeSymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *S=A3->makeSymmetricMatrix();
  S->printOn(cout);
  delete S; S=0;
  }
  delete A3; A3=0;

  *A2=scalar;
  cout << "\nafter A2=scalar" << endl;
  cout << "A2:" << endl;
  A2->printOn(cout);

  A1->resize(n,1);
  cout << "\nafter A1->resize(n,nsub)" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);
  delete A1; A1=0;

  A1=OPERATOR_NEW SymmetricBandMatrix<F,Z>();
  A1->resize(*A2);
  cout << "\nafter A1->resize(*A2)" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);

  A1->copy(*A2);
  cout << "\nafter A1->copy(A2)" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);
  delete A1; A1=0;

  cout << "\nA1+A2:" << endl;
  A1=OPERATOR_NEW SymmetricBandMatrix<F,Z>(n,1,scalar);;
  SymmetricBandMatrix<F,Z> *T=(*A1)+(*A2);
  T->printOn(cout);
  delete T; T=0;

  {
  cout << "\nSymmetricBandMatrix + BandMatrix" << endl;
  BandMatrix<F,Z> *B=OPERATOR_NEW BandMatrix<F,Z>(n,1,2,scalar);
  BandMatrix<F,Z> *S=(*A1)+(*B);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nBandMatrix + SymmetricBandMatrix" << endl;
  S=(*B)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete B; B=0;
  }

  {
  cout << "\nSymmetricBandMatrix + UpperHessenbergMatrix" << endl;
  UpperHessenbergMatrix<F,Z> *H=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*H);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperHessenbergMatrix + SymmetricBandMatrix" << endl;
  S=(*H)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete H; H=0;
  }

  {
  cout << "\nSymmetricBandMatrix + DiagonalMatrix" << endl;
  DiagonalMatrix<F,Z> *D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  SymmetricBandMatrix<F,Z> *B=(*A1)+(*D);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nDiagonalMatrix + SymmetricBandMatrix" << endl;
  B=(*D)+(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete D; D=0;
  }

  {
  cout << "\nSymmetricBandMatrix + SymmetricPositiveTridiagonalMatrix"
       << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *SPT=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  SymmetricBandMatrix<F,Z> *B=(*A1)+(*SPT);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix + SymmetricBandMatrix"
       << endl;
  B=(*SPT)+(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete SPT; SPT=0;
  }

  {
  cout << "\nSymmetricBandMatrix + SymmetricTridiagonalMatrix" << endl;
  SymmetricTridiagonalMatrix<F,Z> *ST=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  SymmetricBandMatrix<F,Z> *B=(*A1)+(*ST);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricTridiagonalMatrix + SymmetricBandMatrix" << endl;
  B=(*ST)+(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete ST; ST=0;
  }

  {
  cout << "\nSymmetricBandMatrix + TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *TM=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)+(*TM);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nTridiagonalMatrix + SymmetricBandMatrix" << endl;
  B=(*TM)+(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete TM; TM=0;
  }

  {
  cout << "\nSymmetricBandMatrix + SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *SP=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *S=(*A1)+(*SP);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSymmetricPositiveMatrix + SymmetricBandMatrix" << endl;
  S=(*SP)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SP; SP=0;
  }

  {
  cout << "\nSymmetricBandMatrix + SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *S=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *Sq=(*A1)+(*S);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix + SymmetricBandMatrix" << endl;
  Sq=(*S)+(*A2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete S; S=0;
  }

  {
  cout << "\nSymmetricBandMatrix + OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  SquareMatrix<F,Z> *S=(*A1)+(*Q);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nOrthogonalMatrix + SymmetricBandMatrix" << endl;
  S=(*Q)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete Q; Q=0;
  }

  {
  cout << "\nSymmetricBandMatrix + UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTriangularMatrix + SymmetricBandMatrix" << endl;
  S=(*U)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix + UnitUpperTrapezoidalMatrix" << endl;
  UnitUpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTrapezoidalMatrix + SymmetricBandMatrix" << endl;
  S=(*U)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix + UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTriangularMatrix + SymmetricBandMatrix" << endl;
  S=(*U)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix + UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTrapezoidalMatrix + SymmetricBandMatrix" << endl;
  S=(*U)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix + UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTriangularMatrix + SymmetricBandMatrix" << endl;
  S=(*L)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix + UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTrapezoidalMatrix + SymmetricBandMatrix" << endl;
  S=(*L)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix + LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTriangularMatrix + SymmetricBandMatrix" << endl;
  S=(*L)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix + LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTrapezoidalMatrix + SymmetricBandMatrix" << endl;
  S=(*L)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix + SquareMatrix" << endl;
  SquareMatrix<F,Z> *SM=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*SM);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSquareMatrix + SymmetricBandMatrix" << endl;
  S=(*SM)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SM; SM=0;
  }

  {
  cout << "\nSymmetricBandMatrix + Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*M);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nMatrix + SymmetricBandMatrix" << endl;
  S=(*M)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete M; M=0;
  }

  cout << "\nA1-A2:" << endl;
  T=(*A1)-(*A2);
  T->printOn(cout);
  delete T; T=0;

  {
  cout << "\nSymmetricBandMatrix - BandMatrix" << endl;
  BandMatrix<F,Z> *B=OPERATOR_NEW BandMatrix<F,Z>(n,1,2,scalar);
  BandMatrix<F,Z> *S=(*A1)-(*B);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nBandMatrix - SymmetricBandMatrix" << endl;
  S=(*B)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete B; B=0;
  }

  {
  cout << "\nSymmetricBandMatrix - UpperHessenbergMatrix" << endl;
  UpperHessenbergMatrix<F,Z> *H=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*H);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperHessenbergMatrix - SymmetricBandMatrix" << endl;
  S=(*H)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete H; H=0;
  }

  {
  cout << "\nSymmetricBandMatrix - DiagonalMatrix" << endl;
  DiagonalMatrix<F,Z> *D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  SymmetricBandMatrix<F,Z> *B=(*A1)-(*D);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nDiagonalMatrix - SymmetricBandMatrix" << endl;
  B=(*D)-(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete D; D=0;
  }

  {
  cout << "\nSymmetricBandMatrix - SymmetricPositiveTridiagonalMatrix"
       << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *SPT=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  SymmetricBandMatrix<F,Z> *B=(*A1)-(*SPT);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix - SymmetricBandMatrix"
       << endl;
  B=(*SPT)-(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete SPT; SPT=0;
  }

  {
  cout << "\nSymmetricBandMatrix - SymmetricTridiagonalMatrix" << endl;
  SymmetricTridiagonalMatrix<F,Z> *ST=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  SymmetricBandMatrix<F,Z> *B=(*A1)-(*ST);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricTridiagonalMatrix - SymmetricBandMatrix" << endl;
  B=(*ST)-(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete ST; ST=0;
  }

  {
  cout << "\nSymmetricBandMatrix - TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *TM=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)-(*TM);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nTridiagonalMatrix - SymmetricBandMatrix" << endl;
  B=(*TM)-(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete TM; TM=0;
  }

  {
  cout << "\nSymmetricBandMatrix - SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *SP=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *S=(*A1)-(*SP);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSymmetricPositiveMatrix - SymmetricBandMatrix" << endl;
  S=(*SP)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SP; SP=0;
  }

  {
  cout << "\nSymmetricBandMatrix - SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *S=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *Sq=(*A1)-(*S);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix - SymmetricBandMatrix" << endl;
  Sq=(*S)-(*A2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete S; S=0;
  }

  {
  cout << "\nSymmetricBandMatrix - OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  SquareMatrix<F,Z> *S=(*A1)-(*Q);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nOrthogonalMatrix - SymmetricBandMatrix" << endl;
  S=(*Q)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete Q; Q=0;
  }

  {
  cout << "\nSymmetricBandMatrix - UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTriangularMatrix - SymmetricBandMatrix" << endl;
  S=(*U)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix - UnitUpperTrapezoidalMatrix" << endl;
  UnitUpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTrapezoidalMatrix - SymmetricBandMatrix" << endl;
  S=(*U)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix - UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTriangularMatrix - SymmetricBandMatrix" << endl;
  S=(*U)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix - UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTrapezoidalMatrix - SymmetricBandMatrix" << endl;
  S=(*U)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix - UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTriangularMatrix - SymmetricBandMatrix" << endl;
  S=(*L)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix - UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTrapezoidalMatrix - SymmetricBandMatrix" << endl;
  S=(*L)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix - LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTriangularMatrix - SymmetricBandMatrix" << endl;
  S=(*L)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix - LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTrapezoidalMatrix - SymmetricBandMatrix" << endl;
  S=(*L)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix - SquareMatrix" << endl;
  SquareMatrix<F,Z> *SM=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*SM);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSquareMatrix - SymmetricBandMatrix" << endl;
  S=(*SM)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SM; SM=0;
  }

  {
  cout << "\nSymmetricBandMatrix - Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)-(*M);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nMatrix - SymmetricBandMatrix" << endl;
  S=(*M)-(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete M; M=0;
  }

  {
  cout << "\nA1*scalar :" << endl;
  BandMatrix<F,Z> *T=(*A1)*scalar;
  T->printOn(cout);
  delete T; T=0;

  cout << "\nA2/scalar.:" << endl;
  T=(*A2)/scalar;
  T->printOn(cout);
  delete T; T=0;
  }

  {
  cout << "\nA1 * A2 = " << endl;
  BandMatrix<F,Z> *T=(*A1)*(*A2);
  T->printOn(cout);
  delete T; T=0;

  cout << "\nA2 * A1 = " << endl;
  T=(*A2)*(*A1);
  T->printOn(cout);
  delete T; T=0;
  }

  {
  cout << "\nSymmetricBandMatrix * BandMatrix" << endl;
  BandMatrix<F,Z> *B=OPERATOR_NEW BandMatrix<F,Z>(n,2,1,scalar);
  BandMatrix<F,Z> *S=(*A1)*(*B);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nBandMatrix * SymmetricBandMatrix" << endl;
  S=(*B)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete B; B=0;
  }

  {
  cout << "\nSymmetricBandMatrix * UpperHessenbergMatrix" << endl;
  UpperHessenbergMatrix<F,Z> *H=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)*(*H);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperHessenbergMatrix * SymmetricBandMatrix" << endl;
  S=(*H)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete H; H=0;
  }

  {
  cout << "\nSymmetricBandMatrix * DiagonalMatrix" << endl;
  DiagonalMatrix<F,Z> *D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)*(*D);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nDiagonalMatrix * SymmetricBandMatrix" << endl;
  B=(*D)*(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete D; D=0;
  }

  {
  cout << "\nSymmetricBandMatrix * SymmetricPositiveTridiagonalMatrix"
       << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *SPT=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)*(*SPT);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricPositiveTridiagonalMatrix * SymmetricBandMatrix"
       << endl;
  B=(*SPT)*(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete SPT; SPT=0;
  }

  {
  cout << "\nSymmetricBandMatrix * SymmetricTridiagonalMatrix" << endl;
  SymmetricTridiagonalMatrix<F,Z> *ST=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)*(*ST);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricTridiagonalMatrix * SymmetricBandMatrix" << endl;
  B=(*ST)*(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete ST; ST=0;
  }

  {
  cout << "\nSymmetricBandMatrix * TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *TM=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)*(*TM);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nTridiagonalMatrix * SymmetricBandMatrix" << endl;
  B=(*TM)*(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete TM; TM=0;
  }

  {
  cout << "\nSymmetricBandMatrix * SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *SP=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)*(*SP);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSymmetricPositiveMatrix * SymmetricBandMatrix" << endl;
  S=(*SP)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SP; SP=0;
  }

  {
  cout << "\nSymmetricBandMatrix * SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *S=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *Sq=(*A1)*(*S);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix * SymmetricBandMatrix" << endl;
  Sq=(*S)*(*A2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete S; S=0;
  }

  {
  cout << "\nSymmetricBandMatrix * OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  Matrix<F,Z> *M=(*A1)*(*Q);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nOrthogonalMatrix * SymmetricBandMatrix" << endl;
  M=(*Q)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete Q; Q=0;
  }

  {
  cout << "\nSymmetricBandMatrix * UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  Matrix<F,Z> *M=(*A1)*(*U);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nUnitUpperTriangularMatrix * SymmetricBandMatrix" << endl;
  M=(*U)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix * UnitUpperTrapezoidalMatrix" << endl;
  UnitUpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n+1,scalar);
  Matrix<F,Z> *M=(*A1)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  cout << "\nUnitUpperTrapezoidalMatrix * SymmetricBandMatrix" << endl;
  U=OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n-1,n,scalar);
  M=(*U)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix * UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  Matrix<F,Z> *M=(*A1)*(*U);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nUpperTriangularMatrix * SymmetricBandMatrix" << endl;
  M=(*U)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix * UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n+1,scalar);
  Matrix<F,Z> *M=(*A1)*(*U);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;

  cout << "\nUpperTrapezoidalMatrix * SymmetricBandMatrix" << endl;
  U=OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n-1,n,scalar);
  M=(*U)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricBandMatrix * UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  Matrix<F,Z> *M=(*A1)*(*L);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nUnitLowerTriangularMatrix * SymmetricBandMatrix" << endl;
  M=(*L)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix * UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n-1,scalar);
  Matrix<F,Z> *M=(*A1)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  cout << "\nUnitLowerTrapezoidalMatrix * SymmetricBandMatrix" << endl;
  L=OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n+1,n,scalar);
  M=(*L)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix * LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  Matrix<F,Z> *M=(*A1)*(*L);
  M->printOn(cout);
  delete M; M=0;

  cout << "\nLowerTriangularMatrix * SymmetricBandMatrix" << endl;
  M=(*L)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix * LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n-1,scalar);
  Matrix<F,Z> *M=(*A1)*(*L);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;

  cout << "\nLowerTrapezoidalMatrix * SymmetricBandMatrix" << endl;
  L=OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n+1,n,scalar);
  M=(*L)*(*A2);
  M->printOn(cout);
  delete M; M=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricBandMatrix * SquareMatrix" << endl;
  SquareMatrix<F,Z> *SM=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)*(*SM);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSquareMatrix * SymmetricBandMatrix" << endl;
  S=(*SM)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SM; SM=0;
  }

  {
  cout << "\nSymmetricBandMatrix * Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  Matrix<F,Z> *S=(*A1)*(*M);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nMatrix * SymmetricBandMatrix" << endl;
  S=(*M)*(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete M; M=0;
  }

  Vector<F,Z> *v=OPERATOR_NEW Vector<F,Z>(n,scalar);
  cout << "\nSymmetricBandMatrix * Vector" << endl;
  Vector<F,Z> *w=(*A1)*(*v);
  w->printOn(cout);
  delete w; w=0;

  w=(*A2)*(*v);
  w->printOn(cout);
  delete w; w=0;
  delete v; v=0;
  delete A1; A1=0;
  delete A2; A2=0;

  SymmetricBandMatrix<F,Z> *AA=
    OPERATOR_NEW SymmetricBandMatrix<F,Z>(3,1,1);
  (*AA)(0,0)=static_cast<F>(2.)*scalar;
  (*AA)(1,0)=static_cast<F>(-1.)*scalar;
    (*AA)(1,1)=static_cast<F>(2.)*scalar;
  (*AA)(2,1)=static_cast<F>(-1.)*scalar;
    (*AA)(2,2)=static_cast<F>(2.)*scalar;
//cout << "\nequilibrate" << endl;
//AA->printOn(cout);
//Vector<F,F> *rr=OPERATOR_NEW Vector<F,F>(3);
//F rowcnd,colcnd;
//Vector<F,F> *c=OPERATOR_NEW Vector<F,F>(3);
//F anorm=AA->equilibrate(*rr,*c,rowcnd,colcnd);
//cout << "rr = " << endl;
//rr->printOn(cout);
//cout << "c = " << endl;
//c->printOn(cout);
//cout << "rowcnd,colcnd = " << rowcnd << " " << colcnd << endl;
//delete rr; rr=0;
//delete c; c=0;

  cout << "\nnormFrobenius = " << AA->normFrobenius() << endl;
  cout << "normInfinity = " << AA->normInfinity() << endl;
  cout << "normMaxEntry = " << AA->normMaxEntry() << endl;
  cout << "normOne = " << AA->normOne() << endl;
//cout << "reciprocalConditionNumber('I') = "
//     << AA->reciprocalConditionNumber('I') << endl;
//cout << "reciprocalConditionNumber('O') = "
//     << AA->reciprocalConditionNumber('O') << endl;

//cout << "\ntranspose" << endl;
//T=AA->transpose();
//T->printOn(cout);
//delete T; T=0;

/*
  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(3);
  (*b)[0]=0.;
  (*b)[1]=6.;
  (*b)[2]=0.;
  cout << "\nsquare system of linear equations" << endl;
  cout << "A:" << endl;
  AA->printOn(cout);
  cout << "b:" << endl;
  b->printOn(cout);
  cout << "solve('N'):" << endl;
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(3);
  AA->solve(*b,*x);
  cout << "x = " << endl;
  x->printOn(cout);
  Vector<F,Z> *r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  AA->gbmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'N');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;

  cout << "\nsolve('T')" << endl;
  AA->solve(*b,*x,'T');
  x->printOn(cout);
  r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  AA->gbmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'T');
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "\nsolve('C')" << endl;
    AA->solve(*b,*x,'C');
    x->printOn(cout);
    r=OPERATOR_NEW Vector<F,Z>(b->size());
    r->copy(*b);
    AA->gbmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r,'C');
    cout << "r = " << endl;
    r->printOn(cout);
    delete r; r=0;
  }
  delete x; x=0;
  delete b; b=0;

  Matrix<F,Z> *B=OPERATOR_NEW Matrix<F,Z>(3,1);
  (*B)(0,0)=0.;
  (*B)(1,0)=6.;
  (*B)(2,0)=0.;
  cout << "\nsolve('L','N'):" << endl;
  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(3,1);
  AA->solve(*B,*X);
  X->printOn(cout);
  Matrix<F,Z> *R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*B);
  AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L','N');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  cout << "\nsolve('L','T'):" << endl;
  AA->solve(*B,*X,'L','T');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*B);
  AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L','T');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "\nsolve('L','C'):" << endl;
    AA->solve(*B,*X,'L','C');
    X->printOn(cout);
    R=OPERATOR_NEW Matrix<F,Z>(3,1);
    R->copy(*B);
    AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L','C');
    cout << "R = " << endl;
    R->printOn(cout);
    delete R; R=0;
  }
  delete X; X=0;
  delete B; B=0;

  B=OPERATOR_NEW Matrix<F,Z>(1,3);
  (*B)(0,0)=0.;
  (*B)(0,1)=6.;
  (*B)(0,2)=0.;
  cout << "\nsolve('R','N'):" << endl;
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  AA->solve(*B,*X,'R','N');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*B);
  AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R','N');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  cout << "\nsolve('R','T'):" << endl;
  AA->solve(*B,*X,'R','T');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*B);
  AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R','T');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;

  if (sizeof(F)!=sizeof(Z)) {
    cout << "\nsolve('R','C'):" << endl;
    AA->solve(*B,*X,'R','C');
    X->printOn(cout);
    R=OPERATOR_NEW Matrix<F,Z>(1,3);
    R->copy(*B);
    AA->gbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R','C');
    cout << "R = " << endl;
    R->printOn(cout);
    delete R; R=0;
  }
  delete X; X=0;
  delete B; B=0;
*/
  delete AA; AA=0;

  SymmetricBandMatrix<F,Z> *TTT=
    OPERATOR_NEW SymmetricBandMatrix<F,Z>(4,2);
  (*TTT)(0,0)=2.*fscalar;
  (*TTT)(1,0)=static_cast<F>(-1.)*scalar;
  (*TTT)(2,0)=static_cast<F>(-2.)*scalar;
    (*TTT)(1,1)=3.*fscalar;
    (*TTT)(2,1)=-2.*fscalar;
    (*TTT)(3,1)=-3.*fscalar;
      (*TTT)(2,2)=4.*fscalar;
      (*TTT)(3,2)=-3.*fscalar;
        (*TTT)(3,3)=5.*fscalar;
  OrthogonalMatrix<F,Z> *Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(4,4);
  Vector<F,F> *lambda=TTT->eigenvalues(Q);
  cout << "\nafter eigenvalues" << endl;
  lambda->printOn(cout);
  Q->printOn(cout);
  Matrix<F,Z> *RR=(*TTT)*(*Q);
  for (int j=0;j<4;j++) {
    for (int i=0;i<4;i++) (*RR)(i,j)-=(*Q)(i,j)*(*lambda)[j];
  }
  cout << "\nTTT * Q - Q * lambda" << endl;
  RR->printOn(cout);
  RR->gemm(Vector<F,Z>::one_,*Q,*Q,Vector<F,Z>::zero_,'C','N');
  cout << "\nQ^H * Q" << endl;
  RR->printOn(cout);
  delete RR; RR=0;
  delete Q; Q=0;
  delete lambda; lambda=0;
  delete TTT; TTT=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
SymmetricPositiveBandMatrix<F,Z>::printOn(ostream& s) const {
  s << "SymmetricPositiveBandMatrix(" << this->size(0) << " x "
    << this->size(1) << ")\n" ;
//s << "subDiags = " << nsub << endl;
  for (int i=0; i< this->size(0); i++) {
    for (int j=0; j< this->size(1); j++) s << (*this)(i,j) << "  ";
    s << "\n";
  }
}

template<typename F,typename Z> void testSymmetricPositiveBandMatrix(
F fscalar,Z scalar) {
  SymmetricPositiveBandMatrix<F,Z> *A1=
    OPERATOR_NEW SymmetricPositiveBandMatrix<F,Z>;
  cout << "\nafter SymmetricPositiveBandMatrix()" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);

  int n=4,nsub=2;
  cout << "\nn, nsub = " << n << " " << nsub << endl;
  SymmetricPositiveBandMatrix<F,Z> *A2=
    OPERATOR_NEW SymmetricPositiveBandMatrix<F,Z>(n,nsub);;
  cout << "\nafter SymmetricPositiveBandMatrix(n,nsub)" << endl;
  cout << "A2:" << endl;
  A2->printOn(cout);

  SymmetricPositiveBandMatrix<F,Z> *A3=
    OPERATOR_NEW SymmetricPositiveBandMatrix<F,Z>(n,nsub,scalar);
  cout << "\nafter SymmetricPositiveBandMatrix(n,nsub,scalar)" << endl;
  cout << "A3:" << endl;
  A3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=A3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }
  delete A3; A3=0;

  *A2=scalar;
  cout << "\nafter A2=scalar" << endl;
  cout << "A2:" << endl;
  A2->printOn(cout);

  A1->resize(n,1);
  cout << "\nafter A1->resize(n,nsub)" << endl;
  A1->printOn(cout);
  delete A1; A1=0;

  A1=OPERATOR_NEW SymmetricPositiveBandMatrix<F,Z>();
  A1->resize(*A2);
  cout << "\nafter A1->resize(A2)" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);

  A1->copy(*A2);
  cout << "\nafter A1->copy(A2)" << endl;
  cout << "A1:" << endl;
  A1->printOn(cout);
  delete A1; A1=0;

  cout << "\nA1+A2:" << endl;
  A1=OPERATOR_NEW SymmetricPositiveBandMatrix<F,Z>(n,1,scalar);
  SymmetricPositiveBandMatrix<F,Z> *T=(*A1)+(*A2);
  T->printOn(cout);
  delete T; T=0;

  {
  cout << "\nSymmetricPositiveBandMatrix + SymmetricBandMatrix" << endl;
  SymmetricBandMatrix<F,Z> *B=
    OPERATOR_NEW SymmetricBandMatrix<F,Z>(n,2,scalar);
  SymmetricBandMatrix<F,Z> *S=(*A1)+(*B);
  S->printOn(cout);
  delete S; S=0;
  delete B; B=0;

  cout << "\nSymmetricBandMatrix + SymmetricPositiveBandMatrix" << endl;
  B=OPERATOR_NEW SymmetricBandMatrix<F,Z>(n,1,scalar);
  S=(*B)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete B; B=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + BandMatrix" << endl;
  BandMatrix<F,Z> *B=OPERATOR_NEW BandMatrix<F,Z>(n,1,2,scalar);
  BandMatrix<F,Z> *S=(*A1)+(*B);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nBandMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*B)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete B; B=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + UpperHessenbergMatrix" << endl;
  UpperHessenbergMatrix<F,Z> *H=
    OPERATOR_NEW UpperHessenbergMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*H);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperHessenbergMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*H)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete H; H=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + DiagonalMatrix" << endl;
  DiagonalMatrix<F,Z> *D=OPERATOR_NEW DiagonalMatrix<F,Z>(n,scalar);
  SymmetricBandMatrix<F,Z> *B=(*A1)+(*D);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nDiagonalMatrix + SymmetricPositiveBandMatrix" << endl;
  B=(*D)+(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete D; D=0;
  }

  {
  cout
    << "\nSymmetricPositiveBandMatrix + SymmetricPositiveTridiagonalMatrix"
    << endl;
  SymmetricPositiveTridiagonalMatrix<F,Z> *SPT=
    OPERATOR_NEW SymmetricPositiveTridiagonalMatrix<F,Z>(n,scalar);
  SymmetricPositiveBandMatrix<F,Z> *B=(*A1)+(*SPT);
  B->printOn(cout);
  delete B; B=0;

  cout
    << "\nSymmetricPositiveTridiagonalMatrix + SymmetricPositiveBandMatrix"
    << endl;
  B=(*SPT)+(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete SPT; SPT=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + SymmetricTridiagonalMatrix"
       << endl;
  SymmetricTridiagonalMatrix<F,Z> *ST=
    OPERATOR_NEW SymmetricTridiagonalMatrix<F,Z>(n,scalar);
  SymmetricBandMatrix<F,Z> *B=(*A1)+(*ST);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nSymmetricTridiagonalMatrix + SymmetricPositiveBandMatrix"
       << endl;
  B=(*ST)+(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete ST; ST=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + TridiagonalMatrix" << endl;
  TridiagonalMatrix<F,Z> *TM=OPERATOR_NEW TridiagonalMatrix<F,Z>(n,scalar);
  BandMatrix<F,Z> *B=(*A1)+(*TM);
  B->printOn(cout);
  delete B; B=0;

  cout << "\nTridiagonalMatrix + SymmetricPositiveBandMatrix" << endl;
  B=(*TM)+(*A2);
  B->printOn(cout);
  delete B; B=0;
  delete TM; TM=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + SymmetricPositiveMatrix" << endl;
  SymmetricPositiveMatrix<F,Z> *SP=
    OPERATOR_NEW SymmetricPositiveMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *S=(*A1)+(*SP);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSymmetricPositiveMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*SP)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SP; SP=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + SymmetricMatrix" << endl;
  SymmetricMatrix<F,Z> *S=OPERATOR_NEW SymmetricMatrix<F,Z>(n,scalar);
  SymmetricMatrix<F,Z> *Sq=(*A1)+(*S);
  Sq->printOn(cout);
  delete Sq; Sq=0;

  cout << "\nSymmetricMatrix + SymmetricPositiveBandMatrix" << endl;
  Sq=(*S)+(*A2);
  Sq->printOn(cout);
  delete Sq; Sq=0;
  delete S; S=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + OrthogonalMatrix" << endl;
  OrthogonalMatrix<F,Z> *Q=OPERATOR_NEW OrthogonalMatrix<F,Z>(n,n);
  SquareMatrix<F,Z> *S=(*A1)+(*Q);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nOrthogonalMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*Q)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete Q; Q=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + UnitUpperTriangularMatrix" << endl;
  UnitUpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTriangularMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*U)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + UnitUpperTrapezoidalMatrix" << endl;
  UnitUpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UnitUpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitUpperTrapezoidalMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*U)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + UpperTriangularMatrix" << endl;
  UpperTriangularMatrix<F,Z> *U=
    OPERATOR_NEW UpperTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTriangularMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*U)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + UpperTrapezoidalMatrix" << endl;
  UpperTrapezoidalMatrix<F,Z> *U=
    OPERATOR_NEW UpperTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*U);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUpperTrapezoidalMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*U)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete U; U=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + UnitLowerTriangularMatrix" << endl;
  UnitLowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTriangularMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*L)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + UnitLowerTrapezoidalMatrix" << endl;
  UnitLowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nUnitLowerTrapezoidalMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*L)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + LowerTriangularMatrix" << endl;
  LowerTriangularMatrix<F,Z> *L=
    OPERATOR_NEW LowerTriangularMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTriangularMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*L)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + LowerTrapezoidalMatrix" << endl;
  LowerTrapezoidalMatrix<F,Z> *L=
    OPERATOR_NEW LowerTrapezoidalMatrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*L);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nLowerTrapezoidalMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*L)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete L; L=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + SquareMatrix" << endl;
  SquareMatrix<F,Z> *SM=OPERATOR_NEW SquareMatrix<F,Z>(n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*SM);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nSquareMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*SM)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete SM; SM=0;
  }

  {
  cout << "\nSymmetricPositiveBandMatrix + Matrix" << endl;
  Matrix<F,Z> *M=OPERATOR_NEW Matrix<F,Z>(n,n,scalar);
  SquareMatrix<F,Z> *S=(*A1)+(*M);
  S->printOn(cout);
  delete S; S=0;

  cout << "\nMatrix + SymmetricPositiveBandMatrix" << endl;
  S=(*M)+(*A2);
  S->printOn(cout);
  delete S; S=0;
  delete M; M=0;
  }
  delete A2; A2=0;
  delete A1; A1=0;

  SymmetricPositiveBandMatrix<F,Z> *AA=
    OPERATOR_NEW SymmetricPositiveBandMatrix<F,Z>(3,1);
  (*AA)(0,0)=2.;
  (*AA)(1,0)=-1.; (*AA)(1,1)=2.;
  (*AA)(2,1)=-1.; (*AA)(2,2)=2.;
  cout << "\nAA = " << endl;
  AA->printOn(cout);

  cout << "\nequilibrate" << endl;
  Vector<F,F> *s=OPERATOR_NEW Vector<F,F>(3);
  F scond;
  F amax=AA->equilibrate(*s,scond);
  cout << "amax,scond = " << amax << " " << scond << endl;
  s->printOn(cout);
  delete s; s=0;

  cout << "\nreciprocalConditionNumber = "
       << AA->reciprocalConditionNumber() << endl;

//OrthogonalMatrix<F,Z> *O=OPERATOR_NEW OrthogonalMatrix<F,Z>(3,3);
//cout << "\neigenvalues" << endl;
//Vector<F,F> *lambda=AA->eigenvalues(O);
//lambda->printOn(cout);
//O->printOn(cout);
//delete O; O=0;
//delete lambda; lambda=0;

  Vector<F,Z> *b=OPERATOR_NEW Vector<F,Z>(3);
  (*b)[0]=0.;
  (*b)[1]=6.;
  (*b)[2]=0.;
  cout << "\nsystem of linear equations" << endl;
  cout << "A:" << endl;
  AA->printOn(cout);
  cout << "b:" << endl;
  b->printOn(cout);
  Vector<F,Z> *x=OPERATOR_NEW Vector<F,Z>(3);
  cout << "AA->solve(b,x):" << endl;
  AA->solve(*b,*x);
  x->printOn(cout);
  Vector<F,Z> *r=OPERATOR_NEW Vector<F,Z>(b->size());
  r->copy(*b);
  AA->sbmv(Vector<F,Z>::one_,*x,Vector<F,Z>::mone_,*r);
  cout << "r = " << endl;
  r->printOn(cout);
  delete r; r=0;
  delete x; x=0;
  delete b; b=0;

  Matrix<F,Z> *B=OPERATOR_NEW Matrix<F,Z>(3,1);
  (*B)(0,0)=0.;
  (*B)(1,0)=6.;
  (*B)(2,0)=0.;
  cout << "\nsystem of linear equations" << endl;
  cout << "B:" << endl;
  B->printOn(cout);
  Matrix<F,Z> *X=OPERATOR_NEW Matrix<F,Z>(3,1);
  cout << "AA->solve(B,X,'L'):" << endl;
  AA->solve(*B,*X);
  X->printOn(cout);
  Matrix<F,Z> *R=OPERATOR_NEW Matrix<F,Z>(3,1);
  R->copy(*B);
  AA->sbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'L');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete B; B=0;

  B=OPERATOR_NEW Matrix<F,Z>(1,3);
  (*B)(0,0)=0.;
  (*B)(0,1)=6.;
  (*B)(0,2)=0.;
  cout << "\nsystem of linear equations" << endl;
  cout << "B:" << endl;
  B->printOn(cout);
  X=OPERATOR_NEW Matrix<F,Z>(1,3);
  cout << "AA->solve(B,X,'R'):" << endl;
  AA->solve(*B,*X,'R');
  X->printOn(cout);
  R=OPERATOR_NEW Matrix<F,Z>(1,3);
  R->copy(*B);
  AA->sbmm(Vector<F,Z>::one_,*X,Vector<F,Z>::mone_,*R,'R');
  cout << "R = " << endl;
  R->printOn(cout);
  delete R; R=0;
  delete X; X=0;
  delete B; B=0;
  delete AA; AA=0;
}

// Modified from LaPack++ gmc.C by John Trangenstein, 11/7/96

//
//              LAPACK++ 1.1 Linear Algebra Package 1.1
//               University of Tennessee, Knoxvilee, TN.
//            Oak Ridge National Laboratory, Oak Ridge, TN.
//        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
//                 (C) 1992-1996 All Rights Reserved
//
//                             NOTICE
//
// Permission to use, copy, modify, and distribute this software and
// its documentation for any purpose and without fee is hereby granted
// provided that the above copyright notice appear in all copies and
// that both the copyright notice and this permission notice appear in
// supporting documentation.
//
// Neither the Institutions (University of Tennessee, and Oak Ridge 
// National Laboratory) nor the Authors make any representations about 
// the suitability of this software for any purpose.  This software is 
// provided ``as is'' without express or implied warranty.
//
// LAPACK++ was funded in part by the U.S. Department of Energy, the
// National Science Foundation and the State of Tennessee.
