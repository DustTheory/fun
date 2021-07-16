#include "OrthogonalMatrix.H"

template<typename F,typename Z>
OrthogonalMatrix<F,Z>::OrthogonalMatrix(int m,int n) : 
Matrix<F,Z>(m,n,OrthogonalMatrix<double,double>::zero_) {
  for (int i=0;i<min(m,n);i++) {
    this->operator()(i,i)=OrthogonalMatrix<double,double>::one_;
  }
}

template<typename F,typename Z> void
OrthogonalMatrix<F,Z>::resize(int m,int n) {
  Matrix<F,Z>::resize(m,n);
  fillWith(OrthogonalMatrix<F,Z>::zero_);
  for (int i=0;i<min(m,n);i++) { 
    this->operator()(i,i)=OrthogonalMatrix<F,Z>::one_;
  }
}

template<typename F,typename Z> void OrthogonalMatrix<F,Z>::printOn(ostream &s) const {
  s << "OrthogonalMatrix(" << this->size(0) << " x " << this->size(1) 
    << ")\n" ;
  int m=this->size(0), g=this->size(1);
  for (int i=0; i< m; i++) {
    for (int j=0; j< g; j++) s << this->operator()(i,j) << "  ";
    s << "\n";
  }
}

template<typename F,typename Z> void testOrthogonalMatrix(F,Z scalar) {
  {
    int m=3,n=2;
    cout << "\nm,n = " << m << " " << n << endl;
    OrthogonalMatrix<F,Z> Q2(m,n);;
    cout << "\nafter OrthogonalMatrix(m,n)" << endl;
    cout << "Q2:" << endl;
    Q2.printOn(cout);

    {
      OrthogonalMatrix<F,Z> Q1;
      cout << "\nafter OrthogonalMatrix()" << endl;
      cout << "Q1:" << endl;
      Q1.printOn(cout);


      cout << "\nafter resize(m,n)" << endl;
      Q1.resize(m,n);
      cout << "Q1:" << endl;
      Q1.printOn(cout);
    }

    {
      OrthogonalMatrix<F,Z> Q1;
      Q1.resize(Q2);
      cout << "\nafter Q1.resize(Q2)" << endl;
      cout << "Q1:" << endl;
      Q1.printOn(cout);
    }

    {
      OrthogonalMatrix<F,Z> Q1(m,n);
      Q1.copy(Q2);
      cout << "\nafter Q1.copy(Q2)" << endl;
      cout << "Q1:" << endl;
      Q1.printOn(cout);
    }

    cout << "\nafter Q2 *= -1" << endl;
    Q2*=OrthogonalMatrix<F,Z>::mone_;
    cout << "Q2:" << endl;
    Q2.printOn(cout);

    cout << "\nafter Q2 /= -1" << endl;
    Q2/=OrthogonalMatrix<F,Z>::mone_;
    cout << "Q2:" << endl;
    Q2.printOn(cout);

    cout << "\noperator*(scalar)" << endl;
    Matrix<F,Z> *M=Q2*scalar;
    M->printOn(cout);
    delete M; M=0;

    cout << "\noperator/(scalar)" << endl;
    M=Q2/scalar;
    M->printOn(cout);
    delete M; M=0;

//  cout << "\nafter transpose" << endl;
//  Q1=Q2.transpose();
//  Q1->printOn(cout);
//  delete Q1; Q1=0;

//  if (sizeof(F)!=sizeof(Z)) {
//    cout << "\nafter conjugateTranspose" << endl;
//    Q1=Q2.conjugateTranspose();
//    Q1->printOn(cout);
//    delete Q1; Q1=0;
//  }
  }

  OrthogonalMatrix<F,Z> Q(4,3);
  Q(0,0)=Q(1,0)=Q(2,0)=Q(3,0)=Q(1,1)=Q(3,1)=Q(2,2)=Q(3,2)=0.5;
  Q(0,1)=Q(2,1)=Q(0,2)=Q(1,2)=-0.5;
  cout << "\nQ = " << endl;
  Q.printOn(cout);

  {
    Vector<F,Z> b(Q.size(0),scalar);
    cout << "\nsolve(Vector,'N')" << endl;
    Vector<F,Z> x(Q.size(1));
    Q.solve(b,x);
    x.printOn(cout);
    Vector<F,Z> r(b.size());
    r.copy(b);
    Q.gemv(Vector<F,Z>::one_,x,Vector<F,Z>::mone_,r,'N');
    cout << "residual = " << endl;
    r.printOn(cout);
    Vector<F,Z> qr(Q.size(1),Vector<F,Z>::zero_);
    Q.gemv(Vector<F,Z>::one_,r,Vector<F,Z>::zero_,qr,'C');
    cout << "Q^H * residual = " << endl;
    qr.printOn(cout);
  }

  if (sizeof(F)!=sizeof(Z)) {
    Vector<F,Z> b(Q.size(1),scalar);
    cout << "\nsolve(Vector,'T')" << endl;
    Vector<F,Z> x(Q.size(0));
    Q.solve(b,x,'T');
    x.printOn(cout);
    Vector<F,Z> r(b.size());
    r.copy(b);
    Q.gemv(Vector<F,Z>::one_,x,Vector<F,Z>::mone_,r,'T');
    cout << "residual = " << endl;
    r.printOn(cout);
  }

  {
    Vector<F,Z> b(Q.size(1),scalar);
    cout << "\nsolve(Vector,'C')" << endl;
    Vector<F,Z> x(Q.size(0));
    Q.solve(b,x,'C');
    x.printOn(cout);
    {
      Vector<F,Z> r(b.size());
      r.copy(b);
      Q.gemv(Vector<F,Z>::one_,x,Vector<F,Z>::mone_,r,'C');
      cout << "residual = " << endl;
      r.printOn(cout);
    }
    {
      Vector<F,Z> s(Q.size(1),Vector<F,Z>::zero_);
      Q.solve(x,s,'N');
      Vector<F,Z> r(x.size());
      r.copy(x);
      Q.gemv(Vector<F,Z>::one_,s,Vector<F,Z>::mone_,r,'N');
      cout << "Q*s-x = " << endl;
      r.printOn(cout);
    }
  }

/*
  {
    UpperTrapezoidalMatrix<F,Z> U(Q.size(0),Q.size(0)+1,scalar);
    cout << "\nsolve(UnitUpperTrapezoidalMatrix,'L','N')" << endl;
    Matrix<F,Z> X(Q.size(1),U.size(1));
    Q.solve(U,X);
    X.printOn(cout);
    Matrix<F,Z> R(U.size(0),U.size(1),Matrix<F,Z>::zero_);
    for (int i=0;i<U.size(0);i++) {
      R(i,i)=Matrix<F,Z>::one_;
      for (int j=i+1;j<U.size(1);j++) R(i,j)=U(i,j);
    }
    R.gemm('N','N',Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_);
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(Q.size(1),R.size(1),Matrix<F,Z>::zero_);
    AR.gemm('C','N',Matrix<F,Z>::one_,Q,R,Matrix<F,Z>::zero_);
    cout << "Q^H * residual = " << endl;
    AR.printOn(cout);
  }

  {
    UnitUpperTrapezoidalMatrix<F,Z> U(Q.size(1),Q.size(1)+1,scalar);
    cout << "\nsolve(UnitUpperTrapezoidalMatrix,'L','T')" << endl;
    Matrix<F,Z> X(Q.size(0),U.size(1));
    Q.solve(U,X,'L','T');
    X.printOn(cout);
    {
      Matrix<F,Z> R(U.size(0),U.size(1),Matrix<F,Z>::zero_);
      for (int i=0;i<U.size(0);i++) {
        R(i,i)=Matrix<F,Z>::one_;
        for (int j=i+1;j<U.size(1);j++) R(i,j)=U(i,j);
      }
      R.gemm('C','N',Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_);
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(Q.size(1),U.size(1),Matrix<F,Z>::zero_);
      Q.solve(X,S,'L','N');
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm('N','N',Matrix<F,Z>::one_,Q,S,Matrix<F,Z>::mone_);
      cout << "Q*s-x = " << endl;
      R.printOn(cout);
    }
  }

  {
    UnitUpperTrapezoidalMatrix<F,Z> U(Q.size(1)-1,Q.size(1),scalar);
    cout << "\nsolve(UnitUpperTrapezoidalMatrix,'R','N')" << endl;
    Matrix<F,Z> X(U.size(0),Q.size(0));
    Q.solve(U,X,'R','N');
    X.printOn(cout);
    {
      Matrix<F,Z> R(U.size(0),U.size(1),Matrix<F,Z>::zero_);
      for (int i=0;i<U.size(0);i++) {
        R(i,i)=Matrix<F,Z>::one_;
        for (int j=i+1;j<U.size(1);j++) R(i,j)=U(i,j);
      }
      R.gemm('N','N',Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_);
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(U.size(0),Q.size(1),Matrix<F,Z>::zero_);
      if (sizeof(F)==sizeof(Z)) {
        Q.solve(X,S,'R','T');
      } else {
        Q.solve(X,S,'R','C');
      }
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm('N','C',Matrix<F,Z>::one_,S,Q,Matrix<F,Z>::mone_);
      cout << "S*Q^T - X = " << endl;
      R.printOn(cout);
    }
  }

  {
    UnitUpperTrapezoidalMatrix<F,Z> U(Q.size(0)-1,Q.size(0),scalar);
    cout << "\nsolve(UnitUpperTrapezoidalMatrix,'R','T')" << endl;
    Matrix<F,Z> X(U.size(0),Q.size(1),scalar);
    Q.solve(U,X,'R','T');
    X.printOn(cout);
    Matrix<F,Z> R(U.size(0),U.size(1),Matrix<F,Z>::zero_);
    for (int i=0;i<U.size(0);i++) {
      R(i,i)=Matrix<F,Z>::one_;
      for (int j=i+1;j<U.size(1);j++) R(i,j)=U(i,j);
    }
    R.gemm('N','C',Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_);
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(R.size(0),Q.size(1),Matrix<F,Z>::zero_);
    AR.gemm('N','N',Matrix<F,Z>::one_,R,Q,Matrix<F,Z>::zero_);
    cout << "R * Q = " << endl;
    AR.printOn(cout);
  }
*/

/*
  {
    UpperTrapezoidalMatrix<F,Z> U(Q.size(0),Q.size(0)+1,scalar);
    cout << "\nsolve(UpperTrapezoidalMatrix,'L','N')" << endl;
    Matrix<F,Z> X(Q.size(1),U.size(1));
    Q.solve(U,X);
    X.printOn(cout);
    Matrix<F,Z> R(U.size(0),U.size(1),Matrix<F,Z>::zero_);
    for (int i=0;i<U.size(0);i++) {
      for (int j=i;j<U.size(1);j++) R(i,j)=U(i,j);
    }
    R.gemm('N','N',Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_);
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(Q.size(1),R.size(1),Matrix<F,Z>::zero_);
    AR.gemm('C','N',Matrix<F,Z>::one_,Q,R,Matrix<F,Z>::zero_);
    cout << "Q^H * residual = " << endl;
    AR.printOn(cout);
  }

  {
    UpperTrapezoidalMatrix<F,Z> U(Q.size(1),Q.size(1)+1,scalar);
    cout << "\nsolve(UpperTrapezoidalMatrix,'L','T')" << endl;
    Matrix<F,Z> X(Q.size(0),U.size(1));
    Q.solve(U,X,'L','T');
    X.printOn(cout);
    {
      Matrix<F,Z> R(U.size(0),U.size(1),Matrix<F,Z>::zero_);
      for (int i=0;i<U.size(0);i++) {
        for (int j=i;j<U.size(1);j++) R(i,j)=U(i,j);
      }
      R.gemm('C','N',Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_);
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(Q.size(1),U.size(1),Matrix<F,Z>::zero_);
      Q.solve(X,S,'L','N');
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm('N','N',Matrix<F,Z>::one_,Q,S,Matrix<F,Z>::mone_);
      cout << "Q*s-x = " << endl;
      R.printOn(cout);
    }
  }

  {
    UpperTrapezoidalMatrix<F,Z> U(Q.size(1)-1,Q.size(1),scalar);
    cout << "\nsolve(UpperTrapezoidalMatrix,'R','N')" << endl;
    Matrix<F,Z> X(U.size(0),Q.size(0));
    Q.solve(U,X,'R','N');
    X.printOn(cout);
    {
      Matrix<F,Z> R(U.size(0),U.size(1),Matrix<F,Z>::zero_);
      for (int i=0;i<U.size(0);i++) {
        for (int j=i;j<U.size(1);j++) R(i,j)=U(i,j);
      }
      R.gemm('N','N',Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_);
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(U.size(0),Q.size(1),Matrix<F,Z>::zero_);
      if (sizeof(F)==sizeof(Z)) {
      Q.solve(X,S,'R','T');
      } else {
      Q.solve(X,S,'R','C');
      }
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm('N','C',Matrix<F,Z>::one_,S,Q,Matrix<F,Z>::mone_);
      cout << "S*Q^T - X = " << endl;
      R.printOn(cout);
    }
  }

  {
    UpperTrapezoidalMatrix<F,Z> U(Q.size(0)-1,Q.size(0),scalar);
    cout << "\nsolve(UpperTrapezoidalMatrix,'R','T')" << endl;
    Matrix<F,Z> X(U.size(0),Q.size(1),scalar);
    Q.solve(U,X,'R','T');
    X.printOn(cout);
    Matrix<F,Z> R(U.size(0),U.size(1),Matrix<F,Z>::zero_);
    for (int i=0;i<U.size(0);i++) {
      for (int j=i;j<U.size(1);j++) R(i,j)=U(i,j);
    }
    R.gemm('N','C',Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_);
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(R.size(0),Q.size(1),Matrix<F,Z>::zero_);
    AR.gemm('N','N',Matrix<F,Z>::one_,R,Q,Matrix<F,Z>::zero_);
    cout << "R * Q = " << endl;
    AR.printOn(cout);
  }
*/

/*
  {
    LowerTrapezoidalMatrix<F,Z> L(Q.size(0),Q.size(0)-1,scalar);
    cout << "\nsolve(UnitLowerTrapezoidalMatrix,'L','N')" << endl;
    Matrix<F,Z> X(Q.size(1),L.size(1));
    Q.solve(L,X);
    X.printOn(cout);
    Matrix<F,Z> R(L.size(0),L.size(1),Matrix<F,Z>::zero_);
    for (int j=0;j<L.size(1);j++) {
      for (int i=j+1;i<L.size(0);i++) R(i,j)=L(i,j);
      R(j,j)=Matrix<F,Z>::one_;
    }
    R.gemm('N','N',Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_);
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(Q.size(1),R.size(1),Matrix<F,Z>::zero_);
    AR.gemm('C','N',Matrix<F,Z>::one_,Q,R,Matrix<F,Z>::zero_);
    cout << "Q^H * residual = " << endl;
    AR.printOn(cout);
  }

  {
    UnitLowerTrapezoidalMatrix<F,Z> L(Q.size(1),Q.size(1)-1,scalar);
    cout << "\nsolve(UnitLowerTrapezoidalMatrix,'L','T')" << endl;
    Matrix<F,Z> X(Q.size(0),L.size(1));
    Q.solve(L,X,'L','T');
    X.printOn(cout);
    {
      Matrix<F,Z> R(L.size(0),L.size(1),Matrix<F,Z>::zero_);
      for (int j=0;j<L.size(1);j++) {
        for (int i=j+1;i<L.size(0);i++) R(i,j)=L(i,j);
        R(j,j)=Matrix<F,Z>::one_;
      }
      R.gemm('C','N',Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_);
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(Q.size(1),L.size(1),Matrix<F,Z>::zero_);
      Q.solve(X,S,'L','N');
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm('N','N',Matrix<F,Z>::one_,Q,S,Matrix<F,Z>::mone_);
      cout << "Q*s-x = " << endl;
      R.printOn(cout);
    }
  }

  {
    UnitLowerTrapezoidalMatrix<F,Z> L(Q.size(1)+1,Q.size(1),scalar);
    cout << "\nsolve(UnitLowerTrapezoidalMatrix,'R','N')" << endl;
    Matrix<F,Z> X(L.size(0),Q.size(0));
    Q.solve(L,X,'R','N');
    X.printOn(cout);
    {
      Matrix<F,Z> R(L.size(0),L.size(1),Matrix<F,Z>::zero_);
      for (int j=0;j<L.size(1);j++) {
        for (int i=j+1;i<L.size(0);i++) R(i,j)=L(i,j);
        R(j,j)=Matrix<F,Z>::one_;
      }
      R.gemm('N','N',Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_);
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(L.size(0),Q.size(1),Matrix<F,Z>::zero_);
      if (sizeof(F)==sizeof(Z)) {
      Q.solve(X,S,'R','T');
      } else {
      Q.solve(X,S,'R','C');
      }
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm('N','C',Matrix<F,Z>::one_,S,Q,Matrix<F,Z>::mone_);
      cout << "S*Q^T - X = " << endl;
      R.printOn(cout);
    }
  }

  {
    UnitLowerTrapezoidalMatrix<F,Z> L(Q.size(0)+1,Q.size(0),scalar);
    cout << "\nsolve(UnitLowerTrapezoidalMatrix,'R','T')" << endl;
    Matrix<F,Z> X(L.size(0),Q.size(1),scalar);
    Q.solve(L,X,'R','T');
    X.printOn(cout);
    Matrix<F,Z> R(L.size(0),L.size(1),Matrix<F,Z>::zero_);
    for (int j=0;j<L.size(1);j++) {
      for (int i=j+1;i<L.size(0);i++) R(i,j)=L(i,j);
      R(j,j)=Matrix<F,Z>::one_;
    }
    R.gemm('N','C',Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_);
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(R.size(0),Q.size(1),Matrix<F,Z>::zero_);
    AR.gemm('N','N',Matrix<F,Z>::one_,R,Q,Matrix<F,Z>::zero_);
    cout << "R * Q = " << endl;
    AR.printOn(cout);
  }
*/

/*
  {
    LowerTrapezoidalMatrix<F,Z> L(Q.size(0),Q.size(0)-1,scalar);
    cout << "\nsolve(LowerTrapezoidalMatrix,'L','N')" << endl;
    Matrix<F,Z> X(Q.size(1),L.size(1));
    Q.solve(L,X);
    X.printOn(cout);
    Matrix<F,Z> R(L.size(0),L.size(1),Matrix<F,Z>::zero_);
    for (int j=0;j<L.size(1);j++) {
      for (int i=j;i<L.size(0);i++) R(i,j)=L(i,j);
    }
    R.gemm('N','N',Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_);
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(Q.size(1),R.size(1),Matrix<F,Z>::zero_);
    AR.gemm('C','N',Matrix<F,Z>::one_,Q,R,Matrix<F,Z>::zero_);
    cout << "Q^H * residual = " << endl;
    AR.printOn(cout);
  }

  {
    LowerTrapezoidalMatrix<F,Z> L(Q.size(1),Q.size(1)-1,scalar);
    cout << "\nsolve(LowerTrapezoidalMatrix,'L','T')" << endl;
    Matrix<F,Z> X(Q.size(0),L.size(1));
    Q.solve(L,X,'L','T');
    X.printOn(cout);
    {
      Matrix<F,Z> R(L.size(0),L.size(1),Matrix<F,Z>::zero_);
      for (int j=0;j<L.size(1);j++) {
        for (int i=j;i<L.size(0);i++) R(i,j)=L(i,j);
      }
      R.gemm('C','N',Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_);
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(Q.size(1),L.size(1),Matrix<F,Z>::zero_);
      Q.solve(X,S,'L','N');
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm('N','N',Matrix<F,Z>::one_,Q,S,Matrix<F,Z>::mone_);
      cout << "Q*s-x = " << endl;
      R.printOn(cout);
    }
  }

  {
    LowerTrapezoidalMatrix<F,Z> L(Q.size(1)+1,Q.size(1),scalar);
    cout << "\nsolve(LowerTrapezoidalMatrix,'R','N')" << endl;
    Matrix<F,Z> X(L.size(0),Q.size(0));
    Q.solve(L,X,'R','N');
    X.printOn(cout);
    {
      Matrix<F,Z> R(L.size(0),L.size(1),Matrix<F,Z>::zero_);
      for (int j=0;j<L.size(1);j++) {
        for (int i=j;i<L.size(0);i++) R(i,j)=L(i,j);
      }
      R.gemm('N','N',Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_);
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(L.size(0),Q.size(1),Matrix<F,Z>::zero_);
      if (sizeof(F)==sizeof(Z)) {
      Q.solve(X,S,'R','T');
      } else {
      Q.solve(X,S,'R','C');
      }
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm('N','C',Matrix<F,Z>::one_,S,Q,Matrix<F,Z>::mone_);
      cout << "S*Q^T - X = " << endl;
      R.printOn(cout);
    }
  }

  {
    LowerTrapezoidalMatrix<F,Z> L(Q.size(0)+1,Q.size(0),scalar);
    cout << "\nsolve(LowerTrapezoidalMatrix,'R','T')" << endl;
    Matrix<F,Z> X(L.size(0),Q.size(1),scalar);
    Q.solve(L,X,'R','T');
    X.printOn(cout);
    Matrix<F,Z> R(L.size(0),L.size(1),Matrix<F,Z>::zero_);
    for (int j=0;j<L.size(1);j++) {
      for (int i=j;i<L.size(0);i++) R(i,j)=L(i,j);
    }
    R.gemm('N','C',Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_);
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(R.size(0),Q.size(1),Matrix<F,Z>::zero_);
    AR.gemm('N','N',Matrix<F,Z>::one_,R,Q,Matrix<F,Z>::zero_);
    cout << "R * Q = " << endl;
    AR.printOn(cout);
  }
*/

  {
    Matrix<F,Z> B(Q.size(0),1,scalar);
    cout << "\nsolve(Matrix,'L','N')" << endl;
    Matrix<F,Z> X(Q.size(1),1);
    Q.solve(B,X);
    X.printOn(cout);
    Matrix<F,Z> R(B.size(0),B.size(1));
    R.copy(B);
    R.gemm(Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_,'N','N');
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(Q.size(1),R.size(1),Matrix<F,Z>::zero_);
    AR.gemm(Matrix<F,Z>::one_,Q,R,Matrix<F,Z>::zero_,'C','N');
    cout << "Q^H * residual = " << endl;
    AR.printOn(cout);
  }

  if (sizeof(F)!=sizeof(Z)) {
    Matrix<F,Z> B(Q.size(1),1,scalar);
    cout << "\nsolve(Matrix,'L','C')" << endl;
    Matrix<F,Z> X(Q.size(0),1);
    Q.solve(B,X,'L','C');
    X.printOn(cout);
    {
      Matrix<F,Z> R(B.size(0),B.size(1));
      R.copy(B);
      R.gemm(Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_,'C','N');
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(Q.size(1),B.size(1),Matrix<F,Z>::zero_);
      Q.solve(X,S,'L','N');
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm(Matrix<F,Z>::one_,Q,S,Matrix<F,Z>::mone_,'N','N');
      cout << "Q*s-x = " << endl;
      R.printOn(cout);
    }
  }

  {
    Matrix<F,Z> B(Q.size(1),1,scalar);
    cout << "\nsolve(Matrix,'L','T')" << endl;
    Matrix<F,Z> X(Q.size(0),1);
    Q.solve(B,X,'L','T');
    X.printOn(cout);
    {
      Matrix<F,Z> R(B.size(0),B.size(1));
      R.copy(B);
      R.gemm(Matrix<F,Z>::one_,Q,X,Matrix<F,Z>::mone_,'T','N');
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    if (sizeof(F)!=sizeof(Z)) {
      Matrix<F,Z> S(Q.size(1),B.size(1),Matrix<F,Z>::zero_);
      Q.solve(X,S,'L','N');
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm(Matrix<F,Z>::one_,Q,S,Matrix<F,Z>::mone_,'N','N');
      cout << "Q*s-x = " << endl;
      R.printOn(cout);
    }
  }

  {
    Matrix<F,Z> B(1,Q.size(1),scalar);
    cout << "\nsolve(Matrix,'R','N')" << endl;
    Matrix<F,Z> X(1,Q.size(0));
    Q.solve(B,X,'R','N');
    X.printOn(cout);
    {
      Matrix<F,Z> R(B.size(0),B.size(1));
      R.copy(B);
      R.gemm(Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_,'N','N');
      cout << "residual = " << endl;
      R.printOn(cout);
    }
    {
      Matrix<F,Z> S(B.size(0),Q.size(1),Matrix<F,Z>::zero_);
      if (sizeof(F)==sizeof(Z)) {
      Q.solve(X,S,'R','T');
      } else {
      Q.solve(X,S,'R','C');
      }
      Matrix<F,Z> R(X.size(0),X.size(1));
      R.copy(X);
      R.gemm(Matrix<F,Z>::one_,S,Q,Matrix<F,Z>::mone_,'N','C');
      cout << "S*Q^T - X = " << endl;
      R.printOn(cout);
    }
  }

  if (sizeof(F)!=sizeof(Z)) {
    Matrix<F,Z> B(1,Q.size(0),scalar);
    cout << "\nsolve(Matrix,'R','C')" << endl;
    Matrix<F,Z> X(1,Q.size(1),scalar);
    Q.solve(B,X,'R','C');
    X.printOn(cout);
    Matrix<F,Z> R(B.size(0),B.size(1));
    R.copy(B);
    R.gemm(Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_,'N','C');
    cout << "residual = " << endl;
    R.printOn(cout);
    Matrix<F,Z> AR(R.size(0),Q.size(1),Matrix<F,Z>::zero_);
    AR.gemm(Matrix<F,Z>::one_,R,Q,Matrix<F,Z>::zero_,'N','N');
    cout << "R * Q = " << endl;
    AR.printOn(cout);
  }

  {
    Matrix<F,Z> B(1,Q.size(0),scalar);
    cout << "\nsolve(Matrix,'R','T')" << endl;
    Matrix<F,Z> X(1,Q.size(1),scalar);
    Q.solve(B,X,'R','T');
    X.printOn(cout);
    Matrix<F,Z> R(B.size(0),B.size(1));
    R.copy(B);
    R.gemm(Matrix<F,Z>::one_,X,Q,Matrix<F,Z>::mone_,'N','T');
    cout << "residual = " << endl;
    R.printOn(cout);
    if (sizeof(F)!=sizeof(Z)) {
      Matrix<F,Z> AR(R.size(0),Q.size(1),Matrix<F,Z>::zero_);
      AR.gemm(Matrix<F,Z>::one_,R,Q,Matrix<F,Z>::zero_,'N','N');
      cout << "R * Q = " << endl;
      AR.printOn(cout);
    }
  }
}
