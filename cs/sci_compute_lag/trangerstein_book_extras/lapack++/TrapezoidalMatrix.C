#include "TrapezoidalMatrix.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
LowerTrapezoidalMatrix<F,Z>::printOn(ostream &s) const {
//a LowerTrapezoidalMatrix could be a reference to a factorization of a 
//  Matrix
  int m=this->size(0),n=this->size(1);
  s << "LowerTrapezoidalMatrix(" << m << " x " << n << ")\n";
  for (int i=0;i<m;i++) {
    for (int j=0;j<n;j++) s << this->operator()(i,j) << " ";
    s << "\n";
  }
}

template<typename F,typename Z> void testLowerTrapezoidalMatrix(F fscalar,
Z scalar) {
  LowerTrapezoidalMatrix<F,Z> *temp=0;
  {
    LowerTrapezoidalMatrix<F,Z> L1;
    cout << "after LowerTrapezoidalMatrix()" << endl;
    cout << "L1:" << endl;
    L1.printOn(cout);

    int m=3,n=2;
    cout << "\nm,n = " << m << " " << n << endl;
    {
      LowerTrapezoidalMatrix<F,Z> L2(m,n);
      cout << "after LowerTrapezoidalMatrix(m,n)" << endl;
      cout << "L2:" << endl;
      L2.printOn(cout);

      {
        LowerTrapezoidalMatrix<F,Z> L3(m,n,scalar);
        cout << "\nafter LowerTrapezoidalMatrix(m,n,scalar)" << endl;
        cout << "L3:" << endl;
        L3.printOn(cout);

        {
          cout << "\nmakeMatrix" << endl;
          Matrix<F,Z> *M=L3.makeMatrix();
          M->printOn(cout);
          delete M; M=0;
        }
      }

      L2=Vector<F,Z>::one_;
      cout << "\nafter L2=one_" << endl;
      cout << "L2:" << endl;
      L2.printOn(cout);

      {
        Matrix<F,Z> M(m,n,scalar);
        L2.copy(M);
        cout << "\nafter L2.copy(M)" << endl;
        cout << "L2:" << endl;
        L2.printOn(cout);

        {
          LowerTrapezoidalMatrix<F,Z> L3;
          L3.resize(m,n);
          cout << "\nafter L3.resize(m,n)" << endl;
          cout << "L3:" << endl;
          L3.printOn(cout);
        }

        L1.resize(L2);
        L1=scalar;
        cout << "\nL1+L2:" << endl;
        temp=L1+L2;
        temp->printOn(cout);
        delete temp; temp=0;

        cout << "\nL1+M:" << endl;
        Matrix<F,Z> *MM=0;
        MM=L1+M;
        MM->printOn(cout);
        delete MM; MM=0;

        cout << "\nM+L2:" << endl;
        MM=M+L2;
        MM->printOn(cout);
        delete MM; MM=0;

        cout << "\nL1-L2:" << endl;
        temp=L1-L2;
        temp->printOn(cout);
        delete temp; temp=0;

        cout << "\nL1-M:" << endl;
        MM=L1-M;
        MM->printOn(cout);
        delete MM; MM=0;

        cout << "\nM-L2:" << endl;
        MM=M-L2;
        MM->printOn(cout);
        delete MM; MM=0;
      }
    }

    cout << "\nL1*scalar:" << endl;
    temp=L1*scalar;
    temp->printOn(cout);
    delete temp; temp=0;

    cout << "\nL1/scalar:" << endl;
    temp=L1/scalar;
    temp->printOn(cout);
    delete temp; temp=0;

//  cout << "\ntranspose(L1)" << endl;
//  UpperTrapezoidalMatrix<F,Z> *U=L1.transpose();
//  U->printOn(cout);
//  delete U; U=0;

//  if (sizeof(F)!=sizeof(Z)) {
//    cout << "\nconjugateTranspose(L1)" << endl;
//    UpperTrapezoidalMatrix<F,Z> *U=L1.conjugateTranspose();
//    U->printOn(cout);
//    delete U; U=0;
//  }

    Vector<F,Z> *t=0;
    {
      Vector<F,Z> v(n,scalar);
      cout << "\ntrmv(v,'N'):" << endl;
      t=L1.trmv(v);
      t->printOn(cout);
      delete t; t=0;
    }

    {
      Vector<F,Z> v(m,scalar);
      cout << "\ntrmv(v,'T'):" << endl;
      t=L1.trmv(v,'T');
      t->printOn(cout);
      delete t; t=0;
    }

    Matrix<F,Z> *T=0;
    {
      Matrix<F,Z> M(n,2,scalar);
      cout << "\ntrmm(M,'L','N'):" << endl;
      T=L1.trmm(M,'L','N');
      T->printOn(cout);
      delete T; T=0;
    }

    {
      Matrix<F,Z> M(m,2,scalar);
      cout << "\ntrmm(M,'L','T'):" << endl;
      T=L1.trmm(M,'L','T');
      T->printOn(cout);
      delete T; T=0;
    }

    {
      Matrix<F,Z> M(2,m,scalar);
      cout << "\ntrmm(M,'R','N'):" << endl;
      T=L1.trmm(M,'R','N');
      T->printOn(cout);
      delete T; T=0;
    }

    {
      Matrix<F,Z> M(2,n,scalar);
      cout << "\ntrmm(M,'R','T'):" << endl;
      T=L1.trmm(M,'R','T');
      T->printOn(cout);
      delete T; T=0;
    }

    {
      Matrix<F,Z> B(2,2,Vector<F,Z>::one_);
      cout << "\nafter copyFrom, L1 = " << endl;
      L1.copyFrom(2,2,B);
      L1.printOn(cout);
    }

    cout << "\nafter scale, L1 = " << endl;
    L1.scale(LowerTrapezoidalMatrix<F,F>::one_,fscalar);
    L1.printOn(cout);

    L1.set(LowerTrapezoidalMatrix<F,Z>::one_,scalar);
    cout << "after set, L1 = " << endl;
    L1.printOn(cout);

    {
      LowerTrapezoidalMatrix<F,Z> LL(3,2);
      LL(0,0)=scalar; 
      LL(1,0)=static_cast<F>(2.)*scalar;
        LL(1,1)=static_cast<F>(4.)*scalar;
      LL(2,0)=static_cast<F>(3.)*scalar;
        LL(2,1)=static_cast<F>(4.)*scalar;
      cout << "\nLL = " << endl;
      LL.printOn(cout);

      cout << "\nLL.normFrobenius = " << LL.normFrobenius() << endl;
      cout << "LL.normInfinity = " << LL.normInfinity() << endl;
      cout << "LL.normMaxEntry = " << LL.normMaxEntry() << endl;
      cout << "LL.normOne = " << LL.normOne() << endl;

      cout << "\nLowerTrapezoidalMatrix * LowerTrapezoidalMatrix = "
           << endl;
      L1=scalar;
      {
        LowerTrapezoidalMatrix<F,Z> L2(4,3,scalar);
        LowerTrapezoidalMatrix<F,Z> *P=L2*L1;
        P->printOn(cout);
        delete P; P=0;
      }

      {
        cout << "\nLowerTrapezoidalMatrix * Matrix = " << endl;
        Matrix<F,Z> M(n,2,scalar);
        T=L1*M;
        T->printOn(cout);
        delete T; T=0;
      }

      {
        cout << "\nMatrix * LowerTrapezoidalMatrix = " << endl;
        Matrix<F,Z> M(2,m,scalar);
        T=M*L1;
        T->printOn(cout);
        delete T; T=0;
      }

      {
        cout << "\nLowerTrapezoidalMatrix * Vector = " << endl;
        Vector<F,Z>v(n,scalar);
        t=L1*v;
        t->printOn(cout);
        delete t; t=0;
      }

      {
        Vector<F,Z> w(L1.size(0),scalar);
        Vector<F,Z> t(L1.size(1));
        cout << "\nsolve('N'):" << endl;
        L1.solve(w,t);
        t.printOn(cout);
        Vector<F,Z> *Lx=L1.trmv(t,'N');
        Vector<F,Z> *r=(*Lx)-w;
        cout << "residual = " << endl;
        r->printOn(cout);
        delete r; r=0;
        delete Lx; Lx=0;
      }

      {
        Vector<F,Z> w(L1.size(1),scalar);
        Vector<F,Z> t(L1.size(0),scalar);
        cout << "\nsolve('T'):" << endl;
        L1.solve(w,t,'T');
        t.printOn(cout);
        Vector<F,Z> *Lx=L1.trmv(t,'T');
        Vector<F,Z> *r=(*Lx)-w;
        cout << "residual = " << endl;
        r->printOn(cout);
        delete r; r=0;
        delete Lx; Lx=0;
      }

      if (sizeof(F)!=sizeof(Z)) {
        Vector<F,Z> w(L1.size(1),scalar);
        Vector<F,Z> t(L1.size(0),scalar);
        cout << "\nsolve('C'):" << endl;
        L1.solve(w,t,'C');
        t.printOn(cout);
        Vector<F,Z> *Lx=L1.trmv(t,'C');
        Vector<F,Z> *r=(*Lx)-w;
        cout << "residual = " << endl;
        r->printOn(cout);
        delete r; r=0;
        delete Lx; Lx=0;
      }

      {
        Matrix<F,Z> B(LL.size(0),2);
        B(0,0)=   scalar;
          B(0,1)= static_cast<F>(2.)*scalar;
        B(1,0)=static_cast<F>(6.)*scalar;
          B(1,1)=static_cast<F>(12.)*scalar;
        B(2,0)=static_cast<F>(7.)*scalar;
          B(2,1)=static_cast<F>(14.)*scalar;
        Matrix<F,Z> X(LL.size(1),B.size(1));
        cout << "\nsolve('L','N',...)" << endl;
        LL.solve(B,X,'L','N');
        cout << "B = " << endl;
        B.printOn(cout);
        cout << "X = " << endl;
        X.printOn(cout);
        Matrix<F,Z> *LX=LL.trmm(X,'L','N');
        Matrix<F,Z> *R=(*LX)-B;
        cout << "residual = " << endl;
        R->printOn(cout);
        delete R; R=0;
        delete LX; LX=0;
      }

      {
        Matrix<F,Z> B(LL.size(1),1);
        B(0,0)=static_cast<F>(6.)*scalar;
        B(1,0)=static_cast<F>(4.)*scalar;
        {
          Matrix<F,Z> X(LL.size(0),1);
          X(2,0)=scalar;
          cout << "\nsolve('L','T',...)" << endl;
          LL.solve(B,X,'L','T');
          cout << "B = " << endl;
          B.printOn(cout);
          cout << "X = " << endl;
          X.printOn(cout);
          Matrix<F,Z> *LX=LL.trmm(X,'L','T');
          Matrix<F,Z> *R=(*LX)-B;
          cout << "residual = " << endl;
          R->printOn(cout);
          delete R; R=0;
          delete LX; LX=0;
        }

        if (sizeof(F)!=sizeof(Z)) {
          Matrix<F,Z> X(LL.size(0),1);
          X(2,0)=scalar;
          cout << "\nsolve('L','C',...)" << endl;
          LL.solve(B,X,'L','C');
          cout << "B = " << endl;
          B.printOn(cout);
          cout << "X = " << endl;
          X.printOn(cout);
          Matrix<F,Z> *LX=LL.trmm(X,'L','C');
          Matrix<F,Z> *R=(*LX)-B;
          cout << "residual = " << endl;
          R->printOn(cout);
          delete R; R=0;
          delete LX; LX=0;
        }
      }

      {
        Matrix<F,Z> B(1,LL.size(1));
        B(0,0)=static_cast<F>(12.)*scalar;
          B(0,1)=static_cast<F>(16.)*scalar;
        Matrix<F,Z> X(B.size(0),LL.size(0));
        X(0,2)=scalar;
        cout << "\nsolve('R','N',...)" << endl;
        LL.solve(B,X,'R','N');
        cout << "B = " << endl;
        B.printOn(cout);
        cout << "X = " << endl;
        X.printOn(cout);
        Matrix<F,Z> *LX=LL.trmm(X,'R','N');
        Matrix<F,Z> *R=(*LX)-B;
        cout << "residual = " << endl;
        R->printOn(cout);
        delete R; R=0;
        delete LX; LX=0;
      }

      {
        Matrix<F,Z> B(2,LL.size(0));
        B(0,0)=scalar;    
          B(1,0)=static_cast<F>(2.)*scalar;
        B(0,1)=static_cast<F>(6.)*scalar;
          B(1,1)=static_cast<F>(12.)*scalar;
        B(0,2)=static_cast<F>(7.)*scalar;
          B(1,2)=static_cast<F>(14.)*scalar;
        {
          Matrix<F,Z> X(B.size(0),LL.size(1));
          cout << "\nsolve('R','T',...)" << endl;
          LL.solve(B,X,'R','T');
          cout << "B = " << endl;
          B.printOn(cout);
          cout << "X = " << endl;
          X.printOn(cout);
          Matrix<F,Z> *LX=LL.trmm(X,'R','T');
          Matrix<F,Z> *R=(*LX)-B;
          cout << "residual = " << endl;
          R->printOn(cout);
          delete R; R=0;
          delete LX; LX=0;
        }

        if (sizeof(F)!=sizeof(Z)) {
          Matrix<F,Z> X(B.size(0),LL.size(1));
          cout << "\nsolve('R','C',...)" << endl;
          LL.solve(B,X,'R','C');
          cout << "B = " << endl;
          B.printOn(cout);
          cout << "X = " << endl;
          X.printOn(cout);
          Matrix<F,Z> *LX=LL.trmm(X,'R','C');
          Matrix<F,Z> *R=(*LX)-B;
          cout << "residual = " << endl;
          R->printOn(cout);
          delete R; R=0;
          delete LX; LX=0;
        }
      }
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
LowerTriangularMatrix<F,Z>::printOn(ostream &s) const {
  int n=this->size(0);
  s << "LowerTriangularMatrix(" << n << " x " << n << ")\n";
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) s << this->operator()(i,j) << " ";
    s << "\n";
  }
}

template<typename F,typename Z> void testLowerTriangularMatrix(
F fscalar,Z scalar) {
  {
    LowerTriangularMatrix<F,Z> L1;
    cout << "\nafter LowerTriangularMatrix()" << endl;
    cout << "L1:" << endl;
    L1.printOn(cout);
  }

  int n=2;
  cout << "n = " << n << endl;
  {
    LowerTriangularMatrix<F,Z> L2(n);
    cout << "\nafter LowerTriangularMatrix(n)" << endl;
    cout << "L2:" << endl;
    L2.printOn(cout);

    L2=scalar;
    cout << "\nafter L2=scalar:" << endl;
    L2.printOn(cout);
  }

  {
    LowerTriangularMatrix<F,Z> L3(n,scalar);
    cout << "\nafter LowerTriangularMatrix(n,scalar)" << endl;
    cout << "L3:" << endl;
    L3.printOn(cout);

    {
      cout << "\nmakeMatrix" << endl;
      Matrix<F,Z> *M=L3.makeMatrix();
      M->printOn(cout);
      delete M; M=0;
    }

    cout << "L3.reciprocalConditionNumber('I') = "
      << L3.reciprocalConditionNumber('I') << endl;
    cout << "L3.reciprocalConditionNumber('O') = "
      << L3.reciprocalConditionNumber('O') << endl;

/*
    LowerTriangularMatrix<F,Z> *T=L3.inverse();
    cout << "L3.inverse = " << endl;
    T->printOn(cout);
    delete T;
*/
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z>
void UnitLowerTrapezoidalMatrix<F,Z>::fillWith(Z d) {
  int m=this->size(0),n=this->size(1);
  for (int j=0;j<min(n,m-1);j++) {
    Z *column_j=addr(j+1,j);
    for (int i=0;i<m-j-1;i++) column_j[i]=d;
  }
}

template<typename F,typename Z> void
UnitLowerTrapezoidalMatrix<F,Z>::printOn(ostream &s) const {
  int m=this->size(0),n=this->size(1);
  s << "UnitLowerTrapezoidalMatrix(" << m << " x " << n << ")\n";
  for (int i=0;i<m;i++) {
    for (int j=0;j<n;j++) s << this->operator()(i,j) << " ";
    s << "\n";
  }
}

template<typename F,typename Z> void testUnitLowerTrapezoidalMatrix(
F fscalar,Z scalar) {
  int m=3,n=2;
  cout << "\nm,n = " << m << " " << n << endl;
  Matrix<F,Z> *T=0;
  Vector<F,Z> *t=0;
  {
    UnitLowerTrapezoidalMatrix<F,Z> L2(m,n);;
    cout << "after UnitLowerTrapezoidalMatrix(m,n)" << endl;
    cout << "L2:" << endl;
    L2.printOn(cout);
    {
      UnitLowerTrapezoidalMatrix<F,Z> L1;
      cout << "\nafter UnitLowerTrapezoidalMatrix()" << endl;
      cout << "L1:" << endl;
      L1.printOn(cout);

      UnitLowerTrapezoidalMatrix<F,Z> *L3=
        OPERATOR_NEW UnitLowerTrapezoidalMatrix<F,Z>(m,n,scalar);
      cout << "\nafter UnitLowerTrapezoidalMatrix(m,n,scalar)" << endl;
      cout << "L3:" << endl;
      L3->printOn(cout);

      {
        cout << "\nmakeMatrix" << endl;
        Matrix<F,Z> *M=L3->makeMatrix();
        M->printOn(cout);
        delete M; M=0;
      }
      delete L3; L3=0;

      L2=scalar;
      cout << "\nafter L2=scalar" << endl;
      cout << "L2:" << endl;
      L2.printOn(cout);

      LowerTrapezoidalMatrix<F,Z> *temp=0;
      {
        Matrix<F,Z> M(m,n,scalar);
        L2.copy(M);
        cout << "\nafter L2.copy(M)" << endl;
        cout << "L2:" << endl;
        L2.printOn(cout);

        L1.resize(m,n);
        cout << "\nafter L1.resize(m,n)" << endl;
        cout << "L1:" << endl;
        L1.printOn(cout);

        L1.resize(L2);
        L1=scalar;
        cout << "\nL1+L2:" << endl;
        temp=L1+L2;
        temp->printOn(cout);
        delete temp; temp=0;

        Matrix<F,Z> *MM=0;
        {
          LowerTrapezoidalMatrix<F,Z> LT(L1.size(0),L1.size(1),scalar);
          cout << "\nUnitLowerTrapezoidalMatrix + LowerTrapezoidalMatrix"
               << endl;
          temp=L1+LT;
          temp->printOn(cout);
          delete temp; temp=0;

          cout << "\nLowerTrapezoidalMatrix + UnitLowerTrapezoidalMatrix"
               << endl;
          temp=LT+L1;
          temp->printOn(cout);
          delete temp; temp=0;

          cout << "\nL1+M:" << endl;
          MM=L1+M;
          MM->printOn(cout);
          delete MM; MM=0;

          cout << "\nM+L2:" << endl;
          MM=M+L2;
          MM->printOn(cout);
          delete MM; MM=0;

          cout << "\nL1-L2:" << endl;
          temp=L1-L2;
          temp->printOn(cout);
          delete temp; temp=0;

          cout << "\nUnitLowerTrapezoidalMatrix - LowerTrapezoidalMatrix"
               << endl;
          temp=L1-LT;
          temp->printOn(cout);
          delete temp; temp=0;

          cout << "\nLowerTrapezoidalMatrix - UnitLowerTrapezoidalMatrix"
               << endl;
          temp=LT-L1;
          temp->printOn(cout);
          delete temp; temp=0;
        }

        cout << "\nL1-M:" << endl;
        MM=L1-M;
        MM->printOn(cout);
        delete MM; MM=0;

        cout << "\nM-L2:" << endl;
        MM=M-L2;
        MM->printOn(cout);
        delete MM; MM=0;
      }

      cout << "\nL1*scalar:" << endl;
      temp=L1*scalar;
      temp->printOn(cout);
      delete temp; temp=0;

      cout << "\nL1/scalar:" << endl;
      temp=L1/scalar;
      temp->printOn(cout);
      delete temp; temp=0;
    }

//  cout << "\ntranspose" << endl;
//  UnitUpperTrapezoidalMatrix<F,Z> *U=L2.transpose();
//  U->printOn(cout);
//  delete U; U=0;

//  if (sizeof(F)!=sizeof(Z)) {
//    cout << "\nconjugateTranspose" << endl;
//    UnitUpperTrapezoidalMatrix<F,Z> *U=L2.conjugateTranspose();
//    U->printOn(cout);
//    delete U; U=0;
//  }

    {
      Vector<F,Z> v(n,scalar);
      t=L2.trmv(v);
      cout << "\ntrmv('N'):" << endl;
      t->printOn(cout);
      delete t; t=0;
    }

    {
      Vector<F,Z> v(m,scalar);
      t=L2.trmv(v,'T');
      cout << "\ntrmv('T'):" << endl;
      t->printOn(cout);
      delete t; t=0;
    }

    {
      Matrix<F,Z> M(n,2,scalar);
      T=L2.trmm(M,'L','N');
      cout << "\ntrmm('L','N'):" << endl;
      T->printOn(cout);
      delete T; T=0;
    }
  
    {
      Matrix<F,Z> M(m,2,scalar);
      T=L2.trmm(M,'L','T');
      cout << "\ntrmm('L','T'):" << endl;
      T->printOn(cout);
      delete T; T=0;
    }
  
    {
      Matrix<F,Z> M(2,m,scalar);
      T=L2.trmm(M,'R','N');
      cout << "\ntrmm('R','N'):" << endl;
      T->printOn(cout);
      delete T; T=0;
    }
  
    {
      Matrix<F,Z> M(2,n,scalar);
      T=L2.trmm(M,'R','T');
      cout << "\ntrmm('R','T'):" << endl;
      T->printOn(cout);
      delete T; T=0;
    }
  
    {
      Matrix<F,Z> B(2,2,Vector<F,Z>::one_);
      L2.copyFrom(2,2,B);
      cout << "\nafter copyFrom, L2 = " << endl;
      L2.printOn(cout);
    }

    Vector<F,Z> *Lx=0;
    Vector<F,Z> *r=0;
    {
      Vector<F,Z> w(m,scalar);
      Vector<F,Z> t(n);
      cout << "\nsolve('N'):" << endl;
      L2.solve(w,t);
      t.printOn(cout);
      Lx=L2.trmv(t,'N');
      r=(*Lx)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Lx; Lx=0;
    }

    {
      Vector<F,Z> w(n,scalar);
      Vector<F,Z> t(m,scalar);
      cout << "\nsolve('T'):" << endl;
      L2.solve(w,t,'T');
      t.printOn(cout);
      Lx=L2.trmv(t,'T');
      r=(*Lx)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Lx; Lx=0;
    }

    if (sizeof(F)!=sizeof(Z)) {
      Vector<F,Z> w(n,scalar);
      Vector<F,Z> t(m,scalar);
      cout << "\nsolve('C'):" << endl;
      L2.solve(w,t,'C');
      t.printOn(cout);
      Lx=L2.trmv(t,'C');
      r=(*Lx)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Lx; Lx=0;
    }
  }

  {
    UnitLowerTrapezoidalMatrix<F,Z> LL(3,2);
    LL(1,0)=static_cast<F>(2.)*scalar;
    LL(2,0)=static_cast<F>(3.)*scalar; LL(2,1)=static_cast<F>(4.)*scalar;
    cout << "\nLL = " << endl;
    LL.printOn(cout);
  
    Matrix<F,Z> *LX=0;
    Matrix<F,Z> *R=0;
    {
      Matrix<F,Z> B(3,2);
      B(0,0)=scalar;                   B(0,1)=static_cast<F>(2.)*scalar;
      B(1,0)=static_cast<F>(6.)*scalar;B(1,1)=static_cast<F>(12.)*scalar;
      B(2,0)=static_cast<F>(7.)*scalar;B(2,1)=static_cast<F>(14.)*scalar;
      cout << "B = " << endl;
      B.printOn(cout);
      cout << "solve('L','N',...)" << endl;
      Matrix<F,Z> X(2,2);
      LL.solve(B,X,'L','N');
      cout << "X = " << endl;
      X.printOn(cout);
      LX=LL.trmm(X,'L','N');
      R=(*LX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete LX; LX=0;
    }
  
    {
      Matrix<F,Z> B(2,1);
      B(0,0)=static_cast<F>(6.)*scalar;
      B(1,0)=static_cast<F>(4.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      cout << "solve('L','T',...)" << endl;
      Matrix<F,Z> X(3,1);
      X(2,0)=scalar;
      LL.solve(B,X,'L','T');
      cout << "X = " << endl;
      X.printOn(cout);
      LX=LL.trmm(X,'L','T');
      R=(*LX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete LX; LX=0;

      if (sizeof(F)!=sizeof(Z)) {
        Matrix<F,Z> X(LL.size(0),1);
        X(2,0)=scalar;
        cout << "\nsolve('L','C',...)" << endl;
        LL.solve(B,X,'L','C');
        cout << "B = " << endl;
        B.printOn(cout);
        cout << "X = " << endl;
        X.printOn(cout);
        LX=LL.trmm(X,'L','C');
        R=(*LX)-B;
        cout << "residual = " << endl;
        R->printOn(cout);
        delete R; R=0;
        delete LX; LX=0;
      }
    }
    
    {
      Matrix<F,Z> B(1,2);
      B(0,0)=static_cast<F>(12.)*scalar;
        B(0,1)=static_cast<F>(16.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      cout << "solve('R','N',...)" << endl;
      Matrix<F,Z> X(1,3);
      X(0,2)=scalar;
      LL.solve(B,X,'R','N');
      cout << "X = " << endl;
      X.printOn(cout);
      LX=LL.trmm(X,'R','N');
      R=(*LX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete LX; LX=0;
    }

    {
      Matrix<F,Z> B(2,3);
      B(0,0)=scalar;                   B(1,0)=static_cast<F>(2.)*scalar;
      B(0,1)=static_cast<F>(6.)*scalar;B(1,1)=static_cast<F>(12.)*scalar;
      B(0,2)=static_cast<F>(7.)*scalar;B(1,2)=static_cast<F>(14.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      Matrix<F,Z> X(2,2);
      LL.solve(B,X,'R','T');
      cout << "solve('R','T',...)" << endl;
      cout << "X = " << endl;
      X.printOn(cout);
      LX=LL.trmm(X,'R','T');
      R=(*LX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete LX; LX=0;

      if (sizeof(F)!=sizeof(Z)) {
        Matrix<F,Z> X(B.size(0),LL.size(1));
        cout << "\nsolve('R','C',...)" << endl;
        LL.solve(B,X,'R','C');
        cout << "B = " << endl;
        B.printOn(cout);
        cout << "X = " << endl;
        X.printOn(cout);
        LX=LL.trmm(X,'R','C');
        R=(*LX)-B;
        cout << "residual = " << endl;
        R->printOn(cout);
        delete R; R=0;
        delete LX; LX=0;
      }
    }

    cout << "\nLL.normFrobenius = " << LL.normFrobenius() << endl;
    cout << "LL.normInfinity = " << LL.normInfinity() << endl;
    cout << "LL.normMaxEntry = " << LL.normMaxEntry() << endl;
    cout << "LL.normOne = " << LL.normOne() << endl;
  }
  
  {
    cout << "\nUnitLowerTrapezoidalMatrix * UnitLowerTrapezoidalMatrix = "
         << endl;
    UnitLowerTrapezoidalMatrix<F,Z> L1(m,n,scalar);
    LowerTrapezoidalMatrix<F,Z> *Q=0;
    {
      UnitLowerTrapezoidalMatrix<F,Z> L2(4,3,scalar);
      LowerTrapezoidalMatrix<F,Z> *P=L2*L1;
      P->printOn(cout);
      delete P; P=0;

      cout << "\nUnitLowerTrapezoidalMatrix * LowerTrapezoidalMatrix = "
           << endl;
      {
        LowerTrapezoidalMatrix<F,Z> N(3,2,scalar);
        Q=L2 * N;
        Q->printOn(cout);
        delete Q; Q=0;
      }
    }

    {
      cout << "\nLowerTrapezoidalMatrix * UnitLowerTrapezoidalMatrix = "
           << endl;
      LowerTrapezoidalMatrix<F,Z> N(4,3,scalar);
      Q=N * L1;
      Q->printOn(cout);
      delete Q; Q=0;
    }
    
    {
      cout << "\nUnitLowerTrapezoidalMatrix * Matrix = " << endl;
      Matrix<F,Z> M(n,2,scalar);
      T=L1*M;
      T->printOn(cout);
      delete T; T=0;
    }
    
    {
      cout << "\nMatrix * UnitLowerTrapezoidalMatrix = " << endl;
      Matrix<F,Z> M(2,m,scalar);
      T=M*L1;
      T->printOn(cout);
      delete T; T=0;
    }
    
    {
      cout << "\nUnitLowerTrapezoidalMatrix * Vector = " << endl;
      Vector<F,Z> v(n,scalar);
      t=L1*v;
      t->printOn(cout);
      delete t; t=0;
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
UnitLowerTriangularMatrix<F,Z>::printOn(ostream &s) const {
  int n=this->size(0);
  s << "UnitLowerTriangularMatrix(" << n << " x " << n << ")\n";
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) s << this->operator()(i,j) << " ";
    s << "\n";
  }
}

template<typename F,typename Z> void testUnitLowerTriangularMatrix(
F fscalar,Z scalar) {
  {
    UnitLowerTriangularMatrix<F,Z> L1;
    cout << "\nafter UnitLowerTriangularMatrix()" << endl;
    cout << "L1:" << endl;
    L1.printOn(cout);
  }

  int n=2;
  cout << "n = " << n << endl;
  {
    UnitLowerTriangularMatrix<F,Z> L2(n);
    cout << "\nafter UnitLowerTriangularMatrix(n)" << endl;
    cout << "L2:" << endl;
    L2.printOn(cout);

    L2=scalar;
    cout << "\nafter L2=scalar:" << endl;
    L2.printOn(cout);
  }

  {
    UnitLowerTriangularMatrix<F,Z> L3(n,scalar);
    cout << "\nafter UnitLowerTriangularMatrix(n,scalar)" << endl;
    cout << "L3:" << endl;
    L3.printOn(cout);

    {
      cout << "\nmakeMatrix" << endl;
      Matrix<F,Z> *M=L3.makeMatrix();
      M->printOn(cout);
      delete M; M=0;
    }

    cout << "L3.reciprocalConditionNumber('I') = "
      << L3.reciprocalConditionNumber('I') << endl;
    cout << "L3.reciprocalConditionNumber('O') = "
      << L3.reciprocalConditionNumber('O') << endl;

/*
    UnitLowerTriangularMatrix<F,Z> *T=L3.inverse();
    cout << "L3.inverse = " << endl;
    T->printOn(cout);
    delete T;
*/
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
UpperTrapezoidalMatrix<F,Z>::printOn(ostream &s) const {
//a upperTrapezoidalMatrix could be a reference to a factorization of a 
//  Matrix
  int m=this->size(0),n=this->size(1);
  s << "UpperTrapezoidalMatrix(" << m << " x " << n << ")\n";
  for (int i=0;i<m;i++) {
    for (int j=0;j<n;j++) s << this->operator()(i,j) << " ";
    s << "\n";
  }
}

template<typename F,typename Z> void testUpperTrapezoidalMatrix(F fscalar,
Z scalar) {
  int m=2,n=3;
  cout << "\nm,n = " << m << " " << n << endl;
  UpperTrapezoidalMatrix<F,Z> U2(m,n);;
  cout << "after UpperTrapezoidalMatrix(m,n)" << endl;
  cout << "U2:" << endl;
  U2.printOn(cout);

  {
    UpperTrapezoidalMatrix<F,Z> *T=0;
    UpperTrapezoidalMatrix<F,Z> U1;
    cout << "\nafter UpperTrapezoidalMatrix()" << endl;
    cout << "U1:" << endl;
    U1.printOn(cout);


    {
      UpperTrapezoidalMatrix<F,Z> U3(m,n,scalar);
      cout << "\nafter UpperTrapezoidalMatrix(m,n,scalar)" << endl;
      cout << "U3:" << endl;
      U3.printOn(cout);

      {
        cout << "\nmakeMatrix" << endl;
        Matrix<F,Z> *M=U3.makeMatrix();
        M->printOn(cout);
        delete M; M=0;
      }
    }

    U2=scalar;
    cout << "\nafter U2=scalar" << endl;
    cout << "U2:" << endl;
    U2.printOn(cout);

    {
      Matrix<F,Z> M(m,n,Vector<F,Z>::one_);
      U2.copy(M);
      cout << "\nafter U2.copy(M)" << endl;
      cout << "U2:" << endl;
      U2.printOn(cout);

      U1.resize(m,n);
      cout << "\nafter U1.resize(m,n)" << endl;
      cout << "U1:" << endl;
      U1.printOn(cout);

      U1=scalar;
      cout << "\nU1+U2:" << endl;
      T=U1+U2;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nU1-U2:" << endl;
      T=U1-U2;
      T->printOn(cout);
      delete T; T=0;

      {
        UpperTrapezoidalMatrix<F,Z> U3(n,n,scalar);
        SquareMatrix<F,Z> *SQ=0;
        {
          LowerTrapezoidalMatrix<F,Z> LL(n,n,scalar);
          cout << "\nUpperTrapezoidalMatrix + LowerTrapezoidalMatrix:"
               << endl;
          SQ=U3+LL;
          SQ->printOn(cout);
          delete SQ; SQ=0;

          cout << "\nLowerTrapezoidalMatrix + UpperTrapezoidalMatrix:"
               << endl;
          SQ=LL+U3;
          SQ->printOn(cout);
          delete SQ; SQ=0;

          cout << "\nUpperTrapezoidalMatrix - LowerTrapezoidalMatrix:"
               << endl;
          SQ=U3-LL;
          SQ->printOn(cout);
          delete SQ; SQ=0;

          cout << "\nLowerTrapezoidalMatrix - UpperTrapezoidalMatrix:"
               << endl;
          SQ=LL-U3;
          SQ->printOn(cout);
          delete SQ; SQ=0;
        }

        {
          UnitLowerTrapezoidalMatrix<F,Z> LL(n,n,scalar);
          cout << "\nUpperTrapezoidalMatrix + UnitLowerTrapezoidalMatrix:"
               << endl;
          SQ=U3+LL;
          SQ->printOn(cout);
          delete SQ; SQ=0;

          cout << "\nUnitLowerTrapezoidalMatrix + UpperTrapezoidalMatrix:"
               << endl;
          SQ=LL+U3;
          SQ->printOn(cout);
          delete SQ; SQ=0;

          cout << "\nUpperTrapezoidalMatrix - UnitLowerTrapezoidalMatrix:"
               << endl;
          SQ=U3-LL;
          SQ->printOn(cout);
          delete SQ; SQ=0;

          cout << "\nUnitLowerTrapezoidalMatrix - UpperTrapezoidalMatrix:"
               << endl;
          SQ=LL-U3;
          SQ->printOn(cout);
          delete SQ; SQ=0;
        }
      }

      Matrix<F,Z> *MM=0;
      cout << "\nUpperTrapezoidalMatrix + Matrix:" << endl;
      MM=U1+M;
      MM->printOn(cout);
      delete MM; MM=0;

      cout << "\nMatrix + UpperTrapezoidalMatrix:" << endl;
      MM=M+U1;
      MM->printOn(cout);
      delete MM; MM=0;

      cout << "\nUpperTrapezoidalMatrix - Matrix:" << endl;
      MM=U1-M;
      MM->printOn(cout);
      delete MM; MM=0;

      cout << "\nMatrix - UpperTrapezoidalMatrix:" << endl;
      MM=M-U1;
      MM->printOn(cout);
      delete MM; MM=0;
    }

    cout << "\nU1*scalar:" << endl;
    T=U1*scalar;
    T->printOn(cout);
    delete T; T=0;

    cout << "\nU1/scalar:" << endl;
    T=U1/scalar;
    T->printOn(cout);
    delete T; T=0;
  }

//cout << "\ntranspose(U2)" << endl;
//LowerTrapezoidalMatrix<F,Z> *L=U2.transpose();
//L->printOn(cout);
//delete L; L=0;

//if (sizeof(F)!=sizeof(Z)) {
//  cout << "\nconjugateTranspose(U2)" << endl;
//  LowerTrapezoidalMatrix<F,Z> *L=U2.conjugateTranspose();
//  L->printOn(cout);
//  delete L; L=0;
//}

  Vector<F,Z> *t=0;
  {
    Vector<F,Z> v(n,scalar);
    t=U2.trmv(v);
    cout << "\ntrmv('N'):" << endl;
    t->printOn(cout);
    delete t; t=0;
  }

  {
    Vector<F,Z> v(m,scalar);
    t=U2.trmv(v,'T');
    cout << "\ntrmv('T'):" << endl;
    t->printOn(cout);
    delete t; t=0;
  }

  if (sizeof(F)!=sizeof(Z)) {
    Vector<F,Z> v(m,scalar);
    t=U2.trmv(v,'C');
    cout << "\ntrmv('C'):" << endl;
    t->printOn(cout);
    delete t; t=0;
  }

  Matrix<F,Z> *M2=0;
  {
    Matrix<F,Z> M(n,2,scalar);
    M2=U2.trmm(M,'L','N');
    cout << "\ntrmm('L','N'):" << endl;
    M2->printOn(cout);
    delete M2; M2=0;
  }

  {
    Matrix<F,Z> M(m,2,scalar);
    M2=U2.trmm(M,'L','T');
    cout << "\ntrmm('L','T'):" << endl;
    M2->printOn(cout);
    delete M2; M2=0;
  }

  if (sizeof(F)!=sizeof(Z)) {
    Matrix<F,Z> M(m,2,scalar);
    M2=U2.trmm(M,'L','C');
    cout << "\ntrmm('L','C'):" << endl;
    M2->printOn(cout);
    delete M2; M2=0;
  }

  {
    Matrix<F,Z> M(2,m,scalar);
    M2=U2.trmm(M,'R','N');
    cout << "\ntrmm('R','N'):" << endl;
    M2->printOn(cout);
    delete M2; M2=0;
  }

  {
    Matrix<F,Z> M(2,n,scalar);
    M2=U2.trmm(M,'R','T');
    cout << "\ntrmm('R','T'):" << endl;
    M2->printOn(cout);
    delete M2; M2=0;
  }

  if (sizeof(F)!=sizeof(Z)) {
    Matrix<F,Z> M(2,n,scalar);
    M2=U2.trmm(M,'R','C');
    cout << "\ntrmm('R','C'):" << endl;
    M2->printOn(cout);
    delete M2; M2=0;
  }

  {
    Matrix<F,Z> M(m,2,Vector<F,Z>::one_);
    U2.copyFrom(m,2,M);
    cout << "\nU2.copyFrom(M)" << endl;
    U2.printOn(cout);
  }

  U2.scale(Vector<F,F>::one_,fscalar);
  cout << "\nafter scale, U2 = " << endl;
  U2.printOn(cout);

  U2.set(Vector<F,Z>::one_,scalar);
  cout << "\nafter set, U2 = " << endl;
  U2.printOn(cout);

  {
    UpperTrapezoidalMatrix<F,Z> UU(2,3);
    UU(0,0)=scalar;
      UU(0,1)=static_cast<F>(2.)*scalar;
      UU(0,2)=static_cast<F>(2.)*scalar;
    UU(1,1)=static_cast<F>(4.)*scalar;
      UU(1,2)=Vector<F,Z>::zero_;
    cout << "\nUU = " << endl;
    UU.printOn(cout);
    cout << "UU.normFrobenius = " << UU.normFrobenius() << endl;
    cout << "UU.normInfinity = " << UU.normInfinity() << endl;
    cout << "UU.normMaxEntry = " << UU.normMaxEntry() << endl;
    cout << "UU.normOne = " << UU.normOne() << endl;

    cout << "\nU2 = " << endl;
    U2.printOn(cout);
    Vector<F,Z> *Ux=0;
    Vector<F,Z> *r=0;
    {
      Vector<F,Z> w(m,scalar);
      cout << "solve('N'):" << endl;
      Vector<F,Z> t(n,scalar);
      U2.solve(w,t);
      t.printOn(cout);
      Ux=U2.trmv(t,'N');
      r=(*Ux)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Ux; Ux=0;
    }
    
    {
      Vector<F,Z> w(n,scalar);
      cout << "\nsolve('T'):" << endl;
      Vector<F,Z> t(m,scalar);
      U2.solve(w,t,'T');
      t.printOn(cout);
      Ux=U2.trmv(t,'T');
      r=(*Ux)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Ux; Ux=0;
    }

    if (sizeof(F)!=sizeof(Z)) {
      Vector<F,Z> w(n,scalar);
      cout << "\nsolve('C'):" << endl;
      Vector<F,Z> t(m,scalar);
      U2.solve(w,t,'C');
      t.printOn(cout);
      Ux=U2.trmv(t,'C');
      r=(*Ux)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Ux; Ux=0;
    }

    Matrix<F,Z> *UX=0;
    Matrix<F,Z> *R=0;
    {
      Matrix<F,Z> B(2,1);
      B(0,0)=static_cast<F>(7.)*scalar;
      B(1,0)=static_cast<F>(4.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      cout << "solve('L','N',...)" << endl;
      Matrix<F,Z> X(3,1);
      X(2,0)=scalar;
      UU.solve(B,X,'L','N');
      cout << "X = " << endl;
      X.printOn(cout);
      UX=UU.trmm(X,'L','N');
      R=(*UX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete UX; UX=0;
    }

    {
      Matrix<F,Z> B(3,1);
      B(0,0)=static_cast<F>(2.)*scalar;
      B(1,0)=static_cast<F>(8.)*scalar;
      B(2,0)=static_cast<F>(4.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      cout << "solve('L','T',...)" << endl;
      Matrix<F,Z> X(2,1);
      UU.solve(B,X,'L','T');
      cout << "X = " << endl;
      X.printOn(cout);
      UX=UU.trmm(X,'L','T');
      R=(*UX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete UX; UX=0;
    }

    if (sizeof(F)!=sizeof(Z)) {
      Matrix<F,Z> B(3,1);
      B(0,0)=static_cast<F>(2.)*scalar;
      B(1,0)=static_cast<F>(8.)*scalar;
      B(2,0)=static_cast<F>(4.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      cout << "solve('L','C',...)" << endl;
      Matrix<F,Z> X(2,1);
      UU.solve(B,X,'L','C');
      cout << "X = " << endl;
      X.printOn(cout);
      Matrix<F,Z> *UX=UU.trmm(X,'L','C');
      Matrix<F,Z> *R=(*UX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete UX; UX=0;
    }

    {
      Matrix<F,Z> B(1,3);
      B(0,0)=scalar;
        B(0,1)=static_cast<F>(6.)*scalar;
        B(0,2)=static_cast<F>(2.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      cout << "solve('R','N',...)" << endl;
      Matrix<F,Z> X(1,2);
      UU.solve(B,X,'R','N');
      cout << "X = " << endl;
      X.printOn(cout);
      UX=UU.trmm(X,'R','N');
      R=(*UX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete UX; UX=0;
    }

    {
      Matrix<F,Z> B(1,2);
      B(0,0)=static_cast<F>(3.)*scalar;B(0,1)=static_cast<F>(-4.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      cout << "solve('R','T',...)" << endl;
      Matrix<F,Z> X(1,3);
      X(0,2)=scalar;
      UU.solve(B,X,'R','T');
      cout << "X = " << endl;
      X.printOn(cout);
      UX=UU.trmm(X,'R','T');
      R=(*UX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete UX; UX=0;
    }

    if (sizeof(F)!=sizeof(Z)) {
      Matrix<F,Z> B(1,2);
      B(0,0)=static_cast<F>(3.)*scalar;B(0,1)=static_cast<F>(-4.)*scalar;
      cout << "\nB = " << endl;
      B.printOn(cout);
      cout << "solve('R','C',...)" << endl;
      Matrix<F,Z> X(1,3);
      X(0,2)=scalar;
      UU.solve(B,X,'R','C');
      cout << "X = " << endl;
      X.printOn(cout);
      UX=UU.trmm(X,'R','C');
      R=(*UX)-B;
      cout << "residual = " << endl;
      R->printOn(cout);
      delete R; R=0;
      delete UX; UX=0;
    }
  }

  Matrix<F,Z> *M=0;
  {
    UpperTrapezoidalMatrix<F,Z> U1(3,4,scalar);
    U2=scalar;
    cout << "\nUpperTrapezoidalMatrix * UpperTrapezoidalMatrix = " << endl;
    Matrix<F,Z> *T=U2 * U1;
    T->printOn(cout);
    delete T; T=0;

    LowerTrapezoidalMatrix<F,Z> L(4,2,scalar);
    cout << "\nUpperTrapezoidalMatrix * LowerTrapezoidalMatrix = " << endl;
    M=U1 * L;
    M->printOn(cout);
    delete M; M=0;
  }

  {
    UpperTrapezoidalMatrix<F,Z> U1(2,4,scalar);
    {
      LowerTrapezoidalMatrix<F,Z> L(4,3,scalar);
      cout << "\nUpperTrapezoidalMatrix * LowerTrapezoidalMatrix = "
           << endl;
      M=U1 * L;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      LowerTrapezoidalMatrix<F,Z> L(4,2,scalar);
      cout << "\nLowerTrapezoidalMatrix * UpperTrapezoidalMatrix = "
           << endl;
      M=L * U2;
      M->printOn(cout);
      delete M; M=0;
    }
  }

  {
    UpperTrapezoidalMatrix<F,Z> U1(3,4,scalar);
    {
      LowerTrapezoidalMatrix<F,Z> L(4,3,scalar);
      cout << "\nLowerTrapezoidalMatrix * UpperTrapezoidalMatrix = "
           << endl;
      M=L * U1;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      LowerTrapezoidalMatrix<F,Z> L(3,3,scalar);
      cout << "\nLowerTrapezoidalMatrix * UpperTrapezoidalMatrix = "
           << endl;
      M=L * U1;
      M->printOn(cout);
      delete M; M=0;
    }

//  {
//    LowerTrapezoidalMatrix<F,Z> L(4,2,scalar);
//    cout << "\nLowerTrapezoidalMatrix * UpperTrapezoidalMatrix = "
//         << endl;
//    M=L * U2;
//    M->printOn(cout);
//    delete M; M=0;
//  }

    {
      UnitLowerTrapezoidalMatrix<F,Z> N(4,2,scalar);
      cout << "\nUpperTrapezoidalMatrix * UnitLowerTrapezoidalMatrix = "
           << endl;
      M=U1 * N;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      UnitLowerTrapezoidalMatrix<F,Z> L(4,2,scalar);
      cout << "\nUpperTrapezoidalMatrix * UnitLowerTrapezoidalMatrix = "
           << endl;
      M=U1 * L;
      M->printOn(cout);
      delete M; M=0;
    }
  }

  {
    UpperTrapezoidalMatrix<F,Z> U1(2,4,scalar);
    {
      UnitLowerTrapezoidalMatrix<F,Z> N(4,3,scalar);
      cout << "\nUpperTrapezoidalMatrix * UnitLowerTrapezoidalMatrix = "
           << endl;
      M=U1 * N;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      UnitLowerTrapezoidalMatrix<F,Z> L(4,3,scalar);
      cout << "\nUpperTrapezoidalMatrix * UnitLowerTrapezoidalMatrix = "
           << endl;
      M=U1 * L;
      M->printOn(cout);
      delete M; M=0;
    }
  }
                
  {
    LowerTrapezoidalMatrix<F,Z> L(4,2,scalar);
    cout << "\nLowerTrapezoidalMatrix * UpperTrapezoidalMatrix = "
         << endl;
    M=L * U2;
    M->printOn(cout);
    delete M; M=0;
  }

  {
    UnitLowerTrapezoidalMatrix<F,Z> L(4,2,scalar);
    cout << "\nUnitLowerTrapezoidalMatrix * UpperTrapezoidalMatrix = "
         << endl;
    M=L * U2;
    M->printOn(cout);
    delete M; M=0;
  }

  {
    UnitLowerTrapezoidalMatrix<F,Z> N(4,2,scalar);
    cout << "\nUnitLowerTrapezoidalMatrix * UpperTrapezoidalMatrix = "
         << endl;
    M=N * U2;
    M->printOn(cout);
    delete M; M=0;
  }

  {
    Matrix<F,Z> B(3,2,scalar);
    cout << "\nUpperTrapezoidalMatrix * Matrix = " << endl;
    M=U2 * B;
    M->printOn(cout);
    delete M; M=0;

    cout << "\nMatrix * UpperTrapezoidalMatrix = " << endl;
    M=B * U2;
    M->printOn(cout);
    delete M; M=0;
  }

  {
    Vector<F,Z> v(3,scalar);
    cout << "\nUpperTrapezoidalMatrix * Vector = " << endl;
    t=U2 * v;
    t->printOn(cout);
    delete t; t=0;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
UpperTriangularMatrix<F,Z>::printOn(ostream &s) const {
  int n=this->size(0);
  s << "UpperTriangularMatrix(" << n << " x " << n << ")\n";
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) s << this->operator()(i,j) << " ";
    s << "\n";
  }
}

template<typename F,typename Z> void testUpperTriangularMatrix(
F fscalar,Z scalar) {
  {
    UpperTriangularMatrix<F,Z> U1;
    cout << "\nafter UpperTriangularMatrix()" << endl;
    cout << "U1:" << endl;
    U1.printOn(cout);
  }

  int n=2;
  cout << "n = " << n << endl;
  {
    UpperTriangularMatrix<F,Z> U2(n);
    cout << "\nafter UpperTriangularMatrix(n)" << endl;
    cout << "U2:" << endl;
    U2.printOn(cout);

    U2=scalar;
    cout << "\nafter U2=scalar:" << endl;
    U2.printOn(cout);
  }

  UpperTriangularMatrix<F,Z> U3(n,scalar);
  cout << "\nafter UpperTriangularMatrix(n,scalar)" << endl;
  cout << "U3:" << endl;
  U3.printOn(cout);

  {
    cout << "\nmakeMatrix" << endl;
    Matrix<F,Z> *M=U3.makeMatrix();
    M->printOn(cout);
    delete M; M=0;
  }

  cout << "U3.reciprocalConditionNumber('I') = "
    << U3.reciprocalConditionNumber('I') << endl;
  cout << "U3.reciprocalConditionNumber('O') = "
    << U3.reciprocalConditionNumber('O') << endl;

//UpperTriangularMatrix<F,Z> *T=U3.inverse();
//cout << "U3.inverse = " << endl;
//T->printOn(cout);
//delete T;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z>
void UnitUpperTrapezoidalMatrix<F,Z>::fillWith(Z d) {
  int m=this->size(0),n=this->size(1);
  for (int j=1;j<n;j++) {
    Z *column_j=addr(0,j);
    for (int i=0;i<min(j,m);i++) column_j[i]=d;
  }
}

template<typename F,typename Z> void
UnitUpperTrapezoidalMatrix<F,Z>::printOn(ostream &s) const {
  int m=this->size(0),n=this->size(1);
  s << "UnitUpperTrapezoidalMatrix(" << m << " x " << n << ")\n";
  for (int i=0;i<m;i++) {
    for (int j=0;j<n;j++) s << this->operator()(i,j) << " ";
    s << "\n";
  }
}

template<typename F,typename Z> void testUnitUpperTrapezoidalMatrix(
F fscalar,Z scalar) {
  int m=2,n=3;
  cout << "\nm,n = " << m << " " << n << endl;
  UpperTrapezoidalMatrix<F,Z> *T=0;
  {
    UnitUpperTrapezoidalMatrix<F,Z> U1;
    cout << "\nafter UnitUpperTrapezoidalMatrix()" << endl;
    cout << "U1:" << endl;
    U1.printOn(cout);

    {
      UnitUpperTrapezoidalMatrix<F,Z> U2(m,n);;
      cout << "after UnitUpperTrapezoidalMatrix(m,n)" << endl;
      cout << "U2:" << endl;
      U2.printOn(cout);

      {
        UnitUpperTrapezoidalMatrix<F,Z> U3(m,n,scalar);
        cout << "\nafter UnitUpperTrapezoidalMatrix(m,n,scalar)" << endl;
        cout << "U3:" << endl;
        U3.printOn(cout);

        {
          cout << "\nmakeMatrix" << endl;
          Matrix<F,Z> *M=U3.makeMatrix();
          M->printOn(cout);
          delete M; M=0;
        }
      }

      cout << "\nU2=scalar" << endl;
      U2=scalar;
      cout << "U2:" << endl;
      U2.printOn(cout);

      cout << "\nU1.copy(U2)" << endl;
      U1.resize(U2);
      U1.copy(U2);
      cout << "U1:" << endl;
      U1.printOn(cout);

      {
        cout << "\nU3.resize(m,n)" << endl;
        UnitUpperTrapezoidalMatrix<F,Z> U3;
        U3.resize(m,n);
        cout << "U3:" << endl;
        U3.printOn(cout);
      }

      {
        cout << "\nU1.copyFrom(M)" << endl;
        Matrix<F,Z> M(m,n,scalar);
        U1.copyFrom(m,n,M);
        U1.printOn(cout);
      }

      cout << "\nUnitUpperTrapezoidalMatrix + UnitUpperTrapezoidalMatrix:"
           << endl;
      U1=scalar;
      T=U1+U2;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nUnitUpperTrapezoidalMatrix - UnitUpperTrapezoidalMatrix:"
           << endl;
      T=U1-U2;
      T->printOn(cout);
      delete T; T=0;
    }

    cout << "\nUnitUpperTrapezoidalMatrix + UpperTrapezoidalMatrix:"
         << endl;
    {
      UpperTrapezoidalMatrix<F,Z> UT(m,n,scalar);
      T=U1+UT;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nUnitUpperTrapezoidalMatrix - UpperTrapezoidalMatrix:"
           << endl;
      T=U1-UT;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nUpperTrapezoidalMatrix + UnitUpperTrapezoidalMatrix:"
           << endl;
      T=UT+U1;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nUpperTrapezoidalMatrix - UnitUpperTrapezoidalMatrix:"
           << endl;
      T=U1-UT;
      T->printOn(cout);
      delete T; T=0;
    }

    {
      cout << "\nUnitUpperTrapezoidalMatrix + LowerTrapezoidalMatrix:"
           << endl;
      UnitUpperTrapezoidalMatrix<F,Z> U3(n,n,scalar);
      SquareMatrix<F,Z> *SQ=0;
      {
        LowerTrapezoidalMatrix<F,Z> LL(n,n,scalar);
        SQ=U3+LL;
        SQ->printOn(cout);
        delete SQ; SQ=0;

        cout << "\nLowerTrapezoidalMatrix + UnitUpperTrapezoidalMatrix:"
             << endl;
        SQ=LL+U3;
        SQ->printOn(cout);
        delete SQ; SQ=0;

        cout << "\nUnitUpperTrapezoidalMatrix - LowerTrapezoidalMatrix:"
             << endl;
        SQ=U3-LL;
        SQ->printOn(cout);
        delete SQ; SQ=0;

        cout << "\nLowerTrapezoidalMatrix - UnitUpperTrapezoidalMatrix:"
             << endl;
        SQ=LL-U3;
        SQ->printOn(cout);
        delete SQ; SQ=0;
      }

      cout << "\nUnitUpperTrapezoidalMatrix + UnitLowerTrapezoidalMatrix:"
           << endl;
      {
        UnitLowerTrapezoidalMatrix<F,Z> LT(n,n,scalar);
        SQ=U3+LT;
        SQ->printOn(cout);
        delete SQ; SQ=0;

        cout
          << "\nUnitLowerTrapezoidalMatrix + UnitUpperTrapezoidalMatrix:"
          << endl;
        SQ=LT+U3;
        SQ->printOn(cout);
        delete SQ; SQ=0;

        cout
          << "\nUnitUpperTrapezoidalMatrix - UnitLowerTrapezoidalMatrix:"
          << endl;
        SQ=U3-LT;
        SQ->printOn(cout);
        delete SQ; SQ=0;

        cout
          << "\nUnitLowerTrapezoidalMatrix - UnitUpperTrapezoidalMatrix:"
          << endl;
        SQ=LT-U3;
        SQ->printOn(cout);
        delete SQ; SQ=0;
      }
    }

    {
      cout << "\nUnitUpperTrapezoidalMatrix + Matrix:" << endl;
      Matrix<F,Z> M(m,n,scalar);
      Matrix<F,Z> *MM=U1+M;
      MM->printOn(cout);
      delete MM; MM=0;

      cout << "\nMatrix + UnitUpperTrapezoidalMatrix:" << endl;
      MM=M+U1;
      MM->printOn(cout);
      delete MM; MM=0;

      cout << "\nUnitUpperTrapezoidalMatrix - Matrix:" << endl;
      MM=U1-M;
      MM->printOn(cout);
      delete MM; MM=0;

      cout << "\nMatrix - UnitUpperTrapezoidalMatrix:" << endl;
      MM=M-U1;
      MM->printOn(cout);
      delete MM; MM=0;
    }

    cout << "\nU1*scalar:" << endl;
    T=U1*scalar;
    T->printOn(cout);
    delete T; T=0;

    cout << "\nU1/scalar:" << endl;
    T=U1/scalar;
    T->printOn(cout);
    delete T; T=0;
  }

  {
    UnitUpperTrapezoidalMatrix<F,Z> U2(n,n+1,scalar);
    UpperTrapezoidalMatrix<F,Z> *V=0;
    {
      cout
        << "\nUnitUpperTrapezoidalMatrix * UnitUpperTrapezoidalMatrix = "
        << endl;
      UnitUpperTrapezoidalMatrix<F,Z> U1(m,n,scalar);
      T=U1*U2;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nUnitUpperTrapezoidalMatrix * UpperTrapezoidalMatrix = "
           << endl;
      UpperTrapezoidalMatrix<F,Z> W(n,n+1,scalar);
      V=U1*W;
      V->printOn(cout);
      delete V; V=0;
    }

    cout << "\nUpperTrapezoidalMatrix * UnitUpperTrapezoidalMatrix = "
         << endl;
    UpperTrapezoidalMatrix<F,Z> W(m,n,scalar);
    V=W*U2;
    V->printOn(cout);
    delete V; V=0;
  }

  {
    cout << "\nUnitUpperTrapezoidalMatrix * UnitLowerTrapezoidalMatrix = "
         << endl;
    UnitUpperTrapezoidalMatrix<F,Z> U1(n,n+2,scalar);
    Matrix<F,Z> *M=0;
    {
      UnitLowerTrapezoidalMatrix<F,Z> LT(n+2,m,scalar);
      M=U1*LT;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      cout
        << "\nUnitLowerTrapezoidalMatrix * UnitUpperTrapezoidalMatrix = "
        << endl;
      UnitUpperTrapezoidalMatrix<F,Z> U2(m,m+2,scalar);
      {
        UnitLowerTrapezoidalMatrix<F,Z> LT(m+1,m,scalar);
        M=LT*U2;
        M->printOn(cout);
        delete M; M=0;
      }

      {
        cout << "\nUnitUpperTrapezoidalMatrix * LowerTrapezoidalMatrix = "
             << endl;
        LowerTrapezoidalMatrix<F,Z> K(n+2,m,scalar);
        M=U1*K;
        M->printOn(cout);
        delete M; M=0;
      }

      cout << "\nLowerTrapezoidalMatrix * UnitUpperTrapezoidalMatrix = "
           << endl;
      LowerTrapezoidalMatrix<F,Z> K(n,m,scalar);
      M=K*U2;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      cout << "\nUnitUpperTrapezoidalMatrix * Matrix:" << endl;
      Matrix<F,Z> N(n+2,m+1,scalar);
      M=U1*N;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      cout << "\nMatrix * UnitUpperTrapezoidalMatrix:" << endl;
      Matrix<F,Z> N(m,n,scalar);
      M=N*U1;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      cout << "\nUnitUpperTrapezoidalMatrix * Vector:" << endl;
      Vector<F,Z> b(n+2,scalar);
      Vector<F,Z> *t=U1*b;
      t->printOn(cout);
      delete t; t=0;
    }
  }

//{
//  cout << "\nU2.transpose()" << endl;
//  UnitLowerTrapezoidalMatrix<F,Z> *L=U2.transpose();
//  L->printOn(cout);
//  delete L; L=0;
//}

  {
    UnitUpperTrapezoidalMatrix<F,Z> U2(m,n,scalar);
    Vector<F,Z> *t=0;
    {
      Vector<F,Z> v(n,scalar);
      t=U2.trmv(v);
      cout << "\ntrmv('N'):" << endl;
      t->printOn(cout);
      delete t; t=0;
    }

    {
      cout << "\ntrmv('T'):" << endl;
      Vector<F,Z> v(m,scalar);
      t=U2.trmv(v,'T');
      t->printOn(cout);
      delete t; t=0;
    }

    Matrix<F,Z> *M2=0;
    {
      cout << "\ntrmm('L','N'):" << endl;
      Matrix<F,Z> M(n,2,scalar);
      M2=U2.trmm(M,'L','N');
      M2->printOn(cout);
      delete M2; M2=0;
    }

    {
      cout << "\ntrmm('L','T'):" << endl;
      Matrix<F,Z> M(m,2,scalar);
      M2=U2.trmm(M,'L','T');
      M2->printOn(cout);
      delete M2; M2=0;
    }

    {
      cout << "\ntrmm('R','N'):" << endl;
      Matrix<F,Z> M(2,m,scalar);
      M2=U2.trmm(M,'R','N');
      M2->printOn(cout);
      delete M2; M2=0;
    }

    {
      cout << "\ntrmm('R','T'):" << endl;
      Matrix<F,Z> M(2,n,scalar);
      M2=U2.trmm(M,'R','T');
      M2->printOn(cout);
      delete M2; M2=0;
    }

    Vector<F,Z> *Ux=0;
    Vector<F,Z> *r=0;
    {
      cout << "\nsolve('N'):" << endl;
      Vector<F,Z> w(m,scalar);
      Vector<F,Z> t(n,scalar);
      U2.solve(w,t);
      t.printOn(cout);
      Ux=U2.trmv(t,'N');
      r=(*Ux)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Ux; Ux=0;
    }
    
    {
      cout << "\nsolve('T'):" << endl;
      Vector<F,Z> w(n,scalar);
      Vector<F,Z> t(m,scalar);
      U2.solve(w,t,'T');
      t.printOn(cout);
      Ux=U2.trmv(t,'T');
      r=(*Ux)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Ux; Ux=0;
    }

    if (sizeof(F)!=sizeof(Z)) {
      Vector<F,Z> w(n,scalar);
      Vector<F,Z> t(m,scalar);
      U2.solve(w,t,'C');
      cout << "\nsolve('C'):" << endl;
      t.printOn(cout);
      Ux=U2.trmv(t,'C');
      r=(*Ux)-w;
      cout << "residual = " << endl;
      r->printOn(cout);
      delete r; r=0;
      delete Ux; Ux=0;
    }
  }

  UnitUpperTrapezoidalMatrix<F,Z> UU(2,3);
  UU(0,1)=static_cast<F>(2.)*scalar; UU(0,2)=static_cast<F>(3.)*scalar;
                                     UU(1,2)=static_cast<F>(4.)*scalar;
  cout << "\nUU = " << endl;
  UU.printOn(cout);

  cout << "\nUU = " << endl;
  UU.printOn(cout);
  cout << "UU.normFrobenius = " << UU.normFrobenius() << endl;
  cout << "UU.normInfinity = " << UU.normInfinity() << endl;
  cout << "UU.normMaxEntry = " << UU.normMaxEntry() << endl;
  cout << "UU.normOne = " << UU.normOne() << endl;
  
  Matrix<F,Z> *UX=0;
  Matrix<F,Z> *R=0;
  {
    Matrix<F,Z> B(2,1);
    B(0,0)=static_cast<F>(6.)*scalar; B(1,0)=static_cast<F>(4.)*scalar;
    cout << "B = " << endl;
    B.printOn(cout);
    cout << "\nsolve('L','N')" << endl;
    Matrix<F,Z> X(3,1);
    X(2,0)=scalar;
    UU.solve(B,X,'L','N');
    cout << "X = " << endl;
    X.printOn(cout);
    UX=UU.trmm(X,'L','N');
    R=(*UX)-B;
    cout << "residual = " << endl;
    R->printOn(cout);
    delete R; R=0;
    delete UX; UX=0;
  }
  
  {
    Matrix<F,Z> B(3,2);
    B(0,0)=scalar;
      B(1,0)=static_cast<F>(6.)*scalar;
      B(2,0)=static_cast<F>(7.)*scalar;
    B(0,1)=static_cast<F>(2.)*scalar;
      B(1,1)=static_cast<F>(12.)*scalar;
      B(2,1)=static_cast<F>(14.)*scalar;
    cout << "\nB = " << endl;
    B.printOn(cout);
    cout << "solve('L','T')" << endl;
    Matrix<F,Z> X(2,2);
    UU.solve(B,X,'L','T');
    cout << "X = " << endl;
    X.printOn(cout);
    UX=UU.trmm(X,'L','T');
    R=(*UX)-B;
    cout << "residual = " << endl;
    R->printOn(cout);
    delete R; R=0;
    delete UX; UX=0;
  }

  if (sizeof(F)!=sizeof(Z)) {
    Matrix<F,Z> B(3,2);
    B(0,0)=scalar;
      B(1,0)=static_cast<F>(6.)*scalar;
      B(2,0)=static_cast<F>(7.)*scalar;
    B(0,1)=static_cast<F>(2.)*scalar;
      B(1,1)=static_cast<F>(12.)*scalar;
      B(2,1)=static_cast<F>(14.)*scalar;
    cout << "\nB = " << endl;
    B.printOn(cout);
    cout << "solve('L','C')" << endl;
    Matrix<F,Z> X(2,2);
    UU.solve(B,X,'L','C');
    cout << "X = " << endl;
    X.printOn(cout);
    UX=UU.trmm(X,'L','C');
    R=(*UX)-B;
    cout << "residual = " << endl;
    R->printOn(cout);
    delete R; R=0;
    delete UX; UX=0;
  }

  {
    Matrix<F,Z> B(2,3);
    B(0,0)=scalar;                    B(1,0)=static_cast<F>(2.)*scalar;
    B(0,1)=static_cast<F>(6.)*scalar; B(1,1)=static_cast<F>(12.)*scalar;
    B(0,2)=static_cast<F>(7.)*scalar; B(1,2)=static_cast<F>(14.)*scalar;
    cout << "\nB = " << endl;
    B.printOn(cout);
    cout << "solve('R','N')" << endl;
    Matrix<F,Z> X(2,2);
    UU.solve(B,X,'R','N');
    cout << "X = " << endl;
    X.printOn(cout);
    UX=UU.trmm(X,'R','N');
    R=(*UX)-B;
    cout << "residual = " << endl;
    R->printOn(cout);
    delete R; R=0;
    delete UX; UX=0;
  }
  
  {
    Matrix<F,Z> B(1,2);
    B(0,0)=static_cast<F>(12.)*scalar; B(0,1)=static_cast<F>(16.)*scalar;
    cout << "\nB = " << endl;
    B.printOn(cout);
    cout << "solve('R','T')" << endl;
    Matrix<F,Z> X(1,3);
    X(0,2)=scalar;
    UU.solve(B,X,'R','T');
    cout << "X = " << endl;
    X.printOn(cout);
    UX=UU.trmm(X,'R','T');
    R=(*UX)-B;
    cout << "residual = " << endl;
    R->printOn(cout);
    delete R; R=0;
    delete UX; UX=0;
  }

  if (sizeof(F)!=sizeof(Z)) {
    Matrix<F,Z> B(1,2);
    B(0,0)=static_cast<F>(12.)*scalar; B(0,1)=static_cast<F>(16.)*scalar;
    cout << "\nB = " << endl;
    B.printOn(cout);
    cout << "solve('R','C')" << endl;
    Matrix<F,Z> X(1,3);
    X(0,2)=scalar;
    UU.solve(B,X,'R','C');
    cout << "X = " << endl;
    X.printOn(cout);
    UX=UU.trmm(X,'R','C');
    R=(*UX)-B;
    cout << "residual = " << endl;
    R->printOn(cout);
    delete R; R=0;
    delete UX; UX=0;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
UnitUpperTriangularMatrix<F,Z>::printOn(ostream &s) const {
  int n=this->size(0);
  s << "UnitUpperTriangularMatrix(" << n << " x " << n << ")\n";
  for (int i=0;i<n;i++) {
    for (int j=0;j<n;j++) s << this->operator()(i,j) << " ";
    s << "\n";
  }
}

template<typename F,typename Z> void testUnitUpperTriangularMatrix(
F fscalar,Z scalar) {
  UnitUpperTriangularMatrix<F,Z> *U1=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>;
  cout << "\nafter UnitUpperTriangularMatrix()" << endl;
  cout << "U1:" << endl;
  U1->printOn(cout);
  delete U1; U1=0;

  int n=2;
  cout << "n = " << n << endl;
  UnitUpperTriangularMatrix<F,Z> *U2=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n);
  cout << "\nafter UnitUpperTriangularMatrix(n)" << endl;
  cout << "U2:" << endl;
  U2->printOn(cout);

  UnitUpperTriangularMatrix<F,Z> *U3=
    OPERATOR_NEW UnitUpperTriangularMatrix<F,Z>(n,scalar);
  cout << "\nafter UnitUpperTriangularMatrix(n,scalar)" << endl;
  cout << "U3:" << endl;
  U3->printOn(cout);

  {
  cout << "\nmakeMatrix" << endl;
  Matrix<F,Z> *M=U3->makeMatrix();
  M->printOn(cout);
  delete M; M=0;
  }

  *U2=scalar;
  cout << "\nafter U2=scalar:" << endl;
  U2->printOn(cout);
  delete U2; U2=0;

  cout << "U3->reciprocalConditionNumber('I') = "
    << U3->reciprocalConditionNumber('I') << endl;
  cout << "U3->reciprocalConditionNumber('O') = "
    << U3->reciprocalConditionNumber('O') << endl;

//UnitUpperTriangularMatrix<F,Z> *T=U3->inverse();
//cout << "U3->inverse = " << endl;
//T->printOn(cout);
//delete T;
  delete U3;
}
     
// Modified from ltgmd.C by John Trangenstein, 11/8/96
//      LAPACK++ (V. 1.1 Beta)
//      (C) 1992-1996 All Rights Reserved.
