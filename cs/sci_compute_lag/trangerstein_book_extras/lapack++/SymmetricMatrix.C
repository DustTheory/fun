#include "SymmetricMatrix.H"

template<typename F,typename Z> void SymmetricMatrix<F,Z>::printOn(ostream &s) const {
  s << "SymmetricMatrix(" << this->size(0) << " x " << this->size(1) 
    << ")\n" ;
  for (int i=0; i< this->size(0); i++) {
    for (int j=0; j< this->size(1); j++) {
      s << this->operator()(i,j) << "  ";
    }
    s << "\n";
  }
}

template<typename F,typename Z> void testSymmetricMatrix(F fscalar,
Z scalar) {
  cout << "\nscalar = " << scalar << endl;
  int n=3;
  cout << "\nn = " << n << endl;

  SquareMatrix<F,Z> *S=0;
  {
    SymmetricMatrix<F,Z> A1;
    cout << "\nafter SymmetricMatrix()" << endl;
    cout << "A1:" << endl;
    A1.printOn(cout);

    {
      SymmetricMatrix<F,Z> A2(n);;
      cout << "after SymmetricMatrix(n)" << endl;
      cout << "A2:" << endl;
      A2.printOn(cout);

      {
        SymmetricMatrix<F,Z> A3(n,scalar);
        cout << "\nafter SymmetricMatrix(n,scalar)" << endl;
        cout << "A3:" << endl;
        A3.printOn(cout);

        {
          cout << "\nmakeMatrix" << endl;
          Matrix<F,Z> *M=A3.makeMatrix();
          M->printOn(cout);
          delete M; M=0;
        }
      }

      A2=scalar;
      cout << "\nafter A2=scalar" << endl;
      cout << "A2:" << endl;
      A2.printOn(cout);

      A1.resize(n);
      A1.copy(A2);
      cout << "\nafter A1.copy(A2)" << endl;
      cout << "A1:" << endl;
      A1.printOn(cout);

      SymmetricMatrix<F,Z> *T=0;
      cout << "\nSymmetricMatrix + SymmetricMatrix:" << endl;
      T=A1+A2;
      T->printOn(cout);
      delete T; T=0;

      {
        cout << "\nSymmetricMatrix + UnitUpperTrapezoidalMatrix" << endl;
        UpperTrapezoidalMatrix<F,Z> U(n,n,scalar);
        S=A1+U;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUnitUpperTrapezoidalMatrix + SymmetricMatrix" << endl;
        S=U+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix + UpperTrapezoidalMatrix" << endl;
        UpperTrapezoidalMatrix<F,Z> U(n,n,scalar);
        S=A1+U;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUpperTrapezoidalMatrix + SymmetricMatrix" << endl;
        S=U+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix + UnitLowerTrapezoidalMatrix" << endl;
        LowerTrapezoidalMatrix<F,Z> L(n,n,scalar);
        S=A1+L;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUnitLowerTrapezoidalMatrix + SymmetricMatrix" << endl;
        S=L+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix + LowerTrapezoidalMatrix" << endl;
        LowerTrapezoidalMatrix<F,Z> L(n,n,scalar);
        S=A1+L;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nLowerTrapezoidalMatrix + SymmetricMatrix" << endl;
        S=L+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix + SquareMatrix" << endl;
        SquareMatrix<F,Z> SQ(n,scalar);
        S=A1+SQ;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nSquareMatrix + SymmetricMatrix" << endl;
        S=SQ+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix + Matrix" << endl;
        Matrix<F,Z> M(n,n,scalar);
        S=A1+M;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nMatrix + SymmetricMatrix" << endl;
        S=M+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      cout << "\nSymmetricMatrix - SymmetricMatrix:" << endl;
      T=A1-A2;
      T->printOn(cout);
      delete T; T=0;

      {
        cout << "\nSymmetricMatrix - UnitUpperTrapezoidalMatrix" << endl;
        UnitUpperTrapezoidalMatrix<F,Z> U(n,n,scalar);
        S=A1-U;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUnitUpperTrapezoidalMatrix - SymmetricMatrix" << endl;
        S=U-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix - UpperTrapezoidalMatrix" << endl;
        UpperTrapezoidalMatrix<F,Z> U(n,n,scalar);
        S=A1-U;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUpperTrapezoidalMatrix - SymmetricMatrix" << endl;
        S=U-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix - UnitLowerTrapezoidalMatrix" << endl;
        UnitLowerTrapezoidalMatrix<F,Z> L(n,n,scalar);
        S=A1-L;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUnitLowerTrapezoidalMatrix - SymmetricMatrix" << endl;
        S=L-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix - LowerTrapezoidalMatrix" << endl;
        LowerTrapezoidalMatrix<F,Z> L(n,n,scalar);
        S=A1-L;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nLowerTrapezoidalMatrix - SymmetricMatrix" << endl;
        S=L-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix - SquareMatrix" << endl;
        SquareMatrix<F,Z> SQ(n,scalar);
        S=A1-SQ;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nSquareMatrix - SymmetricMatrix" << endl;
        S=SQ-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricMatrix - Matrix" << endl;
        Matrix<F,Z> M(n,n,scalar);
        S=A1-M;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nMatrix - SymmetricMatrix" << endl;
        S=M-A2;
        S->printOn(cout);
        delete S; S=0;
      }
    }

    cout << "\nA1*scalar:" << endl;
    SquareMatrix<F,Z> *SQ=A1*scalar;
    SQ->printOn(cout);
    delete SQ; SQ=0;

    cout << "\nA1/scalar:" << endl;
    SQ=A1/scalar;
    SQ->printOn(cout);
    delete SQ; SQ=0;
  }

  {
    SymmetricMatrix<F,Z> A2(n,scalar);
    Matrix<F,Z> *M=0;
    {
      SymmetricMatrix<F,Z> A1(n,scalar);
      cout << "\nSymmetricMatrix * SymmetricMatrix" << endl;
      S=A1*A2;
      S->printOn(cout);
      delete S; S=0;
      
      {
        UnitUpperTrapezoidalMatrix<F,Z> U(n,n+1,scalar);
        cout << "\nSymmetricMatrix * UnitUpperTrapezoidalMatrix" << endl;
        M=A1*U;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UnitUpperTrapezoidalMatrix<F,Z> U(n-1,n,scalar);
        cout << "\nUnitUpperTrapezoidalMatrix * SymmetricMatrix" << endl;
        M=U*A1;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UpperTrapezoidalMatrix<F,Z> U(n,n+1,scalar);
        cout << "\nSymmetricMatrix * UpperTrapezoidalMatrix" << endl;
        M=A1*U;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UpperTrapezoidalMatrix<F,Z> U(n-1,n,scalar);
        cout << "\nUpperTrapezoidalMatrix * SymmetricMatrix" << endl;
        M=U*A1;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UnitLowerTrapezoidalMatrix<F,Z> L(n,n-1,scalar);
        cout << "\nSymmetricMatrix * UnitLowerTrapezoidalMatrix" << endl;
        M=A1*L;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UnitLowerTrapezoidalMatrix<F,Z> L(n+1,n,scalar);
        cout << "\nUnitLowerTrapezoidalMatrix * SymmetricMatrix" << endl;
        M=L*A1;
        M->printOn(cout);
        delete M; M=0;
      }
  
      {
        LowerTrapezoidalMatrix<F,Z> L(n,n-1,scalar);
        cout << "\nSymmetricMatrix * LowerTrapezoidalMatrix" << endl;
        M=A1*L;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        LowerTrapezoidalMatrix<F,Z> L(n+1,n,scalar);
        cout << "\nLowerTrapezoidalMatrix * SymmetricMatrix" << endl;
        M=L*A1;
        M->printOn(cout);
        delete M; M=0;
      }

      {
        SquareMatrix<F,Z> S(n);
        for (int j=0;j<n;j++) {
          for (int i=0;i<n;i++) S(i,j)=static_cast<F>(i+2*j)*scalar;
        }
        cout << "\nSymmetricMatrix * SquareMatrix" << endl;
        SquareMatrix<F,Z> *S2=A1*S;
        S2->printOn(cout);
        delete S2; S2=0;

        cout << "\nSquareMatrix * SymmetricMatrix" << endl;
        S2=S*A1;
        S2->printOn(cout);
        delete S2; S2=0;
      }

      {
        Matrix<F,Z> N(n,n+1,scalar);
        cout << "\nSymmetricMatrix *  Matrix" << endl;
        M=A1*N;
        M->printOn(cout);
        delete M; M=0;
      }

      {
        Matrix<F,Z> N(n+1,n,scalar);
        cout << "\nMatrix * SymmetricMatrix" << endl;
        M=N*A1;
        M->printOn(cout);
        delete M; M=0;
      }

      {
        Vector<F,Z> v(n,scalar);
        cout << "\nSymmetricMatrix * Vector" << endl;
        Vector<F,Z> *t=A1*v;
        t->printOn(cout);
        delete t; t=0;
      }
    }

    {
      Vector<F,Z> x(n,scalar);
      {
        Vector<F,Z> y(n,scalar);
        cout << "\nsymv:" << endl;
        A2.symv(scalar,x,SymmetricMatrix<F,Z>::one_,y);
        y.printOn(cout);
      }

      {
        SymmetricMatrix<F,Z> A1(n,scalar);
        cout << "\nsyr:" << endl;
        A1.syr(scalar,x);
        A1.printOn(cout);
      }

      {
        Vector<F,Z> y(n,scalar);
        SymmetricMatrix<F,Z> A1(n,scalar);
        cout << "\nsyr2:" << endl;
        A1.syr2(scalar,x,y);
        A1.printOn(cout);
      }
    }

    {
      Matrix<F,Z> B(n,1,scalar);
      Matrix<F,Z> C(n,1,scalar);
      cout << "\nsymm('L')" << endl;
      A2.symm(scalar,B,SymmetricMatrix<F,Z>::one_,C);
      C.printOn(cout);
    }

    {
      Matrix<F,Z> B(1,n,scalar);
      Matrix<F,Z> C(1,n,scalar);
      cout << "\nsymm('R')" << endl;
      A2.symm(scalar,B,SymmetricMatrix<F,Z>::one_,C,'R');
      C.printOn(cout);
    }

    {
      Matrix<F,Z> B(n,1,scalar);
      cout << "\nsyrk('N')" << endl;
      A2.syrk(scalar,B,SymmetricMatrix<F,Z>::one_);
      A2.printOn(cout);
    }
  }

  {
    SymmetricMatrix<F,Z> A2(n,scalar);
    Matrix<F,Z> B(1,n,scalar);
    cout << "\nsyrk('C')" << endl;
    A2.syrk(scalar,B,SymmetricMatrix<F,Z>::one_,'C');
    A2.printOn(cout);
  }

  {
    SymmetricMatrix<F,Z> A2(n,scalar);
    Matrix<F,Z> B(n,1,scalar);
    Matrix<F,Z> C(n,1,scalar);
    cout << "\nsyr2k('N')" << endl;
    A2.syr2k(scalar,B,C,SymmetricMatrix<F,Z>::one_);
    A2.printOn(cout);
  }

  {
    SymmetricMatrix<F,Z> A2(n,scalar);
    Matrix<F,Z> B(1,n,scalar);
    Matrix<F,Z> C(1,n,scalar);
    cout << "\nsyr2k('C')" << endl;
    A2.syr2k(scalar,B,C,SymmetricMatrix<F,Z>::one_,'C');
    A2.printOn(cout);
  }

  {
    SymmetricMatrix<F,Z> A1(3);
    A1(0,0)=fscalar;
    A1(1,0)=static_cast<F>(2.)*scalar;
      A1(1,1)=static_cast<F>(3.)*fscalar;
    A1(2,0)=static_cast<F>(4.)*scalar;
      A1(2,1)=static_cast<F>(5.)*scalar;
      A1(2,2)=static_cast<F>(6.)*fscalar;
    cout << "\nequilibrate" << endl;
    A1.printOn(cout);
    {
      Vector<F,F> s(3);
      F scond;
      F amax=A1.equilibrate(s,scond);
      cout << "amax,scond = " << amax << " " << scond << endl;
      s.printOn(cout);
    }

    cout << "\nset" << endl;
    A1.set(SymmetricMatrix<F,Z>::one_,scalar);
    A1.printOn(cout);

    cout << "\nnormFrobenius = " << A1.normFrobenius() << endl;
    cout << "normInfinity = " << A1.normInfinity() << endl;
    cout << "normMaxEntry = " << A1.normMaxEntry() << endl;
    cout << "normOne = " << A1.normOne() << endl;

    cout << "\nreciprocalConditionNumber = "
      << A1.reciprocalConditionNumber() << endl;

//  cout << "\ninverse" << endl;
//  A3=A1.inverse();
//  A3->printOn(cout);
//  delete A3; A3=0;

    A1(0,0)=fscalar;
    A1(1,0)=static_cast<F>(2.)*scalar;
      A1(1,1)=static_cast<F>(3.)*fscalar;
    A1(2,0)=static_cast<F>(4.)*scalar;
      A1(2,1)=static_cast<F>(5.)*scalar;
      A1(2,2)=static_cast<F>(6.)*fscalar;
    cout << "\nA1" << endl;
    A1.printOn(cout);
    OrthogonalMatrix<F,Z> *V=OPERATOR_NEW OrthogonalMatrix<F,Z>(3,3);
    cout << "\neigenvalues" << endl;
    Vector<F,F> *lambda=A1.eigenvalues(V);
    lambda->printOn(cout);
    V->printOn(cout);
    SquareMatrix<F,Z> RR(3);
    for (int j=0;j<3;j++) {
      for (int i=0;i<3;i++) {
        RR(i,j)=-(*V)(i,j)*(*lambda)[j];
        for (int k=0;k<3;k++) {
          RR(i,j)+=(*static_cast<const SymmetricMatrix<F,Z>*>(&A1))(i,k)
            *(*V)(k,j);
        }
      }
    }
    cout << "\nAA * V - V * lambda" << endl;
    RR.printOn(cout);
    delete lambda; lambda=0;
    delete V; V=0;
  }

  {
    SymmetricMatrix<F,Z> AA(2);
    {
      Vector<F,Z> bb(2);
      AA(0,0)=static_cast<F>(2.)*scalar; bb[0]=static_cast<F>(0.)*scalar;
      AA(1,0)=static_cast<F>(3.)*scalar;
        AA(1,1)=static_cast<F>(4.)*scalar;
        bb[1]=static_cast<F>(6.)*scalar;
      cout << "\nsymmetric indefinite system of linear equations" << endl;
      cout << "A:" << endl;
      AA.printOn(cout);
      cout << "b:" << endl;
      bb.printOn(cout);
      cout << "solve(b,x):" << endl;
      Vector<F,Z> x(2);
      AA.solve(bb,x);
      x.printOn(cout);
      Vector<F,Z> r(bb.size());
      r.copy(bb);
      AA.symv(Vector<F,Z>::one_,x,Vector<F,Z>::mone_,r);
      cout << "residual = " << endl;
      r.printOn(cout);
    }

    {
      Matrix<F,Z> B(2,1);
      B(0,0)=static_cast<F>(0.)*scalar;
      B(1,0)=static_cast<F>(6.)*scalar;
      cout << "\nsolve('L')" << endl;
      Matrix<F,Z> X(2,1);
      AA.solve(B,X);
      X.printOn(cout);
      Matrix<F,Z> R(B.size(0),B.size(1));
      R.copy(B);
      AA.symm(Matrix<F,Z>::one_,X,Matrix<F,Z>::mone_,R,'L');
      cout << "residual = " << endl;
      R.printOn(cout);
    }

    {
      Matrix<F,Z> B(1,2);
      B(0,0)=static_cast<F>(0.)*scalar;
      B(0,1)=static_cast<F>(6.)*scalar;
      cout << "\nsolve('R')" << endl;
      Matrix<F,Z> X(1,2);
      AA.solve(B,X,'R');
      X.printOn(cout);
      Matrix<F,Z> R(B.size(0),B.size(1));
      R.copy(B);
      AA.symm(Matrix<F,Z>::one_,X,Matrix<F,Z>::mone_,R,'R');
      cout << "residual = " << endl;
      R.printOn(cout);
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
template<typename F,typename Z> void
SymmetricPositiveMatrix<F,Z>::printOn(ostream &s) const {
  s << "SymmetricPositiveMatrix(" << this->size(0) << " x " 
    << this->size(1) << ")\n" ;
  for (int i=0; i< this->size(0); i++) {
    for (int j=0; j< this->size(1); j++) {
      s << this->operator()(i,j) << "  ";
    }
    s << "\n";
  }
}

template<typename F,typename Z> void testSymmetricPositiveMatrix(F fscalar,
Z scalar) {
  cout << "\nscalar = " << scalar << endl;
  int n=3;
  cout << "\nn = " << n << endl;

  SquareMatrix<F,Z> *S=0;
  {
    SymmetricPositiveMatrix<F,Z> A1;
    cout << "\nafter SymmetricPositiveMatrix()" << endl;
    cout << "A1:" << endl;
    A1.printOn(cout);

    {
      SymmetricPositiveMatrix<F,Z> A2(n);;
      cout << "after SymmetricPositiveMatrix(n)" << endl;
      cout << "A2:" << endl;
      A2.printOn(cout);

      {
        SymmetricPositiveMatrix<F,Z> A3(n,scalar);
        cout << "\nafter SymmetricPositiveMatrix(n,scalar)" << endl;
        cout << "A3:" << endl;
        A3.printOn(cout);

        {
          cout << "\nmakeMatrix" << endl;
          Matrix<F,Z> *M=A3.makeMatrix();
          M->printOn(cout);
          delete M; M=0;
        }
      }

      A2=scalar;
      cout << "\nafter A2=scalar" << endl;
      cout << "A2:" << endl;
      A2.printOn(cout);

      A1.resize(n);
      A1.copy(A2);
      cout << "\nafter A1.copy(A2)" << endl;
      cout << "A1:" << endl;
      A1.printOn(cout);

      {
        cout << "\nSymmetricPositiveMatrix + SymmetricPositiveMatrix:"
             << endl;
        SymmetricPositiveMatrix<F,Z> *T=A1+A2;
        T->printOn(cout);
        delete T; T=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix + UnitUpperTrapezoidalMatrix"
             << endl;
        UpperTrapezoidalMatrix<F,Z> U(n,n,scalar);
        S=A1+U;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUnitUpperTrapezoidalMatrix + SymmetricPositiveMatrix"
             << endl;
        S=U+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix + UpperTrapezoidalMatrix"
             << endl;
        UpperTrapezoidalMatrix<F,Z> U(n,n,scalar);
        S=A1+U;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUpperTrapezoidalMatrix + SymmetricPositiveMatrix"
             << endl;
        S=U+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix + UnitLowerTrapezoidalMatrix"
             << endl;
        LowerTrapezoidalMatrix<F,Z> L(n,n,scalar);
        S=A1+L;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUnitLowerTrapezoidalMatrix + SymmetricPositiveMatrix"
             << endl;
        S=L+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix + LowerTrapezoidalMatrix"
             << endl;
        LowerTrapezoidalMatrix<F,Z> L(n,n,scalar);
        S=A1+L;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nLowerTrapezoidalMatrix + SymmetricPositiveMatrix"
             << endl;
        S=L+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix + SquareMatrix" << endl;
        SquareMatrix<F,Z> SQ(n,scalar);
        S=A1+SQ;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nSquareMatrix + SymmetricPositiveMatrix" << endl;
        S=SQ+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix + Matrix" << endl;
        Matrix<F,Z> M(n,n,scalar);
        S=A1+M;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nMatrix + SymmetricPositiveMatrix" << endl;
        S=M+A2;
        S->printOn(cout);
        delete S; S=0;
      }

      cout << "\nSymmetricPositiveMatrix - SymmetricPositiveMatrix:"
           << endl;
      S=A1-A2;
      S->printOn(cout);
      delete S; S=0;

      {
        cout << "\nSymmetricPositiveMatrix - UnitUpperTrapezoidalMatrix"
             << endl;
        UnitUpperTrapezoidalMatrix<F,Z> U(n,n,scalar);
        S=A1-U;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUnitUpperTrapezoidalMatrix - SymmetricPositiveMatrix"
             << endl;
        S=U-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix - UpperTrapezoidalMatrix"
             << endl;
        UpperTrapezoidalMatrix<F,Z> U(n,n,scalar);
        S=A1-U;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUpperTrapezoidalMatrix - SymmetricPositiveMatrix"
             << endl;
        S=U-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix - UnitLowerTrapezoidalMatrix"
             << endl;
        UnitLowerTrapezoidalMatrix<F,Z> L(n,n,scalar);
        S=A1-L;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nUnitLowerTrapezoidalMatrix - SymmetricPositiveMatrix"
             << endl;
        S=L-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix - LowerTrapezoidalMatrix"
             << endl;
        LowerTrapezoidalMatrix<F,Z> L(n,n,scalar);
        S=A1-L;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nLowerTrapezoidalMatrix - SymmetricPositiveMatrix"
             << endl;
        S=L-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix - SquareMatrix" << endl;
        SquareMatrix<F,Z> SQ(n,scalar);
        S=A1-SQ;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nSquareMatrix - SymmetricPositiveMatrix" << endl;
        S=SQ-A2;
        S->printOn(cout);
        delete S; S=0;
      }

      {
        cout << "\nSymmetricPositiveMatrix - Matrix" << endl;
        Matrix<F,Z> M(n,n,scalar);
        S=A1-M;
        S->printOn(cout);
        delete S; S=0;

        cout << "\nMatrix - SymmetricPositiveMatrix" << endl;
        S=M-A2;
        S->printOn(cout);
        delete S; S=0;
      }
    }

    cout << "\nA1*scalar:" << endl;
    SquareMatrix<F,Z> *SQ=A1*scalar;
    SQ->printOn(cout);
    delete SQ; SQ=0;

    cout << "\nA1/scalar:" << endl;
    SQ=A1/scalar;
    SQ->printOn(cout);
    delete SQ; SQ=0;
  }

  {
    SymmetricPositiveMatrix<F,Z> A2(n,scalar);
    {
      SymmetricPositiveMatrix<F,Z> A1(n,scalar);
      cout << "\nSymmetricPositiveMatrix * SymmetricPositiveMatrix"
           << endl;
      S=A1*A2;
      S->printOn(cout);
      delete S; S=0;
      
      Matrix<F,Z> *M=0;
      {
        UnitUpperTrapezoidalMatrix<F,Z> U(n,n+1,scalar);
        cout << "\nSymmetricPositiveMatrix * UnitUpperTrapezoidalMatrix"
             << endl;
        M=A1*U;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UnitUpperTrapezoidalMatrix<F,Z> U(n-1,n,scalar);
        cout << "\nUnitUpperTrapezoidalMatrix * SymmetricPositiveMatrix"
             << endl;
        M=U*A1;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UpperTrapezoidalMatrix<F,Z> U(n,n+1,scalar);
        cout << "\nSymmetricPositiveMatrix * UpperTrapezoidalMatrix"
             << endl;
        M=A1*U;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UpperTrapezoidalMatrix<F,Z> U(n-1,n,scalar);
        cout << "\nUpperTrapezoidalMatrix * SymmetricPositiveMatrix"
             << endl;
        M=U*A1;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UnitLowerTrapezoidalMatrix<F,Z> L(n,n-1,scalar);
        cout << "\nSymmetricPositiveMatrix * UnitLowerTrapezoidalMatrix"
             << endl;
        M=A1*L;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        UnitLowerTrapezoidalMatrix<F,Z> L(n+1,n,scalar);
        cout << "\nUnitLowerTrapezoidalMatrix * SymmetricPositiveMatrix"
             << endl;
        M=L*A1;
        M->printOn(cout);
        delete M; M=0;
      }
  
      {
        LowerTrapezoidalMatrix<F,Z> L(n,n-1,scalar);
        cout << "\nSymmetricPositiveMatrix * LowerTrapezoidalMatrix"
             << endl;
        M=A1*L;
        M->printOn(cout);
        delete M; M=0;
      }
      
      {
        LowerTrapezoidalMatrix<F,Z> L(n+1,n,scalar);
        cout << "\nLowerTrapezoidalMatrix * SymmetricPositiveMatrix"
             << endl;
        M=L*A1;
        M->printOn(cout);
        delete M; M=0;
      }

      {
        SquareMatrix<F,Z> S(n);
        for (int j=0;j<n;j++) {
          for (int i=0;i<n;i++) S(i,j)=static_cast<F>(i+2*j)*scalar;
        }
        cout << "\nSymmetricPositiveMatrix * SquareMatrix" << endl;
        SquareMatrix<F,Z> *S2=A1*S;
        S2->printOn(cout);
        delete S2; S2=0;

        cout << "\nSquareMatrix * SymmetricPositiveMatrix" << endl;
        S2=S*A1;
        S2->printOn(cout);
        delete S2; S2=0;
      }

      {
        Matrix<F,Z> N(n,n+1,scalar);
        cout << "\nSymmetricPositiveMatrix *  Matrix" << endl;
        M=A1*N;
        M->printOn(cout);
        delete M; M=0;
      }

      {
        Matrix<F,Z> N(n+1,n,scalar);
        cout << "\nMatrix * SymmetricPositiveMatrix" << endl;
        M=N*A1;
        M->printOn(cout);
        delete M; M=0;
      }

      {
        Vector<F,Z> v(n,scalar);
        cout << "\nSymmetricPositiveMatrix * Vector" << endl;
        Vector<F,Z> *t=A1*v;
        t->printOn(cout);
        delete t; t=0;
      }
    }

    {
      Vector<F,Z> x(n,scalar);
      {
        Vector<F,Z> y(n,scalar);
        cout << "\nsymv:" << endl;
        A2.symv(scalar,x,SymmetricPositiveMatrix<F,Z>::one_,y);
        y.printOn(cout);
      }

      {
        SymmetricPositiveMatrix<F,Z> A1(n,scalar);
        cout << "\nsyr:" << endl;
        A1.syr(scalar,x);
        A1.printOn(cout);
      }

      {
        Vector<F,Z> y(n,scalar);
        SymmetricPositiveMatrix<F,Z> A1(n,scalar);
        cout << "\nsyr2:" << endl;
        A1.syr2(scalar,x,y);
        A1.printOn(cout);
      }
    }

    {
      Matrix<F,Z> B(n,1,scalar);
      Matrix<F,Z> C(n,1,scalar);
      cout << "\nsymm('L')" << endl;
      A2.symm(scalar,B,SymmetricPositiveMatrix<F,Z>::one_,C);
      C.printOn(cout);
    }

    {
      Matrix<F,Z> B(1,n,scalar);
      Matrix<F,Z> C(1,n,scalar);
      cout << "\nsymm('R')" << endl;
      A2.symm(scalar,B,SymmetricPositiveMatrix<F,Z>::one_,C,'R');
      C.printOn(cout);
    }

    {
      Matrix<F,Z> B(n,1,scalar);
      cout << "\nsyrk('N')" << endl;
      A2.syrk(scalar,B,SymmetricPositiveMatrix<F,Z>::one_);
      A2.printOn(cout);
    }
  }

  {
    SymmetricPositiveMatrix<F,Z> A2(n,scalar);
    Matrix<F,Z> B(1,n,scalar);
    cout << "\nsyrk('C')" << endl;
    A2.syrk(scalar,B,SymmetricPositiveMatrix<F,Z>::one_,'C');
    A2.printOn(cout);
  }

  {
    SymmetricPositiveMatrix<F,Z> A2(n,scalar);
    Matrix<F,Z> B(n,1,scalar);
    Matrix<F,Z> C(n,1,scalar);
    cout << "\nsyr2k('N')" << endl;
    A2.syr2k(scalar,B,C,SymmetricPositiveMatrix<F,Z>::one_);
    A2.printOn(cout);
  }

  {
    SymmetricPositiveMatrix<F,Z> A2(n,scalar);
    Matrix<F,Z> B(1,n,scalar);
    Matrix<F,Z> C(1,n,scalar);
    cout << "\nsyr2k('C')" << endl;
    A2.syr2k(scalar,B,C,SymmetricPositiveMatrix<F,Z>::one_,'C');
    A2.printOn(cout);
  }

  {
    SymmetricPositiveMatrix<F,Z> A1(3);
    A1(0,0)=fscalar;
    A1(1,0)=static_cast<F>(2.)*scalar;
      A1(1,1)=static_cast<F>(3.)*fscalar;
    A1(2,0)=static_cast<F>(4.)*scalar;
      A1(2,1)=static_cast<F>(5.)*scalar;
      A1(2,2)=static_cast<F>(6.)*fscalar;
    cout << "\nequilibrate" << endl;
    A1.printOn(cout);
    {
      Vector<F,F> s(3);
      F scond;
      F amax=A1.equilibrate(s,scond);
      cout << "amax,scond = " << amax << " " << scond << endl;
      s.printOn(cout);
    }

    cout << "\nset" << endl;
    A1.set(SymmetricPositiveMatrix<F,Z>::one_,scalar);
    A1.printOn(cout);

    cout << "\nnormFrobenius = " << A1.normFrobenius() << endl;
    cout << "normInfinity = " << A1.normInfinity() << endl;
    cout << "normMaxEntry = " << A1.normMaxEntry() << endl;
    cout << "normOne = " << A1.normOne() << endl;

    cout << "\nreciprocalConditionNumber = "
      << A1.reciprocalConditionNumber() << endl;

//  cout << "\ninverse" << endl;
//  A3=A1.inverse();
//  A3->printOn(cout);
//  delete A3; A3=0;

    A1(0,0)=fscalar;
    A1(1,0)=static_cast<F>(2.)*scalar;
      A1(1,1)=static_cast<F>(3.)*fscalar;
    A1(2,0)=static_cast<F>(4.)*scalar;
      A1(2,1)=static_cast<F>(5.)*scalar;
      A1(2,2)=static_cast<F>(6.)*fscalar;
    OrthogonalMatrix<F,Z> *V=OPERATOR_NEW OrthogonalMatrix<F,Z>(3,3);
    cout << "\neigenvalues" << endl;
    Vector<F,F> *lambda=A1.eigenvalues(V);
    lambda->printOn(cout);
    V->printOn(cout);
    SquareMatrix<F,Z> RR(3);
    for (int j=0;j<3;j++) {
      for (int i=0;i<3;i++) {
        RR(i,j)=-(*V)(i,j)*(*lambda)[j];
        for (int k=0;k<3;k++) {
          RR(i,j)+=(*static_cast<const SymmetricMatrix<F,Z>*>(&A1))(i,k)
            *(*V)(k,j);
        }
      }
    }
    cout << "\nAA * V - V * lambda" << endl;
    RR.printOn(cout);
    delete lambda; lambda=0;
    delete V; V=0;
  }

  {
    SymmetricPositiveMatrix<F,Z> AA(2);
    {
      Vector<F,Z> bb(2);
      AA(0,0)=static_cast<F>(2.)*abs(scalar);
        bb[0]=static_cast<F>(0.)*scalar;
      AA(1,0)=static_cast<F>(3.)*scalar;
        AA(1,1)=static_cast<F>(5.)*abs(scalar);
        bb[1]=static_cast<F>(6.)*scalar;
      cout << "\nsymmetric system of linear equations" << endl;
      cout << "A:" << endl;
      AA.printOn(cout);
      cout << "b:" << endl;
      bb.printOn(cout);
      cout << "solve(b,x):" << endl;
      Vector<F,Z> x(2);
      AA.solve(bb,x);
      x.printOn(cout);
      Vector<F,Z> r(bb.size());
      r.copy(bb);
      AA.symv(Vector<F,Z>::one_,x,Vector<F,Z>::mone_,r);
      cout << "residual = " << endl;
      r.printOn(cout);
    }

    {
      Matrix<F,Z> B(2,1);
      B(0,0)=static_cast<F>(0.)*scalar;
      B(1,0)=static_cast<F>(6.)*scalar;
      cout << "\nsolve('L')" << endl;
      Matrix<F,Z> X(2,1);
      AA.solve(B,X);
      X.printOn(cout);
      Matrix<F,Z> R(B.size(0),B.size(1));
      R.copy(B);
      AA.symm(Matrix<F,Z>::one_,X,Matrix<F,Z>::mone_,R,'L');
      cout << "residual = " << endl;
      R.printOn(cout);
    }

    {
      Matrix<F,Z> B(1,2);
      B(0,0)=static_cast<F>(0.)*scalar;
      B(0,1)=static_cast<F>(6.)*scalar;
      cout << "\nsolve('R')" << endl;
      Matrix<F,Z> X(1,2);
      AA.solve(B,X,'R');
      X.printOn(cout);
      Matrix<F,Z> R(B.size(0),B.size(1));
      R.copy(B);
      AA.symm(Matrix<F,Z>::one_,X,Matrix<F,Z>::mone_,R,'R');
      cout << "residual = " << endl;
      R.printOn(cout);
    }
  }
}

// Modified from LaPack++ gmc.C by John Trangenstein, 11/7/96 and 5/27/12

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
