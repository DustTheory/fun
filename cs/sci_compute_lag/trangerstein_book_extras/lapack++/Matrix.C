#include "Matrix.H"

template<typename F,typename Z> void Matrix<F,Z>::printOn(ostream& s)
const {
  s << "Matrix(" << size(0) << " x " << size(1) << ")\n" ;
  for (int i=0; i< size(0); i++) {
    for (int j=0; j< size(1); j++) s << operator()(i,j) << "  ";
    s << "\n";
  }
}

template<typename F,typename Z> void testMatrix(F f,Z scalar) {
  {
    int m=3,n=2;
    Matrix<F,Z> A2(m,n);;
    cout << "\nafter Matrix(m,n)" << endl;
    cout << "A2:" << endl;
    A2.printOn(cout);
    Matrix<F,Z> *T=0;
    {
      Matrix<F,Z> A1;
      cout << "after Matrix()" << endl;
      cout << "A1:" << endl;
      A1.printOn(cout);

      {
        Matrix<F,Z> A3(m,n,scalar);
        cout << "\nafter Matrix(m,n,scalar)" << endl;
        cout << "A3:" << endl;
        A3.printOn(cout);
      }

      A2=scalar;
      cout << "\nafter A2=scalar" << endl;
      cout << "A2:" << endl;
      A2.printOn(cout);

      {
        Matrix<F,Z> A4;
        A4.resize(m,n);
        cout << "\nafter A4.resize(m,n)" << endl;
        cout << "A4:" << endl;
        A4.printOn(cout);
      }

      A1.resize(A2);
      cout << "\nafter A1.resize(A2)" << endl;
      cout << "A1:" << endl;
      A1.printOn(cout);

      cout << "\nA1.copy(A2)" << endl;
      A1.copy(A2);
      cout << "A1:" << endl;
      A1.printOn(cout);

      A2=scalar*static_cast<F>(2.);
      cout << "\nA1+A2:" << endl;
      T=A1+A2;
      T->printOn(cout);
      delete T; T=0;

//    Vector<F,Z> a(6,scalar);
//    cout << "\na+A1:" << endl;
//    Vector<F,Z> *b=a+A1; // run-time error
//    b->printOn(cout);
//    delete b; b=0;

//    cout << "\nA1+a:" << endl;
//    b=A1+a; // compiler error
//    b->printOn(cout);
//    delete b; b=0;

//    cout << "\na-A1:" << endl;
//    b=a-A1; // run-time error
//    b->printOn(cout);
//    delete b; b=0;

//    cout << "\nA1-a:" << endl;
//    b=A1-a; // compiler error
//    b->printOn(cout);
//    delete b; b=0;

      cout << "\nA1-A2:" << endl;
      T=A1-A2;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nA1*scalar :" << endl;
      T=A1*scalar;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nA1/scalar :" << endl;
      T=A1/scalar;
      T->printOn(cout);
      delete T; T=0;
    }

//  cout << "\ntranspose(A2)" << endl;
//  T=A2.transpose();
//  T->printOn(cout);
//  delete T; T=0;

//  if (sizeof(F)!=sizeof(Z)) {
//    cout << "\nconjugateTranspose(A2)" << endl;
//    T=A2.conjugateTranspose();
//    T->printOn(cout);
//    delete T; T=0;
//  }

    {
      Matrix<F,Z> AAAA(2,3);
      AAAA(0,0)=scalar;                   
        AAAA(0,1)=scalar*static_cast<F>(3.); 
        AAAA(0,2)=scalar*static_cast<F>(5.);
      AAAA(1,0)=scalar*static_cast<F>(2.);
        AAAA(1,1)=scalar*static_cast<F>(4.);
        AAAA(1,2)=scalar*static_cast<F>(6.);
      cout << "\ninterchange columns:" << endl;
      cout << "before" << endl;
      AAAA.printOn(cout);
      AAAA.interchangeColumns(0,2);
      cout << "after" << endl;
      AAAA.printOn(cout);
      cout << "\ninterchange rows:" << endl;
      AAAA.interchangeRows(0,1);
      cout << "after" << endl;
      AAAA.printOn(cout);
    }

    {
      Vector<F,Z> x(m,scalar);
      Vector<F,Z> y(n,static_cast<F>(2.)*scalar);
      cout << "\nbefore gemv, A2 = " << endl;
      A2.printOn(cout);
      cout << "x = " << endl;
      x.printOn(cout);
      cout << "y = " << endl;
      y.printOn(cout);
      A2.gemv(Vector<F,Z>::one_,x,Vector<F,Z>::one_,y,'T');
      cout << "after gemv, y = " << endl;
      y.printOn(cout);

      if (sizeof(F)==sizeof(Z)) {
        A2.ger(scalar,x,y);
        cout << "\nA2.ger(scalar,x,y)" << endl;
        A2.printOn(cout);
      } else {
        A2.gerc(scalar,x,y);
        cout << "\nA2.gerc(scalar,x,y)" << endl;
        A2.printOn(cout);
        A2.geru(scalar,x,y);
        cout << "\nA2.geru(scalar,x,y)" << endl;
        A2.printOn(cout);
      }
    }

    {
      Matrix<F,Z> C(n,n);
      C=scalar;
      C.gemm(Vector<F,Z>::one_,A2,A2,Vector<F,Z>::zero_,'T','N');
      cout << "\nafter gemm, C = " << endl;
      C.printOn(cout);
    }

    {
      Vector<F,F> r(m);
      Vector<F,F> c(n);
      F colcnd,rowcnd;
      F amax=A2.equilibrate(r,c,rowcnd,colcnd);
      cout << "\nafter equilibrate, rowcnd,colcnd,amax = "
        << rowcnd << " " << colcnd << " " << amax << endl;
      cout << "r = " << endl;
      r.printOn(cout);
      cout << "c = " << endl;
      c.printOn(cout);
    }

    {
      Matrix<F,Z> B(n,2);
      B.copyFrom('A',m,n,A2);
      cout << "\nafter copyFrom('A',...,A2), B = " << endl;
      B.printOn(cout);

      B=scalar;
      A2.copyFrom('L',m,n,B);
      cout << "\nafter copyFrom('L',...,B), A2 = " << endl;
      A2.printOn(cout);

      A2.scale('U',1,1,Vector<F,F>::one_,f);
      cout << "\nafter scale, A2 = " << endl;
      A2.printOn(cout);

      A2.set('U',static_cast<F>(2.)*scalar,scalar);
      cout << "\nafter set, A2 = " << endl;
      A2.printOn(cout);

      cout << "\nA2.lastNonzeroColumn = " << A2.lastNonzeroColumn()
           << endl;
      cout << "A2.lastNonzeroRow = " << A2.lastNonzeroRow() << endl;
      cout << "A2.normFrobenius = " << A2.normFrobenius() << endl;
      cout << "A2.normInfinity = " << A2.normInfinity() << endl;
      cout << "A2.normMaxEntry = " << A2.normMaxEntry() << endl;
      cout << "A2.normOne = " << A2.normOne() << endl;

      B=scalar;
      T=A2 * B;
      cout << "\nA2 * B:" << endl;
      T->printOn(cout);
      delete T; T=0;
    }

    Vector<F,Z> v(n,scalar);
    Vector<F,Z> *t = A2 * v;
    cout << "\nA2 * v:" << endl;
    t->printOn(cout);
    delete t; t=0;
  }

  {
    Matrix<F,Z> AA(2,2);
    AA(0,0)=scalar;                    AA(0,1)=static_cast<F>(3.)*scalar;
    AA(1,0)=static_cast<F>(2.)*scalar; AA(1,1)=static_cast<F>(4.)*scalar;
    {
      Vector<F,Z> b(AA.size(0));
      b[0]=static_cast<F>(0.)*scalar;
      b[1]=static_cast<F>(6.)*scalar;
      cout << "\nsquare system of linear equations" << endl;
      cout << "A:" << endl;
      AA.printOn(cout);
      cout << "b:" << endl;
      b.printOn(cout);
      {
        Vector<F,Z> t(AA.size(1));
        cout << "A->solve(b,x,'N'):" << endl;
        AA.solve(b,t);
        t.printOn(cout);
        {
          Vector<F,Z> res(b.size());
          res.copy(b);
          AA.gemv(Vector<F,Z>::one_,t,Vector<F,Z>::mone_,res,'N');
          cout << "residual = " << endl;
          res.printOn(cout);

          cout << "\nA->solve(b,x,'T'):" << endl;
          AA.solve(b,t,'T');
          t.printOn(cout);
          res.copy(b);
          AA.gemv(Vector<F,Z>::one_,t,Vector<F,Z>::mone_,res,'T');
          cout << "residual = " << endl;
          res.printOn(cout);

          if (sizeof(F)!=sizeof(Z)) {
            cout << "\nA->solve(b,x,'C'):" << endl;
            AA.solve(b,t,'C');
            t.printOn(cout);
            res.copy(b);
            AA.gemv(Vector<F,Z>::one_,t,Vector<F,Z>::mone_,res,'C');
            cout << "residual = " << endl;
            res.printOn(cout);
          }
        }
      }
    }

    {
      Matrix<F,Z> BB(AA.size(0),1);
      BB(0,0)=static_cast<F>(0.)*scalar;
      BB(1,0)=static_cast<F>(6.)*scalar;
      cout << "\nB:" << endl;
      BB.printOn(cout);
      Matrix<F,Z> T(AA.size(1),1);
      cout << "A->solve(B,X,'L','N'):" << endl;
      AA.solve(BB,T);
      T.printOn(cout);
      {
        Matrix<F,Z> R(BB.size(0),BB.size(1));
        R.copy(BB);
        R.gemm(Matrix<F,Z>::one_,AA,T,Matrix<F,Z>::mone_,'N','N');
        cout << "residual = " << endl;
        R.printOn(cout);

        cout << "A->solve(B,X,'L','T'):" << endl;
        AA.solve(BB,T,'L','T');
        T.printOn(cout);
        R.copy(BB);
        R.gemm(Matrix<F,Z>::one_,AA,T,Matrix<F,Z>::mone_,'T','N');
        cout << "residual = " << endl;
        R.printOn(cout);

        if (sizeof(F)!=sizeof(Z)) {
          cout << "A->solve(B,X,'L','C'):" << endl;
          AA.solve(BB,T,'L','C');
          T.printOn(cout);
          R.copy(BB);
          R.gemm(Matrix<F,Z>::one_,AA,T,Matrix<F,Z>::mone_,'C','N');
          cout << "residual = " << endl;
          R.printOn(cout);
        }
      }
    }

    {
      Matrix<F,Z> BB(1,AA.size(1));
      BB(0,0)=static_cast<F>(0.)*scalar;
      BB(0,1)=static_cast<F>(6.)*scalar;
      cout << "\nB:" << endl;
      BB.printOn(cout);
      cout << "A->solve(B,X,'R','N'):" << endl;
      Matrix<F,Z> T(BB.size(0),AA.size(0));
      AA.solve(BB,T,'R','N');
      T.printOn(cout);
      {
        Matrix<F,Z> R(BB.size(0),BB.size(1));
        R.copy(BB);
        R.gemm(Matrix<F,Z>::one_,T,AA,Matrix<F,Z>::mone_,'N','N');
        cout << "residual = " << endl;
        R.printOn(cout);

        cout << "A->solve(B,X,'R','T'):" << endl;
        Matrix<F,Z> TT(BB.size(0),AA.size(1));
        AA.solve(BB,TT,'R','T');
        TT.printOn(cout);
        R.copy(BB);
        R.gemm(Matrix<F,Z>::one_,TT,AA,Matrix<F,Z>::mone_,'N','T');
        cout << "residual = " << endl;
        R.printOn(cout);

        if (sizeof(F)!=sizeof(Z)) {
          cout << "A->solve(B,X,'R','C'):" << endl;
          AA.solve(BB,TT,'R','C');
          TT.printOn(cout);
          R.copy(BB);
          R.gemm(Matrix<F,Z>::one_,T,AA,Matrix<F,Z>::mone_,'N','C');
          cout << "residual = " << endl;
          R.printOn(cout);
        }
      }
    }
  }

  {
    Matrix<F,Z> A(3,2);
    A(0,0)=A(0,1)=scalar;
    A(1,0)=A(2,1)=static_cast<F>(1.e-9)*scalar;
    A(2,0)=A(1,1)=static_cast<F>(0.)*scalar;
    {
      Vector<F,Z> b(A.size(0));
      b[0]=scalar;
      b[1]=b[2]=static_cast<F>(0.)*scalar;
      cout << "\noverdetermined least squares" << endl;
      cout << "A:" << endl;
      A.printOn(cout);
      cout << "b:" << endl;
      b.printOn(cout);
      {
        Vector<F,Z> t(A.size(1));
        cout << "A.solve(b,x,'N'):" << endl;
        A.solve(b,t);
        t.printOn(cout);
        {
          Vector<F,Z> res(A.size(0));
          res.copy(b);
          A.gemv(Vector<F,Z>::one_,t,Vector<F,Z>::mone_,res,'N');
          cout << "residual = " << endl;
          res.printOn(cout);
          {
            Vector<F,Z> ar(A.size(1),Vector<F,Z>::zero_);
            A.gemv(Vector<F,Z>::one_,res,Vector<F,Z>::zero_,ar,'C');
            cout << "A^H * residual = " << endl;
            ar.printOn(cout);
          }
        }
      }
    }

    {
      Vector<F,Z> b(A.size(1));
      b[0]=scalar;
      b[1]=static_cast<F>(2.)*scalar;
      cout << "\nb:" << endl;
      b.printOn(cout);
      Vector<F,Z> t(A.size(0));
      cout << "A.solve(b,x,'T'):" << endl;
      if (sizeof(F)==sizeof(Z)) {
        A.solve(b,t,'T');
      } else {
        A.solve(b,t,'C');
      }
      t.printOn(cout);
      {
        Vector<F,Z> res(A.size(1));
        res.copy(b);
        A.gemv(Vector<F,Z>::one_,t,Vector<F,Z>::mone_,res,'C');
        cout << "residual = " << endl;
        res.printOn(cout);
      }
      {
        Vector<F,Z> s(A.size(1),Vector<F,Z>::zero_);
        A.solve(t,s,'N');
        Vector<F,Z> res(A.size(0));
        res.copy(t);
        A.gemv(Vector<F,Z>::one_,s,Vector<F,Z>::mone_,res,'N');
        cout << "A*s-x = " << endl;
        res.printOn(cout);
      }
    }

    {
      Matrix<F,Z> B(A.size(0),1);
      B(0,0)=scalar;
      B(1,0)=B(2,0)=static_cast<F>(0.)*scalar;
      cout << "\nB:" << endl;
      B.printOn(cout);
      cout << "A.solve(B,X,'L','N'):" << endl;
      Matrix<F,Z> T(A.size(1),B.size(1));
      A.solve(B,T);
      T.printOn(cout);
      {
        Matrix<F,Z> R(B.size(0),B.size(1));
        R.copy(B);
        R.gemm(Matrix<F,Z>::one_,A,T,Matrix<F,Z>::mone_,'N','N');
        cout << "residual = " << endl;
        R.printOn(cout);
        {
          Matrix<F,Z> AR(A.size(1),R.size(1),Matrix<F,Z>::zero_);
          AR.gemm(Matrix<F,Z>::one_,A,R,Matrix<F,Z>::zero_,'C','N');
          cout << "A^H * residual = " << endl;
          AR.printOn(cout);
        }
      }
    }

    {
      Matrix<F,Z> B(A.size(1),1);
      B(0,0)=scalar;
      B(1,0)=static_cast<F>(2.)*scalar;
      cout << "\nB:" << endl;
      B.printOn(cout);
      Matrix<F,Z> T(A.size(0),B.size(1));
      cout << "A.solve(B,X,'L','T'):" << endl;
      if (sizeof(F)==sizeof(Z)) {
        A.solve(B,T,'L','T');
      } else {
        A.solve(B,T,'L','C');
      }
      T.printOn(cout);
      {
        Matrix<F,Z> R(B.size(0),B.size(1));
        R.copy(B);
        R.gemm(Matrix<F,Z>::one_,A,T,Matrix<F,Z>::mone_,'C','N');
        cout << "residual = " << endl;
        R.printOn(cout);
      }
      {
        Matrix<F,Z> S(A.size(1),B.size(1),Matrix<F,Z>::zero_);
        A.solve(T,S,'L','N');
        {
          Matrix<F,Z> R(T.size(0),T.size(1));
          R.copy(T);
          R.gemm(Matrix<F,Z>::one_,A,S,Matrix<F,Z>::mone_,'N','N');
          cout << "A*s-x = " << endl;
          R.printOn(cout);
        }
      }
    }

    {
      Matrix<F,Z> B(1,A.size(1));
      B(0,0)=scalar;
      B(0,1)=static_cast<F>(0.)*scalar;
      cout << "\nB:" << endl;
      B.printOn(cout);
      Matrix<F,Z> T(B.size(0),A.size(0));
      cout << "A.solve(B,X,'R','N'):" << endl;
      A.solve(B,T,'R','N');
      T.printOn(cout);
      {
        Matrix<F,Z> R(B.size(0),B.size(1));
        R.copy(B);
        R.gemm(Matrix<F,Z>::one_,T,A,Matrix<F,Z>::mone_,'N','N');
        cout << "residual = " << endl;
        R.printOn(cout);
      }
      {
        Matrix<F,Z> S(B.size(0),A.size(1),Matrix<F,Z>::zero_);
        if (sizeof(F)==sizeof(Z)) {
          A.solve(T,S,'R','T');
        } else {
          A.solve(T,S,'R','C');
        }
        {
          Matrix<F,Z> R(T.size(0),T.size(1));
          R.copy(T);
          R.gemm(Matrix<F,Z>::one_,S,A,Matrix<F,Z>::mone_,'N','C');
          cout << "S*A^T - X = " << endl;
          R.printOn(cout);
        }
      }
    }

    {
      Matrix<F,Z> B(1,A.size(0));
      B(0,0)=scalar;
      B(0,1)=static_cast<F>(2.)*scalar;
      B(0,2)=static_cast<F>(4.)*scalar;
      cout << "\nB:" << endl;
      B.printOn(cout);
      Matrix<F,Z> T(B.size(0),A.size(1));
      cout << "A.solve(B,X,'R','T'):" << endl;
      if (sizeof(F)==sizeof(Z)) {
        A.solve(B,T,'R','T');
      } else {
        A.solve(B,T,'R','C');
      }
      T.printOn(cout);
      {
        Matrix<F,Z> R(B.size(0),B.size(1));
        R.copy(B);
        R.gemm(Matrix<F,Z>::one_,T,A,Matrix<F,Z>::mone_,'N','C');
        cout << "residual = " << endl;
        R.printOn(cout);
        {
          Matrix<F,Z> AR(R.size(0),A.size(1),Matrix<F,Z>::zero_);
          AR.gemm(Matrix<F,Z>::one_,R,A,Matrix<F,Z>::zero_,'N','N');
          cout << "R * A = " << endl;
          AR.printOn(cout);
        }
      }
    }
  }

  {
    Matrix<F,Z> AAA(2,3);
    AAA(0,0)=scalar;     
      AAA(0,1)=static_cast<F>(3.)*scalar;
      AAA(0,2)=static_cast<F>(4.)*scalar;
    AAA(1,0)=static_cast<F>(-6.)*scalar;
      AAA(1,1)=static_cast<F>(6.)*scalar;
      AAA(1,2)=static_cast<F>(8.)*scalar;
    {
      Vector<F,Z> b(AAA.size(0));
      b[0]=scalar;
      b[1]=static_cast<F>(2.)*scalar;
      cout << "\nunderdetermined least squares" << endl;
      cout << "A:" << endl;
      AAA.printOn(cout);
      cout << "b:" << endl;
      b.printOn(cout);
      Vector<F,Z> t(AAA.size(1));
      cout << "A->solve(b,x,'N'):" << endl;
      AAA.solve(b,t);
      t.printOn(cout);
      {
        Vector<F,Z> res(b.size());
        res.copy(b);
        AAA.gemv(Vector<F,Z>::one_,t,Vector<F,Z>::mone_,res,'N');
        cout << "residual = " << endl;
        res.printOn(cout);
      }
      {
        Vector<F,Z> s(AAA.size(0),Vector<F,Z>::zero_);
        if (sizeof(F)==sizeof(Z)) {
          AAA.solve(t,s,'T');
        } else {
          AAA.solve(t,s,'C');
        }
        {
          Vector<F,Z> res(t.size());
          res.copy(t);
          AAA.gemv(Vector<F,Z>::one_,s,Vector<F,Z>::mone_,res,'C');
          cout << "A^T*s-x = " << endl;
          res.printOn(cout);
        }
      }
    }

    {
      Vector<F,Z> b(AAA.size(1));
      b[0]=scalar;
      b[1]=b[2]=static_cast<F>(0.)*scalar;
      cout << "\nb:" << endl;
      b.printOn(cout);
      Vector<F,Z> t(AAA.size(0));
      cout << "A->solve(b,x,'T'):" << endl;
      if (sizeof(F)==sizeof(Z)) {
        AAA.solve(b,t,'T');
      } else {
        AAA.solve(b,t,'C');
      }
      t.printOn(cout);
      {
        Vector<F,Z> res(b.size());
        res.copy(b);
        AAA.gemv(Vector<F,Z>::one_,t,Vector<F,Z>::mone_,res,'C');
        cout << "residual = " << endl;
        res.printOn(cout);
        {
          Vector<F,Z> ar(AAA.size(0),Vector<F,Z>::zero_);
          AAA.gemv(Vector<F,Z>::one_,res,Vector<F,Z>::zero_,ar,'N');
          cout << "A * residual = " << endl;
          ar.printOn(cout);
        }
      }
    }

    {
      Matrix<F,Z> BBB(AAA.size(0),1);
      BBB(0,0)=scalar;
      BBB(1,0)=static_cast<F>(2.)*scalar;
      cout << "\nB:" << endl;
      BBB.printOn(cout);
      cout << "A->solve(B,X,'L','N'):" << endl;
      Matrix<F,Z> T(AAA.size(1),BBB.size(1));
      AAA.solve(BBB,T);
      T.printOn(cout);
      {
        Matrix<F,Z> R(BBB.size(0),BBB.size(1));
        R.copy(BBB);
        R.gemm(Matrix<F,Z>::one_,AAA,T,Matrix<F,Z>::mone_,'N','N');
        cout << "residual = " << endl;
        R.printOn(cout);
      }
      {
        Matrix<F,Z> S(AAA.size(0),BBB.size(1),Matrix<F,Z>::zero_);
        if (sizeof(F)==sizeof(Z)) {
          AAA.solve(T,S,'L','T');
        } else {
          AAA.solve(T,S,'L','C');
        }
        {
          Matrix<F,Z> R(T.size(0),T.size(1));
          R.copy(T);
          R.gemm(Matrix<F,Z>::one_,AAA,S,Matrix<F,Z>::mone_,'C','N');
          cout << "A^T*s-x = " << endl;
          R.printOn(cout);
        }
      }
    }

    {
      Matrix<F,Z> BBB(AAA.size(1),1);
      BBB(0,0)=scalar;
      BBB(1,0)=BBB(2,0)=static_cast<F>(0.)*scalar;
      cout << "\nB:" << endl;
      BBB.printOn(cout);
      cout << "A->solve(B,X,'L','T'):" << endl;
      Matrix<F,Z> T(AAA.size(0),BBB.size(1));
      if (sizeof(F)==sizeof(Z)) {
        AAA.solve(BBB,T,'L','T');
      } else {
        AAA.solve(BBB,T,'L','C');
      }
      T.printOn(cout);
      {
        Matrix<F,Z> R(BBB.size(0),BBB.size(1));
        R.copy(BBB);
        R.gemm(Matrix<F,Z>::one_,AAA,T,Matrix<F,Z>::mone_,'C','N');
        cout << "residual = " << endl;
        R.printOn(cout);
        {
          Matrix<F,Z> AR(AAA.size(0),R.size(1),Matrix<F,Z>::zero_);
          AR.gemm(Matrix<F,Z>::one_,AAA,R,Matrix<F,Z>::zero_,'N','N');
          cout << "A * residual = " << endl;
          AR.printOn(cout);
        }
      }
    }

    {
      Matrix<F,Z> BBB(1,3);
      BBB(0,0)=scalar;
      BBB(0,1)=BBB(0,2)=static_cast<F>(0.)*scalar;
      cout << "\nB:" << endl;
      BBB.printOn(cout);
      cout << "A->solve(B,X,'R','N'):" << endl;
      Matrix<F,Z> T(1,2);
      AAA.solve(BBB,T,'R','N');
      T.printOn(cout);
      {
        Matrix<F,Z> R(BBB.size(0),BBB.size(1));
        R.copy(BBB);
        R.gemm(Matrix<F,Z>::one_,T,AAA,Matrix<F,Z>::mone_,'N','N');
        cout << "residual = " << endl;
        R.printOn(cout);
        {
          Matrix<F,Z> AR(R.size(0),AAA.size(0),Matrix<F,Z>::zero_);
          AR.gemm(Matrix<F,Z>::one_,R,AAA,Matrix<F,Z>::zero_,'N','C');
          cout << "residual * A^T = " << endl;
          AR.printOn(cout);
        }
      }
    }

    {
      Matrix<F,Z> BBB(1,2);
      BBB(0,0)=scalar;
      BBB(0,1)=static_cast<F>(0.)*scalar;
      cout << "\nB:" << endl;
      BBB.printOn(cout);
      cout << "A->solve(B,X,'R','T'):" << endl;
      Matrix<F,Z> T(1,3);
      if (sizeof(F)==sizeof(Z)) {
        AAA.solve(BBB,T,'R','T');
      } else {
        AAA.solve(BBB,T,'R','C');
      }
      T.printOn(cout);
      {
        Matrix<F,Z> R(BBB.size(0),BBB.size(1));
        R.copy(BBB);
        R.gemm(Matrix<F,Z>::one_,T,AAA,Matrix<F,Z>::mone_,'N','C');
        cout << "residual = " << endl;
        R.printOn(cout);
      }
      {
        Matrix<F,Z> S(BBB.size(0),AAA.size(0),Matrix<F,Z>::zero_);
        AAA.solve(T,S,'R','N');
        {
          Matrix<F,Z> R(T.size(0),T.size(1));
          R.copy(T);
          R.gemm(Matrix<F,Z>::one_,S,AAA,Matrix<F,Z>::mone_,'N','N');
          cout << "S * A - X = " << endl;
          R.printOn(cout);
        }
      }
    }
  }
}
