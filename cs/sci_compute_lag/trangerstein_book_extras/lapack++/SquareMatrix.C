#include "SquareMatrix.H"

template<typename F,typename Z> void SquareMatrix<F,Z>::printOn(ostream &s)
const {
  s << "SquareMatrix(" << this->size(0) << " x " << this->size(1) 
    << ")\n" ;
  for (int i=0; i< this->size(0); i++) {
    for (int j=0; j< this->size(1); j++) {
      s << this->operator()(i,j) << "  ";
    }
    s << "\n";
  }
}

template<typename F,typename Z> void testSquareMatrix(F fscalar,Z scalar) {
  cout << "\nscalar = " << scalar << endl;

  int n=3;
  SquareMatrix<F,Z> A2(n);;
  cout << "\nafter SquareMatrix(n)" << endl;
  cout << "A2:" << endl;
  A2.printOn(cout);

  {
    SquareMatrix<F,Z> A1;
    cout << "\nafter SquareMatrix()" << endl;
    cout << "A1:" << endl;
    A1.printOn(cout);

    {
      SquareMatrix<F,Z> A3(n,scalar);
      cout << "\nafter SquareMatrix(n,d)" << endl;
      cout << "A3:" << endl;
      A3.printOn(cout);

      A2=scalar;
      cout << "\nafter A2=scalar" << endl;
      cout << "A2:" << endl;
      A2.printOn(cout);

      A3.copy(A2);
      cout << "\nafter A3.copy(A2)" << endl;
      cout << "A3:" << endl;
      A3.printOn(cout);
    }
    cout << "\nafter delete A3" << endl;
    cout << "A2:" << endl;
    A2.printOn(cout);
  }

  SquareMatrix<F,Z> *T=0;
  {
    SquareMatrix<F,Z> A1;
    {
      SquareMatrix<F,Z> A3;
      A3.resize(n);
      cout << "\nafter A3.resize(n)" << endl;
      cout << "A3:" << endl;
      A3.printOn(cout);

      A1.resize(A3);
      A1=SquareMatrix<F,Z>::one_;
      cout << "\nafter A1.resize(A3)" << endl;
      cout << "A1:" << endl;
      A1.printOn(cout);
    }

    cout << "\nA1+A2:" << endl;
    T=A1+A2;
    T->printOn(cout);
    delete T; T=0;

    {
      Matrix<F,Z> M(n,n,scalar);
      cout << "\nA1+M:" << endl;
      T=A1+M;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nM+A1:" << endl;
      T=M+A1;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nA1-A2:" << endl;
      T=A1-A2;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nA1-M:" << endl;
      T=A1-M;
      T->printOn(cout);
      delete T; T=0;

      cout << "\nM-A1:" << endl;
      T=M-A1;
      T->printOn(cout);
      delete T; T=0;
    }

    cout << "\nA1*scalar:" << endl;
    T=A1*scalar;
    T->printOn(cout);
    delete T; T=0;

    cout << "\nA1/scalar:" << endl;
    T=A1/scalar;
    T->printOn(cout);
    delete T; T=0;
  }

  {
    SquareMatrix<F,Z> AA(2);
    AA(0,0)=scalar;        AA(0,1)=scalar+scalar+scalar;
    AA(1,0)=scalar+scalar; AA(1,1)=scalar+scalar+scalar+scalar;
    cout << "\nAA = " << endl;
    AA.printOn(cout);
//  cout << "\ntranspose(AA)" << endl;
//  T=AA.transpose();
//  T->printOn(cout);
//  delete T; T=0;

//  if (sizeof(F)!=sizeof(Z)) {
//    cout << "\nconjugateTranspose(A2)" << endl;
//    T=AA.conjugateTranspose();
//    T->printOn(cout);
//    delete T; T=0;
//  }

    cout << "\nreciprocalConditionNumber(I) = "
      << AA.reciprocalConditionNumber('I') << endl;
    cout << "reciprocalConditionNumber(O) = "
      << AA.reciprocalConditionNumber('O') << endl;

//  cout << "\ninverse:" << endl;
//  T=AA.inverse();
//  T->printOn(cout);
//  delete T; T=0;

    Vector<F,complex<F> > *lambda=0;
    {
      SquareMatrix<F,complex<F> > *U=
        OPERATOR_NEW SquareMatrix<F,complex<F> >(2);
      SquareMatrix<F,complex<F> > *V=
        OPERATOR_NEW SquareMatrix<F,complex<F> >(2);
      lambda=AA.eigenvalues(V,U);
      cout << "\neigenvalues:" << endl;
      lambda->printOn(cout);
      cout << "V = left eigenvectors = " << endl;
      V->printOn(cout);
      cout << "U = right eigenvectors = " << endl;
      U->printOn(cout);
      SquareMatrix<F,complex<F> > *R=
        OPERATOR_NEW SquareMatrix<F,complex<F> >(2);
      for (int j=0;j<2;j++) {
        for (int i=0;i<2;i++) {
          (*R)(i,j)=-(*U)(i,j)*(*lambda)[j];
          for (int k=0;k<2;k++) {
            (*R)(i,j)+=AA(i,k)*(*U)(k,j);
          }
        }
      }
      cout << "\nAA * U - U * lambda" << endl;
      R->printOn(cout);
      for (int j=0;j<2;j++) {
        for (int i=0;i<2;i++) {
          (*R)(i,j)=-(*lambda)[i]*conj((*V)(j,i));
          for (int k=0;k<2;k++) {
            (*R)(i,j)+=conj((*V)(k,i))*AA(k,j);
          }
        }
      }
      cout << "\nV^H * AA - lambda * V^H" << endl;
      R->printOn(cout);
      delete R; R=0;
      delete lambda; lambda=0;
      delete U; U=0;
      delete V; V=0;
    }

    {
      SquareMatrix<F,Z> BB(2);
      BB(0,0)=scalar;  BB(0,1)=scalar;
      BB(1,0)=-scalar; BB(1,1)=scalar;
      cout << "\nBB = " << endl;
      BB.printOn(cout);
      cout << "reciprocalConditionNumber('I') = "
        << BB.reciprocalConditionNumber('I') << endl;
      cout << "reciprocalConditionNumber('O') = "
        << BB.reciprocalConditionNumber('O') << endl;

//    T=BB.inverse();
//    cout << "\ninverse:" << endl;
//    T->printOn(cout);
//    delete T; T=0;

      {
        SquareMatrix<F,complex<F> > *V=
          OPERATOR_NEW SquareMatrix<F,complex<F> >(2);
        SquareMatrix<F,complex<F> > *U=
          OPERATOR_NEW SquareMatrix<F,complex<F> >(2);
        lambda=BB.eigenvalues(V,U);
        cout << "\neigenvalues:" << endl;
        lambda->printOn(cout);
        cout << "V = left eigenvectors = " << endl;
        V->printOn(cout);
        cout << "U = right eigenvectors = " << endl;
        U->printOn(cout);
        {
          SquareMatrix<F,complex<F> > R(2);
          for (int j=0;j<2;j++) {
            for (int i=0;i<2;i++) {
              R(i,j)=-(*U)(i,j)*(*lambda)[j];
              for (int k=0;k<2;k++) {
                R(i,j)+=BB(i,k)*(*U)(k,j);
              }
            }
          }
          cout << "\nBB * U - U * lambda" << endl;
          R.printOn(cout);
          for (int j=0;j<2;j++) {
            for (int i=0;i<2;i++) {
              R(i,j)=-(*lambda)[i]*conj((*V)(j,i));
              for (int k=0;k<2;k++) {
                R(i,j)+=conj((*V)(k,i))*BB(k,j);
              }
            }
          }
          cout << "\nV^H * BB - lambda * V^H" << endl;
          R.printOn(cout);
        }
        delete lambda; lambda=0;
        delete U; U=0;
        delete V; V=0;
      }

      T=AA * BB;
      cout << "SquareMatrix * SquareMatrix:" << endl;
      T->printOn(cout);
      delete T; T=0;
    }

    Matrix<F,Z> *M=0;
    {
      Matrix<F,Z> MM(n,1,scalar);
      cout << "\nSquareMatrix * Matrix:" << endl;
      M=A2 * MM;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      Matrix<F,Z> MM(n,n,scalar);
      cout << "\nSquareMatrix * Matrix: should be SquareMatrix" << endl;
      M=A2 * MM;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      Matrix<F,Z> MM(1,n,scalar);
      cout << "\nMatrix * SquareMatrix:" << endl;
      M=MM * A2;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      Matrix<F,Z> MM(n,n,scalar);
      cout << "\nMatrix * SquareMatrix: should be SquareMatrix" << endl;
      M=MM * A2;
      M->printOn(cout);
      delete M; M=0;
    }

    {
      Vector<F,Z> v(n);
      v=scalar;
      Vector<F,Z> *w=A2 * v;
      cout << "\nSquareMatrix * Vector:" << endl;
      w->printOn(cout);
      delete w; w=0;
    }

/*
    {
      Matrix<F,Z> B(2,1);
      B(0,0)=SquareMatrix<F,Z>::zero_;
      B(1,0)=6.*scalar;
      cout << "\nsquare system of linear equations" << endl;
      cout << "A:" << endl;
      AA.printOn(cout);
      cout << "B:" << endl;
      B.printOn(cout);
      cout << "solve(A,B,'L','N'):" << endl;
      {
        Matrix<F,Z> X(2,1);
        AA.solve(B,X);
        X.printOn(cout);
      }

      cout << "\nsolve(A,B,'L','T'):" << endl;
      Matrix<F,Z> X(2,1);
      AA.solve(B,X,'L','T');
      X.printOn(cout);
    }

    {
      Matrix<F,Z> B(1,2);
      B(0,0)=SquareMatrix<F,Z>::zero_;
      B(0,1)=6.*scalar;
      cout << "\nsolve(A,B,'R','N'):" << endl;
      {
        Matrix<F,Z> X(1,2);
        AA.solve(B,X,'R','N');
        X.printOn(cout);
      }

      cout << "\nsolve(A,B,'R','N'):" << endl;
      Matrix<F,Z> X(1,2);
      AA.solve(B,X,'R','T');
      X.printOn(cout);
    }
*/
  }
}
