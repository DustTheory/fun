#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "GaussianFactorization.H"
#include "MemoryDebugger.H"
#include "SquareMatrix.H"
#include "TimedObject.H"
#include "Tracer.H"

using namespace std;

int main(int /*argc*/,char** /*argv*/) {
  MemoryDebugger md(1);
  {
    ifstream in_file;
    in_file.open("s1rmq4m1.mtx",ios::in);
    int m,n,nonzero;
    in_file >> m >> n >> nonzero;

    SquareMatrix<double,double> AA(m,0.);
    for (int k=0;k<nonzero;k++) {
      int i,j;
      double aij;
      in_file >> i >> j >> aij;
      i--; j--; // Fortran array indices converted to C indices
      AA(i,j)=aij;
      if (i != j) {
        AA(j,i)=aij;
      }
    }
    in_file.close();

    SquareMatrix<float,float> A(m,0.);
    for (int j=0;j<m;j++) {
      for (int i=0;i<m;i++) A(i,j)=AA(i,j);
    }

    Vector<double,double> xx_exact(m);
    for (int j=0;j<m;j++) xx_exact[j]=rand();
    Vector<double,double> *bb=AA*xx_exact; // matrix-vector multiply

    Vector<float,float> x_exact(m);
    for (int j=0;j<m;j++) x_exact[j]=xx_exact[j];
    Vector<float,float> *b=OPERATOR_NEW Vector<float,float>(m);
    for (int i=0;i<m;i++) (*b)[i]=(*bb)[i];
    float bnorm = abs((*b)[b->amax()]);

    Vector<float,float> x_gauss(m);
//  copy A, call sgesv to factor A copy and compute x_gauss, forget factors:
//  A.solve(b,x_gauss);

    { TRACER_CALL(t,"GaussianFactorization with row pivoting");
      TimedObject z("GaussianFactorization with row pivoting");
      { Timer tz(&z);
      GaussianFactorization<float,float,SquareMatrix<float,float> >
        gf(A,Factorization::PIVOT_ROWS); //copy A and factor with sgetrf
      gf.solve(*b,x_gauss); // call sgetrs to compute x_gauss from factors

      Vector<float,float> r(m);
      r.copy(*b);
      A.gemv(-1.,x_gauss,1.,r); // r=b-A*x
      float rnorm = abs(r[r.amax()]);
      float rcond=
        gf.reciprocalConditionNumber(Factorization::INFINITY_NORM);
      cout << "\testimated condition number = " << 1./rcond << endl;
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      float relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      Vector<float,float> *dx=x_gauss-x_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(x_exact[x_exact.amax()])
           << endl;
      delete dx; dx=0;

      float berr,ferr;
      gf.improve(*b,x_gauss,berr,ferr);
      cout << "\n\tafter improve in float, berr,ferr = "
           << berr << " " << ferr << endl;
/*
      r.copy(*b);
      A.gemv(-1.,x_gauss,1.,r); // r=b-A*x
      rnorm = abs(r[r.amax()]);
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      dx=x_gauss-x_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(x_exact[x_exact.amax()])
           << endl;
      delete dx; dx=0;

      Vector<double,double> xx_gauss(m);
      for (int i=0;i<m;i++) xx_gauss[i]=x_gauss[i];
      Vector<double,double> rr(m);
      rr.copy(*bb);
      AA.gemv(-1.,xx_gauss,1.,rr);
      for (int i=0;i<m;i++) r[i]=rr[i];
      Vector<float,float> d(m);
      gf.solve(r,d);
      x_gauss+=d;
//    gf.improve(*b,x_gauss,berr,ferr);
//    cout << "\n\tafter improve in double, berr,ferr = "
//         << berr << " " << ferr << endl;
      for (int i=0;i<m;i++) xx_gauss[i]=x_gauss[i];
      rr.copy(*bb);
      AA.gemv(-1.,xx_gauss,1.,rr);
      rnorm = abs(rr[rr.amax()]);
      cout << "\n\tafter improve in double, rnorm, bnorm = "
           << rnorm << " " << bnorm << endl;
      relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      dx=x_gauss-x_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(x_exact[x_exact.amax()])
           << endl;
      delete dx; dx=0;
*/

/*
      Vector<float,float> bone(m,1.);
      Vector<float,float> xone_exact(m);
      gf.solve(bone,xone_exact);
      Vector<float,float> *A_xone=A*xone_exact;
      Vector<float,float> xone(m);
      gf.solve(*A_xone,xone);

      r.copy(*A_xone);
      A.gemv(1.,xone,-1.,r);
      rnorm=abs(r[r.amax()]);
      cout << "\n\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      float relative_error_in_xone=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in xone = "
           << relative_error_in_xone << endl;
      Vector<float,float> *dx2=xone_exact-xone;
      cout << "\ttrue relative error in xone = "
           << abs((*dx2)[dx2->amax()]) / abs(xone_exact[xone_exact.amax()])
           << endl;
      delete dx2; dx2=0;

      gf.improve(*A_xone,xone,berr,ferr);
      cout << "\tafter improve, berr,ferr = " << berr << " " << ferr
           << endl;
      r.copy(*A_xone);
      A.gemv(1.,xone,-1.,r); // r=A*x-b
      rnorm = abs(r[r.amax()]);
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      dx=x_gauss-x_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(x_exact[x_exact.amax()])
           << endl;
      delete dx; dx=0;
      delete A_xone; A_xone=0;
*/
      }
      z.printOn(cout);
    }
/*
    { TRACER_CALL(t,"GaussianFactorization with row and column pivoting");
      TimedObject z("GaussianFactorization with row and column pivoting");
      { Timer tz(&z);
      GaussianFactorization<float,float,SquareMatrix<float,float> >
        gf(A,Factorization::PIVOT_ROWS_AND_COLUMNS);
      gf.solve(*b,x_gauss);
      float rcond=
        gf.reciprocalConditionNumber(Factorization::INFINITY_NORM);
      cout << "\testimated condition number = " << 1./rcond << endl;
      Vector<float,float> r(m);
      r.copy(*b);
      A.gemv(1.,x_gauss,-1.,r);
      float rnorm = abs(r[r.amax()]);
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      float relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      Vector<float,float> *dx=x_gauss-x_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(x_exact[x_exact.amax()])
           << endl;
      delete dx; dx=0;

      Vector<float,float> bone(m,1.);
      Vector<float,float> xone_exact(m);
      gf.solve(bone,xone_exact);
      Vector<float,float> *A_xone=A*xone_exact;
      Vector<float,float> xone(m);
      gf.solve(*A_xone,xone);
      r.copy(*A_xone);
      A.gemv(1.,xone,-1.,r);
      rnorm=abs(r[r.amax()]);
      bnorm=abs((*A_xone)[A_xone->amax()]);
      cout << "\n\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      float relative_error_in_xone=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in xone = "
           << relative_error_in_xone << endl;
      Vector<float,float> *dx2=xone_exact-xone;
      cout << "\ttrue relative error in xone = "
           << abs((*dx2)[dx2->amax()]) / abs(xone_exact[xone_exact.amax()])
           << endl;
      delete dx2; dx2=0;
      delete A_xone; A_xone=0;
      }
      z.printOn(cout);
    }
*/
/*
    { TRACER_CALL(t,"GaussianFactorization with row pivoting and equilibration");
      TimedObject z("GaussianFactorization with row pivoting and equilibration");
      { Timer tz(&z);
      GaussianFactorization<float,float,SquareMatrix<float,float> >
        gf2(A,Factorization::PIVOT_ROWS,
        Factorization::EQUILIBRATE_ROWS_AND_COLUMNS);
      gf2.solve(*b,x_gauss);

      Vector<float,float> r(m);
      r.copy(*b);
      A.gemv(1.,x_gauss,-1.,r); // r=A*x-b
      float rnorm = abs(r[r.amax()]);
      float rcond=
        gf2.reciprocalConditionNumber(Factorization::INFINITY_NORM);
      cout << "\testimated condition number = " << 1./rcond << endl;
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      float relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      Vector<float,float> *dx=x_gauss-x_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(x_exact[x_exact.amax()])
           << endl;
      delete dx; dx=0;

      float berr,ferr;
      gf2.improve(*b,x_gauss,berr,ferr);
      cout << "\tafter improve, berr,ferr = " << berr << " " << ferr
           << endl;
      r.copy(*b);
      A.gemv(1.,x_gauss,-1.,r);
      rnorm = abs(r[r.amax()]);
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      dx=x_gauss-x_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(x_exact[x_exact.amax()])
           << endl;
      delete dx; dx=0;

      Vector<float,float> bone(m,1.);
      Vector<float,float> xone_exact(m);
      gf2.solve(bone,xone_exact);
      Vector<float,float> *A_xone=A*xone_exact;
      Vector<float,float> xone(m);
      gf2.solve(*A_xone,xone);
      gf2.improve(*A_xone,xone,berr,ferr);
      cout << "\n\tafter improve, berr,ferr = " << berr << " " << ferr
           << endl;
      r.copy(*A_xone);
      A.gemv(1.,xone,-1.,r);
      rnorm=abs(r[r.amax()]);
      bnorm=abs((*A_xone)[A_xone->amax()]);
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      float relative_error_in_xone=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in xone = "
           << relative_error_in_xone << endl;
      Vector<float,float> *dx2=xone_exact-xone;
      cout << "\ttrue relative error in xone = "
           << abs((*dx2)[dx2->amax()]) / abs(xone_exact[xone_exact.amax()])
           << endl;
      delete dx2; dx2=0;
      delete A_xone; A_xone=0;
      }
      z.printOn(cout);
    }
*/
    delete b; b=0;

    Vector<double,double> xx_gauss(m);
//copy AA, call dgesv to factor AA copy, compute xx_gauss, forget factors:
//  AA.solve(bb,xx_gauss);

    { TRACER_CALL(t,"GaussianFactorization with row pivoting");
      TimedObject z("GaussianFactorization with row pivoting");
      { Timer tz(&z);
      GaussianFactorization<double,double,SquareMatrix<double,double> >
        gf(AA,Factorization::PIVOT_ROWS);
      gf.solve(*bb,xx_gauss);

      double bnorm = abs((*bb)[bb->amax()]);
      double rcond=
        gf.reciprocalConditionNumber(Factorization::INFINITY_NORM);
      cout << "\testimated condition number = " << 1./rcond << endl;
      Vector<double,double> r(m);
      r.copy(*bb);
      AA.gemv(1.,xx_gauss,-1.,r); // r=AA*x-bb
      double rnorm = abs(r[r.amax()]);
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      double relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      Vector<double,double> *dx=xx_gauss-xx_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(xx_exact[xx_exact.amax()])
           << endl;
      delete dx; dx=0;


      double berr,ferr;
      gf.improve(*bb,xx_gauss,berr,ferr);
      cout << "\tafter improve, berr,ferr = " << berr << " " << ferr
           << endl;
/*
      r.copy(*bb);
      AA.gemv(1.,xx_gauss,-1.,r); // r=AA*x-bb
      rnorm = abs(r[r.amax()]);
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      dx=xx_gauss-xx_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(xx_exact[xx_exact.amax()])
           << endl;
      delete dx; dx=0;
*/

/*
      Vector<double,double> bone(m,1.);
      Vector<double,double> xone_exact(m);
      gf.solve(bone,xone_exact);
      Vector<double,double> *A_xone=AA*xone_exact;
      Vector<double,double> xone(m);
      gf.solve(*A_xone,xone);
      gf.improve(*bb,xx_gauss,berr,ferr);
      cout << "\tafter improve, berr,ferr = " << berr << " " << ferr
           << endl;
      r.copy(*A_xone);
      AA.gemv(1.,xone,-1.,r);
      rnorm=abs(r[r.amax()]);
      bnorm=abs((*A_xone)[A_xone->amax()]);
      cout << "\n\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      double relative_error_in_xone=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in xone = "
           << relative_error_in_xone << endl;
      Vector<double,double> *dx2=xone_exact-xone;
      cout << "\ttrue relative error in xone = "
           << abs((*dx2)[dx2->amax()]) / abs(xone_exact[xone_exact.amax()])
           << endl;
      delete dx2; dx2=0;
      delete A_xone; A_xone=0;
*/
      }
      z.printOn(cout);
    }
    delete bb; bb=0;
  }
}
