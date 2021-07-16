#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "GaussianFactorization.H"
#include "MemoryDebugger.H"
#include "SpecializedMatrix.H"
#include "TimedObject.H"
#include "Tracer.H"

using namespace std;

int main(int /*argc*/,char** /*argv*/) {
  MemoryDebugger md(1);
  {
    int n=4;
    Vector<double,double> v(n);
    for (int i=0;i<n;i++) {
      v[i]=0.1*static_cast<double>(i)/static_cast<double>(n-1);
    }
    VandermondeMatrix<double,double> AA(v);

    SquareMatrix<float,float> A(n);
    for (int j=0;j<n;j++) {
      for (int i=0;i<n;i++) A(i,j)=AA(i,j);
    }

    Vector<double,double> xx_exact(n);
    for (int j=0;j<n;j++) xx_exact[j]=rand();
    Vector<double,double> *bb=AA*xx_exact; // matrix-vector multiply

    Vector<float,float> x_exact(n);
    for (int j=0;j<n;j++) x_exact[j]=xx_exact[j];
    Vector<float,float> *b=OPERATOR_NEW Vector<float,float>(n);
    for (int i=0;i<n;i++) (*b)[i]=(*bb)[i];
    float bnorm = abs((*b)[b->amax()]);

    Vector<float,float> x_gauss(n);

    { TRACER_CALL(t,"GaussianFactorization");
      GaussianFactorization<float,float,SquareMatrix<float,float> >
        gf(A,Factorization::PIVOT_ROWS); //copy A and factor with sgetrf
      gf.solve(*b,x_gauss); // call sgetrs to compute x_gauss from factors

      float rcond=
        gf.reciprocalConditionNumber(Factorization::INFINITY_NORM);
      cout << "\testimated condition number = " << 1./rcond << endl;

      Vector<float,float> r(n);
      r.copy(*b);
      A.gemv(-1.,x_gauss,1.,r); // r=b-A*x
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

      float berr,ferr;
      gf.improve(*b,x_gauss,berr,ferr);
      cout << "\n\tberr,ferr = " << berr << " " << ferr << endl;
      r.copy(*b);
      A.gemv(-1.,x_gauss,1.,r); // r=b-A*x
      rnorm = abs(r[r.amax()]);
      cout << "\tafter improve in float, rnorm, bnorm = "
           << rnorm << " " << bnorm << endl;
      relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      dx=x_gauss-x_exact;
      cout << "\ttrue relative error in x = "
           << abs((*dx)[dx->amax()]) / abs(x_exact[x_exact.amax()])
           << endl;
      gf.improve(*b,x_gauss,berr,ferr);
      cout << "tberr,ferr = " << berr << " " << ferr << endl;
      delete dx; dx=0;
    }

    { TRACER_CALL(t,"GaussianFactorization");
      GaussianFactorization<float,float,SquareMatrix<float,float> >
        gf(A,Factorization::PIVOT_ROWS); //copy A and factor with sgetrf
      gf.solve(*b,x_gauss); // call sgetrs to compute x_gauss from factors

      float rcond=
        gf.reciprocalConditionNumber(Factorization::INFINITY_NORM);

      Vector<double,double> xx_gauss(n);
      for (int i=0;i<n;i++) xx_gauss[i]=x_gauss[i];
      Vector<double,double> rr(n);
      rr.copy(*bb);
      AA.gemv(-1.,xx_gauss,1.,rr);

      Vector<float,float> r(n);
      for (int i=0;i<n;i++) r[i]=rr[i];
      Vector<float,float> d(n);
      gf.solve(r,d);
      x_gauss+=d;

      float berr,ferr;
      gf.improve(*b,x_gauss,berr,ferr);
      cout << "\tafter improve in double, berr,ferr = "
           << berr << " " << ferr << endl;

      r.copy(*b);
      A.gemv(-1.,x_gauss,1.,r); // r=b-A*x
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
    }

    { TRACER_CALL(t,"GaussianFactorization equilibration");
      GaussianFactorization<float,float,SquareMatrix<float,float> >
        gf2(A,Factorization::PIVOT_ROWS,
        Factorization::EQUILIBRATE_ROWS_AND_COLUMNS);
      gf2.solve(*b,x_gauss);

      float rcond=
        gf2.reciprocalConditionNumber(Factorization::INFINITY_NORM);
      cout << "\testimated condition number = " << 1./rcond << endl;

      float berr,ferr;
      gf2.improve(*b,x_gauss,berr,ferr);
      cout << "\tberr,ferr = " << berr << " " << ferr << endl;

      Vector<float,float> r(n);
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
      gf2.improve(*b,x_gauss,berr,ferr);
      cout << "\tafter improve, berr,ferr = "
           << berr << " " << ferr << endl;
      delete dx; dx=0;
    }
    delete b; b=0;

    Vector<double,double> xx_gauss(n);
//copy AA, call dgesv to factor AA copy, compute xx_gauss, forget factors:
//  AA.solve(bb,xx_gauss);

    { TRACER_CALL(t,"GaussianFactorization");
      GaussianFactorization<double,double,SquareMatrix<double,double> >
        gf(AA,Factorization::PIVOT_ROWS);
      gf.solve(*bb,xx_gauss);

      double bnorm = abs((*bb)[bb->amax()]);
      double rcond=
        gf.reciprocalConditionNumber(Factorization::INFINITY_NORM);
      cout << "\testimated condition number = " << 1./rcond << endl;
      Vector<double,double> r(n);
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
    }
    delete bb; bb=0;

    { TRACER_CALL(t,"GaussianFactorization");
      int m=8192;
      SquareMatrix<double,double> S(m);
      for (int j=0;j<m;j++) {
        for (int i=0;i<m;i++) S(i,j)=rand();
      }
      Vector<double,double> y(m);
      for (int i=0;i<m;i++) y[i]=rand();
      Vector<double,double> *rhs=S*y;
      GaussianFactorization<double,double,SquareMatrix<double,double> >
        gf(S,Factorization::PIVOT_ROWS);
      Vector<double,double> y_gauss(m);
      gf.solve(*rhs,y_gauss);

      double rcond=
        gf.reciprocalConditionNumber(Factorization::INFINITY_NORM);
      cout << "\testimated condition number = " << 1./rcond << endl;

      double berr,ferr;
      gf.improve(*rhs,y_gauss,berr,ferr);
      cout << "\tberr,ferr = " << berr << " " << ferr << endl;
      gf.improve(*rhs,y_gauss,berr,ferr);
      cout << "\tafter improve, berr,ferr = "
           << berr << " " << ferr << endl;
    }
  }
}
