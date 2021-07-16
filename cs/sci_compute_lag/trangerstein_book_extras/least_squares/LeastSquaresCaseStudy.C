#include <fstream>
#include <iostream>
#include <limits>
#include <stdlib.h>
#include "CompleteOrthogonalDecomposition.H"
#include "GramSchmidtQRFactorization.H"
#include "HouseholderQRFactorization.H"
#include "MemoryDebugger.H"
#include "SingularValueDecomposition.H"
#include "TimedObject.H"
#include "Tracer.H"

using namespace std;

int main(int /*argc*/,char** /*argv*/) {
  MemoryDebugger md(1);
  {
/*
    int m=3,n=2;
    Matrix<double,double> AA(m,n,0.);
    AA(0,0)=1.;
    AA(0,1)=2.;
    double eps=sqrt(numeric_limits<double>::epsilon()*0.5);
    AA(1,0)=eps;
    AA(2,1)=2.*eps;
*/
//
    ifstream in_file;
    in_file.open("illc1033.mtx",ios::in);
    int m,n,nonzero;
    in_file >> m >> n >> nonzero;

    Matrix<double,double> AA(m,n,0.);
    for (int k=0;k<nonzero;k++) {
      int i,j;
      double aij;
      in_file >> i >> j >> aij;
      i--; j--; // Fortran array indices converted to C indices
      AA(i,j)=aij;
    }
    in_file.close();
//

    Matrix<float,float> A(m,n);
    for (int j=0;j<n;j++) {
      for (int i=0;i<m;i++) A(i,j)=AA(i,j); // lower precision
    }
#ifdef DEBUG
//  cout << "\tA = " << endl;
//  A.printOn(cout);
#endif

    Vector<double,double> xx_exact(n);
    for (int j=0;j<n;j++) xx_exact[j]=drand48();
//  xx_exact[0]=0.5;
//  xx_exact[1]=1.;
    Vector<double,double> *bb=AA*xx_exact; // matrix-vector multiply

    Vector<float,float> x_exact(n);
    for (int j=0;j<n;j++) x_exact[j]=xx_exact[j]; // lower precision
    Vector<float,float> *b=OPERATOR_NEW Vector<float,float>(m);
    for (int i=0;i<m;i++) (*b)[i]=(*bb)[i]; // lower precision
    float bnorm = b->nrm2();
    double bbnorm = bb->nrm2();
#ifdef DEBUG
//  cout << "\tx_exact = " << endl;
//  x_exact.printOn(cout);
//  cout << "\tb = " << endl;
//  b->printOn(cout);
#endif

    Vector<float,float> x_approx(n);
    Vector<double,double> xx_approx(n);
//  copy A, call sgesv to factor A copy and compute x_approx, forget factors:
//  A.solve(b,x_approx);

    double AA_norm,AA_inv_norm;
    double AA_posteriori_norm,AA_posteriori_inv_norm;
    double AA_condition,AA_posteriori_condition;
    { TRACER_CALL(t,"SingularValueDecomposition");
      TimedObject z("SingularValueDecomposition");
      { Timer tz(&z);
      SingularValueDecomposition<double,double> svd(AA);
      Vector<double,double> *sv=svd.singularValues();
      AA_norm=(*sv)[0];
      AA_inv_norm=1. / (*sv)[n-1];
      AA_posteriori_norm=0.5*(1.+sqrt(1.+4.*pow(AA_norm,2)));
      AA_posteriori_inv_norm=max(1.,
        0.5*AA_inv_norm*(AA_inv_norm+sqrt(4.+pow(AA_inv_norm,2))));
      AA_condition=AA_norm*AA_inv_norm;
      AA_posteriori_condition=AA_posteriori_norm * AA_posteriori_inv_norm;
      cout << "\t|| A ||, || A^+ || = " << AA_norm << " "
        << AA_inv_norm << endl;
      cout << "\ta posteriori error estimates: norm, norm of inverse = "
        << AA_posteriori_norm << " " << AA_posteriori_inv_norm << endl;
      cout << "\tcondition number = " << AA_condition << endl;
      cout << "\ta posteriori error condition number = "
        << AA_posteriori_condition << endl;
/*
      svd.solve(*b,x_approx,rcond,Factorization::NO_TRANSPOSE);

      Vector<float,float> r(m);
      r.copy(*b);
      A.gemv(1.,x_approx,-1.,r); // r=A*x-b
      float rnorm = r.nrm2();
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      float relative_error_in_x=rnorm/(bnorm*rcond);
      cout << "\testimated relative error in x = " << relative_error_in_x
           << endl;
      Vector<float,float> *dx=x_approx-x_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / x_exact.nrm2() << endl;
      delete dx; dx=0;
*/
      }
      z.printOn(cout);
    }

    { TRACER_CALL(t,"HouseholderQRFactorization(double)");
      TimedObject z("HouseholderQRFactorization(double)");
      { Timer tz(&z);
      HouseholderQRFactorization<double,double> hf(AA);
        //copy AA and factor with dgeqrf
      hf.solve(*bb,xx_approx,Factorization::NO_TRANSPOSE);
        // call dgeqrs to compute x_approx from factors

      Vector<double,double> rr(m);
      rr.copy(*bb);
      AA.gemv(-1.,xx_approx,1.,rr); // r=b-A*x
      double rnorm = rr.nrm2();
      cout << "\trnorm, bnorm = " << rnorm << " " << bbnorm << endl;
      Vector<double,double> *dx=xx_approx-xx_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / xx_exact.nrm2() << endl;
      delete dx; dx=0;

      Vector<double,double> Atranspose_r(n,0.);
      AA.gemv(1.,rr,0.,Atranspose_r,'C');
      cout << "\ta priori estimate for relative error = "
        << AA_posteriori_condition * Atranspose_r.nrm2() / bbnorm
        << endl; 
      }
      z.printOn(cout);
    }
    { TRACER_CALL(t,"HouseholderQRFactorization(float)");
      TimedObject z("HouseholderQRFactorization(float)");
      { Timer tz(&z);
      HouseholderQRFactorization<float,float> hf(A);
        //copy A and factor with sgeqrf
      hf.solve(*b,x_approx,Factorization::NO_TRANSPOSE);
        // call sgeqrs to compute x_approx from factors

      Vector<float,float> r(m);
      r.copy(*b);
      A.gemv(-1.,x_approx,1.,r); // r=b-A*x
      float rnorm = r.nrm2();
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      Vector<float,float> *dx=x_approx-x_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / x_exact.nrm2() << endl;
      delete dx; dx=0;

      Vector<float,float> Atranspose_r(n,0.);
      A.gemv(1.,r,0.,Atranspose_r,'C');
      cout << "\ta priori estimate for relative error = "
        << AA_posteriori_condition * Atranspose_r.nrm2() / bnorm
        << endl; 

/*
      float berr,ferr;
      hf.improve(*b,x_approx,Factorization::NO_TRANSPOSE);
      r.copy(*b);
      A.gemv(-1.,x_approx,1.,r); // r=b-A*x
      rnorm = r.nrm2();
      cout << "\tafter improve, rnorm, bnorm = " << rnorm << " " << bnorm
           << endl;
      dx=x_approx-x_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / x_exact.nrm2() << endl;
      delete dx; dx=0;

      Vector<double,double> xx_approx(n);
      for (int i=0;i<n;i++) xx_approx[i]=x_approx[i]; // higher precision
      Vector<double,double> rr(m);
      rr.copy(*bb);
      AA.gemv(-1.,xx_approx,1.,rr);
      for (int i=0;i<m;i++) r[i]=rr[i]; // lower precision
      Vector<float,float> d(n);
      hf.solve(r,d);
      x_approx+=d;
      for (int i=0;i<n;i++) xx_approx[i]=x_approx[i];
      rr.copy(*bb);
      AA.gemv(-1.,xx_approx,1.,rr);
      rnorm = rr.nrm2();
      cout << "\n\tafter improve in double, rnorm, bnorm = "
           << rnorm << " " << bnorm << endl;
      dx=x_approx-x_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / x_exact.nrm2() << endl;
      delete dx; dx=0;
*/

      }
      z.printOn(cout);
    }
    { TRACER_CALL(t,"HouseholderQRFactorization with pivoting");
      TimedObject z("HouseholderQRFactorization with pivoting");
      { Timer tz(&z);
      HouseholderQRFactorization<float,float>
        hf(A,Factorization::PIVOT_COLUMNS);
        //copy A and factor with sgeqrf
      hf.solve(*b,x_approx,Factorization::NO_TRANSPOSE);
        // call sgeqrs to compute x_approx from factors

      Vector<float,float> r(m);
      r.copy(*b);
      A.gemv(-1.,x_approx,1.,r); // r=b-A*x
      float rnorm = r.nrm2();
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      Vector<float,float> *dx=x_approx-x_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / x_exact.nrm2() << endl;
      delete dx; dx=0;
      }
      z.printOn(cout);
    }
/*
    { TRACER_CALL(t,"CompleteOrthogonalDecomposition");
      TimedObject z("CompleteOrthogonalDecomposition");
      { Timer tz(&z);
      float rcond=numeric_limits<float>::epsilon();
      CompleteOrthogonalDecomposition<float,float> cod(A,rcond);
#ifdef DEBUG
      cout << "\tcod = " << endl;
      cod.printOn(cout);
#endif
      cout << "\trank = " << cod.rank() << endl;
      float rnorm=cod.solve(*b,x_approx,Factorization::NO_TRANSPOSE);
      cout << "\trnorm, bnorm = " << rnorm << " " << bnorm << endl;
      Vector<float,float> *dx=x_approx-x_exact;
#ifdef DEBUG
      cout << "\tx_approx = " << endl;
      x_approx.printOn(cout);
      cout << "\tdx = " << endl;
      dx->printOn(cout);
#endif
      cout << "\trelative error in x = "
           << dx->nrm2() / x_exact.nrm2() << endl;
      delete dx; dx=0;
      }
      z.printOn(cout);
    }
*/
    { TRACER_CALL(t,"GramSchmidtQRFactorization(double)");
      TimedObject z("GramSchmidtQRFactorization(double)");
      { Timer tz(&z);
      GramSchmidtQRFactorization<double,double> gs(AA);
      Vector<double,double> rr(m);
      gs.solveOverdetermined(*bb,xx_approx,rr);
      cout << "\trnorm, bnorm = " << rr.nrm2() << " " << bbnorm << endl;
      Vector<double,double> *dx=xx_approx-xx_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / xx_exact.nrm2() << endl;
      delete dx; dx=0;

      Vector<double,double> Axmb(m);
      Axmb.copy(*bb);
      AA.gemv(-1.,xx_approx,1.,Axmb); // b-A*x
      Vector<double,double> *dr=rr-Axmb;
      double drnorm = dr->nrm2();
      delete dr; dr=0;

      Vector<double,double> Atranspose_r(n,0.);
      AA.gemv(1.,rr,0.,Atranspose_r,'C');
      cout << "\ta priori estimate for relative error = "
        << AA_posteriori_condition * sqrt( drnorm*drnorm
          + pow(Atranspose_r.nrm2(),2) ) / bbnorm
        << endl; 
      }
      z.printOn(cout);
    }
    { TRACER_CALL(t,"GramSchmidtQRFactorization(float)");
      TimedObject z("GramSchmidtQRFactorization(float)");
      { Timer tz(&z);
      GramSchmidtQRFactorization<float,float> gs(A);
      Vector<float,float> r(m);
      gs.solveOverdetermined(*b,x_approx,r);
      cout << "\trnorm, bnorm = " << r.nrm2() << " " << bnorm << endl;
      Vector<float,float> *dx=x_approx-x_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / x_exact.nrm2() << endl;
      delete dx; dx=0;
/*
      gs.improveOverdetermined(*b,x_approx,r);
      cout << "\tafter improve, rnorm, bnorm = " << r.nrm2() << " "
           << bnorm << endl;
      dx=x_approx-x_exact;
      cout << "\trelative error in x = "
           << dx->nrm2() / x_exact.nrm2() << endl;
      delete dx; dx=0;
*/
      Vector<float,float> Axmb(m);
      Axmb.copy(*b);
      A.gemv(-1.,x_approx,1.,Axmb); // b-A*x
      Vector<float,float> *dr=r-Axmb;
      float drnorm = dr->nrm2();
      delete dr; dr=0;

      Vector<float,float> Atranspose_r(n,0.);
      A.gemv(1.,r,0.,Atranspose_r,'C');
      cout << "\ta priori estimate for relative error = "
        << AA_posteriori_condition * sqrt( drnorm*drnorm
          + pow(Atranspose_r.nrm2(),2) ) / bbnorm
        << endl; 
      }
      z.printOn(cout);
    }
    delete b; b=0;
    delete bb; bb=0;
  }
}
