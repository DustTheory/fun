#include <cmath>
#include <stdlib.h>
#include "NonlinearSystemSolver.H"

Vector<double,double>* Rosenbrock(const Vector<double,double> &x) {
  CHECK_SAME(2,x.size());
  Vector<double,double> *f=OPERATOR_NEW Vector<double,double>(2);
  (*f)[0]=10.*(x[1]-x[0]*x[0]);
  (*f)[1]=1.-x[0];
  return f;
}

Vector<double,double>* FreudensteinRoth(const Vector<double,double> &x) {
  CHECK_SAME(2,x.size());
  Vector<double,double> *f=OPERATOR_NEW Vector<double,double>(2);
  (*f)[0]=-13.+x[0]+((5.-x[1])*x[1]-2.)*x[1];
  (*f)[1]=-29.+x[0]+((x[1]+1.)*x[1]-14.)*x[1];
  return f;
}

Vector<double,double>* PowellBadlyScaled(const Vector<double,double> &x) {
  CHECK_SAME(2,x.size());
  Vector<double,double> *f=OPERATOR_NEW Vector<double,double>(2);
  (*f)[0]=10000.*x[0]*x[1]-1.;
  (*f)[1]=exp(-x[0])+exp(-x[1])-1.0001;
  return f;
}

double helicalValleyTheta(double a,double b) {
  double val=atan(b/a)/(2.*M_PI);
  if (a<0.) val+=.5;
  return val;
}
Vector<double,double>* HelicalValley(const Vector<double,double> &x) {
  CHECK_SAME(3,x.size());
  Vector<double,double> *f=OPERATOR_NEW Vector<double,double>(3);
  (*f)[0]=10*(x[2]-10.*helicalValleyTheta(x[0],x[1]));
  (*f)[1]=10*(sqrt(x[0]*x[0]+x[1]*x[1])-1.);
  (*f)[2]=x[2];
  return f;
}

Vector<double,double>* PowellSingular(const Vector<double,double> &x) {
  CHECK_SAME(4,x.size());
  Vector<double,double> *f=OPERATOR_NEW Vector<double,double>(4);
  (*f)[0]=x[0]+10.*x[1];
  (*f)[1]=sqrt(5.)*(x[2]-x[3]);
  (*f)[2]=pow(x[1]-2.*x[2],2);
  (*f)[3]=sqrt(10.)*pow(x[0]-x[3],2);
  return f;
}

int main( int /* argc */, char* /* (argv[] */) {
  MinpackSolver ms;
  double tol=1.e-12;

  Vector<double,double> x2(2);
  x2[0]=-1.2;
  x2[1]=1.;
  ms.solve(tol,x2,Rosenbrock);
  cout << "Rosenbrock root = " << x2[0] << " " << x2[1] << endl;

//minpack gets close to a local min:
  x2[0]=.5;
  x2[1]=-2.;
  ms.solve(tol,x2,FreudensteinRoth);
  cout << "FreudensteinRoth root = " << x2[0] << " " << x2[1] << endl;

  x2[0]=0.;
  x2[1]=1.;
  ms.solve(tol,x2,PowellBadlyScaled);
  cout << "PowellBadlyScaled root = " << x2[0] << " " << x2[1] << endl;

  Vector<double,double> x3(3);
  x3[0]=-1;
  x3[1]=0.;
  x3[2]=0.;
  ms.solve(tol,x3,HelicalValley);
  cout << "HelicalValley root = " << x3[0] << " " << x3[1] << " " << x3[2]
       << endl;

  Vector<double,double> x4(4);
  x4[0]=3;
  x4[1]=-1.;
  x4[2]=0.;
  x4[3]=1.;
  ms.solve(tol,x4,PowellSingular);
  cout << "PowellSingular root = " << x4[0] << " " << x4[1] << " " << x4[2]
       << " " << x4[3] << endl;

  return EXIT_SUCCESS;
}
