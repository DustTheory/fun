#include <float.h>
#include <string.h>
#include "MemoryDebugger.H"
#include "ScalingFunction.H"
#include "SquareMatrix.H"
#include "Tracer.H"
//#include "Vector.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ScalingFunction::ScalingFunction(const Signal &sig) : lowpass(&sig) {
//TRACER_CALL(t,"ScalingFunction::ScalingFunction");
  int lp_lo=lowpass->firstIndex();
  int lp_hi=lowpass->lastIndex();
  int len=lowpass->length();
#ifdef DEBUG
//for (int n=lowpass->firstIndex();n<=lowpass->lastIndex();n++) {
//  cout << "\ts[ " << n << " ] = " << lowpass->value(n) << endl;
//}

  double summ=0.;
  double sumr=0.;
  for (int n=lowpass->firstIndex();n<=lowpass->lastIndex();n++) {
    double svn=lowpass->value(n);
    summ+=svn;
    if (n%2!=0) svn=-svn;
    sumr+= svn;
  }
  CHECK_TEST(abs(summ-1.)<512.*DBL_EPSILON);
  CHECK_TEST(abs(sumr)<512.*DBL_EPSILON);
#endif
  SquareMatrix<double,double> M(len-1,0.);
  Vector<double,double> b(len-1,1.);
  double eigvalue=1.+DBL_EPSILON;
  for (int i=0;i<len-1;i++) {
    int jlo=max(2*i-len+1,0);
    int jhi=min(2*i,len-2);
    for (int j=jlo;j<=jhi;j++) M(i,j)=2.*lowpass->value(2*i-j+lp_lo);
    M(i,i)-=eigvalue;
  }
#ifdef DEBUG
//for (int i=0;i<len-1;i++) {
//  cout << "\tM[ " << i << " , * ] = ";
//  for (int j=0;j<len-1;j++) {
//    cout << M(i,j) << " ";
//  }
//  cout << endl;
//}
#endif
  Vector<double,double> vai(len-1,HUGE_VAL);
  M.solve(b,vai);
  double sum=0.;
  double big=-HUGE_VAL;
  for (int n=0;n<len-1;n++) {
    double vain=(vai)[n];
    sum+=vai[n];
    big=max(big,abs(vain));
  }
  if (abs(sum)>=sqrt(DBL_EPSILON)*big) vai.scal(1./sum);
  else vai.scal(1./big);
  values_at_integers=OPERATOR_NEW Signal(lp_lo,lp_hi-1);
  for (int n=lp_lo;n<lp_hi;n++) {
    values_at_integers->value(n)=vai[n-lp_lo];
  }
#ifdef DEBUG
//for (int n=lp_lo;n<lp_hi;n++) {
//  cout << "\tvalues_at_integers[ " << n << " ] = "
//       << values_at_integers->value(n) << endl;
//}
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ScalingFunction::~ScalingFunction() {
  if (values_at_integers) delete values_at_integers;
  values_at_integers=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Signal* ScalingFunction::values(int scale) const {
//TRACER_CALL(t,"ScalingFunction::values");
#ifdef DEBUG
//for (int n=lowpass->firstIndex();n<=lowpass->lastIndex();n++) {
//  cout << "\ts[ " << n << "] = " << lowpass->value(n) << endl;
//}
//for (int n=values_at_integers->firstIndex();
//n<=values_at_integers->lastIndex();n++) {
//  cout << "\tvalues_at_integers[ " << n << "] = "
//       << values_at_integers->value(n) << endl;
//}
#endif
  CHECK_TEST(scale>=1);
  int twotoscale=pow(2,scale);
  int vai_lo=values_at_integers->firstIndex();
  int vai_hi=values_at_integers->lastIndex();
  int v_lo=vai_lo*twotoscale;
  int v_hi=(vai_hi+1)*twotoscale-1;
  Signal *v=OPERATOR_NEW Signal(v_lo,v_hi);
  for (int n=vai_lo;n<=vai_hi;n++) {
    v->value(n*twotoscale)=values_at_integers->value(n);
  }
#ifdef DEBUG
//for (int n=v_lo;n<=v_hi;n+=twotoscale) {
//  cout << "\tv[ " << n << " ] = " << v->value(n) << endl;
//}
//cout << endl;
#endif
  int lp_lo=lowpass->firstIndex();
  int lp_hi=lowpass->lastIndex();
  int inc=twotoscale;
  for (int j=1;j<=scale;j++) {
#ifdef DEBUG
//  cout << "\n\n\tj,inc = " << j << " " << inc << endl;
#endif
    for (int n=v_lo+inc/2;n<=v_hi;n+=inc) {
      double sum=0.;
      int twon=2*n;
#ifdef DEBUG
//    cout << "\n\tn = " << n << endl;
#endif
      for (int m=lp_lo;m<=lp_hi;m++) {
        int mm=twon-m*twotoscale;
        if (mm>=v_lo && mm<=v_hi) {
#ifdef DEBUG
//        cout << "\tm,lowpass[" << m << "],v[" << mm << "] = " << m << " "
//             << lowpass->value(m) << " " << v->value(mm) << endl;
#endif
          sum+=lowpass->value(m)*v->value(mm);
        }
      }
      v->value(n)=2.*sum;
    }
#ifdef DEBUG
//  for (int n=lp_lo*twotoscale+inc/2;n<lp_hi*twotoscale;n+=inc) {
//    cout << "\tv[ " << n << " ] = " << v->value(n) << endl;
//  }
//  cout << endl;
#endif
    inc /= 2;
  }
  return v;
}
