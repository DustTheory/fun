#include "Wavelet.H"
#include <cmath>
//#include <float.h>
//#include <string.h>
#include "MemoryDebugger.H"
#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Signal* Wavelet::values(int scale) const {
//TRACER_CALL(t,"Wavelet::values");
  CHECK_TEST(scale>=0);
  int twotoscale=pow(2,scale);
  int hp_lo=highpass->firstIndex();
  int hp_hi=highpass->lastIndex();
#ifdef DEBUG
//for (int n=hp_lo;n<=highpass->lastIndex();n++) {
//  cout << "\ts[ " << n << "] = " << highpass->value(n) << endl;
//}
#endif

  Signal *sf_v=sf->values(scale);
  int sfv_lo=sf_v->firstIndex();
  int sfv_hi=sf_v->lastIndex();
#ifdef DEBUG
//for (int n=sf_v->firstIndex();n<=sf_v->lastIndex();n++) {
//  cout << "\tsf_v[ " << n << "] = " << sf_v->value(n) << endl;
//}
#endif
  int v_lo=(sfv_lo+hp_lo*twotoscale)/2;
  int v_hi=(sfv_hi+hp_hi*twotoscale)/2-1;
  Signal *v=OPERATOR_NEW Signal(v_lo,v_hi);
#ifdef DEBUG
//for (int n=0;n< sf_v->getNumber();n++) {
//  cout << "\tsf_v[ " << n << "] = " << (*sf_v)[n] << endl;
//}
#endif
  for (int n=v->firstIndex();n<=v->lastIndex();n++) {
    int twon=2*n;
    double sum=0.;
    for (int m=hp_lo;m<=hp_hi;m++) {
      int mm=twon-m*twotoscale;
      if (mm>=sfv_lo && mm<=sfv_hi) {
#ifdef DEBUG
//      cout << "\ts[" << m << "],v[" << mm << "] = " << highpass->value(m)
//           << " " << sf_v->value(mm) << endl;
#endif
        sum+=highpass->value(m)*sf_v->value(mm);
      }
    }
    v->value(n)=2.*sum;
  }
  delete sf_v; sf_v=0;
  return v;
}
