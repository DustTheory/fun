#include "SplineFilterBank.H"
#include <cmath>
#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//L = product filter order
//N = lowpass filter order
SplineFilterBank::SplineFilterBank(int N,int L) {
//TRACER_CALL(t,"SplineFilterBank::SplineFilterBank");
  CHECK_TEST(N>=1);
  CHECK_TEST(L>=1);
  if ((N+L)%2==1) L++;
  int M=(N+L)/2;

  lowpass=OPERATOR_NEW Filter(0,N);

  BinomialCoefficient bc(max(max(N,L),2*M-2));
  double two_power=pow(.5,N);
  for (int n=0;n<=N;n++) lowpass->value(n)=two_power*bc.value(N,n);
#ifdef DEBUG
//for (int n=lowpass->firstIndex();n<=lowpass->lastIndex();n++) {
//  cout << "\tlowpass[" << n << "] = " << lowpass->value(n) << endl;
//}
#endif

  Signal alpha(0,L);
  two_power=pow(.5,L);
  for (int n=0;n<=L;n++) alpha.value(n)=two_power*bc.value(L,n);
  alpha.shift((L-N)/2);

  Filter beta_filter(1-M,M-1);
  for (int n=0;n<M;n++) {
    double sum=0.;
    for (int m=n;m<M;m++) {
      int twom=2*m;
      sum+=bc.value(M+m-1,m)*bc.value(twom,m-n)*pow(.5,twom);
    }
    if (n%2 != 0) sum=-sum;
    beta_filter.value(n)=sum;
    beta_filter.value(-n)=sum;
  }
#ifdef DEBUG
//for (int n=beta_filter.firstIndex();n<=beta_filter.lastIndex();n++) {
//  cout << "\tbeta[" << n << "] = " << beta_filter.value(n) << endl;
//}
#endif

  dual_lowpass=OPERATOR_NEW Filter(
    alpha.firstIndex()+beta_filter.firstIndex(),
    alpha.lastIndex()+beta_filter.lastIndex());
  beta_filter.operateOn(alpha,dual_lowpass->impulseResponse());
//dual_lowpass->replaceWithAdjoint();
#ifdef DEBUG
//for (int n=dual_lowpass->firstIndex();n<=dual_lowpass->lastIndex();n++) {
//  cout << "\tdual_lowpass[" << n << "] = " << dual_lowpass->value(n)
//    << endl;
//}
#endif

  highpass=OPERATOR_NEW Filter(*dual_lowpass);
  highpass->replaceWithAdjoint();
  int shif=highpass->firstIndex();
  highpass->impulseResponse().shift(shif);
  highpass->impulseResponse().reverse();
#ifdef DEBUG
//for (int n=highpass->firstIndex();n<=highpass->lastIndex();n++) {
//  cout << "\thighpass[" << n << "] = " << highpass->value(n)
//    << endl;
//}
#endif

  dual_highpass=OPERATOR_NEW Filter(*lowpass);
  dual_highpass->impulseResponse().shift(-shif);
  dual_highpass->impulseResponse().reverse();
  dual_highpass->replaceWithAdjoint();
#ifdef DEBUG
//for (int n=dual_highpass->firstIndex();n<=dual_highpass->lastIndex();n++)
//{
//  cout << "\tdual_highpass[" << n << "] = " << dual_highpass->value(n)
//    << endl;
//}
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SplineFilterBank::analyze(const Signal &s,Signal &low_out,
Signal &high_out) const {
//TRACER_CALL(t,"SplineFilterBank::analyze");
#ifdef DEBUG
//for (int n=s.firstIndex();n<=s.lastIndex();n++) {
//  cout << "\ts[ " << n << " ] = " << s.value(n) << endl;
//}
#endif
  dual_lowpass->operateAdjointOnAndDownsample(s,low_out);
  dual_highpass->operateAdjointOnAndDownsample(s,high_out);
#ifdef DEBUG
//for (int n=low_out.firstIndex();n<=low_out.lastIndex();n++) {
//  cout << "\tlow_out[ " << n << " ] = " << low_out.value(n) << endl;
//}
//for (int n=high_out.firstIndex();n<=high_out.lastIndex();n++) {
//  cout << "\thigh_out[ " << n << " ] = " << high_out.value(n) << endl;
//}
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void SplineFilterBank::synthesize(const Signal &low_in,
const Signal &high_in,Signal &out) const {
//TRACER_CALL(t,"SplineFilterBank::synthesize");
#ifdef DEBUG
//for (int n=low_in.firstIndex();n<=low_in.lastIndex();n++) {
//  cout << "\tlow_in[ " << n << " ] = " << low_in.value(n) << endl;
//}
//for (int n=high_in.firstIndex();n<=high_in.lastIndex();n++) {
//  cout << "\thigh_in[ " << n << " ] = " << high_in.value(n) << endl;
//}
#endif
  out=0.;
  lowpass->addUpsampleAndOperateOn(low_in,out);
  highpass->addUpsampleAndOperateOn(high_in,out);
  out*=2.;
#ifdef DEBUG
//for (int n=out.firstIndex();n<=out.lastIndex();n++) {
//  cout << "\tout[ " << n << " ] = " << out.value(n) << endl;
//}
#endif
}
