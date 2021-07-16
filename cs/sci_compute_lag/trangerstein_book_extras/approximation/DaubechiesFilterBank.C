#include "BinomialCoefficient.H"
#include "DaubechiesFilterBank.H"
#include "Errors.H"
#include "MemoryDebugger.H"
#include "SquareMatrix.H"
#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DaubechiesFilterBank::DaubechiesFilterBank(int k) {
//TRACER_CALL(t,"DaubechiesFilterBank::DaubechiesFilterBank");
  if (k==1) {
    lowpass=OPERATOR_NEW Filter(0,1);
    lowpass->value(0)=0.5;
    lowpass->value(1)=0.5;
    highpass=OPERATOR_NEW Filter(0,1);
    highpass->value(0)=0.5;
    highpass->value(1)=-0.5;
    return;
  }
  CHECK_TEST(k>=1);
  Filter spline_lowpass(0,k);

  BinomialCoefficient bc(max(k,2*k-2));
  double two_power=pow(.5,k);
  for (int n=0;n<=k;n++) spline_lowpass.value(n)=two_power*bc.value(k,n);
#ifdef DEBUG
//for (int n=spline_lowpass.firstIndex();n<=spline_lowpass.lastIndex();n++)
//{
//  cout << "\tspline_lowpass[" << n << "] = " << spline_lowpass.value(n)
//       << endl;
//}
#endif

  int cms=2*k-1;
  double *beta_array_data=OPERATOR_NEW_BRACKET(double,cms);
  double *beta_array=beta_array_data+(k-1);
  for (int n=0;n<k;n++) {
    double sum=0.;
    for (int m=n;m<k;m++) {
      int twom=2*m;
      sum+=bc.value(k+m-1,m)*bc.value(twom,m-n)*pow(.5,twom);
    }
    if (n%2 != 0) sum=-sum;
    beta_array[n]=sum;
    beta_array[-n]=sum;
  }
#ifdef DEBUG
//for (int n=1-k;n<=k-1;n++) {
//  cout << "\tbeta_array[" << n << "] = " << beta_array[n] << endl;
//}
#endif
  SquareMatrix<double,double> companion_matrix(cms-1,0.);
  for (int i=1;i<cms-1;i++) companion_matrix(i,i-1)=1.;
  double den=1./beta_array_data[cms-1];
  for (int j=0;j<cms-1;j++) {
    companion_matrix(0,j)=-beta_array_data[cms-2-j]*den;
  }
#ifdef DEBUG
//for (int i=0;i<companion_matrix.size(0);i++) {
//  cout << "\tcompanion_matrix[ " << i << " , * ] = ";
//  for (int j = 0;j<companion_matrix.size(1);j++) {
//    cout << companion_matrix(i,j) << " ";
//  }
//  cout << endl;
//}
#endif

  SquareMatrix<double,complex<double> > *V=0;
  SquareMatrix<double,complex<double> > *U=0;
  Vector<double,complex<double> > *roots=companion_matrix.eigenvalues(V,U);
#ifdef DEBUG
//for (int i=0;i<roots->size();i++) {
//  cout << "\troots[ " << i << " ] = " << (*roots)[i] << endl;
//}
#endif
  complex<double> *wanted_roots=OPERATOR_NEW_BRACKET(complex<double>,k-1);
  int count=0;
  for (int n=0;n<cms-1;n++) {
    if (abs((*roots)[n]) > 1.) {
      wanted_roots[count]=(*roots)[n];
      count++;
    }
  }
  delete roots; roots=0;
#ifdef DEBUG
//for (int i=0;i<k-1;i++) {
//  cout << "\twanted_roots[ " << i << " ] = " << wanted_roots[i] << endl;
//}
#endif
  CHECK_TEST(count==k-1);

  complex<double> *second=OPERATOR_NEW_BRACKET(complex<double>,k);
  second[0]=-wanted_roots[0];
  second[1]=1.;
#ifdef DEBUG
//cout << "\tk,second[0,1] = " << k << " " << second[0] << " "
//     << second[1] << endl;
#endif
  for (int n=1;n<k-1;n++) {
    complex<double> old=second[0];
    second[0]*=-wanted_roots[n];
#ifdef DEBUG
//cout << "\n\tn,second[0] = " << n << " " << second[0] << endl;
#endif
    for (int m=1;m<=n;m++) {
      complex<double> current=second[m];
      second[m]=old-current*wanted_roots[n];
      old=current;
#ifdef DEBUG
//    cout << "\tsecond[ " << m << " ] = " << second[m] << endl;
#endif
    }
    second[n+1]=old;
  }
  double factor=sqrt( abs(beta_array_data[0]/second[0].real() ) );
  if (second[0].real()<0.) factor=-factor;
  for (int n=0;n<k;n++) second[n]*=factor;
  delete [] beta_array_data; beta_array_data=0;
  delete [] wanted_roots; wanted_roots=0;
#ifdef DEBUG
//for (int i=0;i<=k-1;i++) {
//  cout << "\tsecond[ " << i << " ] = " << second[i] << endl;
//}
#endif
  lowpass=OPERATOR_NEW Filter(0,2*k-1);
  Filter second_filter(0,k-1);
  for (int n=0;n<k;n++) second_filter.value(n)=second[n].real();
  delete [] second;

  second_filter.operateOn(spline_lowpass.impulseResponse(),
    lowpass->impulseResponse());
#ifdef DEBUG
//for (int n=lowpass->firstIndex();n<=lowpass->lastIndex();n++) {
//  cout << "\tlowpass[" << n << "] = " << lowpass->value(n)
//    << endl;
//}
#endif

  highpass=OPERATOR_NEW Filter(*lowpass);
  highpass->impulseResponse().replaceWithAdjoint();
  highpass->impulseResponse().shift(1-2*k);
  highpass->impulseResponse().reverse();
#ifdef DEBUG
//for (int n=highpass->firstIndex();n<=highpass->lastIndex();n++) {
//  cout << "\thighpass[" << n << "] = " << highpass->value(n)
//    << endl;
//}
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void DaubechiesFilterBank::analyze(const Signal &s,Signal &low_out,
Signal &high_out) const {
//TRACER_CALL(t,"DaubechiesFilterBank::analyze");
#ifdef DEBUG
//for (int n=s.firstIndex();n<=s.lastIndex();n++) {
//  cout << "\ts[ " << n << " ] = " << s.value(n) << endl;
//}
#endif
  lowpass->operateAdjointOnAndDownsample(s,low_out);
  highpass->operateAdjointOnAndDownsample(s,high_out);
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
void DaubechiesFilterBank::synthesize(const Signal &low_in,
const Signal &high_in,Signal &out) const {
//TRACER_CALL(t,"DaubechiesFilterBank::synthesize");
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
