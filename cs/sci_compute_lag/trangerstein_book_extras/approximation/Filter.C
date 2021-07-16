#include "Filter.H"
#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Filter::operateOn(const Signal &s,Signal &r) const {
//TRACER_CALL(t,"Filter::operateOn");
  int signal_start=s.firstIndex();
  int signal_finish=s.lastIndex();
  int impulse_start=impulse_response->firstIndex();
  int impulse_finish=impulse_response->lastIndex();
  int result_start=signal_start+impulse_start;
  int result_finish=signal_finish+impulse_finish;
  CHECK_TEST(result_start>=r.firstIndex());
  CHECK_TEST(result_finish<=r.lastIndex());

  for (int n=result_start;n<=result_finish;n++) {
    double sum=0.;
    int m_start=max(signal_start,n-impulse_finish);
    int m_finish=min(signal_finish,n-impulse_start);
    for (int m=m_start;m<=m_finish;m++) {
      sum+=s.value(m)*impulse_response->value(n-m);
    }
    r.value(n)=sum;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Filter::operateAdjointOn(const Signal &s,Signal &r) const {
//TRACER_CALL(t,"Filter::operateAdjointOn");
  int signal_start=s.firstIndex();
  int signal_finish=s.lastIndex();
  int impulse_start=impulse_response->firstIndex();
  int impulse_finish=impulse_response->lastIndex();
  int result_start=signal_start+impulse_start;
  int result_finish=signal_finish+impulse_finish;
  CHECK_TEST(result_start>=r.firstIndex());
  CHECK_TEST(result_finish<=r.lastIndex());

  for (int n=result_start;n<=result_finish;n++) {
    double sum=0.;
    int m_start=max(signal_start,n+impulse_start);
    int m_finish=min(signal_finish,n+impulse_finish);
    for (int m=m_start;m<=m_finish;m++) {
      sum+=s.value(m)*impulse_response->value(m-n);
    }
    r.value(n)=sum;
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Filter::operateAdjointOnAndDownsample(const Signal &s,Signal &r) const
{
//TRACER_CALL(t,"Filter::operateAdjointOnAndDownsample");
  int signal_start=s.firstIndex();
  int signal_finish=s.lastIndex();
  int impulse_start=impulse_response->firstIndex();
  int impulse_finish=impulse_response->lastIndex();
  int result_start=(signal_start-impulse_finish)/2;
  int result_finish=(signal_finish-impulse_start)/2;
  CHECK_TEST(result_start>=r.firstIndex());
  CHECK_TEST(result_finish<=r.lastIndex());
#ifdef DEBUG
//for (int n=impulse_response->firstIndex();
//n<=impulse_response->lastIndex();n++) {
//  cout << "\timpulse_response[ " << n << "] = "
//       << impulse_response->value(n) << endl;
//}
//for (int n=s.firstIndex();n<=s.lastIndex();n++) {
//  cout << "\ts[ " << n << "] = " << s.value(n) << endl;
//}
#endif

  for (int n=result_start;n<=result_finish;n++) {
    double sum=0.;
    int twon=2*n;
    int m_start=max(signal_start,twon+impulse_start);
    int m_finish=min(signal_finish,twon+impulse_finish);
#ifdef DEBUG
//  cout << "\tn,m_start,m_finish = " << n << " " << m_start << " "
//       << m_finish << endl;
#endif
    for (int m=m_start;m<=m_finish;m++) {
      sum+=s.value(m)*impulse_response->value(m-twon);
    }
    r.value(n)=sum;
  }
#ifdef DEBUG
//for (int n=r.firstIndex();n<=r.lastIndex();n++) {
//  cout << "\tr[ " << n << "] = " << r.value(n) << endl;
//}
#endif
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Filter::addUpsampleAndOperateOn(const Signal &s,Signal &r) const {
//TRACER_CALL(t,"Filter::addUpsampleAndOperateOn");
#ifdef DEBUG
//cout << "\ts = " << s.firstIndex() << " " << s.lastIndex() << endl;
//cout << "\tr = " << r.firstIndex() << " " << r.lastIndex() << endl;
#endif
  int signal_start=s.firstIndex();
  int signal_finish=s.lastIndex();
  int impulse_start=impulse_response->firstIndex();
  int impulse_finish=impulse_response->lastIndex();
  int result_start=impulse_start+2*signal_start;
  int result_finish=impulse_finish+2*signal_finish;
#ifdef DEBUG
//cout << "\timpulse = " << impulse_start << " " << impulse_finish << endl;
//cout << "\tresult = " << result_start << " " << result_finish << endl;
#endif
  CHECK_TEST(result_start>=r.firstIndex());
  CHECK_TEST(result_finish<=r.lastIndex());
#ifdef DEBUG
//for (int n=impulse_response->firstIndex();
//n<=impulse_response->lastIndex();n++) {
//  cout << "\timpulse_response[ " << n << "] = "
//       << impulse_response->value(n) << endl;
//}
//for (int n=s.firstIndex();n<=s.lastIndex();n++) {
//  cout << "\ts[ " << n << "] = " << s.value(n) << endl;
//  cout << "\tn = " << n << endl;
//}
//for (int n=r.firstIndex();n<=r.lastIndex();n++) {
//  cout << "\tr[ " << n << "] = " << r.value(n) << endl;
//}
#endif

  for (int n=result_start;n<=result_finish;n++) {
    double sum=0.;
    int m_start=max(signal_start,(n-impulse_finish)/2);
    if (2*m_start<n-impulse_finish) m_start++;
    int m_finish=min(signal_finish,(n-impulse_start)/2);
    if (2*m_finish>n-impulse_start) m_finish--;
    for (int m=m_start;m<=m_finish;m++) {
      sum+=s.value(m)*impulse_response->value(n-2*m);
    }
    r.value(n)+=sum;
  }
#ifdef DEBUG
//for (int n=r.firstIndex();n<=r.lastIndex();n++) {
//  cout << "\tr[ " << n << "] = " << r.value(n) << endl;
//}
#endif
}
