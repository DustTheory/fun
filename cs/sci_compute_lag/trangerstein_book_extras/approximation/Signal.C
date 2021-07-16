#include <string.h>
#include "MemoryDebugger.H"
#include "Signal.H"
#include "Tracer.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Signal::Signal(int f,int l) : data(0),first(f),last(l) {
//TRACER_CALL(t,"Signal::Signal");
  int len=last-first+1;
  CHECK_POSITIVE(len);
  data=OPERATOR_NEW_BRACKET(double,len);
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Signal::Signal(const Signal &b) : first(b.first), last(b.last) {
  int len=last-first+1;
  data=OPERATOR_NEW_BRACKET(double,len);
  memcpy(data,b.data,len*sizeof(double));
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Signal::~Signal() {
//TRACER_CALL(t,"Signal::~Signal");
  if (data) delete [] data; data=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Signal& Signal::operator=(double v) {
  int len=last-first+1;
  for (int n=0;n<len;n++) data[n]=v;
  return *this;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Signal& Signal::operator*=(double v) {
//TRACER_CALL(t,"Signal::operator*=");
#ifdef DEBUG
//cout << "\tv = " << v << endl;
#endif
  int len=last-first+1;
  for (int n=0;n<len;n++) data[n]*=v;
  return *this;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Signal::reverse() {
//TRACER_CALL(t,"Signal::reverse");
  int start=first;
  if (start%2==0) start++;
  for (int n=start;n<=last;n+=2) data[n-first]=-data[n-first];
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
void Signal::replaceWithAdjoint() {
//TRACER_CALL(t,"Signal::replaceWithAdjoint");
  int len=last-first+1;
  for (int n=0;n<len/2;n++) {
    int m=len-1-n;
    double dn=data[n];
    data[n]=data[m];
    data[m]=dn;
  }
  int n=first;
  first=-last;
  last=-n;
}
