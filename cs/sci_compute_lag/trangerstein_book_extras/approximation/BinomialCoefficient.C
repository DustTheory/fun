#include "BinomialCoefficient.H"
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
BinomialCoefficient::BinomialCoefficient(int ord) : data(0),coefs(0),
order(ord) {
//TRACER_CALL(t,"BinomialCoefficient::BinomialCoefficient");
  if (order>1) {
    data=OPERATOR_NEW_BRACKET(int,numberStoredForAll(order));
    coefs=OPERATOR_NEW_BRACKET(int*,order-1);
    int pos=0;
    for (int n=2; n<=order; n++) {
      coefs[n-2]=data+pos;
      pos+=numberStoredForOne(n);
      coefs[n-2][0] = n;
      if (n==2) continue;
      for (int k=2;k<n-k;k++) {
        coefs[n-2][k-1]=coefs[n-3][k-2]+coefs[n-3][k-1];
      }
      int k=n/2;
      if (2*k==n) coefs[n-2][k-1]=2*coefs[n-3][k-2];
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
BinomialCoefficient::BinomialCoefficient(const BinomialCoefficient &b) :
data(0),coefs(0),order(b.order) {
//TRACER_CALL(t,"BinomialCoefficient::BinomialCoefficient");
  if (order>=2) {
    int n_coefs=numberStoredForAll(order);
    data=OPERATOR_NEW_BRACKET(int,n_coefs);
    memcpy(data,b.data,n_coefs*sizeof(int));
    coefs=OPERATOR_NEW_BRACKET(int*,order-1);
    int pos=0;
    for (int n=2;n<=order;n++) {
      coefs[n-2]=data+pos;
      pos+=numberStoredForOne(n);
    }
  }
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
BinomialCoefficient::~BinomialCoefficient() {
//TRACER_CALL(t,"BinomialCoefficient::~BinomialCoefficient");
  if (coefs) delete [] coefs; coefs=0;
  if (data) delete [] data; data=0;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
BinomialCoefficient& BinomialCoefficient::operator=(
const BinomialCoefficient &b) {
//TRACER_CALL(t,"BinomialCoefficient::operator=");
  if (order==b.order) return *this;
  if (order>=2) {
    if (coefs) delete [] coefs; coefs=0;
    if (data) delete [] data; data=0;
  } 
  order=b.order;
  if (order>=2) {
    int n_coefs=numberStoredForAll(order);
    data=OPERATOR_NEW_BRACKET(int,n_coefs);
    memcpy(data,b.data,n_coefs*sizeof(int));
    coefs=OPERATOR_NEW_BRACKET(int*,order-1);
    int pos=0;
    for (int n=2;n<=order;n++) {
      coefs[n-2]=data+pos;
      pos+=numberStoredForOne(n);
    }
  } 
  return *this;
}
