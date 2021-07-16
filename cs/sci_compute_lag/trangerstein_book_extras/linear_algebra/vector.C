#include "vector.H"
#include <math.h>
Vector::Vector(const Vector &x) : sz(x.sz), data(0) {
  data=new float[sz]; copy(x);
}
Vector& Vector::operator+=(const Vector &y) {
  ASSERT(sz==y.sz);
  for (int i=0;i<sz;i++) data[i]+=y.data[i]; return *this;
}
Vector& Vector::operator-=(const Vector &y) {
  ASSERT(sz==y.sz);
  for (int i=0;i<sz;i++) data[i]-=y.data[i]; return *this;
}
Vector& Vector::operator*=(float alpha) {
  for (int i=0;i<sz;i++) data[i]*=alpha; return *this;
}
Vector& Vector::operator/=(float alpha) {
  ASSERT(alpha != 0.);
  float temp=1./alpha; for (int i=0;i<sz;i++) data[i]*=temp;
  return *this;
}
Vector& Vector::operator+(const Vector &y) const {
  ASSERT(sz==y.sz);
  Vector *sum=new Vector(sz);
  for (int i=0;i<sz;i++) sum->data[i]=data[i]+y.data[i]; return *sum;
}
Vector& Vector::operator-(const Vector &y) const {
  ASSERT(sz==y.sz);
  Vector *dif=new Vector(sz);
  for (int i=0;i<sz;i++) dif->data[i]=data[i]-y.data[i]; return *dif;
}
Vector& Vector::operator*(float alpha) const {
  Vector *prod=new Vector(sz);
  for (int i=0;i<sz;i++) prod->data[i]=data[i]*alpha; return *prod;
}
Vector& Vector::operator/(float alpha) const {
  ASSERT(alpha != 0.);
  float temp=1./alpha; Vector *quot=new Vector(sz);
  for (int i=0;i<sz;i++) quot->data[i]=data[i]*temp; return *quot;
}
int Vector::amax() const {
  float amx=fabs(data[0]);
  int imx=0;
  for (int i=1;i<sz;i++) {
    float ai=fabs(data[i]); if (ai>amx) { imx=i; amx=ai; }
  }
  return imx;
}
float Vector::asum() const {
  float sum=0.; for (int i=0;i<sz;i++) sum+=fabs(data[i]); return sum;
}
void Vector::axpy(float alpha,const Vector &x) {
  ASSERT(sz==x.sz);
  for (int i=0;i<sz;i++) data[i]+=x.data[i]*alpha;
}
float Vector::dot(const Vector &x) const {
  ASSERT(sz==x.sz);
  float sum=0.; for (int i=0;i<sz;i++) sum += data[i]*x.data[i];
  return sum;
}
float Vector::nrm2() const {
  float sum=0; for (int i=0;i<sz;i++) sum += data[i]*data[i];
  return sqrt(sum);
}
void Vector::copy(const Vector &x) {
  if (sz != x.sz) { delete data; sz=x.sz; data=new float[sz]; }
  for (int i=0;i<sz;i++) data[i]=x.data[i];
}
void Vector::copyFrom(int n,Vector &x) const {
  int m=min(n,min(sz,x.sz));
  for (int i=0;i<m;i++) x.data[i]=data[i];
}
void Vector::fillWith(float alpha) {
  for (int i=0;i<sz;i++) data[i]=alpha;
}
void Vector::swap(Vector &x) {
  ASSERT(sz==x.sz);
  for (int i=0;i<sz;i++) {
    float temp=data[i]; data[i]=x.data[i]; x.data[i]=temp;
  }
}
void Vector::printOn(ostream &os) const {
  os << "[";
  for (int i=0;i<sz-1;i++) os << data[i] << ",";
  os << data[sz-1] << "]";
}
