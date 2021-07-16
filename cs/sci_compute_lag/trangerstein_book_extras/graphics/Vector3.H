//**********************************************************************
// Copyright 2006 John A. Trangenstein
//
// This software is made available for research and instructional use 
// only. 
// You may copy and use this software without charge for these 
// non-commercial purposes, provided that the copyright notice and 
// associated text is reproduced on all copies.  
// For all other uses (including distribution of modified versions), 
// please contact the author at
//   John A. Trangenstein
//   Department of Mathematics
//   Duke University
//   Durham, NC 27708-0320
//   USA
// or
//   johnt@math.duke.edu
// 
// This software is made available "as is" without any assurance that it
// is completely correct, or that it will work for your purposes.  
// Use the software at your own risk.
//**********************************************************************
#ifndef vECTOR3_H
#define vECTOR3_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/Vector3.H,v 1.1 2009/08/20 17:31:47 johnt Exp $"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include "Const.H"
#include "DDecl.H"
//#include "Errors.H"
using namespace std;
# ifndef ASSERT
# define ASSERT(s) { if (!(s)) { cerr << "assertion " << __STRING(s) << " failed in file " << __FILE__ << " at line " << __LINE__ << endl; abort(); } }
# endif
# ifndef OBSOLETE
# define OBSOLETE(s) { cerr << s << " in file " << __FILE__ << " at line " << __LINE__ << endl; abort(); }
# endif
# ifdef DEBUG
#   define CHECK_TEST(s) { ASSERT(s) }
#   define CHECK_POINTER(p) { ASSERT(p!=0) }
#   define CHECK_NONNEGATIVE(r) { ASSERT(r>=0) }
#   define CHECK_POSITIVE(r) { ASSERT(r>0) }
#   define CHECK_BOUNDS(v,l,u) { ASSERT(v>=l && v<u) }
#   define CHECK_NONZERO(q) { ASSERT(q!=0) }
#   define CHECK_SAME(v1,v2) { ASSERT((v1)==(v2)) }
# else
#   define CHECK_TEST(s) ;
#   define CHECK_POINTER(p) ;
#   define CHECK_NONNEGATIVE(r) ;
#   define CHECK_POSITIVE(r) ;
#   define CHECK_BOUNDS(v,l,u) ;
#   define CHECK_NONZERO(q) ;
#   define CHECK_SAME(v1,v2) ;
# endif

class Vector3 {
  friend ifstream& operator>>(ifstream&,Vector3&);
  private:
     double vect[3];
  public:
    Vector3() { vect[0]=ZERO;vect[1]=ZERO;vect[2]=ZERO; }
    Vector3(const double *a) { vect[0]=a[0];vect[1]=a[1];vect[2]=a[2]; }
    Vector3(double x,double y,double z) { vect[0]=x;vect[1]=y;vect[2]=z; }
    void getVect(double *v) const { 
      v[0]=vect[0];v[1]=vect[1];v[2]=vect[2];
    }
    double& operator[](int d) { 
      CHECK_BOUNDS(d,0,3)
      return vect[d];
    }
    double operator[](int d) const { 
      CHECK_BOUNDS(d,0,3)
      return vect[d];
    }
    double min() const {
      return std::min(vect[0],std::min(vect[1],vect[2]));
    }
    double max() const {
      return std::max(vect[0],std::max(vect[1],vect[2]));
    }
    double min(int &d) const { 
      d = ( (vect[0]<=vect[1]) ? ((vect[0]<=vect[2])?0:2)
                               : ((vect[1]<=vect[2])?1:2) );
      return vect[d];
    }
    double max(int &d) const { 
      d = ( (vect[0]>=vect[1]) ? ((vect[0]>=vect[2])?0:2)
                               : ((vect[1]>=vect[2])?1:2) );
      return vect[d];
    }

    Vector3 operator+(double p) const { 
      return Vector3(vect[0]+p,vect[1]+p,vect[2]+p);
    }
    Vector3 operator-(double p) const { 
      return Vector3(vect[0]-p,vect[1]-p,vect[2]-p);
    }
    Vector3 operator*(double p) const { 
      return Vector3(vect[0]*p,vect[1]*p,vect[2]*p); 
    }
    Vector3 operator/(double p) const { 
      CHECK_NONZERO(p)
      p=ONE/p; return Vector3(vect[0]*p,vect[1]*p,vect[2]*p); 
    }
    Vector3& operator+=(double p) { 
      vect[0] += p;vect[1] += p;vect[2] += p; return *this; 
    }
    Vector3& operator-=(double p) { 
      vect[0] -= p;vect[1] -= p;vect[2] -= p; return *this; 
    }
    Vector3& operator*=(double p) { 
      vect[0] *= p;vect[1] *= p;vect[2] *= p; return *this; 
    }
    Vector3& operator/=(double p) { 
      CHECK_NONZERO(p)
      p=ONE/p; vect[0] *= p;vect[1] *= p;vect[2] *= p; return *this; 
    }
    bool operator==(const Vector3& p) const {
      return vect[0]==p.vect[0]&&vect[1]==p.vect[1]
                               &&vect[2]==p.vect[2];
    }
    bool operator!=(const Vector3& p) const {
      return vect[0]!=p.vect[0]||vect[1]!=p.vect[1]
                               ||vect[2]!=p.vect[2];
    }
    bool operator<(const Vector3& p) const {
      return vect[0]<p.vect[0]&&vect[1]<p.vect[1]
                              &&vect[2]<p.vect[2];
    }
    bool operator<=(const Vector3& p) const {
      return vect[0]<=p.vect[0]&&vect[1]<=p.vect[1]
                               &&vect[2]<=p.vect[2];
    }
    bool operator>(const Vector3& p) const {
      return vect[0]>p.vect[0]&&vect[1]>p.vect[1]
                              &&vect[2]>p.vect[2];
    }
    bool operator>=(const Vector3& p) const {
      return vect[0]>=p.vect[0]&&vect[1]>=p.vect[1]
                               &&vect[2]>=p.vect[2];
    }
    Vector3 operator+(const Vector3 &p) const {
      return Vector3(vect[0]+p.vect[0],vect[1]+p.vect[1],
                                       vect[2]+p.vect[2] );
    }
    Vector3 operator-(const Vector3 &p) const {
      return Vector3(vect[0]-p.vect[0],vect[1]-p.vect[1],
                                       vect[2]-p.vect[2] );
    }
    Vector3 operator*(const Vector3 &p) const {
      return Vector3(vect[0]*p.vect[0],vect[1]*p.vect[1],
                                       vect[2]*p.vect[2]);
    }
    Vector3 operator/(const Vector3 &p) const {
      CHECK_NONZERO(p.vect[0])
      CHECK_NONZERO(p.vect[1])
      CHECK_NONZERO(p.vect[2])
      return Vector3(vect[0]/p.vect[0],vect[1]/p.vect[1],
                                       vect[2]/p.vect[2]);
    }
    Vector3& operator+=(const Vector3 &p) {
      vect[0] += p.vect[0];vect[1] += p.vect[1]; vect[2] += p.vect[2];
      return *this;
    }
    Vector3& operator-=(const Vector3 &p) {
      vect[0] -= p.vect[0];vect[1] -= p.vect[1]; vect[2] -= p.vect[2];
      return *this;
    }
    Vector3& operator*=(const Vector3 &p) {
      vect[0] *= p.vect[0];vect[1] *= p.vect[1];vect[2] *= p.vect[2]; 
      return *this;
    }
    Vector3& operator/=(const Vector3 &p) {
      CHECK_NONZERO(p.vect[0])
      CHECK_NONZERO(p.vect[1])
      CHECK_NONZERO(p.vect[2])
      vect[0] /= p.vect[0];vect[1] /= p.vect[1];vect[2] /= p.vect[2]; 
      return *this;
    }
    Vector3& min(const Vector3 &p) {
      vect[0] = std::min(vect[0],p.vect[0]);
      vect[1] = std::min(vect[1],p.vect[1]);
      vect[2] = std::min(vect[2],p.vect[2]); 
      return *this;
    }
    Vector3& max(const Vector3 &p) {
      vect[0] = std::max(vect[0],p.vect[0]);
      vect[1] = std::max(vect[1],p.vect[1]);
      vect[2] = std::max(vect[2],p.vect[2]); 
      return *this;
    }

    void printOn(ostream& = cerr) const;
};
inline ostream& operator<<(ostream &os,const Vector3 &p) {
  os << '<' << p[0] << ',' << p[1] << ',' << p[2] << '>';
  return os;
}
#define STRIP(c) while(is.get() != c)
inline istream& operator>>(istream &is,Vector3 &p) {
  STRIP('<');
  is >> p[0]; STRIP(','); is >> p[1]; STRIP(','); is >> p[2];
  STRIP('>');
  return is;
}
#undef STRIP
inline ofstream& operator<<(ofstream &os,const Vector3 &p) {
  double v[3]; p.getVect(v);
  os.write(reinterpret_cast<char*>(v), 3*sizeof(double));
  return os;
}
inline ifstream& operator>>(ifstream &is,Vector3 &p) {
  is.read(reinterpret_cast<char*>(p.vect), 3*sizeof(double));
  return is;
}
inline void Vector3::printOn(ostream &os) const {
  os << "Vector3: " << *this << '\n';
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline double dot(const Vector3 &v1,const Vector3 &v2) {
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}
inline double sumSquares(const Vector3 &v) {
  return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
}
inline double norm(const Vector3 &v) { return sqrt(sumSquares(v)); }
inline Vector3 operator-( const Vector3 &v ) { 
  return Vector3(-v[0],-v[1],-v[2]);
}
inline Vector3 min2(const Vector3 &p1, const Vector3 &p2) {
  Vector3 p(p1); p.min(p2); return p;
}
inline Vector3 max2(const Vector3 &p1, const Vector3 &p2) {
  Vector3 p(p1); p.max(p2); return p;
}
inline Vector3 cross(const Vector3 &v1,const Vector3 &v2) {
  Vector3 c; 
  c[0]=v1[1]*v2[2]-v1[2]*v2[1];
  c[1]=v1[2]*v2[0]-v1[0]*v2[2]; 
  c[2]=v1[0]*v2[1]-v1[1]*v2[0];
  return c;
}
#endif
