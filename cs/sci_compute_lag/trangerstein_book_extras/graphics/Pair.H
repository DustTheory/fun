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
#ifndef PAIR_H
#define PAIR_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/Pair.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#include <iostream>
#include <fstream>
#include <math.h>
#include "Const.H"
//#include "Errors.H"
#include "Vector3.H"

template<class T> class Pair {
  private:
    T vect[2];
  public:
    inline Pair() { vect[0]=ZERO; vect[1]=ZERO; }
    inline Pair(const T *a) { vect[0]=a[0]; vect[1]=a[1]; }
    inline Pair(T x,T y) { vect[0]=x; vect[1]=y; }
    inline virtual ~Pair() {;}
    inline T& operator[](int d) { 
      CHECK_BOUNDS(d,0,2)
      return vect[d];
    }
    inline T operator[](int d) const { 
      CHECK_BOUNDS(d,0,2)
      return vect[d];
    }
    inline Pair& operator/=(T t) { 
      vect[0]/=t; 
      vect[1]/=t; 
      return *this; 
    }
    inline bool operator==(const Pair& p) const {
      return vect[0]==p.vect[0] && vect[1]==p.vect[1];
    }
    inline T sumSquares() const { 
      return vect[0]*vect[0]+vect[1]*vect[1]; 
    }
    inline T norm() const { return sqrt(sumSquares()); }
    inline Vector3 projectToSphere(T radius) const {
      T d2=sumSquares(); T d=sqrt(d2); T t=radius*M_SQRT1_2;
      if (d<t) return Vector3(vect[0],vect[1],sqrt(radius*radius-d2));
      else return Vector3(vect[0],vect[1],t*t/d); 
    }
    virtual void printOn(ostream &os) const {
      os << "Pair:" << endl;
    }
};
template<class T> Pair<T> operator-(const Pair<T> &p,const Pair<T> &q) {
  Pair<T> dif=p; dif[0]-=q[0]; dif[1]-=q[1]; return dif;
}
template<class T> T innerProduct(const Pair<T> &p,const Pair<T> &q) {
  return p[0]*q[0]+p[1]*q[1];
}
template<class T> ostream& operator<<(ostream &os,const Pair<T> &p) {
  os << "<" << p[0] << "," << p[1] << ">";
  return os;
}
template<class T> Pair<T> closestInteriorPosition(
const Pair<T> p0,const Pair<T> p1,const Pair<T> p2,const Pair<T> &p){
  Pair<T> dif1=p1-p0;
  Pair<T> dif2=p2-p0;
  Pair<T> dif=p-p0;
  T area=dif1[0]*dif2[1]-dif1[1]*dif2[0];
  T t1=ZERO,t2=ZERO;
  if (area!=ZERO) {
    t1=(dif2[1]*dif[0]-dif2[0]*dif[1])/area;
    t2=(dif[1]*dif1[0]-dif[0]*dif1[1])/area;
  } else {
    T n1=dif1.norm();
    T n2=dif2.norm();
    if (n1>n2) t1=(dif1[0]*dif[0]+dif1[1]*dif[1])/n1;
    else t2=(dif2[0]*dif[0]+dif2[1]*dif[1])/n1;
  } 
  return Pair<T>(max(ZERO,min(ONE,t1)) , max(ZERO,min(ONE,t2)) );
}

template<class T> class HomoCoord {
  private:
    T coord[4];
  public:
    HomoCoord(T x,T y,T z,T w) { 
      coord[0]=x; coord[1]=y; coord[2]=z; coord[3]=w;
    }
    T& operator[](int i) {
      CHECK_BOUNDS(i,0,4)
      return coord[i];
    }
    T operator[](int i) const {
      CHECK_BOUNDS(i,0,4)
      return coord[i];
    }
    void negate() { 
      coord[0]=-coord[0]; coord[1]=-coord[1];
      coord[2]=-coord[2]; coord[3]=-coord[3];
    }
    void printOn(ostream &os) const {
      os << "HomoCoord = <" << coord[0] << " " << coord[1] << " "
         << coord[2] << " " << coord[3] << ">" << endl;
    }
};
#endif
