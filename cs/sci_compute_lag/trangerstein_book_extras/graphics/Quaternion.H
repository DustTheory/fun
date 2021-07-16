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
#ifndef QUATERNION_H
#define QUATERNION_H "$Header: /home/faculty/johnt/cvs/deal_new/graphics/Quaternion.H,v 1.1 2009/08/20 17:31:46 johnt Exp $"
#include <iostream>
#include <fstream>
#include <math.h>
#include "Const.H"
//#include "Errors.H"
#include "Pair.H"
#include "Vector3.H"

#if (SPACEDIM>1)
/*
 * Quaternions have the form
 *   q = a*i + b*j + c*k +d
 * where
 *   i*i = j*j = k*k = -1
 *   i*j = k = -j*i
 *   j*k = i = -k*j
 *   k*i = j = -i*k
 * Unit Quaternions always obey:  a^2 + b^2 + c^2 + d^2 = 1.0
 *
 * Note: See the following for more information on quaternions:
 *
 * - Shoemake, K., Animating rotation with quaternion curves, Computer
 *   Graphics 19, No 3 (Proc. SIGGRAPH'85), 245-254, 1985.
 * - Pletinckx, D., Quaternion calculus as a basic tool in computer
 *   graphics, The Visual Computer 5, 2-13, 1989.
 */
template<class T> class Quaternion {
  private:
    T vect[4];
  public:
    Quaternion() { vect[0]=ZERO;vect[1]=ZERO;vect[2]=ZERO;vect[3]=ONE; }
    T sumSquares() const {
      return vect[0]*vect[0]+vect[1]*vect[1]
            +vect[2]*vect[2]+vect[3]*vect[3];
    }
    T norm() const { return sqrt(sumSquares()); }
    void normalize() { 
      T val=1./norm(); 
      vect[0] *= val; vect[1] *= val; vect[2] *= val; vect[3] *= val;
    }
    Quaternion(const T *a) { 
      vect[0]=a[0];vect[1]=a[1];vect[2]=a[2];vect[3]=a[3]; 
      normalize();
    }
    Quaternion(T x,T y,T z,T a) { 
      vect[0]=x; vect[1]=y; vect[2]=z; vect[3]=a;
      normalize();
    }
    Quaternion(Vector3 axis,T angle) { 
      angle *= 0.5; axis *= (sin(angle)/::norm(axis));
      vect[0]=axis[0]; vect[1]=axis[1]; vect[2]=axis[2]; 
      vect[3]=cos(angle);
    }

    T& operator[](int d) { 
      CHECK_BOUNDS(d,0,4)
      return vect[d];
    }
    T operator[](int d) const { 
      CHECK_BOUNDS(d,0,4)
      return vect[d];
    }
//  multiplication is not commutative
    Quaternion operator*(const Quaternion &p) const {
      Vector3 t1(vect); 
      Vector3 t2(p.vect); 
      Vector3 t(cross(t1,t2)+(t1*p.vect[3]+t2*vect[3]));
      return Quaternion(t[0],t[1],t[2],vect[3]*p.vect[3]-dot(t1,t2));
    }
    Quaternion inverse() const {
      return Quaternion(-vect[0],-vect[1],-vect[2],vect[3]);
    }
//  rotationMatrix related to Quaterion * Vector3 * Quaternion.inverse
    void rotationMatrix(GLfloat m[4][4]) const {
      m[0][0] = 1.0 - 2.0 * (vect[1] * vect[1] + vect[2] * vect[2]);
      m[0][1] = 2.0 * (vect[0] * vect[1] - vect[2] * vect[3]);
      m[0][2] = 2.0 * (vect[2] * vect[0] + vect[1] * vect[3]);
      m[0][3] = 0.0;

      m[1][0] = 2.0 * (vect[0] * vect[1] + vect[2] * vect[3]);
      m[1][1]= 1.0 - 2.0 * (vect[2] * vect[2] + vect[0] * vect[0]);
      m[1][2] = 2.0 * (vect[1] * vect[2] - vect[0] * vect[3]);
      m[1][3] = 0.0;

      m[2][0] = 2.0 * (vect[2] * vect[0] - vect[1] * vect[3]);
      m[2][1] = 2.0 * (vect[1] * vect[2] + vect[0] * vect[3]);
      m[2][2] = 1.0 - 2.0 * (vect[1] * vect[1] + vect[0] * vect[0]);
      m[2][3] = 0.0;

      m[3][0] = 0.0;
      m[3][1] = 0.0;
      m[3][2] = 0.0;
      m[3][3] = 1.0;
    }
    void printOn(ostream& os = cerr) const {
      os << "Quaternion: ";
      os << *this;
    }
};
template<class T> ostream& operator<<(ostream &os,
const Quaternion<T> &q) {
  os << "<" << q[0] << "," << q[1] << "," << q[2] << "," << q[3] << ">"
     << endl;
  return os;
}

#define TRACKBALLSIZE  (0.8)
template<class T> Quaternion<T> trackBall(Pair<T> p1,Pair<T> p2) {
  if (p1 == p2) return Quaternion<T>();

  Vector3 vp1=p1.projectToSphere(TRACKBALLSIZE);
  Vector3 vp2=p2.projectToSphere(TRACKBALLSIZE);
  float t = norm(vp1-vp2) / (2.0*TRACKBALLSIZE);
  if (t > 1.0) t = 1.0;
  if (t < -1.0) t = -1.0;
  return Quaternion<T>(cross(vp2,vp1),2.*asin(t));
}

#endif
#endif
