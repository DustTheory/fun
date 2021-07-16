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
#ifndef DVECTOR_H
#define DVECTOR_H "$Header: /home/faculty/johnt/cvs/graphics/DVector.H,v 1.7 2007/01/10 13:18:54 johnt Exp $"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include "Const.H"
#include "DDecl.H"
#include "Errors.H"

class DVector {
  friend ifstream& operator>>(ifstream&,DVector&);
  private:
     REAL vect[SPACEDIM];
  public:
    DVector() { D_ASSIGN(vect[0] = ZERO,vect[1] = ZERO,vect[2] = ZERO) }
    DVector(const REAL *a) { 
      D_ASSIGN(vect[0] = a[0],vect[1] = a[1],vect[2]=a[2])
    }
    DVector(D_DECL(REAL x,REAL y,REAL z)) { 
      D_ASSIGN(vect[0]=x,vect[1]=y,vect[2]=z)
    }
    void getVect(REAL *v) const { 
      D_ASSIGN(v[0]=vect[0],v[1]=vect[1],v[2]=vect[2])
    }
    REAL& operator[](int d) { 
      CHECK_DIRECTION(d)
      return vect[d];
    }
    REAL operator[](int d) const { 
      CHECK_DIRECTION(d)
      return vect[d];
    }
#if (SPACEDIM == 1)
    REAL min() const { return vect[0]; }
    REAL max() const { return vect[0]; }
    REAL min(int &d) const { d=0; return vect[d]; }
    REAL max(int &d) const { d=0; return vect[d]; }
#endif
#if (SPACEDIM == 2)
    REAL min() const { return std::min(vect[0],vect[1]); }
    REAL max() const { return std::max(vect[0],vect[1]); }
    REAL min(int &d) const { 
      d = ( (vect[0]<=vect[1]) ? 0 : 1 );
      return vect[d];
    }
    REAL max(int &d) const { 
      d = ( (vect[0]>=vect[1]) ? 0 : 1 );
      return vect[d];
    }
#endif
#if (SPACEDIM == 3)
    REAL min() const { 
      return std::min(vect[0],std::min(vect[1],vect[2])); 
    }
    REAL max() const { 
      return std::max(vect[0],std::max(vect[1],vect[2])); 
    }
    REAL min(int &d) const { 
      d = ( (vect[0]<=vect[1]) ? ((vect[0]<=vect[2])?0:2)
                               : ((vect[1]<=vect[2])?1:2) );
      return vect[d];
    }
    REAL max(int &d) const { 
      d = ( (vect[0]>=vect[1]) ? ((vect[0]>=vect[2])?0:2)
                               : ((vect[1]>=vect[2])?1:2) );
      return vect[d];
    }
#endif


    DVector operator+(REAL p) const { 
      return DVector(D_DECL(vect[0]+p,vect[1]+p,vect[2]+p));
    }
    DVector operator-(REAL p) const { 
      return DVector(D_DECL(vect[0]-p,vect[1]-p,vect[2]-p));
    }
    DVector operator*(REAL p) const { 
      return DVector(D_DECL(vect[0]*p,vect[1]*p,vect[2]*p)); 
    }
    DVector operator/(REAL p) const { 
      CHECK_NONZERO(p)
      return DVector(D_DECL(vect[0]/p,vect[1]/p,vect[2]/p)); 
    }
    DVector& operator+=(REAL p) { 
      D_ASSIGN(vect[0] += p,vect[1] += p,vect[2] += p) 
      return *this; 
    }
    DVector& operator-=(REAL p) { 
      D_ASSIGN(vect[0] -= p,vect[1] -= p,vect[2] -= p) 
      return *this; 
    }
    DVector& operator*=(REAL p) { 
      D_ASSIGN(vect[0] *= p,vect[1] *= p,vect[2] *= p) 
      return *this; 
    }
    DVector& operator/=(REAL p) { 
      CHECK_NONZERO(p)
      D_ASSIGN(vect[0] /= p,vect[1] /= p,vect[2] /= p) 
      return *this; 
    }
    BOOLEAN operator==(const DVector& p) const {
      return D_TERM(vect[0]==p.vect[0],&&vect[1]==p.vect[1],\
                                       &&vect[2]==p.vect[2]);
    }
    BOOLEAN operator!=(const DVector& p) const {
      return D_TERM(vect[0]!=p.vect[0],||vect[1]!=p.vect[1],\
                                       ||vect[2]!=p.vect[2]);
    }
    BOOLEAN operator<(const DVector& p) const {
      return D_TERM(vect[0]<p.vect[0],&&vect[1]<p.vect[1],\
                                      &&vect[2]<p.vect[2]);
    }
    BOOLEAN operator<=(const DVector& p) const {
      return D_TERM(vect[0]<=p.vect[0],&&vect[1]<=p.vect[1],\
                                       &&vect[2]<=p.vect[2]);
    }
    BOOLEAN operator>(const DVector& p) const {
      return D_TERM(vect[0]>p.vect[0],&&vect[1]>p.vect[1],\
                                      &&vect[2]>p.vect[2]);
    }
    BOOLEAN operator>=(const DVector& p) const {
      return D_TERM(vect[0]>=p.vect[0],&&vect[1]>=p.vect[1],\
                                       &&vect[2]>=p.vect[2]);
    }
    DVector operator+(const DVector &p) const {
      return DVector( D_DECL(vect[0]+p.vect[0],vect[1]+p.vect[1],\
                                             vect[2]+p.vect[2]) );
    }
    DVector operator-(const DVector &p) const {
      return DVector( D_DECL(vect[0]-p.vect[0],vect[1]-p.vect[1],\
                                             vect[2]-p.vect[2]) );
    }
    DVector operator*(const DVector &p) const {
      return DVector(D_DECL(vect[0]*p.vect[0],vect[1]*p.vect[1],\
                                             vect[2]*p.vect[2]) );
    }
    DVector operator/(const DVector &p) const {
      D_TERM(CHECK_NONZERO(p.vect[0]),CHECK_NONZERO(p.vect[1]),\
                                      CHECK_NONZERO(p.vect[2]))
      return DVector(D_DECL(vect[0]/p.vect[0],vect[1]/p.vect[1],\
                                              vect[2]/p.vect[2]) );
    }
    DVector& operator+=(const DVector &p) {
      D_ASSIGN(vect[0] += p.vect[0],vect[1] += p.vect[1],\
                                    vect[2] += p.vect[2]) 
      return *this;
    }
    DVector& operator-=(const DVector &p) {
      D_ASSIGN(vect[0] -= p.vect[0],vect[1] -= p.vect[1],\
                                    vect[2] -= p.vect[2]) 
      return *this;
    }
    DVector& operator*=(const DVector &p) {
      D_ASSIGN(vect[0] *= p.vect[0],vect[1] *= p.vect[1],\
                                    vect[2] *= p.vect[2]) 
      return *this;
    }
    DVector& operator/=(const DVector &p) {
      D_TERM(CHECK_NONZERO(p.vect[0]),CHECK_NONZERO(p.vect[1]),\
                                      CHECK_NONZERO(p.vect[2]))
      D_ASSIGN(vect[0] /= p.vect[0],vect[1] /= p.vect[1],\
                                    vect[2] /= p.vect[2]) 
      return *this;
    }
    DVector& min(const DVector &p) {
      D_ASSIGN(vect[0] = std::min(vect[0],p.vect[0]),\
               vect[1] = std::min(vect[1],p.vect[1]),\
               vect[2] = std::min(vect[2],p.vect[2])) 
      return *this;
    }
    DVector& max(const DVector &p) {
      D_ASSIGN(vect[0] = std::max(vect[0],p.vect[0]),\
               vect[1] = std::max(vect[1],p.vect[1]),\
               vect[2] = std::max(vect[2],p.vect[2])) 
      return *this;
    }

    void printOn(ostream& = cerr) const;
};
inline ostream& operator<<(ostream &os,const DVector &p) {
#if (SPACEDIM==1)
  os << '<' << p[0] << '>'; 
#endif
#if (SPACEDIM==2)
  os << '<' << p[0] << ',' << p[1] << '>'; 
#endif
#if (SPACEDIM==3)
  os << '<' << p[0] << ',' << p[1] << ',' << p[2] << '>'; 
#endif
  return os;
}
#define STRIP(c) while(is.get() != c)
inline istream& operator>>(istream &is,DVector &p) {
  STRIP('<'); 
#if (SPACEDIM==1)
  is >> p[0]; 
#endif
#if (SPACEDIM==2)
  is >> p[0]; STRIP(','); is >> p[1]; 
#endif
#if (SPACEDIM==3)
  is >> p[0]; STRIP(','); is >> p[1]; STRIP(','); is >> p[2]; 
#endif
  STRIP('>');
  return is;
}
#undef STRIP
inline ofstream& operator<<(ofstream &os,const DVector &p) {
  REAL v[SPACEDIM]; p.getVect(v);
  os.write(reinterpret_cast<char*>(v), SPACEDIM*sizeof(REAL));
  return os;
}
inline ifstream& operator>>(ifstream &is,DVector &p) {
  is.read(reinterpret_cast<char*>(p.vect), SPACEDIM*sizeof(REAL));
  return is;
}
inline void DVector::printOn(ostream &os) const {
  os << "DVector: " << *this << '\n';
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
inline REAL dot(const DVector &v1,const DVector &v2) {
  return D_TERM(v1[0]*v2[0],+v1[1]*v2[1],+v1[2]*v2[2]);
}
inline REAL sumSquares(const DVector &v) {
  return D_TERM(v[0]*v[0],+v[1]*v[1],+v[2]*v[2]);
}
inline REAL norm(const DVector &v) { return sqrt(sumSquares(v)); }
inline DVector operator-( const DVector &v ) { 
  return DVector( D_DECL(-v[0],-v[1],-v[2]) );
}
#ifndef __DECCXX
inline DVector min(const DVector &p1, const DVector &p2) {
  OBSOLETE("use min2 to avoid conflict with min in <vector>");
  return p1;
}
inline DVector max(const DVector &p1, const DVector &p2) {
  OBSOLETE("use max2 to avoid conflict with max in <vector>");
  return p1;
}
#endif
inline DVector min2(const DVector &p1, const DVector &p2) {
  return DVector(D_DECL(std::min(p1[0],p2[0]),std::min(p1[1],p2[1]),\
                       std::min(p1[2],p2[2])));
}
inline DVector max2(const DVector &p1, const DVector &p2) {
  return DVector(D_DECL(std::max(p1[0],p2[0]),std::max(p1[1],p2[1]),\
                       std::max(p1[2],p2[2])));
}
#if (SPACEDIM==3)
  inline DVector cross(const DVector &v1,const DVector &v2) {
    DVector c; 
    c[0]=v1[1]*v2[2]-v1[2]*v2[1];
    c[1]=v1[2]*v2[0]-v1[0]*v2[2]; 
    c[2]=v1[0]*v2[1]-v1[1]*v2[0];
    return c;
  }
#endif

#endif
